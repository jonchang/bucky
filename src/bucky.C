// BUCKy 1.2 - Bayesian Untangling of Concordance Knots (applied to yeast and other organisms)
// Copyright (C) 2006-2007 by Bret Larget and Cecile Ane

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License, version 2, as
// published by the Free Software Foundation.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License (the file gpl.txt included with this
// distribution or http://www.gnu.org/licenses/gpl.txt) for more
// details.

// Version History

// First alpha version 7 March, 2005
// Version 1.1 released 30 October, 2006
// Version 1.2a 21 August, 2007
// Version 1.2b 17 January, 2008

// File:     bucky.C

// Input:    [options] <summary files>
//
// Options:  
//            [-a alpha-value]                --- Dirichlet process hyper-parameter
//            [-n number-of-MCMC-updates]     --- number of MCMC updates per chain
//            [-c number-of-chains]           --- number of separate chains for MCMCMC
//            [-r MCMCMC-rate]                --- rate in updates to propose MCMCMC swap
//            [-m alpha-multiplier]           --- multiple used to determine alpha for heated chains
//            [-s subsample-rate]             --- subsample rate for output file (summaries do not subsample)
//            [-o output-file-root]           --- root name of all output files
//            [-s1 seed1]                     --- 0 < seed1 < 4294967295
//            [-s2 seed2]                     --- 0 < seed1 < 4294967295
//            [--use-independence-prior]
//            [--calculate-pairs]
//            [--create-sample-file]
//            [--create-single-file]
//            [--create-joint-file]
//            [-h] prints the brief usage message
//            [--help] prints the brief usage message
//            [--version] prints brief intro and version number
//
// Output:   
//           
//  fileRoot.cluster       --- summary of cluster size (number of groups)
//  fileRoot.concordance   --- split-by-split summary of concordance factor
//  fileRoot.gene          --- gene by gene summary of joint posterior information
//  fileRoot.joint         --- table of raw second stage MCMC counts for joint posterior
//  fileRoot.input         --- names of input files
//  fileRoot.out           --- output file
//  fileRoot.pairs         --- output if pairwise analysis is desired (very slow!)
//  fileRoot.sample        --- output of raw sampled GTMs (very big!)
//  fileRoot.single        --- summary of table of single gene posterior probabilities

// Changes in version 1.1
// --- added probabilities to concordance file output
// --- fixed bug with split index that assumed no more than 8 taxa
// --- implemented new and better updateGroups; use it as the default
// --- added runtime argument ---no-use-update-groups

// Changes in version 1.2
// --- Added more information to the --help output including default values
// --- Used Huffman coding idea to create a binary tree for proposing trees that minimizes
//     the average number of comparisons necessary to map the uniform number to the tree
//     which should increase the speed of the program, especially for data sets with many poorly resolved trees.
// --- Added probabilities to the .cluster file.
// --- Changed function getProb which speeds up retrieving probabilities for genes (by a lot!) by using a lookup table.
// --- Fixed error in writing clades in primary concordance tree
// --- Concordance summary includes primary concordance tree
// --- Eliminated some output files, and do not print some output by default. Add controls to print more output.

// Changes in version 1.2b
// --- Added the distribution of genome-wide concordance factors for splits in the primary concordance tree
// --- Fixed error in updatePairCounts, count at the bottom-right
// --- Fixed uninitialized default values for createXxxFile, fixed error in setFileNames

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <limits>

#include "bucky.h"

using namespace std;
string VERSION = "1.2";
string DATE = "17 January 2008";
string COPYRIGHT = "Copyright (C) 2006-2008 by Bret Larget and Cecile Ane";

//
// readFile
//
// input file is expected to have each line contain two fields separated by white space
//   the first field is a topology string (parenthetic representation)
//   the second field is a positive weight (usually a count from an MCMC sample)
// add read topologies to row i of the table
// when new topologies are discovered, add a new column to the table

bool readFile(string filename, int i, vector<string> &topologies, vector<vector<double> > &table, int &max,
	     vector<int>& taxid, vector<string>& translateTable)
{
  ifstream f(filename.c_str());
  if(f.fail()) {
    cerr <<"Error: Cannot open file " << filename << "." << endl;
    exit(1);
  }

  // search for the gene's translate table
  bool hasTtable = false;
  string keyword;
  f >> keyword;
  if(keyword=="translate" || keyword=="TRANSLATE" || keyword=="Translate")
    hasTtable=true;
  f.close();

  int numTaxa=-1;

  f.open(filename.c_str());
  if (hasTtable){ 
  // read the gene's translate table
    string taxonname;
    int taxonnumber;
    f >> keyword; // take out 'translate'
    bool done=false;  // will be done when semicolon is reached
    while(f.good() && !f.eof() && !done){
      f >> taxonnumber >> taxonname;
      if(f.fail()){
	cerr << "Error in reading the translate table for gene "<<i<<endl;
	exit(1);
      }
      size_t pos=taxonname.rfind(",");
      if (pos!=string::npos)
	taxonname.replace(pos,1,"");
      else {
	pos = taxonname.rfind(";");
	if (pos!=string::npos){
	  taxonname.replace(pos,1,"");
	  done=true;
	}
	else {
	  f >> keyword; // take out the comma
	  if (keyword==";") done=true;
	}
      }
      int taxonid = 0;
      for (int i=0; i<translateTable.size(); i++)
	if (taxonname == translateTable[i]){
	  taxonid = i+1;
	  taxid.push_back(i+1);
	  break;
	}
      if (taxonid==0){ // new taxon
	translateTable.push_back(taxonname);
	taxid.push_back(translateTable.size());
      }
    }
    numTaxa=taxid.size();
  }

  // set zeros for ith row of table
  table[i].resize(topologies.size(),0);
  while(f.good() && !f.eof()) {
    string top;
    double weight;
    bool done;
    f >> top >> weight;
    if(f.fail())
      break;
    else
      done = false;

    // count taxa
    int countTaxa=0;
    if(!done) {
      istringstream s(top);
      char c;
      int z;
      do {
	c = s.peek();
	if(!isdigit(c))
	  s >> c;
	else {
	  s >> z;
	  countTaxa++;
	}
      } while(c != ';' && c != EOF && s.good());
      if(numTaxa == -1){ // no translate table was found and first tree was read
	numTaxa = countTaxa;
	taxid.resize(countTaxa);
	for (int i=0; i<countTaxa; i++)
	  taxid[i]=i+1; // assumes IDs are 1,2,...,Ntax.
	if (translateTable.size()==0){ // then build the overall translate table
	  for (int i=0; i<countTaxa; i++){
	    ostringstream s;
	    s << "taxon" << i+1;
	    translateTable.push_back(s.str());
	  }
	}
      }
      else if(countTaxa != numTaxa) {
	cerr << "Error: File " << filename << " with line " << top << " " << weight << " has " << countTaxa
			    << " taxa where " << numTaxa << " are expected." << endl;
	exit(1);
      }
      done = true;
    }

    // determine longest length of topology string
    if(top.length()>max)
      max = top.length();

    // replace taxon numbers by taxid in the topology
    if (hasTtable){
      istringstream itop(top);
      ostringstream otop;
      char c;
      int z;
      do {
	c = itop.peek();
	if(!isdigit(c)){
	  itop >> c;
	  otop << c;
	}
	else {
	  itop >> z;
	  otop << taxid[z-1];
	}
      } while(c != ';' && c != EOF && itop.good());
      top = otop.str();
    }

    // check if topology has already been read
    int n = find(topologies.begin(),topologies.end(),top) - topologies.begin();
    if(n==topologies.size()) {
      topologies.push_back(top);
      for(int k=0;k<=i;k++)
	table[k].push_back(0);
    }
    table[i][n] = weight;
  }
  f.close();
  return hasTtable;
}

double State::calculateLogPriorProb()
{
  logPriorProb=0;
  for(int i=0;i<numGroups;i++)
    logPriorProb += logA(alphaOverTop,counts[indices[i]]);
  logPriorProb -= logA(alpha,genes.size());
  return logPriorProb;
  
}

int State::updateOne(int i,Rand& rand) {
  int oldTop = tops[i];
  int newTop;
  if(i < genes.size())
    //    newTop = genes[i]->pickTree(rand);
    newTop = genes[i]->pickTreeFast(rand);
  else {
    cerr << "Internal Error: i = " << i << " too large in updateOne." << endl;
    exit(1);
  }

  if(newTop==oldTop)
    return 1;
  int n1 = counts[oldTop];
  int n2 = counts[newTop];
  double logHR=0;
  if(!useIndependencePrior) { // logHR=0 if independence prior
    logHR += log(n2 + alphaOverTop);
    logHR -= log(n1-1 + alphaOverTop);
  }
  if(log(rand.runif()) > logHR) // reject
      return 0;
  else { // accept
    tops[i] = newTop;
    counts[oldTop]--;
    counts[newTop]++;
    if(n2==0) {
      numGroups++;
      indices.push_back(newTop);
      sort(indices.begin(),indices.end());
    }
    if(n1==1) {
      numGroups--;
      int j=0;
      while(indices[j]!=oldTop)
	j++;
      indices.erase(indices.begin()+j,indices.begin()+j+1);
    }
    logPriorProb += logHR;
    if(i < genes.size())
      logPosteriorProbProduct += (log(genes[i]->getProb(newTop)) - log(genes[i]->getProb(oldTop)));
    return 1;
  }
}

/********************************************************************************

This is the old very slow algorithm.

void State::updateOneGroup(int group,Rand& rand) {
  // group is index in indices of the group to update
  // find list of topologies that:
  //  (1) are not topologies for a different group; and
  //  (2) have positive single gene prior for all genes in the group
  // current topology is in this set
  // pick one with prob. proportional to prod_{set} prob(top)

  int oldTop = indices[group];

  // find set of genes currently with oldTop
  vector<int> set(counts[oldTop]); // indices of genes with oldTop
  int j=0;
  for(int i=0;i<genes.size();i++)
    if(tops[i] == oldTop)
      set[j++] = i;

  // loop through topologies to find possible tops and weights
  vector<int> possibleTops(0);
  vector<double> logWeights(0);
  for(int top=0;top<numTreesSampled;top++) {
    if(top == oldTop) {
      possibleTops.push_back(top);
      double w = 0.0;
      for(int i=0;i<set.size();i++) {
	double p = genes[set[i]]->getProb(top);
	if(p>0)
	  w += log(p);
	else {
	  cerr << "Error: Current top has zero probability!" << endl;
	  exit(1);
	}
      }
      logWeights.push_back(w);
    }
    else {
      if(!(counts[top]>0)) { // top not already used by another group
	double w = 0.0;
	bool noZeroProbs = true;
	for(int i=0;i<set.size() && noZeroProbs;i++) {
	  double p = genes[set[i]]->getProb(top);
	  if(p>0)
	    w += log(p);
	  else
	    noZeroProbs = false;
	}
	if(noZeroProbs) {
	  possibleTops.push_back(top);
	  logWeights.push_back(w);
	}
      }
    }
  }

  if(possibleTops.size()==1) // oldTop is the only possible choice
    return;

  // subtract max from each log weight for numerical accuracy
  double m = logWeights[0];
  for(int i=1;i<logWeights.size();i++)
    if(logWeights[i] > m)
      m = logWeights[i];
  double totalWeight=0.0;
  for(int i=0;i<logWeights.size();i++) {
    logWeights[i] = exp(logWeights[i] - m); // set logWeights[i] to weight
    totalWeight += logWeights[i];
  }
  double x = rand.runif() * totalWeight;
  int newTopIndex = 0;
  x -= logWeights[newTopIndex];
  while(x > 0) {
    x -= logWeights[++newTopIndex];
  }
  int newTop = possibleTops[newTopIndex];
  if(newTop == oldTop) // picked oldTop
    return;
  // make changes to newTop
  for(int i=0;i<set.size();i++)
    tops[set[i]] = newTop;
  counts[newTop] = counts[oldTop];
  counts[oldTop] = 0;
  j=0;
  while(indices[j]!=oldTop)
    j++;
  indices.erase(indices.begin()+j,indices.begin()+j+1);
  indices.push_back(newTop);
  sort(indices.begin(),indices.end());
}

*******************************************************************************/

// New, faster update group algorithm
int State::updateOneGroup(int gene,Rand& rand) {
  // pick a random topology for the gene
  // check that it:
  //  (1) is not a topology for another group
  //  (2) has positive single gene prior for all genes in the group
  // return 1 if change accepted, 0 otherwise.

  int oldTop = tops[gene];
  int newTop = genes[gene]->pickTreeFast(rand);
  if(newTop==oldTop)
    return 0;

  // check if newTop is already in use
  if(counts[newTop]>0)
    return 0;

  // possible new topology for group
  vector<int> set(counts[oldTop]); // indices of genes with oldTop
  int j=0;
  double oldLogProb = 0.0;
  double newLogProb = 0.0;
  for(int i=0;i<genes.size();i++) {
    if(tops[i] == oldTop) {
      oldLogProb += log( genes[i]->getProb(oldTop) );
      set[j++] = i;
      double p = genes[i]->getProb(newTop);
	  if(p>0)
	    newLogProb += log(p);
	  else // acceptance probability is 0
	    return 0;
    }
  }

  if(log(rand.runif()) < newLogProb - oldLogProb + genes[gene]->getProb(oldTop) - genes[gene]->getProb(newTop)) {
    // make changes to newTop
    for(int i=0;i<set.size();i++)
      tops[set[i]] = newTop;
    counts[newTop] = counts[oldTop];
    counts[oldTop] = 0;
    j=0;
    while(indices[j]!=oldTop)
      j++;
    indices.erase(indices.begin()+j,indices.begin()+j+1);
    indices.push_back(newTop);
    sort(indices.begin(),indices.end());
    return 1;
  }
  else
    return 0;
}

int State::update(Rand& rand) {
  int sum=0;
  for(int i=0;i<genes.size();i++)
    sum += updateOne(i,rand);
  return sum;
}

void State::print(ostream& f) {
  f << "State:" << endl;
  f << "  alpha = " << alpha << endl;
  f << "  logAlpha = " << logAlpha << endl;
  f << "  numGroups = " << numGroups << endl;
  f << "  length of indices = " << indices.size() << endl; // numGroups
  f << "  logPriorProb = " << logPriorProb << endl;
  f << "  logPosteriorProbProduct = " << logPosteriorProbProduct << endl;
  f << "  logH = " << getLogH() << endl;
  f << "  counts:" << endl;
  f << "    Top Count" << endl;
  for(int i=0;i<numGroups;i++)
    f << setw(7) << indices[i] << setw(6) << counts[indices[i]] << endl;
  f << "  tops:" << endl;
  f << "   gene   top" << endl;
  for(int i=0;i<genes.size();i++)
    f << setw(7) << i << setw(6) << tops[i] << endl;
    
}

void State::updateSplits(vector<vector<int> >& totals,vector<vector<int> >& topSplitsIndex) {
  int numSplitsPerTree = topSplitsIndex[0].size();
  int numSplits = totals.size();
  vector<int> counts(numSplits);
  for(int i=0;i<numSplits;i++)
    counts[i] = 0;
  for(int i=0;i<genes.size();i++)
    for(int j=0;j<numSplitsPerTree;j++) {
      counts[topSplitsIndex[tops[i]][j]]++;
    }
  for(int i=0;i<numSplits;i++)
    totals[i][counts[i]]++;
}

void State::updatePairCounts(vector<vector<int> >& counts)
{
  // slow way to do this.  maybe we should have State maintain a vector of groups....
  // only fill in upper triangle (i <= j)
  for(int i=0;i<genes.size();i++)
    for(int j=i;j<genes.size();j++)
      if( tops[i] == tops[j] )
	counts[i][j]++;
}

Node* Node::getNeighbor(int i) const { return edges[i]->getOtherNode(this); }

int Node::setAllTaxa(Node* parent,unsigned int all) {
  int x,ep;
  Edge* e;
  if(leaf) {
    x = 1;
    for(int i=0;i<number;i++)
      x *= 2;
    taxa[0] = all - x;
    e = edges[0];
  }
  else {
    x = 0;
    for(int i=0;i<3;i++) {
      Node* n = getNeighbor(i); 
      if(n != parent) {
	int y = n->setAllTaxa(this,all);
	setTaxa(i,y);
	x += y;
      }
      else {
	e = edges[i];
	ep = i;
      }
    }
  }
  int s = (x%2 == 1 ? x : all - x);
  if(!leaf)
    setTaxa(ep,all-x);
  e->setSplit(s);
  return x;
}

void Tree::setAllTaxa() {
  Node* root = nodes[0]->getNeighbor(0);
  for(int i=0;i<3;i++)
    root->setTaxa(i,root->getNeighbor(i)->setAllTaxa(root,all));
}

void Node::print(ostream& f) const
{
  f << setw(3) << number << ":";
  f << " taxa = (" << taxa[0];
  if(!leaf)
    for(int i=1;i<3;i++)
      f << "," << taxa[i];
  f << "), ";
  for(int i=0;i<3;i++) {
    if(i>0 && !leaf)
      f << ",";
    if(i==0 || !leaf)
      f << " e[" << edges[i]->getNumber() << "] -> " << "n[" << getNeighbor(i)->getNumber() << "]";
  }
  f << endl;
}

void Edge::print(ostream& f,int numTaxa) const
{
  f << setw(3) << number << ":";
  f << " n[" << nodes[0]->getNumber() << "] <-> " << "n[" << nodes[1]->getNumber() << "], split = ";
  Split s(split,numTaxa);
  s.print(f);
  f << endl;
}

void Tree::print(ostream& f) const
{
  f << top << endl;
  f << "numTaxa = " << numTaxa << ", numNodes" << numNodes << ", numEdges = " << numEdges << endl;
  f << "Nodes:" << endl;
  for(int i=0;i<numNodes;i++)
    nodes[i]->print(f);
  f << "Edges:" << endl;
  for(int i=0;i<numEdges;i++)
    edges[i]->print(f,numTaxa);
}

void Tree::connect(string top)
{
  int currentInternalNode=numTaxa,currentLeafEdge=0,currentInternalEdge=numTaxa;
  Node *node1,*node2;
  char ch;
  int n;

  istringstream topStr(top);
  topStr >> ch;
  if(ch!='(') {
    topStr.putback(ch);
    topStr >> n;
  }
  else {
    node1 = connectInt(topStr,currentInternalNode,currentLeafEdge,currentInternalEdge);
    topStr >> ch;
    if(ch!=',') {
      cerr << "Error: Cannot parse topology string, expected ',', found " << ch << endl; 
     exit(1);
    }
    node2 = connectInt(topStr,currentInternalNode,currentLeafEdge,currentInternalEdge);
    topStr >> ch;
    if(ch!=')') {
      cerr << "Error: Cannot parse topology string, expected ')', found " << ch << endl;
      exit(1);
    }
    connectTwoNodes(node1,node2,currentLeafEdge,currentInternalEdge);
  }
}

Node* Tree::connectInt(istream& topStr,int& currentInternalNode,int& currentLeafEdge,int& currentInternalEdge)
{
  char ch;
  int n;
  topStr >> ch;
  if(ch!='(') {
    topStr.putback(ch);
    topStr >> n;
    return nodes[n-1];
  }
  else {
    Node *node1,*node2;
    node1 = connectInt(topStr,currentInternalNode,currentLeafEdge,currentInternalEdge);
    topStr >> ch;
    if(ch!=',') {
      cerr << "Error: Cannot parse topology string, expected ',', found " << ch << endl;
      exit(1);
    }
    node2 = connectInt(topStr,currentInternalNode,currentLeafEdge,currentInternalEdge);
    topStr >> ch;
    if(ch!=')') {
      cerr << "Error: Cannot parse topology string, expected ')', found " << ch << endl;
      exit(1);
    }
    return connectThreeNodes(node1,node2,nodes[currentInternalNode++],currentLeafEdge,currentInternalEdge);
  }
}

Node *Tree::connectThreeNodes(Node* node1,Node* node2,Node* node3,int& currentLeafEdge,int& currentInternalEdge)
{
  // node3 must be an internal node.
  Edge *edge1,*edge2;
  if(node1->isLeaf())
    edge1 = edges[currentLeafEdge++];
  else
    edge1 = edges[currentInternalEdge++];
  if(node2->isLeaf())
    edge2 = edges[currentLeafEdge++];
  else
    edge2 = edges[currentInternalEdge++];
  node1->setEdge(0,edge1);
  node2->setEdge(0,edge2);
  node3->setEdge(1,edge1);
  node3->setEdge(2,edge2);

  edge1->setNode(0,node1);
  edge1->setNode(1,node3);
  edge2->setNode(0,node2);
  edge2->setNode(1,node3);
  return node3;
}

void Tree::connectTwoNodes(Node* node1, Node* node2, int& currentLeafEdge, int& currentInternalEdge)
{
  Edge *edge;
  if(node1->isLeaf() || node2->isLeaf())
    edge = edges[currentLeafEdge++];
  else
    edge = edges[currentInternalEdge++];
  node1->setEdge(0,edge);
  node2->setEdge(0,edge);
  edge->setNode(0,node1);
  edge->setNode(1,node2);
}

void Tree::getSplits(SplitSet& s) {
    s.resize(numSplits);
    s.setAll(all);
    s.setNumTaxa(numTaxa);
    for(int i=0,j=numTaxa;i<numSplits;i++)
      s.setSplit(i,edges[j++]->getSplit());
  }

// assumes that numChains > 1 and useIndependencePrior == false
void mcmcmc(vector<State*>& states,vector<int>& index,vector<double>& alphas,Rand& rand,vector<int>& accepts, vector<int>& proposals)
{
  int numChains = accepts.size();
  int a = (int) ( (numChains-1)*rand.runif() );
  proposals[a]++;
  int b = a+1;
  int ia = index[a];
  int ib = index[b];
  double logAcceptProb = - states[ia]->getLogPriorProb() - states[ib]->getLogPriorProb();
  states[ia]->setAlpha(alphas[b]);
  states[ib]->setAlpha(alphas[a]);
  logAcceptProb += states[ia]->getLogPriorProb() + states[ib]->getLogPriorProb();
  if(log(rand.runif()) < logAcceptProb) { // accept
    index[a] = ib;
    index[b] = ia;
    accepts[a]++;
  }
  else { // reject
    states[ia]->setAlpha(alphas[a]);
    states[ib]->setAlpha(alphas[b]);
  }
}

void usage(Defaults defaults)
{
  cerr << "Usage: bucky [options] <input files>" << endl << endl;
  cerr << "  Options:" << endl;
  cerr << "  Parameter              | Usage                      | Default Value" << endl;
  cerr << "  -------------------------------------------------------------------" << endl;
  cerr << "  alpha                  | -a number                  | " << defaults.getAlpha() << endl; 
  cerr << "  # of runs              | -k integer                 | " << defaults.getNumRuns() << endl;
  cerr << "  # of MCMC updates      | -n integer                 | " << defaults.getNumUpdates() << endl;
  cerr << "  # of chains            | -c integer                 | " << defaults.getNumChains() << endl;
  cerr << "  MCMCMC Rate            | -r integer                 | " << defaults.getMCMCMCRate() << endl;
  cerr << "  alpha multiplier       | -m number                  | " << defaults.getAlphaMultiplier() << endl;
  cerr << "  subsample rate         | -s integer                 | " << defaults.getSubsampleRate() << endl;
  cerr << "  output root file name  | -o name                    | " << defaults.getRootFileName() << endl;
  cerr << "  input file list file   | -i filename                | " << defaults.getInputListFileName() << endl;
  cerr << "  random seed 1          | -s1 integer                | " << defaults.getSeed1() << endl;
  cerr << "  random seed 2          | -s2 integer                | " << defaults.getSeed2() << endl;
  cerr << "  create sample file     | --create-sample-file       | " << (defaults.getCreateSampleFile() == true ? "true" : "false") << endl;
  cerr << "  create joint file      | --create-joint-file        | " << (defaults.getCreateJointFile() == true ? "true" : "false") << endl;
  cerr << "  create single file     | --create-single-file       | " << (defaults.getCreateSingleFile() == true ? "true" : "false") << endl;
  cerr << "  use independence prior | --use-independence-prior   | " << (defaults.getUseIndependencePrior() == true ? "true" : "false") << endl;
  cerr << "  calculate pairs        | --calculate-pairs          | " << (defaults.getCalculatePairs() == true ? "true" : "false") << endl;
  cerr << "  use update groups      | --use-update-groups        | " << (defaults.getUseUpdateGroups() == true ? "true" : "false") << endl;
  cerr << "  use update groups      | --do-not-use-update-groups |" << endl;
  cerr << "  help                   | -h OR --help               |" << endl;
  cerr << "  version                | --version                  |" << endl;
  cerr << "  -------------------------------------------------------------------" << endl << endl;
  exit(1);
}

void showParameters(ostream& f,FileNames& fn,Defaults defaults,ModelParameters& mp,RunParameters& rp)
{
  f << "  Parameter              | Usage                    | Default Value | Value Used" << endl;
  f << "  ------------------------------------------------------------------------------" << endl;
  f << "  alpha                  | -a number                | " << left << setw(14) << defaults.getAlpha()                                             << "| " << mp.getAlpha() << endl;
  f << "  # of runs              | -k integer               | " << left << setw(14) << defaults.getNumRuns()                                           << "| " << rp.getNumRuns() << endl;
  f << "  # of MCMC updates      | -n integer               | " << left << setw(14) << defaults.getNumUpdates()                                        << "| " << rp.getNumUpdates() << endl;
  f << "  # of chains            | -c integer               | " << left << setw(14) << defaults.getNumChains()                                         << "| " << rp.getNumChains() << endl;
  f << "  MCMCMC Rate            | -r integer               | " << left << setw(14) << defaults.getMCMCMCRate()                                        << "| " << rp.getMCMCMCRate() << endl;
  f << "  alpha multiplier       | -m number                | " << left << setw(14) << defaults.getAlphaMultiplier()                                   << "| " << rp.getAlphaMultiplier() << endl;
  f << "  subsample rate         | -s integer               | " << left << setw(14) << defaults.getSubsampleRate()                                     << "| " << rp.getSubsampleRate() << endl;
  f << "  output root file name  | -o name                  | " << left << setw(14) << defaults.getRootFileName()                                      << "| " << fn.getRootFileName() << endl;
  f << "  input file list file   | -i filename              | " << left << setw(14) << defaults.getInputListFileName()                                      << "| " << fn.getInputListFileName() << endl;
  f << "  random seed 1          | -s1 integer              | " << left << setw(14) << defaults.getSeed1()                                             << "| " << rp.getSeed1() << endl;
  f << "  random seed 2          | -s2 integer              | " << left << setw(14) << defaults.getSeed2()                                             << "| " << rp.getSeed2() << endl;
  f << "  create sample file     | --create-sample-file     | " << left << setw(14) << (defaults.getCreateSampleFile() == true ? "true" : "false")     << "| " << (rp.getCreateSampleFile() ? "true" : "false") << endl;
  f << "  create joint file      | --create-joint-file      | " << left << setw(14) << (defaults.getCreateJointFile() == true ? "true" : "false")      << "| " << (rp.getCreateJointFile() ? "true" : "false") << endl;
  f << "  create single file     | --create-single-file     | " << left << setw(14) << (defaults.getCreateSingleFile() == true ? "true" : "false")     << "| " << (rp.getCreateSingleFile() ? "true" : "false") << endl;
  f << "  use independence prior | --use-independence-prior | " << left << setw(14) << (defaults.getUseIndependencePrior() == true ? "true" : "false") << "| " << (mp.getUseIndependencePrior() ? "true" : "false") << endl;
  f << "  calculate pairs        | --calculate-pairs        | " << left << setw(14) << (defaults.getCalculatePairs() == true ? "true" : "false")       << "| " << (rp.getCalculatePairs() ? "true" : "false" ) << endl;
  f << "  use update groups      | --use-update-groups      | " << left << setw(14) << (defaults.getUseUpdateGroups() == true ? "true" : "false")      << "| " << (rp.getUseUpdateGroups() ? "true" : "false")<< endl;
  f << "  ------------------------------------------------------------------------------" << endl;
}

void intro(ostream& f) {
  f << "Bayesian Untangling of Concordance Knots (applied to yeast and other organisms)" << endl;
  f << "BUCKy version " << VERSION << ", " << DATE << endl;
  f << COPYRIGHT << endl << endl;
  f << "This is free software; see the source for copying conditions.  There is NO" << endl;
  f << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl << endl;
}

int readArguments(int argc, char *argv[],FileNames& fn,ModelParameters& mp,RunParameters& rp,Defaults& defaults)
{
  int k=1;
  bool done = (argc>1 ? false : true);
  while(!done && k<argc) {
    string flag=argv[k];
    if(flag=="--version")
      exit(0);
    if(flag=="-h" || flag=="--help")
      usage(defaults);
    else if(flag=="-a") {
      double alpha;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> alpha) )
	usage(defaults);
      if(alpha <= 0.0)
	cerr << "Warning: parameter alpha must be positive.  Ignoring argument -a " << alpha << "." << endl;
      else
	mp.setAlpha(alpha);
      k++;
    }
    else if(flag=="-n") {
      unsigned int numUpdates;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> numUpdates) )
	usage(defaults);
      rp.setNumUpdates(numUpdates);
      k++;
    }
    else if(flag=="-k") {
      unsigned int numRuns;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> numRuns) )
	usage(defaults);
      if (numRuns==0){
	cerr << "Warning: parameter number-of-runs must be at least one. Setting to default: -k 2." << endl;
	numRuns=2;
      }
      rp.setNumRuns(numRuns);
      k++;
    }
    else if(flag=="-c") {
      unsigned int numChains;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> numChains) )
	usage(defaults);
      if(numChains==0) {
	cerr << "Warning: parameter number-of-chains must be at least one. Setting to default -c 1." << endl;
	numChains=1;
      }
      if(mp.getUseIndependencePrior()){
	cerr << "Warning: parameter number-of-chains can only be used with non-independence prior. Setting to default -c 1." << endl;
	numChains=1;
      }
      rp.setNumChains(numChains);
      k++;
    }
    else if(flag=="-m") {
      double alphaMultiplier;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> alphaMultiplier) )
	usage(defaults);
      if(alphaMultiplier<=0)
	cerr << "Warning: parameter alpha-multiplier must be positive. Ignoring argument -m " << alphaMultiplier << "." << endl;
      else
	rp.setAlphaMultiplier(alphaMultiplier);
      k++;
    }
    else if(flag=="-r") {
      unsigned int mcmcmcRate;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> mcmcmcRate) )
	usage(defaults);
      if(mcmcmcRate==0) {
	cerr << "Warning: parameter MCMCMC-rate must be at least one. Setting to default -r 1." << endl;
	mcmcmcRate=1;
      }
      rp.setMCMCMCRate(mcmcmcRate);
      k++;
    }
    else if(flag=="-s") {
      unsigned int subsampleRate;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> subsampleRate) )
	usage(defaults);
      if(subsampleRate==0) {
	cerr << "Warning: parameter subsample-rate must be at least one. Setting to default -s 1." << endl;
	subsampleRate=1;
      }
      rp.setSubsampleRate(subsampleRate);
      k++;
    }
    else if(flag=="-o") {
      string root = argv[++k];
      fn.setFileNames(root);
      k++;
    }
    else if(flag=="-i") {
      string filename = argv[++k];
      fn.setInputListFileName(filename);
      k++;
    }
    else if(flag=="-s1") {
      unsigned int seed1;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> seed1) )
	usage(defaults);
      if(seed1==0)
	cerr << "Warning: parameter seed1 must be at least one. Ignorning command -s1 " << seed1 << "." << endl;
      else
	rp.setSeed1(seed1);
      k++;
    }
    else if(flag=="-s2") {
      unsigned int seed2;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> seed2) )
	usage(defaults);
      if(seed2==0)
	cerr << "Warning: parameter seed2 must be at least one. Ignorning command -s2 " << seed2 << "." << endl;
      else
	rp.setSeed2(seed2);
      k++;
    }
    else if(flag=="--create-sample-file") {
      rp.setCreateSampleFile(true);
      k++;
    }
    else if(flag=="--create-joint-file") {
      rp.setCreateJointFile(true);
      k++;
    }
    else if(flag=="--create-single-file") {
      rp.setCreateSingleFile(true);
      k++;
    }
    else if(flag=="--use-independence-prior") {
      mp.setUseIndependencePrior(true);
      // to make sure alpha=Inf and is not used...
      mp.setAlpha(numeric_limits<double>::infinity());
      cerr << "Using independence prior: setting alpha=" << mp.getAlpha() << endl;
      if(rp.getNumChains()>1){
	cerr << "Warning: Independence prior can only be used with 1 chain. Setting number of chains to default -c 1." << endl;
	rp.setNumChains(1);
      }
      k++;
    }
    else if(flag=="--calculate-pairs") {
      rp.setCalculatePairs(true);
      k++;
    }
    else if(flag=="--use-update-groups") {
      rp.setUseUpdateGroups(true);
      k++;
    }
    else if(flag=="--do-not-use-update-groups") {
      rp.setUseUpdateGroups(false);
      k++;
    }
    else
      done = true;
  }

  if(fn.getInputListFileName().empty() && argc-k<1)
    usage(defaults);
  return k;
}

static const string WHITESPACE_CHARS(" \t\n\r\f\v");

// Returns the string defined by S, with all characters in CHARS
// stripped from its righthand side
string stripRight(const string& s,
                  const string& chars = WHITESPACE_CHARS) {
  string::size_type pos = s.find_last_not_of(chars);
  return (pos == string::npos ? "" : s.substr(0, pos + 1));
}

void readInputFileList(string inputListFileName, vector<string>& inputFiles) {
  ifstream f(inputListFileName.c_str());
  string filename;
  while (getline(f, filename)) {
    inputFiles.push_back(stripRight(filename));
  }
}

bool readInputFiles(vector<string>& inputFiles,vector<vector<double> >& table, vector<string>& translateTable,
		    vector<string>& topologies,int& max, ostream& fout, vector<vector<int> >& taxid)
{
  vector<bool> hasTtable(inputFiles.size());
  bool allno = true;
  bool oneno = false;
  for (size_t i=0;i<inputFiles.size();i++){
    hasTtable[i] = readFile(inputFiles[i],i,topologies,table,max,taxid[i],translateTable);
    if (hasTtable[i]==0)
      oneno=true;
    else
      allno = false;
  }
  
  // warning if genes do not all have translate tables
  if (oneno && !allno){ // some genes have a table, but not all
      cout << "\nWarning! the following files had no translate table:"<<endl;
      for (size_t i=0;i<inputFiles.size();i++)
	if (hasTtable[i]==0)
	  cout << inputFiles[i] << endl;
      cout <<endl;
  }
  if (allno){
      cout << "\nWarning! none of the loci had a translate table."<<endl;
  }
  return(oneno);
}

void normalize(vector<vector<double> >& table) {
  for(int i=0;i<table.size();i++) {
    double sum = 0;
    for(int j=0;j<table[i].size();j++)
      sum += table[i][j];
    if(sum > 0)
      for(int j=0;j<table[i].size();j++)
	table[i][j] /= sum;
  }
}

void sortTrees(vector<vector<double> >& table,vector<string>& topologies)
{
  vector<TreeWeight> tw(topologies.size());
  for(int j=0;j<topologies.size();j++) {
    tw[j].setWeight(0);
    tw[j].setIndex(j);
  }
  for(int i=0;i<table.size();i++)
    for(int j=0;j<table[i].size();j++)
      tw[j].addWeight(table[i][j]);
  sort(tw.begin(),tw.end(),cmpTreeWeights);
  vector<string> topCopy = topologies;
  for(int j=0;j<topologies.size();j++)
    topologies[j] = topCopy[tw[j].getIndex()];
  for(int i=0;i<table.size();i++) {
    vector<double> weightCopy = table[i];
    for(int j=0;j<table[i].size();j++)
      table[i][j] = weightCopy[tw[j].getIndex()];
  }
}

void Defaults::print(ostream& f) {
  f << "BUCKy default values:" << endl;
  f << "-a alpha-value = " << alpha << endl;
  f << "-k number-of-chains = " << numRuns << endl;
  f << "-n number-of-MCMC-updates = " << numUpdates << endl;
  f << "-c number-of-chains = " << numChains << endl;
  f << "-r MCMCMC-rate = " << mcmcmcRate << endl;
  f << "-m alpha-multiplier = " << alphaMultiplier << endl;
  f << "-s subsample-rate = " << subsampleRate << endl;
  f << "-f output-file-root = " << rootFileName << endl;
  f << "-s1 seed1 = " << seed1 << endl;
  f << "-s2 seed2 = " << seed2 << endl;
  f << "--create-sample-file = " << (createSampleFile ? "true" : "false") << endl;
  f << "--create-joint-file = " << (createJointFile ? "true" : "false") << endl;
  f << "--create-single-file = " << (createSingleFile ? "true" : "false") << endl;
  f << "--use-independence-prior = " << (useIndependencePrior ? "true" : "false") << endl;
  f << "--calculate-pairs = " << (calculatePairs ? "true" : "false" ) << endl;
  f << "--use-update-groups = " << (useUpdateGroups ? "true" : "false") << endl;
}

void printParameters(ostream& f,FileNames& fn,ModelParameters& mp,RunParameters& rp) {
//this function seems to be obsolete and not used.
  f << "Parameters used:" << endl;
  f << "-a alpha-value =" << mp.getAlpha() << endl;
  f << "-k number-of-chains = " << rp.getNumRuns() << endl;
  f << "-n number-of-MCMC-updates = " << rp.getNumUpdates() << endl;
  f << "-c number-of-chains = " << rp.getNumChains() << endl;
  f << "-r MCMCMC-rate = " << rp.getMCMCMCRate() << endl;
  f << "-m alpha-multiplier = " << rp.getAlphaMultiplier() << endl;
  f << "-s subsample-rate = " << rp.getSubsampleRate() << endl;
  f << "-f output-file-root = " << fn.getRootFileName() << endl;
  f << "-s1 seed1 = " << rp.getSeed1() << endl;
  f << "-s2 seed2 = " << rp.getSeed2() << endl;
  f << "--create-sample-file = " << (rp.getCreateSampleFile() ? "true" : "false") << endl;
  f << "--create-joint-file = " << (rp.getCreateJointFile() ? "true" : "false") << endl;
  f << "--create-single-file = " << (rp.getCreateSingleFile() ? "true" : "false") << endl;
  f << "--use-independence-prior = " << (mp.getUseIndependencePrior() ? "true" : "false") << endl;
  f << "--calculate-pairs = " << (rp.getCalculatePairs() ? "true" : "false" ) << endl;
  f << "--use-update-groups = " << (rp.getUseUpdateGroups() ? "true" : "false") << endl;
}

// Fix error in finding concordance tree
// Two sets are compatible if one is the subset of the other or a subset of the complement of the other or vice versa
bool isCompatible(unsigned int x,unsigned int y,unsigned int all) {
  if( x == (x & y) )
    return true;
  if( y == (x & y) )
    return true;
  unsigned int xc = x^all;
  if( xc == (xc & y) )
    return true;
  if( y == (xc & y) )
    return true;
  return false;
}
void GenomewideDistribution::updateSamplewide(vector<int> &splitGeneCount, unsigned int numUpdates) {
  vector<double> splitGenePP(splitGeneCount.size());
  for (int j=0;j<splitGeneCount.size();j++)
    splitGenePP[j] = (double)(splitGeneCount[j])/numUpdates;
  this->updateSamplewide(splitGenePP);
}

void GenomewideDistribution::updateSamplewide(vector<double> &splitGenePP){
  if (splitGenePP.size() != samplewide.size()){
      cerr << "InternalError: wrong number of genes in GenomewideDistribution::updateSamplewide;" << endl;
      exit(1);
  }
  int ngenes=samplewide.size()-1;
  vector<double> csum(samplewide.size());
  for(int j=0;j<splitGenePP.size();j++) {
    samplewide[j] = splitGenePP[j];
    if (j>0) csum[j] =  csum[j-1]+samplewide[j];
    else csum[j] = samplewide[j];
  }
  int j=0;
  while(csum[j] < 0.005) j++;
  samplewideCredibilityInterval[0]= ((double) j)/ngenes;
  while(csum[j] < 0.025) j++;
  samplewideCredibilityInterval[1]= ((double) j)/ngenes;
  while(csum[j] < 0.05) j++;
  samplewideCredibilityInterval[2]= ((double) j)/ngenes;
  while(csum[j] < 0.95) j++;
  samplewideCredibilityInterval[3]= ((double) j)/ngenes;
  while(csum[j] < 0.975) j++;
  samplewideCredibilityInterval[4]= ((double) j)/ngenes;
  while(csum[j] < 0.995 && j<ngenes) j++;
  samplewideCredibilityInterval[5]= ((double) j)/ngenes;
}

void GenomewideDistribution::updateConvolutionWeight(double alpha){
  int ngenes = samplewide.size()-1;
  int ngrid  = genomewide.size()-1;
  double priorQ = 1 - priorProbability;
  vector<vector<double> > logweight;
  logweight.resize(ngrid+1);
  for (int i=0; i<ngrid+1; i++) logweight[i].resize(ngenes+1);
  vector<double> logalphap(ngenes+1);
  for (int j=0; j<ngenes+1; j++) logalphap[j] = log(alpha*priorProbability + j);
  vector<double> logalphaq(ngenes+1);
  for (int j=0; j<ngenes+1; j++) logalphaq[j] = log(alpha*priorQ + j);
  vector<double> xinside(ngrid-1);
  vector<double> u(ngrid-1);
  for (int i=1; i<ngrid; i++){
    xinside[i-1] = ((double) (i))/ngrid;
    u[i-1] = log(xinside[i-1])-log(1-xinside[i-1]);
  }
  for (int i=1; i<ngrid; i++){
    logweight[i][1]=alpha*priorProbability*log(xinside[i-1])+(alpha*priorQ+ngenes-2)*log(1-xinside[i-1]);
    logweight[i][0]   = logweight[i][1] - u[i-1] + logalphap[0] - logalphaq[ngenes-1];
  }
  for (int j=2; j<ngenes+1; j++){
    for (int i=1; i<ngrid; i++)
      logweight[i][j] = logweight[i][j-1]+u[i-1] - logalphap[j-1]+logalphaq[ngenes-j];
  }
  // Next: calculate the normalizing constant Z=Beta(alpha*priorP+1, alpha*priorQ+ngenes-1).
  // logZ ~ log of sum(exp(logweight[,j]))/ngrid, for any j, avoiding the need for
  // a library for the Beta or the Gamma function. Checked with R: very accurate
  // approximation as soon as ngrid > 100 and j well chosen (so that weight[,j] not too skewed).
  int bestj = (int) (ngenes/2+alpha*(0.5-priorProbability));
  if (bestj>ngenes-1) bestj=ngenes-1;
  else if (bestj<1) bestj=1;
  double logZ = 0; 
  for (int i=1; i<ngrid; i++) logZ += exp(logweight[i][bestj]);
  logZ = log( logZ/ngrid );
  for (int i=1; i<ngrid; i++)
    for (int j=0; j<ngenes+1; j++)
      convolutionWeight[i][j] = exp(logweight[i][j] -logZ);
  for (int j=0; j<ngenes+1; j++){
    convolutionWeight[0][j] = 0.0;
    convolutionWeight[ngrid][j] = 0.0;
  }
  convolutionWeight[0][0]=(double) ngrid;
  for (int i=1; i<ngrid+1; i++)
    convolutionWeight[0][0] -= convolutionWeight[i][0];
  convolutionWeight[ngrid][ngenes]=(double) ngrid;
  for (int i=0; i<ngrid; i++)
    convolutionWeight[ngrid][ngenes] -= convolutionWeight[i][ngenes];
}

void GenomewideDistribution::updateGenomewide(double alpha) {
  int ngenes = samplewide.size()-1;
  int ngrid  = genomewide.size()-1;
  genomewidePosteriorMean = alpha/(alpha+ngenes)*priorProbability + ngenes/(alpha+ngenes)*samplewidePosteriorMean;

  vector<double> csum(genomewide.size());
  for(int i=0;i<ngrid+1;i++) { 
    // matrix multiplication: genomewide = convolutionWeights * samplewide (density)
    double sumprod=0;
    for (int j=0; j<ngenes+1; j++)
      sumprod += convolutionWeight[i][j]*samplewide[j]/ngrid;
    genomewide[i]=sumprod;
    if (i>0) csum[i] += csum[i-1]+genomewide[i];
    else csum[i] = genomewide[i];
  }
  int i=0;
  while(csum[i] < 0.005) i++;
  genomewideCredibilityInterval[0]= ((double) i)/ngrid;
  while(csum[i] < 0.025) i++;
  genomewideCredibilityInterval[1]= ((double) i)/ngrid;
  while(csum[i] < 0.05) i++;
  genomewideCredibilityInterval[2]= ((double) i)/ngrid;
  while(csum[i] < 0.95) i++;
  genomewideCredibilityInterval[3]= ((double) i)/ngrid;
  while(csum[i] < 0.975) i++;
  genomewideCredibilityInterval[4]= ((double) i)/ngrid;
  while(csum[i] < 0.995 && i<ngrid) i++;
  genomewideCredibilityInterval[5]= ((double) i)/ngrid;
}


// Currently one function to write all output
// Should have separate functions for each type of output
void writeOutput(ostream& fout,FileNames& fileNames,int max,int numTrees,int numTaxa,vector<string> topologies,
		 int numGenes,RunParameters& rp,ModelParameters& mp,vector<vector<int> >& newTable,
		 vector<vector<int> >& clusterCount, vector<int>& splits,  vector<vector<vector<int> > >& splitsGeneMatrix,
		 vector<vector<int> >& pairCounts,   vector<Gene*>& genes, vector<double>& alphas,
		 vector<vector<int> >& mcmcmcAccepts,vector<vector<int> >& mcmcmcProposals)
{
  // .joint
  if(rp.getCreateJointFile()) {
    cout << "Writing joint posterior table to file " << fileNames.getJointFile() << "...." << flush;
    fout << "Writing joint posterior table to file " << fileNames.getJointFile() << "...." << flush;
    ofstream jointStr(fileNames.getJointFile().c_str());
    jointStr.setf(ios::fixed, ios::floatfield);
    jointStr.setf(ios::showpoint);
    for(int i=0;i<numTrees;i++) {
      jointStr << setw(4) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right;
      for(int j=0;j<numGenes;j++)
	jointStr << " " << setw(8) << setprecision(5) << newTable[i][j] / ((double)rp.getNumRuns()*rp.getNumUpdates());
      jointStr << endl;
    }
    cout << "done." << endl;
    fout << "done." << endl;
  }

  // .cluster
  cout << "Writing cluster summary to file " << fileNames.getClusterFile() << "...." << flush;
  fout << "Writing cluster summary to file " << fileNames.getClusterFile() << "...." << flush;
  ofstream clusterStr(fileNames.getClusterFile().c_str());
  clusterStr.setf(ios::fixed, ios::floatfield);
  clusterStr.setf(ios::showpoint);

  vector<double> clusterPP(numGenes+1);
  vector<double> wsum(rp.getNumRuns());
  double wsumAvg=0.0, wsumSD=0.0;
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
    wsum[irun]=0.0;

  for (int i=0; i<numGenes+1; i++){
    clusterPP[i]=0.0;
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
      clusterPP[i] += (double)clusterCount[irun][i];
    wsumAvg += i * clusterPP[i];
    clusterPP[i] /= ((double)rp.getNumUpdates() * rp.getNumRuns());
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
      wsum[irun] += i* (double)clusterCount[irun][i];
  }
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
    wsum[irun] /= (double)rp.getNumUpdates();
  wsumAvg /= ((double)rp.getNumUpdates() * rp.getNumRuns());

  clusterStr << "mean #groups = " << setprecision(3) << wsumAvg << endl ;
  if (rp.getNumRuns()>1){
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
      wsumSD += (wsum[irun]-wsumAvg)*(wsum[irun]-wsumAvg);
    wsumSD = sqrt(wsumSD / (rp.getNumRuns()-1));
    clusterStr << "SD across runs = " << wsumSD <<endl;
  }
  clusterStr << endl ;

  clusterStr << "credible regions for # of groups" << endl;
  clusterStr << "probability region" << endl;
  clusterStr << "------------------" << endl;

  int a=numGenes,b=0;
  for(int i=0;i<numGenes+1;i++)
    if(clusterPP[i]>0) {
      if(a>i)	a = i;
      if(b<i)	b = i;
    }

  int lo=a,hi;
  double sum = clusterPP[lo];
  while(sum < .005)    sum += clusterPP[++lo];
  hi=lo;
  while(sum < .995)    sum += clusterPP[++hi];
  clusterStr << "  0.99      (" << lo << "," << hi << ")" << endl;

  lo=a;  sum = clusterPP[lo];
  while(sum < .025)    sum += clusterPP[++lo];
  hi=lo;
  while(sum < .975)    sum += clusterPP[++hi];
  clusterStr << "  0.95      (" << lo << "," << hi << ")" << endl;

  lo=a;  sum = clusterPP[lo];
  while(sum < .050)    sum += clusterPP[++lo];
  hi=lo;
  while(sum < .950)    sum += clusterPP[++hi];
  clusterStr << "  0.90      (" << lo << "," << hi << ")" << endl;
  clusterStr << "------------------" << endl << endl;

  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
    clusterStr << "Distribution of cluster number in run " << irun+1 << ":" << endl;
    clusterStr << " # of    raw    posterior" << endl;
    clusterStr << "groups  counts probability" << endl;
    clusterStr << "--------------------------" << endl;
    for(int i=a;i<b+1;i++)
      clusterStr << setw(3) << i << " " << setw(10) << clusterCount[irun][i]
		 << setw(12) << setprecision(8) << clusterCount[irun][i] / (double)rp.getNumUpdates() << endl;
    clusterStr << "--------------------------" << endl << endl;
  }
  clusterStr.close();
  cout << "done." << endl;
  fout << "done." << endl;


  // .concordance
  cout << "Writing concordance factors to " << fileNames.getConcordanceFile() << "...." << flush;
  fout << "Writing concordance factors to " << fileNames.getConcordanceFile() << "...." << flush;
  ofstream concordanceStr(fileNames.getConcordanceFile().c_str());
  concordanceStr.setf(ios::fixed, ios::floatfield);
  concordanceStr.setf(ios::showpoint);

  vector<vector<double> > splitsGeneMatrixPP(splits.size());
  for(int i=0;i<splits.size();i++) {
    splitsGeneMatrixPP[i].resize(numGenes+1);
    for(int j=0;j<numGenes+1;j++){
      splitsGeneMatrixPP[i][j]=0.0;
      for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
	splitsGeneMatrixPP[i][j] += (double)splitsGeneMatrix[irun][i][j];
      splitsGeneMatrixPP[i][j] /= ((double)rp.getNumUpdates() * rp.getNumRuns());
    }
  }

  // use TreeWeight class for sorting splits by mean number of genes
  vector<TreeWeight> sw(splits.size());
  for(int i=0;i<splits.size();i++) {
    sw[i].setWeight(0);
    for(int j=0;j<numGenes+1;j++)
      if(splitsGeneMatrixPP[i][j]>0)
	sw[i].addWeight(j* splitsGeneMatrixPP[i][j]);
    sw[i].setIndex(i);
  }
  sort(sw.begin(),sw.end(),cmpTreeWeights);

  unsigned int all=0; //fixit: What is this for??
  for(int i=0;i<numTaxa;i++) {
    all <<= 1;
    all += 1;
  }

  // ctree is the vector of splits in the primary concordance tree
  // gwDistr is the vector of GenomewideDistribution for splits in the primary concordance tree
  vector<Split> ctree;
  vector<GenomewideDistribution*> gwDistr;
  // otherClade is the vector of splits with CF>.10 but not in the concordance tree
  // otherGwDistr is the vector of GenomewideDistribution for those splits
  vector<Split> otherClade;
  vector<GenomewideDistribution*> otherGwDistr;

  for(int w=0;w<splits.size();w++) {
    int i = sw[w].getIndex();
    Split s(splits[i],numTaxa);
    s.setWeight(sw[w].getWeight());
    unsigned int y = s.getClade();
    bool keep=true;
    if (sw[w].getWeight() < .05 * numGenes) // arbitrary cutoff: keep if sample wide CF >= .05
      keep=false;
    bool add=true;
    for(int j=0;j<ctree.size();j++) {
      unsigned int z = ctree[j].getClade();
      if(!isCompatible(y,z,all)) {
	add = false;
	break;
      }
    }
    if(add || keep){
      s.updatePriorProbability();
      if (add)
	ctree.push_back(s);
      else
	otherClade.push_back(s);
      GenomewideDistribution* g;
      g = new GenomewideDistribution(numGenes, rp.getNumGenomewideGrid());
      g->updatePriorProbability(s);
      g->updateSamplewide(splitsGeneMatrixPP[i]);
      if (!mp.getUseIndependencePrior()){ 
	// under independence: genomewide CF is concentrated on priorProbability.
	g->updateConvolutionWeight(alphas[0]);
	g->updateGenomewide(alphas[0]);
      }
      double meanpostCFsd=0.0; // SD of estimated sample-wide CF, across runs
      if (rp.getNumRuns()>1){
	for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
	  double meanpostCF =0.0;
	  for(int j=0;j<numGenes+1;j++)
	    meanpostCF += j * (double)splitsGeneMatrix[irun][i][j];
	  meanpostCF /= (double)rp.getNumUpdates();
	  meanpostCFsd += (meanpostCF - s.getWeight()) * (meanpostCF - s.getWeight());
	}
	meanpostCFsd = sqrt(meanpostCFsd / (rp.getNumRuns()-1));  // number of genes
	g->setSamplewidePosteriorMeanSD(meanpostCFsd / numGenes); // proportion of genes
      }
      if (add)
	gwDistr.push_back(g);
      else
	otherGwDistr.push_back(g);
    }
  }

  // Now, use ctree to actually build the tree
  // Begin by creating a TaxonSet for each split and then partially ordering them so that subsets precede supersets

  double AvgOfSDacrossSplits = 0.0;

  vector<TaxonSet> tset;
  for(vector<Split>::iterator p=ctree.begin();p!=ctree.end();p++)
    tset.push_back( *(new TaxonSet(*p)) );

  sort(tset.begin(),tset.end());
  TaxonSet allTaxa(numTaxa);
  allTaxa.setAll();
  tset.push_back( *(new TaxonSet(allTaxa)) );

  ConcordanceTree z(tset,numGenes);
  concordanceStr << "Primary Concordance Tree Topology:" << endl;
  z.printTopology(concordanceStr);

  concordanceStr << "Primary Concordance Tree with Sample Concordance Factors:" << endl;
  z.print(concordanceStr);

  concordanceStr << "Splits in the Primary Concordance Tree: sample-wide ";
  if (!mp.getUseIndependencePrior())
    concordanceStr << "and genome-wide ";
  concordanceStr << "mean CF (95% credibility)";
  if (rp.getNumRuns()>1)
    concordanceStr << ", SD of mean sample-wide CF across runs";
  concordanceStr << endl;

  for(int w=0;w<ctree.size();w++) {
    ctree[w].print(concordanceStr); 
    gwDistr[w]->printSampleCF(concordanceStr);
    if (!mp.getUseIndependencePrior())
      gwDistr[w]->printGenomeCF(concordanceStr);
    if (rp.getNumRuns()>1){
      concordanceStr << "\t" << gwDistr[w]->getSamplewidePosteriorMeanSD();
      AvgOfSDacrossSplits += gwDistr[w]->getSamplewidePosteriorMeanSD();
    }
    concordanceStr << endl;
  }
  concordanceStr << endl;
  concordanceStr << "Splits NOT in the Primary Concordance Tree but with estimated CF>.10:"<<endl;
  for(int w=0;w<otherClade.size();w++) {
    otherClade[w].print(concordanceStr); 
    otherGwDistr[w]->printSampleCF(concordanceStr);
    if (!mp.getUseIndependencePrior())
      otherGwDistr[w]->printGenomeCF(concordanceStr);
    if (rp.getNumRuns()>1){
      concordanceStr << "\t" << otherGwDistr[w]->getSamplewidePosteriorMeanSD();
      AvgOfSDacrossSplits += otherGwDistr[w]->getSamplewidePosteriorMeanSD();
    }
    concordanceStr << endl;
  }
  concordanceStr << endl;

  AvgOfSDacrossSplits /= (gwDistr.size()+otherGwDistr.size());
  concordanceStr<<"Average SD of mean sample-wide CF: " << AvgOfSDacrossSplits<<endl<<endl;

  concordanceStr << "All Splits:" << endl;

  for(int w=0;w<splits.size();w++) {
    int i = sw[w].getIndex();
    Split s(splits[i],numTaxa);
    s.print(concordanceStr);
    concordanceStr << endl;

    int a=numGenes,b=0;
    for(int j=0;j<numGenes+1;j++)
      if(splitsGeneMatrixPP[i][j]>0) {
	if(a>j)	a = j;
	if(b<j)	b = j;
      }

    concordanceStr << "#Genes      count in run(s) 1";
    if (rp.getNumRuns()>1)
      concordanceStr << " through "<< rp.getNumRuns();
    concordanceStr << ", Overall probability, Overall cumulative probability" << endl;
    double csum=0;
    for(int j=a;j<b+1;j++) {
      concordanceStr << setw(6) << j ;
      for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
	concordanceStr << " " << setw(10) << splitsGeneMatrix[irun][i][j];
      csum += splitsGeneMatrixPP[i][j];
      concordanceStr << " " << setw(12) << setprecision(6) << splitsGeneMatrixPP[i][j]
		     << " " << setw(11) << setprecision(6) << csum << endl;
    }
    concordanceStr << endl;
    concordanceStr << "mean CF = " << setw(7) << setprecision(3) << sw[w].getWeight() / numGenes << " (proportion of loci)" << endl;
    concordanceStr << "        = " << setw(7) << setprecision(3) << sw[w].getWeight()            << " (number of loci)" << endl;

    int lo=a,hi;
    sum = splitsGeneMatrixPP[i][lo];
    while(sum < .005)  sum += splitsGeneMatrixPP[i][++lo];
    hi=lo;
    while(sum < .995)  sum += splitsGeneMatrixPP[i][++hi];
    concordanceStr << "99% CI for CF = (" << lo << "," << hi << ")" << endl;

    lo=a;
    sum = splitsGeneMatrixPP[i][lo];
    while(sum < .025)  sum += splitsGeneMatrixPP[i][++lo];
    hi=lo;
    while(sum < .975)  sum += splitsGeneMatrixPP[i][++hi];
    concordanceStr << "95% CI for CF = (" << lo << "," << hi << ")" << endl;

    lo=a;
    sum = splitsGeneMatrixPP[i][lo];
    while(sum < .050)  sum += splitsGeneMatrixPP[i][++lo];
    hi=lo;
    while(sum < .950)  sum += splitsGeneMatrixPP[i][++hi];
    concordanceStr << "90% CI for CF = (" << lo << "," << hi << ")" << endl << endl;
  }
  concordanceStr.close();
  cout << "done." << endl;    
  fout << "done." << endl;    

  cout << "Average SD of mean sample-wide CF: " << AvgOfSDacrossSplits<<endl;
  fout << "Average SD of mean sample-wide CF: " << AvgOfSDacrossSplits<<endl;

  // .pairs
  if(rp.getCalculatePairs()) {
    cout << "Writing tree pair data to " << fileNames.getPairTreeFile() << "...." << flush;
    fout << "Writing tree pair data to " << fileNames.getPairTreeFile() << "...." << flush;
    ofstream pairsStr(fileNames.getPairTreeFile().c_str());
    pairsStr.setf(ios::fixed, ios::floatfield);
    pairsStr.setf(ios::showpoint);
    double total = pairCounts[0][0];
    for(int i=0;i<numGenes;i++) {
      pairsStr << setw(4) << i << " ";
      for(int j=0;j<numGenes;j++)
	pairsStr << setw(7) << setprecision(4) << (double) (i < j ? pairCounts[i][j] : pairCounts[j][i]) / total;
      pairsStr << endl;
    }
    pairsStr.close();
    cout << "done." << endl;
    fout << "done." << endl;
  }

  // .topologies --- eliminated
//  cout << "Writing topology single and joint posterior distribution to " << fileNames.getTreePosteriorFile() << "...." << flush;
//  fout << "Writing topology single and joint posterior distribution to " << fileNames.getTreePosteriorFile() << "...." << flush;
//  ofstream treePosteriorStr(fileNames.getTreePosteriorFile().c_str());
//  treePosteriorStr.setf(ios::fixed, ios::floatfield);
//  treePosteriorStr.setf(ios::showpoint);
//  treePosteriorStr << "index "
//		   << setw(max) << left << "topology" << right
//		   << setw(10) << "single"
//		   << setw(10) << "joint" << endl;
//  for(int i=0;i<numTrees;i++) {
//    double single=0.0;
//    int joint=0;
//    for(int j=0;j<numGenes;j++) {
//      single += genes[j]->getProb(i);
//      joint += newTable[i][j];
//    }
//    treePosteriorStr << setw(5) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right
//		     << setw(10) << setprecision(6) << (double) single / (double) numGenes
//		     << setw(10) << setprecision(6) << (double) joint / (double) rp.getNumUpdates() / (double) numGenes
//		     << endl;
//  }
//  treePosteriorStr.close();
//  cout << "done." << endl;
//  fout << "done." << endl;

// .gene
  cout << "Writing single and joint gene posteriors to " << fileNames.getGenePosteriorFile() << "...." << flush;
  fout << "Writing single and joint gene posteriors to " << fileNames.getGenePosteriorFile() << "...." << flush;
  ofstream genePostStr(fileNames.getGenePosteriorFile().c_str());
  for(int i=0;i<numGenes;i++)
    genes[i]->print(genePostStr,newTable,rp.getNumUpdates()*rp.getNumRuns(),topologies,max);
  genePostStr.close();
  cout << "done." << endl;
  fout << "done." << endl;

  if(rp.getNumChains()>1) {
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    fout.setf(ios::fixed, ios::floatfield);
    fout.setf(ios::showpoint);
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
      cout << "\nMCMCMC acceptance statistics in run " << irun+1 <<":" << endl;
      cout << "         alpha1 <-->         alpha2 accepted proposed proportion" << endl;
      fout << "\nMCMCMC acceptance statisticsin run " << irun+1 <<":" << endl;
      fout << "alpha1          <--> alpha2         accepted proposed proportion" << endl;
      for(int i=0;i<rp.getNumChains()-1;i++) {
	cout << setw(15) << alphas[i] << " <-->" << setw(15) << alphas[i+1];
	cout << setw(9) << mcmcmcAccepts[irun][i] << setw(9) << mcmcmcProposals[irun][i];
	fout << setw(15) << alphas[i] << " <--> " << setw(15) << alphas[i+1];
	fout << setw(9) << mcmcmcAccepts[irun][i] << setw(9) << mcmcmcProposals[irun][i];
	if(mcmcmcProposals[irun][i] > 0) {
	  cout << setw(11) << setprecision(6) << (double) (mcmcmcAccepts[irun][i]) / (double) (mcmcmcProposals[irun][i]);
	  fout << setw(11) << setprecision(6) << (double) (mcmcmcAccepts[irun][i]) / (double) (mcmcmcProposals[irun][i]);
	}
	cout << endl;
	fout << endl;
      }
    }
  }
}


int main(int argc, char *argv[])
{

  time_t beginTime;
  time(&beginTime);
  intro(cout);
  // Set Default Parameters
  int debug=0;
  FileNames fileNames("run1");
  ModelParameters mp;
  RunParameters rp;
  Defaults defaults(fileNames,mp,rp);

  // Read arguments from the command line
  // return value k is the position of the first input file
  //   i.e. argv[k] is the first input file
  int k=readArguments(argc,argv,fileNames,mp,rp,defaults);

  // Populate input files vector from remaining command line arguments
  vector<string> inputFiles(argv + k, argv + argc);

  // If input list filename is given as argument, read input filenames
  if (!(fileNames.getInputListFileName().empty())) {
    readInputFileList(fileNames.getInputListFileName(), inputFiles);
  }

  int numGenes = inputFiles.size();

  // table has one row for each input file and one column for each observed topology
  // weights for each are doubles rather than integers to allow for non-sample-based input
  vector<vector<double> > table(numGenes);
  vector<string> topologies;
  int max=0;
  int numTaxa;
  vector<string> translateTable(0);
  vector<vector<int> > taxid(numGenes);

  cout << "Reading in summary files...." << flush;
  ofstream fout(fileNames.getOutFile().c_str());
  fout << "Reading in summary files...." << flush;
  bool missingTtable = readInputFiles(inputFiles,table,translateTable,topologies,max,fout,taxid);
  numTaxa = translateTable.size();
  cout << "done." << endl;

  // Open outfile to save all window output
  cout << "Screen output written to file " << fileNames.getOutFile() << endl;
  intro(fout);

  fout << "Program initiated at " << ctime(&beginTime) << endl << endl;

  fout << "Program invocation:";
  for(int i=0;i<argc;i++)
    fout << " " << argv[i];
  fout << endl << endl;
  //  defaults.print(fout);
  //  fout << endl;
  //  printParameters(fout,fileNames,mp,rp);
  //  fout << endl;
  showParameters(fout,fileNames,defaults,mp,rp);
  fout << endl;

  if (missingTtable) {
    cout << "Warning: at least one locus has no 'translate' block." << endl;
    fout << "Warning: at least one locus has no 'translate' block." << endl;
  }

  // Check for missing taxa.
  for(int i=0; i<inputFiles.size(); i++) {
    bool problem = false;
    if (taxid[i].size() != translateTable.size())
      problem=true;
    if (problem) {
      cerr <<"\nError: Expecting " << translateTable.size() << " taxa."<<endl;
      fout <<"\nError: Expecting " << translateTable.size() << " taxa."<<endl;
      cerr <<  "       File " << inputFiles[i] << " contains trees with the following " << taxid[i].size() << " taxa:" <<endl;
      fout <<  "       File " << inputFiles[i] << " contains trees with the following " << taxid[i].size() << " taxa:" <<endl;
      for (int j=0; j<taxid[i].size(); j++){
	cerr << setw(4) << j+1 << " " << setw(4) << taxid[i][j];
	fout << setw(4) << j+1 << " " << setw(4) << taxid[i][j];
	if (taxid[i][j] <= translateTable.size()){
	  cerr << " " << translateTable[taxid[i][j]-1] << endl;
	  fout << " " << translateTable[taxid[i][j]-1] << endl;}
	else {
	  cerr << " " << "(outside range of translate table)" <<endl;
	  fout << " " << "(outside range of translate table)" <<endl; }
      }
      exit(1);
    }
  }

  // Find taxa sequenced in all genes
  TaxonList allTaxList(translateTable);
  allTaxList.setNumberGenes(taxid);
  vector<int> taxlist = allTaxList.getNonMissing();
  allTaxList.setInclude(taxlist);
  allTaxList.print(cout);
  allTaxList.print(fout);
  // fixit:
  // prune the missing taxa from all topologies, which is a vector<string>
  // find the unique topologies
  // collapse the table to only unique topologies, with updated counts
  // ??update the translate table and the taxid vectors??

  if(numTaxa > 32) { // are you sure it works with exactly 32 taxa?
    cerr << "Error: BUCKy version " << VERSION << " only works with 32 or fewer taxa." << endl;
    fout << "Error: BUCKy version " << VERSION << " only works with 32 or fewer taxa." << endl << flush;
    exit(1);
  }

  int numTrees = topologies.size();
  cout << "Read " << numGenes << " genes with a total of " << numTrees << " different sampled tree topologies" << endl;
  fout << "Read " << numGenes << " genes with a total of " << numTrees << " different sampled tree topologies" << endl;

  cout << "Writing input file names to file " << fileNames.getInputFile() << "...." << flush;
  fout << "Writing input file names to file " << fileNames.getInputFile() << "...." << flush;
  ofstream inputStr(fileNames.getInputFile().c_str());
  inputStr << "Gene Filename" << endl;
  inputStr << "============================================" << endl;
  for(int i=0;i<inputFiles.size();i++)
    inputStr << setw(4) << i << " " << inputFiles[i] << endl;
  inputStr << "============================================" << endl;
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Sorting trees by average posterior probability...." << flush;
  fout << "Sorting trees by average posterior probability...." << flush;
  normalize(table);
  sortTrees(table,topologies);
  cout << "done." << endl;
  fout << "done." << endl;

  // Print table of individual gene posteriors

  // .single
  if(rp.getCreateSingleFile()) {
    cout << "Writing single gene posterior distribution table to file " << fileNames.getSingleFile() << "...." << flush;
    fout << "Writing single gene posterior distribution table to file " << fileNames.getSingleFile() << "...." << flush;
    ofstream singleStr(fileNames.getSingleFile().c_str());
    singleStr.setf(ios::fixed, ios::floatfield);
    singleStr.setf(ios::showpoint);
    
    for(int i=0;i<numTrees;i++) {
      singleStr << setw(4) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right;
      for(int j=0;j<numGenes;j++)
	singleStr << " " << setw(8) << table[j][i];
      singleStr << endl;
    }
    singleStr.close();
    cout << "done." << endl;
    fout << "done." << endl;
  }

  // Initialize random number generator.

  cout << "Initializing random number generator...." << flush;
  fout << "Initializing random number generator...." << flush;
  Rand rand(rp.getSeed1(),rp.getSeed2());
  cout << "done." << endl;
  fout << "done." << endl;

  // Create genes

  cout << "Initializing gene information...." << flush;
  fout << "Initializing gene information...." << flush;
  vector<Gene*> genes(numGenes);
  vector<double> counts(numTrees);
  for(int i=0;i<numGenes;i++) {
    for(int j=0;j<numTrees;j++)
      counts[j] = table[i][j];
    genes[i] = new Gene(i,counts,taxid[i]);
    //genes[i]->printTaxa(cerr,translateTable);
  }
  cout << "done." << endl << flush;
  fout << "done." << endl << flush;

  // Finished with table and taxid, so delete them.

  for(int i=0;i<table.size();i++){
    table[i].clear();
    taxid[i].clear();
  }
  table.clear();
  taxid.clear();

  // old .gene --- eliminated
//  cout << "Writing summary of gene information to file " << fileNames.getGeneFile() << "...." << flush;
//  fout << "Writing summary of gene information to file " << fileNames.getGeneFile() << "...." << flush;
//  ofstream geneStr(fileNames.getGeneFile().c_str());
//  for(int i=0;i<numGenes;i++)
//    genes[i]->print(geneStr,inputFiles[i]);
//  geneStr.close();
//  cout << "done." << endl;
//  fout << "done." << endl;

// .top --- eliminated
//  cout << "Writing topology summary to file " << fileNames.getTopologyFile() << "...." << flush;
//  fout << "Writing topology summary to file " << fileNames.getTopologyFile() << "...." << flush;
  vector<vector<int> > topologySplitsMatrix(topologies.size());
  if(numTaxa>3)
    for(int i=0;i<topologies.size();i++)
      topologySplitsMatrix[i].resize(numTaxa-3);
  vector<int> splits; // sorted list of realized splits
//  ofstream topologyStr(fileNames.getTopologyFile().c_str());
  for(int i=0;i<topologies.size();i++) {
//    topologyStr << setw(4) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right;
    Tree t(numTaxa,topologies[i]);
    SplitSet s;
    t.getSplits(s);
//    s.printShort(topologyStr);
//    topologyStr << endl;

    for(int j=0;j<numTaxa-3;j++) {
      int x = s.getSplit(j);
      topologySplitsMatrix[i][j] = x;
      int n = find(splits.begin(),splits.end(),x) - splits.begin();
      if(n==splits.size())
	splits.push_back(x);
    }
  }
//  topologyStr.close();
//  cout << "done." << endl;
//  fout << "done." << endl;

  cout << "Initializing splits counts (found " << splits.size() << " distinct splits)...." << flush;
  fout << "Initializing splits counts (found " << splits.size() << " distinct splits)...." << flush;
  sort(splits.begin(),splits.end());

  // topologySplitsIndexMatrix[i][j] = index into splits vector of the jth split of topology i
  vector<vector<int> > topologySplitsIndexMatrix(topologies.size());
  if(numTaxa>3) {
    for(int i=0;i<topologies.size();i++) {
      topologySplitsIndexMatrix[i].resize(numTaxa-3);
      for(int j=0;j<numTaxa-3;j++)
	topologySplitsIndexMatrix[i][j] = find(splits.begin(),splits.end(),topologySplitsMatrix[i][j]) - splits.begin();
    }
  }

  // splitsGeneMatrix[irun][i][j] = the number of times that split i is in exactly j genes, in run 'irun'.
  vector<vector<vector<int> > > splitsGeneMatrix(rp.getNumRuns());
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
    splitsGeneMatrix[irun].resize(splits.size());
    for(int i=0;i<splits.size();i++) {
      splitsGeneMatrix[irun][i].resize(numGenes+1);
      for(int j=0;j<numGenes+1;j++)
	splitsGeneMatrix[irun][i][j] = 0;
    }
  }
  cout << "done." << endl;
  fout << "done." << endl;

  // .split --- eliminated
//  cout << "Writing splits key to file " << fileNames.getSplitsFile() << "...." << flush;
//  fout << "Writing splits key to file " << fileNames.getSplitsFile() << "...." << flush;
//  ofstream splitsStr(fileNames.getSplitsFile().c_str());
//  for(int i=0;i<splits.size();i++) {
//    splitsStr << setw(3) << splits[i] << " ";
//    Split s(splits[i],numTaxa);
//    s.print(splitsStr);
//    splitsStr << endl;
//  }
//  splitsStr.close();
//  cout << "done." << endl;
//  fout << "done." << endl;

  if(rp.getNumUpdates()==0)
    return 0;

  cout << "Setting initial MCMC state...." << flush;
  fout << "Setting initial MCMC state...." << flush;

  vector<double> alphas(rp.getNumChains());
  alphas[0] = mp.getAlpha();
  for(int i=1;i<rp.getNumChains();i++)
    alphas[i] = rp.getAlphaMultiplier()*alphas[i-1];

  vector<vector<State*> > states(rp.getNumRuns());
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
    states[irun].resize(rp.getNumChains());

  vector<vector<int> > index(rp.getNumRuns());
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
    index[irun].resize(rp.getNumChains());

  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
    for(int i=0;i<rp.getNumChains();i++) {
      states[irun][i] = new State(alphas[i],numTaxa,numTrees,genes,mp.getUseIndependencePrior(),rand);
      index[irun][i] = i;
  }
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Initializing MCMCMC acceptance counters and pairwise counters...." << flush;
  fout << "Initializing MCMCMC acceptance counters and pairwise counters...." << flush;
  vector<vector<int> > mcmcmcAccepts(rp.getNumRuns());
  vector<vector<int> > mcmcmcProposals(rp.getNumRuns());
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
    mcmcmcAccepts[irun].resize(rp.getNumChains());
    mcmcmcProposals[irun].resize(rp.getNumChains());
    for (int i=0;i<rp.getNumChains();i++)
      mcmcmcAccepts[irun][i] = mcmcmcProposals[irun][i] = 0;
  }

  vector<vector<int> > pairCounts(numGenes); // not adapted to several runs. Runs will be pooled.
  for(int i=0;i<numGenes;i++) {
    pairCounts[i].resize(numGenes);
    for(int j=0;j<numGenes;j++)
      pairCounts[i][j] = 0;
  }
  cout << "done." << endl;
  fout << "done." << endl;

  int numBurn = rp.getNumUpdates()/10;
  cout << "Beginning burn-in with " << numBurn << " updates (10% extra of desired updates)..." << endl;
  cout << "0   10   20   30   40   50   60   70   80   90   100" << endl;
  cout << "+----+----+----+----+----+----+----+----+----+----+" << endl << flush;
  fout << "Beginning burn-in with " << numBurn << " updates (10% extra of desired updates)..." << endl;
  fout << "0   10   20   30   40   50   60   70   80   90   100" << endl;
  fout << "+----+----+----+----+----+----+----+----+----+----+" << endl << flush;
  int part = numBurn / 50;
  for(int cycle=0;cycle<numBurn;cycle++) {
    if( (part > 0)  && (cycle % part == 0) ) {
      cout << "*" << flush;
      fout << "*" << flush;
    }
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
      for(int i=0;i<rp.getNumChains();i++) {
	states[irun][i]->update(rand);
	if(rp.getUseUpdateGroups()) {
	  int gene = (int)(rand.runif()*genes.size());
	  states[irun][i]->updateOneGroup(gene,rand);
	}
      }

      if(cycle % rp.getMCMCMCRate() == 0 && rp.getNumChains()>1)
	mcmcmc(states[irun],index[irun],alphas,rand,mcmcmcAccepts[irun],mcmcmcProposals[irun]);
    }
  }
  cout << "*" << endl;
  cout << " ....done." << endl << flush;
  fout << "*" << endl;
  fout << " ....done." << endl << flush;

  cout << "Initializing summary tables..." << flush;
  fout << "Initializing summary tables..." << flush;

  vector<vector<int> > clusterCount(rp.getNumRuns());
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
    clusterCount[irun].resize(numGenes+1);
    for(int i=0;i<clusterCount[irun].size();i++)
      clusterCount[irun][i] = 0;
  }

  vector<vector<int> > newTable(numTrees); // not adapted to several runs.
  for(int i=0;i<numTrees;i++) {
    newTable[i].resize(numGenes);
    for(int j=0;j<numGenes;j++)
      newTable[i][j] = 0;
  }

  cout << "done." << endl << flush;
  fout << "done." << endl << flush;

  if(rp.getCreateSampleFile()) {
    cout << "Sampled topologies will be in file(s) " << fileNames.getSampleFile(1);
    fout << "Sampled topologies will be in file(s) " << fileNames.getSampleFile(1);
    if(rp.getNumRuns()==1){
      cout << "." << endl;
      fout << "." << endl;
    } else {
      cout << " to " << fileNames.getSampleFile(rp.getNumRuns()) << "." << endl;
      fout << " to " << fileNames.getSampleFile(rp.getNumRuns()) << "." << endl;
    }
  }

  cout << "Beginning " << rp.getNumUpdates() << " MCMC updates..." << endl;
  cout << "0   10   20   30   40   50   60   70   80   90   100" << endl;
  cout << "+----+----+----+----+----+----+----+----+----+----+" << endl << flush;
  fout << "Beginning " << rp.getNumUpdates() << " MCMC updates..." << endl;
  fout << "0   10   20   30   40   50   60   70   80   90   100" << endl;
  fout << "+----+----+----+----+----+----+----+----+----+----+" << endl << flush;
  part = rp.getNumUpdates() / 50;

  vector<ofstream*> sampleFileStr(rp.getNumRuns());
  if(rp.getCreateSampleFile()) {
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
      sampleFileStr[irun] = new ofstream(fileNames.getSampleFile(irun+1).c_str());
      if (sampleFileStr[irun]->fail()){
	cerr << "Error: could not open file " << fileNames.getSampleFile(irun+1) << endl;
	exit(1);
      }
      sampleFileStr[irun]->setf(ios::fixed, ios::floatfield);
      sampleFileStr[irun]->setf(ios::showpoint);
    }
  }

  vector<vector<int> > accept(rp.getNumRuns());
  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
    accept[irun].resize(rp.getNumChains());
    for(int i=0;i<rp.getNumChains();i++)
      mcmcmcAccepts[irun][i] = mcmcmcProposals[irun][i] = accept[irun][i] = 0;
  }

  for(int cycle=0;cycle<rp.getNumUpdates();cycle++) {
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
      for(int i=0;i<rp.getNumChains();i++) {
	accept[irun][index[irun][i]] += states[irun][i]->update(rand);
	if(rp.getUseUpdateGroups()) {
	  int gene = (int)(rand.runif()*genes.size());
	  accept[irun][index[irun][i]] += states[irun][i]->updateOneGroup(gene,rand);
	}
      }
      if(cycle % rp.getMCMCMCRate() == 0 && rp.getNumChains()>1)
	mcmcmc(states[irun],index[irun],alphas,rand,mcmcmcAccepts[irun],mcmcmcProposals[irun]);
      // update counts
      int i0 = index[irun][0];
      states[irun][i0]->updateTable(newTable);
      states[irun][i0]->updateSplits(splitsGeneMatrix[irun],topologySplitsIndexMatrix);
      clusterCount[irun][states[irun][i0]->getNumGroups()]++;
      if( rp.getCalculatePairs() && cycle % rp.getSubsampleRate() == 0)
	states[irun][i0]->updatePairCounts(pairCounts);
      if( rp.getCreateSampleFile() && cycle % rp.getSubsampleRate() == 0) {
	*sampleFileStr[irun] << setw(8) << accept[irun][0];
	accept[irun][0] = 0;
	states[irun][i0]->sample(*sampleFileStr[irun]);
      }
    }
    if( cycle % part == 0) {
      cout << "*" << flush;
      fout << "*" << flush;
    }
  }

  cout << "*" << endl;
  cout << " ....done." << endl << flush;
  fout << "*" << endl;
  fout << " ....done." << endl << flush;

  if(rp.getCreateSampleFile())
    for (unsigned int irun=0; irun<rp.getNumRuns(); irun++){
	 sampleFileStr[irun]->close();
	 delete sampleFileStr[irun];
    }

  writeOutput(fout,fileNames,max,numTrees,numTaxa,topologies,numGenes,rp,mp,
	      newTable,clusterCount,splits,splitsGeneMatrix,
	      pairCounts,genes,alphas,mcmcmcAccepts,mcmcmcProposals);

  for(int i=0;i<numGenes;i++)
    delete genes[i];

  for (unsigned int irun=0; irun<rp.getNumRuns(); irun++)
    for(int i=0;i<rp.getNumChains();i++)
      delete states[irun][i];
  

  time_t endTime;
  time(&endTime);

  fout << "Program ended at " << ctime(&endTime) << endl << endl;
  int diff=endTime-beginTime,days=diff/(24*60*60),hours=diff%(24*60*60)/(60*60);
  int minutes=diff%(60*60)/60,seconds=diff%60;
  fout << "Elapsed time: ";
  if(days>0)
    fout << days << (days==1 ? " day, " : " days, ");
  if(days>0 || hours>0)
    fout << hours << (hours==1 ? " hour, " : " hours, ");
  if(days>0 || hours>0 || minutes>0)
    fout << minutes << (minutes==1 ? " minute, " : " minutes, ");
  fout << seconds << (seconds==1 ? " second." : " seconds.") << endl << flush;

  fout.close();

  return 0;
}
