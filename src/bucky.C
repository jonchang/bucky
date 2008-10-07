// BUCKy 1.1 - Bayesian Untangling of Concordance Knots (applied to yeast and other organisms)
// Copyright (C) 2006 by Bret Larget

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License, version 2, as
// published by the Free Software Foundation.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License (the file gpl.txt included with this
// distribution or http://www.gnu.org/licenses/gpl.txt) for more
// details.

// First alpha version March 7, 2005
// Last modified: October 30, 2006

// File:     bucky.C

// Input:    [Options] <summary files>
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
//            [--use-independence-prior]
//            [--calculate-pairs]
//            [--create-sample-file]
//            [-h] prints the brief usage message
//            [--help] prints the brief usage message
//
// Output:   
//           
//  fileRoot.cluster       --- summary of cluster size (number of groups)
//  fileRoot.concordance   --- split-by-split summary of concordance factor
//  fileRoot.gene          --- gene by gene summary of gene input information
//  fileRoot.genepost      --- gene by gene summary of joint posterior information
//  fileRoot.joint         --- table of raw second stage MCMC counts for joint posterior
//  fileRoot.input         --- names of input files
//  fileRoot.out           --- output file
//  fileRoot.pairs         --- output if pairwise analysis is desired (very slow!)
//  fileRoot.sample        --- output of raw sampled GTMs (very big!)
//  fileRoot.single        --- summary of table of single gene posterior probabilities
//  fileRoot.splits        --- numerical code for splits (based on binary representation of half with taxon 1)
//  fileRoot.top           --- list of sampled trees and corresponding splits
//  fileRoot.topologies    --- topology by topology summary of average single-gene and joint posterior probabilities

// Changes:
// --- added probabilities to concordance file output
// --- fixed bug with split index that assumed no more than 8 taxa
// --- implemented new and better updateGroups; use it as the default
// --- added runtime argument ---no-use-update-groups

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

#include "bucky.h"

using namespace std;
string VERSION = "1.1";
string DATE = "30 October 2006";
string COPYRIGHT = "Copyright (C) 2006 by Bret Larget";

//
// readFile
//
// input file is expected to have each line cotain two fields separated by white space
//   the first field is a topology string (parenthetic representation)
//   the second field is a positive weight (usually a count from an MCMC sample)
// add read topologies to row i of the table
// when new topologies are discovered, add a new column to the table

int readFile(string filename, int i, vector<string> &topologies, vector<vector<double> > &table, int &max)
{
  ifstream f(filename.c_str());
  if(f.fail()) {
    cerr <<"Error: Cannot open file " << filename << "." << endl;
    exit(1);
  }

  // set zeros for ith row of table
  table[i].resize(topologies.size(),0);
  int numTaxa=-1;
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
      if(numTaxa == -1)
	numTaxa = countTaxa;
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
  return numTaxa;
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
    newTop = genes[i]->pickTree(rand);
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

void State::updateOneGroupNew(int gene,Rand& rand) {
  // pick a random topology for the gene
  // check that it:
  //  (1) is not a topology for another group
  //  (2) have positive single gene prior for all genes in the group

  int oldTop = tops[gene];
  int newTop = genes[gene]->pickTree(rand);
  if(newTop==oldTop)
    return;

  // check if newTop is already in use
  if(counts[newTop]>0)
    return;

  // possible new topology for group
  vector<int> set(counts[oldTop]); // indices of genes with oldTop
  int j=0;
  double oldLogProb = 0.0;
  double newLogProb = 0.0;
  for(int i=0;i<genes.size();i++) {
    if(tops[i] == oldTop) {
      oldLogProb += genes[i]->getProb(oldTop);
      set[j++] = i;
      double p = genes[i]->getProb(newTop);
	  if(p>0)
	    newLogProb += log(p);
	  else // acceptance probability is 0
	    return;
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
  }
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
  // only fill in upper triangle.... (i < j)
  for(int i=0;i<genes.size()-1;i++)
    for(int j=0;j<genes.size();j++)
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

// assumes numChains > 1
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

void usage()
{
  cerr << "Usage: bucky \\" << endl;
  cerr << "[-a alpha-value] \\" << endl;
  cerr << "[-n number-of-MCMC-updates] \\" << endl;
  cerr << "[-c number-of-chains] \\" << endl;
  cerr << "[-r MCMCMC-rate] \\" << endl;
  cerr << "[-m alpha-multiplier] \\" << endl;
  cerr << "[-s subsample-rate] \\" << endl;
  cerr << "[-o output-file-root] \\" << endl;
  cerr << "[-s1 seed1] \\" << endl;
  cerr << "[-s2 seed2] \\" << endl;
  cerr << "[--create-sample-file] \\" << endl;
  cerr << "[--use-independence-prior] \\" << endl;
  cerr << "[--calculate-pairs] \\" << endl;
  cerr << "[--use-update-groups] \\" << endl;
  cerr << "[--no-use-update-groups] \\" << endl;
  cerr << "[-h] \\" << endl;
  cerr << "[--help] \\" << endl;
  cerr << "<summary files>" << endl;
  exit(1);
}

void intro(ostream& f) {
  f << "Bayesian Untangling of Concordance Knots (applied to yeast and other organisms)" << endl;
  f << "BUCKy version " << VERSION << ", " << DATE << endl;
  f << COPYRIGHT << endl << endl;
}

int readArguments(int argc, char *argv[],FileNames& fn,ModelParameters& mp,RunParameters& rp)
{
  int k=1;
  bool done = (argc>1 ? false : true);
  while(!done && k<argc) {
    string flag=argv[k];
    if(flag=="-h" || flag=="--help")
      usage();
    else if(flag=="-a") {
      double alpha;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> alpha) )
	usage();
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
	usage();
      rp.setNumUpdates(numUpdates);
      k++;
    }
    else if(flag=="-c") {
      unsigned int numChains;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> numChains) )
	usage();
      if(numChains==0) {
	cerr << "Warning: parameter number-of-chains must be at least one. Setting to default -c 1." << endl;
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
	usage();
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
	usage();
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
	usage();
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
    else if(flag=="-s1") {
      unsigned int seed1;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> seed1) )
	usage();
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
	usage();
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
    else if(flag=="--use-independence-prior") {
      mp.setUseIndependencePrior(true);
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
    else if(flag=="--no-use-update-groups") {
      rp.setUseUpdateGroups(false);
      k++;
    }
    else
      done = true;
  }

  if(argc-k<1)
    usage();
  return k;
}

void readInputFiles(int argc,char *argv[],int k,vector<vector<double> >& table,
		    vector<string>& topologies,vector<string>& inputFiles,int& max,int& numTaxa,ostream& fout)
{
  for(int i=k;i<argc;i++) {
    inputFiles.push_back((string)(argv[i])); 
    int numTaxaInFile = readFile((string)(argv[i]),i-k,topologies,table,max);
    if(i==k)
      numTaxa = numTaxaInFile;
    else
      if(numTaxaInFile != numTaxa) {
	cerr << endl << "Error: file " << argv[i] << " contains trees with " << numTaxaInFile << " taxa.";
	cerr << "Expected " << numTaxa << ", the same as in file " << argv[k] << "." << endl;
	fout << endl << "Error: file " << argv[i] << " contains trees with " << numTaxaInFile << " taxa.";
	fout << "Expected " << numTaxa << ", the same as in file " << argv[k] << "." << endl << flush;
	exit(1);
      }
  }
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
  f << "-n number-of-MCMC-updates = " << numUpdates << endl;
  f << "-c number-of-chains = " << numChains << endl;
  f << "-r MCMCMC-rate = " << mcmcmcRate << endl;
  f << "-m alpha-multiplier = " << alphaMultiplier << endl;
  f << "-s subsample-rate = " << subsampleRate << endl;
  f << "-f output-file-root = " << rootFileName << endl;
  f << "-s1 seed1 = " << seed1 << endl;
  f << "-s2 seed2 = " << seed2 << endl;
  f << "--create-sample-file = " << (createSampleFile ? "true" : "false") << endl;
  f << "--use-independence-prior = " << (useIndependencePrior ? "true" : "false") << endl;
  f << "--calculate-pairs = " << (calculatePairs ? "true" : "false" ) << endl;
  f << "--use-update-groups = " << (useUpdateGroups ? "true" : "false") << endl;
}

void printParameters(ostream& f,FileNames& fn,ModelParameters& mp,RunParameters& rp) {
  f << "Parameters used:" << endl;
  f << "-a alpha-value =" << mp.getAlpha() << endl;
  f << "-n number-of-MCMC-updates = " << rp.getNumUpdates() << endl;
  f << "-c number-of-chains = " << rp.getNumChains() << endl;
  f << "-r MCMCMC-rate = " << rp.getMCMCMCRate() << endl;
  f << "-m alpha-multiplier = " << rp.getAlphaMultiplier() << endl;
  f << "-s subsample-rate = " << rp.getSubsampleRate() << endl;
  f << "-f output-file-root = " << fn.getRootFileName() << endl;
  f << "-s1 seed1 = " << rp.getSeed1() << endl;
  f << "-s2 seed2 = " << rp.getSeed2() << endl;
  f << "--create-sample-file = " << (rp.getCreateSampleFile() ? "true" : "false") << endl;
  f << "--use-independence-prior = " << (mp.getUseIndependencePrior() ? "true" : "false") << endl;
  f << "--calculate-pairs = " << (rp.getCalculatePairs() ? "true" : "false" ) << endl;
  f << "--use-update-groups = " << (rp.getUseUpdateGroups() ? "true" : "false") << endl;
}

void writeOutput(ostream& fout,FileNames& fileNames,int max,int numTrees,int numTaxa,vector<string> topologies,
		 int numGenes,RunParameters& rp,vector<vector<int> >& newTable,vector<int>& clusterCount,
		 vector<int>& splits,vector<vector<int> >& splitsGeneMatrix,
		 vector<vector<int> >& pairCounts,vector<Gene*>& genes,vector<double>& alphas,
		 vector<int>& mcmcmcAccepts,vector<int>& mcmcmcProposals)
{
  cout << "Writing joint posterior table to file " << fileNames.getJointFile() << "....";
  fout << "Writing joint posterior table to file " << fileNames.getJointFile() << "....";
  ofstream jointStr(fileNames.getJointFile().c_str());
  jointStr.setf(ios::fixed, ios::floatfield);
  jointStr.setf(ios::showpoint);
  for(int i=0;i<numTrees;i++) {
    jointStr << setw(4) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right;
    for(int j=0;j<numGenes;j++)
      jointStr << " " << setw(8) << newTable[i][j];
    jointStr << endl;
  }
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Writing cluster summary to file " << fileNames.getClusterFile() << "....";
  cout << "done." << endl;
  fout << "Writing cluster summary to file " << fileNames.getClusterFile() << "....";
  fout << "done." << endl;
  ofstream clusterStr(fileNames.getClusterFile().c_str());
  clusterStr.setf(ios::fixed, ios::floatfield);
  clusterStr.setf(ios::showpoint);
  int a=numGenes,b=0;
  for(int i=0;i<numGenes+1;i++)
    if(clusterCount[i]>0) {
      if(a>i)
	a = i;
      if(b<i)
	b = i;
    }
  int q005 = (int)(rp.getNumUpdates()*0.005),q995 = rp.getNumUpdates() - q005+1;
  int q025 = (int)(rp.getNumUpdates()*0.025),q975 = rp.getNumUpdates() - q025+1;
  int q050 = (int)(rp.getNumUpdates()*0.050),q950 = rp.getNumUpdates() - q050+1;
  clusterStr << "  k      count" << endl;
  for(int i=a;i<b+1;i++)
    clusterStr << setw(3) << i << " " << setw(10) << clusterCount[i] << endl;
  int sum;
  double wsum=0;
  for(int j=a;j<b+1;j++)
    wsum += j*clusterCount[j];
  clusterStr << "mean #groups = " << setprecision(3) << wsum / rp.getNumUpdates() << endl;
  int lo,hi;
  lo=a;
  sum = clusterCount[a];
  while(sum < q005)
    sum += clusterCount[++lo];
  hi=a;
  sum = clusterCount[a];
  while(sum < q995)
    sum += clusterCount[++hi];
  clusterStr << "99% CI for #groups = (" << lo << "," << hi << ")" << endl;
  lo=a;
  sum = clusterCount[a];
  while(sum < q025)
    sum += clusterCount[++lo];
  hi=a;
  sum = clusterCount[a];
  while(sum < q975)
    sum += clusterCount[++hi];
  clusterStr << "95% CI for #groups = (" << lo << "," << hi << ")" << endl;
  lo=a;
  sum = clusterCount[a];
  while(sum < q050)
    sum += clusterCount[++lo];
  hi=a;
  sum = clusterCount[a];
  while(sum < q950)
    sum += clusterCount[++hi];
  clusterStr << "90% CI for #groups = (" << lo << "," << hi << ")" << endl;
  clusterStr.close();

  cout << "Writing concordance factors to " << fileNames.getConcordanceFile() << "....";
  fout << "Writing concordance factors to " << fileNames.getConcordanceFile() << "....";
  ofstream concordanceStr(fileNames.getConcordanceFile().c_str());
  concordanceStr.setf(ios::fixed, ios::floatfield);
  concordanceStr.setf(ios::showpoint);

  // use TreeWeight class for sorting splits by mean number of genes
  vector<TreeWeight> sw(splits.size());
  for(int i=0;i<splits.size();i++) {
    sw[i].setWeight(0);
    for(int j=0;j<numGenes+1;j++)
      if(splitsGeneMatrix[i][j]>0)
	sw[i].addWeight(j*splitsGeneMatrix[i][j]);
    sw[i].setIndex(i);
  }
  sort(sw.begin(),sw.end(),cmpTreeWeights);

  vector<Split> ctree;
  for(int w=0;w<splits.size();w++) {
    int i = sw[w].getIndex();
    Split s(splits[i],numTaxa);
    int y = s.getClade();
    bool add=true;
    for(int j=0;j<ctree.size();j++) {
      int z = ctree[j].getClade();
      int yz = y & z;
      if( yz != z  && yz != y ) {
	add=false;
	break;
      }
    }
    if(add)
      ctree.push_back(s);
  }

  concordanceStr << "Concordance Tree Splits:" << endl;
  for(int w=0;w<ctree.size();w++) {
    ctree[w].print(concordanceStr);
    concordanceStr << endl;
  }
  concordanceStr << endl;

  for(int w=0;w<splits.size();w++) {
    int i = sw[w].getIndex();
    Split s(splits[i],numTaxa);
    s.print(concordanceStr);
    concordanceStr << endl;
    int a=numGenes,b=0;
    for(int j=0;j<numGenes+1;j++)
      if(splitsGeneMatrix[i][j]>0) {
	if(a>j)
	  a = j;
	if(b<j)
	  b = j;
      }
    concordanceStr << "#Genes      count probability cumulative" << endl;
    double csum=0;
    for(int j=a;j<b+1;j++) {
      double prob = (double)(splitsGeneMatrix[i][j])/rp.getNumUpdates();
      csum += prob;
      concordanceStr << setw(6) << j << " " << setw(10) << splitsGeneMatrix[i][j]
		     << setw(12) << setprecision(6) << prob
		     << setw(11) << setprecision(6) << csum << endl;
    }
    concordanceStr << endl;
    wsum=0;
    for(int j=a;j<b+1;j++)
      wsum += j*splitsGeneMatrix[i][j];
    concordanceStr << "mean CF = " << setprecision(3) << wsum / rp.getNumUpdates() << endl;
    int lo,hi;
    lo=a;
    sum = splitsGeneMatrix[i][a];
    while(sum < q005)
      sum += splitsGeneMatrix[i][++lo];
    hi=a;
    sum = splitsGeneMatrix[i][a];
    while(sum < q995)
      sum += splitsGeneMatrix[i][++hi];
    concordanceStr << "99% CI for CF = (" << lo << "," << hi << ")" << endl;
    lo=a;
    sum = splitsGeneMatrix[i][a];
    while(sum < q025)
      sum += splitsGeneMatrix[i][++lo];
    hi=a;
    sum = splitsGeneMatrix[i][a];
    while(sum < q975)
      sum += splitsGeneMatrix[i][++hi];
    concordanceStr << "95% CI for CF = (" << lo << "," << hi << ")" << endl;
    lo=a;
    sum = splitsGeneMatrix[i][a];
    while(sum < q050)
      sum += splitsGeneMatrix[i][++lo];
    hi=a;
    sum = splitsGeneMatrix[i][a];
    while(sum < q950)
      sum += splitsGeneMatrix[i][++hi];
    concordanceStr << "90% CI for CF = (" << lo << "," << hi << ")" << endl << endl;
  }
  concordanceStr.close();
  cout << "done." << endl;    
  fout << "done." << endl;    

  if(rp.getCalculatePairs()) {
    cout << "Writing tree pair data to " << fileNames.getPairTreeFile() << "....";
    fout << "Writing tree pair data to " << fileNames.getPairTreeFile() << "....";
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

  cout << "Writing topology single and joint posterior distribution to " << fileNames.getTreePosteriorFile() << "....";
  fout << "Writing topology single and joint posterior distribution to " << fileNames.getTreePosteriorFile() << "....";
  ofstream treePosteriorStr(fileNames.getTreePosteriorFile().c_str());
  treePosteriorStr.setf(ios::fixed, ios::floatfield);
  treePosteriorStr.setf(ios::showpoint);
  treePosteriorStr << "index "
		   << setw(max) << left << "topology" << right
		   << setw(10) << "single"
		   << setw(10) << "joint" << endl;
  for(int i=0;i<numTrees;i++) {
    double single=0.0;
    int joint=0;
    for(int j=0;j<numGenes;j++) {
      single += genes[j]->getProb(i);
      joint += newTable[i][j];
    }
    treePosteriorStr << setw(5) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right
		     << setw(10) << setprecision(6) << (double) single / (double) numGenes
		     << setw(10) << setprecision(6) << (double) joint / (double) rp.getNumUpdates() / (double) numGenes
		     << endl;
  }
  treePosteriorStr.close();
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Writing single and joint gene posteriors to " << fileNames.getGenePosteriorFile() << "....";
  fout << "Writing single and joint gene posteriors to " << fileNames.getGenePosteriorFile() << "....";
  ofstream genePostStr(fileNames.getGenePosteriorFile().c_str());
  for(int i=0;i<numGenes;i++)
    genes[i]->print(genePostStr,newTable,rp.getNumUpdates(),topologies,max);
  genePostStr.close();
  cout << "done." << endl;
  fout << "done." << endl;

  if(rp.getNumChains()>1) {
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout << "MCMCMC acceptance statistics:" << endl;
    cout << " alpha1 <--> alpha2 accepted proposed proportion" << endl;
    fout.setf(ios::fixed, ios::floatfield);
    fout.setf(ios::showpoint);
    fout << "MCMCMC acceptance statistics:" << endl;
    fout << " alpha1 <--> alpha2 accepted proposed proportion" << endl;
    for(int i=0;i<rp.getNumChains()-1;i++) {
      cout << setw(15) << alphas[i] << " <-->" << setw(15) << alphas[i+1];
      cout << setw(9) << mcmcmcAccepts[i] << setw(9) << mcmcmcProposals[i];
      fout << setw(15) << alphas[i] << " <-->" << setw(15) << alphas[i+1];
      fout << setw(9) << mcmcmcAccepts[i] << setw(9) << mcmcmcProposals[i];
      if(mcmcmcProposals[i] > 0) {
	cout << setw(11) << setprecision(6) << (double) (mcmcmcAccepts[i]) / (double) (mcmcmcProposals[i]);
	fout << setw(11) << setprecision(6) << (double) (mcmcmcAccepts[i]) / (double) (mcmcmcProposals[i]);
      }
      cout << endl;
      fout << endl;
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
  int k=readArguments(argc,argv,fileNames,mp,rp);
  int numGenes = argc-k;

  // table has one row for each input file and one column for each observed topology
  // weights for each are doubles rather than integers to allow for non-sample-based input
  vector<vector<double> > table(numGenes);
  vector<string> topologies;
  vector<string> inputFiles;
  int max=0;
  int numTaxa;

  cout << "Reading in summary files....";
  ofstream fout(fileNames.getOutFile().c_str());
  readInputFiles(argc,argv,k,table,topologies,inputFiles,max,numTaxa,fout);
  cout << "done." << endl;

  // Open outfile to save all window output
  cout << "Screen output written to file " << fileNames.getOutFile() << endl;
  intro(fout);

  fout << "Program initiated at " << ctime(&beginTime) << endl << endl;

  fout << "Program invocation:";
  for(int i=0;i<argc;i++)
    fout << " " << argv[i];
  fout << endl << endl;
  defaults.print(fout);
  fout << endl;
  printParameters(fout,fileNames,mp,rp);
  fout << endl << flush;
  fout << "Reading in summary files....";
  fout << "done." << endl;

  if(numTaxa > 32) {
    cerr << "Error: BUCKy version " << VERSION << " only works with 32 or fewer taxa." << endl;
    fout << "Error: BUCKy version " << VERSION << " only works with 32 or fewer taxa." << endl << flush;
    exit(1);
  }

  int numTrees = topologies.size();
  cout << "Read " << numGenes << " genes with a total of " << numTrees << " different sampled tree topologies" << endl;
  fout << "Read " << numGenes << " genes with a total of " << numTrees << " different sampled tree topologies" << endl;

  cout << "Writing input file names to file " << fileNames.getInputFile() << "....";
  fout << "Writing input file names to file " << fileNames.getInputFile() << "....";
  ofstream inputStr(fileNames.getInputFile().c_str());
  for(int i=0;i<inputFiles.size();i++)
    inputStr << setw(4) << i << " " << inputFiles[i] << endl;
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Sorting trees by average posterior probability....";
  fout << "Sorting trees by average posterior probability....";
  normalize(table);
  sortTrees(table,topologies);
  cout << "done." << endl;
  fout << "done." << endl;

  // Print table of individual gene posteriors

  cout << "Writing single gene posterior distribution table to file " << fileNames.getSingleFile() << "....";
  fout << "Writing single gene posterior distribution table to file " << fileNames.getSingleFile() << "....";
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

  // Initialize random number generator.

  cout << "Initializing random number generator....";
  fout << "Initializing random number generator....";
  Rand rand(rp.getSeed1(),rp.getSeed2());
  cout << "done." << endl;
  fout << "done." << endl;

  // Create genes

  cout << "Initializing gene information....";
  fout << "Initializing gene information....";
  vector<Gene*> genes(numGenes);
  vector<double> counts(numTrees);
  for(int i=0;i<numGenes;i++) {
    for(int j=0;j<numTrees;j++)
      counts[j] = table[i][j];
    genes[i] = new Gene(i,counts);
  }
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Writing summary of gene information to file " << fileNames.getGeneFile() << "....";
  fout << "Writing summary of gene information to file " << fileNames.getGeneFile() << "....";
  ofstream geneStr(fileNames.getGeneFile().c_str());
  for(int i=0;i<numGenes;i++)
    genes[i]->print(geneStr,inputFiles[i]);
  geneStr.close();
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Writing topology summary to file " << fileNames.getTopologyFile() << "....";
  fout << "Writing topology summary to file " << fileNames.getTopologyFile() << "....";
  vector<vector<int> > topologySplitsMatrix(topologies.size());
  if(numTaxa>3)
    for(int i=0;i<topologies.size();i++)
      topologySplitsMatrix[i].resize(numTaxa-3);
  vector<int> splits; // sorted list of realized splits
  ofstream topologyStr(fileNames.getTopologyFile().c_str());
  for(int i=0;i<topologies.size();i++) {
    topologyStr << setw(4) << i << " " << setw(max) << left << topologies[i].substr(0,max) << right;
    Tree t(numTaxa,topologies[i]);
    SplitSet s;
    t.getSplits(s);
    s.printShort(topologyStr);
    topologyStr << endl;

    for(int j=0;j<numTaxa-3;j++) {
      int x = s.getSplit(j);
      topologySplitsMatrix[i][j] = x;
      int n = find(splits.begin(),splits.end(),x) - splits.begin();
      if(n==splits.size())
	splits.push_back(x);
    }
  }
  topologyStr.close();
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Initializing splits counts (found " << splits.size() << " distinct splits)....";
  fout << "Initializing splits counts (found " << splits.size() << " distinct splits)....";
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

  // splitsGeneMatrix[i][j] = the number of times that split i is in exactly j genes
  vector<vector<int> > splitsGeneMatrix(splits.size());
  for(int i=0;i<splits.size();i++) {
    splitsGeneMatrix[i].resize(numGenes+1);
    for(int j=0;j<numGenes+1;j++)
      splitsGeneMatrix[i][j] = 0;
  }
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Writing splits key to file " << fileNames.getSplitsFile() << "....";
  fout << "Writing splits key to file " << fileNames.getSplitsFile() << "....";
  ofstream splitsStr(fileNames.getSplitsFile().c_str());
  for(int i=0;i<splits.size();i++) {
    splitsStr << setw(3) << splits[i] << " ";
    Split s(splits[i],numTaxa);
    s.print(splitsStr);
    splitsStr << endl;
  }
  splitsStr.close();
  cout << "done." << endl;
  fout << "done." << endl;

  if(rp.getNumUpdates()==0)
    return 0;

  cout << "Setting initial MCMC state....";
  fout << "Setting initial MCMC state....";

  vector<double> alphas(rp.getNumChains());
  alphas[0] = mp.getAlpha();
  for(int i=1;i<rp.getNumChains();i++)
    alphas[i] = rp.getAlphaMultiplier()*alphas[i-1];

  vector<State*> states(rp.getNumChains());
  vector<int> index(rp.getNumChains());
  for(int i=0;i<rp.getNumChains();i++) {
    states[i] = new State(alphas[i],numTaxa,numTrees,genes,mp.getUseIndependencePrior(),rand);
    index[i] = i;
  }
  cout << "done." << endl;
  fout << "done." << endl;

  cout << "Initializing MCMCMC acceptance counters and pairwise counters....";
  fout << "Initializing MCMCMC acceptance counters and pairwise counters....";
  vector<int> mcmcmcAccepts(rp.getNumChains());
  vector<int> mcmcmcProposals(rp.getNumChains());
  for(int i=0;i<rp.getNumChains();i++)
    mcmcmcAccepts[i] = mcmcmcProposals[i] = 0;

  vector<vector<int> > pairCounts(numGenes);
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
  cout << "+----+----+----+----+----+----+----+----+----+----+" << endl;
  fout << "Beginning burn-in with " << numBurn << " updates (10% extra of desired updates)..." << endl;
  fout << "0   10   20   30   40   50   60   70   80   90   100" << endl;
  fout << "+----+----+----+----+----+----+----+----+----+----+" << endl;
  int part = numBurn / 50;
  for(int cycle=0;cycle<numBurn;cycle++) {
    if( cycle % part == 0) {
      cout << "*" << flush;
      fout << "*" << flush;
    }
    for(int i=0;i<rp.getNumChains();i++) {
      states[i]->update(rand);
      //      if(rp.getUseUpdateGroups()) {
      //	int group = (int)(rand.runif()*states[i]->getNumGroups());
      //	states[i]->updateOneGroup(group,rand);
      //      }
      if(rp.getUseUpdateGroups()) {
	int gene = (int)(rand.runif()*genes.size());
	states[i]->updateOneGroupNew(gene,rand);
      }
      
    }

    if(cycle % rp.getMCMCMCRate() == 0 && rp.getNumChains()>1) {
      mcmcmc(states,index,alphas,rand,mcmcmcAccepts,mcmcmcProposals);
    }
  }
  cout << "*" << endl;
  cout << " ....done." << endl;
  fout << "*" << endl;
  fout << " ....done." << endl;

  cout << "Initializing summary tables...";
  fout << "Initializing summary tables...";

  vector<int> clusterCount(numGenes+1);
  for(int i=0;i<clusterCount.size();i++)
    clusterCount[i] = 0;

  vector<vector<int> > newTable(numTrees);
  for(int i=0;i<numTrees;i++) {
    newTable[i].resize(numGenes);
    for(int j=0;j<numGenes;j++)
      newTable[i][j] = 0;
  }

  cout << "done." << endl;
  fout << "done." << endl;

  if(rp.getCreateSampleFile()) {
    cout << "Sampled topologies will be in file " << fileNames.getSampleFile() << "." << endl;
    fout << "Sampled topologies will be in file " << fileNames.getSampleFile() << "." << endl;
  }

  cout << "Beginning " << rp.getNumUpdates() << " MCMC updates..." << endl;
  cout << "0   10   20   30   40   50   60   70   80   90   100" << endl;
  cout << "+----+----+----+----+----+----+----+----+----+----+" << endl;
  fout << "Beginning " << rp.getNumUpdates() << " MCMC updates..." << endl;
  fout << "0   10   20   30   40   50   60   70   80   90   100" << endl;
  fout << "+----+----+----+----+----+----+----+----+----+----+" << endl;
  part = rp.getNumUpdates() / 50;
  ofstream sampleFileStr(fileNames.getSampleFile().c_str());
  sampleFileStr.setf(ios::fixed, ios::floatfield);
  sampleFileStr.setf(ios::showpoint);
  vector<int> accept(rp.getNumChains());
  for(int i=0;i<rp.getNumChains();i++)
    mcmcmcAccepts[i] = mcmcmcProposals[i] = accept[i] = 0;

  for(int cycle=0;cycle<rp.getNumUpdates();cycle++) {
    for(int i=0;i<rp.getNumChains();i++) {
      accept[index[i]] += states[i]->update(rand);
      if(rp.getUseUpdateGroups()) {
	int group = (int)(rand.runif()*states[i]->getNumGroups());
	states[i]->updateOneGroup(group,rand);
      }
    }
    if(cycle % rp.getMCMCMCRate() == 0 && rp.getNumChains()>1)
      mcmcmc(states,index,alphas,rand,mcmcmcAccepts,mcmcmcProposals);
    // update counts
    int i0 = index[0];
    states[i0]->updateTable(newTable);
    states[i0]->updateSplits(splitsGeneMatrix,topologySplitsIndexMatrix);
    clusterCount[states[i0]->getNumGroups()]++;
    if( cycle % part == 0) {
      cout << "*";
      cout.flush();
      fout << "*";
      fout.flush();
    }
    if( rp.getCalculatePairs() && cycle % rp.getSubsampleRate() == 0)
      states[i0]->updatePairCounts(pairCounts);
    if( rp.getCreateSampleFile() && cycle % rp.getSubsampleRate() == 0) {
      sampleFileStr << setw(8) << accept[0];
      accept[0] = 0;
      states[i0]->sample(sampleFileStr);
    }
  }
  cout << "*" << endl;
  cout << " ....done." << endl;
  fout << "*" << endl;
  fout << " ....done." << endl;

  writeOutput(fout,fileNames,max,numTrees,numTaxa,topologies,numGenes,rp,newTable,clusterCount,splits,splitsGeneMatrix,
	      pairCounts,genes,alphas,mcmcmcAccepts,mcmcmcProposals);

  for(int i=0;i<numGenes;i++)
    delete genes[i];

  for(int i=0;i<rp.getNumChains();i++)
    delete states[i];

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
  fout << seconds << (seconds==1 ? " second." : " seconds.") << endl;

  fout.close();

  return 0;
}
