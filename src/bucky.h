#ifndef MCMCHDR
#define MCMCHDR

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>

using namespace std;

// General function to return log( prod_{i=0}^{n-1} (i+x) )
//   which is used in the new model repeatedly.
//   Assume that x > 0.
double logA(double x,int n)
{
  if(n<0) {
    cerr << "Error: logA called with negative n = " << n << endl;
    exit(1);
  }
  if(n==0)
    return 0.0;
  double answer = log(x);
  for(int i=1;i<n;i++)
    answer += log(i + x);
  return answer;
}

class Rand {
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000, 2003  The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* A version of Marsaglia-MultiCarry */

 public:

  Rand(unsigned int seed1=1234, unsigned int seed2=5678) : I1(seed1), I2(seed2) {}

  void setSeed(unsigned int i1=1234, unsigned int i2=5678) { I1 = i1; I2 = i2; }

  void getSeed(unsigned int& i1, unsigned int& i2) { i1 = I1; i2 = I2; }

  double runif() {
    // Returns a pseudo-random number between 0 and 1.

    I1= 36969*(I1 & 0177777) + (I1>>16);
    I2= 18000*(I2 & 0177777) + (I2>>16);
    return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
  }
 private:
  unsigned int I1, I2;
};

class Gene {
public:
  Gene(int n,vector<double> x) : number(n) {
    numTrees = 0;
    total = 0.0;
    for(int i=0;i<x.size();i++)
      if(x[i]>0) {
	numTrees++;
	total += x[i];
	indices.push_back(i);
	counts.push_back(x[i]);
      }
  }
  int getNumber() const { return number; }
  int getNumTrees() const { return numTrees; }
  double getTotal() const { return total; }
  int getIndex(int i) const { return indices[i]; }
  double getCount(int i) const { return counts[i]; }
  int pickTree(Rand& r) const {
    double x = total * r.runif();
    int i=0;
    while(x>counts[i])
      x -= counts[i++];
    return indices[i];
  }
  double getProb(int top) const {
    int i=0;
    while(indices[i]!=top && i<indices.size())
      i++;
    if(i==indices.size())
      return 0.0;
    else
      return (double)(counts[i])/(double)(total);
  }
  void print(ostream& f,string file) {
    f.setf(ios::fixed, ios::floatfield);
    f.setf(ios::showpoint);

    f << "Gene " << number << " (from file " << file << "):" << endl;
    f << "  numTrees = " << numTrees << endl;
    f << "    i indices probability" << endl;
    for(int i=0;i<numTrees;i++)
      f << setw(5) << i << ":" << setw(7) << indices[i] << setw(12) << setprecision(8) << counts[i] << endl;
    f << endl;
  }
  void print(ostream& f,const vector<vector<int> > &newTable,int numCycles,const vector<string> &topologies,int max) {
    f.setf(ios::fixed, ios::floatfield);
    f.setf(ios::showpoint);
    f << "Gene " << number << ":" << endl;
    f << "  numTrees = " << numTrees << endl;
    f << "  index "
      << setw(max) << left << "topology" << right
      << setw(10) << "single"
      << setw(10) << "joint" << endl;
    for(int i=0;i<numTrees;i++)
      f << setw(7) << indices[i] << " "
	<< setw(max) << left << topologies[indices[i]].substr(0,max) << right
	<< setw(10) << setprecision(6) << counts[i] / total
	<< setw(10) << setprecision(6) << (double) newTable[indices[i]][number] / (double) numCycles << endl;
    f << endl;
  }
private:
  int number;
  int numTrees;
  double total;
  vector<int> indices;
  vector<double> counts;
};

class State {
public:
  State(double a,int numTaxa,int numTrees,vector<Gene*> g,bool useIP,Rand& rand) {
    alpha = a;
    logAlpha = log(alpha);
    numTreesSampled = numTrees;
    logNumTopologies = 0;
    for(int i=4;i<=numTaxa;i++)
      logNumTopologies += log(2.0*i-5);
    genes = g;
    alphaOverTop = exp(logAlpha - logNumTopologies);
    useIndependencePrior = useIP;
    numGroups = 0;
    counts.resize(numTrees);
    for(int i=0;i<numTrees;i++)
      counts[i] = 0;
    for(int i=0;i<genes.size();i++) {
      int top;
      top = genes[i]->pickTree(rand);
      if(counts[top]==0) {
	numGroups++;
	indices.push_back(top);
      }
      counts[top]++;
      tops.push_back(top);
    }
    int numDataGenes = genes.size();

    sort(indices.begin(),indices.end());
    calculateLogPriorProb();
    calculateLogPosteriorProbProduct();
  }
  double getAlpha() { return alpha; }
  int getNumGroups() { return numGroups; }
  double getLogNumTopologies() { return logNumTopologies; }
  double getAlphaOverTop() { return alphaOverTop; }
  int getCounts(int i) { return counts[i]; }
  int getTops(int i) { return tops[i]; }
  void setAlpha(double a) {
    alpha = a;
    logAlpha = log(a);
    alphaOverTop = exp(logAlpha - logNumTopologies);
    calculateLogPriorProb();
  }
  void setTop(int i,int newTop) {
    int oldTop = tops[i];
    int n1 = counts[oldTop];
    int n2 = counts[newTop];
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
    double logHR=0;
    logHR += log(n2 + alphaOverTop);
    logHR -= log(n1-1 + alphaOverTop);
    logPriorProb += logHR;
    if(i < genes.size())
      logPosteriorProbProduct += (log(genes[i]->getProb(newTop)) - log(genes[i]->getProb(oldTop)));

  }
  double getLogPriorProb() { return logPriorProb; }
  double calculateLogPriorProb();
  double getLogPosteriorProbProduct() { return logPosteriorProbProduct; }
  double calculateLogPosteriorProbProduct() {
    logPosteriorProbProduct = 0;
    for(int i=0;i<genes.size();i++)
      logPosteriorProbProduct += log(genes[i]->getProb(tops[i]));
    return logPosteriorProbProduct;
  }
  double getLogH() { return logPriorProb + logPosteriorProbProduct; }
  int updateOne(int,Rand&);
  int update(Rand&);
  void sample(ostream& f) {
    f.setf(ios::fixed, ios::floatfield);
    f.setf(ios::showpoint);
    f << setw(6) << numGroups << setw(12) << setprecision(6) << getLogH();
    for(int i=0;i<genes.size();i++)
      f << setw(4) << tops[i];
    f << endl;
  }
  void updateTable(vector<vector<int> >& table) {
    for(int j=0;j<genes.size();j++)
      table[tops[j]][j]++;
  }
  void print(ostream&);
  void updateSplits(vector<vector<int> >&,vector<vector<int> >&);
  void updateOneGroup(int,Rand&);
  void updateOneGroupNew(int,Rand&);
  void updatePairCounts(vector<vector<int> >&);
 private:
  double alpha,logAlpha,alphaOverTop;
  bool useIndependencePrior;
  vector<Gene*> genes; // vector of pointers to genes
  int numGroups,numTreesSampled;
  double logNumTopologies;
  vector<int> indices; // indices[i] = index of ith positive entry in counts
  vector<int> counts; // counts[i] = number of genes with topology i
  vector<int> tops; // tops[i] = topology for gene i
  double logPriorProb;
  double logPosteriorProbProduct;
};

class Edge;

class Node {
public:
  Node(int n,int lf) : number(n) {
    taxon = 1;
    for(int i=1;i<n;i++)
      taxon *= 2;
    for(int i=0;i<3;i++)
      taxa[i] = 0;
    if(lf)
      leaf = true;
    else
      leaf = false;
  }
  int getNumber() const { return number; }
  Edge* getEdge(int n) const { return edges[n]; }
  bool isLeaf() const { return leaf; } 
  int getTaxa(int i) const { return taxa[i]; }
  void setEdge(int n,Edge *e) { edges[n]=e; }
  void setTaxa(int i,int x) { taxa[i] = x; }
  int setAllTaxa(Node*,unsigned int);
  void print(ostream&) const;
  Node* getNeighbor(int) const;
private:
  int number;
  int taxon; // 2^number for leaves, used as a bit vector
  int taxa[3]; // integer representing bit vector of taxa in each of three subtrees for internal nodes
  bool leaf;
  Edge* edges[3];
};

class Edge {
public:
  Edge(int n) : number(n) {}
  int getNumber() const { return number; }
  Node* getNode(int i) { return nodes[i]; }
  Node* getOtherNode(const Node *n) const { return (n==nodes[0] ? nodes[1] : nodes[0]); }
  int getSplit() { return split; }
  void setSplit(int s) {split = s; }
  void setNode(int i,Node *n) { nodes[i] = n; }
  void print(ostream&,int) const;
private:
  int number;
  int split;
  Node* nodes[2];
};

class Split;
class SplitSet;

class Tree {
public:
  Tree(int n,string top) : numTaxa(n), numNodes(2*n - 2), numEdges(2*n - 3), numSplits(n - 3) {
    int k=1;
    all = 0;
    for(int i=0;i<numTaxa;i++) {
      all += k;
      k *= 2;
    }
    nodes.resize(numNodes);
    edges.resize(numEdges);
    for(int i=0;i<numTaxa;i++)
      nodes[i] = new Node(i,1);
    for(int i=numTaxa;i<numNodes;i++)
      nodes[i] = new Node(i,0);
    for(int i=0;i<numEdges;i++)
      edges[i] = new Edge(i);
    connect(top);
    setAllTaxa();
  }
  ~Tree() {
    for(int i=0;i<numNodes;i++)
      delete nodes[i];
    for(int i=0;i<numEdges;i++)
      delete edges[i];
  }    
  void connect(string);
  Node* connectInt(istream&,int&,int&,int&);
  Node* connectThreeNodes(Node*,Node*,Node*,int&,int&);
  void connectTwoNodes(Node*,Node*,int&,int&);
  void print(ostream&) const;
  void setAllTaxa();
  void getSplits(SplitSet&);
  int getNumTaxa() const { return numTaxa; }
  int getNumNodes() const { return numNodes; }
  int getNumEdges() const { return numEdges; }
  int getNumSplits() const { return numSplits; }
  int getAll() const { return all; }
  Node* getNode(int i) const { return nodes[i]; }
  Edge* getEdge(int i) const { return edges[i]; }
  string getTop() const { return top; }
private:
  int numTaxa,numNodes,numEdges,numSplits,all;
  vector<Node *> nodes;
  vector<Edge *> edges;
  string top;
};

class Split {
public:
  Split() {}
  Split(int c,int x) : clade(c), numTaxa(x) {}
  void setClade(int x) { clade = x; }
  int getClade() { return clade; }
  void setNumTaxa(int x) { numTaxa = x; }
  int getNumTaxa() { return numTaxa; }
  void print(ostream& f) {
    vector<int> left;
    vector<int> right;
    int c = clade;
    for(int i=1;i<=numTaxa;i++) {
      if(c % 2 == 1)
	left.push_back(i);
      else
	right.push_back(i);
      c /= 2;
    }
    f << "{";
    if(left.size()>0)
      f << left[0];
    for(int i=1;i<left.size();i++)
      f << "," << left[i];
    f << "|";
    if(right.size()>0)
      f << right[0];
    for(int i=1;i<right.size();i++)
      f << "," << right[i];
    f << "}";
  }
private:
  int clade,numTaxa,index;
};

class SplitSet {
 public:
  SplitSet() {}
  void print(ostream& f) {
    for(int i=0;i<splits.size();i++) {
      splits[i].print(f);
      f << endl;
    }
  }
  void printShort(ostream& f) {
    for(int i=0;i<splits.size();i++)
      f << setw(5) << splits[i].getClade() << " ";
  }
  int getNumTaxa() { return numTaxa; }
  int getSplit(int i) { return splits[i].getClade(); }
  void setNumTaxa(int x) { numTaxa = x; }
  void setAll(int a) { all = a; }
  void setSplit(int i,int c) {
    splits[i].setClade(c);
    splits[i].setNumTaxa(numTaxa);
  }
  void copySplits(Tree& t) { t.getSplits(*this); }
  void resize(int n) { splits.resize(n); }
 private:
  int all,numTaxa;
  vector<Split> splits;
};

class FileNames {
public:
  FileNames(string fileRoot) :
    rootFileName(fileRoot),
    clusterFile(fileRoot + ".cluster"),
    concordanceFile(fileRoot + ".concordance"),
    geneFile(fileRoot + ".gene"),
    genePosteriorFile(fileRoot + ".genepost"),
    jointFile(fileRoot + ".joint"),
    inputFile(fileRoot + ".input"),
    outFile(fileRoot + ".out"),
    pairTreeFile(fileRoot + ".pairs"),
    sampleFile(fileRoot + ".sample"),
    singleFile(fileRoot + ".single"),
    splitsFile(fileRoot + ".splits"),
    topologyFile(fileRoot + ".top"),
    treePosteriorFile(fileRoot + ".topologies")
  {}
  void setFileNames(string fileRoot) {
    clusterFile = fileRoot + ".cluster";
    concordanceFile = fileRoot + ".concordance";
    geneFile = fileRoot + ".gene";
    genePosteriorFile = fileRoot + ".genepost";
    jointFile = fileRoot + ".joint";
    inputFile = fileRoot + ".input";
    outFile = fileRoot + ".out";
    pairTreeFile = fileRoot + ".pairs";
    sampleFile = fileRoot + ".sample";
    singleFile = fileRoot + ".single";
    splitsFile = fileRoot + ".splits";
    topologyFile = fileRoot + ".top";
    treePosteriorFile = fileRoot + ".topologies";
  }
  string getRootFileName() { return rootFileName; }
  string getClusterFile() { return clusterFile; }
  string getConcordanceFile() { return concordanceFile; }
  string getGeneFile() { return geneFile; }
  string getGenePosteriorFile() { return genePosteriorFile; }
  string getJointFile() { return jointFile; }
  string getInputFile() { return inputFile; }
  string getOutFile() { return outFile; }
  string getPairTreeFile() { return pairTreeFile; }
  string getSampleFile() { return sampleFile; }
  string getSingleFile() { return singleFile; }
  string getSplitsFile() { return splitsFile; }
  string getTopologyFile() { return topologyFile; }
  string getTreePosteriorFile() { return treePosteriorFile; }
private:
  string rootFileName;
  string clusterFile;
  string concordanceFile;
  string geneFile;
  string genePosteriorFile;
  string jointFile;
  string inputFile;
  string outFile;
  string pairTreeFile;
  string sampleFile;
  string singleFile;
  string splitsFile;
  string topologyFile;
  string treePosteriorFile;
};

class ModelParameters {
public:
  ModelParameters() : alpha(1.0), useIndependencePrior(false) {}
  double getAlpha() { return alpha; }
  void setAlpha(double a) { alpha = a; }
  bool getUseIndependencePrior() { return useIndependencePrior; }
  void setUseIndependencePrior(bool iprior) { useIndependencePrior = iprior; }
private:
  double alpha;
  bool useIndependencePrior;
};

class RunParameters {
public:
  RunParameters() : 
    alphaMultiplier(10.0),
    seed1(1234),
    seed2(5678),
    numUpdates(100000),
    subsampleRate(1),
    numChains(1),
    mcmcmcRate(100),
    calculatePairs(false),
    useUpdateGroups(true),
    createSampleFile(false)
  {}
  double getAlphaMultiplier() { return alphaMultiplier; }
  void setAlphaMultiplier(double x) { alphaMultiplier = x; }
  unsigned int getSeed1() { return seed1; }
  void setSeed1(unsigned int x) { seed1 = x; }
  unsigned int getSeed2() { return seed2; }
  void setSeed2(unsigned int x) { seed2 = x; }
  unsigned int getNumUpdates() { return numUpdates; }
  void setNumUpdates(unsigned int x) { numUpdates = x; }
  unsigned int getSubsampleRate() { return subsampleRate; }
  void setSubsampleRate(unsigned int x) { subsampleRate = x; }
  unsigned int getNumChains() { return numChains; }
  void setNumChains(unsigned int x) { numChains = x; }
  unsigned int getMCMCMCRate() { return mcmcmcRate; }
  void setMCMCMCRate(unsigned int x) { mcmcmcRate = x; }
  bool getCalculatePairs() { return calculatePairs; }
  void setCalculatePairs(bool x) { calculatePairs = x; }
  bool getUseUpdateGroups() { return useUpdateGroups; }
  void setUseUpdateGroups(bool x) { useUpdateGroups = x; }
  bool getCreateSampleFile() { return createSampleFile; }
  void setCreateSampleFile(bool x) { createSampleFile = x; }
private:
  double alphaMultiplier;
  unsigned int seed1,seed2,numUpdates,subsampleRate,numChains,mcmcmcRate;
  bool calculatePairs;
  bool useUpdateGroups;
  bool createSampleFile;
};

class Defaults {
public:
  Defaults(FileNames &fn,ModelParameters &mp,RunParameters &rp) :
    alpha(mp.getAlpha()),
    numUpdates(rp.getNumUpdates()),
    numChains(rp.getNumChains()),
    mcmcmcRate(rp.getMCMCMCRate()),
    alphaMultiplier(rp.getAlphaMultiplier()),
    subsampleRate(rp.getSubsampleRate()),
    rootFileName(fn.getRootFileName()),
    seed1(rp.getSeed1()),
    seed2(rp.getSeed2()),
    useIndependencePrior(mp.getUseIndependencePrior()),
    calculatePairs(rp.getCalculatePairs()),
    useUpdateGroups(rp.getUseUpdateGroups()),
    createSampleFile(rp.getCreateSampleFile())
  {}
  void print(ostream&);
private:
  double alpha;
  unsigned int numUpdates;
  unsigned int numChains;
  unsigned int mcmcmcRate;
  double alphaMultiplier;
  unsigned int subsampleRate;
  string rootFileName;
  unsigned int seed1;
  unsigned int seed2;
  bool useIndependencePrior;
  bool calculatePairs;
  bool useUpdateGroups;
  bool createSampleFile;
};

class TreeWeight {
public:
  TreeWeight() {}
  void setWeight(double x) { weight=x; }
  void addWeight(double x) { weight += x; }
  double getWeight() { return weight; }
  void setIndex(int x) { index=x; }
  int getIndex() { return index; }
private:
  double weight;
  int index;
};

bool cmpTreeWeights(TreeWeight x,TreeWeight y) { return x.getWeight() > y.getWeight(); }

#endif

