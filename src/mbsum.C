// BUCKy - Bayesian Untangling of Concordance Knots (applied to yeast and other organisms)
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

// version 1.2b: multiple files from the same gene can be combined.
// version 1.3.0: Translate table, if present, are copied verbatim from
//              the first file, assumed the same in other files.

// mbsum.C

// Input:    one or more MrBayes .t files
//           Lines with the format `tree <name> = <parenthetic tree representation>'
//             are read as trees.  Other lines are ignored.
//
// Options:  [-n number-of-skipped-trees] will skip over the first number-of-skipped-trees trees in each input file
//           [--skip number-of-skipped-trees] will skip over the first number-of-skipped-trees trees in each input file
//           [-o outfile] name of file for writing output
//           [-h] prints the brief usage message
//           [--help] prints the brief usage message
//           [--version] prints the version number
//
// Output:   a single file with the name <filename - .extension> + .in (by default)
//             (Output file name is the original file name, stripped of the last extension if there is one, plus .in .)
//           The output file contains one line per sampled tree topology with two fields,
//             the topology (parenthetic representation) and the count.
//           Tree topologies are sorted by count, from highest to lowest.
//           The output file is in the appropriate input format for BUCKy.
//
// Usage:    mbsum [--help || -h] [<--skip || -n> number-of-skipped-trees] [<--out || -o> output-file] [--version] [input filename(s)]

#define VERSION "1.3.0"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

class Node;

class Edge {
public:
  Edge() {}
  ~Edge() {
    for(int i=0;i<2;i++)
      nodes[i] = NULL;
  }
  Edge(int n,double x) :
    number(n),
    length(x)
  {}
  int getNumber() { return number; }
  void setLength(double x) { length=x; }
  double getLength() { return length; }
  Node* getOtherNode(Node *n) { return (nodes[0] == n ? nodes[1] : nodes[0]); }
  void setNode(Node *n,int i) { nodes[i]=n; }
  void setNodes(Node *m,Node *n) { nodes[0]=m; nodes[1]=n; }
  Node* getNode(int n) { return nodes[n]; }
  void swapNodes() {
    Node *temp=nodes[0];
    nodes[0] = nodes[1];
    nodes[1] = temp;
  }
private:
  int number;
  double length;
  Node* nodes[2];
};

class Node {
 public:
  Node() {}
  ~Node() { edges.clear(); }
  Node(int n,bool l) : number(n), leaf(l) {}
  int getNumber() const { return number; }
  void setNumber(int x) { number = x; }
  bool isLeaf() { return leaf; }
  void setLeaf(bool x) { leaf = x; }
  int getNumEdges() { return edges.size(); }
  Edge* getEdge(int n) { return edges[n]; }
  void setEdge(Edge* e) { edges.push_back(e); }
  void setEdge(Edge *e,int n) { edges[n] = e; }
  void reorder(Edge*);
  void rootReorder();
  int getMinTaxa() { return minTaxa; }
  int setMinTaxa(int x) { minTaxa = x; return x; }
  int setMinTaxa(Edge*);
  void setNumbers(int&,Edge*);
  void printTop(ostream&,Edge*);
 private:
  int number;
  bool leaf;
  vector<Edge*> edges;
  int minTaxa;
};

class EdgeCompare : public binary_function<Edge*,Edge*,bool> {
public:
  EdgeCompare(Edge* e,Node* n) : parent(e),node(n) {}

  bool operator()(Edge* x,Edge* y) {
    if(x==parent)
      return false;
    if(y==parent)
      return true;
    return (x->getOtherNode(node)->getMinTaxa() < y->getOtherNode(node)->getMinTaxa());
  }
private:
  const Edge* parent;
  Node* node;
};

class RootEdgeCompare : public binary_function<Edge*,Edge*,bool> {
public:
  RootEdgeCompare(Node* r) : root(r) {}

  bool operator()(Edge* x,Edge* y) {
    return (x->getOtherNode(root)->getMinTaxa() < y->getOtherNode(root)->getMinTaxa());
  }
private:
  Node* root;
};

void Node::rootReorder() {
  // only call on the root!
  if(leaf)
    return;

  sort(edges.begin(),edges.end(),RootEdgeCompare(this));
}

void Node::reorder(Edge* par) {
  // sort edges by order of minTaxa of other nodes with parent edge coming last
  // recurse to children
  if(leaf)
    return;
  for(int i=0;i<edges.size();i++) {
    Edge* e = edges[i];
    if(e != par)
      e->getOtherNode(this)->reorder(e);
  }
  sort(edges.begin(),edges.end(),EdgeCompare(par,this));
}

class Tree {
public:
  Tree(string,int);
  ~Tree() {
    for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
      delete *e;
    for(vector<Node*>::iterator n=nodes.begin();n!=nodes.end();n++)
      delete *n;
    nodes.clear();
    edges.clear();
  }
  int getNumTaxa() { return numTaxa; }
  int getNumNodes() { return numNodes; }
  int getNumEdges() { return numEdges; }
  Node* getNode(int n) { return nodes[n]; }
  Edge* getEdge(int n) { return edges[n]; }
  void printTop(ostream&);
  void reorder(Node*);
  void readSubtree(istringstream&,Node*,int&,int);
  void setNumbers();
  void setMinTaxa(Node*);
private:
  int numTaxa;
  int numNodes;
  int numEdges;
  vector<Node*> nodes;
  vector<Edge*> edges;
};

bool cmpNodes(Node* x,Node* y) { return( x->getNumber() < y->getNumber() ); }

void Node::setNumbers(int& nextNumber,Edge* parentEdge) {
  // Set numbers for internal nodes for subtrees along non-parent edges.
  // Numbers for leaves should already be set.
  if(leaf)
    return;
  for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
    if(*e != parentEdge)
      (*e)->getOtherNode(this)->setNumbers(nextNumber,*e);
  setNumber(nextNumber++);
}

void Tree::setNumbers() {
  // Internal nodes are not numbered when tree is read.
  // After reading the tree, we can give numbers to the internal nodes.
  int nextNumber = numTaxa+1;
  for(int i=0;i<nodes[0]->getNumEdges();i++) {
    Edge* e = nodes[0]->getEdge(i);
    e->getOtherNode(nodes[0])->setNumbers(nextNumber,e);
  }
  nodes[0]->setNumber(nextNumber);
}

int Node::setMinTaxa(Edge* parentEdge) {
  // assumes that number is set for node and that it is larger than leaf numbers
  if(leaf) {
    setMinTaxa(number);
    return number;
  }
  int m = number,n;
  for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
    if(*e != parentEdge) {
      n = (*e)->getOtherNode(this)->setMinTaxa(*e);
      if(n < m)
	m = n;
    }
  setMinTaxa(m);
  return m;
}

void Tree::setMinTaxa(Node* root) {
  // assumes that root has the largest number
  int m = root->getNumber(),n;
  for(int i=0;i<root->getNumEdges();i++) {
    Edge* e = root->getEdge(i);
    n = e->getOtherNode(root)->setMinTaxa(e);
    if(n < m)
      m = n;
  }
  root->setMinTaxa(m);
}

void Tree::reorder(Node* root) {
  // reorder edges in all nodes, last edge for root is still the last
  for(int i=0;i<root->getNumEdges();i++) {
    Edge* e = root->getEdge(i);
    e->getOtherNode(root)->reorder(e);
  }
  root->rootReorder();
}

void usageError(char *name)
{
  cerr << "Usage: " << name << " mbsum [--help || -h] [<--skip || -n> number-of-skipped-trees] [<--out || -o> output-file] [--version] [input filename(s)]" << endl;
  exit(1);
}

void advance(istringstream& s,char& c,int numLeft,int lineNumber) {
  c = s.peek();
  if(s.eof() || s.fail()) {
    cerr << "Error: Tree on line " << lineNumber
	 << " ends prematurely while reading subtree after left parenthesis " << numLeft << "." << endl;
    exit(1);
  }
  s >> c;
}

void Tree::readSubtree(istringstream& s,Node* parent,int& numLeft,int numLine)
{
  Node* n = new Node();
  nodes.push_back(n);
  numNodes++;
  Edge* e = new Edge();
  edges.push_back(e);
  numEdges++;
  e->setNodes(n,parent);
  n->setEdge(e);
  parent->setEdge(e);
  char c = s.peek();
  if(c=='(') {
    n->setLeaf(false);
    s >> c;
    numLeft++;
    readSubtree(s,n,numLeft,numLine);
    advance(s,c,numLeft,numLine);
    if(c == ',')
      readSubtree(s,n,numLeft,numLine);
    else {
      cerr << "Error: Tree on line " << numLine << " after left parenthesis " << numLeft
	   << " expected comma, but read " << c << "." << endl;
      string a;
      s >> a;
      cerr << "Remainder of line is /" << a << "/." << endl;
    }
    advance(s,c,numLeft,numLine);
    if(c != ')') {
      cerr << "Error: Tree on line " << numLine << " after left parenthesis " << numLeft
	   << " expected right parenthesis, but read " << c << "." << endl;
      string a;
      s >> a;
      cerr << "Remainder of line is /" << a << "/." << endl;
      exit(1);
    }
  }
  else if(isdigit(c)) {
    int num;
    s >> num;
    n->setNumber(num);
    n->setLeaf(true);
    numTaxa++;
  }
  else {
    cerr << "Error: Expected beginning of subtree, but read '" << c << "' on line " << numLine << "." << endl;
    exit(1);
  }
  c = s.peek();
  double length = 1.0;
  if(c == ':') {
    s >> c;
    s >> length;
    if(s.fail()) {
      cerr << "Error: Tree on line " << numLine << " after left parenthesis " << numLeft
	   << " expected an edge length after reading colon." << endl;
      string a;
      s >> a;
      cerr << "Remainder of line is /" << a << "/." << endl;
      exit(1);
    }
  }
  e->setLength(length);
}

Tree::Tree(string line,int lineNumber) 
{
  // read in the tree from parenthetic representation.
  // Create new nodes and edges on the fly.
  // Then renumber.
  // Then, reorder so that leaves come first.
  
  numEdges = 0;
  numTaxa = 0;

  nodes.push_back(new Node()); // root of the tree, thinks it is a leaf....
  numNodes = 1;

  int numLeft = 0;
  istringstream s(line);
  char c = s.peek();
  if(c =='(') {
    numLeft++;
    // root is not a leaf
    nodes[0]->setLeaf(false);
    s >> c; // read in left parenthesis
    while(true) {
      readSubtree(s,nodes[0],numLeft,lineNumber);
      s >> c;
      if(s.eof() || s.fail()) {
	cerr << "Error: Tree on line " << lineNumber << " ended while reading in subtree beginning with left parenthesis number "
	     << numLeft << "." << endl;
	exit(1);
      }
      else if(c==',') // read in another subtree to the root.
	continue;
      else if(c==')') { // base subtree is finished.
	c = s.peek();
	if(c==';')
	  s >> c;
	break;
      }
    }
  }
  else if(isdigit(c)) { // base tree may just be a single node
    int n;
    s >> n;
    nodes[0]->setLeaf(true);
    nodes[0]->setNumber(n);
    numTaxa++;
    if(c==';')
      s >> c;
    s >> c;
    if(!(s.eof())) {
      cerr << "Error: Tree on line  " << lineNumber << " begins with a number, but is not a one-taxon tree." << endl;
      exit(1);
    }
  }
  // give numbers to internal nodes
  setNumbers();
  // set minTaxa
  setMinTaxa(nodes[0]);
  // sort nodes
  sort(nodes.begin(),nodes.end(),cmpNodes);
  // reorder edges in nodes
  reorder(nodes[numNodes-1]);
}

void Node::printTop(ostream& f,Edge* parent) {
  if(isLeaf()) {
    f << number;
    return;
  }
  int remainingCommas = getNumEdges()-2;
  f << "(";
  for(vector<Edge*>::iterator e=edges.begin();e!=edges.end();e++)
    if(*e != parent) {
      (*e)->getOtherNode(this)->printTop(f,*e);
      if(remainingCommas-- > 0)
	f << ",";
    }
  f << ")";
}

void Tree::printTop(ostream& f) {
  Node* root = nodes[numNodes-1];
  if(numNodes==1) {
    f << root->getNumber() << ";" << endl;
    return;
  }
  f << "(";
  int degree = root->getNumEdges();
  if(degree==2) {
    for(int i=0;i<degree;i++) {
      Edge* e = root->getEdge(i);
      e->getOtherNode(root)->printTop(f,e);
      if(i < degree-1)
	f << ",";
    }
  }
  else { // degree is 3 or more
    // combine last degree-1 subtrees 
    // so as to root the tree --required for bucky
    Edge* first = root->getEdge(0);
    first->getOtherNode(root)->printTop(f,first);
    f << ",(";
    for(int i=1;i<degree;i++) {
      Edge* e = root->getEdge(i);
      e->getOtherNode(root)->printTop(f,e);
      if(i < degree-1)
	f << ",";
    }
    f << ")";
  }
  f << ");";
}

class TopCountNode {
public:
  TopCountNode() {}
  TopCountNode(string t) : count(1), top(t), left(NULL), right(NULL), parent(NULL) {} // create the root node
  TopCountNode(string t,TopCountNode* p) : count(1), top(t), left(NULL), right(NULL), parent(p) {} // create a child node
  ~TopCountNode() {
    if(left != NULL) {
      left->~TopCountNode();
      delete left;
    }
    if(right != NULL) {
      right->~TopCountNode();
      delete right;
    }
  }
  void createLeft(string t) {
    left = new TopCountNode(t);
  } 
  void createRight(string t) {
    right = new TopCountNode(t);
  }
  void pass(string t,TopCountNode* d) {
    d->check(t);
  }
  void check(string t) {
    if(t == top)
      count++;
    else if(t < top) {
      if(left==NULL)
	createLeft(t);
      else
	pass(t,left);
    }
    else { // t > top
      if(right==NULL)
	createRight(t);
      else
	pass(t,right);
    }
  }
  void printSubtree(ostream& f) {
    if(left!=NULL)
      left->printSubtree(f);
    if(right!=NULL)
      right->printSubtree(f);
    print(f);
  }
  void print(ostream& f) {
    f << top << " " << count << endl;
  }
  void collect(vector<TopCountNode*>& cnodes) {
    if(left != NULL)
      left->collect(cnodes);
    if(right != NULL)
      right->collect(cnodes);
    cnodes.push_back(this);
  }
  int getCount() { return count; }
  string getTop() { return top; }
private:
  TopCountNode* left;
  TopCountNode* right;
  TopCountNode* parent;
  int count;
  string top;
};

bool cmpTCNodes(TopCountNode* x,TopCountNode* y) {
  if(x->getCount() != y->getCount())
    return( x->getCount() > y->getCount() );
  else
    return( x->getTop() < y-> getTop() );
}

int main(int argc, char *argv[])
{
  if(argc==1)
    usageError(argv[0]);

  int numSkip = 0;
  string outFile;
  vector<string> fileNames;

  for(int i=1;i<argc;) {
    string a = (string)(argv[i]);
    if(a == "-h" || a == "--help")
      usageError(argv[0]);
    if(a == "--version") {
      cout << "mbsum version " << VERSION << endl;
      exit(1);
    }
    if(a == "-n" || a == "--skip") {
      if(i == argc-1)
	usageError(argv[0]);
      istringstream f((string)(argv[++i]));
      if(f.fail())
	usageError(argv[0]);
      f >> numSkip;
      if(f.fail())
	usageError(argv[0]);
      if(numSkip < 0) {
	cerr << "Error: Cannot skip a negative number of trees." << endl;
	usageError(argv[0]);
      }
    }
    else if(a == "-o" || a == "--out") {
      if(i == argc-1)
	usageError(argv[0]);
      istringstream f((string)(argv[++i]));
      if(f.fail())
	usageError(argv[0]);
      outFile.clear();
      f >> outFile;
      if(f.fail())
	usageError(argv[0]);
    }
    else
      fileNames.push_back(a);
    i++;
  }

  if(outFile.size() == 0) {
    if(fileNames.size() > 0) {
      string fileRoot = (string)(fileNames[0]);
      int ext = fileRoot.rfind('.');
      if(ext!=string::npos) // '.' in filename, erase from last '.'
	fileRoot.erase(ext,fileRoot.length() - ext);
      // erase directory
      ext = fileRoot.rfind('/');
      if(ext!=string::npos)
	fileRoot.erase(0,ext+1);
      outFile = fileRoot + ".in";
    }
    else
      outFile = "bucky.in";
  }

  ofstream sumOut(outFile.c_str());
  if(sumOut.fail()){
    cerr << "Error: Could not open file " << outFile << " ." << endl;
    exit(1);
  }

  int numFiles = fileNames.size();
  vector<int> numTrees(numFiles,0);
  vector<int> numTaxa(numFiles,0);
  TopCountNode* root;
  for(int i=0;i<numFiles;i++) {
    ifstream f(fileNames[i].c_str());
    if(f.fail())
      cerr << "Error: Could not open file " << fileNames[i] << " ." << endl;
    else
      cout << "Reading file " << fileNames[i] << ": " << flush;



    string line;
    int lineNumber=0;
    int numActuallySkipped=0;
    int skipInFile = numSkip;
    bool readTtable = false;
    while(getline(f,line)) {
      lineNumber++;
      istringstream s(line);
      string keyTree,name,equalSign;

      if (readTtable){
	char ch;
	s >> noskipws >> ch;
	while (ch != ';' && ch != EOF && s.good()) {
	  sumOut << ch;
	  s >> noskipws >> ch;
	}
	if (ch == ';'){
	  sumOut << ";\n";
	  readTtable = false;
	} else {
	  sumOut << "\n";
	  continue;
	}
      }

      s >> keyTree;
      if (i==0){
	// read translate table only in the first file, assuming the same in other files.
	if (keyTree=="translate" || keyTree=="TRANSLATE" || keyTree=="Translate"){
	  readTtable = true;
	  sumOut << "translate";
	  char ch;
	  s >> noskipws >> ch;
	  while (ch != ';' && ch != EOF && s.good()) {
	    sumOut << ch;
	    s >> noskipws >> ch;
	  }
	  if (ch == ';'){
	    sumOut << ";\n";
	    readTtable = false;
	  } else {
	    sumOut << "\n";
	    continue;
	  }
	}
      }

      // skip if the remaining line is not in format "  tree name = treeRep"
      if(keyTree != "tree")
	continue;

      s >> name >> equalSign;
      if(equalSign != "=")
	continue;
      // The rest should be a parenthetic representation of a tree.  Assume no spaces!
      // Skip the first numSkip trees
      if(skipInFile-- > 0) {
	numActuallySkipped++;
	continue;
      }

      string treeString;
      s >> treeString;
      Tree tree(treeString,lineNumber);
      numTrees[i]++;

      ostringstream g;
      tree.printTop(g); // changed from mb2badger to not include an endline
      string top = g.str();

      if((i==0) && numTrees[0] == 1) // first tree
	root = new TopCountNode(top);
      else
	root->check(top);
    }
    cout << "Skipped " << numActuallySkipped << " trees, read " << numTrees[i] << " trees." << endl;
  }

  vector<TopCountNode*> cnodes;
  root->collect(cnodes);
  sort(cnodes.begin(),cnodes.end(),cmpTCNodes);

  for(vector<TopCountNode*>::iterator n=cnodes.begin();n!=cnodes.end();n++)
    (*n)->print(sumOut);

  int total = 0;
  for(int i=0;i<numFiles;i++)
    total += numTrees[i];
  cout << "Read " << total << " total trees of which " << cnodes.size() << " are distinct." << endl;
  cout << "Output written to file " << outFile << endl;
  return(0);
}
