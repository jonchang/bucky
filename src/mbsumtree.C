#include<iostream>
#include<cmath>
#include<vector>
#include<map>
#include<sstream>
#include"mbsumtree.h"

using namespace std;

namespace mbsumtree {

bool cmpNodes(Node* x,Node* y) { return( x->getNumber() < y->getNumber() ); }

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
  int nextNumber = maxTaxNumber + 1;
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

void advance(istringstream& s,char& c,int numLeft,int lineNumber) {
  c = s.peek();
  if(s.eof() || s.fail()) {
    cerr << "Error: Tree on line " << lineNumber
    << " ends prematurely while reading subtree after left parenthesis " << numLeft << "." << endl;
    exit(1);
  }
  s >> c;
}

bool Tree::readSubtree(istringstream& s,Node* parent,int& numLeft,int numLine)
{
  Node* n;
  Edge* e;
  char c = s.peek();
  if(c=='(') {
    n = new Node();
    s >> c;
    numLeft++;
    bool bl = readSubtree(s,n,numLeft,numLine);
    advance(s,c,numLeft,numLine);
    bool br = false;
    if(c == ',') {
      br = readSubtree(s,n,numLeft,numLine);
    }
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

    if (bl && br) {
      nodes.push_back(n);
      numNodes++;
      e = new Edge();
      edges.push_back(e);
      numEdges++;
      e->setNodes(n,parent);
      n->setEdge(e);
      parent->setEdge(e);
      n->setLeaf(false);
    }
    else if (!(bl || br)){
      // none of the childs have nodes.
      delete n;
      return false;
    }
    else {
      if (n->getNumEdges() == 1) {
          Edge *edge = n->getEdge(0);
          Node *child = edge->getOtherNode(n);
          edge->setNodes(child, parent);
          parent->setEdge(edge);
      }

      delete n;
      return true;
    }
  }
  else if(isdigit(c)) {
    int num;
    s >> num;
    if (!thePruner->prune(num)) {
      n = new Node();
      nodes.push_back(n);
      numNodes++;
      e = new Edge();
      edges.push_back(e);
      numEdges++;
      e->setNodes(n,parent);
      n->setEdge(e);
      parent->setEdge(e);
      n->setNumber(num);
      if (num > maxTaxNumber)
        maxTaxNumber = num;
      n->setLeaf(true);
      numTaxa++;
    }
    else {
      return false;
    }
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

  if (e == NULL) {
    return false;
  }
  e->setLength(length);
  return true;
}

Tree::Tree(string line,int lineNumber, Pruner* aPruner):thePruner(aPruner)
{
  // read in the tree from parenthetic representation.
  // Create new nodes and edges on the fly.
  // Then renumber.
  // Then, reorder so that leaves come first.

  numEdges = 0;
  numTaxa = 0;
  maxTaxNumber = 0;
  numNodes = 0;

  int numLeft = 0;
  istringstream s(line);
  char c = s.peek();
  if(c =='(') {
    Node *n = new Node(); // root of the tree
    numLeft++;
    // root is not a leaf
    n->setLeaf(false);
    s >> c; // read in left parenthesis
    bool b = false;
    while(true) {
      b |= readSubtree(s,n,numLeft,lineNumber);
      s >> c;
      if(s.eof() || s.fail()) {
        cerr << "Error: Tree on line " << lineNumber << " ended while reading in subtree beginning with left parenthesis number "
        << numLeft << "." << endl;
        exit(0);
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
    if (n->getNumEdges() == 1) { // if root only has one child, then remove root, make child new root
      Edge *edge = n->getEdge(0);
      Node *child = edge->getOtherNode(n);
      for (vector<Node *>::iterator itr = nodes.begin(); itr != nodes.end();) {
        if (*itr == child) {
          nodes.erase(itr);
        }
        else
          itr++;
      }
      nodes.insert(nodes.begin(), child);
      child->removeEdge(edge);
      numEdges--;
      delete n;
      for (vector<Edge *>::iterator itr = edges.begin(); itr != edges.end();) {
        if (*itr == edge) {
          edges.erase(itr);
        }
        else
          itr++;
      }
      delete edge;
    }
    else {
      nodes.insert(nodes.begin(), n);
      numNodes++;
    }
  }
  else if(isdigit(c)) { // base tree may just be a single node
    int n;
    s >> n;
    if (!thePruner->prune(n)) {
      nodes.push_back(new Node()); // root of the tree, thinks it is a leaf....
      numNodes = 1;
      nodes[0]->setLeaf(true);
      nodes[0]->setNumber(n);
      if (n>maxTaxNumber)
        maxTaxNumber = n;
      numTaxa++;
      if(c==';')
        s >> c;
      s >> c;
      if(!(s.eof())) {
        cerr << "Error: Tree on line  " << lineNumber << " begins with a number, but is not a one-taxon tree." << endl;
        exit(1);
      }
    }
  }

  if (numNodes == 0) //pruning can make a void tree
    return;

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
  if (numNodes == 0) {
    return;
  }

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
  f<< endl;
}
}
