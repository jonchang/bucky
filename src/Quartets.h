#ifndef QUARTETS_H_
#define QUARTETS_H_

class Quartets {
public:
    Quartets(string top);
    virtual ~Quartets();
private:
    vector<int> qs;
};

class TaxonSet {
public:
    TaxonSet() {
    }

    void add(int taxon) {
        tset.push_back(taxon);
    }

    void merge(TaxonSet*& other) {
        for (vector<int>::iterator itr = other->tset.begin(); itr != other->tset.end(); itr++) {
            tset.push_back(*itr);
        }
    }

    vector<int>& getTSet() {
        return tset;
    }

    friend ostream& operator<<(ostream& f, TaxonSet& t) {
        for (vector<int>::iterator itr = t.tset.begin(); itr != t.tset.end(); itr++) {
            f << *itr << ",";
        }
        return f;
    }
private:
    vector<int> tset;
};

class Edge;

class Node {
public:
  Node(int n,int lf) : number(n) {
    t = new TaxonSet();
    if(lf)
      leaf = true;
    else
      leaf = false;
  }
  ~Node() { delete t;}
  int getNumber() const { return number; }
  Edge* getEdge(int n) const { return edges[n]; }
  bool isLeaf() const { return leaf; }
  void setEdge(int n,Edge *e) { edges[n]=e; }
  void print(ostream&) const;
  Node* getNeighbor(int) const;
  void mergeTaxa(Node *n) {
      t->merge(n->t);
  }
  void addTaxa(int taxon) {
      t->add(taxon);
  }
  vector<int>& getTaxa() {
      return t->getTSet();
  }
private:
  int number;
  bool leaf;
  Edge* edges[3];
  TaxonSet* t;
};

class Edge {
public:
  Edge(int n) : number(n) {}
  int getNumber() const { return number; }
  Node* getNode(int i) { return nodes[i]; }
  Node* getOtherNode(const Node *n) const { return (n==nodes[0] ? nodes[1] : nodes[0]); }
  void setNode(int i,Node *n) { nodes[i] = n; }
  void print(ostream&,int) const;
  ~Edge() {}
private:
  int number;
  Node* nodes[2];
};

class Tree {
public:
  Tree(int nTaxa,string top) {
    numTaxa = nTaxa;
    numNodes = 2*numTaxa - 2;
    numEdges = 2 * numTaxa - 3;
    numSplits = numTaxa - 3;
    nodes.resize(numNodes);
    edges.resize(numEdges);
    for(int i=0;i<numTaxa;i++)
      nodes[i] = new Node(i,1);
    for(int i=numTaxa;i<numNodes;i++)
      nodes[i] = new Node(i,0);
    for(int i=0;i<numEdges;i++)
      edges[i] = new Edge(i);
    connect(top);
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
  void getQuartets(vector<vector<int> >& quartets);
  int getNumTaxa() const { return numTaxa; }
  int getNumNodes() const { return numNodes; }
  int getNumEdges() const { return numEdges; }
  int getNumSplits() const { return numSplits; }
  Node* getNode(int i) const { return nodes[i]; }
  Edge* getEdge(int i) const { return edges[i]; }
  string getTop() const { return top; }
private:
  int numTaxa,numNodes,numEdges,numSplits;
  vector<Node *> nodes;
  vector<Edge *> edges;
  string top;
};

#endif /* QUARTETS_H_ */
