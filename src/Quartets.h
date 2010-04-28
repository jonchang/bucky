#ifndef QUARTETS_H_
#define QUARTETS_H_
namespace quartet {

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
        //merge two sorted ranges, result is also sorted.
        int initialSize = tset.size();
        tset.resize(initialSize + other->tset.size());
        vector<int> otherset = other->getTSet();
        copy(otherset.begin(), otherset.begin() + otherset.size(), tset.begin() + initialSize);
        inplace_merge(tset.begin(), tset.begin() + initialSize, tset.end());
    }

    vector<int>& getTSet() {
        return tset;
    }

    void remove(TaxonSet*& other) {
        for (vector<int>::iterator itr = other->tset.begin(); itr != other->tset.end(); itr++) {
            vector<int>::iterator pos = find(tset.begin(), tset.end(), *itr);
            if (pos != tset.end())
                tset.erase(pos);
        }
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
    for (int i = 0; i < 3; i++)
      t[i] = new TaxonSet();

    if(lf)
      leaf = true;
    else
      leaf = false;
  }
  ~Node() { /*delete[] t;*/}
  int getNumber() const { return number; }
  Edge* getEdge(int n) const { return edges[n]; }
  bool isLeaf() const { return leaf; }
  void setEdge(int n,Edge *e) { edges[n]=e; }
  void print(ostream&) const;
  Node* getNeighbor(int) const;
  void mergeTaxa(int index, Node *n, int nIndex) {
      t[index]->merge(n->t[nIndex]);
  }
  void addTaxa(int index, int taxon) {
      t[index]->add(taxon);
  }
  vector<int>& getTaxa(int index) {
      return t[index]->getTSet();
  }

  TaxonSet*& getTset(int index) {
      return t[index];
  }
private:
  int number;
  bool leaf;
  Edge* edges[3];
  TaxonSet* t[3];
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
  Node* connectInt(istream&,int&,int&,int&,int&);
  Node* connectThreeNodes(Node*,Node*,Node*,int&,int&);
  void connectTwoNodes(Node*,Node*,int&,int&);
  void print(ostream&) const;
  void setAllTaxa();
//  void getQuartets(vector<vector<int> >& quartets);
  void getQuartets(vector<int>& rInd, vector<int>& cInd);
  int getNumTaxa() const { return numTaxa; }
  int getNumNodes() const { return numNodes; }
  int getNumEdges() const { return numEdges; }
  int getNumSplits() const { return numSplits; }
  Node* getNode(int i) const { return nodes[i]; }
  Edge* getEdge(int i) const { return edges[i]; }
  string getTop() const { return top; }
  int intersect(vector<int>& t1, vector<int>& t2);
  void getTsets(vector<vector <int> >& tsets, Node *n1, Node *n2);
private:
  int numTaxa,numNodes,numEdges,numSplits;
  vector<Node *> nodes;
  vector<Edge *> edges;
  string top;
  Node *root;
};
}
#endif /* QUARTETS_H_ */
