#include<string>
#include<vector>
#include<iostream>
#include <iomanip>
#include<cstdlib>
#include<sstream>
#include<cassert>
#include<algorithm>
#include <functional>
using namespace std;
#include "Quartets.h"
using namespace quartet;

//compute nCr
int choose(int n, int r) {
    if (n < r) return 0;
    if (r > n - r) // save computations by evaluating nCn-r if n-r is smaller than r
        r = n - r;

    double val = 1.0;
    for (int i = 1; i <= r; i++) {
        val *= (n - i + 1) / (double) i;
    }

    return (int) val;
}

int computeRIndex(int t1, int t2, int t3, int t4) {
    //compute subset rank. assume sorted order is t1, t2, t3, t4
    return choose(t1 - 1, 1) + choose(t2 - 1, 2) + choose(t3 - 1, 3) + choose(t4 - 1, 4);
}

void getQuartetRowColumnIndex(int t1, int t2, int t3, int t4, int& rIndex, int& cIndex) {
    // make t3 smallest of t3, t4
    if (t3 > t4) {
        int temp = t4;
        t4 = t3;
        t3 = temp;
    }

    // make t1 smallest of t1, t2
    if (t1 > t2) {
        int temp = t1;
        t1 = t2;
        t2 = temp;
    }

    if (t1 > t3) {
        // swap t1 with t3, t2 with t4 - makes t1 smallest
        int temp = t1;
        t1 = t3;
        t3 =temp
                ;
        temp = t2;
        t2 = t4;
        t4 = temp;
    }

    if (t2 < t3) {
        //t2 2nd smallest - order t1, t2, t3, t4
        assert (t1 < t2);
        assert (t2< t3);
        assert (t3< t4);
        cIndex = 0;
        rIndex = computeRIndex(t1, t2, t3, t4);
    }
    else {
        //t3 2nd smallest
        if (t2 < t4) {
            //  order t1, t3, t2, t4
            cIndex = 1;
            rIndex = computeRIndex(t1, t3, t2, t4);
        }
        else {
            // order t1, t3, t4, t2
            cIndex = 2;
            rIndex = computeRIndex(t1, t3, t4, t2);
        }
    }
}

Node* Node::getNeighbor(int i) const { return edges[i]->getOtherNode(this); }

void Node::print(ostream& f) const
{
    f << setw(3) << number << ":" << *t[0] << ":" << *t[1] << ":" << *t[2] << ":";
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
    f << endl;
}

void Tree::print(ostream& f) const
{
    f << top << endl;
    f << "numTaxa = " << numTaxa << ", numNodes = " << numNodes << ", numEdges = " << numEdges << endl;
    f << "Nodes:" << endl;
    for(int i=0;i<numNodes;i++)
        nodes[i]->print(f);
    f << "Edges:" << endl;
    for(int i=0;i<numEdges;i++)
        edges[i]->print(f,numTaxa);
}

void Tree::connect(string top)
{
    int currentLeafNode = 0, currentInternalNode=numTaxa,currentLeafEdge=0,currentInternalEdge=numTaxa;
    Node *node1,*node2;
    int n;

    istringstream topStr(top);
    char ch = topStr.peek();
    if(ch!='(') {
        topStr >> n;
    }
    else {
        topStr >> ch;
        node1 = connectInt(topStr,currentLeafNode,currentInternalNode,currentLeafEdge,currentInternalEdge);
        topStr >> ch;
        if(ch!=',') {
            cerr << "Error: Cannot parse topology string, expected ',', found " << ch << endl;
            exit(1);
        }
        node2 = connectInt(topStr,currentLeafNode,currentInternalNode,currentLeafEdge,currentInternalEdge);
        topStr >> ch;
        if(ch!=')') {
            cerr << "Error: Cannot parse topology string, expected ')', found " << ch << endl;
            exit(1);
        }
        connectTwoNodes(node1,node2,currentLeafEdge,currentInternalEdge);
    }

    for (int i = numTaxa; i < numNodes; i++) {
        if (nodes[i]->getTaxa(0).size() == numTaxa) {
            root = nodes[i];
            break;
        }
        else if (nodes[i]->getTaxa(1).size() == numTaxa) {
            root = nodes[i];
            break;
        }
        else if (nodes[i]->getTaxa(2).size() == numTaxa) {
            root = nodes[i];
            break;
        }
    }
//    cerr << root->getNumber() << endl;
    TaxonSet* t0 = root->getTset(0);
    TaxonSet* t1 = root->getTset(1);
    TaxonSet* t2 = root->getTset(2);
    if (t0->getTSet().size() < t1->getTSet().size()) {
        TaxonSet *temp = t0;
        t0 = t1;
        t1 = t0;
    }
    if (t0->getTSet().size() < t2->getTSet().size()) {
        TaxonSet *temp = t0;
        t0 = t2;
        t2 = t0;
    }
    t0->remove(t1);
    t0->remove(t2);
}

Node* Tree::connectInt(istream& topStr,int& currentLeafNode,int& currentInternalNode,int& currentLeafEdge,int& currentInternalEdge)
{
    int n;
    char ch = topStr.peek();
    if(ch!='(') {
        topStr >> n;
        nodes[currentLeafNode]->addTaxa(0, n);
        return nodes[currentLeafNode++];
    }
    else {
        topStr >> ch;
        Node *node1,*node2;
        node1 = connectInt(topStr,currentLeafNode,currentInternalNode,currentLeafEdge,currentInternalEdge);
        topStr >> ch;
        if(ch!=',') {
            cerr << "Error: Cannot parse topology string, expected ',', found " << ch << endl;
            exit(1);
        }
        node2 = connectInt(topStr,currentLeafNode,currentInternalNode,currentLeafEdge,currentInternalEdge);
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

    node3->mergeTaxa(1, node1, 0);
    node3->mergeTaxa(2, node2, 0);
    node3->mergeTaxa(0, node1, 0);
    node3->mergeTaxa(0, node2, 0);
//    cerr << "Connect " << node1->getNumber() << " " << node2->getNumber() << " " << node3->getNumber() << endl;
    return node3;
}

void Tree::connectTwoNodes(Node* node1, Node* node2, int& currentLeafEdge, int& currentInternalEdge)
{
    Edge *edge;
    if(node1->isLeaf() || node2->isLeaf())
        edge = edges[currentLeafEdge++];
    else
        edge = edges[currentInternalEdge++];

    if (node1->isLeaf() && !node2->isLeaf()) {
        node2->addTaxa(0, node1->getTaxa(0)[0]);
    }
    else {
        node1->mergeTaxa(0, node2, 0);
    }

    node1->setEdge(0,edge);
    node2->setEdge(0,edge);
    edge->setNode(0,node1);
    edge->setNode(1,node2);
//    cerr << "Connect " << node1->getNumber() << " " << node2->getNumber() << endl;
}

int Tree::intersect(vector<int>& t1, vector<int>& t2) {
    vector<int> v(t1.size()> t2.size() ? t2.size() : t1.size());
//    sort(t1.begin(), t1.end());
//    sort(t2.begin(), t2.end());
 /*   for (int i = 0; i < t1.size(); i++) {
        cerr << t1[i] << " ";
    }
    cerr << endl;
    for (int i = 0; i < t2.size(); i++) {
        cerr << t2[i] << " ";
    }
    cerr << endl;
    cerr << "Intersection " << endl;
    for (int i = 0; i < v.size(); i++) {
        cerr << v[i] << " ";
    }
    cerr << endl;*/
    vector<int>::iterator it = set_intersection(t1.begin(), t1.end(), t2.begin(), t2.end(), v.begin());
    if (it -v.begin() > 0 && t1.size() == t2.size() && t2.size() == v.size()) {
        return numeric_limits<int>::max();
    }
    return int(it - v.begin());
}

void Tree::getTsets(vector<vector <int> >& tsets, Node *n1, Node *n2) {
//    cerr << " NOdes " << n1->getNumber() << " " << n2->getNumber() << endl;
    int i, j;
    int matchI = 0, matchJ = 0;
    int greatestMatch = 0;
    for (i = 0; i < 3; i++) {
        vector<int> t1 = n1->getTaxa(i);
        for (j = 0; j < 3; j++) {
            vector<int> t2 = n2->getTaxa(j);
            int match = intersect(t1, t2);
//            cerr << " intersection size from " << i << ", " << j << " is " << match << " " <<greatestMatch << endl;
            if (match != 0 && match >= greatestMatch) {
                greatestMatch = match;
                matchI = i;
                matchJ= j;
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        if (i != matchI) {
            vector<int> tset = n1->getTaxa(i);
            if (n1 != root && i == 0) {
                vector<int> compliment;
                for (int i = 1; i <= numTaxa; i++) {
                    vector<int>::iterator it = find(tset.begin(), tset.end(), i);
                    if (it == tset.end()) {
                        compliment.push_back(i);
                    }
                }
                tsets.push_back(compliment);
            }
            else {
                tsets.push_back(tset);
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        if (i != matchJ) {
            vector<int> tset = n2->getTaxa(i);
            if (n2 != root && i == 0) {
                vector<int> compliment;
                for (int i = 1; i <= numTaxa; i++) {
                    vector<int>::iterator it = find(tset.begin(), tset.end(), i);
                    if (it == tset.end()) {
                        compliment.push_back(i);
                    }
                }
                tsets.push_back(compliment);
            }
            else {
                tsets.push_back(tset);
            }
        }
    }
//    cerr << matchI << " " << matchJ << " should be skipped for " << n1->getNumber() << " " << n2->getNumber() << endl;
}

void Tree::getQuartets(vector<int>& rInd, vector<int>& cInd) {
    vector<vector<int> > quartets;
    for (int i = numTaxa; i < numNodes; i++) {
        Node *n1 = nodes[i];
        for (int j = i + 1; j < numNodes; j++) {
            Node *n2 = nodes[j];
            vector<vector<int> > tsets;
            getTsets(tsets, n1, n2);
            assert(tsets.size() == 4);
            for (int j1 = 0; j1 < tsets[0].size(); j1++) {
                vector<int> quartet;
                quartet.push_back(tsets[0][j1]);
                for (int j2 = 0; j2 < tsets[1].size(); j2++) {
                    quartet.push_back(tsets[1][j2]);
                    for (int j3 = 0; j3 < tsets[2].size(); j3++) {
                        quartet.push_back(tsets[2][j3]);
                        for (int j4 = 0; j4 < tsets[3].size(); j4++) {
                            quartet.push_back(tsets[3][j4]);
                            int rIndex, cIndex;
                            getQuartetRowColumnIndex(tsets[0][j1], tsets[1][j2], tsets[2][j3], tsets[3][j4], rIndex, cIndex);
                            rInd.push_back(rIndex);
                            cInd.push_back(cIndex);
                            quartets.push_back(quartet);
                            cerr << tsets[0][j1] << "," << tsets[1][j2] << "|" << tsets[2][j3] << "," << tsets[3][j4] << " Row:" << rIndex << ", Column:" << cIndex << "\n";
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
//    string top = "(1,(2,(3,(4,(5,6)))));";
      string top = "(((1,4),2),(3,(5,6)));";
//    string top = "((1,2),((3,4),(5,6)));";
//    string top = "(((((1,2),3),4),5),6);";
//        string top = "(1,((2,3),4));";
    Tree t(6, top);
    t.print(cerr);
    vector<vector<int> > quartets;
    vector<int> rInd, cInd;
    t.getQuartets(rInd,cInd);
//    computeRIndex(1, 2, 3, 6);
}
