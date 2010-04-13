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

Node* Node::getNeighbor(int i) const { return edges[i]->getOtherNode(this); }

void Node::print(ostream& f) const
{
    f << setw(3) << number << ":" << *t << ":";
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
    int currentInternalNode=numTaxa,currentLeafEdge=0,currentInternalEdge=numTaxa;
    Node *node1,*node2;
    int n;

    istringstream topStr(top);
    char ch = topStr.peek();
    if(ch!='(') {
        topStr >> n;
    }
    else {
        topStr >> ch;
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
    int n;
    char ch = topStr.peek();
    if(ch!='(') {
        topStr >> n;
        nodes[n-1]->addTaxa(n);
        return nodes[n-1];
    }
    else {
        topStr >> ch;
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
    node3->mergeTaxa(node1);
    node3->mergeTaxa(node2);

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

    if (node1->isLeaf() && !node2->isLeaf()) {
        node2->mergeTaxa(node1);
    }
    else {
    node1->mergeTaxa(node2);
    }
    node1->setEdge(0,edge);
    node2->setEdge(0,edge);
    edge->setNode(0,node1);
    edge->setNode(1,node2);

}

void Tree::getQuartets(vector<vector<int> >& quartets) {
    Node *root;
    for (int i = numTaxa; i < numNodes; i++) {
        if (nodes[i]->getTaxa().size() == numTaxa) {
            root = nodes[i];
            break;
        }
    }
    //for each internal edge i, pick 4 taxa from 4 different sides to create a quartet
    for(int i=numTaxa;i<numEdges;i++) {
        vector<vector<int> > tsets; // stores taxonsets on 4 different sides of a internal edge.
        Node *n1 = edges[i]->getNode(0);
        Node *n2 = edges[i]->getNode(1);
//        cerr << n1->getNumber() << " nbr nodes: " <<endl;
        for (int j = 0; j < 3; j++) {
            Node *other = n1->getNeighbor(j);
//            cerr << other->getNumber();
            if (other != n2) {
                vector<int> taxa(other->getTaxa());
                if (n1 != root && j == 0) {//in case other node is parent, remove all taxa below it
                    vector<int> taxa(root->getTaxa());
                    vector<int> ntaxa = n1->getTaxa();
                    for (int i = 0; i < ntaxa.size(); i++) {
                        taxa.erase(find(taxa.begin(), taxa.end(), ntaxa[i]));
                    }
                    tsets.push_back(taxa);
                }
                else {
                    tsets.push_back(other->getTaxa());
                }
//                cerr << "<--Including:";
//                for (int ii = 0; ii < tsets[tsets.size() - 1].size(); ii++) {
//                    cerr << tsets[tsets.size() - 1][ii] <<",";
//                }
            }
//            cerr << endl;
        }

//        cerr << n2->getNumber() << " nbr nodes: " <<endl;
        for (int j = 0; j < 3; j++) {
            Node *other = n2->getNeighbor(j);
//            cerr << other->getNumber();
            if (other != n1) {
                if (n2 != root && j == 0) {//in case other node is parent, remove all taxa below it
                    vector<int> taxa(root->getTaxa());
                    vector<int> ntaxa = n2->getTaxa();
                    for (int i = 0; i < ntaxa.size(); i++) {
                        taxa.erase(find(taxa.begin(), taxa.end(), ntaxa[i]));
                    }
                    tsets.push_back(taxa);
                }
                else {
                    tsets.push_back(other->getTaxa());
                }
//                cerr << "<--Including";
//                for (int ii = 0; ii < tsets[tsets.size() - 1].size(); ii++) {
//                    cerr << tsets[tsets.size() - 1][ii] <<",";
//                }
            }
//            cerr << endl;
        }

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
                        cerr << tsets[0][j1] << "," << tsets[1][j2] << "," << tsets[2][j3] << "," << tsets[3][j4] << "\n";
                        quartets.push_back(quartet);
                    }
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
//    string top = "(1,(2,(3,(4,(5,6)))));";
      string top = "(((1,6),2),((3,4),5));";
//    string top = "((1,2),((3,6),(4,5)));";
    Tree t(6, top);
    t.print(cerr);
    vector<vector<int> > quartets;
    t.getQuartets(quartets);
}
