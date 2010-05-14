#include<string>
#include<vector>
#include<iostream>
#include <iomanip>
#include<cstdlib>
#include<sstream>
#include<cassert>
#include<algorithm>
#include<functional>
#include<set>
using namespace std;
#include "boost/unordered_map.hpp"
using namespace boost;
#include "TGM.h"
#include "Quartets.h"
using namespace quartet;

vector<vector<int> > nCr;
void precomputeNcR(int numTaxa) {
    // compute nC1, nC2, nC3, nC4 for n in [0, numTaxa - 1], this is used for ranking quartets
    int numSuperNodes = 2 * numTaxa - 3;
    for (int i = 0; i < numSuperNodes; i++) {
        double val = 1.0;
        vector<int> iCj;
        int k = 1;
        for (int j = 1; j <= 4; j++) {
            //compute iCj and store it in vector nCr
            for (; k <= j; k++) {
                val *= (i - k + 1) / (double) k;
            }
            iCj.push_back((int) val);
        }
        nCr.push_back(iCj);
    }
}

int computeRIndex(int t1, int t2, int t3, int t4) {
    //compute subset rank. assume sorted order is t1, t2, t3, t4 -
    // rank for a quartet can be obtained by computing t1 Choose 1 + t2 Choose 2 + t3 Choose 3 + t4 Choose 4
    return nCr[t1 - 1][0] + nCr[t2 - 1][1] + nCr[t3 - 1][2] + nCr[t4 - 1][3];
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

string TreeBuilder::getTree(Table* newTable, int numTaxa) {
    vector<string> topologies = newTable->getTopologies();
    precomputeNcR(numTaxa);
    int numNodes = numTaxa + numTaxa - 3;
    int numQuartets = numNodes * (numNodes - 1) * (numNodes - 2) * (numNodes -3) / 24;
    vector<vector<double> > counts;
    counts.resize(numQuartets);//(numTaxa -3) new super nodes are created
    for (int i = 0; i < counts.size(); i++) {
        counts[i].resize(3, 0.0);
    }

    for (int i = 0; i < topologies.size(); i++) {
        double totalCount = newTable->getTotalCounts(i);
        Tree t(numTaxa, topologies[i]);
        vector<int> rind, cind;
        t.getQuartets(rind, cind);
        for (int i = 0; i < rind.size(); i++) {
            counts[rind[i]][cind[i]] += totalCount;
        }
    }

    // normalize counts.
    int numOfCurrentQuartets = nCr[numTaxa - 1][2] + nCr[numTaxa - 1][3]; //C(numTaxa,4)  = C(numTaxa - 1, 3) + C(numTaxa - 1, 4)
    // current version: resolution with biggest count gets probability 1, other two get 0
    // TODO: if three counts are approximately equal, may need to split into 0.3333 each
    for (int i = 0; i < numOfCurrentQuartets; i++) {
//        if (counts[i][0] == counts[i][1] && counts[i][1] == counts[i][2]) {
//            counts[i][0] = .333333;
//            counts[i][1] = .333333;
//            counts[i][2] = .333334;
//        }
        if (counts[i][0] > counts[i][1]) {
            if (counts[i][0] > counts[i][2]) {
                counts[i][0] = 1;
                counts[i][1] = 0;
                counts[i][2] = 0;
            }
            else {
                counts[i][0] = 0;
                counts[i][1] = 0;
                counts[i][2] = 1;
            }
        }
        else if (counts[i][1] > counts[i][2]) {
            counts[i][0] = 0;
            counts[i][1] = 1;
            counts[i][2] = 0;
        }
        else {
            counts[i][0] = 0;
            counts[i][1] = 0;
            counts[i][2] = 1;
        }
    }

    return getTreeFromQuartetCounts(counts, numTaxa);
}

string TreeBuilder::getTreeFromQuartetCounts(vector<vector<double> >& counts, int numTaxa) {
    int numNodes = numTaxa + numTaxa - 3;
    for (int i = 0; i < numNodes; i++) {
        superNodes.push_back(new SuperNode(i + 1, i < numTaxa));
    }

    vector<int> activeNodes;
    for (int i = 0; i < numTaxa; i++) {
        activeNodes.push_back(i + 1);
    }

    vector<vector<double> > confidence(numNodes, vector<double>(numNodes, 0));
    vector<vector<double> > size(numNodes, vector<double>(numNodes, 0));
    vector<vector<double> > support(numNodes, vector<double>(numNodes, 0));
    for (int j = 1; j <= numTaxa; j++) {
        for (int i = 1; i < j; i++) {
            confidence[i - 1][j -1] = computeConfidence(i, j, activeNodes, counts);
            size[i - 1][j - 1] = computeCardinality(i, j, activeNodes);
            support[i -1][j - 1] = confidence[i - 1][j -1]/size[i - 1][j - 1];
        }
    }

    int currentNode = numTaxa;
    while (true) {
        int maxI = activeNodes[0], maxJ = activeNodes[1];
        double maxSupport = support[maxI - 1][maxJ - 1];
        for (int j = 0; j < activeNodes.size(); j++) {
            int node1 = activeNodes[j];
            for (int i = j + 1; i < activeNodes.size(); i++) {
                int node2 = activeNodes[i];
                if (support[node1 - 1][node2 - 1] > maxSupport) {
                    maxSupport = support[node1 - 1][node2 - 1];
                    maxI = node1;
                    maxJ = node2;
                }
            }
        }

        currentNode++;
        superNodes[currentNode - 1]->add(superNodes[maxI - 1], superNodes[maxJ - 1]);
        vector<int>::iterator itr = find(activeNodes.begin(), activeNodes.end(), maxI);
        assert(itr != activeNodes.end());
        activeNodes.erase(itr);
        itr = find(activeNodes.begin(), activeNodes.end(), maxJ);
        assert(itr != activeNodes.end());
        activeNodes.erase(itr);

        if (activeNodes.size() == 2) {
            activeNodes.push_back(currentNode); //quit when there are 3 active nodes(super nodes)
            break;
        }

        for (int i = 0; i < activeNodes.size(); i++) {
            int node1 = activeNodes[i];
            for (int j = i + 1; j < activeNodes.size(); j++) {
                int node2 = activeNodes[j];
                for (int k = j + 1; k < activeNodes.size(); k++) {
                    int node3 = activeNodes[k];
                    //find weights of 3 different resolutions for node1, node2, node3 and currentNode
                    int rInd1, cInd1, rInd2, cInd2;
                    int rIndResult, cIndResult;

                    getQuartetRowColumnIndex(node1, node2, node3, maxI, rInd1, cInd1);
                    getQuartetRowColumnIndex(node1, node2, node3, maxJ, rInd2, cInd2);
                    getQuartetRowColumnIndex(node1, node2, node3, currentNode, rIndResult, cIndResult);
                    counts[rIndResult][cIndResult] = counts[rInd1][cInd1] + counts[rInd2][cInd2];

                    getQuartetRowColumnIndex(node1, node3, node2, maxI, rInd1, cInd1);
                    getQuartetRowColumnIndex(node1, node3, node2, maxJ, rInd2, cInd2);
                    getQuartetRowColumnIndex(node1, node3, node2, currentNode, rIndResult, cIndResult);
                    counts[rIndResult][cIndResult] = counts[rInd1][cInd1] + counts[rInd2][cInd2];

                    getQuartetRowColumnIndex(node3, node2, node1, maxI, rInd1, cInd1);
                    getQuartetRowColumnIndex(node3, node2, node1, maxJ, rInd2, cInd2);
                    getQuartetRowColumnIndex(node3, node2, node1, currentNode, rIndResult, cIndResult);
                    counts[rIndResult][cIndResult] = counts[rInd1][cInd1] + counts[rInd2][cInd2];
                }
            }

        }

        //compute confidence, support for newnode and every other active node
        for (int i = 0; i < activeNodes.size(); i++) {
            int node1 = activeNodes[i];
            confidence[node1 - 1][currentNode - 1] = computeConfidence(node1, currentNode, activeNodes, counts);
            size[node1 - 1][currentNode - 1] = computeCardinality(node1, currentNode, activeNodes);
            support[node1 - 1][currentNode - 1] = confidence[node1 - 1][currentNode - 1] /size[node1 - 1][currentNode - 1];
        }

        // modify confidence, support of active node pairs from previous iteration
        for (int i = 0; i < activeNodes.size(); i++) {
            int node1 = activeNodes[i];
            for (int j = i + 1; j < activeNodes.size(); j++) {
                int node2 = activeNodes[j];
                int rind, cind;
                getQuartetRowColumnIndex(maxI, maxJ, node1, node2, rind, cind);
                confidence[node1 - 1][node2 - 1] -=  counts[rind][cind];
                size[node1 - 1][node2 - 1] -= superNodes[maxI - 1]->getNumNodes()
                        * superNodes[maxJ - 1]->getNumNodes()
                        * superNodes[node1 - 1]->getNumNodes()
                        * superNodes[node2-1]->getNumNodes();

                support[node1 -1][node2 - 1] = confidence[node1 - 1][node2 - 1] / size [node1 - 1][node2 - 1];
            }
        }
        activeNodes.push_back(currentNode);
    }

    SuperNode *node1, *node2, *node3;
    if (superNodes[activeNodes[2] - 1]->getLowestTaxon() < superNodes[activeNodes[1] - 1]->getLowestTaxon()) {
        if (superNodes[activeNodes[2] - 1]->getLowestTaxon() < superNodes[activeNodes[0] - 1]->getLowestTaxon()) {
            node1 = superNodes[activeNodes[2] - 1];
            if (superNodes[activeNodes[0] - 1]->getLowestTaxon() < superNodes[activeNodes[1] - 1]->getLowestTaxon()) {
                node2 =superNodes[activeNodes[0] - 1];
                node3 = superNodes[activeNodes[1] - 1];
            }
            else {
                node2 =superNodes[activeNodes[1] - 1];
                node3 = superNodes[activeNodes[0] - 1];
            }
        }
        else {
            node1 = superNodes[activeNodes[0] - 1];
            node2 =superNodes[activeNodes[2] - 1];
            node3 = superNodes[activeNodes[1] - 1];
        }
    }
    else {
        if (superNodes[activeNodes[1] - 1]->getLowestTaxon() < superNodes[activeNodes[0] - 1]->getLowestTaxon()) {
            node1 = superNodes[activeNodes[1] - 1];
            if (superNodes[activeNodes[2] - 1]->getLowestTaxon() < superNodes[activeNodes[0] - 1]->getLowestTaxon()) {
                node2 =superNodes[activeNodes[2] - 1];
                node3 = superNodes[activeNodes[0] - 1];
            }
            else {
                node2 =superNodes[activeNodes[0] - 1];
                node3 = superNodes[activeNodes[2] - 1];
            }
        }
        else {
            node1 = superNodes[activeNodes[0] - 1];
            node2 = superNodes[activeNodes[1] - 1];
            node3 = superNodes[activeNodes[2] - 1];
        }
    }

    //TODO: choose lowest taxon as outgroup
    stringstream top;
    top << "(";
    node1->print(top);
    top << ",(";
    node2->print(top);
    top << ",";
    node3->print(top);
    top << "));";
    return top.str();
}

double TreeBuilder::computeConfidence(int m, int n, vector<int>& activeNodes, vector<vector<double> >& counts) {
    double confidence = 0.0;
    for (int i = 0; i < activeNodes.size(); i++) {
        int node1 = activeNodes[i];
        if (node1 == m || node1 == n)
            continue;
        for (int j = i + 1; j < activeNodes.size(); j++) {
            int node2 = activeNodes[j];
            if (node2 == m || node2 == n)
                continue;

            int rIndex, cIndex;
            getQuartetRowColumnIndex(m, n, node1, node2, rIndex, cIndex);
            confidence += counts[rIndex][cIndex];
        }
    }

    return confidence;
}

int TreeBuilder::computeCardinality(int m, int n, vector<int>& activeNodes) {
    int cardinality = 0;
    for (int i = 0; i < activeNodes.size(); i++) {
        int node1 = activeNodes[i];
        if (node1 == m || node1 == n)
            continue;
        for (int j = i + 1; j < activeNodes.size(); j++) {
            int node2 = activeNodes[j];
            if (node2 == m || node2 == n)
                continue;

            cardinality += superNodes[node1 - 1]->getNumNodes() * superNodes[node2 - 1]->getNumNodes();
        }
    }

    cardinality *= superNodes[m - 1]->getNumNodes() * superNodes[n - 1]->getNumNodes();
    return cardinality;
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
    int i, j;
    int matchI = 0, matchJ = 0;
    int greatestMatch = 0;
    for (i = 0; i < 3; i++) {
        vector<int> t1 = n1->getTaxa(i);
        for (j = 0; j < 3; j++) {
            vector<int> t2 = n2->getTaxa(j);
            int match = intersect(t1, t2);
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
                        }
                    }
                }
            }
        }
    }
}

/*int main(int argc, char *argv[]) {
    string top = "(1,(2,(3,(4,(5,6)))));";
//      string top = "(((1,4),2),(3,(5,6)));";
//    string top = "((1,2),((3,4),(5,6)));";
//    string top = "(((((1,2),3),4),5),6);";
//        string top = "(1,((2,3),4));";


      vector<string> topologies;
      topologies.push_back(top);
      Table *newTable = new TGM(topologies);
      newTable->addGeneCount(0,0,10);
      newTable->addGeneCount(0,1,10);
      newTable->addGeneCount(0,2,10);
      newTable->addGeneCount(0,3,10);

      TreeBuilder t;
      top = t.getTree(newTable, 6);
      cerr << "Final topology: " << top<< endl;
}*/
