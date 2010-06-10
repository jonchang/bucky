#define NDEBUG
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
            assert (t1 < t3);
            assert (t3 < t2);
            assert (t2< t4);
            cIndex = 1;
            rIndex = computeRIndex(t1, t3, t2, t4);
        }
        else {
            // order t1, t3, t4, t2
            assert (t1 < t3);
            assert (t3 < t4);
            assert (t4 < t2);
            cIndex = 2;
            rIndex = computeRIndex(t1, t3, t4, t2);
        }
    }
}

void TreeBuilder::getTree(Table* newTable, int numTaxa, string& top, string& topWithWts) {
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

    for (int i = 0; i < counts.size(); i++) {
        double total = counts[i][0] + counts[i][1] + counts[i][2];
        if (total != 0.0) {
            counts[i][0] /= total;
            counts[i][1] /= total;
            counts[i][2] /= total;
        }
    }
    top = getTreeFromQuartetCounts(counts, numTaxa);
    Tree t(numTaxa, top);
    for (int i = 0; i < numTaxa; i++) {
        t.getEdge(i)->addWeight(10.0);//set prob 1.0, branch length in #coalescent units as 10 for all leaf edges
    }

    //compute weights for internal edges using formula W(AB|CD)/|A||B||C||D|
    for (int i = numTaxa; i < t.getNumEdges(); i++) {
        Edge *e = t.getEdge(i);
        Node *n1 = e->getNode(0);
        Node *n2 = e->getNode(1);
        vector<int> rInd, cInd;
        int abcd = t.getQuartets(rInd, cInd, n1, n2); // returns |A||B||C||D|
        double wt = 0.0;
        for (int i = 0; i < rInd.size(); i++) {
            wt += counts[rInd[i]][cInd[i]];
        }

        wt /= (double) abcd;
        //convert prob to # of coalescent units.
        if (wt > 0.99997)
            wt = 10.0;
        else if (wt <= .33334)
            wt = 0.0;
        else
            wt = -log(3.0/2.0 * (1-wt));

        e->addWeight(wt);
    }

    t.modifyOutgroup(numTaxa, top, topWithWts);//choose last taxon as outgroup and print topologies to top, topWithWts.
}

// returns smallest number in range [1, numTaxa] that is not present in t1, t2
int getComplement(vector<int>& t1, vector<int>& t2, int numTaxa)
{
    vector<int>::iterator it1 = t1.begin();
    vector<int>::iterator it2 = t2.begin();
    for (int i = 1; i <= numTaxa; i++) {
        if (i == *it1) {
            it1++;
            continue;
        }

        if (i == *it2) {
            it2++;
            continue;
        }

        return i;
    }
}

string TreeBuilder::getTreeFromQuartetCounts(vector<vector<double> > counts, int numTaxa) {
    // normalize counts.
    int numOfCurrentQuartets = nCr[numTaxa - 1][2] + nCr[numTaxa - 1][3]; //C(numTaxa,4)  = C(numTaxa - 1, 3) + C(numTaxa - 1, 4)
    // current version: resolution with biggest count gets probability 1, other two get 0
    // TODO: if three counts are approximately equal, may need to split into 0.3333 each
    for (int i = 0; i < numOfCurrentQuartets; i++) {
//        if (counts[i][0] == counts[i][1] && counts[i][1] == counts[i][2]) {
//            counts[i][0] = .333334;
//            counts[i][1] = .333333;
//            counts[i][2] = .333333;
//        }
        if (counts[i][0] >= counts[i][1]) {
            if (counts[i][0] >= counts[i][2]) {
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
        else if (counts[i][1] >= counts[i][2]) {
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

    int numNodes = numTaxa + numTaxa - 3;
    for (int i = 0; i < numNodes; i++) {
        superNodes.push_back(new SuperNode(i + 1, i < numTaxa));
    }

    vector<int> activeNodes; // keeps track of list of super nodes active in current iteration
    for (int i = 0; i < numTaxa; i++) {
        activeNodes.push_back(i + 1);
    }

    vector<vector<double> > confidence(numNodes, vector<double>(numNodes, 0));
    vector<vector<double> > size(numNodes, vector<double>(numNodes, 0));
    vector<vector<double> > support(numNodes, vector<double>(numNodes, 0));
    for (int j = 1; j <= numTaxa; j++) {
        for (int i = 1; i < j; i++) {
            confidence[i - 1][j -1] = computeConfidence(i, j, activeNodes, counts);
            size[i - 1][j - 1] = nCr[numTaxa - 2][1]; //initial cardinality is always numTaxa-2 choose 2, no need to do summation
            support[i -1][j - 1] = confidence[i - 1][j -1]/size[i - 1][j - 1];
        }
    }

    int currentNode = numTaxa;
    while (true) {
        // find nodes that have maximum support, call them maxI, maxJ
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
                    //find weights of 3 different resolutions for node1, node2, node3 and newly created node(currentNode)
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

        //compute confidence, support for currentNode and every other active node
        for (int i = 0; i < activeNodes.size(); i++) {
            int node1 = activeNodes[i];
            confidence[node1 - 1][currentNode - 1] = computeNewConfidence(maxI, maxJ, node1,activeNodes,counts,confidence);
            size[node1 - 1][currentNode - 1] = computeNewCardinality(maxI, maxJ, node1, numTaxa, activeNodes, size);
            support[node1 - 1][currentNode - 1] = confidence[node1 - 1][currentNode - 1] / size[node1 - 1][currentNode - 1];
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

    SuperNode *node1 = superNodes[activeNodes[0] - 1];
    SuperNode *node2 = superNodes[activeNodes[1] - 1];
    SuperNode *node3 = superNodes[activeNodes[2] - 1];;

    stringstream top;
    top << "(";
    node1->print(top);
    top << ",";
    node2->print(top);
    top << ",";
    node3->print(top);
    top << ");";
    return top.str();
}

//compute initial confidence, can also be used for computing confidence in later steps but inefficient,  use computeNewConfidence instead
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

// compute confidence of pairs that contain new super node, formula(D!=K,B): C(K, B) = C(I,B) + C(J,B) - Sum(counts[I][B][J][D] + counts[J][B][I][D])
double TreeBuilder::computeNewConfidence(int i, int j, int b, vector<int>& activeNodes, vector<vector<double> >& counts, vector<vector<double> >& confidence) {
    double conf = 0.0;
    if (i < b) {
        conf += confidence[i - 1][b - 1];
    }
    else {
        conf += confidence[b - 1][i - 1];
    }

    if (j < b) {
        conf += confidence[j - 1][b - 1];
    }
    else {
        conf += confidence[b - 1][j - 1];
    }

    for (int ii = 0; ii < activeNodes.size(); ii++) {
        int d = activeNodes[ii];
        if (d == b)
            continue;

        int rind, cind;
        getQuartetRowColumnIndex(i, b, j, d, rind, cind);
        conf -= counts[rind][cind];

        getQuartetRowColumnIndex(j, b, i, d, rind, cind);
        conf -= counts[rind][cind];
    }

    return conf;
}

// compute cardinality of pairs that contain new super node, formula T(K, B) = T(I,B) + T(J,B) - 2 * size(I) * size(J) * size(B) * (ntaxa - size(B) - size(K))
double TreeBuilder::computeNewCardinality(int maxI, int maxJ, int node1, int numTaxa, vector<int>& activeNodes, vector<vector<double> >& size) {
    double sz = 0.0;
    if (maxI < node1) {
        sz += size[maxI - 1][node1 - 1];
    }
    else {
        sz += size[node1 - 1][maxI - 1];
    }

    if (maxJ < node1) {
        sz += size[maxJ - 1][node1 - 1];
    }
    else {
        sz += size[node1 - 1][maxJ - 1];
    }

    int sizeI = superNodes[maxI - 1]->getNumNodes();
    int sizeJ = superNodes[maxJ - 1]->getNumNodes();
    int sizeB = superNodes[node1 - 1]->getNumNodes();
    int sizeK = sizeI + sizeJ;
    sz -= 2 * sizeI * sizeJ * sizeB * (numTaxa - sizeB - sizeK);
    return sz;
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

void Node::print(ostream& f, ostream& topStr, const Node *caller, const Node *root, int numTaxa) const
{
    if (leaf) {
        f << t[0]->getTSet()[0];
        topStr << t[0]->getTSet()[0];
        return;
    }

    int rank[2];
    Node *n[2];
    double weight[2];

    // print is called from one of the three neighbors, find other two neighbors
    for (int i = 2, pos = 1; pos >= 0; i--) {
        Node *temp = getNeighbor(i);
        if (temp != caller) {
            n[pos] = temp;
            weight[pos] =edges[i]->getWeight();
            if (i != 0 || root == this) {
                rank[pos--] = t[i]->getTSet()[0];
            }
            else {
                rank[pos--] = getComplement(t[1]->getTSet(), t[2]->getTSet(), numTaxa);
            }
        }
    }

    // check which neighbor has lower taxa and print it first
    f << "(";
    topStr << "(";
    if (rank[0] < rank[1]) {
        n[0]->print(f, topStr, this, root, numTaxa);
        f << ":" << fixed << weight[0] << ",";
        topStr << ",";
        n[1]->print(f, topStr, this, root, numTaxa);
        f << ":" << fixed << weight[1];
    }
    else {
        n[1]->print(f, topStr, this, root, numTaxa);
        f << ":" << fixed << weight[1] << ",";
        topStr << ",";
        n[0]->print(f, topStr, this, root, numTaxa);
        f << ":" << fixed << weight[0];

    }
    f << ")";
    topStr << ")";

}

void Edge::print(ostream& f,int numTaxa) const
{
    f << setw(3) << number << "(" << weight << ")" << ":";
    f << " n[" << nodes[0]->getNumber() << "] <-> " << "n[" << nodes[1]->getNumber() << "], split = ";
    f << endl;
}

// reroots the tree to parent node of node representing taxa
// Note: does not change the tree structure, just returns a 'new' topology string
void Tree::modifyOutgroup(int taxon, string& top, string& topWithWts) const
{
    stringstream f;// stream to print topology with weights
    f.precision(3);
    stringstream topStr;// stream to print topology without weights
    int index = 0;
    //identify index of leaf node representing this taxa.
    for (;index < numTaxa; index++) {
        if (nodes[index]->getTset(0)->getTSet()[0] == taxon) {
            break;
        }
    }

    Node *newRoot = nodes[index]->getNeighbor(0);
    Node *nbrs[3];
    int tset[3];
    for (int i = 0; i < 3; i++) {
        nbrs[i] = newRoot->getNeighbor(i);
        tset[i] = newRoot->getTset(i)->getTSet()[0];
    }

    if (newRoot != root) {
        // since new root is not actually a root, taxonsets on 3 edges will be {a}U{b}, {a}, {b}
        // change the first set to: compliment of {a}U{b}
        tset[0] = getComplement(newRoot->getTset(1)->getTSet(), newRoot->getTset(2)->getTSet(), numTaxa);
    }

    f << "(";
    topStr << "(";
    if (tset[0] < tset[1]) {
        if (tset[0] < tset[2]) {
            nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(0)->getWeight() << ",";
            topStr << ",";
            if (tset[1] < tset[2]) {
                nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
                f << ":" << fixed << newRoot->getEdge(1)->getWeight() << ",";
                topStr << ",";
                nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
                f << ":" << fixed << newRoot->getEdge(2)->getWeight();
            }
            else {
                nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
                f << ":" << fixed << newRoot->getEdge(2)->getWeight() << ",";
                topStr << ",";
                nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
                f << ":" << fixed << newRoot->getEdge(1)->getWeight();
            }
        }
        else {
            nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(2)->getWeight() << ",";
            topStr << ",";
            nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(0)->getWeight() << ",";
            topStr << ",";
            nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(1)->getWeight();
        }
    }
    else if (tset[1] < tset[2]) {
        nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
        f << ":" << fixed << newRoot->getEdge(1)->getWeight() << ",";
        topStr << ",";
        if (tset[0] < tset[2]) {
            nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(0)->getWeight() << ",";
            topStr << ",";
            nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(2)->getWeight();
        }
        else {
            nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(2)->getWeight() << ",";
            topStr << ",";
            nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
            f << ":" << fixed << newRoot->getEdge(0)->getWeight();
        }
    }
    else {
        nbrs[2]->print(f, topStr, newRoot, root, numTaxa);
        f << ":" << fixed << newRoot->getEdge(2)->getWeight() << ",";
        topStr << ",";
        nbrs[1]->print(f, topStr, newRoot, root, numTaxa);
        f << ":" << fixed << newRoot->getEdge(1)->getWeight() << ",";
        topStr << ",";
        nbrs[0]->print(f, topStr, newRoot, root, numTaxa);
        f << ":" << fixed << newRoot->getEdge(0)->getWeight();
    }

    f << ");";
    topStr << ");";

    topWithWts = f.str();
    top = topStr.str();
}

void Tree::construct(string& top) {
    int currentLeafNode = 0, currentInternalNode=numTaxa,currentLeafEdge=0,currentInternalEdge=numTaxa;
    Node *n1,*n2,*n3;
    char ch;
    istringstream topStr(top);
    topStr >> ch;
    n1 = connectInt(topStr, currentLeafNode, currentInternalNode, currentLeafEdge, currentInternalEdge);
    topStr >> ch;
    n2 = connectInt(topStr, currentLeafNode, currentInternalNode, currentLeafEdge, currentInternalEdge);
    topStr >> ch;
    n3 = connectInt(topStr, currentLeafNode, currentInternalNode, currentLeafEdge, currentInternalEdge);
    root = nodes[currentInternalNode++];
    Edge *e1, *e2, *e3;
    if (n1->isLeaf())
        e1 = edges[currentLeafEdge++];
    else
        e1 = edges[currentInternalEdge++];

    if (n2->isLeaf())
        e2 = edges[currentLeafEdge++];
    else
        e2 = edges[currentInternalEdge++];

    if (n3->isLeaf())
        e3 = edges[currentLeafEdge++];
    else
        e3 = edges[currentInternalEdge++];

    root->setEdge(0, e1);
    root->setEdge(1, e2);
    root->setEdge(2, e3);

    n1->setEdge(0, e1);
    n2->setEdge(0, e2);
    n3->setEdge(0, e3);

    e1->setNode(0, root);
    e1->setNode(1, n1);
    e2->setNode(0, root);
    e2->setNode(1, n2);
    e3->setNode(0, root);
    e3->setNode(1, n3);

    root->mergeTaxa(0, n1, 0);
    root->mergeTaxa(1, n2, 0);
    root->mergeTaxa(2, n3, 0);
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

// returns size of set {t1 intersect t2}
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

// returns tsets of size 4. Taxonsets on path from n1 to n2 are skipped
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
    for (int i = numTaxa; i < numNodes; i++) {
        Node *n1 = nodes[i];
        for (int j = i + 1; j < numNodes; j++) {
            Node *n2 = nodes[j];
            getQuartets(rInd, cInd, n1 , n2);
        }
    }
}

int Tree::getQuartets(vector<int>& rInd, vector<int>& cInd, Node *n1, Node *n2) {
//    for (int i = numTaxa; i < numNodes; i++) {

//        for (int j = i + 1; j < numNodes; j++) {

            vector<vector<int> > tsets;
            getTsets(tsets, n1, n2);
            assert(tsets.size() == 4);
            for (int j1 = 0; j1 < tsets[0].size(); j1++) {
                for (int j2 = 0; j2 < tsets[1].size(); j2++) {
                    for (int j3 = 0; j3 < tsets[2].size(); j3++) {
                        for (int j4 = 0; j4 < tsets[3].size(); j4++) {
                            int rIndex, cIndex;
                            getQuartetRowColumnIndex(tsets[0][j1], tsets[1][j2], tsets[2][j3], tsets[3][j4], rIndex, cIndex);
                            rInd.push_back(rIndex);
                            cInd.push_back(cIndex);
                        }
                    }
                }
            }
//        }
//    }

            return tsets[0].size() * tsets[1].size() * tsets[2].size() * tsets[3].size();
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
