#include "Alias.h"

//// ps - probabilities for all topologies
//// is - indices for which probabilities are non-zero
Alias::Alias(vector<double>& ps, vector<int>& is) {
  indices = is;
  nCells = is.size();
  vector<int> lessThanOne;
  vector<int> gtrThanOne;
  int index = 0;
  for(vector<int>::iterator pit = indices.begin(); pit != indices.end(); pit++, index++) {
    double cutoff = ps[*pit] * nCells;
    cutoffs.push_back(cutoff);
    if (cutoff > 1.0) {
      gtrThanOne.push_back(index);
    }
    else {
      lessThanOne.push_back(index);
    }
  }

  aliases.resize(nCells);
  computeTwoPointDistributions(lessThanOne, gtrThanOne);
//  printTwoPointDistributions(ps, is, nCells);
}

int Alias::pick(Rand& rg) {
  double rand = rg.runif(); // rand is in [0, 1)
  double kfr = cutoffs.size() * rand;
  int k = (int) floor(kfr); // k is in [0, n)
  double r = kfr - (double) k;
  if (r > cutoffs[k]) {
    return aliases[k];
  }

  return indices[k];
}

void Alias::printTwoPointDistributions(vector<double>& ps, vector<int>& is, int numCells) {
  cerr << "\t";
  for(vector<int>::iterator itr = is.begin(); itr != is.end(); itr++) {
    cerr << ps[*itr] * numCells << "\t";
  }
  cerr << endl << endl;

  for (int i = 0; i < is.size(); i++) {
    double c1 = cutoffs[i];
    double c2 = 1.00000 - c1;
    int i1 = is[i];
    int i2 = is[aliases[i]];
    cerr << i + 1 << "\t";
    for (int j = 0; j < numCells; j++) {
      if (j == i1) {
        cerr << c1 << "\t";
      }
      else if (j == i2) {
        cerr << c2 << "\t";
      }
      else {
        cerr << "\t";
      }
    }

    cerr << endl;
  }

}

void Alias::computeTwoPointDistributions(vector<int>& lessThanOne, vector<int>& gtrThanOne) {
  vector<int>::iterator gitr = gtrThanOne.begin();
  if (gitr == gtrThanOne.end())
    return;

  for (int i = 0; i < lessThanOne.size(); i++) {
    int j = lessThanOne[i];
    int k = *gitr;
    aliases[j] = k;
    cutoffs[k] -= (1 - cutoffs[j]);
    if (cutoffs[k] < 1.0) {
      lessThanOne.push_back(k);
      gitr++;
      if (gitr == gtrThanOne.end())
        break;
    }
//    cerr << j << " " << k  << " " << cutoffs[j] << " " << 1 - cutoffs[j] << endl;
  }
}

// simple test

/*
int main(int argc, char *argv[]) {
  vector<double> ps;
  ps.push_back(.17);
  ps.push_back(.02);
  ps.push_back(.15);
  ps.push_back(.01);
  ps.push_back(.04);
  ps.push_back(.25);
  ps.push_back(.05);
  ps.push_back(.03);
  ps.push_back(.2);
  ps.push_back(.08);

  vector<int> is;
  for (int i = 0; i < ps.size(); i++) {
    is.push_back(i);
  }

  Alias a(ps, is);
}
*/
