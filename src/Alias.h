/*
 * Alias.h
 *
 *  Created on: Mar 17, 2010
 *      Author: satish
 */

#ifndef ALIAS_H_
#define ALIAS_H_
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

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

class Alias {
public:
  Alias(vector<double>& ps, vector<int>& is);

  virtual ~Alias() {
    cutoffs.clear();
    aliases.clear();
  }

  int pick(Rand& rg);

  void printTwoPointDistributions(vector<double>& ps, vector<int>& is, int numCells);
private:
  vector<double> cutoffs;
  vector<int> aliases;
  vector<int> indices;
  int nCells;

  void computeTwoPointDistributions(vector<int>& lessThanOne, vector<int>& gtrThanOne);
};

#endif /* ALIAS_H_ */
