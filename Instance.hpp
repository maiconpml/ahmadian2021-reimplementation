#pragma once

#include <fstream>
#include <string.h>
#include <vector>
#include <set>
#include <cassert>
#include <climits>
#include <iostream>
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

class Inst {
public:
	void parse(const string& instPath);

	unsigned nJobs;
	unsigned nMach;
	unsigned nOpers;

	vector<unsigned> costs;
	vector<unsigned> deadlines;
	vector<double> earlCoefs;
	vector<double> tardCoefs;

	vector<unsigned> next;
	vector<unsigned> prev;

	vector<unsigned> opToMach;
	vector<unsigned> opToJob;

	vector<unsigned> roots;
	vector<unsigned> leafs;
};