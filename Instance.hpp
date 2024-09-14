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
#include "Settings.hpp"

class Inst {
public:
	void parse(const string& instPath) {

		double buffer;
		ifstream stream;

		stream.open(instPath);

		stream >> nJobs;
		stream >> nMach;

		costs.push_back(0);
		deadlines.push_back(0);
		earlCoefs.push_back(0);
		tardCoefs.push_back(0);
		next.push_back(0);
		prev.push_back(0);
		opToJob.push_back(0);
		opToMach.push_back(0);
		nOpers = 1;

		for (unsigned i = 0; i < nJobs; ++i) {
			prev.push_back(0);
			roots.push_back(nOpers);
			for (unsigned j = 0; j < nMach; ++j) {
				opToJob.push_back(i);
				stream >> buffer;
				opToMach.push_back(buffer);
				stream >> buffer;
				costs.push_back(buffer);
				stream >> buffer;
				deadlines.push_back(buffer);
				stream >> buffer;
				earlCoefs.push_back(buffer);
				stream >> buffer;
				tardCoefs.push_back(buffer);
				if (j > 0) prev.push_back(nOpers - 1);
				++nOpers;
				if (j < (nMach - 1)) next.push_back(nOpers);
			}
			leafs.push_back(nOpers - 1);
			next.push_back(0);
		}

	}

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

Inst inst;