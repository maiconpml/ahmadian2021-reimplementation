#pragma once

#include <string.h>
#include <vector>

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