#pragma once

#include <string.h>
#include <vector>

using namespace std;

class Inst {
public:
	void parse(const string& instPath);

	unsigned nJobs;								// number of jobs
	unsigned nMach;								// number of machines
	unsigned nOpers;							// number of operations

	vector<unsigned> costs;				// costs[o] is cost of operation o
	vector<unsigned> deadlines;		// deadlines[o] is deadline of operation o
	vector<double> earlCoefs;			// earlCoefs[o] is earliness penalty coeficient of operation o
	vector<double> tardCoefs;			// tardCoefs[o] is tardiness penalty coeficient of operation o

	vector<unsigned> next;				// next[o] is successor operation of o. 0 means no successor.
	vector<unsigned> prev;				// prev[o] is predecessor operation of o. 0 means no predecessor.

	vector<unsigned> opToMach;		// opToMach[o] is machine of o.
	vector<unsigned> opToJob;			// opToJob[o] is job of o.

	vector<unsigned> roots;				// roots[j] is first operation of job j.
	vector<unsigned> leafs;				// leafs[j] is last operation of job j.
};