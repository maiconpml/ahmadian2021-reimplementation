#pragma once

#include <vector>

using namespace std;

class Solver {
public:

	void gifflerThompson();

	bool verifySchedule();

	static double scheduler(vector<unsigned> seq, vector<unsigned>& starts);

	static double relax_1(vector<unsigned>& seq, vector<unsigned>& starts);

	static double relax_2(vector<unsigned>& seq, vector<unsigned>& starts);

	static double insert(vector<unsigned>& seq, vector<unsigned>& starts);
	
	static double swap(vector<unsigned>& seq, vector<unsigned>& starts);

	void local_search(vector<unsigned>& seq, vector<unsigned>& starts, double& targetPenalties, unsigned k);

	void vns();
	
	void solve(string instPath);

	static int isFirstInSomeBlock(const vector<vector<unsigned>>& relaxedBlocks, const unsigned oper);

	static unsigned rc();

	static unsigned rr();

	vector<unsigned> sequence;
	vector<unsigned> startTimes;
	double penalties;
};