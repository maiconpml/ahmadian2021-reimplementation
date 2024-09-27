#pragma once

#include <vector>

using namespace std;

class Solver {
public:

	// sets initial solution to current sequence
	void gifflerThompson();

	// current solution is feasible?
	bool verifySchedule();

	// uses cplex to schedule seq and puts this schedule in starts. Returns total penalties
	static double scheduler(vector<unsigned> seq, vector<unsigned>& starts);

	// uses cplex to schedule seq but relaxing constraints of two machines and put the obtained schedule in starts. Returns total penalties
	static double relax_1(vector<unsigned>& seq, vector<unsigned>& starts);

	// uses cplex to schedule seq but relaxing constraints of randomly selected operations and put the obtained schedule in starts. Returns total penalties
	static double relax_2(vector<unsigned>& seq, vector<unsigned>& starts);

	// transform seq in a neighbour by using insert procedure. Put the obtained schedule in starts. Return total penalties.
	static double insert(vector<unsigned>& seq, vector<unsigned>& starts);
	
	// transform seq in a neighbout by using swap procedure. Put the obtained schedule in starts. Return total penalties. 
	static double swap(vector<unsigned>& seq, vector<unsigned>& starts);

	// try to improve solution represented by seq, starts and targetPenalties using neighbourhood Nk
	void local_search(vector<unsigned>& seq, vector<unsigned>& starts, double& targetPenalties, unsigned k);

	// try to improve current solution
	void vns();
	
	// solve JIT-JSS problem in file instPath
	void solve(string instPath);

	// if oper is the first element of some relaxedBlocks[i] returns i. If not return -1
	static int isFirstInSomeBlock(const vector<vector<unsigned>>& relaxedBlocks, const unsigned oper);

	// returns relaxation center
	static unsigned rc();

	// returns relaxation radius
	static unsigned rr();

	vector<unsigned> sequence;
	vector<unsigned> startTimes;
	double penalties;
};