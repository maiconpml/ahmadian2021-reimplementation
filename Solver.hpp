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
#include "Instance.hpp"

class Solver {
public:
	void gifflerThompson() {
		set<unsigned> ready;
		vector<unsigned> ready0;
		vector<unsigned> ready1;
		vector<unsigned> jobStartTimes(inst.nJobs, 0);
		vector<unsigned> machStartTimes(inst.nMach, 0);
		vector<unsigned> jobLeafs(inst.nJobs, 0);
		vector<unsigned> jobLeafs2(inst.nJobs, 0);
		vector<unsigned> machLeafs(inst.nMach, 0);
		vector<unsigned> machLeafs2(inst.nMach, 0);
		vector<unsigned> sequence2;
		vector<unsigned> job;
		vector<unsigned> machV;

		machV.resize(inst.nOpers);

		unsigned completionTime;
		unsigned earlCompletion;
		unsigned mach = 0;
		unsigned auxStartTime;

		sequence.clear();

		for (unsigned op : inst.roots) {
			ready.insert(op);
		}

		while (!ready.empty()) {

			earlCompletion = UINT_MAX;

			ready0.clear();
			ready1.clear();

			for (unsigned op : ready) {
				assert(op < inst.nOpers);
				completionTime = max(jobStartTimes[inst.opToJob[op]], machStartTimes[inst.opToMach[op]]) + inst.costs[op];
				if (completionTime < earlCompletion) {
					earlCompletion = completionTime;
					mach = inst.opToMach[op];
				}
			}

			for (unsigned op : ready) {
				if (inst.opToMach[op] == mach) {
					ready0.push_back(op);
				}
			}

			for (unsigned op : ready0) {
				auxStartTime = max(jobStartTimes[inst.opToJob[op]], machStartTimes[inst.opToMach[op]]);
				if (auxStartTime < earlCompletion) {
					ready1.push_back(op);
				}
			}

			assert(ready1.size() > 0);
			unsigned op = ready1[0];

			for (unsigned i = 1; i < ready1.size(); ++i) {
				if (inst.deadlines[op] > inst.deadlines[ready1[i]]) {
					op = ready1[i];
				}
			}

			if (jobLeafs[inst.opToJob[op]]) {
				//insertJobOper(jobLeafs[inst.opToJob[op]], jobLeafs2[inst.opToJob[op]], op);
				//sequence.push_back(inst.opToJob[op]);
				jobLeafs2[inst.opToJob[op]] = jobLeafs[inst.opToJob[op]];
			}

			jobLeafs[inst.opToJob[op]] = op;

			if (machLeafs[inst.opToMach[op]]) {
				//insertMachOper(machLeafs[inst.opToMach[op]], machLeafs2[inst.opToMach[op]], op);
				machV[machLeafs[inst.opToMach[op]]] = op;
				machLeafs2[inst.opToMach[op]] = machLeafs[inst.opToMach[op]];
			}
			sequence.push_back(inst.opToJob[op]);
			sequence2.push_back(op);

			machLeafs[inst.opToMach[op]] = op;

			ready.erase(op);

			if (inst.next[op]) {
				ready.insert(inst.next[op]);
			}

			auxStartTime = max(jobStartTimes[inst.opToJob[op]], machStartTimes[inst.opToMach[op]]);
			jobStartTimes[inst.opToJob[op]] = auxStartTime + inst.costs[op];
			machStartTimes[inst.opToMach[op]] = auxStartTime + inst.costs[op];
		}
	}

	bool verifySchedule() {

		for (unsigned i = 1; i < inst.nMach; ++i) {
			for (unsigned j = 1; j < inst.nJobs; ++j) {
				if (i != j) {
					if (inst.opToJob[i] == inst.opToJob[j]) {
						if (startTimes[i] < startTimes[j]) {
							if (startTimes[i] + inst.costs[i] > startTimes[j]) {
								cout << "job not ok: op1: " << i << " op2:" << j << endl;
								cout << "starts: op1(start, cost): " << startTimes[i] << " " << inst.costs[i] << " op2(start, cost): " << startTimes[j] << " " << inst.costs[j] << endl;
								return false;
							}
						}
						else {
							if (startTimes[j] + inst.costs[j] > startTimes[i]) {
								cout << "job not ok: op1: " << j << " op2:" << i << endl;
								cout << "starts: op1(start, cost): " << startTimes[j] << " " << inst.costs[j] << " op2(start, cost): " << startTimes[i] << " " << inst.costs[i] << endl;
								return false;
							}
						}
					}

					if (inst.opToMach[i] == inst.opToMach[j]) {
						if (startTimes[i] < startTimes[j]) {
							if (startTimes[i] + inst.costs[i] > startTimes[j]) {
								cout << "mach not ok: op1: " << i << " op2:" << j << endl;
								cout << "starts: op1(start, cost): " << startTimes[i] << " " << inst.costs[i] << " op2(start, cost): " << startTimes[j] << " " << inst.costs[j] << endl;
								return false;
							}
						}
						else {
							if (startTimes[j] + inst.costs[j] > startTimes[i]) {
								cout << "mach not ok: op1: " << j << " op2:" << i << endl;
								cout << "starts: op1(start, cost): " << startTimes[j] << " " << inst.costs[j] << " op2(start, cost): " << startTimes[i] << " " << inst.costs[i] << endl;
								return false;
							}
						}
					}
				}
			}
		}

		return true;
	}

	double scheduler() {

		vector<unsigned> auxJobs = inst.roots;
		vector<vector<unsigned>> machOrder(inst.nMach, vector<unsigned>(0));
		unsigned penalties = UINT32_MAX;

		startTimes.clear();
		startTimes.push_back(0);

		for (unsigned i = 0; i < sequence.size(); ++i) {
			machOrder[inst.opToMach[auxJobs[sequence[i]]]].push_back(auxJobs[sequence[i]]);
			auxJobs[sequence[i]] = inst.next[auxJobs[sequence[i]]];
		}

		IloEnv jitEnv;
		IloModel jitModel(jitEnv);

		IloNumVarArray completionTimes(jitEnv, inst.nOpers - 1, 0, IloInfinity, ILOINT);
		try {
			IloExpr objExpr(jitEnv);
			for (unsigned i = 1; i < inst.nOpers; ++i) {

				objExpr += (IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) * inst.earlCoefs[i]) + (IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) * inst.tardCoefs[i]);

				if (i % inst.nMach) {
					jitModel.add(completionTimes[i - 1] <= completionTimes[i] - inst.costs[i + 1]);
				}

				jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= inst.deadlines[i] - completionTimes[i - 1]);
				jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= completionTimes[i - 1] - inst.deadlines[i]);
				jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= 0);
				jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= 0);
			}

			for (unsigned i = 1; i < inst.nOpers; i += inst.nMach) {
				jitModel.add(completionTimes[i - 1] - inst.costs[i] >= 0);
			}

			for (unsigned i = 0; i < inst.nMach; ++i) {
				for (unsigned j = 0; j < inst.nJobs - 1; ++j) {
					jitModel.add(completionTimes[machOrder[i][j] - 1] <= completionTimes[machOrder[i][j + 1] - 1] - inst.costs[machOrder[i][j + 1]]);
				}
			}

			jitModel.add(IloMinimize(jitEnv, objExpr));

			IloCplex jitCplex(jitEnv);
			jitCplex.setOut(jitEnv.getNullStream());

			jitCplex.extract(jitModel);
			jitCplex.solve();

			for (unsigned i = 1; i < inst.nOpers; ++i) {
				startTimes.push_back(jitCplex.getValue(completionTimes[i - 1]) - inst.costs[i]);
			}

			penalties = jitCplex.getObjValue();
		}
		catch (IloException& ex) {
			cerr << "Error: " << ex << endl;
		}
		catch (...) {
			cerr << "Error" << endl;
		}

		return penalties;
	}

	double relax_1(vector<unsigned>& seq) {
		vector<unsigned> relaxedMachs;
		unsigned mach1 = rand() % inst.nMach;
		unsigned mach2 = rand() % inst.nMach;

		int rcValue = rc()/inst.nMach;
		int rrValue = rr()/2;

		vector<vector<unsigned>> machOrder(inst.nMach, vector<unsigned>(0));
		vector<unsigned> auxJobs = inst.roots;
		for (unsigned i = 0; i < sequence.size(); ++i) {
			machOrder[inst.opToMach[auxJobs[sequence[i]]]].push_back(auxJobs[sequence[i]]);
			auxJobs[sequence[i]] = inst.next[auxJobs[sequence[i]]];
		}

		unsigned start = rcValue - rrValue < 0 ? 0 : rcValue - rrValue;
		unsigned end = rcValue + rrValue > inst.nJobs ? inst.nJobs : rcValue + rrValue;
		vector<unsigned> mach1Relaxed;
		for (unsigned i = start; i < end; ++i) {
			mach1Relaxed.push_back(machOrder[mach1][i]);
		}
		vector<unsigned> mach2Relaxed;
		for (unsigned i = start; i < end; ++i) {
			mach2Relaxed.push_back(machOrder[mach2][i]);
		}

		unsigned penalties;
		IloEnv jitEnv;
		IloModel jitModel(jitEnv);
		IloNumVarArray completionTimes(jitEnv, inst.nOpers - 1, 0, IloInfinity, ILOINT);
		try {
			IloExpr objExpr(jitEnv);
			for (unsigned i = 1; i < inst.nOpers; ++i) {

				objExpr += (IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) * inst.earlCoefs[i]) + (IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) * inst.tardCoefs[i]);

				if (i % inst.nMach) {
					jitModel.add(completionTimes[i - 1] <= completionTimes[i] - inst.costs[i + 1]);
				}

				jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= inst.deadlines[i] - completionTimes[i - 1]);
				jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= completionTimes[i - 1] - inst.deadlines[i]);
				jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= 0);
				jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= 0);
			}

			for (unsigned i = 1; i < inst.nOpers; i += inst.nMach) {
				jitModel.add(completionTimes[i - 1] - inst.costs[i] >= 0);
			}

			for (unsigned i = 0; i < mach2Relaxed.size(); ++i) {
				for (unsigned j = i+1; j < mach2Relaxed.size(); ++j) {
					jitModel.add(completionTimes[mach1Relaxed[i] - 1] <= completionTimes[mach1Relaxed[j] - 1] - inst.costs[mach1Relaxed[j]] || completionTimes[mach1Relaxed[j] - 1] <= completionTimes[mach1Relaxed[i] - 1] - inst.costs[mach1Relaxed[i]]);
				}
			}

			for (unsigned m = 0; m < inst.nMach; ++m) {
				vector<unsigned> relaxedOpers;
				if (m == mach1) {
					relaxedOpers = mach1Relaxed;
				}
				else if(m == mach2) {
					relaxedOpers = mach2Relaxed;
				}
				for (unsigned j = 0; j < inst.nJobs-1; ++j) {
					if (m == mach1 || m == mach2) {
						if (machOrder[m][j + 1] == relaxedOpers[0]) {
							for (unsigned k = 0; k < relaxedOpers.size(); ++k) {
								jitModel.add(completionTimes[machOrder[m][j] - 1] <= completionTimes[relaxedOpers[k] - 1] - inst.costs[relaxedOpers[k]]);
							}
						}
						else if (j && machOrder[m][j - 1] == relaxedOpers[relaxedOpers.size()-1]) {
							for (unsigned k = 0; k < relaxedOpers.size(); ++k) {
								jitModel.add(completionTimes[relaxedOpers[k] - 1] <= completionTimes[machOrder[m][j] - 1] - inst.costs[machOrder[m][j]]);
							}
							jitModel.add(completionTimes[machOrder[m][j] - 1] <= completionTimes[machOrder[m][j + 1] - 1] - inst.costs[machOrder[m][j + 1]]);
						}
						else if(machOrder[m][j] != relaxedOpers[0]) {
							jitModel.add(completionTimes[machOrder[m][j] - 1] <= completionTimes[machOrder[m][j + 1] - 1] - inst.costs[machOrder[m][j + 1]]);
						}
					}
					else {
						jitModel.add(completionTimes[machOrder[m][j] - 1] <= completionTimes[machOrder[m][j + 1] - 1] - inst.costs[machOrder[m][j + 1]]);
					}
				}
			}

			jitModel.add(IloMinimize(jitEnv, objExpr));

			IloCplex jitCplex(jitEnv);
			jitCplex.setOut(jitEnv.getNullStream());

			jitCplex.extract(jitModel);
			jitCplex.solve();

			for (unsigned i = 1; i < inst.nOpers; ++i) {
				startTimes.push_back(jitCplex.getValue(completionTimes[i - 1]) - inst.costs[i]);
			}

			penalties = jitCplex.getObjValue();
		}
		catch (IloException& ex) {
			cerr << "Error: " << ex << endl;
		}
		catch (...) {
			cerr << "Error" << endl;
		}

		return penalties;
	}

	double relax_2(vector<unsigned>& seq) {

		int rrValue = rr();
		int rcValue = rc();

		vector<unsigned> auxJobs = inst.roots;
		vector<vector<unsigned>> machOrder(inst.nMach, vector<unsigned>(0));
		vector<bool> willRelax(seq.size(), false);
		unsigned seqCount = rcValue - rrValue < 0 ? 0 : rcValue - rrValue;
		unsigned seqCount2 = rrValue + rcValue > seq.size() - 1 ? seq.size() - 1 : rrValue + rcValue;
		for (unsigned i = seqCount; i <= seqCount2; ++i) {
			willRelax[i] = true;
		}

		for (unsigned i = 0; i < seq.size(); ++i) {
			{
				machOrder[inst.opToMach[auxJobs[seq[i]]]].push_back(auxJobs[seq[i]]);
				auxJobs[seq[i]] = inst.next[auxJobs[seq[i]]];
			}
		}

		IloEnv jitEnv;
		IloModel jitModel(jitEnv);
		unsigned penalties = UINT32_MAX;

		IloNumVarArray completionTimes(jitEnv, inst.nOpers - 1, 0, IloInfinity, ILOINT);
		try {
			IloExpr objExpr(jitEnv);
			for (unsigned i = 1; i < inst.nOpers; ++i) {

				objExpr += (IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) * inst.earlCoefs[i]) + (IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) * inst.tardCoefs[i]);

				if (i % inst.nMach) jitModel.add(completionTimes[i - 1] <= completionTimes[i] - inst.costs[i + 1]);

				jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= inst.deadlines[i] - completionTimes[i - 1]);
				jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= completionTimes[i - 1] - inst.deadlines[i]);
				jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= 0);
				jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= 0);
			}

			for (unsigned i = 1; i < inst.nOpers; i += inst.nMach) {
				jitModel.add(completionTimes[i - 1] - inst.costs[i] >= 0);
			}

			/*for (unsigned i = 0; i < inst.nMach; ++i) {
				for (unsigned j = 0; j < inst.nJobs - 1; ++j) {
					jitModel.add(completionTimes[machOrder[i][j] - 1] <= completionTimes[machOrder[i][j + 1] - 1] - inst.costs[machOrder[i][j + 1]]);
				}
			}*/

			for (unsigned i = 0; i < inst.nMach; ++i) {
				vector<unsigned> relaxedOps(0);
				unsigned lastFixedOp = 0;
				for (unsigned j = 0; j < inst.nJobs - 1; ++j) {
					if (willRelax[machOrder[i][j]]) {
						relaxedOps.push_back(machOrder[i][j]);
						if (!willRelax[machOrder[i][j + 1]]) {
							jitModel.add(completionTimes[machOrder[i][j] - 1] <= completionTimes[machOrder[i][j + 1] - 1] - inst.costs[machOrder[i][j + 1]]);
						}
					}
					else {

						if (relaxedOps.size() > 1) {
							for (unsigned k = 0; k < relaxedOps.size(); ++k) {
								if (lastFixedOp) jitModel.add(completionTimes[lastFixedOp - 1] <= completionTimes[relaxedOps[k] - 1] - inst.costs[relaxedOps[k]]);
								jitModel.add(completionTimes[relaxedOps[k] - 1] <= completionTimes[machOrder[i][j] - 1] - inst.costs[machOrder[i][j]]);

								for (unsigned l = k + 1; l < relaxedOps.size(); ++l) {
									jitModel.add((completionTimes[relaxedOps[k] - 1] <= completionTimes[relaxedOps[l] - 1] - inst.costs[relaxedOps[l]]) || (completionTimes[relaxedOps[l] - 1] <= completionTimes[relaxedOps[k] - 1] - inst.costs[relaxedOps[k]]));
								}
							}
						}

						
						jitModel.add(completionTimes[machOrder[i][j] - 1] <= completionTimes[machOrder[i][j + 1] - 1] - inst.costs[machOrder[i][j + 1]]);
						

						lastFixedOp = machOrder[i][j];
						relaxedOps.clear();
					}

				}
				if (relaxedOps.size() > 1) {
					for (unsigned k = 0; k < relaxedOps.size(); ++k) {
						if (lastFixedOp) jitModel.add(completionTimes[lastFixedOp - 1] <= completionTimes[relaxedOps[k] - 1] - inst.costs[relaxedOps[k]]);

						for (unsigned l = k; l < relaxedOps.size(); ++l) {
							jitModel.add((completionTimes[relaxedOps[k] - 1] <= completionTimes[relaxedOps[l] - 1] - inst.costs[relaxedOps[l]]) || (completionTimes[relaxedOps[l] - 1] <= completionTimes[relaxedOps[k] - 1] - inst.costs[relaxedOps[k]]));
						}
					}
				}
				relaxedOps.clear();
			}

			jitModel.add(IloMinimize(jitEnv, objExpr));

			IloCplex jitCplex(jitEnv);
			jitCplex.setOut(jitEnv.getNullStream());

			jitCplex.extract(jitModel);
			jitCplex.solve();

			for (unsigned i = 1; i < inst.nOpers; ++i) {
				startTimes.push_back(jitCplex.getValue(completionTimes[i - 1]) - inst.costs[i]);
			}

			penalties = jitCplex.getObjValue();
		}
		catch (IloException& ex) {
			cerr << "Error: " << ex << endl;
		}
		catch (...) {
			cerr << "Error" << endl;
		}

		cout << "rc:" << rcValue << " rr: " << rrValue << endl;
		return penalties;
	}



	vector<unsigned> sequence;
	vector<unsigned> startTimes;
	unsigned penalties;

private:

	unsigned rc() {
		unsigned rcType = rand() % 100;
		vector<unsigned> chances = { 20, 20, 25, 35 };
		vector<unsigned> ranges = { 10, 30, 40, 20 };
		vector<unsigned> sizes = { (inst.nOpers - 1) * 10 / 100, (inst.nOpers - 1) * 30 / 100, (inst.nOpers - 1) * 40 / 100, (inst.nOpers - 1) * 20 / 100 };
		if (sizes[0] == 7) sizes[0]++;
		unsigned base = 0;
		unsigned middle = sizes[1] + sizes[2] + sizes[3];
		for (unsigned i = 0; i < chances.size(); ++i) {
			if (rcType < chances[i]) {
				unsigned rangeIndex = (rand() % sizes[i]);
				//cout << rcType << " " << i << " " << rangeIndex << " " << (inst.nOpers - 1) * ranges[i] / 100 << endl;
				return rangeIndex < sizes[i] / 2 ? base + rangeIndex : base + middle + rangeIndex;
			}
			rcType -= chances[i];
			base += sizes[0] / 2;
			middle -= sizes[i + 1];
		}
		return 0;
	}

	unsigned rr() {
		switch (inst.nJobs)
		{
		case 10:
			return (inst.nOpers - 1) / 3;
		case 15:
			return (inst.nOpers - 1) / 4;
		case 20:
			return (inst.nOpers - 1) / 5;
		default:
			break;
		}

		return 0;
	}
};