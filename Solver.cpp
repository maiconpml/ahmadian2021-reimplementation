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
#include <queue>
#include "Solver.hpp"
#include "Instance.hpp"

Inst inst;

void Inst::parse(const string& instPath) {

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

void Solver::gifflerThompson() {
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

bool Solver::verifySchedule() {

	for (unsigned i = 1; i < inst.nOpers; ++i) {
		for (unsigned j = 1; j < inst.nOpers; ++j) {
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

double Solver::scheduler(vector<unsigned> seq, vector<unsigned>& starts) {

	vector<unsigned> auxJobs = inst.roots;
	vector<vector<unsigned>> machOrder(inst.nMach, vector<unsigned>(0));
	double penalties = UINT32_MAX;

	starts.clear();
	starts.push_back(0);

	for (unsigned i = 0; i < seq.size(); ++i) {
		machOrder[inst.opToMach[auxJobs[seq[i]]]].push_back(auxJobs[seq[i]]);
		auxJobs[seq[i]] = inst.next[auxJobs[seq[i]]];
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

			/*jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= inst.deadlines[i] - completionTimes[i - 1]);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= completionTimes[i - 1] - inst.deadlines[i]);
			jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= 0);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= 0);*/
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
		jitCplex.setParam(IloCplex::Param::TimeLimit, 5);
		jitCplex.setOut(jitEnv.getNullStream());

		jitCplex.extract(jitModel);
		jitCplex.solve();

		for (unsigned i = 1; i < inst.nOpers; ++i) {
			starts.push_back(round(jitCplex.getValue(completionTimes[i - 1]) - inst.costs[i]));
		}

		penalties = jitCplex.getObjValue();

		/*jitCplex.end();
		objExpr.end();*/
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}

	//jitModel.end();
	completionTimes.end();
	jitEnv.end();

	return penalties;
}

double Solver::relax_1(vector<unsigned>& seq, vector<unsigned>& starts) {
	vector<unsigned> relaxedMachs;
	unsigned mach1 = rand() % inst.nMach;
	unsigned mach2 = rand() % inst.nMach;

	int rcValue = rc() / inst.nMach;
	int rrValue = rr();
	rrValue /= 2;

	vector<vector<unsigned>> machOrder(inst.nMach, vector<unsigned>(0));
	vector<unsigned> auxJobs = inst.roots;
	for (unsigned i = 0; i < seq.size(); ++i) {
		machOrder[inst.opToMach[auxJobs[seq[i]]]].push_back(auxJobs[seq[i]]);
		auxJobs[seq[i]] = inst.next[auxJobs[seq[i]]];
	}

	unsigned start = rcValue - rrValue < 0 ? 0 : rcValue - rrValue;
	unsigned end = rcValue + rrValue > (int)inst.nJobs ? inst.nJobs : rcValue + rrValue;
	vector<unsigned> mach1Relaxed;
	for (unsigned i = start; i < end; ++i) {
		mach1Relaxed.push_back(machOrder[mach1][i]);
	}
	vector<unsigned> mach2Relaxed;
	for (unsigned i = start; i < end; ++i) {
		mach2Relaxed.push_back(machOrder[mach2][i]);
	}

	starts.clear();
	starts.push_back(0);

	double penalties;
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

			/*jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= inst.deadlines[i] - completionTimes[i - 1]);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= completionTimes[i - 1] - inst.deadlines[i]);
			jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= 0);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= 0);*/
		}

		for (unsigned i = 1; i < inst.nOpers; i += inst.nMach) {
			jitModel.add(completionTimes[i - 1] - inst.costs[i] >= 0);
		}

		for (unsigned i = 0; i < mach2Relaxed.size(); ++i) {
			for (unsigned j = i + 1; j < mach2Relaxed.size(); ++j) {
				IloOr jitOr1(jitEnv);
				IloOr jitOr2(jitEnv);
				jitOr1.add(completionTimes[mach1Relaxed[i] - 1] <= completionTimes[mach1Relaxed[j] - 1] - inst.costs[mach1Relaxed[j]]);
				jitOr1.add(completionTimes[mach1Relaxed[j] - 1] <= completionTimes[mach1Relaxed[i] - 1] - inst.costs[mach1Relaxed[i]]);
				jitModel.add(jitOr1);
				if (mach2 != mach1) {
					jitOr2.add(completionTimes[mach2Relaxed[i] - 1] <= completionTimes[mach2Relaxed[j] - 1] - inst.costs[mach2Relaxed[j]]);
					jitOr2.add(completionTimes[mach2Relaxed[j] - 1] <= completionTimes[mach2Relaxed[i] - 1] - inst.costs[mach2Relaxed[i]]);
					jitModel.add(jitOr2);
				}
			}
		}

		for (unsigned m = 0; m < inst.nMach; ++m) {
			vector<unsigned> relaxedOpers;
			if (m == mach1) {
				relaxedOpers = mach1Relaxed;
			}
			else if (m == mach2) {
				relaxedOpers = mach2Relaxed;
			}
			for (unsigned j = 0; j < inst.nJobs - 1; ++j) {
				if (m == mach1 || m == mach2) {
					if (machOrder[m][j + 1] == relaxedOpers[0]) {
						for (unsigned k = 0; k < relaxedOpers.size(); ++k) {
							jitModel.add(completionTimes[machOrder[m][j] - 1] <= completionTimes[relaxedOpers[k] - 1] - inst.costs[relaxedOpers[k]]);
						}
					}
					else if (machOrder[m][j] == relaxedOpers[0]) {
						j += relaxedOpers.size();
						if (j < inst.nJobs) {
							for (unsigned k = 0; k < relaxedOpers.size(); ++k) {
								jitModel.add(completionTimes[relaxedOpers[k] - 1] <= completionTimes[machOrder[m][j] - 1] - inst.costs[machOrder[m][j]]);
							}
						}
						--j;
					}
					else {
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
		jitCplex.setParam(IloCplex::Param::TimeLimit, 3);
		jitCplex.setOut(jitEnv.getNullStream());

		jitCplex.extract(jitModel);
		jitCplex.solve();

		for (unsigned i = 1; i < inst.nOpers; ++i) {
			starts.push_back(round(jitCplex.getValue(completionTimes[i - 1]) - inst.costs[i]));
		}

		penalties = jitCplex.getObjValue();
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}

	//completionTimes.end();
	jitEnv.end();

	priority_queue<pair<unsigned, unsigned>, vector<pair<unsigned, unsigned>>, greater<pair<unsigned, unsigned>>> startsQ;

	for (unsigned i = 1; i < inst.nOpers; ++i) {
		startsQ.push(pair<unsigned, unsigned>(starts[i], i));
	}
	seq.clear();
	for (unsigned i = 1; i < inst.nOpers; ++i) {
		unsigned curOp = startsQ.top().second;
		startsQ.pop();

		seq.push_back(inst.opToJob[curOp]);
	}

	return penalties;
}

double Solver::relax_2(vector<unsigned>& seq, vector<unsigned>& starts) {

	int rrValue = rr();
	int rcValue = rc();

	vector<unsigned> auxJobs = inst.roots;
	vector<vector<unsigned>> machOrder(inst.nMach, vector<unsigned>(0));
	vector<bool> willRelax(seq.size(), false);
	vector<vector<unsigned>> relaxedOpersInMach(inst.nMach, vector<unsigned>(0));
	unsigned seqCount = rcValue - rrValue < 0 ? 0 : rcValue - rrValue;
	unsigned seqCount2 = rrValue + rcValue > (int)seq.size() - 1 ? seq.size() - 1 : rrValue + rcValue;
	for (unsigned i = seqCount; i <= seqCount2; ++i) {
		willRelax[i] = true;
	}

	for (unsigned i = 0; i < seq.size(); ++i) {
		machOrder[inst.opToMach[auxJobs[seq[i]]]].push_back(auxJobs[seq[i]]);
		if (willRelax[auxJobs[seq[i]]]) relaxedOpersInMach[inst.opToMach[auxJobs[seq[i]]]].push_back(auxJobs[seq[i]]);
		auxJobs[seq[i]] = inst.next[auxJobs[seq[i]]];
	}
	vector<vector<unsigned>> relaxedBlocks;

	for (unsigned i = 0; i < inst.nMach; ++i) {
		if (relaxedOpersInMach[i].size() > 1) {
			unsigned inMachIndex = 0;
			unsigned curOp = 0;
			while (relaxedOpersInMach[i][0] != machOrder[i][inMachIndex++]);
			vector<unsigned> relaxedBlock;
			relaxedBlock.push_back(relaxedOpersInMach[i][curOp++]);
			while (curOp < relaxedOpersInMach[i].size()) {
				if (relaxedOpersInMach[i][curOp] != machOrder[i][inMachIndex]) {
					if (relaxedBlock.size() > 1) {
						relaxedBlocks.push_back(relaxedBlock);
					}
					relaxedBlock.clear();
					while (relaxedOpersInMach[i][curOp] != machOrder[i][inMachIndex++]);
					--inMachIndex;
				}
				relaxedBlock.push_back(relaxedOpersInMach[i][curOp++]);
				inMachIndex++;
			}
			if (relaxedBlock.size() > 1) {
				relaxedBlocks.push_back(relaxedBlock);
			}
		}
	}


	starts.clear();
	starts.push_back(0);

	IloEnv jitEnv;
	IloModel jitModel(jitEnv);
	double penalties = UINT32_MAX;

	IloNumVarArray completionTimes(jitEnv, inst.nOpers - 1, 0, IloInfinity, ILOINT);
	try {
		IloExpr objExpr(jitEnv);
		for (unsigned i = 1; i < inst.nOpers; ++i) {

			objExpr += (IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) * inst.earlCoefs[i]) + (IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) * inst.tardCoefs[i]);

			if (i % inst.nMach) jitModel.add(completionTimes[i - 1] <= completionTimes[i] - inst.costs[i + 1]);

			/*jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= inst.deadlines[i] - completionTimes[i - 1]);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= completionTimes[i - 1] - inst.deadlines[i]);
			jitModel.add(IloMax(inst.deadlines[i] - completionTimes[i - 1], 0) >= 0);
			jitModel.add(IloMax(completionTimes[i - 1] - inst.deadlines[i], 0) >= 0);*/
		}

		for (unsigned i = 1; i < inst.nOpers; i += inst.nMach) {
			jitModel.add(completionTimes[i - 1] - inst.costs[i] >= 0);
		}

		for (unsigned i = 0; i < relaxedBlocks.size(); ++i) {
			for (unsigned j = 0; j < relaxedBlocks[i].size(); ++j) {
				for (unsigned k = j + 1; k < relaxedBlocks[i].size(); ++k) {
					IloOr jitOr(jitEnv);
					jitOr.add(completionTimes[relaxedBlocks[i][j] - 1] <= completionTimes[relaxedBlocks[i][k] - 1] - inst.costs[relaxedBlocks[i][k]]);
					jitOr.add(completionTimes[relaxedBlocks[i][k] - 1] <= completionTimes[relaxedBlocks[i][j] - 1] - inst.costs[relaxedBlocks[i][j]]);
					jitModel.add(jitOr);
				}
			}
		}

		for (unsigned m = 0; m < inst.nMach; ++m) {
			for (unsigned j = 0; j < inst.nJobs - 1; ++j) {

				int block = isFirstInSomeBlock(relaxedBlocks, machOrder[m][j]);
				int blockNext = isFirstInSomeBlock(relaxedBlocks, machOrder[m][j + 1]);

				if (block == -1 && blockNext == -1) {
					jitModel.add(completionTimes[machOrder[m][j] - 1] <= completionTimes[machOrder[m][j + 1] - 1] - inst.costs[machOrder[m][j + 1]]);
				}
				else if (blockNext != -1) {
					for (unsigned i = 0; i < relaxedBlocks[blockNext].size(); ++i) {
						jitModel.add(completionTimes[machOrder[m][j] - 1] <= completionTimes[relaxedBlocks[blockNext][i] - 1] - inst.costs[relaxedBlocks[blockNext][i]]);
					}
				}
				else if (block != -1) {
					j += relaxedBlocks[block].size();
					if (j < inst.nJobs) {
						for (unsigned i = 0; i < relaxedBlocks[block].size(); ++i) {
							jitModel.add(completionTimes[relaxedBlocks[block][i] - 1] <= completionTimes[machOrder[m][j] - 1] - inst.costs[machOrder[m][j]]);
						}
					}
					--j;
				}
			}
		}

		jitModel.add(IloMinimize(jitEnv, objExpr));

		IloCplex jitCplex(jitEnv);
		jitCplex.setParam(IloCplex::Param::TimeLimit, 3);
		jitCplex.setOut(jitEnv.getNullStream());

		jitCplex.extract(jitModel);
		jitCplex.solve();

		for (unsigned i = 1; i < inst.nOpers; ++i) {
			starts.push_back(round(jitCplex.getValue(completionTimes[i - 1]) - inst.costs[i]));
		}

		penalties = jitCplex.getObjValue();
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}

	completionTimes.end();
	jitEnv.end();

	priority_queue<pair<unsigned, unsigned>, vector<pair<unsigned, unsigned>>, greater<pair<unsigned, unsigned>>> startsQ;

	for (unsigned i = 1; i < inst.nOpers; ++i) {
		startsQ.push(pair<unsigned, unsigned>(starts[i], i));
	}
	seq.clear();
	for (unsigned i = 1; i < inst.nOpers; ++i) {
		unsigned curOp = startsQ.top().second;
		startsQ.pop();

		seq.push_back(inst.opToJob[curOp]);
	}

	return penalties;
}

double Solver::insert(vector<unsigned>& seq, vector<unsigned>& starts) {

	unsigned origin = rand() % seq.size();
	unsigned destination = rand() % seq.size();

	while (origin == destination) {
		origin = rand() % seq.size();
		destination = rand() % seq.size();
	}

	unsigned aux = seq[origin];

	seq.erase(seq.begin() + origin);
	if (origin < destination) --destination;
	seq.insert(seq.begin() + destination, aux);

	return scheduler(seq, starts);
}

double Solver::swap(vector<unsigned>& seq, vector<unsigned>& starts) {

	unsigned origin = rand() % seq.size();
	unsigned destination = rand() % seq.size();

	while (seq[origin] == seq[destination]) {
		origin = rand() % seq.size();
		destination = rand() % seq.size();
	}

	unsigned aux = seq[origin];
	seq[origin] = seq[destination];
	seq[destination] = aux;

	return scheduler(seq, starts);
}

void Solver::local_search(vector<unsigned>& seq, vector<unsigned>& starts, double& targetPenalties, unsigned k) {

	double (*n)(vector<unsigned>&, vector<unsigned>&);

	switch (k) {
	case 1:
		n = &Solver::relax_2;
		break;
	case 2:
		n = &Solver::insert;
		break;
	case 3:
		n = &Solver::swap;
		break;
	}

	vector<unsigned> tempSeq = seq;
	vector<unsigned> tempStarts = starts;
	unsigned iter = 0;
	unsigned iterMax = k == 1 ? 10 : 20;
	while (iter < iterMax) {
		double tempPenalties = n(tempSeq, tempStarts);

		if (tempPenalties < targetPenalties) {
			seq = tempSeq;
			starts = tempStarts;
			targetPenalties = tempPenalties;
		}
		else {
			++iter;
		}
	}
}

void Solver::vns() {

	vector<unsigned> seq = sequence;
	vector<unsigned> starts = startTimes;
	double auxPenalties = penalties;
	unsigned k = 1;
	unsigned iter = 0;
	unsigned iterMax = (inst.nJobs * inst.nMach) / 4;
	while (k <= 3 && iter < iterMax) {
		iter++;
		//shake phase
		auxPenalties = relax_1(seq, starts);

		local_search(seq, starts, auxPenalties, k);

		if (auxPenalties < penalties) {
			penalties = auxPenalties;
			startTimes = starts;
			sequence = seq;
		}
		else {
			k++;
		}
	}
}

void Solver::solve(string instPath) {

	inst.parse(instPath);


	gifflerThompson();

	double result = scheduler(sequence, startTimes);

	penalties = result;

	vns();


	verifySchedule();
}

int Solver::isFirstInSomeBlock(const vector<vector<unsigned>>& relaxedBlocks, const unsigned oper) {
	for (unsigned i = 0; i < relaxedBlocks.size(); ++i) {
		if (oper == relaxedBlocks[i][0]) return i;
	}
	return -1;
}

unsigned Solver::rc() {
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

unsigned Solver::rr() {
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