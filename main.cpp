#include <string.h>
#include <iostream>
#include <stdio.h>
#include <chrono>
#include "Solver.hpp"
#include "Instance.hpp"

using namespace std;

int main(int argc, char** argv) {

	string instPath = argv[1];

	srand(time(NULL));

	Solver jitJss;
	double meanPenalties;
	double bestPenalties = INT32_MAX;
	double bestTime;
	unsigned runNTimes = 1;
	for (unsigned i = 0; i < runNTimes; ++i) {

		chrono::high_resolution_clock::time_point tpStart = chrono::high_resolution_clock::now();

		jitJss.solve(instPath);

		double millisecsFound = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tpStart).count();;

		if (jitJss.verifySchedule()) {
			meanPenalties += jitJss.penalties;
			if (bestPenalties > jitJss.penalties) {
				bestPenalties = jitJss.penalties;
				bestTime = millisecsFound / 1000;
			}
		}
	}

	cout << instPath << " " << bestPenalties << " " << meanPenalties/runNTimes << " " << bestTime << endl;

}

