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
#include "Solver.hpp"
#include "Instance.hpp"

using namespace std;

int main(int argc, char** argv) {

	string instPath = argv[1];

	srand(time(NULL));

	Solver jitJss;

	inst.parse(instPath);

	jitJss.gifflerThompson();

	double result = jitJss.scheduler();

	double resultRelax = jitJss.relax_1(jitJss.sequence);

	if (!jitJss.verifySchedule()) {
		cout << "Solucao nao viavel!!!\n\n";
	}
	else {
		cout << "initial sol: " << result << " relax2: " << resultRelax << "\n Solucao viavel!!!\n\n";
	}
}

