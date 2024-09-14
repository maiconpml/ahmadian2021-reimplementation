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

	jitJss.solve(instPath);

}

