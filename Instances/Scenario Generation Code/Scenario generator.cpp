/*********************************************
 * Concert Model
 * Autor: Aldo
 * Data de criação: 10-04-2023
 * Problema de Localização de hubs p-mediana (UMApHMP)
 *********************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h> 
#include <random>
#include <chrono>

using namespace std;
ILOSTLBEGIN

int main(int argc, char* argv[]) {
	try {
		/**============================
		 *  Data imput
		 *=============================== */
		ifstream arq(argv[1]);
		if (!arq.is_open()) {
			cout << "Error openning file: " << argv[1] << endl;
			arq.close();
			exit(EXIT_FAILURE);
		}
		int n; // number of nodes
		int numcen;
                if (argc >= 2) 
			numcen = atoi(argv[2]); // number of scenarios
		else numcen = 10;
		
		arq >> n;


		vector<vector<double > > w(n, vector<double>(n)); // a vector of vectors to represent an n × n matrix – amount of demand to be sent between nodes i and j
		vector<vector<double > > c(n, vector<double>(n)); // a vector of vectors to represent an n × n matrix storing the cost between nodes i and j
		vector<vector<double > > wl(n, vector<double>(n)); // a vector of vectors to represent an n × n matrix – amount of demand to be sent between nodes i and j under the scenario
		vector<vector<double > > cl(n, vector<double>(n)); // a vector of vectors to represent an n × n matrix storing the cost between nodes i and j under the scenario

			// ========================================================================================
			// (Command-line arguments: cab.txt numcen)
			// ========================================================================================


		for (int i = 0; i < n; i++) {   //Read cost data
			for (int j = 0; j < n; j++) {
				arq >> c[i][j];
			}
		}

		for (int i = 0; i < n; i++) {   //Read flows data
			for (int j = 0; j < n; j++) {
				arq >> w[i][j];
			}
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		// High
		double min = 0.5;
		double max = 1.5;
		// if Low: min = 0.9, max = 1.1
		// if Medium: min = 0.7, max = 1.3



		FILE* re;
		re = fopen("Instance.txt", "w+");

		// Create a uniform distribution over the interval [min, max]
		std::uniform_real_distribution<double> dis(min, max);

		
		for (int m = 0; m < numcen; m++) {


			for (int i = 0; i < n; i++) { 
				for (int j = 0; j < n; j++) {
					wl[i][j] = dis(gen) * w[i][j];
				}
			}

			for (int i = 0; i < n; i++) { 
				for (int j = 0; j < n; j++) {
					cl[i][j] = dis(gen) * c[i][j];
				}
			}

			

			for (int i = 0; i < n; i++) { 
				for (int j = 0; j < n; j++) {
					fprintf(re, "%10.5f \t ", cl[i][j]);
				}
				fprintf(re, "\n");
			}
			for (int i = 0; i < n; i++) { 
				for (int j = 0; j < n; j++) {
					fprintf(re, "%10.5f \t ", wl[i][j]);
				}
				fprintf(re, "\n");
			}
		}
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	return 0;
}
