/*********************************************
 * Concert Model
 * Autor: Aldo
 * Data de criação: 10-04-2023
 * Problema de Localização de hubs p-mediana
 *********************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <chrono>

using namespace std;
ILOSTLBEGIN



int main(int argc, char* argv[]) {
	try {
		/**============================
		 *  (Command-line arguments: CAB-High-10.txt p alpha beta)
		 *=============================== */
		ifstream arq(argv[1]);
		if (!arq.is_open()) {
			cout << "Error openning file: " << argv[1] << endl;
			arq.close();
			exit(EXIT_FAILURE);
		}

		int n; // Number of nodes
		int p; // Number of hubs
		int numcen; // Number of Scenarios
		double alpha; // discount factor
		double beta; // risk level
		
		int eps = 0.5;
		int eps2 = 1 - eps;		
		                
                arq >> n;
                arq >> numcen;

		if (argc >= 2) 
			p = atoi(argv[2]);
		else p = 2;

		alpha = atof(argv[3]);
		beta = atof(argv[4]);
		double m = numcen;
		double prob = 1 / m; // probability of each scenario (uniform)
		double beta2 = 1 / (1 - beta);
		std::vector<std::vector<std::vector<double>>> custo(numcen, std::vector<std::vector<double>>(n, std::vector<double>(n, 0)));
		std::vector<std::vector<std::vector<double>>> fluxo(numcen, std::vector<std::vector<double>>(n, std::vector<double>(n, 0)));
		vector<vector<double > > o(numcen, vector<double>(n));

		for (int m = 0; m < numcen; m++) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					arq >> custo[m][i][j];
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					arq >> fluxo[m][i][j];
				}
			}
		}
		
		for (int m = 0; m < numcen; m++) {
			for (int i = 0; i < n; i++) { 
				o[m][i] = 0.;
				for (int j = 0; j < n; j++) {
					o[m][i] += fluxo[m][i][j];
				}
			}
		}
		/**============================
		*  Model formulation
		*=============================== */

		IloEnv env;
		IloModel mod(env);
		IloCplex cplex(mod);


		IloNumVarArray h(env, n, 0, 1, ILOBOOL); // h[k] is a binary variable equal to 1 if node k is selected as a hub and 0 otherwise

		IloArray<IloArray<IloNumVarArray>> z(env, numcen); // z[m][i][k] represents the flow sent from node i to hub k under scenario m
		for (int m = 0; m < numcen; m++) {
			z[m] = IloArray<IloNumVarArray>(env, n);
			for (int i = 0; i < n; i++) {
				z[m][i] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}
		IloArray<IloArray<IloArray<IloNumVarArray>>> y(env, numcen); // y[m][i][k][l] represents the flow originating from node i and routed through hubs k and l under scenario m
		for (int m = 0; m < numcen; m++) {
			y[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
			for (int i = 0; i < n; i++) {
				y[m][i] = IloArray<IloNumVarArray>(env, n);
				for (int k = 0; k < n; k++) {
					y[m][i][k] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
				}
			}
		}
		IloArray<IloArray<IloArray<IloNumVarArray>>> x(env, numcen); // x[m][i][l][j] represents the flow originating from node i and sent from hub l to node j under scenario m
		for (int m = 0; m < numcen; m++) {
			x[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
			for (int i = 0; i < n; i++) {
				x[m][i] = IloArray<IloNumVarArray>(env, n);
				for (int l = 0; l < n; l++) {
					x[m][i][l] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
				}
			}
		}
		IloNumVarArray custocenario(env, numcen, 0, IloInfinity, ILOFLOAT); // custocenario[m] represents the total cost under scenario m

		IloNumVar eta(env, -IloInfinity, IloInfinity, ILOFLOAT);

		IloNumVarArray v(env, numcen, 0, IloInfinity, ILOFLOAT); // v[m] used to linearize custo[m] in the CVaR formulation


		// ====================================Formulation=============================================================== 
		// minimize (1-eps) * E(C) + eps * (eta + beta_2 * sum (m in 1..n) prob*v[m])
		// ===============================================================================================================

		IloExpr expfo(env);
		expfo += eps * eta;
		for (int i = 0; i < numcen; i++) {
			expfo += eps * beta2 * prob * v[i];
		}
		for (int i = 0; i < numcen; i++) {
			expfo += eps2 * prob * custocenario[i];
		}
		IloAdd(mod, IloMinimize(env, expfo));
		expfo.end();

		// Set of constraints

		//===========================================================================
		// forall (m in 1..numcen) 
		// custocenario[m] = sum(i in 1 .. n, k in 1..n) cen[m].custo[i][k] * z[m][i][k]+
		//                 +  sum(i in 1 .. n, k in 1..n, l in 1..n) alpha * cen[m].custo[k][l] * y[m][i][k][l]
		//                 +  sum(i in 1 .. n, l in 1..n, j in 1..n) cen[m].custo[l][j] * x[m][i][l][j]
		//===========================================================================
		for (int m = 0; m < numcen; m++) {
			IloExpr r1(env);
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < n; k++) {
					r1 += custo[m][i][k] * z[m][i][k];
				}
			}
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < n; k++) {
					for (int l = 0; l < n; l++) {
						r1 += alpha * custo[m][k][l] * y[m][i][k][l];
					}
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					for (int l = 0; l < n; l++) {
						r1 += custo[m][l][j] * x[m][i][l][j];
					}
				}
			}
			mod.add(r1 == custocenario[m]);
			r1.end();
		}

		//===========================================================================
		// sum(k in 1..n) h[i] == p 
		//===========================================================================
		IloExpr r2(env);
		for (int i = 0; i < n; i++) {
			r2 += h[i];
		}
		mod.add(r2 == p);
		r2.end();

		//===========================================================================
		// forall(m in 1..numcen, i in 1..n)  sum(k in 1..n) z[m][i][k] = o[m][i]   
		//===========================================================================
		for (int m = 0; m < numcen; m++) {
			for (int i = 0; i < n; i++) {
				IloExpr r3(env);
				for (int k = 0; k < n; k++) {
					r3 += z[m][i][k];
				}
				mod.add(r3 == o[m][i]);
				r3.end();
			}
		}

		//===========================================================================
		// forall (m in 1..numcen, i in 1..n, j in 1..n)  sum(l in 1..n) x[m][i][j][l] == cen[m].fluxo[i][j]  
		//===========================================================================
		for (int m = 0; m < numcen; m++) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					IloExpr r4(env);
					for (int l = 0; l < n; l++) {
						r4 += x[m][i][l][j];
					}
					mod.add(r4 == fluxo[m][i][j]);
					r4.end();
				}
			}
		}
		//==================================================================================================================================================================================
		// forall(m in 1.numcen, i in 1..n, k in 1..n) 
		// sum(l in 1..n:l!=k) y[m][i][k][l]+sum(j in 1..n)x[m][i][k][j] == sum(l in 1..n:l!=k) y[m][i][l][k]+z[m][i][k]
		// aqui codificado como
		// forall(m in 1..numcen, i in 1..n, k in 1..n) 
		// (sum(j in 1..n) x[m][i][k][j]) + (sum(l in 1..n:l!=k) (y[m][i][k][l]-y[m][i][l][k]))-z[m][i][k] == 0
		//==================================================================================================================================================================================
		for (int m = 0; m < numcen; m++) {
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < n; k++) {
					IloExpr r5(env);
					for (int j = 0; j < n; j++) {
						r5 += x[m][i][k][j];
					}
					for (int l = 0; l < n; l++) {
						if (l != k) {
							r5 += y[m][i][k][l] - y[m][i][l][k];
						}
					}
					r5 -= z[m][i][k];
					mod.add(r5 == 0);
					r5.end();
				}
			}
		}
		//===========================================================================
		// forall (m in 1..numcen, i in 1..n, k in 1..n) z[m][i][k] <= o[m][i]*h[k]
		//===========================================================================
		for (int m = 0; m < numcen; m++) {
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < n; k++) {
					mod.add(z[m][i][k] <= o[m][i] * h[k]);
				}
			}
		}
		//===========================================================================
		// forall (m in 1..numcen, i in 1..n, j in 1..n, l in 1..n) x[m][i][l][j] <= cen[m].fluxo[i][j]*h[l]
		//===========================================================================
		for (int m = 0; m < numcen; m++) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					for (int l = 0; l < n; l++) {
						mod.add(x[m][i][l][j] <= fluxo[m][i][j] * h[l]);
					}
				}
			}
		}

		for (int m = 0; m < numcen; m++) {
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < n; k++) {
					IloExpr r6(env);
					for (int l = 0; l < n; l++) {
						//if (l != k) {
						r6 += y[m][i][k][l];
						//}
					}
					mod.add(r6 <= z[m][i][k]);
					r6.end();
				}
			}
		}

		
		
		//===========================================================================
		// forall (m in 1..numcen) v[m] >= custocenario[m]-eta
		//===========================================================================
		for (int m = 0; m < numcen; m++) {
			mod.add(v[m] >= custocenario[m] - eta);
		}

		/// ==========================
		/// configurações do Cplex
		/// ==========================
		cplex.setParam(IloCplex::EpGap, 0.0000001); 
		cplex.setParam(IloCplex::TiLim, 86400); 
		//cplex.setWarning(env.getNullStream()); 
		//cplex.setOut(mono->env.getNullStream()); 
		ofstream logs("Logs.txt", std::ofstream::app);
		cplex.setOut(logs);
		cplex.setParam(IloCplex::Threads, 1); //cplex.setParam(IloCplex::Param::Benders::Strategy, 3); 
		


		 ///==============================
		 /// Solving the problem
		 ///==============================

		IloTimer crono(env);
		double lb = 0;
		double ub = 10e-10;
		double gap;


		crono.start();
		cplex.solve();
		crono.stop();


		if (cplex.getStatus() == IloAlgorithm::Optimal) {
			lb = cplex.getObjValue();
			ub = cplex.getObjValue();
			gap = 0.0;
		}

		///=====================================
		/// Saving the results
		///=====================================
		/// 
		/// 
		cout << "Instancia:"<< argv[1] << endl;
		cout << "Numero de Cenarios:"<< numcen << endl;
		cout << "O valor de n:"<< n << endl;
		cout << "O valor de p:" << p << endl;
		cout << "O valor de alpha:" << alpha << endl;
		cout << "O valor de beta:" << beta << endl;


		printf("\n Valor de Funcao Objetivo %10.4f", cplex.getObjValue());

		cout << endl;

		printf("\n Tempo de CPU: " "%f\t \n", (double)crono.getTime());
		cout << endl;

		printf("Hubs Instalados:");
		for (int j = 0; j < n; j++) {
			if (cplex.getValue(h[j]) >= 0.1) {
				printf(" %d\t ", j + 1);
			}
		}
		cout << endl;
		cout << endl;

		/**=====================================
		*  Apresenta a configuração final
		* ====================================*/

		FILE *re;
		re = fopen("Results.txt","aw+");
		fprintf(re, "\n Informacoes Gerais: " "%s\t%d\t%d\t%10.2f\t%10.2f\n",  argv[1], n, p, alpha,beta);
		fprintf(re, "\n Valor funcao objetivo: " "%f\t \n", (double) cplex.getObjValue ());
		fprintf(re, "\n Tempo de CPU: " "%f\t \n", (double) crono.getTime());
		fprintf(re, "\nHubs Instalados:\n");
		for (int j = 0; j < n; j++) {
			if (cplex.getValue(h[j]) >= 0.1) {
				fprintf(re, "%d\t ", j + 1);
			}
		}
		fprintf(re, "\n======================================================================");
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	return 0;
}
