
/*********************************************

 * Concert Model
 * Autor: Aldo
 * Data de criação: 16-02-2025
 * Problema de Localização de hubs p-mediana (UMApHMP) com restrição de risco CVaR
 *********************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <chrono>
#include <time.h>
#include <cstdio>


using namespace std;
ILOSTLBEGIN

#define inf 10e10 /// 




double media(int tam, const std::vector<double>& vetor) {
	double soma = 0;
	for (int i = 0; i < tam; i++) {
		soma = soma + vetor[i];
	}
	soma = soma / tam;
	return soma;
}



void DEP_grupo(int n, int p, double alpha, double beta, double beta2, int numcen, int cenarios_por_grupo, int gr, const std::vector <int>& cenarios,
	const std::vector<std::vector<std::vector<double>>>& custo,
	const std::vector<std::vector<std::vector<double>>>& fluxo,
	std::vector<double>& solucao,
	std::vector<vector<int>>& hubs_grupo,
	std::vector<int>& s
	) {
	double c = cenarios_por_grupo;
	double prob = 1 / c; // probabilidade de cada cenário (uniforme)
	vector<vector<double > > o(cenarios_por_grupo, vector<double>(n));
	for (int m = 0; m < cenarios_por_grupo; m++) {
		for (int i = 0; i < n; i++) { 
			o[m][i] = 0.;
			for (int j = 0; j < n; j++) {
				o[m][i] += fluxo[(gr * cenarios_por_grupo)+m][i][j];
			}
		}
	}

	IloEnv env;
	IloModel mod(env);
	IloCplex cplex(mod);

	IloNumVarArray h(env, n, 0, 1, ILOBOOL); // h[k] variável binária que vale 1 se o nó k é um hub e 0 se não é ativado
	
	IloNumArray initialh(env);

	IloNumVar eta(env, 0.0, IloInfinity, ILOFLOAT);


	IloArray<IloArray<IloNumVarArray>> z(env, cenarios_por_grupo); // z[m][i][k] representa Fluxo enviado do nó i para o hub k no cenário m
	for (int m = 0; m < cenarios_por_grupo; m++) {
		z[m] = IloArray<IloNumVarArray>(env, n);
		for (int i = 0; i < n; i++) {
			z[m][i] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
		}
	}
	IloArray<IloArray<IloArray<IloNumVarArray>>> y(env, cenarios_por_grupo); // y[m][i][k][l] representa o fluxo originado do nó i e roteado por meio dos hubs k e l no cenário m
	for (int m = 0; m < cenarios_por_grupo; m++) {
		y[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
		for (int i = 0; i < n; i++) {
			y[m][i] = IloArray<IloNumVarArray>(env, n);
			for (int k = 0; k < n; k++) {
				y[m][i][k] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}
	}
	IloArray<IloArray<IloArray<IloNumVarArray>>> x(env, cenarios_por_grupo); // x[m][i][l][j] representa o fluxo originado do nó i e enviado do hub l para o nó j, no cenário m
	for (int m = 0; m < cenarios_por_grupo; m++) {
		x[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
		for (int i = 0; i < n; i++) {
			x[m][i] = IloArray<IloNumVarArray>(env, n);
			for (int l = 0; l < n; l++) {
				x[m][i][l] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}
	}


	IloNumVarArray custocenario(env, cenarios_por_grupo, 0, IloInfinity, ILOFLOAT); // custo[m] variável é o custo total no cenário m
	IloNumVarArray v(env, cenarios_por_grupo, 0, IloInfinity, ILOFLOAT); // custo[m] variável é o custo total no cenário m

		// ====================================Formulation=============================================================== 
		// minimize (1-eps) * E(C) + eps * (eta + beta_2 * sum (m in 1..n) prob*v[m])
		// ===============================================================================================================

	IloExpr expfo(env);
	expfo += eps * eta;
	for (int i = 0; i < numcen; i++) {
		expfo += 0.5 * beta2 * prob * v[i];
	}
	for (int i = 0; i < numcen; i++) {
		expfo += 0.5 * prob * custocenario[i];
	}
	IloAdd(mod, IloMinimize(env, expfo));
	expfo.end();

	// Restrições

	for (int m = 0; m < cenarios_por_grupo; m++) {
		IloExpr r1(env);
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < n; k++) {
				r1 += custo[(gr * cenarios_por_grupo) + m][i][k] * z[m][i][k];
			}
		}
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < n; k++) {
				for (int l = 0; l < n; l++) {
					r1 += alpha * custo[(gr * cenarios_por_grupo) + m][k][l] * y[m][i][k][l];
				}
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int l = 0; l < n; l++) {
					r1 += custo[(gr * cenarios_por_grupo) + m][l][j] * x[m][i][l][j];
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
	for (int m = 0; m < cenarios_por_grupo; m++) {
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
	for (int m = 0; m < cenarios_por_grupo; m++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				IloExpr r4(env);
				for (int l = 0; l < n; l++) {
					r4 += x[m][i][l][j];
				}
				mod.add(r4 == fluxo[(gr * cenarios_por_grupo) + m][i][j]);
				r4.end();
			}
		}
	}
	//=================================================================================================================================================================================
	// forall(m in 1.numcen, i in 1..n, k in 1..n) 
	// sum(l in 1..n:l!=k) y[m][i][k][l]+sum(j in 1..n)x[m][i][k][j] == sum(l in 1..n:l!=k) y[m][i][l][k]+z[m][i][k]
	// aqui codificado como
	// forall(m in 1..numcen, i in 1..n, k in 1..n) 
	// (sum(j in 1..n) x[m][i][k][j]) + (sum(l in 1..n:l!=k) (y[m][i][k][l]-y[m][i][l][k]))-z[m][i][k] == 0
	//=================================================================================================================================================================================
	for (int m = 0; m < cenarios_por_grupo; m++) {
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
	for (int m = 0; m < cenarios_por_grupo; m++) {
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < n; k++) {
				mod.add(z[m][i][k] <= o[m][i] * h[k]);
			}
		}
	}
	//===========================================================================
	// forall (m in 1..numcen, i in 1..n, j in 1..n, l in 1..n) x[m][i][l][j] <= cen[m].fluxo[i][j]*h[l]
	//===========================================================================
	for (int m = 0; m < cenarios_por_grupo; m++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int l = 0; l < n; l++) {
					mod.add(x[m][i][l][j] <= fluxo[(gr * cenarios_por_grupo) + m][i][j] * h[l]);
				}
			}
		}
	}
	
	for (int m = 0; m < cenarios_por_grupo; m++) {
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


	IloExpr r7(env);
	for (int m = 0; m < cenarios_por_grupo; m++) {
		mod.add(v[m] >= custocenario[m] - eta);
	}
	r7.end();

	/// ==========================
	/// configurações do Cplex
	/// ==========================
	cplex.setParam(IloCplex::EpGap, 0.000000001); 
	cplex.setParam(IloCplex::TiLim, 86400); 
	cplex.setWarning(env.getNullStream()); 
	ofstream logs("LogSolver.txt", std::ofstream::app);
	cplex.setOut(logs);
	cplex.setParam(IloCplex::Threads, 1); 
	
	///==============================
	/// Resolvendo o problema
	///==============================

	IloTimer crono(env);
	double lb = 0;
	double ub = 10e-10;
	double gap;
	
	for (int k = 16; k < 81; k++) {
	    h[s[k]].setBounds(0, 0);
	}
	for (int k = 16; k < 81; k++) {
	    for (int m = 0; m < cenarios_por_grupo; m++) {
	        for (int i = 0; i < 81; i++) {
	            z[m][i][s[k]].setBounds(0, 0);
	            for (int j = 0; j < 81; j++) {
	                x[m][i][s[k]][j].setBounds(0, 0);
	                y[m][i][s[k]][j].setBounds(0, 0);
	                y[m][i][j][s[k]].setBounds(0, 0);
	            }
	        }
	    }
	}


	crono.start();
	cplex.solve();
	crono.stop();


	if (cplex.getStatus() == IloAlgorithm::Optimal) {
		lb = cplex.getObjValue();
		ub = cplex.getObjValue();
		gap = 0.0;
	}


	solucao[gr] = double(cplex.getObjValue());
	for (int j = 0; j < n; j++){
		if (cplex.getValue(h[j]) >= 0.1) {
			hubs_grupo[gr][j] = 1;
		}
		else {
			hubs_grupo[gr][j] = 0;
		}
	}
	env.end();
}



double DEP_todos(int gr, int n, int p, double alpha, double beta, double beta2, int numcen,
	const std::vector<std::vector<std::vector<double>>>& custo,
	const std::vector<std::vector<std::vector<double>>>& fluxo,
	std::vector<vector<int>>& hubs_grupo, 
	std::vector<int>& s) {
	double cnt = numcen;
	double prob = 1 / (cnt); // probabilidade de cada cenário (uniforme)
	vector<vector<double > > o(numcen, vector<double>(n));
	for (int m = 0; m < numcen; m++) {
		for (int i = 0; i < n; i++) { //Inicializa vetor o de demandas
			o[m][i] = 0.;
			for (int j = 0; j < n; j++) {
				o[m][i] += fluxo[m][i][j];
			}
		}
	}

	IloEnv env;
	IloModel mod(env);
	IloCplex cplex(mod);

	IloNumVarArray h(env, n, 0, 1, ILOBOOL); // h[k] variável binária que vale 1 se o nó k é um hub e 0 se não é ativado

	IloNumVar eta(env, 0.0, IloInfinity, ILOFLOAT);


	IloArray<IloArray<IloNumVarArray>> z(env, numcen); // z[m][i][k] representa Fluxo enviado do nó i para o hub k no cenário m
	for (int m = 0; m < numcen; m++) {
		z[m] = IloArray<IloNumVarArray>(env, n);
		for (int i = 0; i < n; i++) {
			z[m][i] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
		}
	}
	IloArray<IloArray<IloArray<IloNumVarArray>>> y(env, numcen); // y[m][i][k][l] representa o fluxo originado do nó i e roteado por meio dos hubs k e l no cenário m
	for (int m = 0; m < numcen; m++) {
		y[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
		for (int i = 0; i < n; i++) {
			y[m][i] = IloArray<IloNumVarArray>(env, n);
			for (int k = 0; k < n; k++) {
				y[m][i][k] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}
	}
	IloArray<IloArray<IloArray<IloNumVarArray>>> x(env, numcen); // x[m][i][l][j] representa o fluxo originado do nó i e enviado do hub l para o nó j, no cenário m
	for (int m = 0; m < numcen; m++) {
		x[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
		for (int i = 0; i < n; i++) {
			x[m][i] = IloArray<IloNumVarArray>(env, n);
			for (int l = 0; l < n; l++) {
				x[m][i][l] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}
	}

	IloNumVarArray custocenario(env, numcen, 0, IloInfinity, ILOFLOAT); // custo[m] variável é o custo total no cenário m
	IloNumVarArray v(env, numcen, 0, IloInfinity, ILOFLOAT); // custo[m] variável é o custo total no cenário m

	// ====================================Formulation=============================================================== 
	// minimize (1-eps) * E(C) + eps * (eta + beta_2 * sum (m in 1..n) prob*v[m])
	// ===============================================================================================================

	IloExpr expfo(env);
	expfo += eps * eta;
	for (int i = 0; i < numcen; i++) {
		expfo += 0.5 * beta2 * prob * v[i];
	}
	for (int i = 0; i < numcen; i++) {
		expfo += 0.5 * prob * custocenario[i];
	}
	IloAdd(mod, IloMinimize(env, expfo));
	expfo.end();

	// Restrições

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
	//==========================================================================================================================================================
	// forall(m in 1.numcen, i in 1..n, k in 1..n) 
	// sum(l in 1..n:l!=k) y[m][i][k][l]+sum(j in 1..n)x[m][i][k][j] == sum(l in 1..n:l!=k) y[m][i][l][k]+z[m][i][k]
	// aqui codificado como
	// forall(m in 1..numcen, i in 1..n, k in 1..n) 
	// (sum(j in 1..n) x[m][i][k][j]) + (sum(l in 1..n:l!=k) (y[m][i][k][l]-y[m][i][l][k]))-z[m][i][k] == 0
	//==========================================================================================================================================================
	
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


	IloExpr r7(env);
	for (int m = 0; m < numcen; m++) {
		mod.add(v[m] >= custocenario[m] - eta);
	}
	r7.end();


	/// ==========================
	/// configurações do Cplex
	/// ==========================
	cplex.setParam(IloCplex::EpGap, 0.000000001); 
	cplex.setParam(IloCplex::TiLim, 86400); 
	cplex.setWarning(env.getNullStream()); 
	ofstream logs("LogSolver.txt", std::ofstream::app);
	cplex.setOut(logs);
	cplex.setParam(IloCplex::Threads, 1); 



	///==============================
	/// Resolvendo o problema
	///==============================

	IloTimer crono(env);// Variável para coletar o tempo
	double lb = 0;
	double ub = 10e-10;
	double gap;

	for (int i = 0; i < n; i++) {
		if (hubs_grupo[gr][i] > 0.1) {
			h[i].setBounds(1, 1);
		}
	}
	
	for (int k = 16; k < 81; k++) {
	    h[s[k]].setBounds(0, 0);
	}
	for (int k = 16; k < 81; k++) {
	    for (int m = 0; m < numcen; m++) {
	        for (int i = 0; i < 81; i++) {
	            z[m][i][s[k]].setBounds(0, 0);
	            for (int j = 0; j < 81; j++) {
	                x[m][i][s[k]][j].setBounds(0, 0);
	                y[m][i][s[k]][j].setBounds(0, 0);
	                y[m][i][j][s[k]].setBounds(0, 0);
	            }
	        }
	    }
	}


	crono.start();
	cplex.solve();
	crono.stop();


	if (cplex.getStatus() == IloAlgorithm::Optimal) {
		lb = cplex.getObjValue();
		ub = cplex.getObjValue();
		gap = 0.0;
	}

	double auxiliar = double(cplex.getObjValue());
	env.end();
	return auxiliar;
}




int main(int argc, char* argv[]) {
	try {
		/**============================
		 *  Leitura dos dados
		 *=============================== */
		clock_t inicio_CPU,       
			fim_CPU;
		inicio_CPU = clock();
		ifstream arq(argv[1]);
		if (!arq.is_open()) {
			cout << "Error openning file: " << argv[1] << endl;
			arq.close();
			exit(EXIT_FAILURE);
		}

		int n; 
		int p;
		int numcen;
		double alpha;
		double beta;
		if (argc >= 2) 
		    n = atoi(argv[2]);
		else n = 20;
		numcen = atoi(argv[3]);
		p = atoi(argv[4]);
		alpha = atof(argv[5]);
		beta = atof(argv[6]);
		int total_grupos;
		total_grupos = atoi(argv[7]);
		int cenarios_por_grupo = numcen/total_grupos;
		double beta2 = 1 / (1 - beta);
		std::vector<std::vector<std::vector<double>>> custo(numcen, std::vector<std::vector<double>>(n, std::vector<double>(n, 0)));
		std::vector<std::vector<std::vector<double>>> fluxo(numcen, std::vector<std::vector<double>>(n, std::vector<double>(n, 0)));
		std::vector<int> cenario(total_grupos);
		std::vector<double> solucao(total_grupos);
		std::vector<vector<int > > hubs_grupo(total_grupos, vector<int>(n)); 
		std::vector<int> s(n);


                
                for (int i = 0; i < n; i++) {
                        arq >> s[i];
                }
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
    		//gr é o contador do grupo
		for (int gr = 0; gr < total_grupos; gr++) {
			DEP_grupo(n, p, alpha, beta, beta2, numcen, cenarios_por_grupo, gr, cenario, custo, fluxo, solucao, hubs_grupo, s);
			printf("\n\nGrupo %d : %10.4f \n\n", gr, solucao[gr]);
		}
		double soma = 0;
		for(int i = 0; i < total_grupos; i++){
		      soma = soma + solucao[i];
		}
		double lb = soma/total_grupos; 
		// Calculo UB
		double ub = inf;
		double aux = 0;
		double pior = 0;
		double auxiliar;
		int indice;
		int indicepior;
		for (int gr = 0; gr < total_grupos; gr++) {
			auxiliar = DEP_todos(gr, n, p, alpha, beta, beta2, numcen, custo, fluxo, hubs_grupo, s);
			if (auxiliar < ub ){
			      ub = auxiliar;
			      indice = gr;
			}
		}

		
		double gap_cen = 100 * (ub - lb) / ub;
		fim_CPU = clock();
		double tempo = double(fim_CPU - inicio_CPU) / (CLOCKS_PER_SEC);
         	

		cout <<"Instancia:"<< argv[1] << endl;
		cout <<"Numero de Cenarios:"<< numcen << endl;
		cout <<"O valor de n:"<< n << endl;
		cout << "O valor de p:" << p << endl;
		cout << "O valor de alpha:" << alpha << endl;
		cout << "O valor de beta:" << beta << endl;
		printf("\nTempo: %10.4f\n", tempo);
		printf("\n *************** Solucao FINAL: **************************");
		printf("\n\n UB igual a \t %10.4f", ub);
		printf("\n\n LB igual a \t %10.4f", lb);
		printf("\n\n GAP igual a \t %10.4f \t por cento", gap_cen);
		printf("\n\nHubs Instalados:\n\n");
		
		int alpha0 = 10 * alpha;
		int beta0 = 10 * beta;
		
		FILE *re2;
		FILE *re3;
		re2 = fopen("Resultados.txt","aw+");
		re3 = fopen("Tabela.txt","aw+");
                fprintf(re2, "\n\n %s - %d - %d - %d - %d - %d \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t", argv[1], n, numcen , p, alpha0, beta0, tempo, ub, lb, gap_cen);
                fprintf(re3, "\n %s - %d - %d - %d - %d - %d \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t ", argv[1], n, numcen , p, alpha0, beta0, tempo, ub, lb, gap_cen);
		for (int j = 0; j < n; j++) {
			if (hubs_grupo[indice][j] >= 0.1) {
				fprintf(re2, " %d\t ", j + 1);
				fprintf(re3, " %d - ", j + 1);
			}
		}
        }
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	return 0;
}
