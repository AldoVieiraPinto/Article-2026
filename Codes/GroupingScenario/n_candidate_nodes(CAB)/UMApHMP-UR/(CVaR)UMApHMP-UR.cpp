/*********************************************
 * Concert Model
 * Autor: Aldo
 * Data de criação: 16-02-2025
 * Problema de Localização de hubs p-mediana 
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

#define inf 10e10

#ifndef CENARIO_H
#define CENARIO_H



struct Cenario{
			double prob;
			double **custo;
			double **fluxo;
		};
#endif // PESSOA_H



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
	std::vector<vector<int>>& hubs_grupo
	) {
	double c = cenarios_por_grupo;
	double prob = 1 / c; 
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

	IloNumVarArray h(env, n, 0, 1, ILOBOOL); 
	
	IloNumArray initialh(env);

	IloNumVar eta(env, 0.0, IloInfinity, ILOFLOAT);


	IloArray<IloArray<IloNumVarArray>> z(env, cenarios_por_grupo); 
	for (int m = 0; m < cenarios_por_grupo; m++) {
		z[m] = IloArray<IloNumVarArray>(env, n);
		for (int i = 0; i < n; i++) {
			z[m][i] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
		}
	}
	IloArray<IloArray<IloArray<IloNumVarArray>>> y(env, cenarios_por_grupo); 
	for (int m = 0; m < cenarios_por_grupo; m++) {
		y[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
		for (int i = 0; i < n; i++) {
			y[m][i] = IloArray<IloNumVarArray>(env, n);
			for (int k = 0; k < n; k++) {
				y[m][i][k] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}
	}
	IloArray<IloArray<IloArray<IloNumVarArray>>> x(env, cenarios_por_grupo); 
	for (int m = 0; m < cenarios_por_grupo; m++) {
		x[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
		for (int i = 0; i < n; i++) {
			x[m][i] = IloArray<IloNumVarArray>(env, n);
			for (int l = 0; l < n; l++) {
				x[m][i][l] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}
	}


	IloNumVarArray custocenario(env, cenarios_por_grupo, 0, IloInfinity, ILOFLOAT); 
	IloNumVarArray v(env, cenarios_por_grupo, 0, IloInfinity, ILOFLOAT); 

		// ====================================Formulation=============================================================== 
		// minimize eta + sum (m in 1..n) prob*v[m][k]
		// ===============================================================================================================

		IloExpr expfo(env);
		for (int m = 0; m < numcen; m++) {
			expfo += beta2 * prob * v[m];
		}
		expfo += eta;
		IloAdd(mod, IloMinimize(env, expfo));
		expfo.end();


	// Set of constraints

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
				mod.add(r6 <= o[m][i]*h[k]);
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
	/// Cplex settings
	/// ==========================
	cplex.setParam(IloCplex::EpGap, 0.000000001); 
	cplex.setParam(IloCplex::TiLim, 86400); 
	cplex.setWarning(env.getNullStream()); 
	ofstream logs("LogSolver.txt", std::ofstream::app);
	cplex.setOut(logs);
	cplex.setParam(IloCplex::Threads, 1); 

	
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
	std::vector<vector<int>>& hubs_grupo) {
	double cnt = numcen;
	double prob = 1 / (cnt); 
	vector<vector<double > > o(numcen, vector<double>(n));
	for (int m = 0; m < numcen; m++) {
		for (int i = 0; i < n; i++) { 
			o[m][i] = 0.;
			for (int j = 0; j < n; j++) {
				o[m][i] += fluxo[m][i][j];
			}
		}
	}

	IloEnv env;
	IloModel mod(env);
	IloCplex cplex(mod);

	IloNumVarArray h(env, n, 0, 1, ILOBOOL); 

	IloNumVar eta(env, 0.0, IloInfinity, ILOFLOAT);


	IloArray<IloArray<IloNumVarArray>> z(env, numcen); 
	for (int m = 0; m < numcen; m++) {
		z[m] = IloArray<IloNumVarArray>(env, n);
		for (int i = 0; i < n; i++) {
			z[m][i] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
		}
	}
	IloArray<IloArray<IloArray<IloNumVarArray>>> y(env, numcen); 
	for (int m = 0; m < numcen; m++) {
		y[m] = IloArray<IloArray<IloNumVarArray>>(env, n);
		for (int i = 0; i < n; i++) {
			y[m][i] = IloArray<IloNumVarArray>(env, n);
			for (int k = 0; k < n; k++) {
				y[m][i][k] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
			}
		}
	}
	IloArray<IloArray<IloArray<IloNumVarArray>>> x(env, numcen); 
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
	// minimize eta + sum (m in 1..n) prob*v[m][k]
	// ===============================================================================================================

	IloExpr expfo(env);
	for (int m = 0; m < numcen; m++) {
		expfo += beta2 * prob * v[m];
	}
	expfo += eta;
	IloAdd(mod, IloMinimize(env, expfo));
	expfo.end();

	// Set of constraints

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
				mod.add(r6 <= o[m][i]*h[k]);
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
	/// Cplex settings
	/// ==========================
	cplex.setParam(IloCplex::EpGap, 0.000000001); 
	cplex.setParam(IloCplex::TiLim, 86400); 
	cplex.setWarning(env.getNullStream()); 
	ofstream logs("LogSolver.txt", std::ofstream::app);
	cplex.setOut(logs);
	cplex.setParam(IloCplex::Threads, 1); 


	///==============================
	/// Solving the problem
	///==============================

	IloTimer crono(env);
	double lb = 0;
	double ub = 10e-10;
	double gap;

        //open hubs
	for (int i = 0; i < n; i++) {
		if (hubs_grupo[gr][i] > 0.1) {
			h[i].setBounds(1, 1);
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
		 * (Command-line arguments: CAB-High-10.txt p alpha beta)
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
		// Entrada via linha de comando: CAB-High10.txt p alpha beta total_grupos
		int n; // number of nodes
		int numcen;// number of scenarios
		arq >> n;
		arq >> numcen;
		int p;
		double alpha;
		double beta;
		if (argc >= 2) //Coletar números de hubs a serem ativados
		    p = atoi(argv[2]);
		else n = 20;
		alpha = atof(argv[3]);
		beta = atof(argv[4]);
		int total_grupos; // number of groups
		total_grupos = atoi(argv[7]);
		int cenarios_por_grupo = numcen/total_grupos;
		double beta2 = 1 / (1 - beta);
		std::vector<std::vector<std::vector<double>>> custo(numcen, std::vector<std::vector<double>>(n, std::vector<double>(n, 0)));
		std::vector<std::vector<std::vector<double>>> fluxo(numcen, std::vector<std::vector<double>>(n, std::vector<double>(n, 0)));
		std::vector<int> cenario(total_grupos);
		std::vector<double> solucao(total_grupos);
		std::vector<vector<int > > hubs_grupo(total_grupos, vector<int>(n)); 
		std::vector<double> s(n);
		std::vector<double> s_aux(n);

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
			DEP_grupo(n, p, alpha, beta, beta2, numcen, cenarios_por_grupo, gr, cenario, custo, fluxo, solucao, hubs_grupo);
		}
		double lb = media(total_grupos, solucao);
		
		// Calculo UB
		double ub = inf;
		double pior = 0;
		double auxiliar;
		int indice;
		for (int gr = 0; gr < total_grupos; gr++) {
			auxiliar = DEP_todos(gr, n, p, alpha, beta, beta2, numcen, custo, fluxo, hubs_grupo);
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
		//printf("\n\nHubs Instalados:\n\n");
		for (int j = 0; j < n; j++) {
			if (hubs_grupo[indice][j] >= 0.1) {
				printf(" %d\t ", j + 1);
			}
		}
		cout << endl;
		cout << endl;
		
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
				fprintf(re3, " %d\t ", j + 1);
			}
		}
		for (int j = 0; j < n; j++) {
			if (hubs_grupo[indice][j] >= 0.1) {
				fprintf(re3, " %d - ", j + 1);
			}
		}
		printf("\n");
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	return 0;
}
