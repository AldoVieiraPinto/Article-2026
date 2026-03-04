#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <vector>
#include <time.h>
#include <string>
#include <chrono>


using namespace std;

#define inf 10e10 /// 

#ifndef CENARIO_H
#define CENARIO_H

struct Cenario{
	double prob;
	double **custo;
	double **fluxo;
};
#endif 

// Auxiliary functions – memory allocation

void atualiza_vetor(int *s_star, int *s, int n)
{
   for (int j=0; j < n; j++) s_star[j] = s[j];
}


double **cria_matriz_double(int nlinhas, int ncolunas)
{
  register int i;
  double **matriz;

  matriz = (double **) malloc(nlinhas*sizeof(double *));
  if (!matriz) {
        printf("Falta memoria para alocar a matriz de ponteiros\n");
        exit(1);
  }
  for (i=0; i < nlinhas; i++) {
    matriz[i] = (double *) malloc(ncolunas*sizeof(double));
    if (!matriz[i]){
      printf("Falta memoria para alocar a matriz de ponteiros.\n");
      exit(1);
    }
  }
  return matriz;
}

/* Allocates memory for a vector with tam positions. */
int *cria_vetor(int tam)
{
  int *vetor;

  vetor = (int *) malloc(tam*sizeof(int));
  if (!vetor){
  	printf("Falta memoria para alocar o vetor de ponteiros");
    exit(1);
  }
  return vetor;
}


void inicializa_vetor_int(int *vetor, int tam)
{
    for (int j=0; j<tam; j++) vetor[j] = j;
}

void libera_vetor(int *vetor)
{
    free(vetor);
}

void libera_matriz_double(double **matriz, int nlinhas)
{
  register int i;

  for (i=nlinhas-1; i >= 0; i--)
    free((double *) matriz[i]);
  free((double *) matriz);
}


// Auxiliary functions – shuffle, sort

void embaralha_vetor(int *vetor, int n)
{
  int aux, j1, j2;
  for (int i=0; i < n; i++){
    j1 = rand() % (n);
    j2 = rand() % (n);
    while (j1 == j2) j2 = rand() % (n);
    aux = vetor[j1];
    vetor[j1] = vetor[j2];
    vetor[j2] = aux;
  }
}


void insertionSort( double* v, int *w, int n ) 
{
	int i= 0;
	int j= 1;
	double aux = 0;
	int aux2 = 0;
	while (j < n)
	{
		aux = v[j];
		aux2 = w[j];
		i = j - 1;
		while ((i >= 0) && (v[i] > aux))
			{
			v[i + 1] = v[i];
			w[i + 1] = w[i];
			i = i - 1;
			}
		v[i + 1] = aux;
		w[i + 1] = aux2;
	j = j + 1;
	}
}

void insertionSortWReal( double* v, double *w, int n ) 
{
	int i= 0;
	int j= 1;
	double aux = 0;
	double aux2 = 0;
	while (j < n)
	{
		aux = v[j];
		aux2 = w[j];
		i = j - 1;
		while ((i >= 0) && (v[i] > aux))
			{
			v[i + 1] = v[i];
			w[i + 1] = w[i];
			i = i - 1;
			}
		v[i + 1] = aux;
		w[i + 1] = aux2;
	        j = j + 1;
	}
}

// Risk functions

double media(double* c, double* pp, int size){
	double soma = 0.;
        for (int j = 0; j < size; j++){
		soma = soma + pp[j]*c[j];
	}
	return soma;
}


double cvar(double* c, double* pp, double beta, int size ){
	insertionSortWReal(c, pp, size);
	double soma = 0.;
	int i = 0;
	int kb;
	while((soma < beta-0.000000000000001)&&(i<size)){
		kb=i;
		soma = soma + pp[i];
	      	i = i+1;
	};
	double var = c[kb];
	double cvar = (soma-beta)*var; 
	for (int j = kb+1; j < size; j++){
		cvar = cvar + pp[j]*c[j];
	}
	cvar = cvar/(1-beta);
	return cvar;
}

double mcvar(double* c, double* pp, double beta, int size ){
	double neutro = media(c, pp, size);
	double risco = cvar(c, pp, beta, size); 
	double mcvar = (neutro + risco)/2;
	return mcvar;
}


// Movements

void embaralha_hub_e_naohub(int *vetor, int n, int p)
{
  int *haux = cria_vetor(p);
  int *nhaux = cria_vetor(n-p);
  for (int i = 0; i < p; i++){
		haux[i]=int(vetor[i]);
  }
  embaralha_vetor(haux, p);
  for (int i = 0; i < n-p; i++){
		nhaux[i]=int(vetor[p+i]);
  }
  for (int i = 0; i < p; i++){
		vetor[i]=haux[i];
  }
  embaralha_vetor(nhaux, n-p);
  for (int i = p; i < n; i++){
		vetor[i]=nhaux[i-p];
  }
  libera_vetor(haux);
  libera_vetor(nhaux);
  }
  
void movimento_troca_hub_naohub(int i, int j, int *s){
	int aux;
	aux = s[i];
	s[i] = s[j];
	s[j] = aux;
}



void movimento_ktrocas_hub_naohub(int p, int n, int *s, int k) {
	int aux;
	embaralha_hub_e_naohub(s, n, p);
  	for (int i = 0; i<k; i++){
		aux = s[i];
		s[i] = s[i+p];
		s[i+p] = aux;
    }    
}

// Greedy construction

void construcao_gulosa (int n, int *s, const Cenario cen[], int size){
  	double v[n];
	for (int k = 0; k < n; k++){
		v[k]=0;
		s[k]=k;
	}
  	for (int m = 0; m < size; m++){
		for (int i = 0; i<n; i++){
	    	for (int j = 0; j<n; j++){
	    		v[i]=v[i]+(cen[m].custo[i][j])*(cen[m].fluxo[i][j]);
	    	}
	    }
	}
	insertionSort(v, s, n);	
}


// Cost function -THR: two-hubs routes
double floyd_thr(int n, int p, double alpha, int* s, int w, const Cenario cen[], int size){
	double** d = cria_matriz_double(n, n);
	double** c = cria_matriz_double(n, n);
	double** e = cria_matriz_double(n, n);
	for(int i=0 ; i < p ; i++){
		for (int j=0; j<n; j++){
			if(j<p){
				c[s[i]][s[j]]=cen[w].custo[s[i]][s[j]];
			} 
			else{
				c[s[i]][s[j]]=cen[w].custo[s[i]][s[j]];
			}
		}
	}
	for(int i=p ; i < n ; i++){
		for (int j=p; j<n; j++){
			if(s[i]==s[j]){
				c[s[i]][s[j]]=0.;
				}else{
					c[s[i]][s[j]]=inf;
				}
			}
		for (int j=0; j<p; j++){
			c[s[i]][s[j]]=cen[w].custo[s[i]][s[j]];
		}	
	}
 	for(int k=0 ; k < p ; k++){
		for (int l=0; l<n; l++){
			d[s[k]][s[l]]=inf;
			for (int m=0; m<p; m++){
				if(d[s[k]][s[l]] > alpha*c[s[k]][s[m]]+c[s[m]][s[l]]){
					d[s[k]][s[l]] = alpha*c[s[k]][s[m]]+c[s[m]][s[l]];
				}
			}
		}
	}
	for(int i=0 ; i < n ; i++){
		for (int j=0; j<n; j++){
			e[s[i]][s[j]]=inf;
			for (int k=0; k<p; k++){
				if(e[s[i]][s[j]] > c[s[i]][s[k]]+d[s[k]][s[j]]){
					e[s[i]][s[j]] = c[s[i]][s[k]]+d[s[k]][s[j]];
				}
			}
		}
	}
	double soma=0;
	for(int i=0 ; i < n ; i++){
		for (int j=0; j<n; j++){
			soma = soma + e[s[i]][s[j]] * (cen[w].fluxo[s[i]][s[j]]);
		}
	}
	libera_matriz_double(c, n);
	libera_matriz_double(d, n);
	libera_matriz_double(e, n);
	return soma;
}

// Objective Function
double FO(int n, int p, double alpha, double beta, int* s, const Cenario cen[], int size){
	double custo_cen[size];
	double pp[size];
	for (int t = 0; t<size; t++){
		custo_cen[t]=floyd_thr(n, p, alpha, s, t, cen, size);
		pp[t]=cen[t].prob;
	}
	double result;
	result = media(custo_cen, pp, size);
	return result;
}


// Local search

double best_neighbor(int n, int p, double alpha, double beta, int *s, const Cenario cen[], int size, double fo, int *melhor_i, int *melhor_j)
{
  int aux;
  double fo_melhor_viz = fo;

  double fo_viz;
  for(int i=0 ; i < p ; i++){
    for(int j=p ; j < n ; j++) {
      
      movimento_troca_hub_naohub(i,j,s);

      // Calculate the new distance
      fo_viz = FO(n, p, alpha, beta, s, cen, size);
      // Store the best move (best swap)
      if(fo_viz < fo_melhor_viz){
        *melhor_i = i;
        *melhor_j = j;
        fo_melhor_viz = fo_viz;
      }

      // Undo the move
      movimento_troca_hub_naohub(i,j,s);
    }
  }
  // Return the distance of the best neighbor
  return fo_melhor_viz;

}//


double best_improvement(int n, int p, double alpha, double beta, int *s, const Cenario cen[], int size)
{
  int aux, melhor_i, melhor_j;
  double fo_viz, fo;
  bool melhorou;

  fo = fo_viz = FO(n, p, alpha, beta, s, cen, size);

  do{
     melhorou = false;
     
     //Select the best neighbor using the city-swap move – explore the entire neighborhood of 
     fo_viz = best_neighbor(n, p, alpha, beta, s, cen, size, fo, &melhor_i, &melhor_j);
     if (fo_viz < fo ){
	 	//Perform the move to the best neighbor      s <- s'
		  movimento_troca_hub_naohub(melhor_j,melhor_i,s);
          
          fo = fo_viz;
          melhorou = true;
     }
  } while (melhorou == true);
  return fo;
}


double first_neighbor(int n, int p, double alpha, double beta, int *s, const Cenario cen[], int size, double fo, int *melhor_i, int *melhor_j)
{
  int aux;
  double fo_melhor_viz = fo;
  bool melhorou = false;
  int *vet;
  vet = cria_vetor(n);
  for (int i=0; i < n; i++){
  	vet[i] = i;
  }
  embaralha_hub_e_naohub(vet, n, p);

  double fo_viz;
  for(int i=0 ; i < p ; i++){
    for(int j=p ; j < n ; j++) {
      
      movimento_troca_hub_naohub(vet[i],vet[j],s); 

      fo_viz = FO(n, p, alpha, beta, s, cen, size);
      // Store the best move (best swap)
      if(fo_viz < fo_melhor_viz){
        *melhor_i = vet[i];
        *melhor_j = vet[j];
        fo_melhor_viz = fo_viz;
        fo_melhor_viz = fo_viz;
        melhorou = true;
      }

      // Undo the move
 		movimento_troca_hub_naohub(vet[i],vet[j],s);	
		if (melhorou) break;
    }
    if (melhorou) break;
  }
  libera_vetor(vet);
  return fo_melhor_viz;
}



double first_improvement(int n, int p, double alpha, double beta, int *s, const Cenario cen[], int size)
{
  int aux, melhor_i, melhor_j;
  double fo_viz, fo;
  bool melhorou;

  fo = fo_viz = FO(n, p, alpha, beta, s, cen, size);

  do{
     melhorou = false;
     
     fo_viz = first_neighbor(n, p, alpha, beta, s, cen, size, fo, &melhor_i, &melhor_j);
     if (fo_viz < fo ){
	 	//Perform the move to the best neighbor      s <- s'
		  movimento_troca_hub_naohub(melhor_j,melhor_i,s);
          
          fo = fo_viz;
          melhorou = true;
     }
  } while (melhorou == true);
  return fo;
}



//ILS with Best Improvement
double ILS(int n, int p, double alpha, double beta, int *s, const Cenario cen[], int size, int ILSmax)
{
	int nivel, iter, MelhorIter;
	int *s_2l, *s_star;
	double fo, fo_2l, fo_star;
	s_star = cria_vetor(n);
	atualiza_vetor(s_star, s, n);
        construcao_gulosa(n, s, cen, size);
	fo_star = FO(n, p, alpha, beta, s_star, cen, size);
	atualiza_vetor(s,s_star,n);
	s_2l = cria_vetor(n);
	fo =  best_improvement(n, p, alpha, beta, s_star, cen, size);
	atualiza_vetor(s,s_star,n);
	iter = MelhorIter = 0;
	nivel = 2;
	int nivelmax;
	if (p == 10){
	    nivelmax = 7;
	}else{
	    nivelmax = p;
	}
	while (iter - MelhorIter < ILSmax && nivel < nivelmax){
	    iter++;
	    // s_2l <- s
	    atualiza_vetor(s_2l,s,n);
	    fo_2l = fo;
	    movimento_ktrocas_hub_naohub(p, n, s_2l, nivel); //Perturbation
	    fo_2l = best_improvement(n, p, alpha, beta, s_2l, cen, size);
	    if (fo_2l < fo){
	          fo = fo_2l;
	          atualiza_vetor(s,s_2l,n);
		  MelhorIter = iter;
		  nivel = 2;
	      }else{
	          nivel++;
	      }
	} 
	libera_vetor(s_star);
	libera_vetor(s_2l);
	return fo;
}

//ILS with First Improvement
double ILS_first(int n, int p, double alpha, double beta, int *s, const Cenario cen[], int size, int ILSmax)
{
	int nivel, iter, MelhorIter;
	int *s_2l, *s_star;
	double fo, fo_2l, fo_star;
	s_star = cria_vetor(n);
	atualiza_vetor(s_star, s, n);
        construcao_gulosa(n, s, cen, size);
	fo_star = FO(n, p, alpha, beta, s_star, cen, size);
	atualiza_vetor(s,s_star,n);
	s_2l = cria_vetor(n);
	fo =  first_improvement(n, p, alpha, beta, s_star, cen, size);
	atualiza_vetor(s,s_star,n);
	iter = MelhorIter = 0;
	nivel = 2;
	int nivelmax;
	if (p == 10){
	    nivelmax = 7;
	}else{
	    nivelmax = p;
	}
	while (iter - MelhorIter < ILSmax && nivel < nivelmax){
	    iter++;
	    // s_2l <- s
	    atualiza_vetor(s_2l,s,n);
	    fo_2l = fo;
	    movimento_ktrocas_hub_naohub(p, n, s_2l, nivel); //perturbation
	    fo_2l = first_improvement(n, p, alpha, beta, s_2l, cen, size);
	    if (fo_2l < fo){
	          fo = fo_2l;
	          atualiza_vetor(s,s_2l,n);
		  MelhorIter = iter;
		  nivel = 2;
	      }else{
	          nivel++;
	      }
	} 
	libera_vetor(s_star);
	libera_vetor(s_2l);
	return fo;
}

// E-ILS with Best Improvement
double EILS_best(int n, int p, double alpha, double beta, int *s, const Cenario cen[], int size, int vezesnivel, int ILSmax)
{
  int nivel, iter, MelhorIter, vezes;
  int *s_2l, *s_star;
  double fo, fo_2l, fo_star;

  s_star = cria_vetor(n);
  inicializa_vetor_int(s_star,n);

  //---SOLUCAO INICIAL---
  
  construcao_gulosa(n, s_star, cen, size);
  fo_star = FO(n, p, alpha, beta, s_star, cen, size);
  atualiza_vetor(s,s_star,n);

  fo =  best_improvement(n, p, alpha, beta, s_star, cen, size);
  atualiza_vetor(s,s_star,n);
  
  iter = MelhorIter = 0;
  nivel = 2;
  vezes = 0;

  s_2l = cria_vetor(n);
  
  while (iter - MelhorIter < ILSmax && nivel < p){
    iter++;
    atualiza_vetor(s_2l,s,n);
    fo_2l = fo;
    movimento_ktrocas_hub_naohub(p, n, s_2l, nivel); //perturbation
    fo_2l = best_improvement(n, p, alpha, beta, s_2l, cen, size);
    if (fo_2l < fo){
        fo = fo_2l;
        atualiza_vetor(s,s_2l,n);
        vezes = 0;
        nivel = 2;
        MelhorIter = iter;
    }else{
        if(vezes < vezesnivel){
            vezes++;
        }else{
            nivel++;
            vezes = 0;
        }
    }
 } // end while

  libera_vetor(s_star);
  libera_vetor(s_2l);
  return fo;
}

// E-ILS with First Improvement
double EILS_first(int n, int p, double alpha, double beta, int *s, const Cenario cen[], int size, int vezesnivel, int ILSmax)
{
  int nivel, iter, MelhorIter, vezes;
  int *s_2l, *s_star;
  double fo, fo_2l, fo_star;

  s_star = cria_vetor(n);
  inicializa_vetor_int(s_star,n);

  //---INITIAL SOLUTION---
  
  construcao_gulosa(n, s_star, cen, size);
  fo_star = FO(n, p, alpha, beta, s_star, cen, size);
  atualiza_vetor(s,s_star,n);

  fo =  first_improvement(n, p, alpha, beta, s_star, cen, size);
  atualiza_vetor(s,s_star,n);
  
  iter = MelhorIter = 0;
  nivel = 2;
  vezes = 0;

  s_2l = cria_vetor(n);
  
  while (iter - MelhorIter < ILSmax && nivel < p){
    iter++;
    atualiza_vetor(s_2l,s,n);
    fo_2l = fo;
    movimento_ktrocas_hub_naohub(p, n, s_2l, nivel); //perturbacao
    fo_2l = first_improvement(n, p, alpha, beta, s_2l, cen, size);
    if (fo_2l < fo){
        fo = fo_2l;
        atualiza_vetor(s,s_2l,n);
        vezes = 0;
        nivel = 2;
        MelhorIter = iter;
    }else{
        if(vezes < vezesnivel){
            vezes++;
        }else{
            nivel++;
            vezes = 0;
        }
    }
 } // end while

  libera_vetor(s_star);
  libera_vetor(s_2l);
  return fo;
}


int main(int argc, char* argv[]) {
	try {
		// ==========================================
		// (Command-line arguments: CAB-High-10.txt p alpha beta)
		// ==========================================
		ifstream arq(argv[1]);
		if (!arq.is_open()) {
			cout << "Error openning file: " << argv[1] << endl;
			arq.close();
			exit(EXIT_FAILURE);
		}
		int n; // Number of nodes
		int p; // Number of hubs to be activated
		int numcen;
		double alpha;
		double beta;
		if (argc >= 2) //Read the number of hubs to be activated
			p = atoi(argv[2]);
		else p = 3;
		alpha = atof(argv[3]);
		beta = atof(argv[4]);
		clock_t inicio_CPU, fim_CPU;  
		inicio_CPU = clock();
		unsigned int seed = static_cast<unsigned int>(time(0));
                std::srand(seed);
		// Instance reading
		arq >> n;
		arq >> numcen;
		double m = numcen;
		int ILSmax = 3;	
		Cenario cen[numcen];
		int *s;
		s = cria_vetor(n);
                inicializa_vetor_int(s, n);
		for (int w = 0; w < numcen; w++){
			cen[w].prob=1/m;
			cen[w].fluxo = cria_matriz_double(n, n);
			cen[w].custo = cria_matriz_double(n, n);
		}
		for (int w = 0; w < numcen; w++){
			for (int i = 0; i < n; i++){
				for (int j = 0; j < n; j++){
					arq >> cen[w].custo[i][j];
				}
			}
			for (int i = 0; i < n; i++){
				for (int j = 0; j < n; j++){
					arq >> cen[w].fluxo[i][j];
				}
			}
		}
		double fotimo;
		double tempo;
		fotimo = ILS(n, p, alpha, beta, s, cen, numcen, ILSmax);
		// if ILS - First Improvement
		// fotimo = ILS_first(n, p, alpha, beta, s, cen, numcen, ILSmax);
		// if E-ILS, add int vezesnivel = ...;
		// if E-ILS - Best Improvement
		// fotimo = EILS_best(n, p, alpha, beta, s, cen, numcen, ILSmax);
		// if E-ILS - First Improvement
		// fotimo = EILS_first(n, p, alpha, beta, s, cen, numcen, ILSmax);
		fim_CPU = clock();
		tempo=double(fim_CPU - inicio_CPU)/(CLOCKS_PER_SEC);
		FILE *re;
		re = fopen("Results.txt","aw+");
		fprintf(re, "\n %10.4f \t \t %10.4f \t \t", tempo, fotimo);
		for(int i = 0; i < p; i++){
			fprintf(re, "%d \t", s[i]+1);
		}
		fprintf(re, "\t %d", seed);
	}
	catch(int erro) {
		cerr << "Error: " << endl;
	}
	return 0;
}
