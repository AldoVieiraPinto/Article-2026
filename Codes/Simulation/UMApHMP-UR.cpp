#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include<iostream>
#include <limits.h>
#include <vector>
#include <time.h>
#include <string>
#include <chrono>
#include <thread>
#include <random>


using namespace std;

#define inf 10e10 

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

void insertionSorte( double* v, double *w, int n ) 
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

/* cria memoria para um vetor de tam posicoes */
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




void imprime_vetor(int *s, int n)
{
    for (int j=0; j < n; j++) 
	   printf("s[%2d]=%d \n",j,s[j]);
}


void inicializa_vetor(int *vetor, int tam)
{
    for (int j=0; j<tam; j++) vetor[j] = 0;
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



double cvar(double* c, double* pp, double beta, int tam ){
		insertionSorte(c, pp, tam);
		double soma = 0.;
		int i = 0;
		int kb;
		while((soma < beta-0.000000000000001)&&(i<tam)){
			kb=i;
			soma = soma + pp[i];
			i = i+1;
		};
		double var = c[kb];
		double cvar = (soma-beta)*var; 
		for (int j = kb+1; j < tam; j++){
			cvar = cvar + pp[j]*c[j];
		}
		cvar = cvar/(1-beta);
	return cvar;
}


double media(double* c, double* pp, int tam){
		double soma = 0;
		for (int t = 0; t < tam; t++){
			soma = soma + pp[t]*c[t];
		}
	return soma;
}

// Cost function -UR: unrestricted routes ;
double FO(int n, int p, double alpha, int* s, double** custo_det, double** fluxo_det){
	double** c = cria_matriz_double(n, n);
 	for(int i=0 ; i < p ; i++){
		for (int j=0; j<n; j++){
			if(j<p){
				c[s[i]][s[j]]=alpha*custo_det[s[i]][s[j]];
			} 
			else{
			c[s[i]][s[j]]=custo_det[s[i]][s[j]];
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
			c[s[i]][s[j]]=custo_det[s[i]][s[j]];
		}	
	}
	for(int k=0 ; k < p ; k++){
		for(int i=0 ; i < n ; i++) {
			for(int j = 0; j < n; j++){
					if(c[s[i]][s[j]] > c[s[i]][s[k]]+c[s[k]][s[j]]){
							c[s[i]][s[j]] = c[s[i]][s[k]]+c[s[k]][s[j]];
				}
			}
		}
	}
	double soma=0;
	for(int i=0 ; i < n ; i++){
		for (int j=0; j<n; j++){
			soma = soma + c[s[i]][s[j]] * fluxo_det[s[i]][s[j]];
		}
	}
	libera_matriz_double(c, n);
	return soma;
}

int main(int argc, char* argv[]) {
	try {
		ifstream arq(argv[1]);
		if (!arq.is_open()) {
			cout << "Error openning file: " << argv[1] << endl;
			arq.close();
			exit(EXIT_FAILURE);
		}
		int n; 
		int p; 
		double alpha;
		double beta;
		double nominal;
		double delta;

		if (argc >= 2) 
			n = atoi(argv[2]);
		else n = 20;
		p = atoi(argv[3]);
		alpha = atof(argv[4]);
		beta = atof(argv[5]);
		delta = atof(argv[6]);
		nominal = atof(argv[7]);
		double** custo_det;
		double** fluxo_det;
		custo_det = cria_matriz_double(n, n);
		fluxo_det = cria_matriz_double(n, n);
		int *s_auxi;
		int *s;
		int *h;
		s_auxi = cria_vetor(p);
		s = cria_vetor(n);
		h = cria_vetor(n);
		for (int i = 0; i < p; i++){
		      s_auxi[i] = atof(argv[i + 8]);
		      s_auxi[i] = s_auxi[i]-1;
		}
		for (int i = 0; i < n; i++){
		      h[i] = 0;
		}
		for (int i = 0; i < p; i++){
		      h[s_auxi[i]] = 1;
		}
		int conth = 0;
		int conthh = n-1;
		for (int j = 0; j < n; j++) {
			if (h[j] >= 0.1) {
				s[conth] = j;
				conth++;
			}else{
			        s[conthh] = j;
			        conthh--;
			}
		}
		for (int i = 0; i < n; i++){
		      for (int j = 0; j < n; j++){
			   arq >> custo_det[i][j];
		      }
		}
		for (int i = 0; i < n; i++){
		      for (int j = 0; j < n; j++){
			   arq >> fluxo_det[i][j];
		      }
		}
	        int tam = 10000;
		double prob = 0.0001;
		double pp[tam];
		double media_sim;
		double cvar_sim;
		double mcvar_sim;
		
		for (int m = 0; m < tam; m++) {
		      pp[m]=prob;
		}

		std::random_device rd;
		std::mt19937 gen(rd());

		// Define o intervalo desejado 
		double min = 1-delta;
		double max = 1+delta;
		
		std::uniform_real_distribution<double> dis(min, max);
                double c[tam];
		for (int m = 0; m < tam; m++) {
		        double** catual;
		        double** watual;
		        catual = cria_matriz_double(n, n);
		        watual = cria_matriz_double(n, n);
			for (int i = 0; i < n; i++) { 
				for (int j = 0; j < n; j++) {
					catual[i][j] = dis(gen) * custo_det[i][j];
				}
			}
		        for (int i = 0; i < n; i++) { 
			          for (int j = 0; j < n; j++) {
				          watual[i][j] = dis(gen) * fluxo_det[i][j];
			          }
		        }
		        c[m] = FO(n, p, alpha, s, catual, watual); 
                        libera_matriz_double(catual, n);
		        libera_matriz_double(watual, n);
                }
                
                media_sim = media(c, pp, tam);
                cvar_sim = cvar(c, pp, beta, tam);
                mcvar_sim = 0.5*media_sim+0.5*cvar_sim;

                FILE *re;
		re = fopen("Resultados.txt","aw+");
		fprintf(re, "\n %d \t %10.4f \t %10.4f \t %10.4f \t%10.4f \t%10.4f \n", p, alpha, beta, media_sim, cvar_sim, mcvar_sim);
 	}
	catch(int erro) {
		cerr << "Error: " << endl;
	}
	return 0;
}



