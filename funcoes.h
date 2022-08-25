#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 2		//tamanho da rede
#define N L*L		//número de spins
#define Nskip 500000	//número de passos para o sistema atingir o equilíbro
int Rede[L][L];		//rede quadrada

double Rand(unsigned int *seed);

void inicia_malha(unsigned int *seed);

void escolhe_pos(int *p, unsigned int *seed);

int hamiltoniana(int *p, int J, int B);

int testa_flip(int *p, int *dE, int J, int B, int k, double T, unsigned int *seed);

void equilibra(int J, int B, int k, double T, unsigned int *seed);

int magnetizacao_total();

int energia_total(int J, int B);
