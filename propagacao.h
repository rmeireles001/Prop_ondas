#ifndef aco_h
#define aco_h

#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <fstream>
#include "parameters.h"


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
#define ponderacao 1.0
#define ro 2700
#define pontos 1001
#define contagem 1000


using namespace std;

class propagacao{
public:
	double *Gexp, *G, *A, *Z, *R, *F, *P, *posicao;
	int *config;
	double **A1, **B1, **Y;
	double E, c;
	propagacao();
	~propagacao();
	double **alocar_mat(int size);
	void liberar_mat(double **mat, int size);
	void inserir(string s);
	void prob_direto(int n, double *G);
	void prob_inverso();
	void prob_inverso(int n);
	double erroG(int i);
	void atribuirA(int n, double v);
	void dados_experimentais(string s, int n, int ini);
	void estimativa_inicial(double n);
	void config_area(double n);
	void attr_config(int k, double var);
	void atualizar_area(int k);
	void zerar(int n);
	double luus_jaakola(int pos);
	void run_lj(string arq_areas, string saida, int inicio, int fim);
	//void run_aco(string arq_areas, string saida, int inicio, int fim);
	void escrever_txt(string saida, string metodo, double tempo, int inicio, int fim);
	double cgrasp(double llimit, double ulimit, double incr, int section);
	double time_cgrasp(string arq_areas, string saida, int inicio, int fim);
	void run_cgrasp(FILE *resposta, int inicio, int fim);
};


class variable{
public:
	double ulimit;
	double llimit;
	double increment;
	int nvalues;
	double *trail, *prob;
	void init(double ulimit, double llimit, double increment);
	void free();
	void find_nvalue();
	double *allocate(int size);
	void reset(double *vetor, int size);
	void ls_nvalue();
};

class values{
public:
	double var, ofn;
	int iter;
	values();
	void copy_from(values *from);
};

class ants{
public:
	int variable;
	double ofn;
	double bestofn;
	int iter;
	ants();
	void init(int itno, double iofn);
	void copy_from(ants *from);
};

struct trash{
	double *help_b;
	values tempval;
	trash();
	~trash();
	void reset(double *vec, int size);
};

class aco : public propagacao{
public:
	variable var;
	values itbestval, itworstval, bestval, *runbestval;
	ants *ant, bestant, itbestant, itworstant, *runbestant; 
	long int nevaluations, *nevalbest;
	trash t;
	long int seed;
	void get_data(double llimit, double ulimit, double increment);
	void end();
	void initialize_ants_variables(int runno);
	void initialize_trail();
	void iteration_init(int iterationno, int nants);
	void find_values(int nants);
	int find_var_valueno();
	int find_no_of_max(int *noofmax);
	void find_prob();
	void decode_varvalues(double *varvalues, int antno);
	void analysis(int iterationno, int nants, int section);
	//Objective function
	double objective_function(double *varvalue, int section);
	void some_stats(int iterationno, int runno);
	void evaporate_trail();
	void update_trail_weighted(ants *antptr, int weight);
	void ras_update(double *help_b);
	void trail();
	void run(int section);
	double get_var();
	double time_aco(string arq_areas, string saida, int inicio, int fim);
	void run_aco(FILE *resposta, int inicio, int fim);
};

class lsaco : public aco{
public:
	void get_data(values *gvalptr, variable *varv);
	void run(values *gvalptr, variable *varv, int section, long int *s);
	void some_stats(int iterationno);
	void evaporate_trail();
	void update_trail(ants *antptr);
	void trail();
	void end();
};


class monomio{
private:
	double coef;
	int exp;
	monomio *prev, *next;
public:
	monomio(double coef, int exp);
	monomio(string s);
	int getExp();
	double getCoef();
	void setExp(int exp);
	void setCoef(double coef);
	monomio* getNext();
	monomio* getPrev();
	void setNext(monomio *m);
	void setPrev(monomio *m);
	void str2monomio(string s);
	int isValid(char s);
	monomio* newMonomio(string s);
	int isNum(char s);
	void derivar();
	string monomio2str();
};

class polinomio{
private:
	monomio* cabeca;
	monomio *cauda;
	string s;
public:
	polinomio();
	polinomio(double coef, int exp);
	polinomio(string s);
	int empty();
	void insert(double coef, int exp);
	void insert(string s);
	void percorrer();
	void percorrer_inv();
	void removerCoef(double coef);
	monomio* getCabeca();
	monomio* getCauda();
	void setPtrs(monomio *cabeca, monomio *cauda);
	void add_polinomio(const polinomio& p);
	polinomio operator+(const polinomio& p);
	polinomio operator-(const polinomio& p);
	polinomio operator*(const polinomio& p);
	//void operator=(const polinomio &p);
	void self_mult(string s);
	void atualizar_string();
	friend ostream &operator<<(ostream &output, const polinomio &p){
		output << p.s;
		return output;
	}
	void escalar(double n);
	double p(double x);
	void insert_string(string s);
	void derivar();
	void zerar();
};

class interpoLagrange{
private:
	int grau;
	double *pts;
	double *y;
	polinomio **L;
public:
	interpoLagrange(int grau);
	~interpoLagrange();
	void inserir_pts(int p, double pto, double valor);
	string binomio(double num);
	void insertLs();
	double objective_function(double x);
	polinomio *interpolar();
};


double drand(double low, double high);
double ran01( long *idum );
double function(double x);
string int2str(int num);

#endif
