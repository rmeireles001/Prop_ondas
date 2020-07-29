#include "propagacao.h"

propagacao::propagacao(){
	Gexp = new double[pontos];
	G = new double[pontos];
	A = new double[pontos];
	Z = new double[pontos];
	R = new double[pontos];
	F = new double[pontos];
	P = new double[pontos];
	config = new int[pontos];
	posicao = new double[pontos];
	A1 = alocar_mat(pontos);
	B1 = alocar_mat(pontos);
	Y = alocar_mat(pontos);
	for(int i=0; i<pontos; i++){
		Gexp[i] = 0.0;
		G[i] = 0.0;
		A[i] = 0.0;
		Z[i] = 0.0;
		R[i] = 0.0;
		F[i] = 0.0;
		P[i] = 0.0;
		posicao[i] = 0.0;
		config[i] = 0.0;
	}
	for(int i=0; i<pontos; i++){
		for(int j=0; j<pontos; j++){
			A1[i][j] = 0.0;
			B1[i][j] = 0.0;
			Y[i][j] = 0.0;
		}
	}
	E = 71*pow(10, 9);
	c = sqrt(E/ro);
	cout << "Alocados\n";
}

propagacao::~propagacao(){
	delete Gexp;
	delete G;
	delete A;
	delete Z;
	delete R;
	delete F;
	delete P;
	delete posicao;
	delete config;
	liberar_mat(A1, pontos);
	liberar_mat(B1, pontos);
	liberar_mat(Y, pontos);
	cout << "liberados\n";
}

double **propagacao::alocar_mat(int size){
	double **nv = new double*[size];
	for(int i=0; i<size; i++){
		nv[i] = new double[size];
	}
	return nv;
}

void propagacao::inserir(string s){
	int k;
	ifstream mf(s.c_str());
	for(int i=1; i<pontos; i++){
		mf >> posicao[i] >> A[i]; 
	}
	posicao[0] = 0;
	A[0] = A[1];
	/*for(int i=0, k=480; i<100; i++, k++){
		A[i] = A[k];
		posicao[i] = posicao[k];
	}*/
}

void propagacao::zerar(int n){
	for(int i=0; i<=n; i++){
		G[i] = 0.0;
		Z[i] = 0.0;
		R[i] = 0.0;
		F[i] = 0.0;
		P[i] = 0.0;
		A[i] = 0.0;
		config[i] = 0;
	}
	for(int i=0; i<=n; i++){
		for(int j=0; j<=n; j++){
			A1[i][j] = 0.0;
			B1[i][j] = 0.0;
			Y[i][j] = 0.0;
		}
	}
}

void propagacao::estimativa_inicial(double n){
	for(int i=0; i<=n; i++){
		A[i] = 300.0;
	}
}

void propagacao::config_area(double n){
	for(int i=0; i<=n; i++){
		config[i] = A[i]/A[0];
	}
}

void propagacao::attr_config(int k, double var){
	config[k] = var;
}

void propagacao::atualizar_area(int k){
	A[k] = config[k]*300;
}

void propagacao::dados_experimentais(string s, int n, int ini){
	int i, j;
	double *aux1, *aux2;
	double aux;
	ifstream mf(s.c_str());
	aux1 = new double[pontos];
	aux2 = new double[pontos];
	for(i=1; i<pontos; i++){
		mf >> aux1[i] >> aux >> aux2[i];
	}
	mf.close();
	for(i=ini, j=3; i<=(ini+n-3); i++, j++){
		Gexp[j] = aux2[i];
		posicao[j] = aux1[i];
	}
	Gexp[1] = aux2[ini-2];
	Gexp[2] = aux2[ini-1];
	posicao[1] = aux1[ini-2];
	posicao[2] = aux1[ini-1];
	delete aux1;
	delete aux2;
}

void propagacao::atribuirA(int n, double v){
	A[n] = v;
}

void propagacao::liberar_mat(double **mat, int size){
	for(int i=0; i<size; i++){
		delete mat[i];
	}
	delete mat;
}

void propagacao::prob_inverso(){
	int n1;
	F[1]=1.0;

	for(int i=0; i<=2; i++){
		Z[i]=300*A[i]*ro*c;
	}

	for(int i=1; i<=2; i++){
		if((Z[i]+Z[i-1])!=0.0)
			R[i] = (Z[i]-Z[i-1])/(Z[i]+Z[i-1]);
	}
	P[1] = 0.0;
	P[2] = 0.0;

	n1 = 1;
	G[1] = ((G[1] + (R[n1+1]+P[n1])*F[1-n1+1]));

	for(n1=1; n1<=2; n1++){
		G[2]=((G[2]+(R[n1+1]+P[n1])*F[2-n1+1]));
	}
}

void propagacao::prob_inverso(int n){
	double somatorio1;
	int n1, p1, k;

	Z[n]=300*A[n]*ro*c;

	if((Z[n]+Z[n-1])!=0.0)
		R[n] = (Z[n]-Z[n-1])/(Z[n]+Z[n-1]);
	for(p1=1; p1<=n-2; p1++){
		if(p1<=(n-3))
			A1[p1][n]=0;
		A1[n-2][n-1]=0;
		B1[p1][n-1] = 0;
		Y[p1][n] = 0;
		P[n] = 0;
	}

	for(p1=1; p1<=n-2; p1++){
		if(p1<=(n-3))
			A1[p1][n]=A1[p1][n-2]-B1[p1][n-2];
		A1[n-2][n-1]=0;
		somatorio1=0.0;
		for(k=1; k<=p1-1; k++)
			somatorio1=somatorio1+Y[k][n-1];
		B1[p1][n-1] = (R[n-p1-1]*(somatorio1+R[n-1]));
		Y[p1][n] = (R[n-p1]*(A1[p1][n-1] - B1[p1][n-1]));
		P[n] = (P[n]+Y[p1][n]);
	}

	G[n] = 0;
	for(n1=1; n1<=n; n1++){
		G[n]=((G[n]+(R[n1]+P[n1]))*F[n-n1+1]);
	}
}

void propagacao::prob_direto(int n, double *G){
	double somatorio1;
	int n1, p1, k;

	zerar(n);

	F[1]=1.0;

	for(int i=0; i<=n; i++){
		Z[i]=300*A[i]*ro*c;
	}

	for(int i=1; i<=n; i++){
		if((Z[i]+Z[i-1])!=0.0)
			R[i] = (Z[i]-Z[i-1])/(Z[i]+Z[i-1]);
	}

	P[1] = 0.0;
	P[2] = 0.0;

	n1 = 1;
	G[1] = ((G[1] + (R[n1+1]+P[n1])*F[1-n1+1]));

	for(n1=1; n1<=2; n1++){
		G[2]=((G[2]+(R[n1+1]+P[n1])*F[2-n1+1]));
	}

	for(n1=3; n1<=n; n1++){
		for(p1=1; p1<=n1-2; p1++){
			if(p1<=(n1-3))
				A1[p1][n1]=A1[p1][n1-2]-B1[p1][n1-2];
			A1[n1-2][n1-1]=0;
			somatorio1=0.0;
			for(k=1; k<=p1-1; k++)
				somatorio1=somatorio1+Y[k][n1-1];
			B1[p1][n1-1] = (R[n1-p1-1]*(somatorio1+R[n1-1]));
			Y[p1][n1] = (R[n1-p1]*(A1[p1][n1-1] - B1[p1][n1-1]));
			P[n1] = (P[n1]+Y[p1][n1]);
		}
	}
	for(int j=3; j<=n; j++){
		for(n1=1; n1<=j; n1++){
			G[j]=((G[j]+(R[n1]+P[n1]))*F[j-n1+1]);
		}
	}
}

double propagacao::erroG(int i){
	return sqrt(pow(G[i]-Gexp[i], 2));//pow(A[i]-0.5, 2)+2;//sqrt(pow(G[i]-Gexp[i], 2));
}

double func(double x){
	return pow(x-2, 2)+2;
}

double propagacao::luus_jaakola(int pos){
	double mini, maxi, Rr, r, eps, aux1, aux2, aux1ant, qbest, t1, condicao, randomico, newconfig, oldconfig;
	int n_in, i, j, k, it, n_out, PAROU=0, aux, atrib, x, custo=0;
	Rr = 0.0;
	//G[] = 0.0;
	n_out = 1;
	n_in = 1;
	eps = 0.05;

	mini = 0.0;
	maxi = 1.0;

	r = maxi-mini;

	randomico = drand(0, 1);
	oldconfig = mini + r*randomico;

	aux1 = 1000;
	qbest = 1000;
	i = 1.0;
	condicao = pow(10, -5);
	while(qbest > condicao && PAROU<500){
		randomico = drand(0, 1);
		Rr = ((-0.5 + randomico)*(r));
		newconfig = oldconfig + Rr;
		if(newconfig < mini){
			newconfig = mini;
		}
		if(newconfig > maxi){
			newconfig = maxi;
		}
		atribuirA(pos, newconfig);
		prob_inverso(pos);
		aux1 = erroG(pos);
		custo++;
		atribuirA(pos, oldconfig);
		prob_inverso(pos);
		aux2 = erroG(pos);
		custo++;
		if(aux1<=aux2){
			qbest = aux1;
			oldconfig = newconfig;
		}

		i++;
		PAROU++;
		r = ((1 - eps)*r);
	}
	config[pos] = custo;
	return oldconfig;
}

void propagacao::run_lj(string arq_areas, string saida, int inicio, int fim){
	clock_t start, end;
	double img;
	inserir(arq_areas);
	prob_direto(1000, Gexp);
	prob_direto(inicio, G);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		img = luus_jaakola(i);
		atribuirA(i, img);
		prob_inverso(i);
	}
	end = clock();
	escrever_txt(saida, "LUUS JAAKOLA", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
}

/*void propagacao::run_aco(string arq_areas, string saida, int inicio, int fim){
	clock_t start, end;
	inserir(arq_areas);
	prob_direto(1000, Gexp);
	prob_direto(inicio, G);
	aco a;
	a.get_data(0,1,0.05);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		a.run(i, this);
		atribuirA(i, a.get_var());
		prob_inverso(i);
	}
	end = clock();
	escrever_txt(saida, "ANT COLONY OPTIMIZATION", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
}*/

//(double)(end-start)/(double)(CLOCKS_PER_SEC)
void propagacao::escrever_txt(string saida, string metodo, double tempo, int inicio, int fim){
	FILE *resultados = fopen(saida.c_str(), "w");
	fprintf(resultados, "Método usado: %s\n\nPOSIÇÃO\tÁREA(mm²)\tCUSTO\n\n", metodo.c_str());
	for(int i=inicio; i<=fim; i++){
		fprintf(resultados, "%d\t%.15e\t%d\n", (int) posicao[i], A[i], (int) config[i]);
	}
	fprintf(resultados, "\n\n");
	fprintf(resultados, "Tempo gasto: %lf\n", tempo);
	fclose(resultados);
}

double function(double x){
	return pow(x-1.5, 2)+3;
}

double propagacao::cgrasp(double llimit, double ulimit, double incr, int section){
	double fbest = 1e10;
	double xbest = 0;
	double xi, xf;
	double f;
	double max=ulimit+0.1;
	int divisor = -2;
	//Busca linear / construção
	for(double x=llimit; x<=max && fbest>1e-5; x+=incr){
		atribuirA(section, x);
		prob_inverso(section);
		f = erroG(section);
		config[section]++;
		if(f<fbest){
			fbest = f;
			xbest = x;
		}
	}
	while(divisor>=-5){
		incr = pow(10, divisor);
		xi = xbest - incr;
		if(xi < llimit)
			xi = llimit;
		xf = xbest + incr;
		if(xf > ulimit)
			xf = ulimit;
		for(double x=xi; x<=xf && fbest>1e-5; x+=incr){
			atribuirA(section, x);
			prob_inverso(section);
			f = erroG(section);
			config[section]++;
			if(f<fbest){
				fbest = f;
				xbest = x;
			}
		}
		divisor--;
	}
	
	return xbest;
}

double propagacao::time_cgrasp(string arq_areas, string saida, int inicio, int fim){
	clock_t start, end;
	double img;
	inserir(arq_areas);
	prob_direto(1000, Gexp);
	prob_direto(inicio, G);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		img = cgrasp(0, 1, 0.05, i);
		atribuirA(i, img);
		prob_inverso(i);
	}
	end = clock();
	escrever_txt(saida, "C-GRASP", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
	zerar(fim);
	return (double)(end-start)/(double)(CLOCKS_PER_SEC);
}

void propagacao::run_cgrasp(FILE *resposta, int inicio, int fim){
	double med=0, var=0;
	double tempos[contagem];
	fprintf(resposta, "C-GRASP ÁREAS %d - %d\n\n", inicio, fim);
	for(int i=0; i<contagem; i++){
		tempos[i] = time_cgrasp("areas.txt", "cgrasp_run#"+int2str(i+1)+".txt", inicio, fim);
		fprintf(resposta, "%lf\t", tempos[i]);
		cout << tempos[i] << endl;
	}
	for(int i=0; i<contagem; i++){
		med = med + tempos[i];
	}
	med=med/contagem;
	for(int i=0; i<contagem; i++){
		var = var + pow(tempos[i] - med, 2);
	}
	var = var/(contagem);
	fprintf(resposta, "\n\nMédia: %lf\tVariância: %lf\tDesvio-padrão: %lf\n\n\n", med, var, sqrt(var));
}

void variable::init(double llimit, double ulimit, double increment){
	this->ulimit = ulimit;
	this->llimit = llimit;
	this->increment = increment;
	find_nvalue();
	trail = allocate(nvalues);
	prob = allocate(nvalues);
}

void variable::free(){
	delete trail;
	delete prob;
}

void variable::find_nvalue(){
	double nv, quotient;
	int divisor;
	nv = (ulimit-llimit)/increment;
	divisor = (int) nv;
	nvalues = 1+divisor;
	quotient = nv/((double)divisor);
	if(quotient != 1.0f){
		if((fabs((quotient - 1.000000000000000) ) > 0.01)){
			nvalues++;
		}
	}
}

double *variable::allocate(int size){
	double *vet = new double[size];
	if(vet == NULL){
		cout << "Memória insuficiente" << endl;
	}
	reset(vet, size);
	return vet;
}

void variable::reset(double *vetor, int size){
	for(int i=0; i<size; i++){
		vetor[i] = 0;
	}
}

void variable::ls_nvalue(){
	
	find_nvalue();
	trail = allocate(nvalues);
	prob = allocate(nvalues);
}

values::values(){
	var = 0;
}

void values::copy_from(values *from){
	ofn = from->ofn;
	iter = from->iter;
	var = from->var;
}

ants::ants(){
	
}

void ants::init(int itno, double iofn){
	this->ofn = iofn;
	this->iter = itno;
	this->bestofn = 1e25;
}

void ants::copy_from(ants *from){
	ofn = from->ofn;
	bestofn = from->bestofn;
	iter = from->iter;
	variable = from->variable;
}

void aco::get_data(double llimit, double ulimit, double increment){
	ant = new ants[ANTS];
	runbestant = new ants[runs];
	nevalbest = new long int[runs];
	runbestval = new values[runs];
	nevaluations=12;
	var.init(llimit, ulimit, increment);
	long int seed = time(NULL);
}

void aco::end(){
	cout << "Destruindo geral" << endl;
	delete ant;
	delete runbestant;
	delete nevalbest;
	delete runbestval;
	var.free();
}

void aco::initialize_ants_variables(int runno){
	nevaluations = 0;
	nevalbest[runno] = 0;
	for(int antno=0; antno<ANTS; antno++){
		ant[antno].init(0, 0.0);
	}
	bestant.init(0, 1e25);
	bestval.ofn=1e25;
}

void aco::iteration_init(int iterationno, int nants){
	for(int i=0; i<nants; i++){
		ant[i].iter = iterationno;
	}
	itbestant.init(iterationno, 1e25);
	itbestval.ofn = 1e25;
	itworstant.init(iterationno, 0.0);
}

//Trail
void aco::initialize_trail(){
	for(int j=0; j<var.nvalues; j++){
		var.trail[j] = tau0;
	}
}

void aco::evaporate_trail(){
	for(int j=0; j<var.nvalues; j++){
		var.trail[j] *= (1-rho);
	}
}

void aco::ras_update(double *help_b){
	int i, antno, target;
	double b;
	for(antno=0; antno<ANTS; antno++)
		help_b[antno] = ant[antno].ofn;
	for(i=0; i<ras_ranks-1; i++){
		b = help_b[0];
		target = 0;
		for(antno=0; antno<ANTS; antno++){
			if(help_b[antno]<b){
				b = help_b[antno];
				target = antno;
			}
		}
		help_b[target] = 1e25;
		update_trail_weighted(&ant[target], ras_ranks-1-i);
	}
	update_trail_weighted(&bestant, ras_ranks);
}

void aco::update_trail_weighted(ants *antptr, int weight){
	int vno, valueno;
	valueno=antptr->variable;
	if(antptr->ofn != 0.0){
		var.trail[valueno] += weight;
	}
	else{
		var.trail[valueno] += weight;
	}
}

void aco::trail(){
	double *help_b;
	evaporate_trail();
	t.reset(t.help_b, ANTS);
	ras_update(t.help_b);
}

void aco::find_values(int nants){
	for(int antno=0; antno<nants; antno++){
		ant[antno].variable = find_var_valueno();
	}
}

//find values for ants
int aco::find_var_valueno(){
	double q, sum, ranno;
	int valueno, nmax;
	find_prob();
	/*q = ran01(&seed);
	ranno = ran01(&seed);*/
	q = drand(0,1);
	ranno = drand(0,1);
	valueno = find_no_of_max(&nmax);

	if((q<q0) && (nmax==1)){
		#if(localupdate==1)
		var.trail[valueno] *= (1-gamma);
		#endif
		return valueno;
	}
	else{
		sum=0.0;
		for(valueno=0; valueno-1<var.nvalues && sum < ranno; valueno++){
			sum += var.prob[valueno];
		}
		#if(localupdate==1)
		if(valueno > 0 && valueno-1 < var.nvalues){
			var.trail[valueno-1] *= (1-rho);
		}
		#endif
		return (valueno-1);
	}
}

int aco::find_no_of_max(int *noofmax){
	int valueno, chooseno, wentinside=0, inloop=0;
	*noofmax=var.nvalues;
	double max=0;
	for(valueno=0; valueno<var.nvalues; valueno++){
		if(max < var.prob[valueno]){
			max = var.prob[valueno];
			chooseno = valueno;
		}
	}
	for(valueno=0; valueno<var.nvalues; valueno++){
		if(max > var.prob[valueno]){
			wentinside++;
			(*noofmax) = (*noofmax) - 1;
		}
	}
	return chooseno;
}

void aco::find_prob(){
	int valueno;
	double denominator=0.0;
	for(valueno=0; valueno<var.nvalues; valueno++){
		denominator += pow((double) var.trail[valueno], (double) alpha);
	}
	for(valueno=0; valueno<var.nvalues; valueno++){
		var.prob[valueno] = pow((double) var.trail[valueno], (double) alpha)/denominator;
	}
}

//Analysis
void aco::decode_varvalues(double *varvalues, int antno){
	int vno, valueno;
	if(ant[antno].variable == (var.nvalues)-1)
		*varvalues = var.ulimit;
	else
		*varvalues = (var.llimit+(ant[antno].variable*var.increment));
}

void aco::analysis(int iterationno, int nants, int section){
	int antno, vno, valueno;
	t.tempval.iter = iterationno;
	for(antno=0; antno<nants; antno++){
		ant[antno].iter = iterationno;
		decode_varvalues(&t.tempval.var, antno);
		ant[antno].ofn = objective_function(&t.tempval.var, section);
		t.tempval.ofn = ant[antno].ofn;

		if(ant[antno].ofn < ant[antno].bestofn){
			ant[antno].bestofn = ant[antno].ofn;
			ant[antno].iter = iterationno;
		}

		if(ant[antno].ofn<itbestant.ofn){
			itbestval.copy_from(&t.tempval);
			itbestant.copy_from(&ant[antno]);
		}
		if(ant[antno].ofn>itworstant.ofn){
			itworstval.copy_from(&t.tempval);
			itworstant.copy_from(&ant[antno]);
		}
	}
}

double aco::objective_function(double *varvalue, int section){
	atribuirA(section, *varvalue);
	prob_inverso(section);
	config[section]++;
	return erroG(section);
}

void aco::some_stats(int iterationno, int runno){
	if(itbestant.ofn<bestant.ofn){
		if(itbestval.ofn<bestval.ofn){
			bestant.copy_from(&itbestant);
			bestval.copy_from(&itbestval);
			nevalbest[runno] = nevaluations;
		}
	}
}

void aco::run(int section){
	#if localsearch
	lsaco *ls;
	#endif
	int runno, itno;
	double condicao = pow(10, -5);
	for(runno=0; runno<runs; runno++){
		initialize_ants_variables(runno);
		initialize_trail();
		for(itno=0; itno<ncmax && bestval.ofn > condicao; itno++){
			iteration_init(itno, ANTS);
			find_values(ANTS);
			analysis(itno, ANTS, section);
			some_stats(itno, runno);
			#if localsearch
			ls = new lsaco;
			ls->run(&bestval, &var, section, &seed);
			#endif
			trail();
		}
		some_stats(itno, runno);
	}
}

double aco::get_var(){
	return bestval.var;
}

double aco::time_aco(string arq_areas, string saida, int inicio, int fim){
	clock_t start, end;
	inserir(arq_areas);
	prob_direto(1000, Gexp);
	prob_direto(inicio, G);
	get_data(0, 1, 0.05);
	start = clock();
	for(int i=inicio; i<=fim; i++){
		this->run(i);
		atribuirA(i, get_var());
		prob_inverso(i);
	}
	end = clock();
	escrever_txt(saida, "ANTCOLONY OPTIMIZATION", (double)(end-start)/(double)(CLOCKS_PER_SEC), inicio, fim);
	zerar(fim);
	this->end();
	return (double)(end-start)/(double)(CLOCKS_PER_SEC);
}

void aco::run_aco(FILE *resposta, int inicio, int fim){
	double med=0, var=0;
	double tempos[contagem];
	fprintf(resposta, "ANTCOLONY OPTIMIZATION ÁREAS %d - %d\n\n", inicio, fim);
	for(int i=0; i<contagem; i++){
		tempos[i] = time_aco("areas.txt", "aco_run#"+int2str(i+1)+".txt", inicio, fim);
		fprintf(resposta, "%lf\t", tempos[i]);
		cout << tempos[i] << endl;
	}
	for(int i=0; i<contagem; i++){
		med = med + tempos[i];
	}
	med=med/contagem;
	for(int i=0; i<contagem; i++){
		var = var + pow(tempos[i] - med, 2);
	}
	var = var/(contagem);
	fprintf(resposta, "\n\nMédia: %lf\tVariância: %lf\tDesvio-padrão: %lf\n\n\n", med, var, sqrt(var));
}

void lsaco::get_data(values *gvalptr, variable *varv){
	double ll, ul, inc;
	ul = gvalptr->var + (lsnvalues * (varv->increment));
	if(ul > varv->ulimit){
		ul = varv->ulimit;
	}
	ll = gvalptr->var - (lsnvalues * (varv->increment));
	if(ll < varv->llimit){
		ll = varv->llimit;
	}
	inc = varv->increment;
	var.init(ll, ul, inc);
	ant = (ants *) calloc(LSANTS, sizeof(ants));
	for(int antno=0; antno<LSANTS; antno++){
		ant[antno].init(0,0);
	}
	bestant.init(0, 1e10);
}

void lsaco::run(values *gvalptr, variable *varv, int section, long int *s){
	seed = *s;
	get_data(gvalptr, varv);
	initialize_trail();
	for(int itno=0; itno<lsncmax; itno++){
		iteration_init(itno, LSANTS);
		find_values(LSANTS);
		analysis(itno, LSANTS, section);
		some_stats(itno);
		trail();
	}
	if(bestant.ofn<gvalptr->ofn){
		gvalptr->copy_from(&bestval);
	}
	end();
}

void lsaco::some_stats(int iterationno){
	if(itbestant.ofn<bestant.ofn){
		bestant.copy_from(&itbestant);
		bestval.copy_from(&itbestval);
	}
}

void lsaco::evaporate_trail(){
	for(int j=0; j<var.nvalues; j++){
		var.trail[j] *= (1-lsrho);
	}
}

void lsaco::update_trail(ants *antptr){
	int vno, valueno;
	valueno=antptr->variable;
	if(antptr->ofn != 0.0){
		var.trail[valueno] += 1;
	}
	else{
		var.trail[valueno] += 1;
	}
}

void lsaco::trail(){
	evaporate_trail();
	for(int antno=0; antno<LSANTS; antno++){
		update_trail(&ant[antno]);
	//cout << "marcador" << endl;
	}
	update_trail(&itbestant);
}

void lsaco::end(){
	var.free();
}

trash::trash(){
	help_b = new double[ANTS];
}

trash::~trash(){
	delete help_b;
}

void trash::reset(double *vec, int size){
	for(int i=0; i<size; i++){
		vec[i] = 0;
	}
}
double drand(double low, double high)
{
	//srand(time(NULL));
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

double ran01( long *idum )
/*    
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

string int2str(int num){
	string ret="";
	int bit;
	while(num>0){
		bit = num%10;
		num = (num-bit)/10;
		ret = (char)(bit+48)+ret;
	}
	return ret;
}



monomio::monomio(double coef, int exp){
	this->coef = coef;
	this->exp = exp;
	this->next = NULL;
	this->prev = NULL;
}

monomio::monomio(string s){
	this->next = NULL;
	this->prev = NULL;
	str2monomio(s);
}

int monomio::getExp(){
	return this->exp;
}
double monomio::getCoef(){
	return this->coef;
}

void monomio::setExp(int exp){
	this->exp = exp;
}

void monomio::setCoef(double coef){
	this->coef = coef;
}

monomio* monomio::getNext(){
	return next;
}
monomio* monomio::getPrev(){
	return prev;
}
void monomio::setNext(monomio *m){
	this->next = m;
}
void monomio::setPrev(monomio *m){
	this->prev = m;
}


void monomio::str2monomio(string s){
	bool bexp=false, bcoef=true, bdecimal=false;
	double coef=0, decimal=0, dec=10;
	int exp=0, exps=1;
	int signal = 1;
	if(s.at(0)=='x'){
		coef=1;
	}
	else if(s.at(0)=='-'){
		signal = -1;
		if(s.length()>1&&s.at(1)=='x'){
			coef = 1;
		}
	}
	for(int i=0; i<s.length(); i++){
		if(s.at(i)=='x'){
			bcoef=false;
		}
		if(s.at(i)=='^'&&!bcoef){
			bexp=true;
		}
		if(s.at(i)=='-'&&bexp){
			exps=-1;
		}
		if(bcoef){
			if(isNum(s.at(i))){
				if(!bdecimal){
					coef = coef*10+(double)(s.at(i)-48);
				}
				else{
					decimal = decimal + ((double)(s.at(i)-48))/dec;
					dec*=10;
				}
			}
			else if(s.at(i)=='.'&&!bdecimal){
				bdecimal=true;
			}
		}
		if(bexp){
			if(isNum(s.at(i))){
				exp += exp*10+(double)(s.at(i)-48);
			}
		}
	}
	if(!bcoef&&!bexp){
		exp = 1;
	}
	coef+=decimal;
	exp*=exps;
	this->coef = signal*coef;
	this->exp = exp;
}

int monomio::isNum(char s){
	return (s>=48&&s<=57);
}

int monomio::isValid(char s){
	if(s>=48&&s<=57||s==94||s==120){
		return 1;
	}
	return 0;
}

monomio* monomio::newMonomio(string s){
	monomio* n = new monomio(1, 1);
	return n;
}

void monomio::derivar(){
	coef = coef*exp;
	exp--;
}

string monomio::monomio2str(){
	stringstream buffer;
	if(!(coef==1&&exp!=0)){
		buffer << coef;
	}
	/*if(coef==-1&&exp!=0){
		buffer << "-";
	}*/
	if(coef!=0&&exp!=0){
		buffer << "x";
		if(exp!=1){
			buffer << "^" << exp;
		}
	}
	return buffer.str();
}

polinomio::polinomio(){
	this->cabeca = NULL;
	this->cauda = NULL;
	s="0";
}

polinomio::polinomio(double coef, int exp){
	this->cabeca = new monomio(coef, exp);
	this->cauda = this->cabeca;
	s=cabeca->monomio2str();
}

polinomio::polinomio(string s){
	this->cabeca = NULL;
	this->cauda = NULL;
	this->insert_string(s);
}

int polinomio::empty(){
	return cabeca==NULL;
}

void polinomio::insert(double coef, int exp){
	monomio *p, *ant;
	p = this->cabeca;
	ant = NULL;
	while(p!=NULL&&p->getExp()<exp){
		ant = p;
		p = p->getNext();
	}
	if(p!=NULL&&p->getExp()==exp){
		double oldCoef = p->getCoef();
		p->setCoef(oldCoef + coef);
	}
	else{
		monomio *novo = new monomio(coef, exp);
		if(ant==NULL){
			if(empty()){
				cabeca = novo;
				cauda = novo;
			}
			else{
				cabeca->setPrev(novo);
				novo->setNext(cabeca);
				cabeca = novo;
			}
		}
		else{
			if(p==NULL){
				cauda->setNext(novo);
				novo->setPrev(cauda);
				cauda = novo;
			}
			else{
				novo->setPrev(ant);
				novo->setNext(ant->getNext());
				ant->getNext()->setPrev(novo);
				ant->setNext(novo);
			}
		}	
	}
}

void polinomio::insert(string s){
	monomio *p, *ant;
	p = this->cabeca;
	ant = NULL;
	monomio aux(s);
	int exp = aux.getExp();
	double coef = aux.getCoef();
	while(p!=NULL&&p->getExp()<exp){
		ant = p;
		p = p->getNext();
	}
	if(p!=NULL&&p->getExp()==exp){
		double oldCoef = p->getCoef();
		p->setCoef(oldCoef + coef);
	}
	else{
		monomio *novo = new monomio(coef, exp);
		if(ant==NULL){
			if(empty()){
				cabeca = novo;
				cauda = novo;
			}
			else{
				cabeca->setPrev(novo);
				novo->setNext(cabeca);
				cabeca = novo;
			}
		}
		else{
			if(p==NULL){
				cauda->setNext(novo);
				novo->setPrev(cauda);
				cauda = novo;
			}
			else{
				novo->setPrev(ant);
				novo->setNext(ant->getNext());
				ant->getNext()->setPrev(novo);
				ant->setNext(novo);
			}
		}	
	}
}

void polinomio::percorrer(){
	cout << "\nDo princípio ao fim:\n";
	for(monomio *p=cabeca; p!=NULL; p=p->getNext()){
		cout << "Coeficiente: " << p->getCoef() << " Expoente: " << p->getExp() << endl;
	}
}

void polinomio::percorrer_inv(){
	cout << "\nDo fim ao princípio:\n";
	for(monomio *p=cauda; p!=NULL; p=p->getPrev()){
		cout << "Coeficiente: " << p->getCoef() << " Expoente: " << p->getExp() << endl;
	}
}

void polinomio::removerCoef(double coef){
	monomio *p, *ant;
	p = this->cabeca;
	ant = NULL;
	while(p!=NULL&&p->getCoef()!=coef){
		ant = p;
		p = p->getNext();
	}
	if(p!=NULL){
		if(ant==NULL){
			cabeca = cabeca->getNext();
			cabeca->setPrev(NULL);
		}
		else{
			ant->setNext(p->getNext());
			if(p->getNext()==NULL){
				cauda = p->getPrev();
			}
			else{
				p->getNext()->setPrev(ant);
			}
		}
		delete p;
	}
}

void polinomio::add_polinomio(const polinomio& p){
	for(monomio *ptr = p.cabeca; ptr!=NULL; ptr=ptr->getNext()){
		int exp = ptr->getExp();
		double coef = ptr->getCoef();
		this->insert(coef, exp);
	}
	this->removerCoef(0);
	this->atualizar_string();
}

polinomio polinomio::operator+(const polinomio& p){
	polinomio poli;
	for(monomio *ptr = p.cabeca; ptr!=NULL; ptr=ptr->getNext()){
		int exp = ptr->getExp();
		double coef = ptr->getCoef();
		poli.insert(coef, exp);
	}
	for(monomio *ptr = this->getCabeca(); ptr!=NULL; ptr=ptr->getNext()){
		int exp = ptr->getExp();
		double coef = ptr->getCoef();
		poli.insert(coef, exp);
	}
	poli.removerCoef(0);
	poli.atualizar_string();
	return poli;
}

polinomio polinomio::operator-(const polinomio& p){
	polinomio poli;
	for(monomio *ptr = p.cabeca; ptr!=NULL; ptr=ptr->getNext()){
		int exp = ptr->getExp();
		double coef = -ptr->getCoef();
		poli.insert(coef, exp);
	}
	for(monomio *ptr = this->getCabeca(); ptr!=NULL; ptr=ptr->getNext()){
		int exp = ptr->getExp();
		double coef = ptr->getCoef();
		poli.insert(coef, exp);
	}
	poli.removerCoef(0);
	poli.atualizar_string();
	return poli;
}

polinomio polinomio::operator*(const polinomio& p){
	polinomio poli;
	for(monomio *ptr = this->getCabeca(); ptr!=NULL; ptr = ptr->getNext()){
		for(monomio *ptr2 = p.cabeca; ptr2!=NULL; ptr2 = ptr2->getNext()){
			double coef = ptr->getCoef()*ptr2->getCoef();
			int exp = ptr->getExp()+ptr2->getExp();
			poli.insert(coef, exp);
		}
	}
	poli.removerCoef(0);
	poli.atualizar_string();
	return poli;
}

/*void polinomio::operator=(const polinomio &p){
	zerar();
	for(monomio *ptr2 = p.cabeca; ptr2!=NULL; ptr2 = ptr2->getNext()){
		double coef = ptr2->getCoef();
		int exp = ptr2->getExp();
		this->insert(coef, exp);
	}
	atualizar_string();
}*/

void polinomio::self_mult(string s){
	polinomio p(s);
	polinomio poli;
	for(monomio *ptr = this->getCabeca(); ptr!=NULL; ptr = ptr->getNext()){
		for(monomio *ptr2 = p.getCabeca(); ptr2!=NULL; ptr2 = ptr2->getNext()){
			double coef = ptr->getCoef()*ptr2->getCoef();
			int exp = ptr->getExp()+ptr2->getExp();
			poli.insert(coef, exp);
		}
	}
	poli.removerCoef(0);
	this->zerar();
	this->setPtrs(poli.getCabeca(), poli.getCauda());
	this->atualizar_string();
	poli.setPtrs(NULL, NULL);
}


void polinomio::setPtrs(monomio *cabeca, monomio *cauda){
	this->cabeca = cabeca;
	this->cauda = cauda;
}
void polinomio::atualizar_string(){
	monomio *ptr;
	ptr = cauda;
	if(empty()){
		s = "0";
	}
	else{
		s = "";
	}
	while(ptr!=NULL){
		s+=ptr->monomio2str();
		ptr = ptr->getPrev();
		if(ptr!=NULL&&ptr->getCoef()>0){
			s+='+';
		}
	}
}

double polinomio::p(double x){
	double value = 0;
	for(monomio *ptr = this->cabeca; ptr!=NULL; ptr=ptr->getNext()){
		value += ptr->getCoef()*pow(x, ptr->getExp());
	}
	return value;
}

void polinomio::escalar(double n){
	for(monomio *ptr = this->cabeca; ptr!=NULL; ptr=ptr->getNext()){
		double coef = ptr->getCoef();
		ptr->setCoef(n*coef);
	}
}

void polinomio::insert_string(string s){
	string aux = "";
	int begin=0;
	if(s.at(0)=='-'){
		aux="-";
		begin=1;
	}
	for(int i=0; i<s.length(); i++){
		if(s.at(i)=='+'){
			this->insert(aux);
			aux="";
		}
		else if(s.at(i)=='-'){
			if(i!=0&&s.at(i-1)=='^'){
				aux+=s.at(i);
				
			}
			else{
				this->insert(aux);
				aux = "-";
			}
		}
		else{
			aux+=s.at(i);
		}
	}
	if(aux.length()>0){
		this->insert(aux);
	}
	this->removerCoef(0);
	this->atualizar_string();
}

void polinomio::derivar(){
	for(monomio *ptr = this->cabeca; ptr!=NULL; ptr=ptr->getNext()){
		ptr->derivar();
	}
	this->removerCoef(0);
	this->atualizar_string();
}

void polinomio::zerar(){
	monomio *ptr, *ant;
	ptr = cabeca;
	ant = NULL;
	while(ptr!=NULL){
		ant = ptr;
		ptr = ptr->getNext();
		delete ant;
	}
	setPtrs(NULL, NULL);
	atualizar_string();
}

monomio* polinomio::getCabeca(){
	return this->cabeca;
}

monomio* polinomio::getCauda(){
	return this->cauda;
}

interpoLagrange::interpoLagrange(int grau){
	this->grau = grau;
	this->pts = new double[grau+1];
	this->y = new double[grau+1];
	this->L = new polinomio*[grau+1];
	for(int i=0; i<=grau; i++){
		L[i] = new polinomio("1");
	}
}
interpoLagrange::~interpoLagrange(){
	delete pts;
	delete y;
	for(int i=0; i<=grau; i++){
		L[i]->zerar();
		delete L[i];
	}
	delete L;
}
void interpoLagrange::inserir_pts(int p, double pto, double valor){
	if(p>=0&&p<=grau){
		this->pts[p] = pto;
		this->y[p] = valor;
	}
}

string interpoLagrange::binomio(double num){
	stringstream buffer;
	buffer << "x"; 
	if(num != 0){
		if(num>0){
			buffer << "+" << num;
		}
		else{
			buffer << num;
		}
	}
	return buffer.str();
}

void interpoLagrange::insertLs(){
	double divisor;
	for(int k=0; k<=grau; k++){
		divisor = 1;
		for(int j=0; j<=grau; j++){
			if(k==j){
				continue;
			}
			L[k]->self_mult(binomio(-pts[j]));
			divisor *= 1/(pts[k]-pts[j]);
		}
		L[k]->escalar(divisor);
		L[k]->escalar(y[k]);
	}
}

polinomio* interpoLagrange::interpolar(){
	polinomio *p = new polinomio;
	insertLs();
	for(int i=0; i<=grau; i++){
		p->add_polinomio(*L[i]);
	}
	return p;
}
