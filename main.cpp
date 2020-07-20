#include <iostream>
#include <cstdio>
#include "propagacao.h"

using namespace std;

int main(){
	FILE *result = fopen("resultados_finais.txt", "w");
	propagacao p;
	p.run_cgrasp(result, 1, 1000);
	aco a;
	a.run_aco(result, 1, 1000);
	system("mkdir resultados\nmv resultados_finais.txt resultados\nmv *#*.txt resultados\n");
	return 0;
}

/*int main(){
    polinomio *p;
    interpoLagrange l(2);
    l.inserir_pts(0, -1, 4);
    l.inserir_pts(1, 0, 1);
    l.inserir_pts(2, 2, -1);
    p = l.interpolar();
    polinomio p1("2x");
    polinomio p2("x^2-2");
    polinomio p3 = p1+p2;
    cout << p3 << endl;
    p3.derivar();
    cout << p3 << endl;
    p3.derivar();
    cout << p3 << endl;
    cout << p->p(-1) << endl;
    cout << p->p(0) << endl;
    cout << p->p(2) << endl;
    p->zerar();
    p1.zerar();
    p2.zerar();
    p3.zerar();
    delete p;
    return 0;
}*/