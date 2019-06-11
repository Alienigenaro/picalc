#include "pi.h"


//Calcula o enézimo termo da serie
void BorweinTermo(unsigned long int n, unsigned long int pr, mpf_t retval){
    //(4n!/(n!)^4)*(26390*n+1103)/(396^(4*n))
	//A*B/C, onde A = (4n!/(n!)^4); B = (26390*n+1103) e C =(396^(4*n))
	mpf_t A, B, C, auxf;
	mpz_t auxi;
	mpf_init2(retval,pr);
	mpf_init2(A,pr);
	mpf_init2(B,pr);
	mpf_init2(C,pr);
	mpf_init2(auxf,pr);
	mpz_init2(auxi,pr);
	//Calculando parcela A = ((4n)!/(n!)^4):
	mpz_fac_ui(auxi,n*4); 	//Calcula 4*n fatorial
	mpf_set_z(A, auxi);		//A = (4*n)!
	mpz_fac_ui(auxi,n);		//auxi = n!
	mpf_set_z(auxf, auxi);	//auxf = n!
	mpf_pow_ui(auxf,auxf,4);//auxf = n!^4
	mpf_div(A, A, auxf);	//A = ((4n)!/(n!)^4)

	//Calculando a parcela B = (26390*n+1103):
	mpf_add_ui(B,B,26390*n+1103); //B = 26390*n+1103
	//Calculando a parcela C =(396^(4*n):
	mpf_add_ui(C,C,396);		//C = 396
	mpf_pow_ui(C,C,4*n); 	//C =(396^(4*n)

	mpf_mul(retval, A,B);
	mpf_div(retval, retval, C);

	mpf_clear (A);
	mpf_clear (B);
	mpf_clear (C);
	mpf_clear (auxf);
	mpz_clear (auxi);

}

//Recebe o número de iterações, número máximo de threads, precisao
void Borwein(unsigned long int it,int th, unsigned long int pr,mpf_t  retval){//https://crypto.stanford.edu/pbc/notes/pi/ramanujan.html
    register unsigned long int i;
    mpf_t parcela;
    mpf_t aux;
    //Inidializa retval com 0
    mpf_init2(retval,pr);//https://gmplib.org/manual/Initializing-Floats.html#Initializing-Floats
    mpf_init2(parcela,pr);
    mpf_init2(aux,pr);

    //https://gmplib.org/manual/Float-Arithmetic.html#Float-Arithmetic
    for(i=0;i<it;i++){
		BorweinTermo(i,pr, parcela);		//Calcula a parcela
        mpf_add(retval,retval, parcela);	//add parcela a retval
    }
	mpf_sqrt_ui(aux,8);
	mpf_mul(retval, retval, aux);
	mpf_ui_div(retval, 9801,retval);
}






