#include <stdio.h>
#include <stdlib.h>
#include "pi.h"


int main(){
	mpf_t  var;

	Borwein(1000,1,1000000,var);

	gmp_printf ("%.*Ff", 10000, var);//https://gmplib.org/manual/Formatted-Output-Strings.html


	return 0;
}


