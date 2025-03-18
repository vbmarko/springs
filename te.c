#include <stdio.h>
#include <stdlib.h>
int main(){

	if(0.1){
		printf("0.1 is true\n");
	}
	char *a = "0.0";
	double val = atof(a);
	if(0.0 == 0){
		printf("%g ==  0\n",atof(a));
	}
	return 0;

}
