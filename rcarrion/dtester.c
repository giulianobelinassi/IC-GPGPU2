#include <stdio.h>
#include <string.h>
#include <complex.h>

#define BUF_SIZE 4096


static char buffer[BUF_SIZE];


void get_dimension_data(int* n,
					   int* nbe,
					   int* ni,
					   FILE* file
					   )
{
	#define MALHA_STR    " N\332MERO DE ELEMENTOS DA MALHA="
	#define CONTORNO_STR " N\332MERO DE ELEMENTOS DE CONTORNO="
	#define INTERNO_STR  " N\332MERO DE PONTOS INTERNOS="

	rewind(file);

	while (fgets(buffer, BUF_SIZE, file) != NULL)
	{
		if (strstr(buffer, MALHA_STR))
		{
			sscanf(buffer, MALHA_STR "%5d\n", n);
		}
		else if (strstr(buffer, CONTORNO_STR))
		{
			sscanf(buffer, CONTORNO_STR "%5d", nbe);
		}
		else if (strstr(buffer, INTERNO_STR))
		{
			sscanf(buffer, INTERNO_STR "%5d", ni);
		}

	}
	rewind(file);
}

void get_nos_contorno(int nbe, float complex nos_contorno[nbe][3], FILE* file)
{
	float reals[3][2];
	int index, i, j;

	rewind(file);

	#define NOS_CONTORNO_STR " N\323S DO CONTORNO"

	while (fgets(buffer, BUF_SIZE, file) != NULL)
	{
		if (strstr(buffer, NOS_CONTORNO_STR))
		{
			break;
		}
	}

	/*Descarta as tres próximas linhas*/
	for (i = 0; i < 4; ++i)
	{
		fgets(buffer, BUF_SIZE, file);
	}

	/*Lê de fato os resultados.*/
	for (i = 0; i < nbe; ++i)
	{
		#define FORMAT_67 "%4d  %f,%f  %f,%f  %f,%f"
		sscanf(buffer, FORMAT_67, &index, &reals[0][0], &reals[0][1], &reals[1][0], &reals[1][1], &reals[2][0], &reals[2][1]);
		
		for (j = 0; j < 3; ++j)
		{
			nos_contorno[index][j] = reals[j][0] + I*reals[j][1];
		}
	}
}


int main(int argc, char* argv[])
{
	FILE* file;
	int n, nbe, ni;


	file = fopen(argv[1], "r");
	if (!file)	
	{
		fputs("Error: Invalid filepath\n", stderr);
		return 1;
	}
	get_dimension_data(&n, &nbe, &ni, file);

	float complex nos_contorno[nbe][3];
	get_nos_contorno(nbe, nos_contorno, file);


	return 0;
}
