#define _GNU_SOURCE
/*#include <stdio.h>*/
#include <sys/resource.h>

/*Subrotina que deve ser chamada pelo programa em fortran para pedir*/
/*um tamanho não limitado da pilha, pois no programa desenvolvido pelo*/
/*Ronaldo, não há alocação no Heap*/
void request_unlimited_stack_(int* ret)
{
	struct rlimit rl;
	int result;

	result = getrlimit(RLIMIT_STACK, &rl);
	
	if (result != 0)
	{
		/*Request denied by the OS*/
		*ret = 1;
		return;
	}
	rl.rlim_cur = RLIM_INFINITY;
	result = setrlimit(RLIMIT_STACK, &rl);
	if (result != 0)
	{
		*ret = 1;
		/*Request denied by the OS*/
		return;
	}
	*ret = 0;
}
/*
void allocate_junk_in_stack(void)
{

	int A[20*1024*1024];
	A[20*1024*1024-1] = 1;
}

int main(void)
{
	int ret;
	request_unlimited_stack_(&ret);
	if (ret != 0)
		puts("Request Failed");
	
	allocate_junk_in_stack();	
	return 0;
}
*/
