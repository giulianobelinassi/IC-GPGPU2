#!/bin/bash
# Este arquivo especifica os testes automatizados.

# Especifica os caracteres de cores.
GREEN='\033[1;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color


function print_ok {
	echo -e "[${GREEN}OK${NC}]\n"
}

function print_failed {
	echo -e "[${RED}FALHOU${NC}]\n"
}

#recompila tudo com as macros que habilitam testes que rodam internamente no programa.
make clean
make EXTRA=-DTEST_GHMATECD_CUDA

#remove os arquivos antigos.
rm -f SSOLO240E_-5+5.DAT SSOLO240D_-5+5.DAT SSOLO2160E_-5+5.DAT SSOLO2160D_-5+5.DAT

# Teste 1: O arquivo SSOLO240E_-5+5.DAT, que é a solução encontrada
# por este programa, deve ser igual ao SSOLO240E_-5+5.SOL, que é o
# arquivo calculado pelo programa compilado no compilador original
# (Compaq Visual Fortran 2003, no Windows) convertido para o formato
# Unix (sem o carriage return).

echo -e "Dado a execução do programa 'main_small'"

./main_small > /dev/null

echo -e "O arquivo gerado SSOLO240E_-5+5.DAT é igual ao SSOLO240E_-5+5.SOL ?"
if cmp -s "SSOLO240E_-5+5.DAT" "SSOLO240E_-5+5.SOL"; then
	print_ok
else 
	print_failed
fi

echo -e "O arquivo gerado SSOLO240D_-5+5.DAT é compatível em precião ao SSOLO240D_-5+5.SOL ?"
./dtester "SSOLO240D_-5+5.DAT" "SSOLO240D_-5+5.SOL"
if [ $? -eq 0 ]; then
	print_ok
else 
	print_failed
fi

echo -e "Dado a execução do programa 'main'"

./main > /dev/null

echo -e "O arquivo gerado SSOLO2160E_-5+5.DAT é igual ao SSOLO2160E_-5+5.SOL ?"
if cmp -s "SSOLO2160E_-5+5.DAT" "SSOLO2160E_-5+5.SOL"; then
	print_ok
else 
	print_failed
fi

echo -e "O arquivo gerado SSOLO2160D_-5+5.DAT é igual ao SSOLO2160D_-5+5.SOL ?"
if cmp -s "SSOLO2160D_-5+5.DAT" "SSOLO2160D_-5+5.SOL"; then
	print_ok
else 
	print_failed
fi

echo -e "Este é o fim dos testes.\n\n"
