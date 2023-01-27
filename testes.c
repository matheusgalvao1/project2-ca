#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    int *ponteiro, **ponteiroDoPonteiro, valor;
    valor = 50;
    ponteiro = &valor;
    *ponteiro = 60;
    ponteiroDoPonteiro = &ponteiro;
    **ponteiroDoPonteiro = 70;
    printf("\nValor de ponteiroDoPonteiro: %p\nEndereço de memoria de ponteiro: %p\nValor de ponteiro: %p\nEndereco de memoria de valor: %p\nValor: %d\n", ponteiroDoPonteiro, &ponteiro, ponteiro, &valor, valor);
    return 0;
}