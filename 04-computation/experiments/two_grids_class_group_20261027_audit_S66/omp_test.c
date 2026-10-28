#include <omp.h>
#include <stdio.h>
int main(void) {
    int n = 0;
#pragma omp parallel reduction(+:n)
    n += 1;
    printf("threads %d\n", n);
    return 0;
}
