/* experiments in reverting back to previous order */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define NR 4
#define NC 12
#define N NC*NR

int main(int argc, char *argv[])
{
    int i, j;
   /* declarations */
   int *a=malloc(N*sizeof(int));
   int *b=malloc(N*sizeof(int));
   int *c=malloc(N*sizeof(int));
   int *d=malloc(N*sizeof(int));
   for(i=0;i<N;++i) 
       a[i]=i;

   for(i=0;i<NR;++i) 
       for(j=0;j<NC;++j) 
           b[NR*j+i] = a[NC*i+j];

   for(i=0;i<NR;++i) {
       for(j=0;j<NC;++j) {
           c[NC*i+j] = b[NR*j+i];
           d[NC*i+j] = NR*j+i;
       }
   }

   for(i=0;i<N;++i) 
       printf("%i ", b[i]); 
   printf("\n"); 
   for(i=0;i<N;++i) 
       printf("%i ", i+b[b[i]]); 
   printf("\n"); 
   for(i=0;i<N;++i) 
       printf("%i ", d[i]); 
   printf("\n"); 
   for(i=0;i<N;++i) 
       printf("%i ", d[i]-b[b[i]]); 
   printf("\n"); 



   // so this finally is how it is done
   for(i=0;i<N;++i) 
       b[i] =i;
   for(i=0;i<N;++i) 
       printf("%i ", b[i]); 
   printf("\n"); 

    free(a);
    free(b);
    free(c);
    free(d);
   return 0;
}
