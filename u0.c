#include <stdio.h>
#include <wchar.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{

   // wchar_t *a="En la venta del Molinillo, que está puesta en los fines de los famosos";
   int i;
   int *a=NULL;
   char *s="En la venta del Molinillo, que está puesta en los fines de los famosos";
   size_t sl=strlen(s);
   printf("%zu\n", sl); 
   char *b=calloc(1+strlen(s), sizeof(char));
   for(i=0;i<strlen(s);++i) 
       b[i]=s[i];

   /* print stuff out */
   printf("%s\n", b);

   return 0;
}
