#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define CBUF 8
#define GBUF 8
#define WABUF 4
#define LBUF 8

typedef unsigned char boole;

typedef struct /* word type */
{
    char *w;
    unsigned b; /* buffer */
    unsigned lp1; /* length */
} w_t;

typedef struct /* wa_t: word array */
{
    w_t **wa;
    unsigned ab;
    unsigned al;
} wa_t;

typedef struct /* wseq_t */
{
    size_t numl; /* number of lines, i.e. rows */
    wa_t **awat; /* an array of pointers to wa_t */
} wseq_t;

