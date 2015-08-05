#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define CBUF 8
#define GBUF 8
#define WBUF 8
#define WABUF 8

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
    size_t *wln; /* Mne on this unfathomable: word lines? what those? wods per line perhaps? */
    size_t wsbuf;
    size_t quan;
    size_t lbuf; /* a buffer for the number of lines */
    size_t numl; /* number of lines, i.e. rows */
    size_t *wpla; /* words per line array: the number of words on each line */
    wa_t **awat; /* an array of pointers to wa_t */
} wseq_t;

