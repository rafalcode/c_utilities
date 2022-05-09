/* it's all too easy to start over-complicating this: for example quotations. Here you woul dneed to check the last 2 characters of everyword, not just the last one, i.e "stop!", that adds new layers. */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h> // for the opts handling and abort()

#ifdef DBG
#define CBUF 2
#define WABUF 1
#define LBUF 2
#else
#define CBUF 12
#define WABUF 20
#define LBUF 32
#endif

typedef unsigned char boole;
typedef unsigned char trival;
typedef struct  /* optstruct, a struct for the options */
{
    int dflag; // show duplicate names in one group per line
} optstruct;

typedef struct /* word type */
{
    char *w;
    unsigned lp1; /* length of word plus one (i.e. includes null char */
} w_c;

typedef struct /* aw_c: array of words container */
{
    w_c **aw;
    unsigned ab;
    unsigned al;
    // only one of the following will be used each time, I just didn't wnat to reset.
    boole ma1; // a match with the first name file?
    boole ma2; // a match with the second name file?
    trival dup;
    unsigned dgrp;
} aw_c;

typedef struct /* aaw_c: array of array of words container */
{
    long long numl; /* number of lines, i.e. rows */
    aw_c **aaw; /* an array of pointers to aw_c */
} aaw_c;

struct strchainodm /* m for map */
{
    aw_c *aw;
    struct strchainodm *n;
    int idx; // the index corresponding to this mp element: it's the sort of thing youex
};
typedef struct strchainodm snodm;

