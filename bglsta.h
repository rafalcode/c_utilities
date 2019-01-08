/* it's all too easy to start over-complicating this: for example quotations. Here you woul dneed to check the last 2 characters of everyword, not just the last one, i.e "stop!", that adds new layers. */
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#ifdef DBG
#define CBUF 2
#define WABUF 1
#define LBUF 2
#else
#define CBUF 12
#define WABUF 20
#define LBUF 32
#endif
#define NUMGTS 20
typedef unsigned char boole;
typedef enum /* gt_t, genotype type */
{
    AA,
    CC,
    GG,
    TT,
    AC, /* anything >= 0x100 is therefore hetzyg */
    AG,
    AT,
    CA,
    CG,
    CT,
    GA,
    GC,
    GT,
    TA,
    TC,
    TG,
    D1, /* only one of the alleles is a I or D ( policy will be to treat as uncalled) */
    D2, /* both alleles is a I or D ( policy will be to treat as uncalled) */
    Z1, /* one SNP uncalled or unrecognised */
    ZZ /* both uncalled, unrecognised, the obligatory catch-all so to speak */
} gt_t;
/* GT Name arrays based on above, "gtna" sets 00 as the unknown 17th GT, while gtna2 sets it to NN */
char gtna0[NUMGTS][3]={"AA", "CC", "GG", "TT", "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG", "D1", "D2", "Z1", "ZZ"};
char gtna2[NUMGTS][3]={"AA", "CC", "GG", "TT", "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG", "NN", "NN", "NN", "NN"};
char gtna[NUMGTS][3]={"AA", "CC", "GG", "TT", "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG", "00", "00", "00", "00"};

typedef struct  /* optstruct, a struct for the options */
{
    int cflag; /* convert, print out output as ped and map, aftetr having resolved resolution events */
    char *iname; /* input file name which should be a tped file of course */
    char *fname; /* optional input tfam file for sanity checking of the tped/tfam pairing */
    char *mname; /* the marker file */
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
    int gd /* whether it's a chrompos or nay */, gdn /* the dupgroup */;
} aw_c;

typedef struct /* aaw_c: array of array of words container */
{
    size_t numl; /* number of lines, i.e. rows */
    aw_c **aaw; /* an array of pointers to aw_c */
} aaw_c;
struct strchainodm /* for the marker file */
{
    aw_c *aw;
    struct strchainodm *n;
    int idx; // the index corresponding to this marker element: it's the sort of thing youex
    boole vis; // whether it's been visited or not ... a double check, if you will.
};
typedef struct strchainodm snodm;
