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
#define CPSTRSZ 16 // Chrom-pos size stringified size
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

typedef enum /* the old t_t typing of general words ... could be deleted */
{
    STRG, /* unknown type, so default to string */
    NUMS, /* NUMberS: but does not include currency. Date, time, float or int, i.e. just a sequence of numbers with maybe some special symbils.*/
    PNI, /* pos or beg int */
    STCP, /* string with closing punctuation attached.. a comma, or a full stop, semicolon, !? etc */
    SCST, /* string with starting capital */
    SCCP, /* string with starting capital AND closing punctuation */
    ALLC /* string with all caps */
} t_t;

typedef struct  /* optstruct, a struct for the options */
{
    int cflag; /* convert, print out output as ped and map, aftetr having resolved resolution events */
    char *iname; /* input file name which should be a tped file of course */
    char *fname; /* optional input tfam file for sanity checking of the tped/tfam pairing */
} optstruct;

typedef struct /* word type */
{
    char *w;
    t_t t; /* number or not */
    unsigned lp1; /* length of word plus one (i.e. includes null char */
} w_c;

typedef struct /* aw_c: array of words container */
{
    w_c **aw;
    unsigned ab;
    unsigned al;
    char *nn; /* this is the name name ... actually the CromPos string */
    int gd /* whether it's a chrompos or nay */, gdn /* the dupgroup */;
} aw_c;

typedef struct /* aaw_c: array of array of words container */
{
    size_t numl; /* number of lines, i.e. rows */
    aw_c **aaw; /* an array of pointers to aw_c */
} aaw_c;

/* checking each character can be comptiue-intensive, so I've offloaded off to MACROS */

/* Macro fo GET Leading Char TYPE */
/* so, this refers to the first character: "+-.0123456789" only these are allowed. These is how we infer
 * a quantity of some sort ... except for currency */
#define GETLCTYPE(c, typ); \
            if(((c) == 0x2B) | ((c) == 0x2D) | ((c) == 0x2E) | (((c) >= 0x30) && ((c) <= 0x39))) { \
                if( ((c) == 0x2B) | ((c) == 0x2D) | (((c) >= 0x30) && ((c) <= 0x39))) \
                    typ = PNI; \
                else \
                    typ = NUMS; \
            } else if(((c) >= 0x41) && ((c) <= 0x5A)) \
                    typ = SCST;

/* Macro for InWord MODify TYPE */
#define IWMODTYPEIF(c, typ); \
            if( ((typ) == NUMS) & (((c) != 0x2E) & (((c) < 0x30) || ((c) > 0x39)))) \
                typ = STRG; \
            else if( ((typ) == PNI) & (c == 0x2E)) \
                typ = NUMS; \
            else if( ((typ) == PNI) & ((c < 0x30) || (c > 0x39)) ) \
                typ = STRG;

/* Macro for SETting CLosing Punctuation TYPE, based on oldc (oc) not c-var */
/* 21=! 29=) 2C=, 2E=. 3B=; 3F=? 5D=] 7D=}*/
#define SETCPTYPE(oc, typ); \
            if( ((oc)==0x21)|((oc)==0x29)|((oc)==0x2C)|((oc)==0x2E)|((oc)==0x3B)|((oc)==0x3F)|((oc)==0x5D)|((oc)==0x7D) ) { \
                if((typ) == STRG) \
                    typ = STCP; \
                else if((typ) == SCST) \
                    typ = SCCP; \
            }
