/* based on genread.h */
#include<stdio.h>
#include<stdlib.h>
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

typedef unsigned char boole;
typedef enum
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
    ZZ /* unrecognised */
} t_t;

typedef struct /* sample genotype array type */
{
    char **iid; /* an array of individ IDs */
    unsigned nsamps, gasz; /* number of samples, uniform length of genotype array */
    t_t** ga; /* an array length nsamps, of genotype arrays each of length gasz */
    boole tseles; /* tab sep allele-pairs, like from Neogen */
} sampga_t; /* sample genotype array */

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
} aw_c;

typedef struct /* aaw_c: array of array of words container */
{
    size_t numl; /* number of lines, i.e. rows */
    aw_c **aaw; /* an array of pointers to aw_c */
} aaw_c;
