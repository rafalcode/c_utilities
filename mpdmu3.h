/* it's all too easy to start over-complicating this: for example quotations. Here you woul dneed to check the last 2 characters of everyword, not just the last one, i.e "stop!", that adds new layers. */
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h> // for the opts handling and abort()
#include<string.h>

#ifdef DBG
#define CBUF 4
#define WABUF 2
#define LBUF 4
#else
#define CBUF 12
#define WABUF 20
#define LBUF 32
#endif

#ifdef DBG
#define GBUF 4
#define WBUF 4
#else
#define GBUF 64
#define WBUF 64
#endif
#define MNCOLS 4 // mandatory number of columns

#define CPSTRSZ 16
#define NDEV 7 // number of duplicate resolution event types
#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
		memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
    }

typedef unsigned char boole;
typedef boole bool3; // not necessary, but makes t celar there's 3 values 0, 1, 2 for notseen, seen once, seen >1
typedef struct  /* optstruct, a struct for the options */
{
    int pflag; // output as ped and map, aftetr having resolved resolution events
    int tflag; // if you want the tech rep list printed out
    int eflag; // if you want to see the techrep resolution categories
    int fflag; // summary of Tech Rep resolutions
    int cflag; // TRDupsets with CN_PN strngs for comfort
    char *nf; // a file with SNP names
} optstruct;

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
    ZZ /* unrecognised, the obligatory catch-all so to speak */
} gt_t;

char gtna[17][3]={"AA", "CC", "GG", "TT", "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG", "00"};

typedef enum /* dv_t, dupevent type */
{
    NNNNNN, /* no event, initial value */
    OGTA00, /* 1 GT for all dups, and it's a 00 */
    OGTN00, /* TR are all same GT (not 00) */
    MGTNWI, /* many GTs not 00 but no winners */
    MGTNWZ, /* many GTs but no winners though there was a 00 among them */
    VGTB1W, /* valid GT competition but 1 winner */
    ZGTB1W /* ZZ GT incompetition but 1 winner */
} dv_t;

char devna[NDEV][7]={"NNNNNN", "OGTA00", "OGTN00", "MGTNWI", "MGTNWZ", "VGTB1W", "ZGTB1W"};

typedef struct /* i2g_t2 */
{
    unsigned **i, sz, bf;
    gt_t **gt; // dynamic .. per sample
    int **dpcou; // a count of the dups for this loc ... necessarily starts as 2.
    bool3 **mif; // after filtering: already mentioned in filter?
} i2g_t2; /* map indices to with gts */

typedef struct /* dgia_t */
{
    unsigned **is /* indices */, bf, sz;
    unsigned **js /* second-level indices ... refers to the global index */;
    gt_t gt; // category label index
} dgia_t; /* dupe index array */

typedef struct /* adgia_t */
{
    dgia_t **dg;
    unsigned bf, sz;
} adgia_t; /* dupe index array */

typedef struct /* dia_t */
{
    unsigned **is /* indices */, bf, sz;
    int lidx; // category label index
    gt_t wgt; // the the winning GT, only after processing
    unsigned wi; // the winning index
    char posstr[17]; /* position string, the thing we're hashing on */
} dia_t; /* dupe index array */

typedef struct /* adia_t */
{
    dia_t **d;
    unsigned bf, sz;
} adia_t; /* dupe index array */

typedef struct /* mp_t, map type, one line in the map file */
{
	char *n;
	char *nn; /* numbered name CXX_XXXXX, etc. will alway sbe 16 chars in length */
	size_t nsz; /* size of the name r ID field */
    char cnu[3]; // first column the chromosome number. 1-32 and Y
    // float cmo; // the centimorgans, it truly should be float, but opt for int
    int cmo; // the centimorgans, it truly should be float, but opt for int
	long pos; /* just the one number */
    bool3 gd;
    int gdn;
    boole ngd; // new versionsof the above for the hash on SNP names
    int ngdn;
    boole keep; // new versionsof the above for the hash on SNP names
} mp_t; /* map type */

typedef struct /* wff_t, word from file type */
{
	char *w;
	size_t wsz; /* size of the name r ID field */
    int gdn; /* will be zero not a duplicate */
    boole nomatch;
    boole hitmdup; // hit master dup
    boole hitddup; // hit discarded dup
    boole ddupwom; // a discarded whose master dup was not called.
} wff_t;

typedef struct /* wseq_t */
{
    size_t *wln;
    size_t wsbuf;
    size_t quan;
    size_t lbuf; /* a buffer for the number of lines */
    size_t numl; /* number of lines, i.e. rows */
    size_t *wpla; /* words per line array: the number of words on each line */
} wseq_t;

struct strchainodw
{
    wff_t *wff;
    struct strchainodw *n;
    int idx; // the index corresponding to this mp element: it's the sort of thing youex
};
typedef struct strchainodw snodw;


struct strchainode
{
    mp_t *mp;
    struct strchainode *n;
    int idx; // the index corresponding to this mp element: it's the sort of thing youex
};

typedef struct strchainode snod;


typedef enum /* t_t categorising words */
{
    STRG, /* unknown type, so default to string */
    NUMS, /* NUMberS: but does not include currency. Date, time, float or int, i.e. just a sequence of numbers with maybe some special symbils.*/
    PNI, /* pos or beg int */
    STCP, /* string with closing punctuation attached.. a comma, or a full stop, semicolon, !? etc */
    SCST, /* string with starting capital */
    SCCP, /* string with starting capital AND closing punctuation */
    ALLC /* string with all caps */
} t_t;

typedef struct /* w_c word type container */
{
    char *w;
    unsigned wbf;
    t_t t; /* number or not */
    unsigned lp1; /* length of word plus one (i.e. includes null char */
} w_c;

typedef struct /* aw_c: array of words container */
{
    w_c **aw;
    unsigned ab;
    unsigned al;
    char *nn;
    boole gd;
    int gdn;
    boole ngd; // new versionsof the above for the hash on SNP names
    int ngdn;
} aw_c;

typedef struct /* aaw_c: array of array of words container */
{
    size_t numl; /* number of lines, i.e. rows */
    aw_c **aaw; /* an array of pointers to aw_c */
} aaw_c;

struct strchainodm /* m for map */
{
    aw_c *mpaw;
    struct strchainodm *n;
    int idx; // the index corresponding to this mp element: it's the sort of thing youex
};
typedef struct strchainodm snodm;

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
