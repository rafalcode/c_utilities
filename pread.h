/* it's all too easy to start over-complicating this: for example quotations. Here you woul dneed to check the last 2 characters of everyword, not just the last one, i.e "stop!", that adds new layers. */
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

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-0)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
    }

typedef unsigned char boole;
typedef enum
{
    STRG, /* unknown type, so default to string */
    NUM, /*number integer*/
    STPU, /* string which ends with punctuation :,;- */
    STES, /* string with closing full stop, uestion mark or exclamation which Ends Sentence */
} t_t;

typedef struct /* word container */
{
    char *w;
    t_t t; /* number or not */
    unsigned lp1; /* length of word plus one (i.e. will include null char, for good and for bad) */
} w_c;

typedef struct /* aw_c: array of words containers: essentially, a line of words */
{
    w_c **aw; /* a pointer to an array of words. */
    unsigned ab; /* array buffer */
    unsigned al; /* array length ... number of words in the array */
    /* you can add characteristics here is you like */
    short stsps; /* number of starting spaces */
    short sttbs; /* number of starting tabs */
} aw_c;

typedef struct /* aaw_c: array of array of words container */
{
    size_t numl; /* number of lines, i.e. rows */
    aw_c **aaw; /* an array of pointers to aw_c */
    int *ppa; /* paragraph point array: will record lines which have double newlines */
    int ppb, ppsz; /* buffer and size for our ppa */
} aaw_c;

/* checking each character can be compute-intensive, so I've offloaded off to MACROS, see below */

/* Macro for GET Leading Char TYPE */
#define GETLCTYPE(c, typ); \
            if(((c) >= 0x30) && ((c) <= 0x39)) \
                typ = NUM; \
            else \
                typ = STRG;

/* Macro for InWord MODify TYPE */
#define IWMODTYPEIF(c, typ); \
            if( ((typ) == NUM) & (((c) != ':') & ((c) != ',') & (((c) < 0x30) || ((c) > 0x39)))) \
                typ = STRG;

/* Macro for SETting CLosing Punctuation TYPE, based on oldc (oc) not c-var */
#define SETCPTYPE(oc, typ); \
            if( ((oc)==':') |((oc)==',')|((oc)==';')|((oc)=='-')) { \
                if((typ) == STRG) \
                    typ = STPU; \
            } else if( ((oc)=='!') | ((oc)=='?') |((oc)=='.')) { \
                if((typ) == STRG) \
                    typ = STES;  \
            }
