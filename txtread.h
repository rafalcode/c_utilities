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

typedef unsigned char boole;
typedef enum /* type of word might be, it is */
{
    STRG, /* unknown type, so default to string */
    NUMS, /* NUMberS: but does not include currency. Date, time, float or int, i.e. just a sequence of numbers with maybe some special symbils.*/
    PNI, /* pos or beg int */
    FLO, /* floating pt. */
    STPP, /* string with closing period symbol */
    STPA, /* string with closing colon: a pause */
    ENMR, /* enumerating ... digits ending with a period */
    STCP, /* string with closing punctuation attached.. a comma, or a full stop, semicolon, !? etc */
    SCST, /* string with starting capital */
    SCCP, /* string with starting capital AND closing punctuation */
    LCPA, /* leading char is a parenthesis */
    LCLT, /* leading char is a less than ... HT/XML prob */
    TNUM, /* time number: type numbers .... mainly a colon inside */
    ALLC /* string with all caps */
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

/* Macro fo GET Leading Char TYPE: 3C is <, 28 is ( */
#define GETLCTYPE(c, typ); \
            if(((c) == ':') | ((c) == '-') | ((c) == '.') | (((c) >= '0') && ((c) <= '9'))) { \
                if( ((c) == '-') | (((c) >= 0x30) && ((c) <= 0x39))) \
                    typ = PNI; \
                else \
                    typ = TNUM; \
            } else if(((c) >= 'A') && ((c) <= 'Z')) \
                    typ = SCST;

/* the 2nd version of this macro was the original, but it makes the a leading + or - quite important, despite the fact that these are quite
 * rare in text couched number text */
#define GETLCTYPE2(c, typ); \
            if(((c) == 0x2B) | ((c) == 0x2D) | ((c) == 0x2E) | (((c) >= 0x30) && ((c) <= 0x39))) { \
                if( ((c) == 0x2B) | ((c) == 0x2D) | (((c) >= 0x30) && ((c) <= 0x39))) \
                    typ = PNI; \
                else \
                    typ = NUMS; \
            } else if(((c) >= 0x41) && ((c) <= 0x5A)) \
                    typ = SCST;

/* Macro for InWord MODify TYPE */
/* 2C is , 3A is : 2E is . Times are to the minute, */
/* if PNI integer,demote to NUMS if : or . is seen. If nums, demote to strng if something else seen */
#define IWMODTYPEIF(c, typ); \
            if( ((typ) == PNI) & ((c == 0x2E))) \
                typ = FLO; \
            else if( ((typ) == PNI) & (c == 0x3A)) \
                typ = TNUM; \
            else if( ((typ) == NUMS) & (((c) != 0x2E) & ((c) != 0x3A) & (((c) < 0x30) || ((c) > 0x39)))) \
                typ = STRG; \
            else if( ((typ) == PNI) & ((c < 0x30) || (c > 0x39)) ) \
                typ = STRG;

/* Macro for SETting CLosing Punctuation TYPE, based on oldc (oc) not c-var */
/* 21=! 29=) 2C=, 2E=. 3A=: 3B=; 3E=> 3F=? 5D=] 7D=}*/
#define SETCPTYPE(oc, typ); \
            if( ((oc)==0x21)|((oc)==0x29)|((oc)==0x2C)|((oc)==0x2E)|((oc)==0x3A)|((oc)==0x3B)|((oc)==0x3F)|((oc)==0x5D)|((oc)==0x7D)|((oc)==0x3E) ) { \
                if( ((oc)==0x2E) & ((typ) == STRG) ) \
                    typ = STPP; \
                else if( ((oc)==0x2E) & ((typ) == NUMS) ) \
                    typ = ENMR; \
                else if( ((oc)==0x3A) & ((typ) == TNUM) ) \
                    typ = STPA; \
                else if((typ) == STRG) \
                    typ = STCP; \
                else if((typ) == SCST) \
                    typ = SCCP; \
            }

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-0)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
    }
