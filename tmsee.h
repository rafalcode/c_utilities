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
typedef enum
{
    STRG, /* unknown type, so default to string */
    NUM, /*number integer*/
    TMNG, /* Timing this is clearly for srt files */
    STCP, /* string with closing punctuation attached.. a comma, or a full stop, semicolon, !? etc */
    SCST, /* string with starting capital */
    SCCP, /* string with starting capital AND closing punctuation */
} t_t;

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

/* checking each character can be compute-intensive, so I've offloaded off to MACROS */

/* Macro fo GET Leading Char TYPE */
/* so, this refers to the first character: "+-.0123456789" only these are allowed. These is how we infer
 * a quantity of some sort ... except for currency */
#define GETLCTYPE(c, typ); \
            if(((c) >= 0x30) && ((c) <= 0x39)) \
                typ = TMNG; \
            else if(((c) >= 0x41) && ((c) <= 0x5A)) \
                typ = SCST; \
            else \
                typ = STRG;

/* Macro for InWord MODify TYPE */
#define IWMODTYPEIF(c, typ); \
            if( ((typ) == NUM) & (((c) == ':') & ((c) == ',') & (((c) >= 0x30) || ((c) <= 0x39)))) \
                typ = TMNG; \
            else if( ((typ) == NUM) & (((c) != ':') & ((c) != ',') & (((c) < 0x30) || ((c) > 0x39)))) \
                typ = STRG;

/* Macro for SETting CLosing Punctuation TYPE, based on oldc (oc) not c-var */
/* 21=! 29=) 2C=, 2E=. 3B=; 3F=? 5D=] 7D=}*/
#define SETCPTYPE(oc, typ); \
            if( ((oc)=='!')|((oc)==')')|((oc)==',')|((oc)=='.')|((oc)==';')|((oc)=='?')|((oc)==']')|((oc)=='}') ) { \
                if((typ) == STRG) \
                    typ = STCP; \
            }
