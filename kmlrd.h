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
    ORDW, // ordinary word
    WHEW, // when word
    GXW // gx word
} t_t;

typedef struct /* hmst_t, hours mins secs thousandths time format */
{
    int h, m, s;
} hmst_t;

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
    boole seengxw; // the line where the first locpoint occurs (for kmls)
    boole seenwhew; // the line where the first when (timepoint) occurs (for kmls)
    size_t flp; // the line where the first locpoint occurs (for kmls)
    size_t ftp; // the line where the first timepoint occurs (a "when" ..for kmls)
    size_t numl; /* number of lines, i.e. rows */
    aw_c **aaw; /* an array of pointers to aw_c */
} aaw_c;

/* checking each character can be comptiue-intensive, so I've offloaded off to MACROS */

/* Macro fo GET Leading Char TYPE */
/* so, this refers to the first character: "+-.0123456789" only these are allowed. These is how we infer
 * a quantity of some sort ... except for currency */
#define GETLCTYPE(c, typ); \
            if((c) == 'g') \
                typ = GXW; \
            else if((c) == 'w') \
                typ = WHEW;

/* Macro for InWord MODify TYPE */
#define IWMODTYPEIF(c, typ); \
            if(((c) != 'x') & ((c) != ':') & ((c) != 'c') & ((c) != 'o') & ((c) != 'r') & ((c) != 'd') & ((typ) == GXW)) \
                typ = ORDW; \
            else if(((c) != 'h') & ((c) != 'e') & ((c) != 'n') & ((typ) == WHEW)) \
                typ = ORDW;

/* Macro for SETting CLosing Punctuation TYPE, based on oldc (oc) not c-var */
/* 21=! 29=) 2C=, 2E=. 3B=; 3F=? 5D=] 7D=}*/
#define SETCPTYPE(oc, typ); \
            if(((oc) != 'd') & ((typ) == GXW)) \
                typ = ORDW; \
            else if(((oc) != 'n') & ((typ) == WHEW)) \
                typ = ORDW;
