/* vttgo3: taken from vtggo2
 * accepts two sub files (say Russian and then English)
 * only works with two autosubs.
 * Manual subs are totally different, do not try with them
 * BEWARE spaces after last word in line mplayer/mpv very fussy about these.*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "vttgo2.h"

w_c *crea_wc(unsigned initsz)
{
    w_c *wc=malloc(sizeof(w_c));
    wc->lp1=initsz;
    wc->t=STRG;
    wc->w=malloc(wc->lp1*sizeof(char));
    return wc;
}

void reall_wc(w_c **wc, unsigned *cbuf)
{
    w_c *twc=*wc;
    unsigned tcbuf=*cbuf;
    tcbuf += CBUF;
    twc->lp1=tcbuf;
    twc->w=realloc(twc->w, tcbuf*sizeof(char));
    *wc=twc; /* realloc can often change the ptr */
    *cbuf=tcbuf;
    return;
}

void norm_wc(w_c **wc)
{
    w_c *twc=*wc;
    twc->w=realloc(twc->w, twc->lp1*sizeof(char));
    *wc=twc; /* realloc can often change the ptr */
    return;
}

void free_wc(w_c **wc)
{
    w_c *twc=*wc;
    free(twc->w);
    free(twc);
    return;
}

aw_c *crea_awc(unsigned initsz)
{
    int i;
    aw_c *awc=malloc(sizeof(aw_c));
    awc->ab=initsz;
    awc->al=awc->ab;
    awc->aw=malloc(awc->ab*sizeof(w_c*));
    for(i=0;i<awc->ab;++i) 
        awc->aw[i]=crea_wc(CBUF);
    return awc;
}

void reall_awc(aw_c **awc, unsigned buf)
{
    int i;
    aw_c *tawc=*awc;
    tawc->ab += buf;
    tawc->al=tawc->ab;
    tawc->aw=realloc(tawc->aw, tawc->ab*sizeof(aw_c*));
    for(i=tawc->ab-buf;i<tawc->ab;++i)
        tawc->aw[i]=crea_wc(CBUF);
    *awc=tawc;
    return;
}

void norm_awc(aw_c **awc)
{
    int i;
    aw_c *tawc=*awc;
    /* free the individual w_c's */
    for(i=tawc->al;i<tawc->ab;++i) 
        free_wc(tawc->aw+i);
    /* now release the pointers to those freed w_c's */
    tawc->aw=realloc(tawc->aw, tawc->al*sizeof(aw_c*));
    *awc=tawc;
    return;
}

void free_awc(aw_c **awc)
{
    int i;
    aw_c *tawc=*awc;
    for(i=0;i<tawc->al;++i) 
        free_wc(tawc->aw+i);
    free(tawc->aw); /* unbelieveable: I left this out, couldn't find where I leaking the memory! */
    free(tawc);
    return;
}

aaw_c *crea_aawc(unsigned initsz)
{
    int i;
    unsigned lbuf=initsz;
    aaw_c *aawc=malloc(sizeof(aaw_c));
    aawc->numl=0;
    aawc->aaw=malloc(lbuf*sizeof(aw_c*));
    for(i=0;i<initsz;++i) 
        aawc->aaw[i]=crea_awc(WABUF);
    return aawc;
}

void free_aawc(aaw_c **aw)
{
    int i;
    aaw_c *taw=*aw;
    for(i=0;i<taw->numl;++i) /* tried to release 1 more, no go */
        free_awc(taw->aaw+i);
    free(taw->aaw);
    free(taw);
}

void prtaawcdbg(aaw_c *aawc)
{
    int i, j, k;
    for(i=0;i<aawc->numl;++i) {
        printf("l.%u(%u): ", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf("w_%u: ", j); 
            if(aawc->aaw[i]->aw[j]->t == NINT) {
                printf("NINT! "); 
                continue;
            } else if(aawc->aaw[i]->aw[j]->t == NFLT) {
                printf("NFLT! "); 
                continue;
            } else if(aawc->aaw[i]->aw[j]->t == NTIM) {
                printf("NTIM! "); 
                continue;
            } else if(aawc->aaw[i]->aw[j]->t == LTW) {
                printf("LTW! "); 
                continue;
            } else if(aawc->aaw[i]->aw[j]->t == STRG) {
                printf("STRG! "); 
                continue;
            }
            for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                putchar(aawc->aaw[i]->aw[j]->w[k]);
            printf("/%u ", aawc->aaw[i]->aw[j]->lp1-1); 
        }
        printf("\n"); 
    }
}

void prtspec(aaw_c *aawc, aaw_c *aawc2)
{
    int i, j, k;
    int ii;
    for(i=0;i<aawc->numl;++i) {
        if((i<8) & (i!=5)) {
            for(j=0;j<aawc->aaw[i]->al;++j) {
                for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                    putchar(aawc->aaw[i]->aw[j]->w[k]);
                if(j!=aawc->aaw[i]->al-1)
                    putchar(' ');
            }
            printf("\n"); 
        } else {
            ii=(i-12)%8;
            if(ii==0) {
                // printf("[i==%i]: ", i);
                for(j=0;j<aawc->aaw[i]->al;++j) {
                    for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                        putchar(aawc->aaw[i]->aw[j]->w[k]);
                    if(j!=aawc->aaw[i]->al-1)
                        putchar(' ');
                }
                printf("\n"); 
            } else if(ii==5) {
                for(j=0;j<aawc->aaw[i]->al;++j) {
                    for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                        putchar(aawc->aaw[i]->aw[j]->w[k]);
                    if(j!=aawc->aaw[i]->al-1)
                        putchar(' ');
                }
                printf("\n"); 
                for(j=0;j<aawc2->aaw[i]->al;++j) {
                    for(k=0;k<aawc2->aaw[i]->aw[j]->lp1-1; k++)
                        putchar(aawc2->aaw[i]->aw[j]->w[k]);
                    if(j!=aawc2->aaw[i]->al-1)
                        putchar(' ');
                }
                printf("\n"); 
            } else if(ii==6)
                printf("\n"); 
        }
    }
}

void prtspec2(aaw_c *aawc, aaw_c *aawc2)
{
    int i, j, k;
    int ii;
    int leastlines=(aawc->numl>aawc2->numl)?aawc2->numl:aawc->numl;
    for(i=0;i<leastlines;++i) {
        if((i<8) & (i!=5)) {
            for(j=0;j<aawc->aaw[i]->al;++j) {
                for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                    putchar(aawc->aaw[i]->aw[j]->w[k]);
                if(j!=aawc->aaw[i]->al-1)
                    putchar(' ');
            }
            printf("\n"); 
            if(i==6) {
                for(j=0;j<aawc2->aaw[i]->al;++j) {
                    for(k=0;k<aawc2->aaw[i]->aw[j]->lp1-1; k++)
                        putchar(aawc2->aaw[i]->aw[j]->w[k]);
                    if(j!=aawc2->aaw[i]->al-1)
                        putchar(' ');
                }
                printf("\n"); 
            }
        } else {
            ii=(i-12)%8;
            if(ii==0) {
                // printf("[i==%i]: ", i);
                for(j=0;j<aawc->aaw[i]->al;++j) {
                    for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                        putchar(aawc->aaw[i]->aw[j]->w[k]);
                    if(j!=aawc->aaw[i]->al-1)
                        putchar(' ');
                }
                printf("\n"); 
            } else if(ii==5) {
                for(j=0;j<aawc->aaw[i]->al;++j) {
                    for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                        putchar(aawc->aaw[i]->aw[j]->w[k]);
                    if(j!=aawc->aaw[i]->al-1)
                        putchar(' ');
                }
                printf("\n"); 
                for(j=0;j<aawc2->aaw[i]->al;++j) {
                    for(k=0;k<aawc2->aaw[i]->aw[j]->lp1-1; k++)
                        putchar(aawc2->aaw[i]->aw[j]->w[k]);
                    if(j!=aawc2->aaw[i]->al-1)
                        putchar(' ');
                }
                printf("\n"); 
            } else if(ii==6)
                printf("\n"); 
        }
    }
}

aaw_c *processinpf(char *fname)
{
    /* declarations */
    FILE *fp=fopen(fname,"r");
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */;
    int c, oldc='\0';
    boole inword=0;
    unsigned lbuf=LBUF /* buffer for number of lines */, cbuf=CBUF /* char buffer for size of w_c's: reused for every word */;
    aaw_c *aawc=crea_aawc(lbuf); /* array of words per line */

    while( (c=fgetc(fp)) != EOF) {
        if( (c== '\n') | (c == ' ') | (c == '\t') ) {
            if( inword==1) { /* cue word-ending procedure */
                aawc->aaw[aawc->numl]->aw[couw]->w[couc++]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=couc;
                IWMODTYPEIF(oldc, aawc->aaw[aawc->numl]->aw[couw]->t);
                norm_wc(aawc->aaw[aawc->numl]->aw+couw);
                couw++; /* verified: this has to be here */
            }
            if(c=='\n') { /* cue line-ending procedure */
                if(aawc->numl ==lbuf-1) {
                    lbuf += LBUF;
                    aawc->aaw=realloc(aawc->aaw, lbuf*sizeof(aw_c*));
                    for(i=lbuf-LBUF; i<lbuf; ++i)
                        aawc->aaw[i]=crea_awc(WABUF);
                }
                aawc->aaw[aawc->numl]->al=couw;
                norm_awc(aawc->aaw+aawc->numl);
                aawc->numl++;
                couw=0;
            }
            inword=0;
        } else if(inword==0) { /* a normal character opens word */
            if(couw ==aawc->aaw[aawc->numl]->ab-1) /* new word opening */
                reall_awc(aawc->aaw+aawc->numl, WABUF);
            couc=0;
            cbuf=CBUF;
            aawc->aaw[aawc->numl]->aw[couw]->w[couc++]=c;
            GETLCTYPE(c, aawc->aaw[aawc->numl]->aw[couw]->t); /* MACRO: the firt character gives a clue */
            inword=1;
        } else if(inword) { /* simply store */
            if(couc == cbuf-1)
                reall_wc(aawc->aaw[aawc->numl]->aw+couw, &cbuf);
            aawc->aaw[aawc->numl]->aw[couw]->w[couc++]=c;
            IWMODTYPEIF(c, aawc->aaw[aawc->numl]->aw[couw]->t);
        }
        oldc=c;
    } /* end of big for statement */
    fclose(fp);

    /* normalization stage */
    for(i=aawc->numl; i<lbuf; ++i)
        free_awc(aawc->aaw+i);
    aawc->aaw=realloc(aawc->aaw, aawc->numl*sizeof(aw_c*));

    return aawc;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=3) {
        printf("Error. Pls supply 2 arguments: first and second subtitle files.\n");
        exit(EXIT_FAILURE);
    }

    aaw_c *aawc=processinpf(argv[1]);
    aaw_c *aawc2=processinpf(argv[2]);
    prtspec2(aawc, aawc2);

    free_aawc(&aawc);
    free_aawc(&aawc2);

    return 0;
}
