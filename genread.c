/* modification of matread but operating on words instead of floats */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "genread.h"

w_t *creawt(unsigned initsz)
{
    w_t *wt=malloc(sizeof(w_t));
    wt->b=initsz;
    wt->lp1=0;
    wt->w=malloc(wt->b*sizeof(char));
    return wt;
}

void reallwt(w_t **wt, unsigned buf)
{
    w_t *twt=*wt;
    twt->b += buf;
    twt->w=realloc(twt->w, twt->b*sizeof(char));
    *wt=twt; /* realloc can often change the ptr */
    return;
}

void normwt(w_t **wt)
{
    w_t *twt=*wt;
    twt->w=realloc(twt->w, twt->lp1*sizeof(char));
    *wt=twt; /* realloc can often change the ptr */
    return;
}

void freewt(w_t **wt)
{
    w_t *twt=*wt;
    free(twt->w);
    free(twt);
    return;
}

wa_t *creawat(unsigned initsz)
{
    int i;
    wa_t *wat=malloc(sizeof(wa_t));
    wat->ab=initsz;
    wat->al=0;
    wat->wa=malloc(wat->ab*sizeof(w_t*));
    for(i=0;i<wat->ab;++i) 
        wat->wa[i]=creawt(CBUF);
    return wat;
}

void reallwat(wa_t **wat, unsigned buf)
{
    int i;
    wa_t *twat=*wat;
    twat->ab += buf;
    twat->wa=realloc(twat->wa, twat->ab*sizeof(wa_t*));
    for(i=twat->ab-buf;i<twat->ab;++i)
        twat->wa[i]=creawt(CBUF);
    *wat=twat;
    return;
}

void normwat(wa_t **wat)
{
    int i;
    wa_t *twat=*wat;
    /* free the individual w_t's */
    for(i=twat->al;i<twat->ab;++i) 
        freewt(twat->wa+i);
    /* now release the pointers to those freed w_t's */
    twat->wa=realloc(twat->wa, twat->al*sizeof(wa_t*));
    *wat=twat;
    return;
}

void freewat(wa_t **wat)
{
    int i;
    wa_t *twat=*wat;
    for(i=0;i<twat->al;++i) 
        freewt(twat->wa+i);
    free(twat);
    return;
}

wseq_t *create_wseq_t(unsigned initsz)
{
    int i;
    unsigned lbuf=initsz;
    wseq_t *awpl=malloc(sizeof(wseq_t));
    awpl->numl=0;
    awpl->awat=malloc(lbuf*sizeof(wa_t*));
    for(i=0;i<initsz;++i) 
        awpl->awat[i]=creawat(WABUF);
    return awpl;
}

void free_wseq(wseq_t **wa)
{
    int i;
    wseq_t *twa=*wa;
    for(i=0;i<twa->numl;++i) 
        freewat(twa->awat+i);
    free(twa->awat);
    free(twa);
}

wseq_t *processinpf(char *fname)
{
    /* declarations */
    FILE *fp=fopen(fname,"r");
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */;
    int c;
    boole inword=0;
    unsigned lbuf=LBUF;
    wseq_t *awpl=create_wseq_t(lbuf); /* array of words per line */

    while( (c=fgetc(fp)) != EOF) {
        if( (c== '\n') | (c == ' ') | (c == '\t') ) {
            if( inword==1) { /* cue word -edning procedure */
                awpl->awat[awpl->numl]->wa[couw]->w[couc++]='\0';
                awpl->awat[awpl->numl]->wa[couw]->lp1=couc;
                normwt(awpl->awat[awpl->numl]->wa+couw);
                couw++;
            }
            if(c=='\n') { /* cue line-ending procedure */
                if(awpl->numl ==lbuf-1) {
                    lbuf += LBUF;
                    awpl->awat=realloc(awpl->awat, lbuf*sizeof(wa_t*));
                    for(i=lbuf-LBUF; i<lbuf; ++i)
                        awpl->awat[i]=creawat(WABUF);
                }
                awpl->awat[awpl->numl]->al=couw;
                normwat(awpl->awat+awpl->numl);
                awpl->numl++;
                couw=0;
            }
            inword=0;
        } else if(inword==0) { /* a normal character opens word */
            if(couw ==awpl->awat[awpl->numl]->ab-1) /* new word opening */
                reallwat(awpl->awat+awpl->numl, WABUF);
            couc=0;
#ifdef DBG
            printf("numl: %zu couw: %zu couc: %zu\n", awpl->numl, couw, couc); 
#endif
            awpl->awat[awpl->numl]->wa[couw]->w[couc++]=c;
            inword=1;
        } else if(inword) { /* simply store */
            if(couc == awpl->awat[awpl->numl]->wa[couw]->b-1)
                reallwt(awpl->awat[awpl->numl]->wa+couw, CBUF);
            awpl->awat[awpl->numl]->wa[couw]->w[couc++]=c;
        }
    } /* end of big for statement */
    fclose(fp);

    /* normalization stage */
    for(i=awpl->numl; i<lbuf; ++i)
        freewat(awpl->awat+i);
    awpl->awat=realloc(awpl->awat, awpl->numl*sizeof(wa_t*));

    return awpl;
}

void prtwseq(wseq_t *awpl)
{
    int i, j, k;
    for(i=0;i<awpl->numl;++i) {
        printf("l.%u(%u): ", i, awpl->awat[i]->al); 
        for(j=0;j<awpl->awat[i]->al;++j) {
            printf("w_%u: ", j); 
            for(k=0;k<awpl->awat[i]->wa[j]->lp1-1; k++)
                putchar(awpl->awat[i]->wa[j]->w[k]);
            printf("/%u ", awpl->awat[i]->wa[j]->lp1-1); 
        }
        printf("\n"); 
    }
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=2) {
        printf("Error. Pls supply argument (name of text file).\n");
        exit(EXIT_FAILURE);
    }
#ifdef DBG
    printf("typeszs: wseq_t: %zu wa_t: %zu w_t: %zu\n", sizeof(wseq_t), sizeof(wa_t), sizeof(w_t));
#endif

    wseq_t *awpl=processinpf(argv[1]);
    prtwseq(awpl);
    printf("Numlines: %zu\n", awpl->numl); 

    free_wseq(&awpl);

    return 0;
}
