#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"csvfocg.h"

int cmpocabyoc(const void *a, const void *b) /* compare aoca by occurence */
{
    oca_t *oca = (oca_t*)a; /* cast our void! */
    oca_t *ocb = (oca_t*)b; /* cast our void! */
    return oca->uo  - ocb->uo; /* integer comparison: returns positive if b > a and nagetive if a > b: i.e. highest values first */
}

int cmpocabyoasz(const void *a, const void *b) /* compare aoca by occurence */
{
    oca_t *oca = (oca_t*)a; /* cast our void! */
    oca_t *ocb = (oca_t*)b; /* cast our void! */
    return ocb->ocisz  - oca->ocisz; /* integer comparison: returns positive if b > a and nagetive if a > b: i.e. highest values first */
}

void prtaoca0(oca_t *aoca, int aocasz)
{
    int i, j;
    for(j=0; j<aocasz;++j) {
        printf("word %s (ocisz=%i): ", aoca[j].uo, aoca[j].ocisz);
        for(i=0;i<aoca[j].ocisz;++i) 
            printf("%u ", aoca[j].ocinda[i]);
        putchar('\n');
    }
}

void prtaoca(aaw_c *aawc, oca_t *aoca, int aocasz)
{
    int i, j, k, ii, iidx, n;
    int scmax; // score max
    int intfield=8;

    // header
    for(n=0;n<aawc->aaw[0]->al;++n)
        printf((n!=aawc->aaw[0]->al-1)?"%s,":"%s\n", aawc->aaw[0]->aw[n]->w);

    // this will aawc in order of highest aocasz, you can correct this in R.
    for(j=0; j<aocasz;++j) {
        if(aoca[j].ocisz>1) {
            scmax=0;
            for(i=0;i<aoca[j].ocisz;++i) {
                k = aoca[j].ocinda[i]; // convenience variable
                ii=atoi(aawc->aaw[k]->aw[intfield]->w);
                if(ii>scmax) {
                    scmax=ii;
                }
            }
            for(i=0;i<aoca[j].ocisz;++i) {
                k = aoca[j].ocinda[i]; // convenience variable
                ii=atoi(aawc->aaw[k]->aw[intfield]->w);
                if(ii==scmax) {
                    for(n=0;n<aawc->aaw[iidx]->al;++n)
                        printf((n!=aawc->aaw[k]->al-1)?"%s,":"%s\n", aawc->aaw[k]->aw[n]->w);
                }
            }
        } else {
            iidx = aoca[j].ocinda[0];
            for(n=0;n<aawc->aaw[iidx]->al;++n)
                printf((n!=aawc->aaw[iidx]->al-1)?"%s,":"%s\n", aawc->aaw[iidx]->aw[n]->w);
        }
    }
}

void prtaoca2(oca_t *aoca, int aocasz)
{
    int i, j;
    for(j=0; j<aocasz;++j) {
        printf("word %s of size %i has ocisz=%i\n ", aoca[j].uo, aoca[j].uosz, aoca[j].ocisz);
    }
}

oca_t *uniquevals(aaw_c *aawc, int *aocasz_, int col)
{
    unsigned char new;
    unsigned i, j;
    unsigned aocabuf=GBUF;
    char *colword=NULL;
    int aocasz=0;
    oca_t *aoca=malloc(aocabuf*sizeof(oca_t));
    for(i=0;i<aocabuf;++i) {
        aoca[i].uo=NULL;
        aoca[i].ocibf=GBUF;
        aoca[i].ocinda=malloc(aoca[i].ocibf*sizeof(int));
        aoca[i].ocisz=0;
    }

    for(i=0; i<aawc->numl;++i) {
        new=1;
        colword=aawc->aaw[i]->aw[col]->w;
        for(j=0; j<aocasz;++j) {
            if(!strcmp(aoca[j].uo, aawc->aaw[i]->aw[col]->w)) {
                CONDREALLOC(aoca[j].ocisz, aoca[j].ocibf, GBUF, aoca[j].ocinda, int);
#ifdef DBG
                printf("APPENDING aocaind %i (ou=%s) @sz %i @buf %i adding vind %i\n", j, aawc->aaw[i]->aw[col]->w, aoca[j].ocisz, aoca[j].ocibf, i);
#endif
                aoca[j].ocinda[aoca[j].ocisz] = i;
                aoca[j].ocisz++;
                new=0;
                break;
            }
        }
        if(new) {
            aocasz++;
            CONDREALLOC2(aocasz, aocabuf, GBUF, aoca);
#ifdef DBG
                printf("NEW aocaind %i (ou=%s) @sz %i @buf %i adding vind %i\n", j, colword, aoca[j].ocisz, aoca[j].ocibf, i);
#endif
            aoca[aocasz-1].uosz = aawc->aaw[i]->aw[col]->lp1;
            aoca[aocasz-1].uo = calloc(aawc->aaw[i]->aw[col]->lp1, sizeof(char));
            strcpy(aoca[aocasz-1].uo, colword);
            aoca[aocasz-1].ocisz++;
            aoca[aocasz-1].ocinda[aoca[j].ocisz-1] = i;
        }
    }

    // normalizing
    for(i=aocasz;i<aocabuf;++i) {
        free(aoca[i].ocinda);
    }
    aoca=realloc(aoca, aocasz*sizeof(oca_t));

    /* order */
    // qsort(aoca, aocasz, sizeof(oca_t), cmpocabyoc);
    qsort(aoca, aocasz, sizeof(oca_t), cmpocabyoasz);
#ifdef DBG
    printf("number of different values: %i\n", aocasz);
    for(j=0; j<aocasz;++j) {
        printf("%i (%u): ", aoca[j].uo, aoca[j].ocisz);
        for(i=0;i<aoca[j].ocisz;++i) 
            printf("%i ", aoca[j].ocinda[i]);
        printf("\n"); 
    }
#endif
    *aocasz_ = aocasz;
    return aoca;
}

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
            if(aawc->aaw[i]->aw[j]->t == NUMS) {
                printf("NUM! "); 
                continue;
            } else if(aawc->aaw[i]->aw[j]->t == PNI) {
                printf("PNI! "); 
                continue;
            } else if(aawc->aaw[i]->aw[j]->t == STCP) {
                printf("STCP! "); 
                continue;
            }
            for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                putchar(aawc->aaw[i]->aw[j]->w[k]);
            printf("/%u ", aawc->aaw[i]->aw[j]->lp1-1); 
        }
        printf("\n"); 
    }
}

void prtaawcdata(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf("l%ut", aawc->aaw[i]->aw[j]->lp1-1);
            switch(aawc->aaw[i]->aw[j]->t) {
                case NUMS: printf("N "); break;
                case PNI: printf("I "); break;
                case STRG: printf("S "); break;
                case STCP: printf("C "); break; /* closing punctuation */
                case SCST: printf("Z "); break; /* starting capital */
                case SCCP: printf("Y "); break; /* starting capital and closing punctuation */
                case ALLC: printf("A "); break; /* horrid! all capitals */
            }
        }
    }
    printf("\n"); 
	printf("L is a line, l is length of word, S is normal string, C closing punct, Z, starting cap, Y Starting cap and closing punct.\n"); 
}

void prtaawcplain(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf("l%u/", aawc->aaw[i]->aw[j]->lp1-1);
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
}

void prtaawcp5(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=3;i<aawc->numl;++i) {
        j=(i-3)%12;
        printf((j!=11)?"%s ":"%s\n", aawc->aaw[i]->aw[4]->w);
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
                SETCPTYPE(oldc, aawc->aaw[aawc->numl]->aw[couw]->t);
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
            /* if word is a candidate for a NUM or PNI (i.e. via its first character), make sure it continues to obey rules: a MACRO */
            IWMODTYPEIF(c, aawc->aaw[aawc->numl]->aw[couw]->t);
        }
        oldc=c;
    } /* end of big for statement */
    fclose(fp);

    /* normalization stage */
    for(i=aawc->numl; i<lbuf; ++i) {
        free_awc(aawc->aaw+i);
    }
    aawc->aaw=realloc(aawc->aaw, aawc->numl*sizeof(aw_c*));

    return aawc;
}

aaw_c *processincsv(char *fname)
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
        if( (c== '\n') | (c == ',') ) {
            if( inword==1) { /* cue word-ending procedure */
                aawc->aaw[aawc->numl]->aw[couw]->w[couc++]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=couc;
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
            inword=1;
        } else if(inword) { /* simply store */
            if(couc == cbuf-1)
                reall_wc(aawc->aaw[aawc->numl]->aw+couw, &cbuf);
            aawc->aaw[aawc->numl]->aw[couw]->w[couc++]=c;
        }
        oldc=c;
    } /* end of big for statement */
    fclose(fp);

    /* normalization stage */
    for(i=aawc->numl; i<lbuf; ++i) {
        free_awc(aawc->aaw+i);
    }
    aawc->aaw=realloc(aawc->aaw, aawc->numl*sizeof(aw_c*));

    return aawc;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=2) {
        printf("Error. Pls supply argument (name of text file).\n");
        exit(EXIT_FAILURE);
    }
#ifdef DBG2
    printf("typeszs: aaw_c: %zu aw_c: %zu w_c: %zu\n", sizeof(aaw_c), sizeof(aw_c), sizeof(w_c));
#endif

    aaw_c *aawc=processincsv(argv[1]);
#ifdef DBG
    prtaawcdbg(aawc);
// #else
    // prtaawcp5(aawc); // printout original text as well as you can.
    prtaawcplain(aawc); // printout original text as well as you can.
#endif
    // printf("Numlines: %zu\n", aawc->numl); 
    int aocasz;
    int focusfield=0; // which field to focus on
    oca_t *aoca = uniquevals(aawc, &aocasz, focusfield);
    prtaoca(aawc, aoca, aocasz);

    for(int i=0;i<aocasz;++i) {
        free(aoca[i].uo);
        free(aoca[i].ocinda);
    }
    free(aoca);
    free_aawc(&aawc);

    return 0;
}
