/* read in a csv of floats, allowing for skipp first few rows or columns (i.e. heaards, rownames etc */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "fcsvrd.h"

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
            printf("s1+%u ", aawc->aaw[i]->aw[j]->lp1);
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
}

double **givemat(aaw_c *aawc, int skiprows, int skipcols) /* print line and word details, but not the words themselves */
{
    char *zerostr=malloc(64*sizeof(char));
    memset(zerostr, '0', 63);
    zerostr[63]='\0';
    zerostr[0]='4';
    zerostr[4]='.';
    zerostr[7]='2';
    printf("zerostr=%s\n", zerostr);
    printf("%6.12f\n", strtod(zerostr, NULL));
    int i, j, k;
    //allocate memory
    int numcols1st=aawc->aaw[skiprows]->al-skipcols; // we take this column to represent the size of all others.
    double **mat=malloc((aawc->numl-skiprows)*sizeof(double*));
    for(i=0;i<aawc->numl-skiprows;++i)
        mat[i]=calloc(numcols1st, sizeof(double));
    //assign values:
    for(i=0;i<aawc->numl-skiprows;++i) {
        for(j=0;j<numcols1st;++j) {
            memset(zerostr, '0', 63);
            for(k=0;k<aawc->aaw[i+skiprows]->aw[j+skipcols]->lp1-1;k++)
                zerostr[k]=aawc->aaw[i+skiprows]->aw[j+skipcols]->w[k];
            // mat[i][j]=strtod(aawc->aaw[i+skiprows]->aw[j+skipcols]->w, NULL);
            mat[i][j]=strtod(zerostr, NULL);
        }
    }
    free(zerostr);
    return(mat);
}

void prtmat(double **mat, int matrows, int matcols) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<matrows;++i)
        for(j=0;j<matcols;++j)
            printf((j!=matcols-1)?"%4.4f ":"%4.4f\n", mat[i][j]);
}

void prtmat0(double **mat, int firstrows, int firstcols) /* print line and word details, but not the words themselves */
{
    int i, j;
    printf("Matrix derived from csv: first %i rows and first %i columns:\n", firstrows, firstcols); 
    for(i=0;i<firstrows;++i)
        for(j=0;j<firstcols;++j)
            printf((j!=firstcols-1)?"%4.6f ":"%4.6f\n", mat[i][j]);
}

void prtmat00(double **mat, int firstrows, int firstcols, int maxprec) /* print line and word details, but not the words themselves */
{
    int i, j;
    printf("max precision: matrix derived from csv: first %i rows and first %i columns:\n", firstrows, firstcols); 
    for(i=0;i<firstrows;++i)
        for(j=0;j<firstcols;++j)
            printf((j!=firstcols-1)?"%6.*f ":"%6.*f\n", maxprec, mat[i][j]);
}

void freemat(double **mat, int matrows)
{
    int i;
    for(i=0;i<matrows;++i)
        free(mat[i]);
    free(mat);
}

void prtaawcplain20(aaw_c *aawc, int skiprows, int skipcols, int firstrows, int firstcols) /* print line and word details, but not the words themselves */
{
    int i, j;
    printf("First %i rows and first %i columns, with nrow=%i and ncol=%i skipped.\n", firstrows, firstcols, skiprows, skipcols); 
    for(i=skiprows;i<firstrows;++i) {
        for(j=skipcols;j<firstcols;++j) {
            printf((j!=firstcols-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
}

void prtaawcplain200(aaw_c *aawc, int skiprows, int skipcols, int firstrows, int firstcols) /* print line and word details, but not the words themselves */
{
    int i, j;
    int lword;
    int prec; //precision
    char *pt;
    printf("Precisions of first %i rows and first %i columns, with nrow=%i and ncol=%i skipped.\n", firstrows, firstcols, skiprows, skipcols); 
    for(i=skiprows;i<firstrows;++i) {
        for(j=skipcols;j<firstcols;++j) {
            lword=aawc->aaw[i]->aw[j]->lp1-2;
            pt=strchr(aawc->aaw[i]->aw[j]->w, '.');
            if(pt!=NULL)
                prec=(int)((aawc->aaw[i]->aw[j]->w+lword) - pt);
            else
                prec=0;
            printf((j!=firstcols-1)?"%s(prec=%i)  ":"%s(prec=%i)\n", aawc->aaw[i]->aw[j]->w, prec);
        }
    }
}

int maxprec(aaw_c *aawc, int skiprows, int skipcols) // workout max precision.
{
    int i, j;
    int lword;
    int prec; //precision
    int maxprec=0;
    char *pt;
    for(i=skiprows;i<aawc->numl;++i) {
        for(j=skipcols;j<aawc->aaw[i]->al;++j) {
            lword=aawc->aaw[i]->aw[j]->lp1-2;
            pt=strchr(aawc->aaw[i]->aw[j]->w, '.');
            if(pt!=NULL) {
                prec=(int)((aawc->aaw[i]->aw[j]->w+lword) - pt);
                if(prec>maxprec)
                    maxprec=prec;
            }
        }
    }
    return(maxprec);
}

void prtaawcplain0(aaw_c *aawc, int firstrows, int firstcols) /* print only first few rows and columns. */
{
    int i, j;
    printf("First %i rows and first %i columns:\n", firstrows, firstcols); 
    for(i=0;i<firstrows;++i) {
        for(j=0;j<firstcols;++j) {
            printf((j!=firstcols-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
    printf("Full csv file specs: nrow= %zu, ncol1st=%i\n", aawc->numl, aawc->aaw[0]->al);
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
            /* if word is a candidate for a NUM or PNI (i.e. via its first character), make sure it continues to obey rules: a MACRO */
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
    if(argc!=4) {
        printf("Error. Pls supply args 1) csvfile name 2) numrows (headers) to skip 3) numcols (rowname columns) to skip.\n");
        exit(EXIT_FAILURE);
    }
    int skiprows=atoi(argv[2]);
    int skipcols=atoi(argv[3]);

    aaw_c *aawc=processincsv(argv[1]);

    prtaawcplain0(aawc, 9, 9);
    printf("\n"); 
    prtaawcplain20(aawc, skiprows, skipcols, 9, 9);
    printf("\n"); 
    // prtaawcplain200(aawc, skiprows, skipcols, 12, 12);
    // printf("\n"); 
    int mxprec = maxprec(aawc, skiprows, skipcols);
    printf("Max precision is %i\n", mxprec);

    // can work out the following right now ... you for ease of viewing.
    int matrows=aawc->numl-skiprows;
    int matcols=aawc->aaw[skiprows]->al-skipcols; //yes the skiprows'th row is used to calculate ncols: all rows must have this numcols.

    // matrows and matcols not passed up up recalculated (perhaps dodgy) in func.
    double **mat0=givemat(aawc, skiprows, skipcols);
    prtmat00(mat0, 8, 8, mxprec);
    freemat(mat0, matrows);

    free_aawc(&aawc);

    return 0;
}
