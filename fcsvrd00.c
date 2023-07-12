/* this should be quoite useful especially hwen examining matrices from R
 * which are output as binary floating representation (to avoid losses).
 *
 * it read sin an calculates romeans, rowvars, comenas, colvrs
 *
 * Also the sum .. which is the depth of the library
 */
/* read in a csv of floats, allowing for skipp first few rows or columns (i.e. heaards, rownames etc */
/* we shall lean on %a, the binary floating point representation for losslessness in handling floats. */

/* fcsvrd00.c actually a fcsvcmp better */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
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

double **givemat0(aaw_c *aawc, int skiprows, int skipcols) /* print line and word details, but not the words themselves */
{
    /* OK what you're trying to do here could be a little tricky ... you're deliving into floatin point representation
     * check out this stov link:
     * https://stackoverflow.com/questions/16839658/printf-width-specifier-to-maintain-precision-of-floating-point-value
     *
     * the roudning func ions in math.h they just round to integers
     *
     * there is also this, on the latest C99 standard fp hex representation
     * https://www.exploringbinary.com/hexadecimal-floating-point-constants/
     *
     */
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
            // for(k=0;k<aawc->aaw[i+skiprows]->aw[j+skipcols]->lp1-1;k++)
            for(k=0;k<8;k++)
                zerostr[k]=aawc->aaw[i+skiprows]->aw[j+skipcols]->w[k];
            // mat[i][j]=strtod(aawc->aaw[i+skiprows]->aw[j+skipcols]->w, NULL);
            mat[i][j]=strtod(zerostr, NULL);
        }
    }
    free(zerostr);
    return(mat);
}

double **givemat(aaw_c *aawc, int skiprows, int skipcols) /* print line and word details, but not the words themselves */
{
    /* OK what you're trying to do here could be a little tricky ... you're deliving into floatin point representation
     * check out this stov link:
     * https://stackoverflow.com/questions/16839658/printf-width-specifier-to-maintain-precision-of-floating-point-value
     *
     * the roudning func ions in math.h they just round to integers
     *
     * there is also this, on the latest C99 standard fp hex representation
     * https://www.exploringbinary.com/hexadecimal-floating-point-constants/
     *
     */
    int i, j;
    //allocate memory
    int numcols1st=aawc->aaw[skiprows]->al-skipcols; // we take this column to represent the size of all others.
    double **mat=malloc((aawc->numl-skiprows)*sizeof(double*));
    for(i=0;i<aawc->numl-skiprows;++i)
        mat[i]=calloc(numcols1st, sizeof(double));

    //assign values:
    for(i=0;i<aawc->numl-skiprows;++i) {
        for(j=0;j<numcols1st;++j) {
            mat[i][j] = strtod(aawc->aaw[i+skiprows]->aw[j+skipcols]->w, NULL);
        }
    }
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

void prtcmpmats(double **mat, double **mat2, int nr, int nc)
{
    int i, j;
    printf("Only printing out differences:\n");
    for(i=0;i<nr;++i)
        for(j=0;j<nc;++j)
            if(mat2[i][j] != mat[i][j]) {
                printf((j!=nc-1)?"diff r%i c%i %6.14f ":"diff r%i c%i %6.14f\n", i, j, mat2[i][j] - mat[i][j]);
                printf((j!=nc-1)?"diff r%i c%i %a ":"diff r%i c%i %a\n", i, j, mat2[i][j] - mat[i][j]);
            }
}

void prtcmpmats2(double **mat, double **mat2, int nr, int nc)
{
    /* print percentage differences */
    int i, j;
    printf("Only printing out differences in percentage of first value:\n");
    for(i=0;i<nr;++i)
        for(j=0;j<nc;++j)
            if(mat2[i][j] != mat[i][j]) {
                printf("diff r%i c%i %4.8f%%\n", i, j, 100*(mat2[i][j] - mat[i][j])/mat[i][j]);
            }
}

void prtcmpmats3(double **mat, double **mat2, int nr, int nc)
{
    /* preparing for many many difference */
    int i, j;
    unsigned quandiffs=0;
    double pdiff, allpdiffs=.0; // all percentage differences
    double minpdiff=1e64, maxpdiff=.0;
    for(i=0;i<nr;++i)
        for(j=0;j<nc;++j)
            if(mat2[i][j] != mat[i][j]) {
                quandiffs++;
                pdiff = 100*(mat2[i][j] - mat[i][j])/mat[i][j];
                if(pdiff<minpdiff)
                    minpdiff=pdiff;
                if(pdiff>maxpdiff)
                    maxpdiff=pdiff;
                allpdiffs += pdiff;
            }
    // printf("Quantity of differences = %u (%2.2f%%) / Minpct diff: %6.6f%%; Maxpct diff: %6.6f%%; Avgpct diff: %6.6f%%\n", quandiffs, 100.*quandiffs/(nc*nr), minpdiff, maxpdiff, allpdiffs/quandiffs);
    printf("Quantity of differences = %u (%2.2f%%) / Minpct diff: %e%%; Maxpct diff: %e%%; Avgpct diff: %e%%\n", quandiffs, 100.*quandiffs/(nc*nr), minpdiff, maxpdiff, allpdiffs/quandiffs);
    printf("Binfloat version - Quantity of differences = %u (%2.2f%%) / Minpct diff: %a%%; Maxpct diff: %a%%; Avgpct diff: %a%%\n", quandiffs, 100.*quandiffs/(nc*nr), minpdiff, maxpdiff, allpdiffs/quandiffs);
}

void prtmat000b(double **mat, int nr, int nc)
{
    int i, j;
    printf("As ibinary (lossless) floating point matrix:\n"); 
    for(i=0;i<nr;++i)
        for(j=0;j<nc;++j)
            printf((j!=nc-1)?"%a ":"%a\n", mat[i][j]);
}

void prtvec(double *vec, int n)
{
    int i;
    for(i=0;i<n;++i)
        printf((i!=n-1)?"%4.6f ":"%4.6f\n", vec[i]);
}

void sumcol(double **mat, double *csu, int nr, int nc) // do your sums on the matrix
{
    int i, j;
    for(i=0;i<nr;++i)
        for(j=0;j<nc;++j)
            csu[j] += mat[i][j];

    printf("Your library size (colSums) is as follows:\n"); 
    for(j=0;j<nc;++j)
        printf((j!=nc-1)?"%4.6f ":"%4.6f\n", csu[j]); 
}

void sumrow(double **mat, double *rsu, int nr, int nc) // do your sums on the matrix
{
    int i, j;
    for(i=0;i<nr;++i)
        for(j=0;j<nc;++j)
            rsu[i] += mat[i][j];

    printf("Your difference rowSums are as follows:\n"); 
    for(i=0;i<nr;++i)
        printf((i!=nr-1)?"%4.6f ":"%4.6f\n", rsu[i]); 
}

void meansmat(double **mat, double *rmeans, double *cmeans, int nr, int nc) // do your sums on the matrix
{
    int i, j;
    double *csu=calloc(nc, sizeof(double));
    double *rsu=calloc(nr, sizeof(double));
    for(i=0;i<nr;++i)
        for(j=0;j<nc;++j) {
            rsu[i] += mat[i][j];
            csu[j] += mat[i][j];
        }
    // get means from the sums
    for(i=0;i<nr;++i)
        rmeans[i] = rsu[i]/nc; // yep watch that nc
    for(j=0;j<nc;++j)
        cmeans[j] = csu[j]/nr;
    free(csu);
    free(rsu);
}

void varsmat(double **mat, double *rvar, double *cvar, double *rmeans, double *cmeans, int nr, int nc) // do your sums on the matrix
{
    int i, j;
    double *csu=calloc(nc, sizeof(double));
    double *rsu=calloc(nr, sizeof(double));
    for(i=0;i<nr;++i)
        for(j=0;j<nc;++j) {
            rsu[i] += pow(mat[i][j]-rmeans[i], 2);
            csu[j] += pow(mat[i][j]-cmeans[j], 2);
        }
    for(i=0;i<nr;++i)
        rvar[i] = rsu[i]/nc; // yep watch that nc
    for(j=0;j<nc;++j)
        cvar[j] = csu[j]/nr;
    free(csu);
    free(rsu);
}

void freemat(double **mat, int matrows)
{
    int i;
    for(i=0;i<matrows;++i)
        free(mat[i]);
    free(mat);
}

void prtaawcplain2(aaw_c *aawc, int skiprows, int skipcols)
{
    int i, j;
    printf("Print aawc with row=%i and ncol=%i skipped.\n", skiprows, skipcols); 
    for(i=skiprows;i<aawc->numl;++i) {
        for(j=skipcols;j<aawc->aaw[i]->al;++j) {
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
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

twoints_t *maxminprec(aaw_c *aawc, int skiprows, int skipcols) // workout max precision.
{
    // zero precision is not counted
    int i, j;
    int lword;
    int prec; //precision
    int maxprec=0;
    int minprec=100000;
    char *pt;
    for(i=skiprows;i<aawc->numl;++i) {
        for(j=skipcols;j<aawc->aaw[i]->al;++j) {
            lword=aawc->aaw[i]->aw[j]->lp1-2;
            pt=strchr(aawc->aaw[i]->aw[j]->w, '.');
            if(pt!=NULL) {
                prec=(int)((aawc->aaw[i]->aw[j]->w+lword) - pt);
                if(prec>maxprec)
                    maxprec=prec;
                if(prec<minprec)
                    minprec=prec;
            }
        }
    }
    twoints_t *mxmn=calloc(1, sizeof(twoints_t));
    mxmn->x=maxprec;
    mxmn->y=minprec;
    return(mxmn);
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
    if(argc!=5) {
        printf("fcsvrd00: reads into two binfloat matrices and compares\n");
        printf("Error. Pls supply args 1) csvfile name 2) second csv 3) numrows (headers) to skip 4) numcols (rowname columns) to skip.\n");
        exit(EXIT_FAILURE);
    }
    int skiprows=atoi(argv[3]);
    int skipcols=atoi(argv[4]);

    // first csv
    aaw_c *aawc=processincsv(argv[1]);

#ifdef DBG
    prtaawcplain(aawc);
    prtaawcplain2(aawc, skiprows, skipcols);
#endif

    // prtaawcplain200(aawc, skiprows, skipcols, 12, 12);
    // printf("\n"); 
    // int mxprec = maxprec(aawc, skiprows, skipcols);
    // printf("Max precision is %i\n", mxprec);
    // twoints_t *mxmn=maxminprec(aawc, skiprows, skipcols);
    // printf("precision max:%i min:%i\n", mxmn->x, mxmn->y); 


    // can work out the following right now ... you for ease of viewing.
    int matrows=aawc->numl-skiprows;
    int matcols=aawc->aaw[skiprows]->al-skipcols; //yes the skiprows'th row is used to calculate ncols: all rows must have this numcols.
    double *cmeans=calloc(matcols, sizeof(double));
    double *rmeans=calloc(matrows, sizeof(double));
    double *cvar=calloc(matcols, sizeof(double));
    double *rvar=calloc(matrows, sizeof(double));


    double **mat0=givemat(aawc, skiprows, skipcols);
    free_aawc(&aawc);

    // OK get second file
    aawc=processincsv(argv[2]);
    double **mat2=givemat(aawc, skiprows, skipcols);
    free_aawc(&aawc);
    double *cmean2=calloc(matcols, sizeof(double));
    double *rmean2=calloc(matrows, sizeof(double));
    double *cvar2=calloc(matcols, sizeof(double));
    double *rvar2=calloc(matrows, sizeof(double));

#ifdef DBG
    prtmat000(mat0, matrows, matcols);
    prtmat000b(mat0, matrows, matcols);
#endif

    double *csu=calloc(matcols, sizeof(double)); // this is for library size
    double *rsu=calloc(matrows, sizeof(double)); // this is for rowSums
    double *csu2=calloc(matcols, sizeof(double)); // this is for library size
    double *rsu2=calloc(matrows, sizeof(double)); // this is for rowSums

#ifdef DBG
    sumcol(mat0, csu, matrows, matcols);
    sumcol(mat2, csu2, matrows, matcols);
    sumrow(mat0, rsu, matrows, matcols);

    meansmat(mat0, rmeans, cmeans, matrows, matcols);
    varsmat(mat0, rvar, cvar, rmeans, cmeans, matrows, matcols);
    meansmat(mat2, rmean2, cmean2, matrows, matcols);
    varsmat(mat2, rvar2, cvar2, rmean2, cmean2, matrows, matcols);
#endif
#ifdef DBG
    printf("%i rmeans:\n", matrows); 
    prtvec(rmeans, matrows);
    printf("%i rvar:\n", matrows); 
    prtvec(rvar, matrows);
    printf("Col-wise tallies:\n"); 
    printf("%i cmeans:\n", matcols); 
    prtvec(cmeans, matcols);
    printf("%i cmean2:\n", matcols); 
    prtvec(cmean2, matcols);
    printf("%i cvar:\n", matcols); 
    prtvec(cvar, matcols);
    printf("%i cvar2:\n", matcols); 
    prtvec(cvar2, matcols);

    prtcmpmats(mat0, mat2, matrows, matcols);
    printf("\n"); 
#endif
    prtcmpmats3(mat0, mat2, matrows, matcols); 

    // matrows and matcols not passed up up recalculated (perhaps dodgy) in func.
    // prtmat00(mat0, 2, 2, mxprec);

    // freeing up section:
    freemat(mat0, matrows);
    freemat(mat2, matrows);
    free(csu);
    free(rsu);
    free(rmeans);
    free(cmeans);
    free(rvar);
    free(cvar);
    free(csu2);
    free(rsu2);
    free(rmean2);
    free(cmean2);
    free(rvar2);
    free(cvar2);
    // free(mxmn);

    return 0;
}
