/* modification of matread but operating on words instead of floats */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define CBUF 8
#define GBUF 8
#define WBUF 8

typedef unsigned char boole;

typedef struct /* word type */
{
    char *w;
    unsigned b; /* buffer */
    unsigned lp1; /* length */
} w_t;

w_t *creawt(void)
{
    w_t *wt=malloc(sizeof(w_t));
    wt->b=CBUF;
    wt->lp1=0;
    wt->w=malloc(wt->b*sizeof(char));
    return wt;
}

void reallwt(w_t **wt)
{
    w_t *twt=*wt;
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
    free(twt->w)
    free(twt)
    return;
}

typedef struct /* wa_t: word array */
{
    w_t **wa;
    unsigned ab;
    unsigned al;
} wa_t;

wa_t *creatwat(void)
{
    int i;
    wa_t *wat=malloc(sizeof(wa_t));
    wat->ab=WABUF;
    wat->al=0;
    for(i=0;i<wat->ab;++i) 
        wat->wa[i]=creawt();
    return wat;
}

void reallwat(wa_t **wat, unsigned buf)
{
    int i;
    wa_t *twat=*wat;
    for(i=twat->ab-buf;i<twat->ab;++i) 
        twat->wa[i]=creawt();
    *wat=twat;
    return;
}

void normwat(wa_t **wat)
{
    int i;
    wa_t *twat=*wat;
    for(i=twat->al;i<twat->ab;++i) 
        freewt(twat->wa[i]);
    *wat=twat;
    return;
}

typedef struct /* wseq_t */
{
    size_t *wln;
    size_t wsbuf;
    size_t quan;
    size_t lbuf; /* a buffer for the number of lines */
    size_t numl; /* number of lines, i.e. rows */
    size_t *wpla; /* words per line array: the number of words on each line */
    wa_t *wat;
} wseq_t;

wseq_t *create_wseq_t(size_t initsz)
{
    wseq_t *words=malloc(sizeof(wseq_t));
    words->wsbuf = initsz;
    words->quan = initsz;
    words->wln=calloc(words->wsbuf, sizeof(size_t));
    words->lbuf=WBUF;
    words->numl=0;
    words->wpla=calloc(words->lbuf, sizeof(size_t));
    return words;
}

void free_wseq(wseq_t *wa)
{
    free(wa->wln);
    free(wa->wpla);
    free(wa);
}

float *processinpf(char *fname, int *m, int *n)
{
    /* declarations */
    FILE *fp=fopen(fname,"r");
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
    int c;
    boole inword=0;
    wseq_t *wa=create_wseq_t(GBUF);
    size_t bwbuf=WBUF;

    wa_t *wa=

    for(i=0;i<n;++i) 
    w_t *bufworda=calloc(bwbuf, sizeof(w_t)); /* this is the string we'll keep overwriting. */

    float *mat=malloc(GBUF*sizeof(float));

    while( (c=fgetc(fp)) != EOF) {
        /*  take care of  */
        if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) {
            if( inword==1) { /* we've been in a word so we have to end it */
                wa->wln[couw]=couc;
                bufword[k][couc++]='\0';
                bufword[k] = realloc(bufword[k], couc*sizeof(char)); /* normalize */
                mat[couw]=atof(bufword);
                couc=0;
                couw++;
            }
            if(c=='\n') {
                if(wa->numl == wa->lbuf-1) {
                    wa->lbuf += WBUF;
                    wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
                    memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
                }
                wa->wpla[wa->numl] = couw-oldcouw;
                oldcouw=couw;
                wa->numl++;
            }
            inword=0;
        } else if(inword==0) {
            if(couw == wa->wsbuf-1) {
                wa->wsbuf += GBUF;
                wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
                mat=realloc(mat, wa->wsbuf*sizeof(float));
                for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
                    wa->wln[i]=0;
            }
            couc=0;
            bwbuf=WBUF;
            bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
            bufword[couc++]=c; /* no need to check here, it's the first character */
            inword=1;
        } else if( (c == 0x2E) | ((c >= 0x30) && (c <= 0x39)) ) {
            if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
                bwbuf += WBUF;
                bufword = realloc(bufword, bwbuf*sizeof(char));
            }
            bufword[couc++]=c;
        } else {
            printf("Error. Non-float character detected. This program is only for reading floats\n"); 
            free_wseq(wa);
            exit(EXIT_FAILURE);
        }

    } /* end of big for statement */
    fclose(fp);
    free(bufword);

    /* normalization stage */
    wa->quan=couw;
    wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
    mat = realloc(mat, wa->quan*sizeof(float)); /* normalize */
    wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

    *m= wa->numl;
    int k=wa->wpla[0];
    for(i=1;i<wa->numl;++i)
        if(k != wa->wpla[i])
            printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
    *n= k; 
    free_wseq(wa);

    return mat;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=2) {
        printf("Error. Pls supply argument (name of text file).\n");
        exit(EXIT_FAILURE);
    }

    int i, j, m, n;
    float *mat=processinpf(argv[1], &m, &n);

    printf("Matrix is %i rows by %i columns and is as follows:\n", m, n); 
    for(i=0;i<m;++i) {
        for(j=0;j<n;++j) 
            printf("%f ", mat[i*n+j]);
        printf("\n"); 
    }

    free(mat);

    return 0;
}
