/* some statistics from a ped file ... specificaly fro Illumina BeadArray chip target */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "pedsta0.h"

w_c *crea_wc(unsigned initsz)
{
    w_c *wc=malloc(sizeof(w_c));
    wc->lp1=initsz;
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
    printf("Legend: Begin with Line number, numwords in brackets, words uscore number, length of words:\n\n");
    for(i=0;i<aawc->numl;++i) {
        printf("l.%u(%u): ", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf("w_%u: ", j); 
            // if(aawc->aaw[i]->aw[j]->t == NUMS) {
            //     printf("NUM! "); 
            //     continue;
            // } else if(aawc->aaw[i]->aw[j]->t == PNI) {
            //     printf("PNI! "); 
            //     continue;
            // } else if(aawc->aaw[i]->aw[j]->t == STCP) {
            //     printf("STCP! "); 
            //     continue;
            // }
            for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                putchar(aawc->aaw[i]->aw[j]->w[k]);
            printf("/%u ", aawc->aaw[i]->aw[j]->lp1-1); 
        }
        printf("\n"); 
    }
}

void prtaawcdbg2(aaw_c *aawc)
{
    int j, k;
    size_t i;
    printf("Legend: Begin with Line number, numwords in brackets, words uscore number, length of words:\n\n");
    for(i=0;i<aawc->numl;++i) {
        if(aawc->aaw[i]->aw[0]->w[0] == '#')
            continue;
        printf("%zu) IID: %s GT:", i, aawc->aaw[i]->aw[1]->w);
        for(j=6;j<aawc->aaw[i]->al;j+=2) {
            putchar(' ');
            for(k=0;k<aawc->aaw[i]->aw[j]->lp1-1; k++)
                putchar(aawc->aaw[i]->aw[j]->w[k]);
            putchar('/');
            for(k=0;k<aawc->aaw[i]->aw[j+1]->lp1-1; k++)
                putchar(aawc->aaw[i]->aw[j+1]->w[k]);
        }
        printf("\n");
    }
}

void statsaawc(aaw_c *aawc)
{
    char a1, a2;
    int j, j2, called, unc, totgts;
    size_t i;

    int **cous=malloc(aawc->numl*sizeof(int*));
    for(i=0;i<aawc->numl;++i) 
        cous[i]=calloc(NTS, sizeof(int));

    t_t **spotts=malloc(aawc->numl*sizeof(t_t*));
    for(i=0;i<aawc->numl;++i) 
        spotts[i]=malloc(((aawc->aaw[i]->al-6)/2)*sizeof(t_t));

    for(i=0;i<aawc->numl;++i) {
        for(j=6;j<aawc->aaw[i]->al;j+=2) {
            a1 = aawc->aaw[i]->aw[j]->w[0];
            a2 = aawc->aaw[i]->aw[j+1]->w[0];
            j2=(j-6)/2;
            switch(a1) {
                case 'A': case 'C': case 'G': case 'T':
                    if(a1 == a2) {
                        spotts[i][j2]=HMF;
                        cous[i][0]++;
                    } else {
                        switch(a2) {
                            case 'A': case 'C': case 'G': case 'T':
                                spotts[i][j2]=HTF;
                                cous[i][1]++; break;
                            case 'I': case 'D':
                                spotts[i][j2]=IGT;
                                cous[i][2]++; break;
                            case '0': case 'N': default:
                                spotts[i][j2]=ZGT;
                                cous[i][3]++; break;
                        }
                    }
                    break;
                case 'I': case 'D':
                    spotts[i][j2]=IGT;
                    cous[i][2]++; break;
                case '0': case 'N': default:
                    spotts[i][j2]=ZGT;
                    cous[i][3]++; break;
            }
        }
        unc=cous[i][2]+cous[i][3];
        called=cous[i][0]+cous[i][1];
        totgts=unc+called;
        printf("Sample %s rates: HMFs= %4.4f HTFs=%4.4f GT Call rate: %4.4f\n", aawc->aaw[i]->aw[1]->w, (float)cous[i][0]/totgts, (float)cous[i][1]/totgts, (float)called/totgts);
    }
    /*free up */
    for(i=0;i<aawc->numl;++i) {
        free(cous[i]);
        free(spotts[i]);
    }
    free(cous);
    free(spotts);

    return;
}

void statsaawc2(aaw_c *aawc, float *allra1, float *allra2, size_t *sumaall, int *retnumsamps)
{
    char a1, a2;
    int j;
    size_t i;
    size_t a[10]; // counts which are "of interest": 
    size_t suma1all_, suma2all_;
    float ra1, ra2;

    int numsamps = *retnumsamps; // because numl is not always numsamps due to # comment lines.
    for(i=0;i<aawc->numl;++i) {
        if(aawc->aaw[i]->aw[0]->w[0] == '#')
            continue;
        numsamps++;
        for(j=0;j<10;++j) 
            a[j]= 0;
        // OK now so now we have a match on the IIDs
        // check num genos at very least
        printf("%s\t", aawc->aaw[i]->aw[1]->w);
        for(j=6;j<aawc->aaw[i]->al;j+=2) {
            a1 = aawc->aaw[i]->aw[j]->w[0];
            a2 = aawc->aaw[i]->aw[j+1]->w[0];
            switch(a1) {
                case 'A': case 'C': case 'G': case 'T': case 'a': case 'c': case 'g': case 't':
                    a[0]++; break;
                case '0':
                    a[2]++; break;
                case 'N': case 'n':
                    a[4]++; break;
                case 'I': case 'D': case 'i': case 'd':
                    a[6]++; break;
                default:
                    a[8]++;
            }
            switch(a2) {
                case 'A': case 'C': case 'G': case 'T': case 'a': case 'c': case 'g': case 't':
                    a[1]++; break;
                case '0':
                    a[3]++; break;
                case 'N': case 'n':
                    a[5]++; break;
                case 'I': case 'D': case 'i': case 'd':
                    a[7]++; break;
                default:
                    a[9]++;
            }
        }
        for(j=0;j<10;++j) 
            printf("%zu\t", a[j]);
        suma1all_=a[0]+a[2]+a[4]+a[6]+a[8];
        suma2all_=a[1]+a[3]+a[5]+a[7]+a[9];

        ra1=(float)(a[0]+a[6])/suma1all_;
        printf("%4.4f\t", ra1);
        ra2=(float)(a[1]+a[7])/suma2all_;
        printf("%4.4f\n", ra2);
        /* now how we include IDs because hese are properly called. Often we don't used them though */
        *allra1 += ra1;
        *allra2 += ra2;
        *sumaall=suma1all_+suma2all_;
        *retnumsamps=numsamps;
    }
}

void prtaawcdata(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf("l%ut", aawc->aaw[i]->aw[j]->lp1-1);
        }
    }
    printf("\n"); 
}

void prtaawcplain(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j)
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
    }
}

aaw_c *processinpf(FILE *fp)
{
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */;
    int c;
    boole inword=0;
    unsigned lbuf=LBUF /* buffer for number of lines */, cbuf=CBUF /* char buffer for size of w_c's: reused for every word */;
    aaw_c *aawc=crea_aawc(lbuf); /* array of words per line */

    for(;;) {
        c=fgetc(fp);
        if(c==EOF)
            goto uit;
        if( (c== '\n') | (c == ' ') | (c == '\t') ) {
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
    } /* end of big for statement */

uit:

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
        printf("Error. Pls supply one pedfile as argument: stats will be generated.\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp=fopen(argv[1],"r");
    aaw_c *aawc=NULL;

    aawc=processinpf(fp);
    statsaawc(aawc);
    free_aawc(&aawc);
    fclose(fp);

    return 0;
}
