/* ENTREZID to SYMBOL converter */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "namgenes2.h"

int catchopts(optstruct *opstru, int oargc, char **oargv)
{
    int c;
    opterr = 0;
    while ((c = getopt (oargc, oargv, "d")) != -1)
        switch (c) {
            case 'd': // want to see the genotypes of these SNPs
                opstru->dflag = 1; break;
            default:
                abort();
        }
    return 0;
}

unsigned givehtsz(unsigned mnf)
{
    unsigned htsz=2*mnf/3;
    // try to grab a prime ... well just avoid 5-multiples, 3-multiples, and evens
    if(!(htsz%5)) 
        htsz++; // incrment by 1 if multiple of 5
    if(!(htsz%3)) 
        htsz++;
    if(!(htsz%2)) 
        htsz++;
    return htsz;
}

unsigned hashit(char *str, unsigned tsz) /* Dan Bernstein's one */
{
    unsigned long hash = 5381;
    int c;

    char *tstr=str;
    while ((c = *tstr++))
        hash = ((hash << 5) + hash) + c; /*  hash * 33 + c */

    return hash % tsz;
}

snodm **hashnam(aaw_c *aawc, unsigned tsz, int hi/*hash index*/)
{
    unsigned i;

    snodm **stab=malloc(tsz*sizeof(snodm *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL;
    snodm *tsnod0, *tsnod2;

    unsigned tint;

    /* OK, we're going to loop through the map file container: i here follows the global SNP name index */
    for(i=0; i< aawc->numl; ++i) {
        tint=hashit(aawc->aaw[i]->aw[hi]->w, tsz); // hash the entrez id
        if( (stab[tint] == NULL) ) { // nothing in that slot right now.
            stab[tint]=malloc(sizeof(snodm));
            stab[tint]->aw=aawc->aaw[i];
            stab[tint]->idx=i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            if(!strcmp(tsnod2->aw->aw[hi]->w, aawc->aaw[i]->aw[hi]->w)) {
                // printf("Duppair: %s vs %s\n", tsnod2->aw->aw[0]->w, aawc->aaw[i]->aw[0]->w);
                aawc->aaw[i]->dup=2;
                break;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        if(aawc->aaw[i]->dup)
            continue;
        tsnod0->n=malloc(sizeof(snodm));
        tsnod0->n->aw = aawc->aaw[i];
        tsnod0->n->idx=i;
        tsnod0->n->n=NULL;
    }
    return stab;
}

void mu_nam0(aaw_c *aawc, snodm **stam, unsigned tsz, int hi)
{
    /* if there's no match, don't print anything */
    unsigned i, j;
    snodm *tsnod0, *tsnod2;
    unsigned tint;
    boole ma;
    int nummisses=0;

    for(i=0; i< aawc->numl; ++i) {
        printf("%s: ", aawc->aaw[i]->aw[0]->w);
        for(j=1;j<aawc->aaw[i]->al;++j) {
            ma=0;
            tint=hashit(aawc->aaw[i]->aw[j]->w, tsz); // hash the snpname
            if( (stam[tint] == NULL) ) {
                if(j==aawc->aaw[i]->al-1) // if it happens to be the last in hte list
                    putchar('\n');
                nummisses++;
                continue; // blank entry, no match, so skip.
            }

            tsnod2=stam[tint];
            while( (tsnod2 != NULL) ) {
                if(!strcmp(tsnod2->aw->aw[hi]->w, aawc->aaw[i]->aw[j]->w)) {
                    printf((j==aawc->aaw[i]->al-1)?"%s\n":"%s ", tsnod2->aw->aw[1]->w);
                    ma=1;
                    break;
                }
                tsnod0=tsnod2;
                tsnod2=tsnod2->n;
            }
            if(!ma) {
                nummisses++;
                if(j==aawc->aaw[i]->al-1) // if it happens to be the last in hte list
                    putchar('\n'); // nothing else required.
            }
        }
    }
    printf("Nummisses = %i\n", nummisses); 
    return;
}

void mu_nam2(aaw_c *aawc, snodm **stam, unsigned tsz, int hi)
{
    /* if there's no match, don't print anything */
    unsigned i, j;
    snodm *tsnod0, *tsnod2;
    unsigned tint;
    boole ma;
    int nummisses=0;

    for(i=0; i< aawc->numl; ++i) {
        printf("%s: ", aawc->aaw[i]->aw[0]->w);
        for(j=1;j<aawc->aaw[i]->al;++j) {
            ma=0;
            tint=hashit(aawc->aaw[i]->aw[j]->w, tsz); // hash the snpname
            if( (stam[tint] == NULL) ) {
                if(j==aawc->aaw[i]->al-1) // if it happens to be the last in hte list
                    putchar('\n');
                nummisses++;
                continue; // blank entry, no match, so skip.
            }

            tsnod2=stam[tint];
            while( (tsnod2 != NULL) ) {
                if(!strcmp(tsnod2->aw->aw[hi]->w, aawc->aaw[i]->aw[j]->w)) {
                    printf((j==aawc->aaw[i]->al-1)?"%s\n":"%s ", tsnod2->aw->aw[1]->w);
                    ma=1;
                    aawc->aaw[i]->aw[j]->ma2 = 1;
                    break;
                }
                tsnod0=tsnod2;
                tsnod2=tsnod2->n;
            }
            if(!ma) {
                nummisses++;
                if(j==aawc->aaw[i]->al-1) // if it happens to be the last in hte list
                    putchar('\n'); // nothing else required.
            }
        }
    }
    printf("Nummisses = %i\n", nummisses); 
    return;
}

void mu_nam22f(aaw_c *aawc, snodm **stam, unsigned tsz, int hi, snodm **stam1, unsigned tsz1, int hi1, char *fn)
{
    /* outputting to file */
    unsigned i, j;
    snodm *tsnod0, *tsnod2;
    unsigned tint, tint1;
    boole ma0, ma;
    int nummisses=0, nummisses1=0;
    FILE *fp=fopen(fn,"w");

    for(i=0; i< aawc->numl; ++i) {
        // fprintf(fp, "%s: ", aawc->aaw[i]->aw[0]->w);
        for(j=1;j<aawc->aaw[i]->al;++j) {
            ma=0;
            tint=hashit(aawc->aaw[i]->aw[j]->w, tsz); // hash the snpname
            if( (stam[tint] == NULL) ) {
                tint1=hashit(aawc->aaw[i]->aw[j]->w, tsz1); // hash the snpname
                nummisses++;
                if( (stam1[tint1] == NULL) ) {
                    nummisses1++;
                    if(j==aawc->aaw[i]->al-1) // if it happens to be the last in hte list
                        fprintf(fp, "\n");
                }
                continue; // blank entry, no match, so skip.
            }

            tsnod2=stam[tint];
            while( (tsnod2 != NULL) ) {
                if(!strcmp(tsnod2->aw->aw[hi]->w, aawc->aaw[i]->aw[j]->w)) {
                    fprintf(fp, (j==aawc->aaw[i]->al-1)?"%s\n":"%s ", tsnod2->aw->aw[1]->w);
                    ma0=1;
                    aawc->aaw[i]->aw[j]->ma1 = 1;
                    break;
                }
                tsnod0=tsnod2;
                tsnod2=tsnod2->n;
            }
            if(!ma0) {
                nummisses++;
                tsnod2=stam1[tint1];
                while( (tsnod2 != NULL) ) {
                    if(!strcmp(tsnod2->aw->aw[hi1]->w, aawc->aaw[i]->aw[j]->w)) {
                        fprintf(fp, (j==aawc->aaw[i]->al-1)?"%s\n":"%s ", tsnod2->aw->aw[1]->w);
                        ma=1;
                        aawc->aaw[i]->aw[j]->ma2 = 1;
                        break;
                    }
                    tsnod0=tsnod2;
                    tsnod2=tsnod2->n;
                }
                if(!ma) {
                    nummisses1++;
                    if(j==aawc->aaw[i]->al-1) // if it happens to be the last in hte list
                        fprintf(fp, "\n");
                }
            }
        }
    }
    fclose(fp);
    return;
}

void mu_nam22_(aaw_c *aawc, snodm **stam, unsigned tsz, int hi, snodm **stam1, unsigned tsz1, int hi1)
{
    /* if there's no matchin ht efirst has, try the second */
    /* _ no print version */
    unsigned i, j;
    snodm *tsnod0, *tsnod2;
    unsigned tint, tint1;
    boole ma0, ma;
    int nummisses=0, nummisses1=0;

    for(i=0; i< aawc->numl; ++i) {
        for(j=1;j<aawc->aaw[i]->al;++j) {
            ma=0;
            tint=hashit(aawc->aaw[i]->aw[j]->w, tsz); // hash the snpname
            if( (stam[tint] == NULL) ) {
                tint1=hashit(aawc->aaw[i]->aw[j]->w, tsz1); // hash the snpname
                nummisses++;
                printf("No entry in HT1 for: %s\n", aawc->aaw[i]->aw[j]->w); 
                if( (stam1[tint1] == NULL) ) {
                    printf("No entry in HT2 for: %s\n", aawc->aaw[i]->aw[j]->w); 
                    nummisses1++;
                }
                continue; // blank entry, no match, so skip.
            }

            tsnod2=stam[tint];
            while( (tsnod2 != NULL) ) {
                if(!strcmp(tsnod2->aw->aw[hi]->w, aawc->aaw[i]->aw[j]->w)) {
                    ma0=1;
                    aawc->aaw[i]->aw[j]->ma1 = 1;
                    break;
                }
                tsnod0=tsnod2;
                tsnod2=tsnod2->n;
            }
            if(!ma0) {
                printf("HT1 entries exhausted, no match for: %s\n", aawc->aaw[i]->aw[j]->w); 
                nummisses++;
                tsnod2=stam1[tint1];
                while( (tsnod2 != NULL) ) {
                    if(!strcmp(tsnod2->aw->aw[hi1]->w, aawc->aaw[i]->aw[j]->w)) {
                        ma=1;
                        aawc->aaw[i]->aw[j]->ma2 = 1;
                        break;
                    }
                    tsnod0=tsnod2;
                    tsnod2=tsnod2->n;
                }
                if(!ma) {
                    printf("HT2 entries exhausted, no match for: %s\n", aawc->aaw[i]->aw[j]->w); 
                    nummisses1++;
                }
            }
        }
    }
    return;
}

void mu_nam(aaw_c *aawc, snodm **stam, unsigned tsz)
{
    /* where no matched, put NOENSEMBL */
    unsigned i, j;
    snodm *tsnod0, *tsnod2;
    unsigned tint;
    boole ma;

    for(i=0; i< aawc->numl; ++i) {
        printf("%s: ", aawc->aaw[i]->aw[0]->w);
        for(j=1;j<aawc->aaw[i]->al;++j) {
            ma=0;
            tint=hashit(aawc->aaw[i]->aw[j]->w, tsz); // hash the snpname
            if( (stam[tint] == NULL) ) {
                printf((j==aawc->aaw[i]->al-1)? "NOENSEMBL\n": "NOENSEMBL ");
                continue; // ma2 stays at zero ... there is no match for this
            }

            tsnod2=stam[tint];
            while( (tsnod2 != NULL) ) {
                if(!strcmp(tsnod2->aw->aw[2]->w, aawc->aaw[i]->aw[j]->w)) {
                    printf((j==aawc->aaw[i]->al-1)?"%s\n":"%s ", tsnod2->aw->aw[1]->w);
                    ma=1;
                    break;
                }
                tsnod0=tsnod2;
                tsnod2=tsnod2->n;
            }
            if(!ma)
                printf((j==aawc->aaw[i]->al-1)? "NOENSEMBL\n": "NOENSEMBL ");
        }
    }
    return;
}

void freechainharr(snodm **stam, size_t tsz)
{
    int i;
    snodm *tsnod0, *tsnod2;
    for(i=0; i<tsz; ++i) {
        if( (stam[i] != NULL) ) {
            while( (stam[i]->n != NULL) ) {
                tsnod0=stam[i];
                tsnod2=stam[i]->n;
                while((tsnod2->n != NULL) ){
                    tsnod0=tsnod2;
                    tsnod2=tsnod2->n;
                }
                free(tsnod0->n);
                tsnod0->n=NULL;
            }
            free(stam[i]);
        }
    }
    free(stam);
    return;
}

void prtchaharr(snodm **stam, unsigned tsz)
{
    unsigned i;
    snodm *tsnod2;
    for(i=0;i<tsz;++i) {
        printf("Tablepos %i: ", i); 
        tsnod2=stam[i];
        while(tsnod2) {
            printf("%s ", tsnod2->aw->aw[2]->w);
            tsnod2=tsnod2->n;
        }
        putchar('\n');
    }
    return;
}

void rep_dups(snodm **stam, unsigned tsz)
{
    unsigned i;
    snodm *tsnod2;
    boole dupseen;

    for(i=0;i<tsz;++i) {
        dupseen=0;
        tsnod2=stam[i];
        while(tsnod2) {
            if(tsnod2->aw->dup == 1) {
                printf("%s (Headgrp %u) ", tsnod2->aw->aw[0]->w, tsnod2->aw->dgrp);
                dupseen=1;
            } else if(tsnod2->aw->dup == 2) {
                printf("%s (Membgrp %u)", tsnod2->aw->aw[0]->w, tsnod2->aw->dgrp);
                dupseen=1;
            }

            tsnod2=tsnod2->n;
        }
        if(dupseen)
            putchar('\n');
    }
    return;
}

w_c *crea_wc(unsigned initsz)
{
    w_c *wc=malloc(sizeof(w_c));
    wc->lp1=initsz;
    wc->ma1=0;
    wc->ma2=0;
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
    awc->dup=0;
    awc->dgrp=0;
    awc->ma2=0;
    awc->ma1=0;
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

void prtaawcplain(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L)%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j)
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
    }
}

void prtaawcplain2(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L)%u(%uw) word is:", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j)
            printf("%s\n", aawc->aaw[i]->aw[2]->w);
    }
}

aaw_c *processinpf0(char *fname)
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

void prtusage()
{
    printf("Program \"namgenes2\" takes 4 args: 2 name tsv files 1) a 4 column one, where index 2 is entrezid, 2) a two column one where index 0 is entrezid 3) the geneset file and 4) output file\n");
    printf("Example usage: ./namgenes2 genes2.tsv id2sy.tsv eids.txt oo.out\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=5)
        prtusage();

    aaw_c *aawc=processinpf0(argv[1]);
    aaw_c *aawc1=processinpf0(argv[2]);
    aaw_c *aawc2=processinpf0(argv[3]);
    int i, j;

    unsigned htsz=givehtsz(aawc->numl);
    snodm **stam = hashnam(aawc, htsz, 2);
    // mu_nam0(aawc2, stam, htsz, 2);

    unsigned htsz1=givehtsz(aawc1->numl);
    snodm **stam1 = hashnam(aawc1, htsz1, 0);
    // mu_nam0(aawc2, stam1, htsz1, 0);
    mu_nam22f(aawc2, stam, htsz, 2, stam1, htsz1, 0, argv[4]);

    // check which id's did not match.
    // for(i=0; i< aawc2->numl; ++i) {
    //     for(j=1;j<aawc2->aaw[i]->al;++j)
    //         if((!aawc2->aaw[i]->aw[j]->ma1) && (!aawc2->aaw[i]->aw[j]->ma2))
    //             printf("%s ", aawc2->aaw[i]->aw[j]->w);
    // }
    // printf("\n"); 

    freechainharr(stam, htsz);
    freechainharr(stam1, htsz1);
    free_aawc(&aawc);
    free_aawc(&aawc1);
    free_aawc(&aawc2);

    return 0;
}
