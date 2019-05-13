/* modification of matread but operating on words instead of floats */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "namsets.h"

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

snodm **hashnam(aaw_c *aawc, unsigned tsz, unsigned *dgcou, trival tf /* target file is first=1 or second=2 */)
{
    unsigned i;

    snodm **stab=malloc(tsz*sizeof(snodm *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL;
    snodm *tsnod0, *tsnod2;

    unsigned tint;
    unsigned dpgrps=0;
    unsigned dupslo=0; // dups less original, I.e. one will be retained.

    /* OK, we're going to loop through the map file container: i here follows the global SNP name index */
    for(i=0; i< aawc->numl; ++i) {

        tint=hashit(aawc->aaw[i]->aw[0]->w, tsz); // hash the snpname
        if( (stab[tint] == NULL) ) { // nothing in that slot right now.
            stab[tint]=malloc(sizeof(snodm));
            stab[tint]->aw=aawc->aaw[i];
            stab[tint]->idx=i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            if(!strcmp(tsnod2->aw->aw[0]->w, aawc->aaw[i]->aw[0]->w)) {
                // printf("Duppair: %s vs %s\n", tsnod2->aw->aw[0]->w, aawc->aaw[i]->aw[0]->w);
                if(!tsnod2->aw->dup) {
                    dpgrps++;
                    tsnod2->aw->dup=1;
                    tsnod2->aw->dgrp = dpgrps;
                }
                aawc->aaw[i]->dup=2;
                dupslo++;
                aawc->aaw[i]->dgrp=tsnod2->aw->dgrp;
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
    *dgcou = dpgrps;
    return stab;
}

void mu_nam(aaw_c *aawc, snodm **stam, unsigned tsz, trival sf /* subject file, either 1 or 2 */)
{
    unsigned i;
    snodm *tsnod0, *tsnod2;
    unsigned tint;

    for(i=0; i< aawc->numl; ++i) {

        tint=hashit(aawc->aaw[i]->aw[0]->w, tsz); // hash the snpname
        if( (stam[tint] == NULL) )
            continue; // ma2 stays at zero ... there is no match for this

        tsnod2=stam[tint];
        while( (tsnod2 != NULL) ) {
            if(!strcmp(tsnod2->aw->aw[0]->w, aawc->aaw[i]->aw[0]->w)) {
                if(sf==1)
                    aawc->aaw[i]->ma2=1;
                else if(sf==2)
                    aawc->aaw[i]->ma1=1;
                break;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
    }
    return;
}

void rep_mu(aaw_c *aawc, aaw_c *aawc2)
{
    unsigned i;
    unsigned resarray[3] = {0};

    for(i=0; i< aawc->numl; ++i) {
        if(aawc->aaw[i]->ma2)
            resarray[1]++;
        else
            resarray[0]++;
    }

    for(i=0; i< aawc2->numl; ++i)
        if(aawc2->aaw[i]->ma1)
            resarray[2]++;

    printf("TOTF1=%lli; #F1\\F2=%u; #F1nF2=%u; #F2\\F1=%u TOTF2=%lli\n", aawc->numl, resarray[0], resarray[1], resarray[2], aawc2->numl);
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
            printf("%s ", tsnod2->aw->aw[0]->w);
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
            printf((j!=aawc->aaw[i]->al-1)?"%s\\":"%s\n", aawc->aaw[i]->aw[j]->w);
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
    boole intitle=0;
    unsigned lbuf=LBUF /* buffer for number of lines */, cbuf=CBUF /* char buffer for size of w_c's: reused for every word */;
    aaw_c *aawc=crea_aawc(lbuf); /* array of words per line */

    while( (c=fgetc(fp)) != EOF) {
        if(c=='>')
            intitle=1;
        else if((intitle) & (oldc=='\n'))
            intitle=0;
        if((c==' ') & !intitle)
            continue;
        if( ((c==' ') & intitle) | (c=='"') | (c=='=') | (c=='|') | (c== '\n') | (c == '\t') ) {
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

aaw_c *processinpf1l(char *fname)
{
    /* declarations */
    FILE *fp=fopen(fname,"r");
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */;
    int c;
    boole inword=0;
    unsigned lbuf=LBUF /* buffer for number of lines */, cbuf=CBUF /* char buffer for size of w_c's: reused for every word */;
    aaw_c *aawc=crea_aawc(lbuf); /* array of words per line */

    while( (c=fgetc(fp)) != EOF) {
        if(c=='\n') { /* cue line-ending procedure */
            if( inword==1) { /* cue word-ending procedure */
                aawc->aaw[aawc->numl]->aw[couw]->w[couc++]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=couc;
                norm_wc(aawc->aaw[aawc->numl]->aw+couw);
                couw++; /* verified: this has to be here */
            }
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
    printf("Program \"namsets\" takes two name files (one name per line)\n"); 
    printf("And compares them giving short report\n"); 
    printf("Arg1 is file 1, arg2 is file2.\n"); 
    printf("> If any of these files contain duplicates the comparison is not really \"fair\"\n");
    printf("  To see the duplicates include option \"-d\" after file1 and file2 args.\n"); 
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc<3)
        prtusage();

    int argignore=2; 
    int oargc=argc-argignore;
    char **oargv=argv+argignore;
    optstruct opstru={0};
    catchopts(&opstru, oargc, oargv);

    unsigned dgf1, dgf2; // dup groups in file 1, in file 2
    aaw_c *aawc2=processinpf1l(argv[2]); // yes second argument is first to process because it's the target.
    unsigned htsz=givehtsz(aawc2->numl);
    snodm **stam = hashnam(aawc2, htsz, &dgf2, 2);
    prtchaharr(stam, htsz);
    if(dgf2) {
        fprintf(stderr, "File_%i duplicates groups= %u\n", 2, dgf2);
        if(opstru.dflag)
            rep_dups(stam, htsz);
    } else
        fprintf(stderr, "File_%i fine, no duplicates\n", 2);

    // prtchaharr(stam, htsz);
    aaw_c *aawc=processinpf1l(argv[1]);
    mu_nam(aawc, stam, htsz, 1);
    freechainharr(stam, htsz);

    // OK, second run
    htsz=givehtsz(aawc->numl);
    stam = hashnam(aawc, htsz, &dgf1, 1);
    if(dgf1) {
        fprintf(stderr, "File_%i duplicates groups= %u\n", 1, dgf1);
        if(opstru.dflag)
            rep_dups(stam, htsz);
    } else
        fprintf(stderr, "File_%i fine, no duplicates\n", 2);

    mu_nam(aawc2, stam, htsz, 2);
    freechainharr(stam, htsz); // yep I can do this now.

    rep_mu(aawc, aawc2);

    free_aawc(&aawc);
    free_aawc(&aawc2);

    return 0;
}
