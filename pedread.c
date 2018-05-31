/* modification of matread but operating on words instead of floats */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "pedread.h"

#define GRABNUMGENOSONTHISLINE 1

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

void prtaawcdata(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf("l%ut", aawc->aaw[i]->aw[j]->lp1-1);
            switch(aawc->aaw[i]->aw[j]->t) {
                case AA: case CC: case GG: case TT: printf("Hom "); break;
                default: printf("Het ");
            }
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

void prt2dets(aaw_c *aawc, sampga_t *sga) /* print details of the two datastrucs we have running here */
{
    int i;
    printf("Straight printout of words per line:\n");
    for(i=0;i<aawc->numl;++i) {
        printf("%u ", aawc->aaw[i]->al); 
    }
    printf("totlines=%zu\n", aawc->numl);
    printf("uniform genotype number per sample=%u\n", sga->gasz);
}

void rendsga(aaw_c *aawc, sampga_t *sga) /* render the sga from the aawc */
{
    int i, j, k, m, startaawcidx;
    
    if(aawc->aaw[0]->al != aawc->aaw[1]->al) {
        sga->nsamps = aawc->numl-1;
        startaawcidx = 1;
    } else {
        sga->nsamps = aawc->numl;
        startaawcidx = 0;
    }
    sga->iid=malloc(sga->nsamps*sizeof(char*));
    sga->ga=malloc(sga->nsamps*sizeof(t_t*));
    for(i=startaawcidx;i< aawc->numl;++i) {
        j= i-startaawcidx;
        sga->ga[j]=malloc(sizeof(t_t));
        sga->iid[j]=malloc(aawc->aaw[i]->aw[1]->lp1*sizeof(char));
        strcpy(sga->iid[j], aawc->aaw[i]->aw[1]->w);
        if(sga->tseles) {
            for(k=6;k<aawc->aaw[i]->al-1; k+=2) {
                m=(k-6)/2;
                if( (aawc->aaw[i]->aw[k]->w[0] == 'A') & (aawc->aaw[i]->aw[k+1]->w[0] == 'A') )
                    sga->ga[j][m] = AA;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'C') & (aawc->aaw[i]->aw[k+1]->w[0] == 'C') )
                    sga->ga[j][m] = CC;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'G') & (aawc->aaw[i]->aw[k+1]->w[0] == 'G') )
                    sga->ga[j][m] = GG;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'T') & (aawc->aaw[i]->aw[k+1]->w[0] == 'T') )
                    sga->ga[j][m] = TT;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'A') & (aawc->aaw[i]->aw[k+1]->w[0] == 'C') )
                    sga->ga[j][m] = AC;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'A') & (aawc->aaw[i]->aw[k+1]->w[0] == 'G') )
                    sga->ga[j][m] = AG;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'A') & (aawc->aaw[i]->aw[k+1]->w[0] == 'T') )
                    sga->ga[j][m] = AT;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'C') & (aawc->aaw[i]->aw[k+1]->w[0] == 'A') )
                    sga->ga[j][m] = CA;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'C') & (aawc->aaw[i]->aw[k+1]->w[0] == 'G') )
                    sga->ga[j][m] = CG;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'C') & (aawc->aaw[i]->aw[k+1]->w[0] == 'T') )
                    sga->ga[j][m] = CT;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'G') & (aawc->aaw[i]->aw[k+1]->w[0] == 'A') )
                    sga->ga[j][m] = GA;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'G') & (aawc->aaw[i]->aw[k+1]->w[0] == 'C') )
                    sga->ga[j][m] = GC;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'G') & (aawc->aaw[i]->aw[k+1]->w[0] == 'T') )
                    sga->ga[j][m] = GT;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'T') & (aawc->aaw[i]->aw[k+1]->w[0] == 'A') )
                    sga->ga[j][m] = TA;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'T') & (aawc->aaw[i]->aw[k+1]->w[0] == 'C') )
                    sga->ga[j][m] = TC;
                else if( (aawc->aaw[i]->aw[k]->w[0] == 'T') & (aawc->aaw[i]->aw[k+1]->w[0] == 'G') )
                    sga->ga[j][m] = TG;
            }
        } else {
            for(k=6;k<aawc->aaw[i]->al; k++) {
                m=k-6;
                if( !strncmp(aawc->aaw[i]->aw[k]->w, "A A", 3))
                    sga->ga[j][m] = AA;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "C C", 3))
                    sga->ga[j][m] = CC;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "G G", 3))
                    sga->ga[j][m] = GG;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "T T", 3))
                    sga->ga[j][m] = TT;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "A C", 3))
                    sga->ga[j][m] = AC;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "A G", 3))
                    sga->ga[j][m] = AG;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "A T", 3))
                    sga->ga[j][m] = AT;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "C A", 3))
                    sga->ga[j][m] = CA;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "C G", 3))
                    sga->ga[j][m] = CG;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "C T", 3))
                    sga->ga[j][m] = CT;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "G A", 3))
                    sga->ga[j][m] = GA;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "G C", 3))
                    sga->ga[j][m] = GC;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "G T", 3))
                    sga->ga[j][m] = GT;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "T A", 3))
                    sga->ga[j][m] = TA;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "T C", 3))
                    sga->ga[j][m] = TC;
                else if( !strncmp(aawc->aaw[i]->aw[k]->w, "T G", 3))
                    sga->ga[j][m] = TG;
            }
        }
    }
}

void prtsga(sampga_t *sga) /* print the sga as a ped, no header, alleles are tab separated */
{
    int i, j;
    
    for(i=0;i<sga->nsamps;++i) {
        printf("%s\t", sga->iid[i]);
        for(j=0;j<sga->gasz;++j) {
            switch(sga->ga[i][j]) {
                case AA: printf("AA\t"); break;
                case CC: printf("CC\t"); break;
                case GG: printf("GG\t"); break;
                case TT: printf("TT\t"); break;
                case AC: printf("AC\t"); break;
                case AG: printf("AG\t"); break;
                case AT: printf("AT\t"); break;
                case CA: printf("CA\t"); break;
                case CG: printf("CG\t"); break;
                case CT: printf("CT\t"); break;
                case GA: printf("GA\t"); break;
                case GC: printf("GC\t"); break;
                case GT: printf("GT\t"); break;
                case TA: printf("TA\t"); break;
                case TC: printf("TC\t"); break;
                case TG: printf("TG\t"); break;
            }
        }
        printf("\n"); 
    }
}

void prtaawcsum0(aaw_c *aawc) /* a summary, or one way of doing one */
{
    int i;
    printf("Straight printout of words per line:\n");
    for(i=0;i<aawc->numl;++i) {
        printf("%u ", aawc->aaw[i]->al); 
    }
    printf("totlines=%zu\n", aawc->numl);
}

aaw_c *processinpf(char *fname, sampga_t *sga)
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
        // if( (c== '\n') | (c == '\t') ) { /* in ped files, genotypes spearated by tabs and each snp separated by spaces. */
            if( inword==1) { /* cue word-ending procedure */
                if( (aawc->numl == GRABNUMGENOSONTHISLINE) & (couw == 6) ) {
                    if(couc == 1) /* i.e. our tab sep tokens are only 1 char long */
                        sga->tseles=1;
                }
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
                if(aawc->numl == GRABNUMGENOSONTHISLINE) {
                    if(sga->tseles)
                        sga->gasz = (couw-6)/2;
                    else
                        sga->gasz = couw-6;
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
    if(aawc->numl ==1) {
        printf("Warning: this PED file only has one line. This program uses the second\n"); 
        printf("line to test existence of header, so this is not a very suitable target file.\n");
    }

    return aawc;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=2) {
        printf("Error. Pls supply argument (name of text file).\n");
        exit(EXIT_FAILURE);
    }
    int i;
    sampga_t *sga=malloc(sizeof(sampga_t));
#ifdef DBG2
    printf("typeszs: aaw_c: %zu aw_c: %zu w_c: %zu\n", sizeof(aaw_c), sizeof(aw_c), sizeof(w_c));
#endif

    aaw_c *aawc=processinpf(argv[1], sga);
#ifdef DBG
    prtaawcdbg(aawc);
#else
    // prtaawcsum0(aawc);
    prt2dets(aawc, sga);
#endif

    /* let's see if the rendga works properly */
    rendsga(aawc, sga);
    free_aawc(&aawc); /* finished with this chappie now */
//    prtsga(sga);

    /* now free the sga sans function */
    for(i=0;i<sga->nsamps;++i) {
        free(sga->iid[i]);
        free(sga->ga[i]);
    }
    free(sga->iid);
    free(sga->ga);
    free(sga);
    return 0;
}
