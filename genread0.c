/* modification of matread but operating on words instead of floats */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "genread0.h"

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

void convt0(char *tmng, float dt)
{
    int hs, ms, ss, fs; // fractions of secs
    char h[3]={0}; h[0] = tmng[0]; h[1] = tmng[1];
    hs=atoi(h);
    char m[3]={0}; m[0] = tmng[3]; m[1] = tmng[4];
    ms=atoi(m);
    char s[3]={0}; s[0] = tmng[6]; s[1] = tmng[7];
    ss=atoi(s);
    char f[4]={0}; f[0] = tmng[9]; f[1] = tmng[10], f[2] = tmng[11];
    fs=atoi(f);
    // printf("hs=%i ms=%i ss=%i fs=%i\n", hs, ms, ss, fs); 
    // printf("%s --- %i:%i:%i,%i\n", tmng, hs, ms, ss, fs); 
    float t0=hs*3600+ms*60+ss+fs/1000.;
    // printf("%4.4f to %4.4f\n", t0, t0+dt); 
    float newt=t0+dt; 
    hs=(int)newt/3600;
    int hsm=hs*3600;
    ms=(int)(newt-hsm)/60;
    int msm=ms*60;
    ss=(int)(newt-hsm-msm);
    fs=100*(newt-hsm-msm-ss);
    fs*=10;
    printf("%s --- %02i:%02i:%02i,%03i\n", tmng, hs, ms, ss, fs); 
}


void convt(char *tmng, float dt)
{
    int hs, ms, ss, fs; // fractions of secs
    char h[3]={0}; h[0] = tmng[0]; h[1] = tmng[1];
    hs=atoi(h);
    char m[3]={0}; m[0] = tmng[3]; m[1] = tmng[4];
    ms=atoi(m);
    char s[3]={0}; s[0] = tmng[6]; s[1] = tmng[7];
    ss=atoi(s);
    char f[4]={0}; f[0] = tmng[9]; f[1] = tmng[10], f[2] = tmng[11];
    fs=atoi(f);
    // printf("hs=%i ms=%i ss=%i fs=%i\n", hs, ms, ss, fs); 
    float t0=hs*3600+ms*60+ss+fs/1000.;
    printf("%4.4f to %4.4f\n", t0, t0+dt); 
}


void prtaawcd2(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        // printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            // printf("l%ut", aawc->aaw[i]->aw[j]->lp1-1);
            if((aawc->aaw[i]->aw[j]->t == TMNG) & (aawc->aaw[i]->aw[j]->lp1 >10))
                // printf("%s len: %i\n", aawc->aaw[i]->aw[j]->w, aawc->aaw[i]->aw[j]->lp1);
                // printf("%s(%i)\n", aawc->aaw[i]->aw[j]->w, aawc->aaw[i]->aw[j]->lp1);
                printf("%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
    printf("\n"); 
	printf("Only numbers printed\n");
}

void prtaawcd3(aaw_c *aawc, float dt) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        // printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            // printf("l%ut", aawc->aaw[i]->aw[j]->lp1-1);
            if((aawc->aaw[i]->aw[j]->t == TMNG) & (aawc->aaw[i]->aw[j]->lp1 >10))
                convt0(aawc->aaw[i]->aw[j]->w, dt);
        }
    }
	printf("Only timings printed\n");
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

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=3) {
        printf("Error. Pls supply 2 arguments (name of text file) and number of seconds and fractional seconds (as floati, to delay timings, just make it negative).\n");
        exit(EXIT_FAILURE);
    }

    aaw_c *aawc=processinpf(argv[1]);
    float dt=atof(argv[2]);
    prtaawcd3(aawc, dt); // just the metadata
    // prtaawcd2(aawc);
    printf("typeszs: aaw_c: %zu aw_c: %zu w_c: %zu\n", sizeof(aaw_c), sizeof(aw_c), sizeof(w_c));
    printf("Numlines: %zu\n", aawc->numl); 

    free_aawc(&aawc);

    return 0;
}
