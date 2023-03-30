/* modification of csvrde, for DNA MEth many-genes-per-cg csv */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "csvrdm.h"

int catchopts(optstruct *opstru, int oargc, char **oargv)
{
    int c;
    opterr = 0;
    while ((c = getopt (oargc, oargv, "aos:e:")) != -1)
        switch (c) {
            case 'o': // out only rows with regular cellnum
                opstru->oflag = 1;
                break;
            case 's': // skip first rows, like the R csvread, -s 1 will skip one line, -s 2, two lines etc.
                opstru->snum = atoi(optarg);
                break;
            case 'e': // max number of empty cells allowed
                opstru->emp = atoi(optarg);
                break;
			case '?':
				if((optopt == 's') | (optopt == 'e')) {
					fprintf (stderr, "Option -%c requires an integer.\n", optopt);
                    exit(EXIT_FAILURE);
                }
            default:
                printf("Options wrong ... probably a wrong one submitted.\n"); 
                abort();
        }
    return 0;
}

av_c *crea_avc(int vbf)
{
    av_c *avc=malloc(sizeof(av_c));
    avc->vbf=vbf;
    avc->v=malloc(avc->vbf*sizeof(int));
    avc->vsz=0;
    return avc;
}

void condrea_avc(av_c *avc)
{
    /* somewhat trivial, but idea is that, as avc is a container, it can be re-alloced inside a function */
    CONDREALLOCAV(avc->vsz, avc->vbf, GBUF, avc->v, int);
    return;
}

void norm_avc(av_c *avc)
{
    /* somewhat trivial, but idea is that, as avc is a container, it can be re-alloced inside a function */
    avc->v=realloc(avc->v, avc->vsz*sizeof(int));
    return;
}

void free_avc(av_c *avc)
{
    free(avc->v);
    free(avc);
    return;
}

void prtavecnd(char *ca, av_c *avc) /* a debug version */
{
    int i, j;
    char *c=ca;
    printf("chars %i to %i: ", 0, avc->v[0]-1);
    for(i=0;i<avc->v[0];++i)
        putchar(c[i]);
    putchar('\n');

    for(i=1;i<avc->vsz;++i) {
        printf("chars %i to %i: ", avc->v[i-1]+1, avc->v[i]-1);
        for(j=avc->v[i-1]+1;j<avc->v[i];++j) 
            putchar(c[j]);
        putchar('\n');
    }
    return;
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
    awc->avc=NULL;
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
    if(tawc->avc!=NULL)
        free_avc(tawc->avc);
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
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
}

void prtaawcplain00(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i;
    for(i=0;i<aawc->numl;++i) {
        printf("%s\n", aawc->aaw[i]->aw[6]->w);
    }
}

void prtaawcplainc7(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("%s(l.%i): ", aawc->aaw[i]->aw[6]->w, aawc->aaw[i]->aw[6]->lp1);
        if(aawc->aaw[i]->avc !=NULL) {
            printf("c7 semicolon pos: "); 
            for(j=0;j<aawc->aaw[i]->avc->vsz;++j)
                printf("%i ", aawc->aaw[i]->avc->v[j]);
        } else
            printf("No c7 semicolons."); 
        printf("\n"); 
    }
}

void prtaawcplain2(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L%u#C%u) ", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf("sz%u:", aawc->aaw[i]->aw[j]->lp1);
            if(aawc->aaw[i]->aw[j]->lp1==1)
                printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", "EEE");
            else
                printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
}

void prtaawcsum0(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;

    // now let's count number of empty cells
    int *emc=calloc(aawc->numl, sizeof(int));
    for(i=0; i<aawc->numl;i++)
        for(j=0;j<aawc->aaw[i]->al;++j)
            if(aawc->aaw[i]->aw[j]->lp1==1)
                emc[i]++;

    printf("Numlines in CSV = %zu. Let's print out number of cells per line, and num empty cells in sequence.\n", aawc->numl);
    for(i=0; i<aawc->numl;i++)
        printf((i!=aawc->numl-1)?"L%i:C%u:E%i ":"L%i:C%u:E%i\n", i, aawc->aaw[i]->al, emc[i]);

    free(emc);
}

void prtaawcsum2(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i;

    printf("Numlines in CSV = %zu. Details:\n", aawc->numl);
    // printf("C%u ", aawc->aaw[0]->al);
    int timesm1=0;
    int mxt=1;
    int mxqn=0;
    for(i=1; i<aawc->numl+1;i++) {
        if(aawc->aaw[i]->al == aawc->aaw[i-1]->al) {
            timesm1++;
            // cell num repeated ... print nothing.
            if(mxt<timesm1) {
                mxt=timesm1;
                mxqn = aawc->aaw[i-1]->al;
            }

        } else {
            // retrospective print
            printf("C%u:T%i ", aawc->aaw[i-1]->al, timesm1+1);
            timesm1 =0;
        }
    }
    printf("\nMost common cellquan=%i with %i times\n", mxqn, mxt+1);
}

cattrs *retattr3(aaw_c *aawc, int skip, int mxemp) /* print line and word details, but not the words themselves */
{
    /*NOTE: most common cell quantity given a limited calculation here: max number of times this number is repeated "in sequence" */
    int i, j;
    cattrs *atr=calloc(1, sizeof(cattrs));
    atr->nrows = aawc->numl;

    int timesm1=0;
    int mxt=1;
    int mxqn=0;
    // following forloop, analyse current but take action on previous
    for(i=skip+1; i<aawc->numl;i++) {
        if(aawc->aaw[i]->al == aawc->aaw[i-1]->al) {
            timesm1++;
            // cell num repeated ... print nothing.
            if(mxt<timesm1) {
                mxt=timesm1;
                mxqn = aawc->aaw[i-1]->al;
            }

        } else {
            timesm1 =0;
        }
    }
    // last line
    if(aawc->aaw[aawc->numl-1]->al == aawc->aaw[aawc->numl-2]->al) {
        timesm1++;
        // cell num repeated ... print nothing.
        if(mxt<timesm1) {
            mxt=timesm1;
            mxqn = aawc->aaw[aawc->numl-1]->al;
        }
    }

    // right get summary stats together
    atr->regcellnum = mxqn;

    // actually to find the number of rows with that regcellnum needs another run through.
    int numregc=0;
    int nempc;
    int nnempc=0;
    for(i=skip; i<aawc->numl;i++) {
        nempc=0;
        for(j=0;j<aawc->aaw[i]->al;++j)
            if(aawc->aaw[i]->aw[j]->lp1 ==1)
                nempc++;
        if(nempc>=mxemp) {
            nnempc++;
            continue;
        } else if(aawc->aaw[i]->al == atr->regcellnum)
            numregc++;
    }

    atr->nregcell = numregc;
    atr->empdel = nnempc; // number of rows deleted because max num emp cells exceeded.
    return atr;
}

void prtsegunattr(aaw_c *aawc, int skip, int mxemp) /* print line and word details, but not the words themselves */
{
    /*NOTE: most common cell quantity given a limited calculation here: max number of times this number is repeated "in sequence" */
    int i, j;

    int timesm1=0;
    int mxt=1;
    int mxqn=0;
    // following forloop, analyse current but take action on previous
    for(i=skip+1; i<aawc->numl;i++) {
        if(aawc->aaw[i]->al == aawc->aaw[i-1]->al) {
            timesm1++;
            // cell num repeated ... print nothing.
            if(mxt<timesm1) {
                mxt=timesm1;
                mxqn = aawc->aaw[i-1]->al;
            }

        } else {
            timesm1 =0;
        }
    }
    // last line
    if(aawc->aaw[aawc->numl-1]->al == aawc->aaw[aawc->numl-2]->al) {
        timesm1++;
        // cell num repeated ... print nothing.
        if(mxt<timesm1) {
            mxt=timesm1;
            mxqn = aawc->aaw[aawc->numl-1]->al;
        }
    }

    // actually to find the number of rows with that regcellnum needs another run through.
    int numregc=0;
    int nempc;
    int nnempc=0;
    for(i=skip; i<aawc->numl;i++) {
        nempc=0;
        if(aawc->aaw[i]->al != mxqn)
            continue;
        for(j=0;j<aawc->aaw[i]->al;++j)
            if(aawc->aaw[i]->aw[j]->lp1 ==1)
                nempc++;
        if(nempc>=mxemp) {
            nnempc++;
            continue;
        }
        numregc++;
        for(j=0;j<aawc->aaw[i]->al;++j)
            printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
    }
    return;
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

aaw_c *processincsv0(char *fname)
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
        if(c=='"')
            continue;
        if( (c== '\n') | (c == ',') ) {
            if(oldc==',') {
                if(couw ==aawc->aaw[aawc->numl]->ab-1) /* new word opening */
                    reall_awc(aawc->aaw+aawc->numl, WABUF);
                couw++; /* verified: this has to be here */
                aawc->aaw[aawc->numl]->aw[couw]->w[0]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=1;
                norm_wc(aawc->aaw[aawc->numl]->aw+couw);
            } else if(inword==1) {
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
        } else if(inword==0) {
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
            if(oldc==',') {
                if(couw ==aawc->aaw[aawc->numl]->ab-1) /* new word opening */
                    reall_awc(aawc->aaw+aawc->numl, WABUF);
                couw++; /* verified: this has to be here */
                aawc->aaw[aawc->numl]->aw[couw]->w[0]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=1;
                if(couw==6)
                    norm_avc(aawc->aaw[aawc->numl]->avc);
                norm_wc(aawc->aaw[aawc->numl]->aw+couw);
            } else if(inword==1) {
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
        } else if(inword==0) {
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
            // locate semcolons in c7:
            if((couw==6) & (c==';')) {
                if(aawc->aaw[aawc->numl]->avc==NULL) {
                    aawc->aaw[aawc->numl]->avc = crea_avc(GBUF);
                    aawc->aaw[aawc->numl]->avc->v[aawc->aaw[aawc->numl]->avc->vsz++]=couc-1;
                } else {
                    CONDREALLOCAV(aawc->aaw[aawc->numl]->avc->vsz, aawc->aaw[aawc->numl]->avc->vbf, GBUF, aawc->aaw[aawc->numl]->avc->v, int);
                    aawc->aaw[aawc->numl]->avc->v[aawc->aaw[aawc->numl]->avc->vsz++]=couc-1;
                }
            }
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

void prtusage(void)
{
    printf("Usage notes: csvrde\n"); 
    printf("             program to read a (rough i.e. surveymonkey) CSV and regularise it\n");
    printf("             First argument must be the name of csv file to analyse\n");
    printf("             Can be followed by -o flag, which will output the regular-cell-num rows\n");
    printf("             Also can be followed by -s option with integer, this is number of first lines to skip.\n");
    return;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc==1) {
        prtusage();
        exit(EXIT_FAILURE);
    }
    int argignore=1; // first arg is not of option type .. it's the input CSV
    int oargc=argc-argignore;
    char **oargv=argv+argignore;
    optstruct opstru={0};
    catchopts(&opstru, oargc, oargv);

    printf("typeszs: aaw_c: %zu aw_c: %zu w_c: %zu\n", sizeof(aaw_c), sizeof(aw_c), sizeof(w_c));

    aaw_c *aawc=processincsv(argv[1]);
    // prtaawcdbg(aawc);
    // prtaawcplain00(aawc);
    prtaawcplainc7(aawc);
    // prtaawcsum2(aawc); // printout original text as well as you can.
    free_aawc(&aawc);

    return 0;
}
