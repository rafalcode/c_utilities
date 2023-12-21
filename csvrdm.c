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
    awc->av2=NULL;
    awc->av3=NULL;
    awc->av4=NULL;
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
    if(tawc->av2!=NULL)
        free_avc(tawc->av2);
    if(tawc->av3!=NULL)
        free_avc(tawc->av3);
    if(tawc->av4!=NULL)
        free_avc(tawc->av4);
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
    int sta, end;
    for(i=0;i<aawc->numl;++i) {
        printf("%s(l.%i): ", aawc->aaw[i]->aw[6]->w, aawc->aaw[i]->aw[6]->lp1);
        sta=0;
        if(aawc->aaw[i]->avc !=NULL) {
            for(j=0;j<aawc->aaw[i]->avc->vsz;++j) {
                end=aawc->aaw[i]->avc->v[j]-1;
                printf("%i->%i ", sta, end);
                sta=aawc->aaw[i]->avc->v[j]+1;
            }
            end=aawc->aaw[i]->aw[6]->lp1-1;
            printf("%i->%i ", sta, end);
        } else
            printf("No c7 semicolons."); 
        printf("\n"); 
    }
}

void prtaawcplainc78(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    int sta, end;
    for(i=0;i<aawc->numl;++i) {
        printf("%s(l.%i): ", aawc->aaw[i]->aw[6]->w, aawc->aaw[i]->aw[6]->lp1);
        printf("%s(l.%i): ", aawc->aaw[i]->aw[7]->w, aawc->aaw[i]->aw[7]->lp1);
        sta=0;
        if(aawc->aaw[i]->avc !=NULL) {
            for(j=0;j<aawc->aaw[i]->avc->vsz;++j) {
                end=aawc->aaw[i]->avc->v[j]-1;
                printf("%i->%i ", sta, end);
                sta=aawc->aaw[i]->avc->v[j]+1;
            }
            end=aawc->aaw[i]->aw[6]->lp1-1;
            printf("%i->%i ", sta, end);
        } else
            printf("No c7 semicolons."); 
        sta=0;
        if(aawc->aaw[i]->av2 !=NULL) {
            for(j=0;j<aawc->aaw[i]->av2->vsz;++j) {
                end=aawc->aaw[i]->av2->v[j]-1;
                printf("%i->%i ", sta, end);
                sta=aawc->aaw[i]->av2->v[j]+1;
            }
            end=aawc->aaw[i]->aw[7]->lp1-1;
            printf("%i->%i ", sta, end);
        } else
            printf("No c8 semicolons."); 
        printf("\n"); 
    }
}

void prtaawcplainc78_(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i;
    for(i=0;i<aawc->numl;++i) {
        if((aawc->aaw[i]->avc !=NULL) & (aawc->aaw[i]->av2 !=NULL))
            if(aawc->aaw[i]->avc->vsz != aawc->aaw[i]->av2->vsz)
                printf("PROBLEM@l%i: c7sc:%i != c8sc:%i\n", i, aawc->aaw[i]->avc->vsz, aawc->aaw[i]->av2->vsz);
    }
}

void prtaawcplainc910_(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i;
    for(i=0;i<aawc->numl;++i) {
        if((aawc->aaw[i]->avc !=NULL) & (aawc->aaw[i]->av2 !=NULL))
            if(aawc->aaw[i]->avc->vsz != aawc->aaw[i]->av2->vsz)
                printf("PROBLEM@l%i: c7sc:%i != c8sc:%i\n", i, aawc->aaw[i]->avc->vsz, aawc->aaw[i]->av2->vsz);
        if((aawc->aaw[i]->av3 !=NULL) & (aawc->aaw[i]->av4 !=NULL))
            if(aawc->aaw[i]->av3->vsz != aawc->aaw[i]->av4->vsz)
                printf("PROBLEM@l%i: c9sc:%i != c10sc:%i\n", i, aawc->aaw[i]->av3->vsz, aawc->aaw[i]->av4->vsz);
        if((aawc->aaw[i]->avc !=NULL) & (aawc->aaw[i]->av3 !=NULL))
            printf("UCSCRGsc=%i vs. GCODEV12sc=%i\n", aawc->aaw[i]->avc->vsz, aawc->aaw[i]->av3->vsz);
    }
}

void prtaawcplainav0(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j, k, kk;
    int sta, end;
    for(i=0;i<aawc->numl;++i) {
        printf("C7length=%i::",aawc->aaw[i]->aw[6]->lp1-1);
        if(aawc->aaw[i]->avc !=NULL) {
            sta=0;
            for(k=0;k<aawc->aaw[i]->avc->vsz;++k) {
                end=aawc->aaw[i]->avc->v[k];
                for(kk=sta;kk<end;++kk)
                    putchar(aawc->aaw[i]->aw[6]->w[kk]);
                putchar(','); 
                for(j=0;j<aawc->aaw[i]->al;++j) {
                    if(j==6) continue;
                    printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
                }
                sta=aawc->aaw[i]->avc->v[k]+1;
            }
            /* last semicolon */
            end=aawc->aaw[i]->aw[6]->lp1-1;
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[6]->w[kk]);
            putchar(','); 
            for(j=0;j<aawc->aaw[i]->al;++j)
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
        } else if(aawc->aaw[i]->aw[6]->lp1>1) {
            printf("%s,", aawc->aaw[i]->aw[6]->w);
            for(j=0;j<aawc->aaw[i]->al;++j)
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
        } else {
            putchar(','); // no semicolon
            for(j=0;j<aawc->aaw[i]->al;++j)
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
}

void prtaawcplainav(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j, k, kk;
    int sta, end;
    for(i=0;i<aawc->numl;++i) {
        //col7:
        if(aawc->aaw[i]->avc !=NULL) {
            sta=0;
            for(k=0;k<aawc->aaw[i]->avc->vsz;++k) {
                for(j=0;j<6;++j)
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                end=aawc->aaw[i]->avc->v[k];
                for(kk=sta;kk<end;++kk)
                    putchar(aawc->aaw[i]->aw[6]->w[kk]);
                putchar(','); 
                for(j=7;j<aawc->aaw[i]->al;++j)
                    printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
                sta=aawc->aaw[i]->avc->v[k]+1;
            }
            /* last semicolon */
            end=aawc->aaw[i]->aw[6]->lp1-1;
            for(j=0;j<6;++j)
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[6]->w[kk]);
            putchar(','); 
            for(j=7;j<aawc->aaw[i]->al;++j)
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
        } else
            for(j=0;j<aawc->aaw[i]->al;++j)
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
    }
}

void prtaawcplainav2(aaw_c *aawc) /* prints col7 and col8 parsed. the semicolons must be the exact same: they mostly(!) are.*/
{
    int i, j, k, kk;
    int sta, end;
    int sta2, end2;
    for(i=0;i<aawc->numl;++i) {
        //col7:
        if(aawc->aaw[i]->avc !=NULL) {
            sta=0;
            sta2=0;
            for(k=0;k<aawc->aaw[i]->avc->vsz;++k) {
                for(j=0;j<6;++j)
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                end=aawc->aaw[i]->avc->v[k];
                for(kk=sta;kk<end;++kk)
                    putchar(aawc->aaw[i]->aw[6]->w[kk]);
                putchar(','); 
                sta=aawc->aaw[i]->avc->v[k]+1;
                end2=aawc->aaw[i]->av2->v[k];
                for(kk=sta2;kk<end2;++kk)
                    putchar(aawc->aaw[i]->aw[7]->w[kk]);
                putchar(','); 
                for(j=8;j<aawc->aaw[i]->al;++j)
                    printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
                sta2=aawc->aaw[i]->av2->v[k]+1;
            }
            /* last semicolon */
            end=aawc->aaw[i]->aw[6]->lp1-1;
            for(j=0;j<6;++j)
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[6]->w[kk]);
            putchar(','); 
            end2=aawc->aaw[i]->aw[7]->lp1-1;
            for(kk=sta2;kk<end2;++kk)
                putchar(aawc->aaw[i]->aw[7]->w[kk]);
            putchar(','); 
            for(j=8;j<aawc->aaw[i]->al;++j)
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
        } else
            for(j=0;j<aawc->aaw[i]->al;++j)
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
    }
}

void prtaawcplainav30(aaw_c *aawc) /* prints col7 and col8 parsed. the semicolons must be the exact same: they mostly(!) are.*/
{
    /* no gencode printing */
    int i, j, k, kk;
    int sta, end;
    int sta2, end2;

    //header: first empty cell gets skipped .. sort of a bug.
    printf("Genename,Genegrp,Cpgname,"); 
    for(j=0;j<3;++j) {
        if(j==5) continue;
        printf("%s,", aawc->aaw[0]->aw[j]->w);
    }
    printf("%s,", aawc->aaw[0]->aw[4]->w);
    for(j=9;j<aawc->aaw[0]->al;++j) {
        printf((j!=aawc->aaw[0]->al-1)?"%s,":"%s\n", aawc->aaw[0]->aw[j]->w);
    }
    for(i=1;i<aawc->numl;++i) {
        //col7:
        if(aawc->aaw[i]->avc !=NULL) {
            sta=0;
            sta2=0;
            for(k=0;k<aawc->aaw[i]->avc->vsz;++k) {
                end=aawc->aaw[i]->avc->v[k];
                for(kk=sta;kk<end;++kk)
                    putchar(aawc->aaw[i]->aw[6]->w[kk]);
                putchar(','); 
                sta=aawc->aaw[i]->avc->v[k]+1;
                end2=aawc->aaw[i]->av2->v[k];
                for(kk=sta2;kk<end2;++kk)
                    putchar(aawc->aaw[i]->aw[7]->w[kk]);
                putchar(','); 
                for(j=0;j<4;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                printf("%s,", aawc->aaw[i]->aw[5]->w);
                for(j=10;j<aawc->aaw[i]->al;++j) {
                    printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
                }
                sta2=aawc->aaw[i]->av2->v[k]+1;
            }
            /* last semicolon */
            end=aawc->aaw[i]->aw[6]->lp1-1;
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[6]->w[kk]);
            putchar(','); 
            end2=aawc->aaw[i]->aw[7]->lp1-1;
            for(kk=sta2;kk<end2;++kk)
                putchar(aawc->aaw[i]->aw[7]->w[kk]);
            putchar(','); 
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
            }
        } else if(aawc->aaw[i]->aw[6]->lp1>1) /* single USCS refgenes */ {
            printf("%s,", aawc->aaw[i]->aw[6]->w);
            printf("%s,", aawc->aaw[i]->aw[7]->w);
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
            }
        } else if(aawc->aaw[i]->aw[8]->lp1>1) /* single GENCODEV12 genes */ {
            printf("%s,", aawc->aaw[i]->aw[8]->w);
            printf("%s,", aawc->aaw[i]->aw[9]->w);
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf((j!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[j]->w);
            }
        }
    }
}

void prtaawcplainav3(aaw_c *aawc) /* prints col7 and col8 parsed. the semicolons must be the exact same: they mostly(!) are.*/
{
    /* no gencode printing */
    int i, j, k, kk;
    int sta, end;
    int sta2, end2;

    //header: first empty cell gets skipped .. sort of a bug.
    printf("Genename,CpgGrp,Cpgname,"); 
    for(j=0;j<3;++j) {
        if(j==5) continue;
        printf("%s,", aawc->aaw[0]->aw[j]->w);
    }
    printf("%s,", aawc->aaw[0]->aw[4]->w);
    for(j=9;j<aawc->aaw[0]->al;++j) {
        printf("%s,", aawc->aaw[0]->aw[j]->w);
    }
    printf("UCSCRG_GCV12\n");  // the new RF probability of biological impact score.
    int couavnn=0; //count avc not nulls: this first if
    int cousrg=0;  // count UCSC single transcript/splicevar genes
    int couav34nn=0;
    int couav34sg=0;
    int coualloth=0; // count all others.
    for(i=1;i<aawc->numl;++i) {
        if(aawc->aaw[i]->avc !=NULL) {
            couavnn++;
            sta=0;
            sta2=0;
            for(k=0;k<aawc->aaw[i]->avc->vsz;++k) {
                end=aawc->aaw[i]->avc->v[k];
                for(kk=sta;kk<end;++kk)
                    putchar(aawc->aaw[i]->aw[6]->w[kk]);
                putchar(','); 
                sta=aawc->aaw[i]->avc->v[k]+1;
                end2=aawc->aaw[i]->av2->v[k];
                for(kk=sta2;kk<end2;++kk)
                    putchar(aawc->aaw[i]->aw[7]->w[kk]);
                putchar(','); 
                for(j=0;j<4;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                printf("%s,", aawc->aaw[i]->aw[5]->w);
                for(j=10;j<aawc->aaw[i]->al;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                printf("UCSCRG\n"); // Genname coded by UCSC Refgene
                sta2=aawc->aaw[i]->av2->v[k]+1;

            }
            /* last semicolon */
            end=aawc->aaw[i]->aw[6]->lp1-1;
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[6]->w[kk]);
            putchar(','); 
            end2=aawc->aaw[i]->aw[7]->lp1-1;
            for(kk=sta2;kk<end2;++kk)
                putchar(aawc->aaw[i]->aw[7]->w[kk]);
            putchar(','); 
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("UCSCRG\n"); // Genname coded by UCSC Refgene
        } else if(aawc->aaw[i]->aw[6]->lp1>1) /* single USCS refgenes */ {
            cousrg++;
            printf("%s,", aawc->aaw[i]->aw[6]->w);
            printf("%s,", aawc->aaw[i]->aw[7]->w);
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("UCSCRG\n"); // Genname coded by UCSC Refgene
        /* OK, now we[re in gencode territory . they are less consistent than refgen */
        } else if((aawc->aaw[i]->av3 !=NULL) & (aawc->aaw[i]->av4 !=NULL)) {
            // Move to GCv12:
            couav34nn++;
            sta=0;
            sta2=0;
            for(k=0;k<aawc->aaw[i]->av3->vsz;++k) {
                end=aawc->aaw[i]->av3->v[k];
                for(kk=sta;kk<end;++kk)
                    putchar(aawc->aaw[i]->aw[8]->w[kk]);
                putchar(','); 
                sta=aawc->aaw[i]->av3->v[k]+1;
                end2=aawc->aaw[i]->av4->v[k];
                for(kk=sta2;kk<end2;++kk)
                    putchar(aawc->aaw[i]->aw[9]->w[kk]);
                putchar(','); 
                for(j=0;j<4;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                printf("%s,", aawc->aaw[i]->aw[5]->w);
                for(j=10;j<aawc->aaw[i]->al;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                printf("GCV12\n"); // Genname coded by UCSC Refgene
                sta2=aawc->aaw[i]->av4->v[k]+1;
            }
            /* last semicolon */
            end=aawc->aaw[i]->aw[8]->lp1-1;
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[8]->w[kk]);
            putchar(','); 
            // with GCV12 there can be more CpgGrps than Genenames! Sacré GCv12!
            if(aawc->aaw[i]->av3->vsz < aawc->aaw[i]->av4->vsz)
                end2=aawc->aaw[i]->av4->v[k]; // k should be robust to this.
            else
                end2=aawc->aaw[i]->aw[9]->lp1-1;
            for(kk=sta2;kk<end2;++kk)
                putchar(aawc->aaw[i]->aw[9]->w[kk]);
            putchar(','); 
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("GCV12\n"); // Genname coded by UCSC Refgene
        } else if((aawc->aaw[i]->aw[8]->lp1>1) & (aawc->aaw[i]->av3 ==NULL) & (aawc->aaw[i]->av4 !=NULL)) {
            /* GENCODEV12 single gene but several GCV12 groups ... choose the first one */ 
            couav34sg++;
            printf("%s,", aawc->aaw[i]->aw[8]->w);
            /* now get first token of GCV12 group */
            sta=0;
            end=aawc->aaw[i]->av4->v[0];
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[9]->w[kk]);
            putchar(','); 
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("GCV12\n"); // Genname coded by UCSC Refgene
        } else {
            coualloth++;
#ifdef CHECKLOST
            for(jj=0;jj<aawc->aaw[i]->al;++jj) {
                printf((jj!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[jj]->w);
            }
#endif
        }
    }
    fprintf(stderr,"Summcounts: couavnn=%i, cousrg=%i, couav34nn=%i, couav34sg=%i, coualloth=%i, numl=%zu\n", couavnn, cousrg, couav34nn, couav34sg, coualloth, aawc->numl);
}

void prtaawcplainav4(aaw_c *aawc) /* prints col7 and col8 parsed. the semicolons must be the exact same: they mostly(!) are.*/
{
    /* no gencode printing */
    int i, j, jj, k, kk;
    int sta, end;
    int sta2, end2;

    //header: first empty cell gets skipped .. sort of a bug.
    printf("Genename,CpgGrp,Cpgname,"); 
    for(j=0;j<3;++j) {
        if(j==5) continue;
        printf("%s,", aawc->aaw[0]->aw[j]->w);
    }
    printf("%s,", aawc->aaw[0]->aw[4]->w);
    for(j=9;j<aawc->aaw[0]->al;++j) {
        printf("%s,", aawc->aaw[0]->aw[j]->w);
    }
    printf("UCSCRG_GCV12,RFBISC\n");  // the new RF probability of biological impact score.
    int couavnn=0; //count avc not nulls: this first if
    int cousrg=0;  // count UCSC single transcript/splicevar genes
    int rfbisc=0;
    int couav34nn=0;
    int couav34sg=0;
    int coualloth=0; // count all others.
    for(i=1;i<aawc->numl;++i) {
        if(aawc->aaw[i]->avc !=NULL) {
            couavnn++;
            sta=0;
            sta2=0;
            for(k=0;k<aawc->aaw[i]->avc->vsz;++k) {
                rfbisc=0;
                end=aawc->aaw[i]->avc->v[k];
                for(kk=sta;kk<end;++kk)
                    putchar(aawc->aaw[i]->aw[6]->w[kk]);
                putchar(','); 
                sta=aawc->aaw[i]->avc->v[k]+1;
                end2=aawc->aaw[i]->av2->v[k];
                for(kk=sta2;kk<end2;++kk) {
                    putchar(aawc->aaw[i]->aw[7]->w[kk]);
                    if(((kk-sta2)==3) & (aawc->aaw[i]->aw[7]->w[kk]=='2'))
                        rfbisc +=5; // a TSS200
                    else if(((kk-sta2)==3) & (aawc->aaw[i]->aw[7]->w[kk]=='1'))
                        rfbisc +=4; // a TSS1500
                    else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[7]->w[kk]=='5'))
                        rfbisc +=3; // a UTR'5
                    else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[7]->w[kk]=='1'))
                        rfbisc +=2; // a 1st Exon
                //     else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[7]->w[kk]=='3'))
                //         rfbisc +=1; // a 3'UTR
                }
                putchar(','); 
                for(j=0;j<4;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                // next is relation to island
                printf("%s,", aawc->aaw[i]->aw[5]->w);
                if(aawc->aaw[i]->aw[5]->w[0]=='I')
                    rfbisc +=3; // Island
                else if(aawc->aaw[i]->aw[5]->w[4]=='o')
                    rfbisc +=2; // north or south shore
                else if(aawc->aaw[i]->aw[5]->w[4]=='e')
                    rfbisc +=1; // north or south shore
                for(j=10;j<aawc->aaw[i]->al;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                printf("UCSCRG,%i\n", rfbisc); // Genname coded by UCSC Refgene
                sta2=aawc->aaw[i]->av2->v[k]+1;

            }
            /* last semicolon */
            rfbisc=0;
            end=aawc->aaw[i]->aw[6]->lp1-1;
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[6]->w[kk]);
            putchar(','); 
            end2=aawc->aaw[i]->aw[7]->lp1-1;
            for(kk=sta2;kk<end2;++kk) {
                putchar(aawc->aaw[i]->aw[7]->w[kk]);
                if(((kk-sta2)==3) & (aawc->aaw[i]->aw[7]->w[kk]=='2'))
                    rfbisc +=5; // a TSS200
                else if(((kk-sta2)==3) & (aawc->aaw[i]->aw[7]->w[kk]=='1'))
                    rfbisc +=4; // a TSS1500
                else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[7]->w[kk]=='5'))
                    rfbisc +=3; // a UTR'5
                else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[7]->w[kk]=='1'))
                    rfbisc +=2; // a 1st Exon
                // else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[7]->w[kk]=='3'))
                //     rfbisc +=1; // a 3'UTR
            }
            putchar(','); 
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            if(aawc->aaw[i]->aw[5]->w[0]=='I')
                rfbisc +=3; // Island
            else if(aawc->aaw[i]->aw[5]->w[4]=='o')
                rfbisc +=2; // north or south shore
            else if(aawc->aaw[i]->aw[5]->w[4]=='e')
                rfbisc +=1; // north or south shore
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("UCSCRG,%i\n", rfbisc); // Genname coded by UCSC Refgene
        } else if(aawc->aaw[i]->aw[6]->lp1>1) /* single USCS refgenes */ {
            cousrg++;
            printf("%s,", aawc->aaw[i]->aw[6]->w);
            printf("%s,", aawc->aaw[i]->aw[7]->w);
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            rfbisc=0;
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            if(aawc->aaw[i]->aw[5]->w[0]=='I')
                rfbisc +=3; // Island
            else if(aawc->aaw[i]->aw[5]->w[4]=='o')
                rfbisc +=2; // north or south shore
            else if(aawc->aaw[i]->aw[5]->w[4]=='e')
                rfbisc +=1; // north or south shore
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            if(aawc->aaw[i]->aw[7]->w[3]=='2')
                rfbisc +=5; // a TSS200
            else if(aawc->aaw[i]->aw[7]->w[3]=='1')
                rfbisc +=4; // a TSS1500
            else if(aawc->aaw[i]->aw[7]->w[0]=='5')
                rfbisc +=3; // a UTR'5
            else if(aawc->aaw[i]->aw[7]->w[0]=='1')
                rfbisc +=2; // a 1st Exon
            printf("UCSCRG,%i\n", rfbisc); // Genname coded by UCSC Refgene
        /* OK, now we[re in gencode territory . they are less consistent than refgen */
        } else if((aawc->aaw[i]->av3 !=NULL) & (aawc->aaw[i]->av4 !=NULL)) {
            // Move to GCv12:
            couav34nn++;
            sta=0;
            sta2=0;
            for(k=0;k<aawc->aaw[i]->av3->vsz;++k) {
                end=aawc->aaw[i]->av3->v[k];
                for(kk=sta;kk<end;++kk)
                    putchar(aawc->aaw[i]->aw[8]->w[kk]);
                putchar(','); 
                sta=aawc->aaw[i]->av3->v[k]+1;
                end2=aawc->aaw[i]->av4->v[k];
                rfbisc=0;
                for(kk=sta2;kk<end2;++kk) {
                    putchar(aawc->aaw[i]->aw[9]->w[kk]);
                    if(((kk-sta2)==3) & (aawc->aaw[i]->aw[9]->w[kk]=='2'))
                        rfbisc +=5; // a TSS200
                    else if(((kk-sta2)==3) & (aawc->aaw[i]->aw[9]->w[kk]=='1'))
                        rfbisc +=4; // a TSS1500
                    else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[9]->w[kk]=='5'))
                        rfbisc +=3; // a UTR'5
                    else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[9]->w[kk]=='1'))
                        rfbisc +=2; // a 1st Exon
                }
                putchar(','); 
                for(j=0;j<4;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                printf("%s,", aawc->aaw[i]->aw[5]->w);
                for(j=10;j<aawc->aaw[i]->al;++j) {
                    printf("%s,", aawc->aaw[i]->aw[j]->w);
                }
                printf("GCV12,%i\n", rfbisc); // Genname coded by UCSC Refgene
                sta2=aawc->aaw[i]->av4->v[k]+1;
            }
            /* last semicolon */
            end=aawc->aaw[i]->aw[8]->lp1-1;
            for(kk=sta;kk<end;++kk)
                putchar(aawc->aaw[i]->aw[8]->w[kk]);
            putchar(','); 
            // with GCV12 there can be more CpgGrps than Genenames! Sacré GCv12!
            if(aawc->aaw[i]->av3->vsz < aawc->aaw[i]->av4->vsz)
                end2=aawc->aaw[i]->av4->v[k]; // k should be robust to this.
            else
                end2=aawc->aaw[i]->aw[9]->lp1-1;
            rfbisc=0;
            for(kk=sta2;kk<end2;++kk) {
                putchar(aawc->aaw[i]->aw[9]->w[kk]);
                if(((kk-sta2)==3) & (aawc->aaw[i]->aw[9]->w[kk]=='2'))
                    rfbisc +=5; // a TSS200
                else if(((kk-sta2)==3) & (aawc->aaw[i]->aw[9]->w[kk]=='1'))
                    rfbisc +=4; // a TSS1500
                else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[9]->w[kk]=='5'))
                    rfbisc +=3; // a UTR'5
                else if(((kk-sta2)==0) & (aawc->aaw[i]->aw[9]->w[kk]=='1'))
                    rfbisc +=2; // a 1st Exon
            }
            putchar(','); 
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("GCV12,%i\n",rfbisc); // Genname coded by UCSC Refgene
        } else if((aawc->aaw[i]->aw[8]->lp1>1) & (aawc->aaw[i]->av3 ==NULL) & (aawc->aaw[i]->av4 !=NULL)) {
            /* GENCODEV12 single gene but several GCV12 groups ... choose the first one */ 
            couav34sg++;
            printf("%s,", aawc->aaw[i]->aw[8]->w);
            /* now get first token of GCV12 group */
            sta=0;
            end=aawc->aaw[i]->av4->v[0];
            rfbisc=0;
            for(kk=sta;kk<end;++kk) {
                putchar(aawc->aaw[i]->aw[9]->w[kk]);
                if(((kk-sta)==3) & (aawc->aaw[i]->aw[9]->w[kk]=='2'))
                    rfbisc +=5; // a TSS200
                else if(((kk-sta)==3) & (aawc->aaw[i]->aw[9]->w[kk]=='1'))
                    rfbisc +=4; // a TSS1500
                else if(((kk-sta)==0) & (aawc->aaw[i]->aw[9]->w[kk]=='5'))
                    rfbisc +=3; // a UTR'5
                else if(((kk-sta)==0) & (aawc->aaw[i]->aw[9]->w[kk]=='1'))
                    rfbisc +=2; // a 1st Exon
            }
            putchar(','); 
            for(j=0;j<4;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("%s,", aawc->aaw[i]->aw[5]->w);
            for(j=10;j<aawc->aaw[i]->al;++j) {
                printf("%s,", aawc->aaw[i]->aw[j]->w);
            }
            printf("GCV12,%i\n", rfbisc); // Genname coded by UCSC Refgene
        } else {
            coualloth++;
            // Nah, we don't want these gene-less Cpg's in our output:
            // printf("LOST: "); 
            // for(jj=0;jj<aawc->aaw[i]->al;++jj) {
            //     printf((jj!=aawc->aaw[i]->al-1)?"%s,":"%s\n", aawc->aaw[i]->aw[jj]->w);
            // }
        }
    }
    fprintf(stderr,"Summcounts: couavnn=%i, cousrg=%i, couav34nn=%i, couav34sg=%i, LOST:=%i, numl=%zu\n", couavnn, cousrg, couav34nn, couav34sg, coualloth, aawc->numl);
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
                aawc->aaw[aawc->numl]->aw[couw]->w[0]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=1;
                norm_wc(aawc->aaw[aawc->numl]->aw+couw);
                couw++; /* verified: this has to be here */
            } else if(inword==1) {
                aawc->aaw[aawc->numl]->aw[couw]->w[couc++]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=couc;
                if((couw==6) & (aawc->aaw[aawc->numl]->avc!=NULL))
                    norm_avc(aawc->aaw[aawc->numl]->avc);
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

aaw_c *processincsv2(char *fname) // this one handles the 2 avc's.
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
                aawc->aaw[aawc->numl]->aw[couw]->w[0]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=1;
                norm_wc(aawc->aaw[aawc->numl]->aw+couw);
                couw++; /* verified: this has to be here */
            } else if(inword==1) /* empty cell */ {
                aawc->aaw[aawc->numl]->aw[couw]->w[couc++]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=couc;
                if((couw==6) & (aawc->aaw[aawc->numl]->avc!=NULL))
                    norm_avc(aawc->aaw[aawc->numl]->avc);
                else if((couw==7) & (aawc->aaw[aawc->numl]->av2!=NULL))
                    norm_avc(aawc->aaw[aawc->numl]->av2);
                else if((couw==8) & (aawc->aaw[aawc->numl]->av3!=NULL))
                    norm_avc(aawc->aaw[aawc->numl]->av3);
                else if((couw==9) & (aawc->aaw[aawc->numl]->av4!=NULL))
                    norm_avc(aawc->aaw[aawc->numl]->av4);
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
            } else if((couw==7) & (c==';')) /* col8 */ {
                if(aawc->aaw[aawc->numl]->av2==NULL) {
                    aawc->aaw[aawc->numl]->av2 = crea_avc(GBUF);
                    aawc->aaw[aawc->numl]->av2->v[aawc->aaw[aawc->numl]->av2->vsz++]=couc-1;
                } else {
                    CONDREALLOCAV(aawc->aaw[aawc->numl]->av2->vsz, aawc->aaw[aawc->numl]->av2->vbf, GBUF, aawc->aaw[aawc->numl]->av2->v, int);
                    aawc->aaw[aawc->numl]->av2->v[aawc->aaw[aawc->numl]->av2->vsz++]=couc-1;
                }
            } else if((couw==8) & (c==';')) {
                if(aawc->aaw[aawc->numl]->av3==NULL) {
                    aawc->aaw[aawc->numl]->av3 = crea_avc(GBUF);
                    aawc->aaw[aawc->numl]->av3->v[aawc->aaw[aawc->numl]->av3->vsz++]=couc-1;
                } else {
                    CONDREALLOCAV(aawc->aaw[aawc->numl]->av3->vsz, aawc->aaw[aawc->numl]->av3->vbf, GBUF, aawc->aaw[aawc->numl]->av3->v, int);
                    aawc->aaw[aawc->numl]->av3->v[aawc->aaw[aawc->numl]->av3->vsz++]=couc-1;
                }
            } else if((couw==9) & (c==';')) /* col10 */ {
                if(aawc->aaw[aawc->numl]->av4==NULL) {
                    aawc->aaw[aawc->numl]->av4 = crea_avc(GBUF);
                    aawc->aaw[aawc->numl]->av4->v[aawc->aaw[aawc->numl]->av4->vsz++]=couc-1;
                } else {
                    CONDREALLOCAV(aawc->aaw[aawc->numl]->av4->vsz, aawc->aaw[aawc->numl]->av4->vbf, GBUF, aawc->aaw[aawc->numl]->av4->v, int);
                    aawc->aaw[aawc->numl]->av4->v[aawc->aaw[aawc->numl]->av4->vsz++]=couc-1;
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
    printf("Usage notes: csvrdm\n"); 
    printf("             program to read a sigDMPs.csv file (output of DNA Meth pipeline, and creates new gene rows, for different CpG groups (CpgGrp) and Island types (Relation_To_Island)\n");
    printf("             First argument must be the name of csv file to analyse\n");
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

    aaw_c *aawc=processincsv2(argv[1]);
    // prtaawcdbg(aawc);
    // prtaawcplain00(aawc);
    // prtaawcplainc7(aawc);
    // prtaawcplainc910_(aawc);
    prtaawcplainav4(aawc);
    free_aawc(&aawc);

    return 0;
}
