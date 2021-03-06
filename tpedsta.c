/* some statistics from a ped file ... specificaly fro Illumina BeadArray chip target */
#include "tpedsta.h"

int catchopts(optstruct *opstru, int oargc, char **oargv)
{
    int c;
    opterr = 0;
    while ((c = getopt (oargc, oargv, "ci:f:")) != -1)
        switch (c) {
            case 'c': // want to see the genotypes of these SNPs
                opstru->cflag = 1;
                break;
            case 'i': // input file 
                opstru->iname = optarg;
                break;
            case 'f': // input file 
                opstru->fname = optarg; /* the tfam of fam file */
                break;
			case '?':
				if (optopt == 'i') {
					fprintf (stderr, "Option -%c requires an input tped name.\n", optopt);
                    exit(EXIT_FAILURE);
                }
            default:
                abort();
        }
    return 0;
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
    awc->nn=calloc(CPSTRSZ, sizeof(char));
    awc->gd=0;
    awc->gdn=0;
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
    free(tawc->nn);
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

gt_t from2l(char A1, char A2)
{
    /* return a gt type for two letters */
    gt_t tgt; // The temporary GT
    switch(A1) {
        case 'A':
            switch(A2){
                case 'A':
                    tgt=AA; break;
                case 'C':
                    tgt=AC; break;
                case 'G':
                    tgt=AG; break;
                case 'T':
                    tgt=AT; break;
                case 'I': case 'D':
                    tgt=D1; break;
                case 'N': case '0':
                    tgt=Z1; break;
                default:
                    tgt=ZZ;
            }
            break;
        case 'C':
            switch(A2){
                case 'A':
                    tgt=CA; break;
                case 'C':
                    tgt=CC; break;
                case 'G':
                    tgt=CG; break;
                case 'T':
                    tgt=CT; break;
                case 'I': case 'D':
                    tgt=D1; break;
                case 'N': case '0':
                    tgt=Z1; break;
                default:
                    tgt=ZZ;
            }
            break;
        case 'G':
            switch(A2){
                case 'A':
                    tgt=GA; break;
                case 'C':
                    tgt=GC; break;
                case 'G':
                    tgt=GG; break;
                case 'T':
                    tgt=GT; break;
                case 'I': case 'D':
                    tgt=D1; break;
                case 'N': case '0':
                    tgt=Z1; break;
                default:
                    tgt=ZZ;
            }
            break;
        case 'T':
            switch(A2){
                case 'A':
                    tgt=TA; break;
                case 'C':
                    tgt=TC; break;
                case 'G':
                    tgt=TG; break;
                case 'T':
                    tgt=TT; break;
                case 'I': case 'D':
                    tgt=D1; break;
                case 'N': case '0':
                    tgt=Z1; break;
                default:
                    tgt=ZZ;
            }
            break;
        case 'D': case 'I':
            switch(A2){
                case 'A': case 'C': case 'G': case 'T':
                    tgt=D1; break;
                default:
                    tgt=D2;
            }
            break;
        case 'N': case '0':
            switch(A2){
                case 'A': case 'C': case 'G': case 'T':
                    tgt=Z1; break;
                case 'I': case 'D':
                    /* special case where NI or ND occur: an ID was involced, so to speak */
                    tgt=D2; break;
                default:
                    tgt=ZZ;
            }
            break;
        default:
            tgt=ZZ;
    }
    return tgt;
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
    int j;
    size_t i;
    size_t a[12]; // counts which are "of interest"
    size_t suma1good, suma1all, suma2good, suma2all;
    float ra1, ra2, allra1=.0, allra2=.0;

    printf("Samplename\t#A1==A\t#A1==C\t#A1==G\t#A1==T\t#A1==0\tA1==?\t #A2==A\t#A2==C\t#A2==G\t#A2==T\t#A2==0\tA2==?\tA1CR\tA2CR\n");
    int numsamps =0; // because numl is not always numsamps due to # comment lines.
    for(i=0;i<aawc->numl;++i) {
        if(aawc->aaw[i]->aw[0]->w[0] == '#')
            continue;
        numsamps++;
        for(j=0;j<12;++j) 
            a[j]= 0;
        // OK now so now we have a match on the IIDs
        // check num genos at very least
        printf("%s\t", aawc->aaw[i]->aw[1]->w);
        for(j=6;j<aawc->aaw[i]->al;j+=2) {
            a1 = aawc->aaw[i]->aw[j]->w[0];
            a2 = aawc->aaw[i]->aw[j+1]->w[0];
            switch(a1) {
                case 'A':
                    a[0]++; break;
                case 'C':
                    a[1]++; break;
                case 'G':
                    a[2]++; break;
                case 'T':
                    a[3]++; break;
                case '0':
                    a[4]++; break;
                default:
                    a[5]++;
            }
            switch(a2) {
                case 'A':
                    a[6]++; break;
                case 'C':
                    a[7]++; break;
                case 'G':
                    a[8]++; break;
                case 'T':
                    a[9]++; break;
                case '0':
                    a[10]++; break;
                default:
                    a[11]++;
            }
        }
        for(j=0;j<12;++j) 
            printf("%zu\t", a[j]);
        suma1good=a[0]+a[1]+a[2]+a[3];
        suma1all=suma1good+a[4]+a[5];
        suma2good=a[6]+a[7]+a[8]+a[9];
        suma2all=suma2good+a[10]+a[11];

        ra1=(float)suma1good/suma1all;
        ra2=(float)suma2good/suma2all;
        allra1 += ra1;
        allra2 += ra2;
        printf("%4.4f\t%4.4f\n", ra1, ra2);
    }
    printf("Overall: AvgA1CR=%4.4f AvgA2CR=%4.4f\n", allra1/numsamps, allra2/numsamps);
}

void statsaawc2(aaw_c *aawc, float *allra1, float *allra2, int *retnumsamps)
{
    char a1, a2;
    int j;
    size_t i;
    size_t a[10]; // counts which are "of interest": 
    size_t suma1all, suma2all;
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
        suma1all=a[0]+a[2]+a[4]+a[6]+a[8];
        suma2all=a[1]+a[3]+a[5]+a[7]+a[9];

        ra1=(float)(a[0]+a[6])/suma1all;
        printf("%4.4f\t", ra1);
        ra2=(float)(a[1]+a[7])/suma2all;
        printf("%4.4f\n", ra2);
        /* now how we include IDs because hese are properly called. Often we don't used them though */
        *allra1 += ra1;
        *allra2 += ra2;
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
    int ll; /* length of line */
    boole asamell=0; /* all same line length */
    for(i=0;i<aawc->numl;++i) {
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        if(!i)
            ll = aawc->aaw[i]->al;
        else
            if( ll != aawc->aaw[i]->al) {
                asamell++;
                printf("line %i, %i words\n", i, aawc->aaw[i]->al); 
            }

        for(j=0;j<aawc->aaw[i]->al;++j)
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
    }
    if(!asamell)
        printf("All lines had same number of words (%i).\n", ll); 
    else
        printf("Not all lines had same number of wordsi line 0 %i.\n", ll); 
}

void prt_tfc(aaw_c *tfc)
{
    int i, j;
    int ll0=0; /* length of line */
    boole asamell0=0; /* all same line length? */
    for(i=0;i<tfc->numl;++i) {
#ifdef DBG
        printf("L%u(%uw):", i, tfc->aaw[i]->al); 
#endif
        if(!i)
            ll0 = tfc->aaw[i]->al;
        else
            if( ll0 != tfc->aaw[i]->al) {
                asamell0++;
#ifdef DBG
                printf("exception - line %i, %i words\n", i, tfc->aaw[i]->al); 
#endif
            }

#ifdef DBG
        for(j=0;j<tfc->aaw[i]->al;++j)
            printf((j!=tfc->aaw[i]->al-1)?"%s ":"%s\n", tfc->aaw[i]->aw[j]->w);
#endif
    }
    if(!asamell0)
        printf("All lines in .tfam had same number of words (%i).\n", ll0); 
    else
        printf("Not all lines in .tfam had same number of words as line 0 %i.\n", ll0); 

    return;
}

void vertptf(aaw_c *tfc, aaw_c *aawc) /* verify the tped and tfam file pair */
{
    int i, j;
    int ll0; /* length of line */
    boole asamell0=0; /* all same line length? */
    for(i=0;i<tfc->numl;++i) {
#ifdef DBG
        printf("L%u(%uw):", i, tfc->aaw[i]->al); 
#endif
        if(!i)
            ll0 = tfc->aaw[i]->al;
        else
            if( ll0 != tfc->aaw[i]->al) {
                asamell0++;
#ifdef DBG
                printf("exception - line %i, %i words\n", i, tfc->aaw[i]->al); 
#endif
            }

#ifdef DBG
        for(j=0;j<tfc->aaw[i]->al;++j)
            printf((j!=tfc->aaw[i]->al-1)?"%s ":"%s\n", tfc->aaw[i]->aw[j]->w);
#endif
    }
    if(!asamell0)
        printf("All lines in .tfam had same number of words (%i).\n", ll0); 
    else
        printf("Not all lines in .tfam had same number of words as line 0 %i.\n", ll0); 


    /* OK. now the .tped's turn */
    boole asamell=0;
    int ll=0;
    for(i=0;i<aawc->numl;++i) {
#ifdef DBG
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
#endif
        if(!i)
            ll = aawc->aaw[i]->al;
        else
            if( ll != aawc->aaw[i]->al) {
                asamell++;
#ifdef DBG
                printf("line %i, %i words\n", i, aawc->aaw[i]->al); 
#endif
            }

#ifdef DBG
        for(j=0;j<aawc->aaw[i]->al;++j)
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
#endif
    }
    if(!asamell)
        printf("All lines in .tped had same number of words (%i).\n", ll); 
    else
        printf("Not all lines .tped had same number of words as line 0 (%i).\n", ll); 

    /* Finally check that .tped has same number of genotypes as lines in .tfam
     * Note this is hardly a great check, just because the numbers coincide doesn't mean
     * a whole lot, but it's just a sanity check. */
    int ngts=(ll-4)/2;
    if(!asamell & !asamell0) {
        if(ngts != tfc->numl)
            printf("Number of lines in .tfam does not match number of GTs in .tped.\n"); 
        else
            printf("All clear: .tped GT quantity and number of lines in .tfam coincide.\n");
    } else
        printf("Either .tped of .tfam do not have a consistently equal number of words per line.\n"); 
    return;
}

void prt_tpedaawc0f(aaw_c *aawc) /* print GT friendly ... not a valid tped, just convenient for reading */
{
    /* Note: IDN0's printed as they are */
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        /* first four are the non GT fields */
        // printf("%s\t", aawc->aaw[i]->aw[0]->w);
        printf("%s ", aawc->aaw[i]->aw[1]->w);
        // printf("%s\t", aawc->aaw[i]->aw[2]->w);
        // printf("%s\t", aawc->aaw[i]->aw[3]->w);
        for(j=4;j<aawc->aaw[i]->al;j+=2) {
            printf("%c", aawc->aaw[i]->aw[j]->w[0]);
            printf(((j+1)!=aawc->aaw[i]->al-1)?"%c ":"%c\n", aawc->aaw[i]->aw[j+1]->w[0]);
        }
    }
    return;
}

void prt_tpedaawc1f(aaw_c *aawc) /* print GT friendly ... not a valid tped, just convenient for reading */
{
    /* Note: IDN0's converted to 0's */
    int i, j;
    char a1, a2;
    for(i=0;i<aawc->numl;++i) {
        /* first four are the non GT fields */
        // printf("%s\t", aawc->aaw[i]->aw[0]->w);
        printf("%s ", aawc->aaw[i]->aw[1]->w);
        // printf("%s\t", aawc->aaw[i]->aw[2]->w);
        // printf("%s\t", aawc->aaw[i]->aw[3]->w);
        for(j=4;j<aawc->aaw[i]->al;j+=2) {
            a1=aawc->aaw[i]->aw[j]->w[0];
            a2=aawc->aaw[i]->aw[j+1]->w[0];
            switch (a1) {
                case 'N': case '0': case 'D': case 'I':
                    a1='0'; break;
                default:
                    break;
            }
            switch (a2) {
                case 'N': case '0': case 'D': case 'I':
                    a2='0'; break;
                default:
                    break;
            }
            printf("%c", a1);
            printf(((j+1)!=aawc->aaw[i]->al-1)?"%c ":"%c\n", a2);
        }
    }
    return;
}

void prt_tpedaawc1p(aaw_c *aawc) /* print proper tped */
{
    /* Note: IDN0's converted to 0's */
    int i, j;
    char a1, a2;
    for(i=0;i<aawc->numl;++i) {
        /* first four are the non GT fields */
        printf("%s ", aawc->aaw[i]->aw[0]->w);
        printf("%s ", aawc->aaw[i]->aw[1]->w);
        printf("%s ", aawc->aaw[i]->aw[2]->w);
        printf("%s ", aawc->aaw[i]->aw[3]->w);
        for(j=4;j<aawc->aaw[i]->al;j+=2) {
            a1=aawc->aaw[i]->aw[j]->w[0];
            a2=aawc->aaw[i]->aw[j+1]->w[0];
            switch (a1) {
                case 'N': case '0': case 'D': case 'I':
                    a1='0'; break;
                default:
                    break;
            }
            switch (a2) {
                case 'N': case '0': case 'D': case 'I':
                    a2='0'; break;
                default:
                    break;
            }
            printf("%c ", a1);
            printf(((j+1)!=aawc->aaw[i]->al-1)?"%c ":"%c\n", a2);
        }
    }
    return;
}

void prt_tpedaawc0(aaw_c *aawc) /* print GT friendly ... not a valid tped, just convenient for reading */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        /* first four are the non GT fields */
        // printf("%s\t", aawc->aaw[i]->aw[0]->w);
        printf("%s ", aawc->aaw[i]->aw[1]->w);
        // printf("%s\t", aawc->aaw[i]->aw[2]->w);
        // printf("%s\t", aawc->aaw[i]->aw[3]->w);
        for(j=4;j<aawc->aaw[i]->al;j+=2) {
            printf("%c", aawc->aaw[i]->aw[j]->w[0]);
            printf(((j+1)!=aawc->aaw[i]->al-1)?"%c ":"%c\n", aawc->aaw[i]->aw[j+1]->w[0]);
        }
    }
    return;
}

void cougt_tpedaawc(aaw_c *aawc, int **cougt, boole *eqngts) /* just count the gts that come out: any SNP with N or 0, means the GT is NN */
{
    int *cougt_=*cougt;
    int i, j;
    char a1, a2;
    boole eqngts_ = 0;
    int ingts, ngts = (aawc->aaw[0]->al-4)/2;
    printf("First sample numgts=%i\n", ngts); 
    for(i=0;i<aawc->numl;++i) {
        ingts =(aawc->aaw[i]->al-4)/2;
        if(ingts != ngts) {
            eqngts_=1;
            printf("Idx %i had %i GTs!\n", i, ingts); 
        }
        for(j=4;j<aawc->aaw[i]->al;j+=2) {
            a1=aawc->aaw[i]->aw[j]->w[0];
            a2=aawc->aaw[i]->aw[j+1]->w[0];
            cougt_[from2l(a1,a2)]++;
        }
    }
    // *cougt=cougt_; // not nec.
    *eqngts = eqngts_;
    return;
}

aaw_c *processinpf(FILE *fp)
{
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */;
    int c, oldc='\0';
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
                // goto uit; /// this forces dropping out when newline is reached, and that's all, right?
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
    if((argc==1) | (argc>5) ){
        printf("usage: Pls supply one tpedfile as argument with -i option.\n");
        printf("       A .tfam file can be included with the -f option.\n");
        printf("       if you include -c flag, a converted tped (missing=0) sent to STDOUT.\n");
        exit(EXIT_FAILURE);
    }
    int argignore=0; //
    int oargc=argc-argignore;
    char **oargv=argv+argignore;
    optstruct opstru={0};
    catchopts(&opstru, oargc, oargv);

    int i;
    FILE *fp=NULL, *ff=NULL;
    fp=fopen(opstru.iname, "r");
    aaw_c *aawc=NULL, *tfc=NULL;
    aawc=processinpf(fp);
    if(opstru.fname) {
        ff=fopen(opstru.fname, "r");
        tfc=processinpf(ff);
    }

    /* Now this is a bare line-word type struct, which could conceivably be converted into a more tped-friendly structure
     * but, as usual, this woul require copying and holding two tped's in memory, so why bother? Let's just be careful and keep
     * the particularities of the tped structure in mind: */
    // prt_tpedaawc1f(aawc); /* friendly not proper tped print: CHRIST I had kep this in */
    /* First, let's start with counting the 18 categorized genotypes. */
    int *cougt=calloc(NUMGTS, sizeof(int));
    boole eqngts = 0; /* equal number of GTs in our tped file? */

    /* OK, we're sticking in some options */
    if( (opstru.iname != NULL) & !(opstru.cflag) & !(opstru.fname)) {
        cougt_tpedaawc(aawc, &cougt, &eqngts);
        printf("Different GT counts: (Z1: just one uncalled allele, ZZ: both alleles uncalled.\n"); 
        for(i=0;i<NUMGTS;++i) 
            printf((i==NUMGTS-1)?"%5s \n":"%5s  ", gtna0[i]);
        for(i=0;i<NUMGTS;++i) 
            printf((i==NUMGTS-1)?"%6i\n":"%6i ", cougt[i]);
        printf((eqngts)?"Problem: unequal number of GTs across samples.\n":"Checked: yes, an equal num of GTs for all samples.\n");
    } else if( !(opstru.cflag) & !(opstru.iname) & (opstru.fname != NULL) ) {
        prt_tfc(tfc);
    } else if( !(opstru.cflag) & (opstru.iname != NULL) & (opstru.fname != NULL)) {
        vertptf(tfc, aawc);/* verify the tped and tfam file pair */
    } else if(opstru.cflag) {
        /* printing out converted tped: 0NID's all set to 00 */
        prt_tpedaawc1p(aawc);
    }

    // prtaawcplain(aawc);
    // statsaawc2(aawc, &allra1, &allra2, &numsamps);
    free_aawc(&aawc);
    // printf("overall: avga1cr=%4.4f avga2cr=%4.4f\n", allra1/numsamps, allra2/numsamps);
    fclose(fp);
    free(cougt);
    if(opstru.fname != NULL) {
        fclose(ff);
        free_aawc(&tfc);
    }

    return 0;
}
