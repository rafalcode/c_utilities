/* some statistics from a ped file ... specificaly fro Illumina BeadArray chip target */
#include "bglsta.h"

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

void prtchaharr(snodm **stab, unsigned tsz) // print out mp_t on baiss of second (SNP name) hashing
{
    unsigned i;
    snodm *tsnod2;
    printf("SNPNAME\tPOS\tA1\tA2\n");
    for(i=0;i<tsz;++i) {
        tsnod2=stab[i];
        while(tsnod2) {
            printf("%s\t%s\t%c\t%c\n", tsnod2->aw->aw[0]->w, tsnod2->aw->aw[1]->w, tsnod2->aw->aw[2]->w[0], tsnod2->aw->aw[3]->w[0]);
            tsnod2=tsnod2->n;
        }
        printf("\n"); 
    }
    return;
}

void freechainharr(snodm **stab, size_t tsz)
{
    int i;
    snodm *tsnod0, *tsnod2;
    for(i=0; i<tsz; ++i) {
        if( (stab[i] != NULL) ) {
            while( (stab[i]->n != NULL) ) {
                tsnod0=stab[i];
                tsnod2=stab[i]->n;
                while((tsnod2->n != NULL) ){
                    tsnod0=tsnod2;
                    tsnod2=tsnod2->n;
                }
                free(tsnod0->n);
                tsnod0->n=NULL;
            }
            free(stab[i]);
        }
    }
    free(stab);
    return;
}

snodm **hashmrk(aaw_c *maawc, unsigned tsz)
{
    unsigned i;

    snodm **stab=malloc(tsz*sizeof(snodm *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL;
    snodm *tsnod0, *tsnod2;

    unsigned tint;

    /* OK, we're going to loop through the map file container: i here follows the global SNP name index */
    for(i=0; i< maawc->numl; ++i) {

        tint=hashit(maawc->aaw[i]->aw[0]->w, tsz); // hash the snpname
        if( (stab[tint] == NULL) ) { // nothing in that slot right now.
            stab[tint]=malloc(sizeof(snodm));
            stab[tint]->aw=maawc->aaw[i];
            stab[tint]->idx=i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(snodm));
        tsnod0->n->aw = maawc->aaw[i];
        tsnod0->n->idx=i;
        tsnod0->n->n=NULL;
    }
    return stab;
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
    awc->al=awc->ab;
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
                    tgt=HH; break;
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
                    tgt=HH; break;
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
                    tgt=HH; break;
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
                    tgt=HH; break;
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
        }
    }
    printf("\n"); 
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

void prtaawc0(aaw_c *aawc, aaw_c *maawc, snodm **stab, unsigned htsz, char *tc)
{
    int i, j, pos, unseencou=0;
    int u2=0;
    unsigned tint;
    snodm *tsnod2;
    boole seen;

    for(i=2;i<aawc->numl;++i) {
        seen = 0;

        tint=hashit(aawc->aaw[i]->aw[1]->w, htsz); // hash the snpname
        if(stab[tint] == NULL) {
            u2++;
            continue;
        } else {
            tsnod2=stab[tint];
            while( (tsnod2 != NULL) ) {
                if(!(strcmp(tsnod2->aw->aw[0]->w, aawc->aaw[i]->aw[1]->w))) {
                    seen =1;
                    pos=atoi(tsnod2->aw->aw[1]->w);
                    printf("%s ", tc); 
                    printf("%s ", aawc->aaw[i]->aw[1]->w);
                    printf("%c ", '0'); 
                    printf("%i ", pos);
                    for(j=2;j<aawc->aaw[i]->al;++j)
                         printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
                    break;
                }
                tsnod2=tsnod2->n;
            }
        }
        if(!seen) {
            unseencou++;
            // printf("%s ", aawc->aaw[i]->aw[1]->w);
        }
    }
    if( u2 | unseencou)
        printf("WARNING: #u2=%i unseen=%i\n", u2, unseencou); 
    return;
}

void prtaawctp(aaw_c *aawc, aaw_c *maawc, snodm **stab, unsigned htsz, char *tc, char *fn)
{
    int i, j, pos, unseencou=0;
    int u2=0;
    unsigned tint;
    snodm *tsnod2;
    boole seen;
    FILE *fo=fopen(fn, "w");

    for(i=2;i<aawc->numl;++i) {
        seen = 0;

        tint=hashit(aawc->aaw[i]->aw[1]->w, htsz); // hash the snpname
        if(stab[tint] == NULL) {
            u2++;
            continue;
        } else {
            tsnod2=stab[tint];
            while( (tsnod2 != NULL) ) {
                if(!(strcmp(tsnod2->aw->aw[0]->w, aawc->aaw[i]->aw[1]->w))) {
                    seen =1;
                    pos=atoi(tsnod2->aw->aw[1]->w);
                    fprintf(fo, "%s ", tc); 
                    fprintf(fo, "%s ", aawc->aaw[i]->aw[1]->w);
                    fprintf(fo, "%c ", '0'); 
                    fprintf(fo, "%i ", pos);
                    for(j=2;j<aawc->aaw[i]->al;++j)
                         fprintf(fo, (j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
                    break;
                }
                tsnod2=tsnod2->n;
            }
        }
        if(!seen) {
            unseencou++;
            // printf("%s ", aawc->aaw[i]->aw[1]->w);
        }
    }
    if( u2 | unseencou)
        fprintf(stderr, "WARNING: #u2=%i unseen=%i\n", u2, unseencou); 
    fclose(fo);
    return;
}

void prt_tfc(aaw_c *tfc)
{
    int i;
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
    int i;
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
                // goto uit; /// this forces dropping out when newline is reached, and that's all, right?
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

void prtusage()
{
    printf("Usage: Insufficent arguments given.\n");
    printf("       Use options -i and -m to to give filenames of the old beagle format\n");
    printf("       and the markers (for the positions).\n");
    printf("       To output a tped, include a prefix string to -o option.\n");
    printf("       Note tfam's are not handled, due to sex desig being required.\n");
    printf("       -s will not output a tped, but only stats.\n");
    return;
}

int catchopts(optstruct *opstru, int oargc, char **oargv)
{
    int c;
    opterr = 0;
    while ((c = getopt (oargc, oargv, "sci:m:o:")) != -1)
        switch (c) {
            case 'c': // want to see the genotypes of these SNPs
                opstru->cflag = 1;
                break;
            case 's': // stats version, no tped is output. Principally to see how hets are called.
                opstru->sflag = 1;
                break;
            case 'i': // input file 
                opstru->iname = optarg;
                break;
            case 'm': // input file 
                opstru->mname = optarg; /* the tfam of fam file */
                break;
            case 'o': // input file 
                opstru->oname = optarg; /* the tfam of fam file */
                break;
			case '?':
				if ((optopt == 'i') | (optopt =='m') | (optopt == 'o') ) {
					fprintf (stderr, "Option -%c requires a string.\n", optopt);
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                printf("Sorry those options are incompatible.\n"); 
                prtusage();
                exit(EXIT_FAILURE);
        }
    return 0;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if((argc<5) | (argc>7)) {
        prtusage();
        exit(EXIT_FAILURE);
    }
    int argignore=0; //
    int oargc=argc-argignore;
    char **oargv=argv+argignore;
    optstruct opstru={0};
    catchopts(&opstru, oargc, oargv);
    if( (opstru.iname==NULL) | (opstru.mname==NULL)) {
        prtusage();
        exit(EXIT_FAILURE);
    }

    /* behold rough hack for getting chromosome number */
    char *tmp=strrchr(opstru.mname, '.');
    char *tmp2=strrchr(opstru.mname, 'r');
    char tc[3]={0};
    sprintf(tc, "%.*s", (int)(tmp-tmp2-1), tmp2+1);
    char tpedout[64]={0};
    if(opstru.oname != NULL)
        sprintf(tpedout, "%s_chr%s.tped", opstru.oname, tc);


    int i;
    FILE *fp=NULL, *fm=NULL;
    fp=fopen(opstru.iname, "r");
    fm=fopen(opstru.mname, "r");
    aaw_c *aawc=NULL, *maawc=NULL;
    aawc=processinpf(fp);
    maawc=processinpf(fm);

    unsigned htsz=givehtsz(maawc->numl);
    snodm **stab=hashmrk(maawc, htsz);
    // prtchaharr(stab, htsz);
    if((opstru.oname != NULL) & (!opstru.sflag) ) {
        prtaawctp(aawc, maawc, stab, htsz, tc, tpedout);
        goto gone;
    }

    /* Now this is a bare line-word type struct, which could conceivably be converted into a more tped-friendly structure
     * but, as usual, this woul require copying and holding two tped's in memory, so why bother? Let's just be careful and keep
     * the particularities of the tped structure in mind: */
    // prt_tpedaawc1f(aawc); /* friendly not proper tped print: CHRIST I had kep this in */
    /* First, let's start with counting the 18 categorized genotypes. */
    int *cougt=calloc(NUMGTS, sizeof(int));
    boole eqngts = 0; /* equal number of GTs in our tped file? */

    /* OK, we're sticking in some options */
    // if( (opstru.iname != NULL) & !(opstru.cflag) & (opstru.fname==NULL) & (opstru.sflag)) {
    if( opstru.sflag) {
        printf("OK in here\n"); 
        cougt_tpedaawc(aawc, &cougt, &eqngts);
        printf("Different GT counts: (Z1: just one uncalled allele, ZZ: both alleles uncalled.\n"); 
        for(i=0;i<NUMGTS;++i) 
            printf((i==NUMGTS-1)?"%4s \n":"%4s  ", gtna0[i]);
        for(i=0;i<NUMGTS;++i) 
            printf((i==NUMGTS-1)?"%5i\n":"%5i ", cougt[i]);
        printf((eqngts)?"Problem: unequal number of GTs across samples.\n":"Checked: yes, an equal num of GTs for all samples.\n");
    } else if(opstru.cflag) {
        prt_tpedaawc1p(aawc);
    }

    free(cougt);
gone:
    freechainharr(stab, htsz);
    free_aawc(&aawc);
    free_aawc(&maawc);
    // printf("overall: avga1cr=%4.4f avga2cr=%4.4f\n", allra1/numsamps, allra2/numsamps);
    fclose(fp);
    fclose(fm);

    return 0;
}
