#include "bglvset.h"

int catchopts(optstruct *opstru, int oargc, char **oargv)
{
    int c;
    opterr = 0;
    while ((c = getopt (oargc, oargv, "pci:m:o:v:")) != -1)
        switch (c) {
            case 'c': // no printing out of files, please. I.e. _C_ancel printing to file.
                opstru->cflag = 1;
                break;
            case 's': // stats version, no tped is output. Principally to see how hets are called.
                opstru->sflag = 1;
                break;
            case 'p': // confirm you are feeding phased GTs, this will set | instead of / in VCF.
                opstru->pflag = 1;
                break;
            case 'i': // input file 
                opstru->iname = optarg;
                break;
            case 'v': // input file 
                opstru->vname = optarg;
                break;
            case 'm': // input file 
                opstru->mname = optarg; /* the tfam of fam file */
                break;
            case 'o': // input file 
                opstru->oname = optarg; /* the tfam of fam file */
                break;
			case '?':
				if ((optopt == 'i') | (optopt =='m') | (optopt == 'o') | (optopt == 'v')) {
					fprintf (stderr, "Option -%c requires a string.\n", optopt);
                    exit(EXIT_FAILURE);
                }
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

snodm **hashmrk2(aaw_c *maawc, unsigned tsz)
{
    unsigned i;

    snodm **stab=malloc(tsz*sizeof(snodm *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL;
    snodm *tsnod0, *tsnod2;

    unsigned tint;

    /* OK, we're going to loop through the map file container: i here follows the global SNP name index */
    for(i=0; i< maawc->numl; ++i) {
        if(maawc->aaw[i]->gd == 2) // it's a dupe
            continue;

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

snodm **hashmrk_nd(aaw_c *maawc, unsigned tsz) /* throw out positional dups */
{
    unsigned i;

    snodm **stab=malloc(tsz*sizeof(snodm *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL;
    snodm *tsnod0, *tsnod2;
    boole pos_seen; // this position already seen.

    unsigned tint;

    /* OK, we're going to loop through the map file container: i here follows the global SNP name index */
    for(i=0; i< maawc->numl; ++i) {

        pos_seen=0;
        tint=hashit(maawc->aaw[i]->aw[1]->w, tsz); // hash the position
        if( (stab[tint] == NULL) ) { // nothing in that slot right now.
            stab[tint]=malloc(sizeof(snodm));
            stab[tint]->aw=maawc->aaw[i];
            stab[tint]->idx=i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            if( !(strcmp(tsnod2->aw->aw[1]->w, maawc->aaw[i]->aw[1]->w))) {
                tsnod2->aw->gd=1; // means this already hashed line was first but has a dup
                maawc->aaw[i]->gd=2; // means this unhashed line was second (or later) to a first dup.
                pos_seen=1;
                break;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        if(pos_seen)
            continue;
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

void prtaawc00(aaw_c *aawc, aaw_c *maawc, snodm **stab, aaw_c *vaawc, unsigned htsz, char *tc)
{
    /* prints the BGL but with VCF interpolations ... a debugging function */
    int i, j, pos, unseencou=0;
    int u2=0;
    unsigned tint;
    snodm *tsnod2;
    boole seen;
    boole switchvgt;
    int vi, vj;
    char va1, va2; // these are a1 and a2 in the VCF file

    for(i=2;i<aawc->numl;++i) {
        vi=i+5;
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
                    va1 = vaawc->aaw[vi]->aw[3]->w[0];
                    va2 = vaawc->aaw[vi]->aw[4]->w[0];
                    printf("%s ", vaawc->aaw[vi]->aw[3]->w);
                    printf("%s ", vaawc->aaw[vi]->aw[4]->w);
                    for(j=2;j<aawc->aaw[i]->al;++j) {
                        printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
                        vj=j/2+8;
                        switchvgt=0; // default don't switch VFCs GT
                        if(j%2==1) {
                            if(aawc->aaw[i]->aw[j]->w[0] == va1)
                                switchvgt=1;
                            if(switchvgt) 
                                printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s ", "1/0");
                            else
                                printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s ", vaawc->aaw[vi]->aw[vj]->w);
                        }
                    }
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

void prt_vaawc(aaw_c *vaawc)
{
    int i, j;
    char fc; // this is the first char, useful to keep n case of distingusihig comments/vcf metadata.
    for(i=0;i<vaawc->numl;++i) {
        fc = vaawc->aaw[i]->aw[0]->w[0];
        if(fc != '#') {
            for(j=0;j<vaawc->aaw[i]->al;++j)
                printf((j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", vaawc->aaw[i]->aw[j]->w);
        }
    }
    return;
}

void prt_vcf(aaw_c *vaawc, char *fname)
{
    /* print out a new vcf file ... _2 appended to its basename */
    int i, j;
    char fc, sc; /* first char, second char: needed because sometimes VCF comen sneed tabs, othertimes not. */
    FILE *fvp=fopen(fname, "w");
    for(i=0;i<vaawc->numl;++i) {
        fc = vaawc->aaw[i]->aw[0]->w[0];
        sc = vaawc->aaw[i]->aw[0]->w[1];
        for(j=0;j<vaawc->aaw[i]->al;++j) {
            if(sc != '#')
                fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", vaawc->aaw[i]->aw[j]->w);
            else
                fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s ":"%s\n", vaawc->aaw[i]->aw[j]->w);
        }
    }
    fclose(fvp);
    return;
}

void prtvaawcvcf(aaw_c *aawc, aaw_c *maawc, snodm **stab, aaw_c *vaawc, unsigned htsz, char *tc, char *fname)
{
    /* function not robust to different prints the BGL but with VCF interpolations ... a debugging function */
    int i, j;
    int u2=0;
    unsigned tint;
    unsigned commlines=0; // number of comments lines
    snodm *tsnod2;
    int bi, bj, bj2;
    char va1, va2; // these are a1 and a2 in the VCF file
    char ba1, ba2; // these are a1 and a2 in the BGL file
    char fc, sc;

    FILE *fvp=fopen(fname, "w");
    for(i=0;i<vaawc->numl;++i) {
        fc = vaawc->aaw[i]->aw[0]->w[0];
        sc = vaawc->aaw[i]->aw[0]->w[1];
        if(fc == '#') { // a comment, so print VCF as is.
            commlines++;
            for(j=0;j<vaawc->aaw[i]->al;++j) {
                if(sc != '#')
                    fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", vaawc->aaw[i]->aw[j]->w);
                else
                    fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s ":"%s\n", vaawc->aaw[i]->aw[j]->w);
            }
        } else {
            bi=i+2-commlines;
            // seen = 0;
            tint=hashit(aawc->aaw[bi]->aw[1]->w, htsz); // hash the snpname in the BGL file
            // printf("BGL hash on %s\n", aawc->aaw[bi]->aw[1]->w);
            if(stab[tint] == NULL) {
                u2++;
                continue;
            } else {
                tsnod2=stab[tint];
                while( (tsnod2 != NULL) ) {
                    if(!(strcmp(tsnod2->aw->aw[0]->w, aawc->aaw[bi]->aw[1]->w))) {
                        // seen =1;
                        va1 = vaawc->aaw[i]->aw[3]->w[0];
                        va2 = vaawc->aaw[i]->aw[4]->w[0];
                        // fprintf(fvp, "VGT:%c%c ", va1, va2); 
                        for(j=0;j<9;++j) // print straight VCF for first 9
                            fprintf(fvp, "%s\t", vaawc->aaw[i]->aw[j]->w);

                        for(j=9;j<vaawc->aaw[i]->al;++j) {
                            bj=2*(j-8);
                            bj2=bj+1;
                            ba1 = aawc->aaw[bi]->aw[bj]->w[0];
                            ba2 = aawc->aaw[bi]->aw[bj2]->w[0];
                            // fprintf(fvp, "BGT:%c%c ", ba1, ba2); 
                            // if((ba1 == va1) & (ba2 == va2)) // heterozygous in order of va1 and va2
                            //     fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", vaawc->aaw[i]->aw[j]->w);
                            // else if( (ba2 == va1) & (ba1 == va2) ) // heterzygous in opposite order to va1, va2
                            if( (ba2 == va1) & (ba1 == va2) ) // heterzygous in opposite order to va1, va2
                                fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", "1/0");
                            else // must be homozygous or in order heterozygous
                                fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", vaawc->aaw[i]->aw[j]->w);
                        }
                        break;
                    }
                    tsnod2=tsnod2->n;
                } // given an occupied location, go through all located SNPnames.
            } // finished seeing if there's a hash location for marker on that SNP name
        } // finished distinguishing between VCF comments and data
    } // finished groing through all the lines in VCF
    fclose(fvp);
    return;
}

void prtvaawcvcfp(aaw_c *aawc, aaw_c *maawc, snodm **stab, aaw_c *vaawc, unsigned htsz, char *tc, char *fname)
{
    /* phased version */
    int i, j;
    int u2=0;
    unsigned tint;
    unsigned commlines=0; // number of comments lines
    snodm *tsnod2;
    int bi, bj, bj2;
    char va1, va2; // these are a1 and a2 in the VCF file
    char ba1, ba2; // these are a1 and a2 in the BGL file
    char fc, sc;

    FILE *fvp=fopen(fname, "w");
    for(i=0;i<vaawc->numl;++i) {
        fc = vaawc->aaw[i]->aw[0]->w[0];
        sc = vaawc->aaw[i]->aw[0]->w[1];
        if(fc == '#') { // a comment, so print VCF as is.
            commlines++;
            for(j=0;j<vaawc->aaw[i]->al;++j) {
                if(sc != '#')
                    fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", vaawc->aaw[i]->aw[j]->w);
                else
                    fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s ":"%s\n", vaawc->aaw[i]->aw[j]->w);
            }
        } else {
            bi=i+2-commlines;
            // seen = 0;
            tint=hashit(aawc->aaw[bi]->aw[1]->w, htsz); // hash the snpname in the BGL file
            // printf("BGL hash on %s\n", aawc->aaw[bi]->aw[1]->w);
            if(stab[tint] == NULL) {
                u2++;
                continue;
            } else {
                tsnod2=stab[tint];
                while( (tsnod2 != NULL) ) {
                    if(!(strcmp(tsnod2->aw->aw[0]->w, aawc->aaw[bi]->aw[1]->w))) {
                        // seen =1;
                        va1 = vaawc->aaw[i]->aw[3]->w[0];
                        va2 = vaawc->aaw[i]->aw[4]->w[0];
                        // fprintf(fvp, "VGT:%c%c ", va1, va2); 
                        for(j=0;j<9;++j) // print straight VCF for first 9
                            fprintf(fvp, "%s\t", vaawc->aaw[i]->aw[j]->w);

                        for(j=9;j<vaawc->aaw[i]->al;++j) {
                            bj=2*(j-8);
                            bj2=bj+1;
                            ba1 = aawc->aaw[bi]->aw[bj]->w[0];
                            ba2 = aawc->aaw[bi]->aw[bj2]->w[0];
                            // fprintf(fvp, "BGT:%c%c ", ba1, ba2); 
                            // if((ba1 == va1) & (ba2 == va2)) // heterozygous in order of va1 and va2
                            //     fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", vaawc->aaw[i]->aw[j]->w);
                            // else if( (ba2 == va1) & (ba1 == va2) ) // heterzygous in opposite order to va1, va2
                            if( (ba2 == va1) & (ba1 == va2) ) // heterzygous in opposite order to va1, va2
                                fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", "1|0");
                            else if( (ba1 == va1) & (ba2 == va1) ) // homoz on a1
                                fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", "0|0");
                            else if( (ba1 == va2) & (ba2 == va2) ) // homoz on a1
                                fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", "1|1");
                            else if( (ba1 == va1) & (ba2 == va2) ) // heteroz on a1
                                fprintf(fvp, (j!=vaawc->aaw[i]->al-1)?"%s\t":"%s\n", "0|1");
                        }
                        break;
                    }
                    tsnod2=tsnod2->n;
                } // given an occupied location, go through all located SNPnames.
            } // finished seeing if there's a hash location for marker on that SNP name
        } // finished distinguishing between VCF comments and data
    } // finished groing through all the lines in VCF
    fclose(fvp);
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

aaw_c *processinpfv(FILE *fp) /* special for the vcf file */
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
        // if( (!couw) & (!couc) & (c=='#')) { // use this if you cwant to zap comment lines entirely, but you shouldn't
        //       while( (c=fgetc(fp)) != '\n') ;
        //       continue; // c =='\n' reached, go back an pick next c
        // }

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
    printf("       bglvset take a bgl file, a marker file and a vcf file\n");
    printf("       and rewrites the vcf with heterozygotes reflecting the bgl.\n");
    printf("       Note that it is a post-processing step intended to correct an\n");
    printf("       already generated (i.e. converted from bgl via plink1.9) vcf file\n");
    printf("       -i the bgl file.\n");
    printf("       -m the marker file.\n");
    printf("       -v the vcf file.\n");
    printf("       -p set this flag if your GTs are all phased. Other wise unphased VCFs output.\n");
    printf("       -c set this flag if you don't files written out.\n");
    printf("       Example run:\n");
    printf("       ./bglvset -p -i phased.chr24.plink.bgl.phased -m markers.chr24.txt -v pha_chr24.vcf\n");
    return;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc==1) {
        prtusage();
        exit(EXIT_FAILURE);
    }
    int argignore=0; //
    int oargc=argc-argignore;
    char **oargv=argv+argignore;
    optstruct opstru={0};
    catchopts(&opstru, oargc, oargv);
    if( (opstru.iname==NULL) & (opstru.mname==NULL) & (opstru.vname==NULL)) {
        prtusage();
        exit(EXIT_FAILURE);
    }

    /* behold rough hack for getting chromosome number */
    char *tmp=NULL, *tmp2=NULL;
    char tc[3]={0};
    char tpedout[64]={0};
    char vcfout[64]={0};


    int i;
    FILE *fp=NULL, *fm=NULL, *fv=NULL;
    fp=fopen(opstru.iname, "r");
    fm=fopen(opstru.mname, "r");
    aaw_c *aawc=NULL, *maawc=NULL, *vaawc=NULL;
    if((opstru.vname != NULL) & (opstru.iname != NULL) & (opstru.mname != NULL)) {
        printf("Options OK.\n"); 
    } else {
        printf("Sorry current only -i, -m and -v options active run like so:\n"); 
        printf("./bglvset -i phased.chr24.plink.bgl.phased -m markers.chr24.txt -v pha_chr24.vcf\n");
        printf("Use -c (cancel) flag if you don't want files written out.\n");
        exit(EXIT_FAILURE);
    }
    fv=fopen(opstru.vname, "r");
    vaawc=processinpfv(fv);
    // prt_vcf(vaawc);
    tmp=strrchr(opstru.mname, '.');
    tmp2=strrchr(opstru.mname, 'r');
    sprintf(tc, "%.*s", (int)(tmp-tmp2-1), tmp2+1);
    if((opstru.oname != NULL) & (!opstru.cflag))
        sprintf(tpedout, "%s_chr%s.tped", opstru.oname, tc);
    else if(!opstru.cflag) {
        tmp=strrchr(opstru.vname, '.');
        sprintf(vcfout, "%.*s_2.vcf", (int)(tmp-opstru.vname), opstru.vname);
    }
    aawc=processinpf(fp);
    maawc=processinpf(fm);

    unsigned htsz=givehtsz(maawc->numl);

    // we want to clear dupes, so first we has on position
    snodm **stab=hashmrk_nd(maawc, htsz);
    // this will set the gd's on maawc's.
    // // but we no longer need this hash
    freechainharr(stab, htsz);
    // we can reuse stab
    stab=hashmrk2(maawc, htsz); // this second version of function takes gd's into account.
    // prtchaharr(stab, htsz);
    if((opstru.oname != NULL) & !(opstru.cflag)) { // cflag prevents files being written
        prtaawctp(aawc, maawc, stab, htsz, tc, tpedout);
    } else if(opstru.cflag) {
        prtaawc00(aawc, maawc, stab, vaawc, htsz, tc);
    } else {
        // prt_vcf(vaawc, vcfout);
        if(opstru.pflag)
            prtvaawcvcfp(aawc, maawc, stab, vaawc, htsz, tc, vcfout);
        else
            prtvaawcvcf(aawc, maawc, stab, vaawc, htsz, tc, vcfout);
    }

    freechainharr(stab, htsz);
    free_aawc(&vaawc);
    free_aawc(&aawc);
    free_aawc(&maawc);
    // printf("overall: avga1cr=%4.4f avga2cr=%4.4f\n", allra1/numsamps, allra2/numsamps);
    fclose(fp);
    fclose(fm);
    fclose(fv);

    return 0;
}
