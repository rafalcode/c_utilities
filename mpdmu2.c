#include "mpdmu2.h" // has the header stuff, and structs.

int catchopts(optstruct *opstru, int oargc, char **oargv)
{
    int c;
    opterr = 0;
    while ((c = getopt (oargc, oargv, "ptefcn:")) != -1)
        switch (c) {
            case 'p': // want output as ped and map
                opstru->pflag = 1;
                break;
            case 't': // want output as tped and tfam
                opstru->tflag = 1;
                break;
            case 'e': // print per-sample Tech Rep resolutions events
                opstru->eflag = 1;
                break;
            case 'f': // print per-sample Tech Rep resolutions events
                opstru->fflag = 1;
                break;
            case 'c': // "comfort" flag ... similar to -t, but want to see CN_PN's strings, "just to make sure"
                opstru->cflag = 1;
                break;
            case 'n': // n for name file flag .. a filename with SNP names
                opstru->nf=optarg;
                break;
			case '?':
				if (optopt == 'n') {
					fprintf (stderr, "Option -%c requires an file with SNP names argument.\n", optopt);
                    exit(EXIT_FAILURE);
                }
            default:
                abort();
        }
    return 0;
}

char *newna(char *fname)
{
    char *tc=strrchr(fname, '.');
    char *retc=calloc(3+strlen(fname), sizeof(char));
    sprintf(retc, "%.*s_2%s", (int)(tc-fname), fname, tc);
    return retc;

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
                default:
                    tgt=ZZ;
            }
            break;
        case 'T':
            switch(A2){
                case 'A':
                    tgt=CA; break;
                case 'C':
                    tgt=CC; break;
                case 'G':
                    tgt=CG; break;
                case 'T':
                    tgt=CT; break;
                default:
                    tgt=ZZ;
            }
            break;
        default:
            tgt=ZZ;
    }
    return tgt;
}

i2g_t2 *crea_i22(void)
{
    int i;
    i2g_t2 *i22=malloc(sizeof(i2g_t2));
    i22->bf=GBUF;
    i22->sz=0;
    i22->i=malloc(sizeof(unsigned*));
    i22->gt=malloc(sizeof(gt_t*));
    i22->dpcou=malloc(sizeof(int*));
    (*i22->i)=malloc(i22->bf*sizeof(unsigned));
    (*i22->gt)=malloc(i22->bf*sizeof(gt_t));
    (*i22->dpcou)=malloc(i22->bf*sizeof(int));
    for(i=0;i<i22->bf;++i) {
        (*i22->gt)[i]=ZZ;
        (*i22->dpcou)[i]=2;
    }
    return i22;
}

void append_i22(i2g_t2 *i22, unsigned i)
{
    int j;
    if(i22->sz == i22->bf) {
        i22->bf += GBUF;
        (*i22->i)=realloc((*i22->i), i22->bf*sizeof(unsigned));
        (*i22->gt)=realloc((*i22->gt), i22->bf*sizeof(gt_t));
        (*i22->dpcou)=realloc((*i22->dpcou), i22->bf*sizeof(int));
        for(j=i22->bf-GBUF;j<i22->bf;++j) {
            (*i22->gt)[j]=ZZ;
            (*i22->dpcou)[j]=2;
        }
    }
    (*i22->i)[i22->sz]=i; // can'tset gt s yet though.
    i22->sz++;
    return;
}

void norm_i22(i2g_t2 *i22)
{
    (*i22->i)=realloc((*i22->i), i22->sz*sizeof(unsigned));
    (*i22->gt)=realloc((*i22->gt), i22->sz*sizeof(gt_t));
    (*i22->dpcou)=realloc((*i22->dpcou), i22->sz*sizeof(int));
    return;
}

void prt_i22(i2g_t2 *i22) // try not to use this before ->gt's have been set.
{
    // this version is when you are inside a sample, for which the gt setting is unique
    int i;
    for(i=0;i<i22->sz;i++) {
        printf("%u:%s:%i ", (*i22->i)[i], gtna[(*i22->gt)[i]], (*i22->dpcou)[i]);
        printf("\n");
    }
    return;
}

void prt2_i22(i2g_t2 *i22, mp_t *mp) // try not to use this before ->gt's have been set.
{
    // this version is for when you wan tthe global details
    int i;
    printf("Printing of Technical replicates:\nGLBALSNPIdx\tSNPName\tPositionStrg\tNumDups\n\n"); 
    for(i=0;i<i22->sz;i++) {
        printf("%u\t%s\t%s\t%i\n", (*i22->i)[i], mp[(*i22->i)[i]].n, mp[(*i22->i)[i]].nn, (*i22->dpcou)[i]);
    }
    return;
}

void prt3_i22(i2g_t2 *i22, mp_t *mp, adia_t *ad) // try not to use this before ->gt's have been set.
{
    // this version is for when you wan tthe global details
    int i, j;
    unsigned tu;
    printf("Printing of Technical replicates:\n");
    // printf("GLBALSNPIdx\tSNPName\tPositionStrg\tNumDups\n\n"); 
    printf("DuplicateSet\tNumDuplicates\tChromosome\tPosition");
    printf("\tRetainedGlbSNPIdx\tRetainedSNPName\tResolEvent\tListDupIdx\tListDupNames\n");
#ifdef DBG
    printf("ad sz %i mid sz %i\n", ad->sz, i22->sz);
    // these are two are same size, but are they also in sync? They should be.
#endif
    for(j=0;j<ad->sz;++j) {
        printf("%i\t%i\t%s\t%li", (*ad->d)[j].lidx, (*i22->dpcou)[j], mp[(*i22->i)[j]].cnu, mp[(*i22->i)[j]].pos);
        // printf("\t%u\t%s\t%s\t", (*i22->i)[j], mp[(*i22->i)[j]].n, devna[deva[j]]);
        printf("\t%u\t%s\t", (*i22->i)[j], mp[(*i22->i)[j]].n);
        for(i=0;i<(*ad->d)[j].sz;++i) {
            tu = (*(*ad->d)[j].is)[i]; // temp variable to make reading easier.
            printf("%li ", mp[tu].pos);
        }
        printf("\t"); 
        for(i=0;i<(*ad->d)[j].sz;++i) {
            tu = (*(*ad->d)[j].is)[i]; // temp variable to make reading easier.
            printf("%s ", mp[tu].n);
        }
        printf("\n"); 
    }
    return;
}

void prt4_i22(i2g_t2 *i22, mp_t *mp, adia_t *ad) // like prt3_i22, but with CN_PN string also just for comfort-of-mind sake.
{
    // this version is for when you wan tthe global details
    int i, j;
    unsigned tu;
    printf("Printing of Technical replicates:\n");
    // printf("GLBALSNPIdx\tSNPName\tPositionStrg\tNumDups\n\n"); 
    printf("DuplicateSet\tNumDuplicates\tChromosome\tPosition");
    printf("\tRetainedGlbSNPIdx\tRetainedSNPName\tResolEvent\tListDupIdx\tListDupNames\n");
#ifdef DBG
    printf("ad sz %i mid sz %i\n", ad->sz, i22->sz);
    // these are two are same size, but are they also in sync? They should be.
#endif
    for(j=0;j<ad->sz;++j) {
        printf("%i\t%i\t%s\t%li", (*ad->d)[j].lidx, (*i22->dpcou)[j], mp[(*i22->i)[j]].cnu, mp[(*i22->i)[j]].pos);
        // printf("\t%u\t%s\t%s\t", (*i22->i)[j], mp[(*i22->i)[j]].n, devna[deva[j]]);
        printf("\t%u\t%s\t", (*i22->i)[j], mp[(*i22->i)[j]].n);
        for(i=0;i<(*ad->d)[j].sz;++i) {
            tu = (*(*ad->d)[j].is)[i]; // temp variable to make reading easier.
            printf("%li ", mp[tu].pos);
        }
        printf("\t"); 
        for(i=0;i<(*ad->d)[j].sz;++i) {
            tu = (*(*ad->d)[j].is)[i]; // temp variable to make reading easier.
            printf("%s ", mp[tu].n);
        }
        printf("\t"); 
        for(i=0;i<(*ad->d)[j].sz;++i) {
            tu = (*(*ad->d)[j].is)[i]; // temp variable to make reading easier.
            printf("%s ", mp[tu].nn);
        }
        printf("\n"); 
    }
    return;
}

void cleangt_i22(i2g_t2 *i22) // each time we start a new sample, we need to clean just the gt part: i and dpcou are global
{
    int i;
    for(i=0;i<i22->sz;i++)
        (*i22->gt)[i] = ZZ;
    return;
}

void free_i22(i2g_t2 *i22)
{
    free((*i22->i));
    free((*i22->gt));
    free((*i22->dpcou));
    free(i22->i);
    free(i22->dpcou);
    free(i22->gt);
    free(i22);
    return;
}

void assign_dgia(dgia_t *dg, gt_t gt)
{
    int i;
    dg->bf=GBUF;
    dg->sz=0;
    dg->gt=gt;
    dg->is=malloc(sizeof(unsigned*));
    dg->js=malloc(sizeof(unsigned*));
    (*dg->is)=malloc(dg->bf*sizeof(unsigned));
    (*dg->js)=malloc(dg->bf*sizeof(unsigned));
    for(i=0;i<dg->bf;++i) { 
        (*dg->is)[i]=9; // a dbug friendly init value
        (*dg->js)[i]=9; 
    }
    return;
}

void reall_dgia(dgia_t *dg)
{
    dg->bf += GBUF;
    (*dg->is)=realloc((*dg->is), dg->bf*sizeof(unsigned));
    (*dg->js)=realloc((*dg->js), dg->bf*sizeof(unsigned));
    return;
}

void reall_adgia(adgia_t *adg)
{
    int i;
    adg->bf += GBUF;
    (*adg->dg)=realloc((*adg->dg), adg->bf*sizeof(dgia_t));
    for(i=adg->bf-GBUF;i<adg->bf;++i)
        assign_dgia((*adg->dg)+i, 8);
    return;
}

adgia_t *crea_adgia(void)
{
    int i;
    adgia_t *adg=malloc(sizeof(adgia_t));
    adg->bf=GBUF;
    adg->sz=0;
    adg->dg=malloc(sizeof(dgia_t*));
    (*adg->dg)=malloc(adg->bf*sizeof(dgia_t));
    for(i=0;i<adg->bf;i++)
        assign_dgia((*adg->dg)+i, 7);
    return adg;
}

void free_dgia(dgia_t *dg)
{
    free((*dg->is));
    free((*dg->js));
    free(dg->is);
    free(dg->js);
    return;
}

void norm_dgia(dgia_t *dg)
{
    (*dg->is)=realloc((*dg->is), dg->sz*sizeof(unsigned));
    (*dg->js)=realloc((*dg->js), dg->sz*sizeof(unsigned));
    return;
}

void free_adgia(adgia_t *adg)
{
    int i;
    for(i=0;i<adg->sz;i++)
        free_dgia((*adg->dg)+i);
    free((*adg->dg));
    free(adg->dg);
    free(adg);
    return;
}

void norm_adgia(adgia_t *adg)
{
    int i;
    for(i=adg->sz;i<adg->bf;i++)
        free_dgia((*adg->dg)+i);
    (*adg->dg)=realloc((*adg->dg), adg->sz*sizeof(dgia_t));
    // then we want to shorten the individual arrays (dgia_t) to their proper size
    for(i=0;i<adg->sz;i++)
        norm_dgia((*adg->dg)+i);
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

void statsaawc(aaw_c *aawc) // actually printout in table fashion for all
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

void loc_catg(adgia_t *adg, int i, unsigned tu, gt_t gt) /* locate the category */
{
    int j;
    boole seencatgry=0;
    for(j=0;j<adg->sz;++j) {
        if(gt == (*adg->dg)[j].gt) {
            if((*adg->dg)[j].sz == (*adg->dg)[j].bf)
                reall_dgia((*adg->dg)+j);
            (*(*adg->dg)[j].is)[(*adg->dg)[j].sz] = i;
            (*(*adg->dg)[j].js)[(*adg->dg)[j].sz] = tu;
            (*adg->dg)[j].sz++;
            seencatgry=1;
        }
        if(seencatgry)
            break;
    }
    /* have gone through all the available clusters with no luck, time to create a new one.
     * Note, you still need to check for seencatgry ... 
     * not immediately obvious why, given break statement */
    if(!seencatgry) {
        if(adg->sz == adg->bf)
            reall_adgia(adg);
        (*adg->dg)[adg->sz].gt = gt;
        (*(*adg->dg)[adg->sz].is)[(*adg->dg)[adg->sz].sz] = i;
        (*(*adg->dg)[adg->sz].js)[(*adg->dg)[adg->sz].sz] = tu;
        (*adg->dg)[adg->sz].sz++;
        adg->sz++;
    }
    return;
}

void prt_adgia(adgia_t *adg)
{
    int i, j;
    for(i=0;i<adg->sz;++i) {
        printf("\t%s", gtna[(*adg->dg)[i].gt]);
        printf("(#:%u): ", (*adg->dg)[i].sz);
        for(j=0;j<(*adg->dg)[i].sz;++j) 
            printf("%u|%u ", (*(*adg->dg)[i].is)[j], (*(*adg->dg)[i].js)[j]);
        printf("\n"); 
    }
    return;
}

void proc_adgia3(adgia_t *adg, mp_t *mp, aaw_c *aawc, int curri, i2g_t2 *mid, int *deva, int nsidx)
{
    /* an advanced proc adgia ... for dupstats3
     * NB mid is passed through with master indices but no gt's for them */
    int i;
    int mxsz=0;
    int ocmx=0;
    boole zzpresent=0; // this dupsetg included a 00 call.
    gt_t rgt; // the GT to represent.

    if(1==adg->sz) { // the case where the techrep only has one GT variant
        if( (*adg->dg)[0].gt == ZZ) { // however, it could be ZZ all round.
            deva[nsidx*NDEV+OGTA00]++;
        } else {
            deva[nsidx*NDEV+OGTN00]++;
        }
        (*mid->gt)[curri] = (*adg->dg)[0].gt; // outside of if()? even when 00, we want this in here
#ifdef DBG2
            printf("\t\t\tSet dupset_%i (rep by globmapidx %u) to %s for this sample, all TRs give 1GT\n", curri, (*mid->i)[curri], gtna[(*mid->gt)[curri]]);
#endif
        return; //get out early
    }

    // preprocess
    for(i=0;i<adg->sz;++i) {
        // ignore but record the idz of the ZZ GT if it exists.
        if( (*adg->dg)[i].gt == ZZ) {
            zzpresent=1;
            continue;
        }
        if(mxsz < (*adg->dg)[i].sz)
            mxsz = (*adg->dg)[i].sz; 
    }
    // now we have the max appearing GT (excepting 00 of course) we don't know how often max occurs though.
    for(i=0;i<adg->sz;++i)
        if( (mxsz == (*adg->dg)[i].sz) && ( (*adg->dg)[i].gt != ZZ) ) {
            ocmx++;
            // we're only seeing if there is more than one max GT, but as we're here we'll also 
            // record it so it's ready when there's only 1 with max.
            rgt = (*adg->dg)[i].gt;
        }

    if(ocmx!=1) { // so more than one GT is of the max length
        // printf("\t\t\tInconclusive TRs: all following must be set to missing (A1=0, A2=0): ");
        // printf("\t\t\tSet to %s for this sample\n", "ZZ");
#ifdef DBG2
        printf("\t\t\tSet dupset_%i (rep by globmapidx %u) to %s for this sample\n", curri, (*mid->i)[curri], gtna[(*mid->gt)[curri]]);
#endif
        if(zzpresent)
            deva[nsidx*NDEV+MGTNWZ]++; // no winnners though there was also a 00 call.
        else 
            deva[nsidx*NDEV+MGTNWI]++; // no winnners, i.e. >1 different GTs came out on top
        /* NOTE: no need to set to 00, because that is already the default */

        // for(i=0;i<adg->sz;++i) 
        //     for(j=0;j<(*adg->dg)[i].sz;++j) 
        //         printf("%u ", (*(*adg->dg)[i].js)[j]);
        // printf("\n"); 
    // and no appending to i2 is done.
    } else {
        // printf("\t\t\tMajority rules: set idx %u to %s from following: ", sidx, gtna[rgt]);
        // printf("\t\t\tSet to %s for this sample\n", gtna[rgt]);
        // for(i=0;i<adg->sz;++i) 
        //     for(j=0;j<(*adg->dg)[i].sz;++j) 
        //         printf("%u ", (*(*adg->dg)[i].js)[j]);
        // printf("\n");
        (*mid->gt)[curri] = rgt;
#ifdef DBG2
        printf("\t\t\tSet dupset_%i (rep by globmapidx %u) to %s for this sample\n", curri, (*mid->i)[curri], gtna[(*mid->gt)[curri]]);
#endif
        // printf("SNP Idx %i with GT %s is the winning Tech Rep for this SNP.\n", sidx, gtna[uwgt]);
        if(zzpresent)
            deva[nsidx*NDEV+ZGTB1W]++; // ZZ appeared but a valid GT won
        else 
            deva[nsidx*NDEV+VGTB1W]++; // no winnners though there was also a 00 call.
    }
    return;
}

void dupstats4(aaw_c *aawc, mp_t *mp, adia_t *ad, i2g_t2 *mid, int *deva, int numsamps)
{
    char a1, a2;
    int j, jj, k;
    unsigned tu;
    size_t i;

    adgia_t *adg =NULL;
    for(i=0;i<aawc->numl;++i) {
        if(aawc->aaw[i]->aw[0]->w[0] == '#')
            continue;
#ifdef DBF
        printf("TR Sets for sample: %s\n", aawc->aaw[i]->aw[1]->w); // sample name
#endif
        for(j=0;j<ad->sz;++j) { // for each dup category, group better said.
            tu = (*(*ad->d)[j].is)[0]; // temp variable to make reading easier.
            jj=tu*2+6;
            adg=crea_adgia();
            for(k=0;k<(*ad->d)[j].sz;++k) {
                tu = (*(*ad->d)[j].is)[k]; // temp variable to make reading easier.
                jj=tu*2+6; // typical conversion of map to ped indices
                a1 = aawc->aaw[i]->aw[jj]->w[0];
                a2 = aawc->aaw[i]->aw[jj+1]->w[0];
                // you need to get the "tu" in here it's the global coordinate of the TRSNP.
                loc_catg(adg, k, tu, from2l(a1,a2));
            }
            norm_adgia(adg);
#ifdef DBG2
            prt_adgia(adg);
            printf("\tTR#%i (master index %u):\n", j, (*mid->i)[j]);
#endif
            proc_adgia3(adg, mp, aawc, j, mid, deva, numsamps);
            free_adgia(adg);
            cleangt_i22(mid);
        }
    }
    return;
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

void prt_adia(adia_t *ad)
{
    /* gives the raw output of the duplicate table, for debugging */
    int i, j;
    printf("%i duplicate categories were detected.\n", ad->sz);
    printf("What now follows is the arrangement of the elements into their corresponding categories\n\n"); 
    for(i=0;i<ad->sz;++i) {
        printf("%i: ", (*ad->d)[i].lidx);
        printf("sz=%u: ", (*ad->d)[i].sz);
        for(j=0;j<(*ad->d)[i].sz;++j) 
            printf("%u ", (*(*ad->d)[i].is)[j]);
        printf("\n"); 
    }
    return;
}

void prt_adia2(adia_t *ad, mp_t *mp) // m not actually required.
{
    int i, j;
    unsigned tu;
    printf("Duplicate-Pos-in-MAP-file Report:\n");
    for(i=0;i<ad->sz;++i) {
        printf("Dupset %i: ", (*ad->d)[i].lidx);
        for(j=0;j<(*ad->d)[i].sz;++j) {
            tu = (*(*ad->d)[i].is)[j]; // temp variable to make reading easier.
            printf("%s/%s/Idx=%u ", mp[tu].n, mp[tu].nn, tu);
        }
        printf("\n"); 
    }
    return;
}

void assign_dia(dia_t *d, unsigned lidx)
{
    int i;
    d->bf=GBUF;
    d->sz=0;
    d->lidx=lidx;
    d->is=malloc(sizeof(unsigned*));
    (*d->is)=malloc(d->bf*sizeof(unsigned));
    for(i=0;i<d->bf;++i) 
        (*d->is)[i]=9;
    return;
}

void reall_dia(dia_t *d)
{
    d->bf += GBUF;
    (*d->is)=realloc((*d->is), d->bf*sizeof(unsigned));
    return;
}

void reall_adia(adia_t *ad)
{
    int i;
    ad->bf += GBUF;
    (*ad->d)=realloc((*ad->d), ad->bf*sizeof(dia_t));
    for(i=ad->bf-GBUF;i<ad->bf;++i)
        assign_dia((*ad->d)+i, 8);
    return;
}

adia_t *crea_adia(void)
{
    int i;
    adia_t *ad=malloc(sizeof(adia_t));
    ad->bf=GBUF;
    ad->sz=0;
    ad->d=malloc(sizeof(dia_t*));
    (*ad->d)=malloc(ad->bf*sizeof(dia_t));
    for(i=0;i<ad->bf;i++) {
        assign_dia((*ad->d)+i, 7);
        // assign_dia(add->d+i);
        //     td=(*ad->d)+i;
        //     td=crea_dia();
        //     // ((*add->d)+i)=crea_dia();
    }
    return ad;
}

void free_dia(dia_t *d)
{
    free((*d->is));
    free(d->is);
    return;
}

void norm_dia(dia_t *d)
{
    (*d->is)=realloc((*d->is), d->sz*sizeof(unsigned));
    return;
}

void free_adia(adia_t *ad)
{
    int i;
    // for(i=0;i<add->sz;i++) // only if normalized.
    for(i=0;i<ad->sz;i++)
        free_dia((*ad->d)+i);
    free((*ad->d));
    free(ad->d);
    free(ad);
    return;
}

void norm_adia(adia_t *ad)
{
    int i;
    for(i=ad->sz;i<ad->bf;i++)
        free_dia((*ad->d)+i);
    (*ad->d)=realloc((*ad->d), ad->sz*sizeof(dia_t));
    // then we want to shorten the individual arrays (dia_t) to their proper size
    for(i=0;i<ad->sz;i++)
        norm_dia((*ad->d)+i);
    return;
}

void loc_cat(adia_t *ad, int i, int c) /* locate the category */
{
    int j;
    boole seencatgry=0;
    for(j=0;j<ad->sz;++j) {
        if(c == (*ad->d)[j].lidx) {
            if((*ad->d)[j].sz == (*ad->d)[j].bf)
                reall_dia((*ad->d)+j);
            (*(*ad->d)[j].is)[(*ad->d)[j].sz] = i;
#ifdef DBG2
            printf("loc_cat: %i added to dupcat %i\n", i, (*ad->d)[j].lidx); 
#endif
            (*ad->d)[j].sz++;
            seencatgry=1;
        }
        if(seencatgry)
            break;
    }
    /* have gone through all the available clusters with no luck, time to create a new one.
     * Note, you still need to check for seencatgry ... 
     * not immediately obvious why, given break statement */
    if(!seencatgry) {
        if(ad->sz == ad->bf)
            reall_adia(ad);
        (*ad->d)[ad->sz].lidx = c;
        (*(*ad->d)[ad->sz].is)[(*ad->d)[ad->sz].sz] = i;
        (*ad->d)[ad->sz].sz++;
        ad->sz++;
    }
    return;
}

snod **tochainharr2(mp_t *mp, int m, unsigned tsz, adia_t *ad, i2g_t2 *mid)
{
    // this "2" version of the function hashes on nn, the C##_P###etc names
    unsigned i;

    snod **stab=malloc(tsz*sizeof(snod *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    snod *tsnod0, *tsnod2;
    // adia_t *ad = crea_adia();

    unsigned tint;
    int gdk=1; /* the k-index for counting up the genuine duplicates, zero ignored as it means no duplicate. */
    for(i=0; i< m; ++i) {

        tint=hashit(mp[i].nn, tsz); 
        // this hash entry is new.
        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(snod));
            stab[tint]->mp=mp+i;
            stab[tint]->idx=i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            /// this bit is controversial, if strings are rtruly the same, leave the original in
            // however, this time I want to include dupes so set the gd
            if(!strcmp(tsnod2->mp->nn, mp[i].nn)) {
                if(!tsnod2->mp->gd) { // new duplicate type therefore new category
                    tsnod2->mp->gd = 2; // will be the master index for this dupset
                    tsnod2->mp->gdn = gdk;
                    loc_cat(ad, tsnod2->idx, tsnod2->mp->gdn);
                    append_i22(mid, tsnod2->idx);
                    // now take care of mp entry which exposed above as a dup
                    mp[i].gd = 1;
                    mp[i].gdn = gdk;
                    loc_cat(ad, i, mp[i].gdn);
                    gdk++;
                } else if (!mp[i].gd) { // this snod already registered as dup
                    // for several repeats, this could happen a few times.
                    // originally a simple else, but duplicates then were being repeatedly inserted
                    // checking whether gd was already set means loc_cat isn't overly called.
                    mp[i].gd = 1;
                    mp[i].gdn = tsnod2->mp->gdn;
                    tsnod2->mp->gd++; // this will give number of reps
                    (*mid->dpcou)[tsnod2->mp->gdn-2]++; // causes probs.
                    loc_cat(ad, i, mp[i].gdn);
                }
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(snod));
        tsnod0->n->mp=mp+i;
        tsnod0->n->idx=i;
        tsnod0->n->n=NULL;
    }
    return stab;
}

void prtchaharr(snod **stab, unsigned tsz)
{
    unsigned i;
    snod *tsnod2;
    for(i=0;i<tsz;++i) {
        printf("Tablepos %i: ", i); 
        tsnod2=stab[i];
        while(tsnod2) {
            // printf("\"%s\" ", tsnod2->mp->n); 
            printf("\"%s\"(%s:%i) ", tsnod2->mp->n, tsnod2->mp->nn, (int)tsnod2->mp->gd); 
            tsnod2=tsnod2->n;
        }
        printf("\n"); 
    }
    return;
}

void prtchaharr2(snod **stab, unsigned tsz)
{
    unsigned i;
    snod *tsnod2;
    for(i=0;i<tsz;++i) {
        printf("Tablepos %i: ", i); 
        tsnod2=stab[i];
        while(tsnod2) {
            printf("%s(I:%i,C:%i) ", tsnod2->mp->nn, tsnod2->idx, (int)tsnod2->mp->gdn);
            tsnod2=tsnod2->n;
        }
        printf("\n"); 
    }
    printf("\n"); 
    return;
}

void freechainharr(snod **stab, size_t tsz)
{
    int i;
    snod *tsnod0, *tsnod2;
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

wseq_t *create_wseq_t(size_t initsz)
{
    wseq_t *words=malloc(sizeof(wseq_t));
    words->wsbuf = initsz;
    words->quan = initsz;
    words->wln=calloc(words->wsbuf, sizeof(size_t));
    words->lbuf=WBUF;
    words->numl=0;
    words->wpla=calloc(words->lbuf, sizeof(size_t));
    return words;
}

void free_wseq(wseq_t *wa)
{
    free(wa->wln);
    free(wa->wpla);
    free(wa);
}

wff_t *process1c(char *fname, int *m, int *n)
{
    /* just the one wor per line ... into a string */
    FILE *fp=fopen(fname,"r");
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
    int c;
    boole inword=0;
    wseq_t *wa=create_wseq_t(GBUF);
    size_t bwbuf=WBUF;
    char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

    wff_t *wff=malloc(GBUF*sizeof(wff_t));

    while( (c=fgetc(fp)) != EOF) { /* grab a char */
        if( (c== '\n') | (c == ' ') | (c == '\t')) { /* word closing events */
            if( inword==1) { /* first word closing event */
                wa->wln[couw]=couc;
                bufword[couc++]='\0';
                bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
                /* for the struct, we want to know if it's the first word in a line, like so: */
                if(couw==oldcouw) {
                    wff[wa->numl].w=malloc(couc*sizeof(char));
                    wff[wa->numl].wsz=couc; // yes will be one bigger
                    strcpy(wff[wa->numl].w, bufword);
                }
                couc=0;
                couw++;
            }
            if(c=='#') { /* comment case */
                while( (c=fgetc(fp)) != '\n') ;
                continue;
            } else if(c=='\n') { /* end of a line */
                if(wa->numl == wa->lbuf-1) { /* enought space in our current array? */
                    wa->lbuf += WBUF;
                    wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
                    wff=realloc(wff, wa->lbuf*sizeof(wff_t));
                    memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
                }
                wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
                if(couw-oldcouw >1) {
                    printf("Error, each row cannot exceed 1 words: revise your input file\n"); 
                    /* need to release all memory too */
                    free_wseq(wa);
                    exit(EXIT_FAILURE);
                }
                oldcouw=couw; /* restart words per line count */
                wa->numl++; /* brand new line coming up */
            }
            inword=0;
        } else if(inword==0) { /* deal with first character of new word, + and - also allowed */
            if(couw == wa->wsbuf-1) {
                wa->wsbuf += GBUF;
                wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
                wff=realloc(wff, wa->wsbuf*sizeof(wff_t));
                for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
                    wa->wln[i]=0;
            }
            couc=0;
            bwbuf=WBUF;
            bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
            bufword[couc++]=c; /* no need to check here, it's the first character */
            inword=1;
        } else {
            if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
                bwbuf += WBUF;
                bufword = realloc(bufword, bwbuf*sizeof(char));
            }
            bufword[couc++]=c;
        }

    } /* end of big for statement */
    fclose(fp);
    free(bufword);

    /* normalization stage */
    wa->quan=couw;
    wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
    wff = realloc(wff, wa->quan*sizeof(wff_t)); /* normalize */
    wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

    *m= wa->numl;
    int k=wa->wpla[0];
    for(i=1;i<wa->numl;++i)
        if(k != wa->wpla[i])
            printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
    *n= k; 
    free_wseq(wa);

    return wff;
}

mp_t *processinpf(char *fname, int *m, int *n)
{
    /* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
     * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
     * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

    /* declarations */
    FILE *fp=fopen(fname,"r");
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
    int c;
    boole inword=0;
    wseq_t *wa=create_wseq_t(GBUF);
    size_t bwbuf=WBUF;
    char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

    mp_t *mp=malloc(GBUF*sizeof(mp_t));

    while( (c=fgetc(fp)) != EOF) { /* grab a char */
        if( (c== '\n') | (c == ' ') | (c == '\t')) { /* word closing events */
            if( inword==1) { /* first word closing event */
                wa->wln[couw]=couc;
                bufword[couc++]='\0';
                bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
                /* for the struct, we want to know if it's the first word in a line, like so: */
                if(couw==oldcouw) {
                    strncpy(mp[wa->numl].cnu, bufword, 3);
                } else if((couw - oldcouw) ==2) {
                    mp[wa->numl].cmo=atoi(bufword);
                } else if((couw - oldcouw) ==3) {
                    mp[wa->numl].pos=atol(bufword);
                    // we're ready to fill in nn: C%02i_P%09i
                    mp[wa->numl].nn=calloc(16, sizeof(char));
                    mp[wa->numl].gd=0; // default genuine dup category is 0, which means no dup.
                    mp[wa->numl].gdn=0; // default genuine dup category is 0, which means no dup.
                    sprintf(mp[wa->numl].nn, "C%s_P%09li", mp[wa->numl].cnu, mp[wa->numl].pos);
                } else if((couw - oldcouw) ==1) {
                    mp[wa->numl].n=malloc(couc*sizeof(char));
                    mp[wa->numl].nsz=couc;
                    strcpy(mp[wa->numl].n, bufword);
                }
                couc=0;
                couw++;
            }
            if(c=='#') { /* comment case */
                while( (c=fgetc(fp)) != '\n') ;
                continue;
            } else if(c=='\n') { /* end of a line */
                if(wa->numl == wa->lbuf-1) { /* enought space in our current array? */
                    wa->lbuf += WBUF;
                    wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
                    mp=realloc(mp, wa->lbuf*sizeof(mp_t));
                    memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
                }
                wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
                if(couw-oldcouw >4) {
                    printf("Error, each row cannot exceed 4 words: revise your input file\n"); 
                    /* need to release all memory too */
                    free_wseq(wa);
                    exit(EXIT_FAILURE);
                }
                oldcouw=couw; /* restart words per line count */
                wa->numl++; /* brand new line coming up */
            }
            inword=0;
        } else if(inword==0) { /* deal with first character of new word, + and - also allowed */
            if(couw == wa->wsbuf-1) {
                wa->wsbuf += GBUF;
                wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
                mp=realloc(mp, wa->wsbuf*sizeof(mp_t));
                for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
                    wa->wln[i]=0;
            }
            couc=0;
            bwbuf=WBUF;
            bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
            bufword[couc++]=c; /* no need to check here, it's the first character */
            inword=1;
        } else {
            if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
                bwbuf += WBUF;
                bufword = realloc(bufword, bwbuf*sizeof(char));
            }
            bufword[couc++]=c;
        }

    } /* end of big for statement */
    fclose(fp);
    free(bufword);

    /* normalization stage */
    wa->quan=couw;
    wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
    mp = realloc(mp, wa->quan*sizeof(mp_t)); /* normalize */
    wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

    *m= wa->numl;
    int k=wa->wpla[0];
    for(i=1;i<wa->numl;++i)
        if(k != wa->wpla[i])
            printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
    *n= k; 
    free_wseq(wa);

    return mp;
}

aaw_c *processinpf3(FILE *fp, int *lastchar)
{
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */;
    int c, oldc='\0';
    boole inword=0;
    unsigned lbuf=LBUF /* buffer for number of lines */, cbuf=CBUF /* char buffer for size of w_c's: reused for every word */;
    aaw_c *aawc=crea_aawc(lbuf); /* array of words per line */

    for(;;) {
        c=fgetc(fp);
        *lastchar=c;
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
                *lastchar=c;
                goto uit;
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

void prt_mp(mp_t *mp, int m, int n)
{
    int i, j;
    printf("var mp type mp_t is %i rows by %i columns and is as follows:\n", m, n); 
    for(i=0;i<m;++i) {
        for(j=0;j<n;++j) {
            if(j==1)
                printf("%s\t", mp[i].n);
            else if (j==3)
                printf("%li\n", mp[i].pos);
            else if (j==0)
                printf("%s\t", mp[i].cnu);
            else if (j==2)
                printf("%i ", mp[i].cmo);
        }
    }
    return;
}

void prt_wff(wff_t *wff, int m, int n)
{
    int i, j;
    printf("var wff type wff_t is %i rows by %i columns and is as follows:\n", m, n); 
    for(i=0;i<m;++i)
        for(j=0;j<n;++j)
                printf("%s\n", wff[i].w);
    return;
}

void prt_map(mp_t *mp, int m, char *fname)
{
    int i;
    char *ofn=newna(fname);
    FILE *of=fopen(ofn, "w");
    for(i=0;i<m;++i) {
        if(mp[i].gd == 1)
            continue;
        fprintf(of, "%s\t%s\t%i\t%li\n", mp[i].cnu, mp[i].n, mp[i].cmo, mp[i].pos);
    }
    fclose(of);
    free(ofn);
    return;
}

void prt_partped(mp_t *mp, int m, int n, i2g_t2 *mid, aaw_c *aawc, FILE *of)
{
    int i, k=0, ii;
    for(i=0;i<6;++i) 
        fprintf(of, "%s\t", aawc->aaw[0]->aw[i]->w);

    int ml1=m-1; // m less one
    for(i=0;i<ml1;++i) {
        if(mp[i].gd == 1) {
            continue;
        } else if (mp[i].gd == 2) {
            fprintf(of, "%c\t", gtna[(*mid->gt)[k]][0]);
            fprintf(of, "%c\t", gtna[(*mid->gt)[k]][1]);
            k++;
        } else {
            ii=i*2+6;
            fprintf(of, "%c\t", aawc->aaw[0]->aw[ii]->w[0]);
            fprintf(of, "%c\t", aawc->aaw[0]->aw[ii+1]->w[0]);
        }
    }
    // final genotype
    if(mp[ml1].gd == 1) {
        fseek(of, -1, SEEK_CUR);
        // fprintf(of, "\n"); 
        fputc('\n', of);
    } else if (mp[ml1].gd > 1) {
        fprintf(of, "%c\t", gtna[(*mid->gt)[k]][0]);
        fprintf(of, "%c\n", gtna[(*mid->gt)[k]][1]);
        k++;
    } else {
        ii=ml1*2+6;
        fprintf(of, "%c\t", aawc->aaw[0]->aw[ii]->w[0]);
        fprintf(of, "%c\n", aawc->aaw[0]->aw[ii+1]->w[0]);
    }
    return;
}

void prt_trrese(int *deva, int numsamps) // print tech rep resol events
{
    int i, j;

    for(i=0;i<NDEV;++i) 
        printf((i==NDEV-1)?"#%s\n":"#%s\t", devna[i]);
    for(i=0;i<numsamps;++i)
        for(j=0;j<NDEV;++j)
            printf((j==NDEV-1)?"%6i \n":"%6i \t", deva[i*NDEV+j]);

    printf("Legend:\n------\nNNNNNN, no TR event; OGTA00, one GT variant all 00;\n");
    printf("OGTN00, one GT but no 00; MGTNWI: conflicting GTs (no 00 nor winning GT, so 00)\n");
    printf("MGTNWZ: conflicting GTs, among which 00, no winner, so 00)\n");
    printf("VGTB1W, valid GTs but 1 winner; ZGTB1W, a valid GT winner, albeit 00 was among the TRs\n");

    return;
}

void prt_trrese2(int *deva, int numsamps) // print tech rep resol events
{
    int i, j;
    int smry[NDEV]={0};

    for(i=0;i<numsamps;++i)
        for(j=0;j<NDEV;++j)
            smry[j] += deva[i*NDEV+j];

    printf("Summary of Tech Rep Resolutions:\n------------------\n");
    printf("Undetermined(NNNNNN): %i\n", smry[0]);
    printf("All same value of 00(OGTAZZ): %i\n", smry[1]);
    printf("All same value, though not 00 (OGTN00): %i\n", smry[2]);
    printf("Many GTs, but no winner (alternates equal rep)(MGTNWI): %i\n", smry[3]);
    printf("Many GTs, among which 00, no winner (MGTNWZ): %i\n", smry[4]);
    printf("Valid GTs but 1 winner (VGTB1W): %i\n", smry[5]);
    printf("Valid GTs (one of which 00) and 1 winner (ZGTB1W): %i\n", smry[6]);

    return;
}

void prtusage(void)
{
    printf("Usage notice. Pls supply at least 2 arguments: 1) Name of plink map file 2) Name of plink ped file.\n");
    printf("Take good note: map file first, ped file second ... no option required.\n");
    printf("For output as a .tped and .tfam pair, add the \"-t\" option as (NB!) 3rd argument.\n"); 
    return;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc<3) {
        prtusage();
        exit(EXIT_FAILURE);
    }
    optstruct opstru={0};

    int argignore=2; //
    int oargc=argc-argignore;
    char **oargv=argv+argignore;
    catchopts(&opstru, oargc, oargv);

    int i, m, n;
    mp_t *mp=processinpf(argv[1], &m, &n);
    aaw_c *aawc=NULL; // like pedsta, to avoid the huge ped files, we're going to go line by line.

    // hash it all
    unsigned htsz=2*m/3;
    // try to grab a prime ... well just avoid 5-multiples, 3-multiples, and evens
    if(!(htsz%5)) 
        htsz++; // incrment by 1 if multiple of 5
    if(!(htsz%3)) 
        htsz++;
    if(!(htsz%2)) 
        htsz++;
    adia_t *ad=crea_adia();

    i2g_t2 *mid=crea_i22(); // master index of dupsets
    snod **mph = tochainharr2(mp, m, htsz, ad, mid);
    norm_i22(mid);
    // prtchaharr(mph, htsz);

    norm_adia(ad);
    // prt_adia(ad);
    //// got wrong end of stick with this guy
    /// dv_t *deva=calloc(NDEV, sizeof(dv_t)); // type of event in duplicate resolution, for all samples (not split per sample).
    int numsamps=0;
    int dvbuf= GBUF*NDEV;
    int *deva=calloc(dvbuf, sizeof(int)); // type of event in duplicate resolution, for all samples (not split per sample).
    FILE *of=NULL, *fp=fopen(argv[2], "r"); // OK, now handle
    int lastchar=9;
    char *ofn=NULL;

    if(opstru.pflag) {
        ofn=newna(argv[2]);
        of=fopen(ofn, "w");
    }
    for(;;) { // for each sample
        aawc=processinpf3(fp, &lastchar);
        if(lastchar==EOF) {
            free_aawc(&aawc);
            break;
        }
        dupstats4(aawc, mp, ad, mid, deva, numsamps);
        if(opstru.pflag)
            prt_partped(mp, m, n, mid, aawc, of);
        free_aawc(&aawc);
        numsamps++;
        CONDREALLOC(numsamps*NDEV, dvbuf, GBUF*NDEV, deva, int);
    }
    fclose(fp);
    if(opstru.pflag) {
        free(ofn);
        fclose(of);
    }
    deva=realloc(deva, numsamps*NDEV*sizeof(int));

    if(opstru.tflag)
        prt4_i22(mid, mp, ad);

    if(opstru.eflag)
        prt_trrese(deva, numsamps);

    if(opstru.fflag)
        prt_trrese2(deva, numsamps);

    // OK get a list of SNP names
    wff_t *wff=NULL;
    int mnf=0, nnf=0;
    if(opstru.nf) {
        wff=process1c(opstru.nf, &mnf, &nnf);
        prt_wff(wff, mnf, nnf);
    }

    // prt_adia3(ad, mp, mid);
    free_adia(ad);
    free_i22(mid);
    freechainharr(mph, htsz);

    if(opstru.pflag)
        prt_map(mp, m, argv[1]);

    for(i=0;i<m;++i) {
        free(mp[i].n);
        free(mp[i].nn);
    }
    free(mp);
    free(deva);
    if(opstru.nf) {
        for(i=0;i<mnf;++i)
            free(wff[i].w);
        free(wff);
    }

    return 0;
}
