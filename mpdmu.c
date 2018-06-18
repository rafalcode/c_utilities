/* size file read .... read a file with contig names and then sizez */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "mpdmu.h"

#ifdef DBG
#define GBUF 2
#define WBUF 2
#else
#define GBUF 64
#define WBUF 64
#endif
#define MNCOLS 4 // mandatory number of columns

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
		memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
    }

typedef unsigned char boole;

typedef struct /* i2g_t */
{
    unsigned **i, sz, bf;
} i2g_t; /* map indices to go/ filter out */

typedef struct /* dgia_t */
{
    unsigned **is /* indices */, bf, sz;
    gt_t gt; // category label index
} dgia_t; /* dupe index array */

typedef struct /* adgia_t */
{
    dgia_t **dg;
    unsigned bf, sz;
} adgia_t; /* dupe index array */

typedef struct /* dia_t */
{
    unsigned **is /* indices */, bf, sz;
    int lidx; // category label index
    char posstr[17]; /* position string, the thing we're hashing on */
} dia_t; /* dupe index array */

typedef struct /* adia_t */
{
    dia_t **d;
    unsigned bf, sz;
} adia_t; /* dupe index array */

typedef struct /* mp_t, map type, one line in the map file */
{
	char *n;
	char *nn; /* numbered name CXX_XXXXX, etc. will alway sbe 16 chars in length */
	size_t nsz; /* size of the name r ID field */
    char cnu; // first column the chromosome number.
    float cmo; // the centimorgans
	long pos; /* just the one number */
    boole gd;
    int gdn;
} mp_t; /* map type */

typedef struct /* wseq_t */
{
    size_t *wln;
    size_t wsbuf;
    size_t quan;
    size_t lbuf; /* a buffer for the number of lines */
    size_t numl; /* number of lines, i.e. rows */
    size_t *wpla; /* words per line array: the number of words on each line */
} wseq_t;

struct strchainode
{
    mp_t *mp;
    struct strchainode *n;
    int idx; // the index corresponding to this mp element: it's the sort of thing youex
};

typedef struct strchainode snod;


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

i2g_t *crea_i2(void)
{
    i2g_t *i2=malloc(sizeof(i2g_t));
    i2->bf=GBUF;
    i2->sz=0;
    i2->i=malloc(sizeof(unsigned*));
    (*i2->i)=malloc(i2->bf*sizeof(unsigned));
    return i2;
}

void append_i2(i2g_t *i2, unsigned i)
{
    if(i2->sz == i2->bf) {
        i2->bf += GBUF;
        (*i2->i)=realloc((*i2->i), i2->bf*sizeof(unsigned));
    }
    (*i2->i)[i2->sz]=i;
    i2->sz++;
    return;
}

void norm_i2(i2g_t *i2)
{
    (*i2->i)=realloc((*i2->i), i2->sz*sizeof(unsigned));
    return;
}

void prt_i2(i2g_t *i2)
{
    int i;
    for(i=0;i<i2->sz;i++)
        printf("%u ", (*i2->i)[i]);
    printf("\n"); 
    return;
}

void free_i2(i2g_t *i2)
{
    free((*i2->i));
    free(i2->i);
    free(i2);
    return;
}

void prt_adgia(adgia_t *adg)
{
    int i, j;
    printf("What now follows is the arrangement of the elements into their corresponding categories\n\n"); 
    for(i=0;i<adg->sz;++i) {
        printf("%s: ", gtna[(*adg->dg)[i].gt]);
        printf("sz=%u: ", (*adg->dg)[i].sz);
        for(j=0;j<(*adg->dg)[i].sz;++j) 
            printf("%u ", (*(*adg->dg)[i].is)[j]);
        printf("\n"); 
    }
    printf("\n"); 
    return;
}

void assign_dgia(dgia_t *dg, gt_t gt)
{
    int i;
    dg->bf=GBUF;
    dg->sz=0;
    dg->gt=gt;
    dg->is=malloc(sizeof(unsigned*));
    (*dg->is)=malloc(dg->bf*sizeof(unsigned));
    for(i=0;i<dg->bf;++i) 
        (*dg->is)[i]=9; // a dbug friendly init value
    return;
}

void reall_dgia(dgia_t *dg)
{
    dg->bf += GBUF;
    (*dg->is)=realloc((*dg->is), dg->bf*sizeof(unsigned));
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
    free(dg->is);
    return;
}

void norm_dgia(dgia_t *dg)
{
    (*dg->is)=realloc((*dg->is), dg->sz*sizeof(unsigned));
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

void dupstats(aaw_c *aawc, mp_t *mp, adia_t *ad) // find a way to print out genos
{
    char a1, a2;
    int j, jj, k;
    unsigned tu;
    size_t i;

    int numsamps =0; // because numl is not always numsamps due to # comment lines.
    for(i=0;i<aawc->numl;++i) {
        if(aawc->aaw[i]->aw[0]->w[0] == '#')
            continue;
        numsamps++;
        printf("Sample: %s\n", aawc->aaw[i]->aw[1]->w); // sample name
        // for(j=6;j<aawc->aaw[i]->al;j+=2) {
        for(j=0;j<ad->sz;++j) { // for each dup category
            printf("\tDS:%i) ", j); 
            for(k=0;k<(*ad->d)[j].sz;++k) {
                tu = (*(*ad->d)[j].is)[k]; // temp variable to make reading easier.
                jj=tu*2+6;
                a1 = aawc->aaw[i]->aw[jj]->w[0];
                a2 = aawc->aaw[i]->aw[jj+1]->w[0];
                printf("%c%c ", a1, a2);
            }
            printf("\n"); 
        }
    }
    return;
}

void dupstats2(aaw_c *aawc, mp_t *mp, adia_t *ad) // find a way to print out genos
{
    char fa1, fa2, /* first of the allele pair */ a1, a2;
    int j, jj, k;
    unsigned tu;
    size_t i;
    boole samea; /* same alleles */

    int numsamps =0; // because numl is not always numsamps due to # comment lines.
    for(i=0;i<aawc->numl;++i) {
        if(aawc->aaw[i]->aw[0]->w[0] == '#')
            continue;
        numsamps++;
        printf("Dups which conflict for sample: %s\n", aawc->aaw[i]->aw[1]->w); // sample name
        for(j=0;j<ad->sz;++j) { // for each dup category
            tu = (*(*ad->d)[j].is)[0]; // temp variable to make reading easier.
            jj=tu*2+6;
            fa1 = aawc->aaw[i]->aw[jj]->w[0];
            fa2 = aawc->aaw[i]->aw[jj+1]->w[0];
            samea=1;
            for(k=1;k<(*ad->d)[j].sz;++k) {
                tu = (*(*ad->d)[j].is)[k]; // temp variable to make reading easier.
                jj=tu*2+6;
                a1 = aawc->aaw[i]->aw[jj]->w[0];
                a2 = aawc->aaw[i]->aw[jj+1]->w[0];
                if( (fa1 != a1) | (fa2 != a2) ) 
                    samea =0;
            }
            if(samea)
                continue; // go to next dupe set.
            printf("DSNUM:%i) ", j); 
            for(k=0;k<(*ad->d)[j].sz;++k) {
                tu = (*(*ad->d)[j].is)[k]; // temp variable to make reading easier.
                jj=tu*2+6;
                a1 = aawc->aaw[i]->aw[jj]->w[0];
                a2 = aawc->aaw[i]->aw[jj+1]->w[0];
                printf("%c%c ", a1, a2);
            }
            printf("// "); 
        }
        printf("\n"); 
    }
    printf("Note: Only conflicting duplicates mentioned. For those that agree, first can be taken.\n"); 
    return;
}

void loc_catg(adgia_t *adg, int i, gt_t gt) /* locate the category */
{
    int j;
    boole seencatgry=0;
    for(j=0;j<adg->sz;++j) {
        if(gt == (*adg->dg)[j].gt) {
            if((*adg->dg)[j].sz == (*adg->dg)[j].bf)
                reall_dgia((*adg->dg)+j);
            (*(*adg->dg)[j].is)[(*adg->dg)[j].sz] = i;
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
        (*adg->dg)[adg->sz].sz++;
        adg->sz++;
    }
    return;
}

void proc_adgia(adgia_t *adg, i2g_t *i2)
{
    int i;
    int mxsz=0;
    int ocmx=0;
    unsigned sidx; /* If there's a winner TRSNP, this is hte first available SNP index for it */

    if(1==adg->sz) {
        /* if there's only one GT then, irrespective how many times it
         * appears, the first index registering it should be given */
       // if(1==(*adg->dg)[0].sz) { // this check not nec.
       // printf("Not a duplicate, only one version of this SNP AND only one GT for it.\n"); 
        sidx = (*(*adg->dg)[0].is)[0]; /* the first TR of these will do (they're all the same). */
        append_i2(i2, sidx); // the first index.
        return;
    }

    // get the maximum represented GT
    for(i=0;i<adg->sz;++i)
        if(mxsz < (*adg->dg)[i].sz)
            mxsz = (*adg->dg)[i].sz;

    for(i=0;i<adg->sz;++i)
        if(mxsz == (*adg->dg)[i].sz) {
            ocmx++;
            // we're only seeing if there is more than one max GT, but as we're here we'll also 
            // record it so it's ready when there's only 1 with max.
            sidx = (*(*adg->dg)[i].is)[0]; /* the first TR of these will do (they're all the same). */
            // append_i2(i2, sidx + strint);
        }

    if(ocmx!=1)
        printf("Conflicting GT of equal count weight: All Tech Reps for this SNP must be removed.\n"); 
        // and no appending to i2 is done.
    else
        append_i2(i2, sidx); // the first index.
        // printf("SNP Idx %i with GT %s is the winning Tech Rep for this SNP.\n", sidx, gtna[uwgt]);

    return;
}

void dupstats3(aaw_c *aawc, mp_t *mp, adia_t *ad) // grab the best TRs
{
    char a1, a2;
    int j, jj, k;
    unsigned tu;
    size_t i;
    // gt_t *gta=NULL;

    int numsamps =0; // because numl is not always numsamps due to # comment lines.
    adgia_t *adg =NULL;
    i2g_t *i2=NULL;
    for(i=0;i<aawc->numl;++i) {
        if(aawc->aaw[i]->aw[0]->w[0] == '#')
            continue;
        numsamps++;
        i2=crea_i2();
        printf("Dups which conflict for sample: %s\n", aawc->aaw[i]->aw[1]->w); // sample name
        for(j=0;j<ad->sz;++j) { // for each dup category, group better said.
            tu = (*(*ad->d)[j].is)[0]; // temp variable to make reading easier.
            jj=tu*2+6;
            // gta=malloc((*ad->d)[j].sz*sizeof(gt_t));
            adg=crea_adgia();
            for(k=0;k<(*ad->d)[j].sz;++k) {
                tu = (*(*ad->d)[j].is)[k]; // temp variable to make reading easier.
                jj=tu*2+6;
                a1 = aawc->aaw[i]->aw[jj]->w[0];
                a2 = aawc->aaw[i]->aw[jj+1]->w[0];
                loc_catg(adg, k, from2l(a1,a2));
                // gta[k]=from2l(a1, a2);
            }
            norm_adgia(adg);
            prt_adgia(adg);
            proc_adgia(adg, i2);
            free_adgia(adg);
        }
        norm_i2(i2);
        prt_i2(i2);
        free_i2(i2);
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
    printf("\n"); 
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
            printf("%s/Idx=%u ", mp[tu].nn, tu);
        }
        printf("\n"); 
    }
    printf("\n"); 
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

snod **tochainharr(mp_t *mp, int m, unsigned tsz)
{
    // this version of thefunction will has on the name
    // look for the "2" version 
    unsigned i;

    snod **stab=malloc(tsz*sizeof(snod *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    snod *tsnod0, *tsnod2;

    unsigned tint;
    for(i=0; i< m; ++i) {
        tint=hashit(mp[i].n, tsz); 

        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(snod));
            stab[tint]->mp=mp+i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            if(!strcmp(tsnod2->mp->n, mp[i].n)) {
                goto nxt;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(snod));
        tsnod0->n->mp=mp+i;
        tsnod0->n->n=NULL;
nxt:        continue;
    }
    return stab;
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


snod **tochainharr2(mp_t *mp, int m, unsigned tsz, adia_t *ad)
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
                    tsnod2->mp->gd = 1;
                    tsnod2->mp->gdn = gdk;
                    loc_cat(ad, tsnod2->idx, tsnod2->mp->gdn);
                    mp[i].gd = 1;
                    mp[i].gdn = gdk;
                    loc_cat(ad, i, mp[i].gdn);
                    gdk++;
                } else {
                    // for several repeats, this could happen a few times.
                    mp[i].gd = 1;
                    mp[i].gdn = tsnod2->mp->gdn;
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
            printf("\"%s\"(%s) ", tsnod2->mp->n, tsnod2->mp->nn); 
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
                	mp[wa->numl].cnu=(char)atoi(bufword);
                } else if((couw - oldcouw) ==2) {
                	mp[wa->numl].cmo=atof(bufword);
                } else if((couw - oldcouw) ==3) {
                	mp[wa->numl].pos=atol(bufword);
                    // we're ready to fill in nn: C%02i_P%09i
                	mp[wa->numl].nn=calloc(16, sizeof(char));
                	mp[wa->numl].gd=0; // default genuine dup category is 0, which means no dup.
                	mp[wa->numl].gdn=0; // default genuine dup category is 0, which means no dup.
                    sprintf(mp[wa->numl].nn, "C%02i_P%09li", (int)mp[wa->numl].cnu, mp[wa->numl].pos);
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
	printf("couw %zu\n", couw); 
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

aaw_c *processinpf2(char *fname)
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

void prt(mp_t *mp, int m, int n)
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
            	printf("%i\t", (char)mp[i].cnu);
			else if (j==2)
            	printf("%2.1f ", mp[i].cmo);
		}
    }
	return;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=3) {
        printf("Error. Pls supply 2 arguments: 1) Name of plink map file 2) Name of plink ped file.\n");
        exit(EXIT_FAILURE);
    }

    int i, m, n;
    mp_t *mp=processinpf(argv[1], &m, &n);
    aaw_c *aawc=processinpf2(argv[2]);

    if(n != MNCOLS) {
        printf("Error: Mandatory uniform number of columns for plink map files is %i\n", MNCOLS); 
        goto abo;
    }

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
    snod **mph = tochainharr2(mp, m, htsz, ad);

    norm_adia(ad);
#ifdef DBG
    prt_adia(ad);
#else
    dupstats3(aawc, mp, ad);
//    prt_adia2(ad, mp);
#endif

    free_adia(ad);
    freechainharr(mph, htsz);

abo: 
    free_aawc(&aawc);
    for(i=0;i<m;++i) {
		free(mp[i].n);
		free(mp[i].nn);
    }
    free(mp);
    return 0;
}
