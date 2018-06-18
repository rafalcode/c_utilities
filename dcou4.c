/* size file read .... read a file with contig names and then sizez */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "mpdmu.h"

#define GBUF 2

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

void prt_i2(i2g_t *i2)
{
    int i;
    for(i=0;i<i2->sz;i++)
        printf("%u ", (*i2->i)[i]);
    printf("\n"); 
    return;
}

void norm_i2(i2g_t *i2)
{
    (*i2->i)=realloc((*i2->i), i2->sz*sizeof(unsigned));
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

void proc_adgia(adgia_t *adg, i2g_t *i2, unsigned strint, unsigned szint)
{
    int i;
    int mxsz=0;
    int ocmx=0;
    unsigned sidx; /* If there's a winner TRSNP, this is hte first available SNP index for it */
    gt_t uwgt; /* the unique "winning" genotype */

    if(1==adg->sz) {
        /* if there's only one GT then, irrespective how many times it
         * appears, the first index registering it should be given */
       // if(1==(*adg->dg)[0].sz) { // this check not nec.
       // printf("Not a duplicate, only one version of this SNP AND only one GT for it.\n"); 
        sidx = (*(*adg->dg)[0].is)[0]; /* the first TR of these will do (they're all the same). */
        append_i2(i2, strint + sidx); // the first index.
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
            uwgt = (*adg->dg)[i].gt;
            sidx = (*(*adg->dg)[i].is)[0]; /* the first TR of these will do (they're all the same). */
            // append_i2(i2, sidx + strint);
        }

    if(ocmx!=1)
        printf("Conflicting GT of equal count weight: All Tech Reps for this SNP must be removed.\n"); 
        // and no appending to i2 is done.
    else
        append_i2(i2, strint + sidx); // the first index.
        // printf("SNP Idx %i with GT %s is the winning Tech Rep for this SNP.\n", sidx, gtna[uwgt]);

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

void loc_cat(adgia_t *adg, int i, gt_t gt) /* locate the category */
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

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=3) {
        printf("Error. Pls supply 2 arguments: 1) number GT dup-events to run. 2) Max dups per event\n");
        exit(EXIT_FAILURE);
    }
    int i, j;
    int ne=atoi(argv[1]); /* num elements to be categorized, the elements will be letters */
    int mne=atoi(argv[2]); /* num elements to be categorized, the elements will be letters */

    gt_t **gta=malloc(ne*sizeof(gt_t*));
    unsigned *gtaz=malloc(ne*sizeof(unsigned)); /* sizes for each */
    unsigned *gtacz=malloc((ne+1)*sizeof(unsigned)); /* sizes for each, cumulative */
    /* the cumulative fixes a global index for all the GTs and decides which one should be kept. */
    gtacz[0]=0;
    for(i=0;i<ne;++i) {
        gtaz[i] = 1+mne*((float)rand()/RAND_MAX);
        gtacz[i+1] = gtaz[i] + gtacz[i];
        gta[i]=malloc(gtaz[i]*sizeof(size_t));
        for(j=0;j<gtaz[i];++j) 
            gta[i][j]=(gt_t)(17*((float)rand()/RAND_MAX)); // did have 65.5 here but [ kept coming up
    }

    printf("Random GTs generated:\n");
    for(i=0;i<ne;++i) {
        for(j=0;j<gtaz[i];++j) 
            printf("%s ", gtna[gta[i][j]]); 
        printf("\n"); 
    }

    adgia_t *adg =NULL;
    i2g_t *i2=crea_i2();
    unsigned strint;
   
    for(i=0;i<ne;++i) {
        adg=crea_adgia();
        for(j=0;j<gtaz[i];++j)
            loc_cat(adg, j, gta[i][j]);
        norm_adgia(adg);
        prt_adgia(adg);
        strint=gtacz[i];
        proc_adgia(adg, i2, strint, szint);
        free_adgia(adg);
    }

    norm_i2(i2);
    prt_i2(i2);

    free_i2(i2);
    for(i=0;i<ne;++i)
        free(gta[i]);
    free(gta);
    free(gtaz);
    free(gtacz);

    return 0;
}
