/* size file read .... read a file with contig names and then sizez */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define GBUF 2

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
    }

typedef unsigned char boole;

typedef struct /* letcla_t, letters here are just names really */
{
    char l; // a letter
    int c; // a category
} letcla_t;

typedef struct /* dia_t */
{
    unsigned **is /* indices */, bf, sz;
    int lidx; // category label index
} dia_t; /* dupe index array */

typedef struct /* adia_t */
{
    dia_t **d;
    unsigned bf, sz;
} adia_t; /* dupe index array */

void prt_adia(adia_t *ad)
{
    int i, j;
    for(i=0;i<ad->sz;++i) {
        printf("%i: ", (*ad->d)[i].lidx);
        printf("sz=%u: ", (*ad->d)[i].sz);
        for(j=0;j<(*ad->d)[i].sz;++j) 
            printf("%u ", (*(*ad->d)[i].is)[j]);
        printf("\n"); 
    }
    return;
}

dia_t *crea_dia(void)
{
    dia_t *d=malloc(sizeof(dia_t));
    d->bf=GBUF;
    d->sz=0;
    d->is=malloc(sizeof(unsigned*));
    (*d->is)=malloc(d->bf*sizeof(unsigned));
    return d;
}

void assign_dia(dia_t **d)
{
    dia_t *dd = *d;
    dd->bf=GBUF;
    dd->sz=0;
    // dd->is=malloc(sizeof(unsigned*));
    (*dd->is)=malloc(dd->bf*sizeof(unsigned));
    *d=dd;
    return;
}

void assign_dia2(dia_t *d, unsigned lidx)
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

void reall_adia(adia_t **ad)
{
    int i;
    // dia_t *td=NULL;
    adia_t *add = *ad;
    add->bf += GBUF;
    (*add->d)=realloc((*add->d), add->bf*sizeof(dia_t));
    for(i=add->bf-GBUF;i<add->bf;++i) {
        // td=(*add->d)+i;
        // td=crea_dia();
        // (*add->d)[i]=crea_dia();
        // assign_dia((*add->d+i));
        assign_dia(add->d+i);
    }
    *ad=add;
    return;
}

void reall_adia2(adia_t *ad)
{
    int i;
    ad->bf += GBUF;
    (*ad->d)=realloc((*ad->d), ad->bf*sizeof(dia_t));
    for(i=ad->bf-GBUF;i<ad->bf;++i) {
        // td=(*add->d)+i;
        // td=crea_dia();
        // (*add->d)[i]=crea_dia();
        // assign_dia((*add->d+i));
        assign_dia2((*ad->d)+i, 8);
    }
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
        assign_dia2((*ad->d)+i, 7);
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
    // free(d);
    return;
}

void free_adia(adia_t *ad)
{
    int i;
    // for(i=0;i<add->sz;i++) // only if normalized.
    for(i=0;i<ad->bf;i++)
        free_dia((*ad->d)+i);
    free((*ad->d));
    free(ad->d);
    free(ad);
    return;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=3) {
        printf("Error. Pls supply  2 arguments 1) Size of letter array 2) Number of clusters to randomly assign.\n");
        exit(EXIT_FAILURE);
    }
    int i, j;
    int ne=atoi(argv[1]); /* num elements to be categorized, the elements will be letters */
    letcla_t *la=malloc(ne*sizeof(letcla_t));
    int nc=atoi(argv[2]); /* num clusters which will be randomly assigned: may not include all of them */

    for(i=0;i<ne;++i) {
        la[i].l=(char)(65+26*((float)rand()/RAND_MAX)); // did have 65.5 here but [ kept coming up
        la[i].c=(int)(1+nc*((float)rand()/RAND_MAX));
        printf("(%c:%i) ", la[i].l, la[i].c); 
    }
    printf("\n"); 

    adia_t *ad = crea_adia(); // seems to get past this stage OK.

    //     // unfinished:
    //     for(i=0;i<24;++i) 
    //         if(na[i]%2) {
    //             if(ad->sz == ad->bf)
    //                 reall_adia(&ad);
    //             ad->d[ad->sz]->is[ad->d[ad->sz]->sz] =na[i];
    // 
    //
    boole seenclu;
    for(i=0;i<ne;++i) {
        seenclu=0;
        for(j=0;j<ad->sz;++j) {
            if(la[i].c == (*ad->d)[j].lidx) {
                if((*ad->d)[j].sz == (*ad->d)[j].bf)
                    reall_dia((*ad->d)+j);
                (*(*ad->d)[j].is)[(*ad->d)[j].sz] = i;
                (*ad->d)[j].sz++;
                seenclu=1;
                // ad->d[ad->sz]->is[ad->d[ad->sz]->sz] =la[i].l;
            }
            if(seenclu)
                break;
        }
        // have gone through all the available clusters with no luck, time to create a new one.
        if(!seenclu) {
            if(ad->sz == ad->bf)
                reall_adia2(ad);
            (*ad->d)[ad->sz].lidx = la[i].c;
            (*(*ad->d)[ad->sz].is)[(*ad->d)[ad->sz].sz] = i;
            (*ad->d)[ad->sz].sz++;
            ad->sz++;
        }
    }

    prt_adia(ad);

    free_adia(ad);
    free(la);

    return 0;

}
