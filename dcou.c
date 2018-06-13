/* size file read .... read a file with contig names and then sizez */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define GBUF 2
#define WBUF 2

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
    }

typedef unsigned char boole;

typedef struct /* dia_t */
{
    unsigned **is /* indices */, bf, sz;
} dia_t; /* dupe index array */

typedef struct /* adia_t */
{
    dia_t **d;
    unsigned bf, sz;
} adia_t; /* dupe index array */

dia_t *crea_dia(void)
{
    dia_t *d=malloc(sizeof(dia_t));
    d->bf=GBUF;
    d->sz=0;
    d->is=malloc(d->bf*sizeof(unsigned));
    return d;
}

void reall_dia(dia_t **d)
{
    dia_t *dd = *d;
    dd->bf += GBUF;
    dd->is=realloc(dd->is, dd->bf*sizeof(unsigned));
    *d=dd;
    return;
}

void reall_adia(adia_t **ad)
{
    int i;
    adia_t *add = *ad;
    add->bf += GBUF;
    add->d=realloc(add->d, add->bf*sizeof(dia_t*));
    for(i=add->bf-GBUF;i<add->bf;++i) 
        add->d[i]=crea_dia();
    *ad=add;
    return;
}

adia_t *crea_adia(void)
{
    int i;
    adia_t *ad=malloc(sizeof(adia_t));
    ad->bf=GBUF;
    ad->sz=0;
    ad->d=malloc(ad->bf*sizeof(dia_t*));
    for(i=0;i<ad->bf;i++)
        ad->d[i]=crea_dia();
    return ad;
}

void free_dia(dia_t **d)
{
    dia_t *dd = *d;
    free(dd->is);
    free(dd);
    return;
}

void free_adia(adia_t **ad)
{
    int i;
    adia_t *add=*ad;
    // for(i=0;i<add->sz;i++) // only if normalized.
    for(i=0;i<add->bf;i++)
        free_dia(add->d+i);
    free(add->d);
    free(add);
    return;
}

int main(int argc, char *argv[])
{
    int na[12]={3,2 ,5, 2, 6,5, 2, 1,2 ,7, 8, 7};

    adia_t *ad = crea_adia();
    reall_adia(&ad);
    reall_adia(&ad);
    free_adia(&ad);

    return 0;

}
