/* size file read .... read a file with contig names and then sizez */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define GBUF 2
#define WBUF 2
#define MNCOLS 4 // mandatory number of columns

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
    int pr; // index as it appear as in map file. By nature, must aligned with assoc ped file
    /* why is the above necessary? the index of the array of these should be enough?
     * or do we need something explicit for the hash table? */
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
    unsigned char gd; // genuine duplicate ... see tochainharr2
    int gdn; /* the number corresponding to the duplicate category */
    int idx; // the index corresponding to this mp element: it's the sort of thing youex
};

typedef struct strchainode snod;

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
    boole seendup; // this sepearates genuine dupes from hash table collisions.

    snod **stab=malloc(tsz*sizeof(snod *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    snod *tsnod0, *tsnod2;
    // adia_t *ad = crea_adia();

    unsigned tint;
    int gdk=0; /* the k-index for counting up the genuine duplicates */
    for(i=0; i< m; ++i) {
        seendup=0;

        tint=hashit(mp[i].nn, tsz); 

        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(snod));
            stab[tint]->mp=mp+i;
            stab[tint]->idx=i;
            stab[tint]->n=NULL;
            stab[tint]->gd=0; // first entry not a genuine dupe.
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            /// this bit is controversial, if strings are rtruly the same, leave the original in
            // however, this time I want to include dupes so set the gd
            if(!strcmp(tsnod2->mp->nn, mp[i].nn)) {
                if(!tsnod2->gd) { // new duplicate type therefore new category
                    tsnod2->gdn = gdk;
                    gdk++;
                    tsnod2->gd = 1;
                }
                if(!tsnod2->gd) { // new duplicate type therefore new category
                    loc_cat(ad, i, gdk);
                }
                tsnod2->mp->pr = tsnod2->idx;
                mp[i].pr = tsnod2->mp->pr;
                seendup =1;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(snod));
        tsnod0->n->mp=mp+i;
        tsnod0->n->idx=i;
        tsnod0->n->gd = (seendup)?1:0;
        /* not setting gdk here, the loc_cat will catch it */
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
            printf("%s(I:%i,D:%i,PR:%i) ", tsnod2->mp->nn, tsnod2->idx, (int)tsnod2->gd, tsnod2->mp->pr);
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
                    sprintf(mp[wa->numl].nn, "C%02i_P%09li", (int)mp[wa->numl].cnu, mp[wa->numl].pos);
                	mp[wa->numl].pr=-1; // default val -1 is no pairing.
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
    if(argc!=2) {
        printf("Error. Pls supply argument (name of plink map file).\n");
        exit(EXIT_FAILURE);
    }

    int i, m, n;
    mp_t *mp=processinpf(argv[1], &m, &n);

    if(n != MNCOLS) {
        printf("Error: Mandatory uniform number of columsn for plink map files is %i\n", MNCOLS); 
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
    prtchaharr2(mph, htsz);

    norm_adia(ad);
    prt_adia(ad);

    free_adia(ad);

    freechainharr(mph, htsz);

abo: 
    for(i=0;i<m;++i) {
		free(mp[i].n);
		free(mp[i].nn);
    }
    free(mp);
    printf("Run summary: nsnps=%i htsz=%u\n", m, htsz); 
    return 0;
}
