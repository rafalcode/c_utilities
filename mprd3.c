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

typedef struct /* mp_t */
{
	char *n;
	char *nn; /* numbered name CXX_XXXXX, etc. will alway sbe 16 chars in length */
	size_t nsz; /* size of the name r ID field */
    char cnu; // first column the chromosome number.
    float cmo; // the centimorgans
	long pos; /* just the one number */
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

snod **tochainharr(mp_t *mp, int m, unsigned tsz)
{
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
    if(!(htsz%2)) // prime is usually recomended, but let's go for at least odd.
        htsz++;
    snod **mph = tochainharr(mp, m, htsz);
    prtchaharr(mph, htsz);
    freechainharr(mph, htsz);

abo: 
    for(i=0;i<m;++i) {
		free(mp[i].n);
		free(mp[i].nn);
    }
    free(mp);

    return 0;
}
