/* Failed attempt to turn this into a row reader. Based on matread.c */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define GBUF 4
#define WBUF 4

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
    }

typedef unsigned char boole;

typedef struct /*bgr_t*/
{
	char *n;
	long c[3]; /* coords: 1) start 2) end 3) coverage */
} bgr_t; /* bedgraph row type */

typedef struct /*bgra_t*/
{
	bgr_t *bg;
	int bf, sz;
} bgr_t; /* bedgraph row type */

typedef struct /* wseq_t */
{
    size_t *wln;
    size_t wsbuf;
    size_t quan;
    size_t lbuf; /* a buffer for the number of lines */
    size_t numl; /* number of lines, i.e. rows */
    size_t *wpla; /* words per line array: the number of words on each line */
} wseq_t;

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

bgr_t *processinpf(char *fname, int *m, int *n)
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

	int bgrabf=GBUF, cc=0;
    bgra_t *bgra=malloc(bgrabf*sizeof(bgra_t));
	for(i=0;i<bgrabf; i++) {
		bgra[i].bf=GBUF;
		bgra[i].sz=0;
    	bgra[i].bg=malloc(bgra[i].bf*sizeof(bgr_t));
	}

	char *tmpch=NULL; /* temp string holder */
	int prevlc=0; /* previous line count */
    while( (c=fgetc(fp)) != EOF) {
        /*  take care of  */
        if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) {
            if( inword==1) { /* we end a word */
                wa->wln[couw]=couc;
                bufword[couc++]='\0';
                bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				if(strcmp(tmpch, bufword)) { /* a new bgra element required. */
					CONDREALLOC(cc, bgrabf, GBUF, bgra, bgra_t);
					for(i=bgrabf-GBUF;i<bgrabf;++i) {
						bgra[i].bf=GBUF;
						bgra[i].sz=0;
    					bgra[i].bg=malloc(bgra[i].bf*sizeof(bgr_t));
					}
					cc++;
					// CONDREALLOC(bgra[cc].sz, bgra[cc].bf, GBUF, bgra[cc].bg, bgr_t);

					if(couw==oldcouw) {
                		bgra[cc][wa->numl-prevlc].n=malloc(couc*sizeof(char));
                		strcpy(bgra[cc][wa->numl-prevlc].n, bufword);
					} else
                		bgra[cc][wa->numl-prevlc].c[couw-oldcouw-1]=atol(bufword);
                couc=0;
                couw++;
            }
            if(c=='#') {
                while( (c=fgetc(fp)) != '\n') ;
                continue;
            } else if(c=='\n') {
                if(wa->numl == wa->lbuf-1) {
                    wa->lbuf += WBUF;
                    wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
                    memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
                    bgra[cc]=realloc(bgra[cc], (wa->lbuf-prevlc)*sizeof(bgr_t));
                }
                wa->wpla[wa->numl] = couw-oldcouw;
				if(couw-oldcouw >4) {
					printf("Error, each row cannot exceed 4 words: revise your input file\n"); 
					/* need to release all memory too */
            		free_wseq(wa);
					exit(EXIT_FAILURE);
				}
                oldcouw=couw;
                wa->numl++;
				cc++;
                if(cc == bgra->bf->lbuf-1) {
                    bgrow=realloc(bgrow, wa->lbuf*sizeof(bgr_t));
                    memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
                }
            }
            inword=0;
        } else if(inword==0) { /* deal with first character of new word, + and - also allowed */
            if(couw == wa->wsbuf-1) {
                wa->wsbuf += GBUF;
                wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
                bgrow=realloc(bgrow, wa->wsbuf*sizeof(bgr_t));
                for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
                    wa->wln[i]=0;
            }
            couc=0;
            bwbuf=WBUF;
            bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
            bufword[couc++]=c; /* no need to check here, it's the first character */
            inword=1;
        // } else if( (c == 0x2E) | ((c >= 0x30) && (c <= 0x39)) ) {
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
    bgrow = realloc(bgrow, wa->quan*sizeof(bgr_t)); /* normalize */
    wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

    *m= wa->numl;
    int k=wa->wpla[0];
    for(i=1;i<wa->numl;++i)
        if(k != wa->wpla[i])
            printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
    *n= k; 
    free_wseq(wa);

    return bgrow;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=2) {
        printf("Error. Pls supply argument (name of text file).\n");
        exit(EXIT_FAILURE);
    }

    int i, j, m, n;
    // float *mat=processinpf(argv[1], &m, &n);
    bgr_t *bgrow=processinpf(argv[1], &m, &n);

    printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
    for(i=0;i<m;++i) {
        for(j=0;j<n;++j) {
            if(j==0)
				printf("%s ", bgrow[i].n);
			else
            	printf("%li ", bgrow[i].c[j-1]);
		}
        printf("\n"); 
    }
    for(i=0;i<m;++i)
		free(bgrow[i].n);
    free(bgrow);
    return 0;
}
