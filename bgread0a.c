/* Hopelessly broken, please fix condreallloc:  another failed attempted, jst trying to get unqiue strings. */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define GBUF 2
#define WBUF 2

/* ugly hacking goign on here witht he char* */
#define CONDREALLOC(x, b, c, a, d, t, i); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        (d)=realloc((d), (b)*sizeof(char*)); \
        for((i)=(((a)+(b)-(c)); (i)<(((a)+(b)); ((i)++)) { \
                    d[(i)]=NULL;
                    } \
		memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
    }

typedef unsigned char boole;

typedef struct /* bgr_t */
{
	char *n;
	size_t nsz;
	long c[3]; /* coords: 1) start 2) end 3) coverage */
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

    bgr_t *bgrow=malloc(GBUF*sizeof(bgr_t));

    while( (c=fgetc(fp)) != EOF) { /* grab a char */
        if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
            if( inword==1) { /* first word closing event */
                wa->wln[couw]=couc;
                bufword[couc++]='\0';
                bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
                	bgrow[wa->numl].n=malloc(couc*sizeof(char));
                	bgrow[wa->numl].nsz=couc;
                	strcpy(bgrow[wa->numl].n, bufword);
				} else /* it's not the first word */
                	bgrow[wa->numl].c[couw-oldcouw-1]=atol(bufword);
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
                    bgrow=realloc(bgrow, wa->lbuf*sizeof(bgr_t));
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
                bgrow=realloc(bgrow, wa->wsbuf*sizeof(bgr_t));
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

void prtbed(bgr_t *bgrow, int m, int n)
{
	int i, j;
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
	return;
}

void prtbed2(bgr_t **bgra, int *dca, int dcasz, int n) /* the 2D version */
{
	int i, j, k;
	for(i=0;i<dcasz;++i) {
		for(j=0;j<dca[i];++j) {
			for(k=0;k<n;++k) {
            	if(k==0)
					printf("%s ", bgra[i][j].n);
				else
            		printf("%li ", bgra[i][j].c[k-1]);
			}
        	printf("\n"); 
		}
	}
	return;
}

int *difca(bgr_t *bgrow, int m, int *dcasz, char **rdcca) /* find out how many differnt chromosomes there are */
{
	int i;
	/* how many different chromosomes are there? the dc (different chromsosome array */
	int dcbf=GBUF, dci=0;
	int *dca=calloc(dcbf, sizeof(int));
	char **dcca=malloc(dcbf*sizeof(char*)); /* want a list of the names */
	size_t sz1=strlen(bgrow[0].n);
	char *tstr=malloc((1+sz1)*sizeof(char)); /* tmp string */
	/* deal with first outside of loop */
	strcpy(tstr, bgrow[dci].n);
	dca[dci]++;
	dcca[0]=malloc((1+sz1)*sizeof(char)); /* tmp string */
	strcpy(dcca[0], bgrow[dci].n);
    for(i=1;i<m;++i) {
		if(!strcmp(tstr, bgrow[i].n)) /* if they are the same, increment the current count */
			dca[dci]++;
		else {
			CONDREALLOC(dci, dcbf, GBUF, dca, dcca, int);
			dci++;
			dca[dci]++;
			/* new string could be differnt length*/
			sz1=strlen(bgrow[i].n);
			tstr=realloc(tstr, bgrow[i].nsz*sizeof(char)); /* tmp string */
			strcpy(tstr, bgrow[i].n);
	        dcca[dci]=malloc((1+sz1)*sizeof(char)); /* tmp string */
	        strcpy(dcca[dci], bgrow[dci].n);
		}
	}
	dca=realloc(dca, (dci+1)*sizeof(int));
 	printf("Num of different chroms=%i. How many each of the same name?\n", dci+1); 
	printf("dcbf=%i\n", dcbf); 
	for(i=0;i<=dci;++i) 
		printf("%i ",dca[i]); 
	printf("\n"); 
	*dcasz=dci+1;
	free(tstr);
    *rdcca=*dcca;
	return dca;
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
	// prtbed(bgrow, m, n);
	int dcasz, cumsz;
    char **dcca=NULL;
	int *dca=difca(bgrow, m, &dcasz, dcca);

	bgr_t **bgra=malloc(dcasz*sizeof(bgr_t*));
	/* splits */
	bgra[0]=malloc(dca[0]*sizeof(bgr_t)); /* first one is special */
	for(j=0;j<dca[0];++j) {
		bgra[0][j].n=malloc(bgrow[j].nsz*sizeof(char));
		strcpy(bgra[0][j].n, bgrow[j].n);
		memcpy(bgra[0][j].c, bgrow[j].c, 3*sizeof(long));
	}
	cumsz=dca[0];
	for(i=1;i<dcasz;++i) {
		bgra[i]=malloc(dca[i]*sizeof(bgr_t));
		for(j=0;j<dca[i];++j) {
			bgra[i][j].n=malloc(bgrow[cumsz+j].nsz*sizeof(char));
			strcpy(bgra[i][j].n, bgrow[cumsz+j].n);
			memcpy(bgra[i][j].c, bgrow[cumsz+j].c, 3*sizeof(long));
		}
		cumsz += dca[i];
	}

	/* now parsed, we can get rid of bgrow now, now that we have bgra */
    for(i=0;i<m;++i)
		free(bgrow[i].n);
    free(bgrow);

	prtbed2(bgra, dca, dcasz, n);

	/* free bgra */
	for(i=0;i<dcasz;++i) {
		for(j=0;j<dca[i];++j)
			free(bgra[i][j].n);
		free(bgra[i]);
	}
	free(bgra);

	free(dca);
	for(i=0;i<dcasz;++i)
        free(dcca[i]);
    free(dcca);
    return 0;
}
