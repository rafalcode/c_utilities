/* bgmergmcstealth: this takes a bed file, and a minimum coverage and merges all the segments that 1) belong to same contig 2) are contiguous 3) satisfy mincov
 * This program functions correctly but is not memory optimized.. with more time, this can be improved, but should only matter for very big bedgraphs */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#ifdef DBG
#define GBUF 2
#define WBUF 2
#else
#define GBUF 32
#define WBUF 32
#endif

#define CONDREALLOC(x, b, c, a, t); \
	if((x)>=((b)-1)) { \
		(b) += (c); \
		(a)=realloc((a), (b)*sizeof(t)); \
		memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
	}

typedef unsigned char boole;

typedef struct /* i4_t */
{
	int sc; /* number of same chromosomes */
	int mc; /* min coverage associated */
	int b1i; /* index of the 1st bgr_t, which satisfies the conditions */
	int lgbi; /* last good bgr_t index */
} i4_t; /* bedgraph row type */

typedef struct /* bgr_t */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) end */
	int co; /* coverage */
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
				} else if((couw-oldcouw)<4) /* it's not the first word, and it's 1st and second col */
					bgrow[wa->numl].c[couw-oldcouw-1]=atol(bufword);
				else if( ((couw-oldcouw)==4) & (((int)(bufword[0])>47) & ((int)(bufword[0])<58)) ) /* 4th col could be a string, in which case, don't store */
					bgrow[wa->numl].co=atoi(bufword);
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

void prtmbed(bgr_t **bgra, i4_t *dca, int dcasz, int n) /* the 2D version */
{
	int i, j;
	for(i=0;i<dcasz;++i) {
		for(j=0;j<dca[i].sc;++j) { // we're cycling though all of them, though we're really only interested in the first and last.
			if(j==0) { 
				printf("%s ", bgra[i][j].n);
				printf("%li ", bgra[i][j].c[0]);
			}
			if(j==dca[i].sc-1) { // note this cannot be an else if, because if only one line j = 0 = dca[i]-1.
				printf("%li ", bgra[i][j].c[1]);
				printf("%i ", dca[i].mc);
			}
		}
		printf("\n"); 
	}
	return;
}

i4_t *difca(bgr_t *bgrow, int m, int *dcasz, int minco) /* An temmpt to merge bgraph quickly, no hope */
{
	int i, goodi=0 /* the last i at which minco was satisfied */;
	boole seenminco=0;
	/* how many different chromosomes are there? the dc (different chromsosome array */
	int dcbf=GBUF, dci=0;
	i4_t *dca=calloc(dcbf, sizeof(i4_t));
	char *tstr=NULL;
	/* find first bgrow element which is over the minimum coverage */
	for(i=0;i<m;++i)
		if(bgrow[i].co >= minco) {
			tstr=malloc(bgrow[i].nsz*sizeof(char)); /* tmp string */
			strcpy(tstr, bgrow[i].n);
			dca[dci].sc++;
			dca[dci].mc=bgrow[i].co;
			dca[dci].b1i=i;
			dca[dci].lgbi=i;
			seenminco=1;
			goodi=i;
			break;
		}
	if(!seenminco) {
		printf("Error. No bedgraph element was able to satisfy the minimum coverage value that was specified: abandoning ship.\n");
		exit(EXIT_FAILURE);
	}

	for(i=goodi+1;i<m;++i) {
		/* the same now means same name and contiguous */
		if( (!strcmp(tstr, bgrow[i].n)) & (bgrow[i].c[0] - /*stealth*/(bgrow[dca[dci].lgbi].c[1])<9999) & (bgrow[i].co >= minco) ) {
			dca[dci].sc++;
			dca[dci].lgbi=i;
			if(bgrow[i].co<dca[dci].mc)
				dca[dci].mc=bgrow[i].co;
		} else if (bgrow[i].co >= minco) {
			CONDREALLOC(dci, dcbf, GBUF, dca, i4_t);
			dci++;
			dca[dci].sc++;
			dca[dci].mc=bgrow[i].co;
			dca[dci].b1i=i;
			dca[dci].lgbi=i;
			/* new string could be differnt length*/
			tstr=realloc(tstr, bgrow[i].nsz*sizeof(char)); /* tmp string */
			strcpy(tstr, bgrow[i].n);
		}
	}
	dca=realloc(dca, (dci+1)*sizeof(i4_t));
#ifdef DBG
	printf("Num of different chromcontigs=%i. How many of each? Let's see:\n", dci+1); 
	printf("dcbf=%i. 4-tupe is sc/mc/b1i/lgbi\n", dcbf); 
	for(i=0;i<=dci;++i) 
		printf("%i/%i/%i/%i ",dca[i].sc, dca[i].mc, dca[i].b1i, dca[i].lgbi); 
	printf("\n"); 
#endif
	*dcasz=dci+1;
	free(tstr);
	return dca;
}

int main(int argc, char *argv[])
{
	/* argument accounting */
	if(argc!=3) {
		printf("bgmergmc: this takes a bed file, and a minimum coverage value and merges all the segments that:\n");
		printf("1) belong to same contig 2) are contiguous 3) satisfy a min coverage. Output is bedgraph to STDOUT.\n");
		printf("Usage. Pls supply 2 arguments 1) name of bedgraph file 2) minimum coverage.\n");
		exit(EXIT_FAILURE);
	}
	int minco=atoi(argv[2]);
	int i, j, m, n;

	bgr_t *bgrow=processinpf(argv[1], &m, &n);
	// prtbed(bgrow, m, n);
	int dcasz;
	i4_t *dca=difca(bgrow, m, &dcasz, minco);

	bgr_t **bgra=malloc(dcasz*sizeof(bgr_t*));
	/* splits */
	bgra[0]=malloc(dca[0].sc*sizeof(bgr_t)); /* first one is special */
	for(j=0;j<dca[0].sc;++j) {
		bgra[0][j].n=malloc(bgrow[dca[0].b1i+j].nsz*sizeof(char));
		strcpy(bgra[0][j].n, bgrow[dca[0].b1i+j].n);
		memcpy(bgra[0][j].c, bgrow[dca[0].b1i+j].c, 2*sizeof(long));
		bgra[0][j].c[2] = dca[0].mc;
	}
	for(i=1;i<dcasz;++i) {
		bgra[i]=malloc(dca[i].sc*sizeof(bgr_t));
		for(j=0;j<dca[i].sc;++j) {
			bgra[i][j].n=malloc(bgrow[dca[i].b1i+j].nsz*sizeof(char));
			strcpy(bgra[i][j].n, bgrow[dca[i].b1i+j].n);
			memcpy(bgra[i][j].c, bgrow[dca[i].b1i+j].c, 2*sizeof(long));
			bgra[i][j].c[2] = dca[i].mc;
		}
	}

	/* now parsed, we can get rid of bgrow now, now that we have bgra */
	for(i=0;i<m;++i)
		free(bgrow[i].n);
	free(bgrow);

	prtmbed(bgra, dca, dcasz, n);

	/* free bgra */
	for(i=0;i<dcasz;++i) {
		for(j=0;j<dca[i].sc;++j)
			free(bgra[i][j].n);
		free(bgra[i]);
	}
	free(bgra);

	free(dca);
	return 0;
}
