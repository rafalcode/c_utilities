/* macsigf: this takes a mac signal bed file, whose fourth column has a float value and filters out those under a certain signal value */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h> // required for optopt, opterr and optarg.

#ifdef DBG
#define GBUF 2
#define WBUF 2
#else
#define GBUF 32
#define WBUF 32
#endif
#define NUMBUCKETS 10

#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
    }

typedef unsigned char boole;

typedef struct  /* opt_t, a struct for the options */
{
    boole dflg; /* details / information only */
    boole pflg; /* simply print out the bed file */
    char *istr; /* input filename */
    char *fstr; /* the floating point value to low filter off */
    char *hstr; /* the floating point value to high filter off */
} opt_t;

typedef struct /* i4_t */
{
    int sc; /* number of same chromosomes */
    float mc; /* min signal value */
    int b1i; /* index of the 1st bgr_t, which satisfies the conditions */
    int lgbi; /* last good bgr_t index */
} i4_t; /* bedgraph row type */

typedef struct /* bgr_t */
{
    char *n;
    size_t nsz; /* size of the name r ID field */
    long c[2]; /* coords: 1) start 2) end */
    float co; /* signal value */
    char *nf; /* no float */
    char **rc; /* rest of the columns as strings */
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

bgr_t *processinpf(char *fname, int *m, int *n, boole *isnf /* this relates to column 4 ... is it a float or not */)
{
    /* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
     * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
     * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

    /* declarations */
    FILE *fp=fopen(fname,"r");
    int i, nrc=0 /* number of rest of columns */ ;
    size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
    int c;
    boole inword=0;
    (*isnf)=0; // is no float? i.e. column 4 i snot a float ... initialise as zero */
    wseq_t *wa=create_wseq_t(GBUF);
    size_t bwbuf=WBUF;
    char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

    bgr_t *bgrow=malloc(GBUF*sizeof(bgr_t));
    for(i=0;i<GBUF;++i) 
        bgrow[i].rc=NULL;

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
                } else if((couw-oldcouw)<3) { /* it's not the first word, and it's 1st and second col */
                    bgrow[wa->numl].c[couw-oldcouw-1]=atol(bufword);
                } else if( (couw-oldcouw)==3) { // assume float
                    for(i=0;i<couc-1;++i) { // note above how final char is a \0: don't check it!
                        if( (!*isnf) & ((bufword[i] <48) | (bufword[i] > 57)) & (bufword[i] != 46) ) {
                            *isnf=1;
                        }
                    }
                    if(!(*isnf))
                        bgrow[wa->numl].co=atof(bufword);
                    else {
                        bgrow[wa->numl].nf=malloc(couc*sizeof(char));
                        strcpy(bgrow[wa->numl].nf, bufword);
                    }
                } else {
                    nrc++;
                    bgrow[wa->numl].rc=realloc(bgrow[wa->numl].rc, nrc*sizeof(char*));
                    bgrow[wa->numl].rc[nrc-1]=malloc(couc*sizeof(char));
                    strcpy(bgrow[wa->numl].rc[nrc-1], bufword);
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
                    bgrow=realloc(bgrow, wa->lbuf*sizeof(bgr_t));
                    for(i=wa->lbuf-WBUF;i<wa->lbuf;++i) 
                        bgrow[i].rc=NULL;
                    memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
                }
                wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */

                oldcouw=couw; /* restart words per line count */
                nrc=0;
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

void prtobed(bgr_t *bgrow, int m, int n, float minsig) // print over bed ... a value that is over a certain signal
{
    int i, j;
    for(i=0;i<m;++i) {
        if(bgrow[i].co >= minsig) {
            for(j=0;j<n;++j) {
                if(j==0)
                    printf("%s ", bgrow[i].n);
                else if(j==3)
                    printf("%2.6f ", bgrow[i].co);
                else
                    printf("%li ", bgrow[i].c[j-1]);
            }
            printf("\n"); 
        }
    }
    return;
}

void prtobed2(bgr_t *bgrow, int m, int n, float minsig, float maxsig) // print over bed ... a value that is over a certain signal
{
    int i, j;
    for(i=0;i<m;++i) {
        if( (bgrow[i].co >= minsig) & (bgrow[i].co < maxsig) ) {
            for(j=0;j<n;++j) {
                if(j==0)
                    printf("%s ", bgrow[i].n);
                else if(j==3)
                    printf("%2.6f ", bgrow[i].co);
                else
                    printf("%li ", bgrow[i].c[j-1]);
            }
            printf("\n"); 
        }
    }
    return;
}

int *hist_co(bgr_t *bgrow, int m, float mxco, float mnco, int numbuckets)
{
    int i, j;
    float step=(mxco-mnco)/(float)numbuckets;
    float *bucketlimarr=malloc((numbuckets-1)*sizeof(float));
    int *bucketarr=calloc(numbuckets, sizeof(int));
    bucketlimarr[0]=step+mnco;
    for(i=1;i<numbuckets-1;++i) 
        bucketlimarr[i]=bucketlimarr[i-1]+step;

    for(i=0;i<m;++i)
        if(bgrow[i].co>=bucketlimarr[numbuckets-2]) {
            bucketarr[numbuckets-1]++;
            continue;
        } else {
            for(j=0;j<numbuckets-1;++j)
                if(bgrow[i].co < bucketlimarr[j]) {
                    bucketarr[j]++;
                    break;
                }
        }
    free(bucketlimarr);
    return bucketarr;
}

void prthist(char *histname, int *bucketarr, int numbuckets, int m, float mxco, float mnco)
{
    int i;
    // printf("HISTOGRAM: \"%s\" NUMBINS: %d for: %-24.24s (totels=%04i):\n", histname, numbuckets, histname, m); 
    printf("HISTOGRAM: \"%s\" NUMBINS: %d RANGEEACHBIN: %4.4f TOTALELEMENTS: %04i:\n", histname, numbuckets, (mxco-mnco)/numbuckets, m); 
    printf("minval=%4.4f |", mnco);
    for(i=0;i<8*12;++i) 
        putchar(' ');
    printf("| maxval=%4.4f\n", mxco); 
    for(i=0;i<14;++i) 
        putchar(' ');
    for(i=0;i<numbuckets;++i) 
        printf("| %6i ", bucketarr[i]);
    // printf("|->maxval=%4.6f\n", mxco); 
    printf("|\n"); 
    return;
}

void prtmacsig(bgr_t *bgrow, int m, int n)
{
    int i;
    float mxco=.0, mnco=10e20;
    printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
    for(i=0;i<m;++i) {
        if(bgrow[i].co > mxco)
            mxco=bgrow[i].co;
        if(bgrow[i].co < mnco)
            mnco = bgrow[i].co;
    }
    int *hco=hist_co(bgrow, m, mxco, mnco, NUMBUCKETS);
    prthist("Macs2 Signal", hco, NUMBUCKETS, m, mxco, mnco);
    free(hco);
    return;
}

void prtbed(bgr_t *bgrow, int m, int n, boole isnf)
{
    int i, j;
    printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
    for(i=0;i<m;++i) {
        for(j=0;j<n;++j) {
            if(j==0)
                printf("%s\t", bgrow[i].n);
            else if(j==3) {
                if(!isnf)
                    printf("%2.6f ", bgrow[i].co);
                else 
                    printf("%s\t", bgrow[i].nf);
            } else if ( (j==1) | (j==2) ) 
                printf("%li ", bgrow[i].c[j-1]);
            else 
                printf("%s\t", bgrow[i].rc[j-4]);
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
                printf("%2.6f ", dca[i].mc);
            }
        }
        printf("\n"); 
    }
    return;
}

i4_t *difca(bgr_t *bgrow, int m, int *dcasz, float minsig) /* An temmpt to merge bgraph quickly, no hope */
{
    int i, goodi=0 /* the last i at which minsig was satisfied */;
    boole seenminsig=0;
    /* how many different chromosomes are there? the dc (different chromsosome array */
    int dcbf=GBUF, dci=0;
    i4_t *dca=calloc(dcbf, sizeof(i4_t));
    char *tstr=NULL;
    /* find first bgrow element which is over the minimum coverage */
    for(i=0;i<m;++i)
        if(bgrow[i].co >= minsig) {
            tstr=malloc(bgrow[i].nsz*sizeof(char)); /* tmp string */
            strcpy(tstr, bgrow[i].n);
            dca[dci].sc++;
            dca[dci].mc=bgrow[i].co;
            dca[dci].b1i=i;
            dca[dci].lgbi=i;
            seenminsig=1;
            goodi=i;
            break;
        }
    if(!seenminsig) {
        printf("Error. No bedgraph element was able to satisfy the minimum signal value that was specified: abandoning ship.\n");
        exit(EXIT_FAILURE);
    }

    for(i=goodi+1;i<m;++i) {
        /* the same now means same name and contiguous */
        if( (!strcmp(tstr, bgrow[i].n)) & (bgrow[i].c[0] == bgrow[dca[dci].lgbi].c[1]) & (bgrow[i].co >= minsig) ) {
            dca[dci].sc++;
            dca[dci].lgbi=i;
            if(bgrow[i].co<dca[dci].mc)
                dca[dci].mc=bgrow[i].co;
        } else if (bgrow[i].co >= minsig) {
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

void prtusage()
{
    printf("bedsumzr:\n");
    printf("--------\n");
    printf("This is a bed file summarizer: you must specify the input file with -i, bed or bedgraph accepted.\n");
    printf("If you also specify the -d (details) option, general information on the bed file will be given.\n");
    printf("(this also serves as a useful first pass, to see what can be done with the bed file).\n");
    printf("If the 4th column is an intensity signal, you may specify a filter with -f, below which lines will be omitted.\n");
    printf("A second filter specified by -h (higher) is available: only values below it will be printed.\n");
    return;
}

int catchopts(opt_t *opts, int oargc, char **oargv)
{
    int i, c;
    opterr = 0;

    while ((c = getopt (oargc, oargv, "dpi:f:h:")) != -1)
        switch (c) {
            case 'd':
                opts->dflg = 1;
                break;
            case 'p':
                opts->pflg = 1;
                break;
            case 'i':
                opts->istr = optarg;
                break;
            case 'f':
                opts->fstr = optarg;
                break;
            case 'h':
                opts->hstr = optarg;
                break;
            case '?':
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                return 1;
            default:
                fprintf (stderr, "Wrong arguments. Please launch without arguments to see help file.\n");
                exit(EXIT_FAILURE);
        }
    boole nonopt=0;
    for(i=optind;i<oargc;++i) {
        printf("Error: \"%s\"?: Arguments only accepted if specified with \"-<letter>\" option argument.\n", oargv[i]);
        nonopt=1;
    }   
    if(nonopt) {
        printf("\nArgument errors: usage instructions are as follows:\n\n"); 
        prtusage();
    }   

    return 0;
}

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc == 1) {
        prtusage();
        exit(EXIT_FAILURE);
    }
    int i, j, m, n;
    float minsig, maxsig;
    opt_t opts={0};
    catchopts(&opts, argc, argv);
    boole filterval=0; // if we have a filter ... convenience boolean
    if(opts.fstr) {
        minsig=atof(opts.fstr);
        filterval=1;
        if(opts.hstr) {
            filterval=2;
            maxsig=atof(opts.hstr);
        }
    }

    boole isnf;
    bgr_t *bgrow=processinpf(opts.istr, &m, &n, &isnf);
    if((opts.dflg) & (n>3) ) { // bed must have a 4th column with (often) mac signal values 
        prtmacsig(bgrow, m, n);
        goto final;
    }

    if( (!isnf) & (filterval==1) & (n>3) ) { // bed must have a 4th column with (often) mac signal values 
        prtobed(bgrow, m, n, minsig);
        goto final;
    } else if( (!isnf) & (filterval==2) & (n>3) ) { // bed must have a 4th column with (often) mac signal values 
        prtobed2(bgrow, m, n, minsig, maxsig);
        goto final;
    }

    if(opts.pflg)
        prtbed(bgrow, m, n, isnf);

final: for(i=0;i<m;++i) {
           if(n>4) {
               for(j=0;j<n-4;++j) {
                   free(bgrow[i].rc[j]);
               }
               free(bgrow[i].rc);
           }
           if(isnf)
               free(bgrow[i].nf);
           free(bgrow[i].n);
       }
       free(bgrow);

       return 0;
}
