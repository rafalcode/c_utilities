/* kmlrd.c taking after gpxrd1.c 
 * one challenge of kml is that all time points are listed first and then all locpoints.
 * so depending on the length of journey, there is no hard code for the start of the loc points.
 *
 * need to cheack for first occurence of
            <gx:coord>-2.93590 43.31035 221</gx:coord>
   that's it.
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "kmlrd.h"

//haversine courtesy of chatgpt.
#define EARTH_RADIUS_KM 6371.0
#define KMLFIRSTPT 25 // kml first timepoint.
#define MYLINE 246
#define PMLAT 47.8127 // Penmarch latitude
#define PMLON -4.3366 // Penmarch long
#define RVLAT 47.82311428292863 //rocher des victimes lat
#define RVLON -4.379444899501478 // lat of  above (from google maps).

double to_radians(double degree)
{
    return degree * M_PI / 180.0;
}

double haversine(double lat1, double lon1, double lat2, double lon2)
{
    double dLat = to_radians(lat2 - lat1);
    double dLon = to_radians(lon2 - lon1);

    lat1 = to_radians(lat1);
    lat2 = to_radians(lat2);

    double a = sin(dLat / 2) * sin(dLat / 2) +
        sin(dLon / 2) * sin(dLon / 2) * cos(lat1) * cos(lat2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    return EARTH_RADIUS_KM * c;
}

void parsetime(char *time, hmst_t *hmst)
{
    // ignores date part (i.e. up tot T character.
    int j;
    int hrlen, minlen, seclen;
    char *t0, *t1, *t2, *t3;
    char hr[12]={'\0'};
    char min[12]={'\0'};
    char sec[12]={'\0'};
    t0=strchr(time, 'T');
    t1=strchr(time, ':');
    t2=strchr(t1+1, ':');
    t3=strchr(t2+1, 'Z');
    hrlen=(int)(t1-t0);
    for(j=0;j<hrlen-1;j++)
        hr[j]=t0[j+1];
    hmst->h=atoi(hr);
    minlen=(int)(t2-t1);
    for(j=0;j<minlen-1;j++)
        min[j]=t1[j+1];
    hmst->m=atoi(min);
    seclen=(int)(t3-t2);
    for(j=0;j<seclen-1;j++)
        sec[j]=t2[j+1];
    hmst->s=atoi(sec);
}

w_c *crea_wc(unsigned initsz)
{
    w_c *wc=malloc(sizeof(w_c));
    wc->lp1=initsz;
    wc->t=ORDW;
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
    aawc->ftp=0;
    aawc->seenwhew=0;
    aawc->flp=0;
    aawc->seengxw=0;
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

void prtaawcdata(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    for(i=0;i<aawc->numl;++i) {
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf("l%ut", aawc->aaw[i]->aw[j]->lp1-1);
            switch(aawc->aaw[i]->aw[j]->t) {
                case GXW: printf("G "); break;
                case WHEW: printf("W "); break;
                case ORDW: printf("O "); break;
            }
        }
    }
    printf("\n"); 
    printf("L is a line, l is length of word, S is normal string, C closing punct, Z, starting cap, Y Starting cap and closing punct.\n"); 
}

void prtaawcplain(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    printf("Just line %i\n", MYLINE);
    for(i=0;i<aawc->numl;++i) {
        if(i==MYLINE) {
        printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        for(j=0;j<aawc->aaw[i]->al;++j)
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
}

void prtaawcplain2(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    boole seengxw=0;
    for(i=0;i<aawc->numl;++i) {
        for(j=0;j<aawc->aaw[i]->al;++j) {
            if((seengxw==0) & (aawc->aaw[i]->aw[j]->t==GXW)) {
                printf("GXW at line %i, word %i: %s\n", i, j, aawc->aaw[i]->aw[j]->w);
                seengxw=1;
            }
        }
    }
}

void prtaawcplain3(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i, j;
    int extent=aawc->flp - aawc->ftp;
    for(i=aawc->ftp;i<aawc->flp+extent;++i) {
        for(j=0;j<aawc->aaw[i]->al;++j) {
            printf((j!=aawc->aaw[i]->al-1)?"%s ":"%s\n", aawc->aaw[i]->aw[j]->w);
        }
    }
}

void prtaawcpla2(aaw_c *aawc) /* print line and word details, but not the words themselves */
{
    int i;
    double lon, lat;
    float ele;
    char *time;
    printf("%24s%14s%14s%14s\n", "TIME", "LON", "LAT", "ELE"); 
    for(i=9;i<aawc->numl-3;i+=4) {
        // printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        // for(j=0;j<aawc->aaw[i]->al;++j)
        lon=strtod(aawc->aaw[i]->aw[2]->w, NULL);
        lat=strtod(aawc->aaw[i]->aw[4]->w, NULL);
        ele=atof(aawc->aaw[i+1]->aw[1]->w);
        time=aawc->aaw[i+2]->aw[1]->w;
        printf("%24s%14.6f%14.6f%14.6f\n", time, lon, lat, ele);
    }
}

void prtaawcpla3(aaw_c *aawc) /* garmin connect running gpx */
{
    int i, k=0;
    double lon, lat;
    float ele;
    char *time;
    printf("%24s%14s%14s%14s\n", "TIME", "LON", "LAT", "ELE"); 
    for(i=16;i<aawc->numl-8;i+=9) {
        // printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        // for(j=0;j<aawc->aaw[i]->al;++j)
        lon=strtod(aawc->aaw[i]->aw[2]->w, NULL);
        lat=strtod(aawc->aaw[i]->aw[4]->w, NULL);
        ele=atof(aawc->aaw[i+1]->aw[1]->w);
        time=aawc->aaw[i+2]->aw[1]->w;
        printf("%i: %24s%14.6f%14.6f%14.6f\n", k++, time, lon, lat, ele);
    }
    printf("total lines=%i\n", k); 
}

void prtaawcpla300(aaw_c *aawc) /* kml from flightware */
{
    int i, k, kk=0;
    double lon, lat, lon2, lat2;
    double dlon, dlat;
    double hdist;
    int ele, dele, ele2;
    hmst_t *hmst=calloc(1, sizeof(hmst_t));
    hmst_t *hmst2=calloc(1, sizeof(hmst_t));
    char *timestr;

    size_t allsecs, allsecs2, tdiff; // all seconds
    // header:
    // printf("%24s%12s%12s%12s%12s%12s%14s%14s%14s\n", "TIME", "HR", "MIN", "SEC", "THOU", "ASECS", "LON", "LAT", "ELE"); 
    // first trkpt: absolute, i.e starting point.
    printf("First trkpt: absolute, i.e starting point:\n");
    i=aawc->ftp;
    k=aawc->flp;
    // printf("L%u(%uw):", i, aawc->aaw[i]->al); 
    // for(j=0;j<aawc->aaw[i]->al;++j)
    lon=strtod(aawc->aaw[k]->aw[1]->w, NULL);
    lat=strtod(aawc->aaw[k]->aw[2]->w, NULL);
    ele=atoi(aawc->aaw[k]->aw[3]->w);

    timestr=aawc->aaw[i]->aw[1]->w;
    parsetime(aawc->aaw[i]->aw[1]->w, hmst);
    allsecs = 3600*hmst->h + 60*hmst->m + hmst->s;
    printf("%i: %24s%12i%12i%12i%12zu%14.6f%14.6f%14i\n", k++, timestr, hmst->h, hmst->m, hmst->s, allsecs, lon, lat, ele);
    // printf("%i: %24s%12i%12i%12i%14.6f\n", k++, timestr, hmst->h, hmst->m, hmst->s, allsecs);

    // the rest shall all be differences:
    printf("From now on, cumulative differences:\n"); 
    printf("%s  %24s%12s%12s%12s\n", "line", "secselapsed", "havkm", "speed", "eldiff");
    for(i=aawc->ftp+1,k=aawc->flp+1;i<aawc->flp;++i,++k) {
        // printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        // for(j=0;j<aawc->aaw[i]->al;++j)
        lon2=strtod(aawc->aaw[k]->aw[1]->w, NULL);
        lat2=strtod(aawc->aaw[k]->aw[2]->w, NULL);
        ele2=atoi(aawc->aaw[k]->aw[3]->w);
        // dlon=lon2-lon;
        // dlat=lat2-lat;
        hdist=haversine(lat, lon, lat2, lon2);
        dele=ele2-ele;
        // time=aawc->aaw[i+2]->aw[1]->w;
        timestr=aawc->aaw[i]->aw[1]->w;
        parsetime(aawc->aaw[i]->aw[1]->w, hmst2);
        allsecs2=3600*hmst2->h + 60*hmst2->m + hmst2->s;
        tdiff=allsecs2-allsecs;

        // printf("%i: %24s%12i%12i%12i%12zu%14.6f%14.6f%14i\n", kk++, timestr, hmst2->h-hmst->h, hmst2->m -hmst->m, hmst2->s-hmst->s, allsecs2-allsecs, dlon, dlat, dele);
        // printf("%i: %24s%12i%12i%12i%12zu%14.6f%14i\n", kk++, timestr, hmst2->h-hmst->h, hmst2->m -hmst->m, hmst2->s-hmst->s, allsecs2-allsecs, hdist, dele);
        printf("%i: %24zu%14.6f%14.6f%14i\n", kk++, tdiff, hdist, 3600.*hdist/tdiff, dele);
        // printf("%i: %24s%12i%12i%12i%14.6f\n", k++, timestr, hmst2->h-hmst->h, hmst2->m -hmst->m, hmst2->s-hmst->s, allsecs2-allsecs);

        lon=lon2;
        lat=lat2;
        ele=ele2;
        hmst->h=hmst2->h;
        hmst->m=hmst2->m;
        hmst->s=hmst2->s;
        allsecs=allsecs2;
    }
    // printf("total lines=%i\n", k); 
    free(hmst);
    free(hmst2);
}

void prtaawcpla301(aaw_c *aawc) /* kml from flightware */
{
    int i, k, kk=0, totalpts;
    double lon, lat, lon2, lat2;
    double dlon, dlat;
    double hdist;
    int ele, dele, ele2;
    hmst_t *hmst=calloc(1, sizeof(hmst_t));
    hmst_t *hmst2=calloc(1, sizeof(hmst_t));
    char *timestr;

    size_t allsecs, allsecs2, tdiff; // all seconds
    // header:
    // printf("%24s%12s%12s%12s%12s%12s%14s%14s%14s\n", "TIME", "HR", "MIN", "SEC", "THOU", "ASECS", "LON", "LAT", "ELE"); 
    // first trkpt: absolute, i.e starting point.
    printf("First trkpt: absolute:\n");
    i=aawc->ftp;
    k=aawc->flp;
    totalpts=k-i;
    // printf("L%u(%uw):", i, aawc->aaw[i]->al); 
    // for(j=0;j<aawc->aaw[i]->al;++j)
    lon=strtod(aawc->aaw[k]->aw[1]->w, NULL);
    lat=strtod(aawc->aaw[k]->aw[2]->w, NULL);
    ele=atoi(aawc->aaw[k]->aw[3]->w);

    timestr=aawc->aaw[i]->aw[1]->w;
    parsetime(timestr, hmst);
    allsecs = 3600*hmst->h + 60*hmst->m + hmst->s;
    printf("%i/%i: %24s%12i%12i%12i%12zu%14.6f%14.6f%14i\n", 1+kk++, totalpts, timestr, hmst->h, hmst->m, hmst->s, allsecs, lon, lat, ele);
    // printf("%i: %24s%12i%12i%12i%14.6f\n", k++, timestr, hmst->h, hmst->m, hmst->s, allsecs);

    // the rest shall all be differences:
    printf("From now on, cumulative differences:\n"); 
    printf("%s  %24s%12s%12s%12s\n", "line", "secselapsed", "havkm", "speed", "eldiff");
    for(i=aawc->ftp+1,k=aawc->flp+1;i<aawc->flp;++i,++k) {
        // printf("L%u(%uw):", i, aawc->aaw[i]->al); 
        // for(j=0;j<aawc->aaw[i]->al;++j)
        lon2=strtod(aawc->aaw[k]->aw[1]->w, NULL);
        lat2=strtod(aawc->aaw[k]->aw[2]->w, NULL);
        ele2=atoi(aawc->aaw[k]->aw[3]->w);
        // dlon=lon2-lon;
        // dlat=lat2-lat;
        hdist=haversine(lat, lon, lat2, lon2);
        dele=ele2-ele;
        // time=aawc->aaw[i+2]->aw[1]->w;
        timestr=aawc->aaw[i]->aw[1]->w;
        parsetime(aawc->aaw[i]->aw[1]->w, hmst2);
        allsecs2=3600*hmst2->h + 60*hmst2->m + hmst2->s;
        tdiff=allsecs2-allsecs;

        // printf("%i: %24s%12i%12i%12i%12zu%14.6f%14.6f%14i\n", kk++, timestr, hmst2->h-hmst->h, hmst2->m -hmst->m, hmst2->s-hmst->s, allsecs2-allsecs, dlon, dlat, dele);
        // printf("%i: %24s%12i%12i%12i%12zu%14.6f%14i\n", kk++, timestr, hmst2->h-hmst->h, hmst2->m -hmst->m, hmst2->s-hmst->s, allsecs2-allsecs, hdist, dele);
        // printf("%i/%i: %24zu%14.6f%14.6f%14i\n", 1+kk++, totalpts, tdiff, hdist, 3600.*hdist/tdiff, dele);
        printf("%i/%i: %24s%12i%12i%12i%12zu%14.6f%14.6f%14i%24zu%14.6f%14.6f%14i\n", 1+kk++, totalpts, timestr, hmst2->h, hmst2->m, hmst2->s, allsecs, lon, lat, ele, tdiff, hdist, 3600.*hdist/tdiff, dele);
        // printf("%i: %24s%12i%12i%12i%14.6f\n", k++, timestr, hmst2->h-hmst->h, hmst2->m -hmst->m, hmst2->s-hmst->s, allsecs2-allsecs);

        lon=lon2;
        lat=lat2;
        ele=ele2;
        hmst->h=hmst2->h;
        hmst->m=hmst2->m;
        hmst->s=hmst2->s;
        allsecs=allsecs2;
    }
    // printf("total lines=%i\n", k); 
    free(hmst);
    free(hmst2);
}

aaw_c *processinpf(char *fname)
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
        if( (c== '\n') | (c == ' ') | (c == '<') | (c == '>') | (c == '"') ) {
            if( inword==1) { /* cue word-ending procedure */
                aawc->aaw[aawc->numl]->aw[couw]->w[couc++]='\0';
                aawc->aaw[aawc->numl]->aw[couw]->lp1=couc;
                SETCPTYPE(oldc, aawc->aaw[aawc->numl]->aw[couw]->t);
                if((aawc->seengxw==0) & (aawc->aaw[aawc->numl]->aw[couw]->t==GXW)) {
                    aawc->seengxw=1;
                    aawc->flp=aawc->numl;
                } else if((aawc->seenwhew==0) & (aawc->aaw[aawc->numl]->aw[couw]->t==WHEW)) {
                    aawc->seenwhew=1;
                    aawc->ftp=aawc->numl;
                }
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

int main(int argc, char *argv[])
{
    /* argument accounting */
    if(argc!=2) {
        printf("Error. Pls supply argument (name of text file).\n");
        exit(EXIT_FAILURE);
    }

   aaw_c *aawc=processinpf(argv[1]);
   prtaawcpla301(aawc); // printout original text as well as you can.
   // prtaawcplain2(aawc); // printout original text as well as you can.
   // prtaawcplain3(aawc); // printout original text as well as you can.
   printf("Numlines: %zu\n", aawc->numl); 
   printf("ftp: %zu\n", aawc->ftp); 
   printf("flp: %zu\n", aawc->flp); 

    free_aawc(&aawc);

    return 0;
}
