/* allfl.c "ALL Files List" .. a program to list all the files in a (heavily populated) directory and
 * store them is a file. Seem to repeat the function of tother tools, but the idea here is that we want a text file of all the files
 * in the current directory */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <sys/param.h>

#define GBUF 4
#define STRBUF 32

typedef struct /* xncou */
{
    char *xn;
    unsigned cou;
} xncou;

static int cmpxncou(const void *p1, const void *p2)
{
    int ret;
    xncou *x1=(xncou*)p1;
    xncou *x2=(xncou*)p2;
    if(x1->cou < x2->cou)
        ret=1;
    else if(x1->cou == x2->cou)
        ret = 0;
    else
        ret=-1;
    return ret;
}

int main(int argc, char *argv[])
{
    DIR *dirp;
    struct dirent *direntp;
    int i, j;

    if(argc != 2) {
        fprintf( stderr, "Usage: %s <directory-to-analyse>\n", argv[0]);
        exit(1);
    }

    if ((dirp = opendir(argv[1])) == NULL) {
        perror("directory name problem");
        exit(2);
    }

    unsigned xsbuf=GBUF;
    xncou *xsa=malloc(xsbuf*sizeof(xncou));
    for(i=0;i<xsbuf;++i) {
        xsa[i].xn=calloc(STRBUF, sizeof(char));
        xsa[i].cou=0;
    }
    unsigned xcou=0;
    unsigned noextcou=0;
    unsigned nonreg=0;
    unsigned char seenx=0;

    char *ppoi=NULL;
    /* d_types:DT_BLK blockdev; DT_CHR chardev; DT_DIR straight dir; DT_FIFO namedpipe; DT_LNK symlink; DT_REG regular file; DT_SOCK UNIX domain socket; DT_UNKNOWN file type unknown */
    while ((direntp = readdir(dirp)) != NULL) {
        if(direntp->d_type == DT_REG) {
            if((ppoi=strrchr(direntp->d_name, '.')) == NULL)
                noextcou++;
            else {
                for(i=0;i<xcou;++i)
                    if(!strcmp(xsa[i].xn, ppoi)) {
                        seenx=1;
                        xsa[i].cou++;
                        break;
                    }
                if(!seenx) {
                    if(xcou==xsbuf-1) {
                        xsbuf += GBUF;
                        xsa=realloc(xsa, xsbuf*sizeof(xncou));
                        for(j=xsbuf-GBUF;j<xsbuf;++j) {
                            xsa[j].xn=calloc(STRBUF, sizeof(char));
                            xsa[j].cou=0;
                        }
                    }
                    strcpy(xsa[xcou].xn, ppoi);
                    xsa[xcou].cou++;
                    xcou++;
                    seenx=0;
                }
            }
        } else
            nonreg++;
    }
    /*normalize */
    for(j=xcou;j<xsbuf;++j)
        free(xsa[j].xn);
    xsa=realloc(xsa, xcou*sizeof(xncou));
    qsort(xsa, xcou, sizeof(xncou), cmpxncou);

    printf("In the directory \"%s\", there were regular files with %u different extensions, Extension names and counts now follow in descending count order:\n", argv[1], xcou);
    for(i=0;i<xcou;++i)
        printf("%s:%u ", xsa[i].xn, xsa[i].cou);
    printf("\n"); 
    printf("Furthermore, %u files had no extension, and there were %u non-regular files (eg. subdirectories (note: nested included), etc).\n", noextcou, nonreg);

    /* free */
    for(j=0;j<xcou;++j)
        free(xsa[j].xn);
    free(xsa);
    closedir(dirp);
    return 0;
}
