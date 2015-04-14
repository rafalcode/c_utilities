/* This used to be called filtyex ... base don a memonic about filtering by extension
 * It's not likely to be useful normally, and covers the edge case of a directory
 * which has a huge number of files of a certain extension type, and allows a command
 * to be executed on them. Because of this, I've renamed it cmd2ext, i.e. command to extension.
 *
 * Beware that it requires a parameter file with three lines:
 * .
 * c
 * ls -l
 * this will execute a ls -l on all files with the .c extension.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <sys/param.h>

#define SBZ 256 /* standard char buffer size */

int main(int argc, char *argv[])
{
    DIR *dirp;
    struct dirent *direntp;

    if(argc != 2) {
        fprintf( stderr, "Usage: %s paramfilename\n", argv[0]);
        exit(1);
    }
    char prms[3][SBZ]={{0}, {0}, {0}};
    FILE *fptr=fopen(argv[1], "r");
    fscanf(fptr, "%s\n%s\n%[^\t\n]", prms[0], prms[1], prms[2]);
    fclose(fptr);

    if ((dirp = opendir(prms[0])) == NULL) {
        perror(prms[0]);
        exit(2);
    }

    char *ppoi=NULL;
    short extsz=strlen(prms[1]);
    char cmdstrng[SBZ*2]={0};

    while ((direntp = readdir(dirp)) != NULL)
        if(direntp->d_type == DT_REG) {
            if( (ppoi=strrchr(direntp->d_name, '.')) != NULL)
                if(strncmp((ppoi+1), prms[1], extsz)==0) {
                    sprintf(cmdstrng, "%s %s/%s", prms[2], prms[0], direntp->d_name);
                    system(cmdstrng);
                }
        }

    closedir(dirp);
    return 0;
}
