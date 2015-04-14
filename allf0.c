#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <sys/param.h>

int main(int argc, char *argv[])
{
    DIR *dirp;
    struct dirent *direntp;

    if(argc != 3) {
        fprintf( stderr, "Usage: %s dir_name\n", argv[0]);
        exit(1);
    }

    if ((dirp = opendir(argv[1])) == NULL) {
        perror(argv[1]);
        exit(2);
    }

    unsigned regcou=0, cou=0;
    char *ppoi=NULL;
    short extsz=strlen(argv[2]);
    /* d_types:DT_BLK blockdev; DT_CHR chardev; DT_DIR straight dir; DT_FIFO namedpipe; DT_LNK symlink; DT_REG regular file; DT_SOCK UNIX domain socket; DT_UNKNOWN file type unknown */
    while ((direntp = readdir(dirp)) != NULL)
        if(direntp->d_type == DT_REG) {
            regcou++;
            if( (ppoi=strrchr(direntp->d_name, '.')) != NULL)
                if(strncmp((ppoi+1), argv[2], extsz)==0)
                    cou++;
        }

    printf("In the directory \"%s\", of the total %u regular files, %u had the extension \"%s\".\n", argv[1], regcou, cou, argv[2]);

    closedir(dirp);
    return 0;
}
