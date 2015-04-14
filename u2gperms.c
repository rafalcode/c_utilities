/* Normally called with the * wildcard from bash, (so hidden files and dirs will not be affected)
 * and mirrors the file's or dir's user permissions, to its group permissions */
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <bash wildcard>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    struct stat sb;
    uint i, up, gp, op, nwtp;

    for(i=1;i<=argc-1;++i) {

        if (stat(argv[i], &sb) == -1) {
            perror("stat");
            exit(EXIT_FAILURE);
        }

        up=(sb.st_mode & S_IRWXU);
        gp=(sb.st_mode & S_IRWXG);
        op=(sb.st_mode & S_IRWXO);

        if(up != gp) {
            nwtp= (uint)(up | (up >>3) | op);
            chmod(argv[i], nwtp);
        }
    }
    exit(EXIT_SUCCESS);
}
