/* c program to mirror the user permissions of all files and directories
 * to their corresponding group permissions:
 * ie. a file with 750 permission octal code will be 770 afterwards.
 * This part of the program is reflected in the "u2g" part of its name.
 *
 * The "trav" part of its name reflects it manner of traversing (recursively) the directory tree starting from the
 * CWD 
 *
 * Coding notes: there are several aspects of thsi program worth explaining.
 * First the prgram itself travels through the directory structure, so at any one moment the code need only concern itself
 * with the current working directory.
 *
 * It seems that even files are treated as directories, so when you see d_name, it means f_name also.
 *
 * */
#include <dirent.h> 
#include <sys/types.h> 
#include <sys/param.h> 
#include <sys/stat.h> 
#include <unistd.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 

#define MYCWD "."
#define UPADIR ".."

void u2g(char *dfent)
{
    struct stat sb;
    uint up, gp, op, nwtp;

    if (stat(dfent, &sb) == -1) {
        perror("stat");
        exit(EXIT_FAILURE);
    }

    up=(sb.st_mode & S_IRWXU);
    gp=(sb.st_mode & S_IRWXG);
    op=(sb.st_mode & S_IRWXO);

    if(up != gp) {
        nwtp= (uint)(up | (up >>3) | op);
        chmod(dfent, nwtp);
    }
}

int walker(char *cd)
{
    struct dirent *f_or_dir;
    DIR *d;
    if((d = opendir(cd)) == NULL) {
        printf("can't access current directory.\n"); 
        exit(EXIT_FAILURE);
    }

    while( (f_or_dir=readdir(d)) ) {

        /* ignore the "." and ".." entities in the Cur Working Dir */
        if( strcmp( f_or_dir->d_name, MYCWD ) == 0 || strcmp( f_or_dir->d_name, UPADIR ) == 0 )
            continue;

        if( f_or_dir->d_type == DT_DIR ) {
            chdir(f_or_dir->d_name);
            walker(MYCWD);
            chdir(UPADIR);
            u2g(f_or_dir->d_name); 
        } else {
            u2g(f_or_dir->d_name); 
        }
    }
    closedir(d);
    return 0;
}

int main(void)
{
    walker(MYCWD);
    return 0;
}
