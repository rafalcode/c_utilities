/* Normally called with the * wildcard from bash, (so hidden files and dirs will not be affected)
 * and mirrors the file's or dir's user permissions, to its group permissions */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <dirent.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h> 
#include <pwd.h>
#include <grp.h>

#define NUMPROJN 6
#define GSTRBF 4096 /* generous string buffer */
#define SSTRBF 64 /* small string buffer */
#define MYCWD "."
#define UPADIR ".."


// const char *pn[NUMPROJN]={"b2010059", "b2010060", "b2012209", "b2013182", "b2014050", "b2014104"};
const char *pn[NUMPROJN]={"b2010059", "b2010060", "b2012209", "nutria", "b2014050", "b2014104"};

void usage(char *argname)
{
    fprintf(stderr, "Program \"%s\" takes a single directoryname as input, enters it,\n", argname);
    fprintf(stderr, "and sanitises groups ownerships on all files and subdirectories within.\n");
    fprintf(stderr, "Usage: %s dirname\n", argname);
}

const char *getprojn(char *cwd) /* get the project name from the pwd */
{
    int i;
    char *tp=strtok(cwd, "/");
    unsigned depthcount=1;
    for(i=0;i<NUMPROJN;++i) 
        if(!strcmp(tp, pn[i])) {
            // printf("%s dir found at depth %d\n", tp, depthcount);
            return pn[i];
        }
    while((tp=strtok(NULL, "/")) != NULL ) {
        depthcount++;
        for(i=0;i<NUMPROJN;++i)  
            if(!strcmp(tp, pn[i])) 
                // printf("%s dir found at depth %d\n", tp, depthcount);
                return pn[i];
    }
    return NULL;
}

void sanitisegrp(char *path, uid_t relui, gid_t relgi)
{
    struct stat sb;

    if (stat(path, &sb) == -1) {
        perror(path);
        exit(EXIT_FAILURE);
    }

    /* Ensure that current user is also owner of the dir_or_file */
    if(sb.st_uid !=  relui)
        printf("Dir/file %s cannot be chgrp'd nor setgid'd\n", path);

    /* Ensure that current group ownership is not the relevant project's */
    if(sb.st_gid != relgi)
        if( chown(path, relui, relgi))  { /* throw error if chown was not posible */
            printf("There's a problem chown'ing %s. Skipped.\n", path);
            return;
        }

    if(S_ISDIR(sb.st_mode)) { /* we are going to forcibly set read, execution and setgid bits */
        sb.st_mode |= S_IRGRP;
        sb.st_mode |= S_IXGRP;
        sb.st_mode |= S_ISGID;
        chmod(path, sb.st_mode);
    }
}

int walker(char *cd, uid_t relui, gid_t relgi)
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
            if(chdir(f_or_dir->d_name)) {
                printf("Couldn't access directory \"%s\"\n", f_or_dir->d_name); 
                exit(EXIT_FAILURE);
            }
            walker(MYCWD, relui, relgi);
            if(chdir(UPADIR)) {
                printf("Couldn't come up back from a directory\n");
                exit(EXIT_FAILURE);
            }
            sanitisegrp(f_or_dir->d_name, relui, relgi); 
        } else {
            sanitisegrp(f_or_dir->d_name, relui, relgi); 
        }
    }
    closedir(d);
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc != 2) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    char *launchdir=get_current_dir_name();

    DIR *targetdir;
    if((targetdir = opendir(argv[1])) == NULL) {
        printf("can't access the \"%s\" directory. Try absolute path.\n", argv[1]); 
        exit(EXIT_FAILURE);
    }
    closedir(targetdir);

    /* OK. our target directory is accessible at least, we can now work out intended users and groups */
    char *cwd=get_current_dir_name();
#ifdef DBG
    printf("found %s\n", getprojn(cwd));
#endif
    const char *relgrp;
    if((relgrp=getprojn(cwd)) == NULL) { /* beware you can only call this function once! */
        printf("Relevant project group not found in path, or getprojn() function called more than once\n"); 
        exit (EXIT_FAILURE);
    }
    struct group *relg;
    if( (relg=getgrnam(relgrp)) == NULL) {
        printf("relgrp %s is not a valid group in system.\n", relgrp);
        exit(EXIT_FAILURE);
    }
    gid_t relgi=relg->gr_gid; /* gets unsigned version of relgrp, which we want directory to be set to */

    char una[SSTRBF]={0};
    getlogin_r(una, SSTRBF);
#ifdef DBG
    printf("Name is %s\n", una);
#endif

    struct passwd *relu=getpwnam(una);
    if( (relu=getpwnam(una)) == NULL) {
        printf("user %s is not a valid user in system.\n", una);
        exit(EXIT_FAILURE);
    }
    uid_t relui=relu->pw_uid; /*we'll only operate on user-set files anyway, put chown needs this */

    /* OK, so we have relui and relgi, before walking within the target directory, we want to change groups on the target directory itself */
    sanitisegrp(argv[1], relui, relgi); 

    /* the following is almost a second check, but anyhow we now want to chdir into the target directory to begin walking */
    int retchdir=chdir(argv[1]);
#ifdef DBG
    printf("%d\n", retchdir); 
#endif
    if(retchdir) {
        printf("Error. Check permissions on \"%s\"\n", argv[1]); 
        exit(EXIT_FAILURE);
    }

    walker(MYCWD, relui, relgi);
    if(chdir(launchdir)) {
        printf("Couldn't get back to original directory \"%s\"\n", launchdir); 
        exit(EXIT_FAILURE);
    }

    free(cwd);
    free(launchdir);
    return 0;
}
