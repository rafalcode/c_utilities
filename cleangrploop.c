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

void sanitisegrp(int argc, char *argv[])
{
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

    int i;
    struct stat sb;
    struct passwd *pw;
    struct group *gr;
    uint sgu;
    for(i=1;i<=argc-1;++i) {
        if (stat(argv[i], &sb) == -1) {
            perror(argv[i]);
            exit(EXIT_FAILURE);
        }
        pw = getpwuid(sb.st_uid);
        if(strcmp(pw->pw_name, una)) {
            printf("Dir/file %s cannot be chgrp'd nor setgid'd\n", argv[i]);
            continue;
        }
        gr = getgrgid(sb.st_gid);
        if(strcmp(gr->gr_name, relgrp)) {
            // printf("Group needs to be changed\n"); 
            if( chown(argv[i], relui, relgi)) {
                printf("There's a problem chown'ing %s. Skipped.\n", argv[i]);
                continue;
            }
        }

        if(S_ISDIR(sb.st_mode)) { /* we are going to forcibly set read, execution ansd setgid bits */
            sb.st_mode |= S_IRGRP;
            sb.st_mode |= S_IXGRP;
            sb.st_mode |= S_ISGID;
            chmod(argv[i], sb.st_mode);
#ifdef DBG
            printf("Have just setgid'd %s\n", argv[i]);
#endif
        }
    }
    free(cwd);
}

int main(int argc, char *argv[])
{
    if (argc == 1) {
        fprintf(stderr, "Usage: %s dir|filename\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    sanitisegrp(argc, argv);
    return 0;
}
