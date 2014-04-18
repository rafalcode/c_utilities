/* u2gperms: Simple use of linux api to read the owner permissions on a file
 and copy them to the group permissiosn of the same file.

 Copyright (C) Ramón Fallon

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 Further notes:
 - Compile with "gcc -Wall -o u2gperms u2gperms.c"
 - Normally called with the * wildcard from bash, (so hidden files and dirs will not be affected)
 and mirrors the file's or dir's user permissions, to its group permissions 
*/
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
    uint i, up, gp, op, nwtp; /* unsigneds under types.h */

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
