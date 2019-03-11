/* @bkb3, licensed under GNU GPLv3 or later at your
opinion*/

/*this intercepts file write and redirects output 
to stdout stream. compilation by gcc as follows:

gcc -shared -fPIC -ldl -o loader.so loader.c*/


#include <stdio.h>
#include <stdarg.h>


FILE *fopen(const char *path, const char *mode) {
    /* instead of returning a file pointer to disk,*/
    /* we return the pointer to stderr stream      */
    /*technically we could write to stdout directly*/
    /*but this may create problems if the process  */
    /*needs to write to stdout on occasions        */
    return stderr;
}


int fprintf(FILE *stream, const char *format, ...) {
    /* redirect data towards file to stdout        */
    va_list arg;
    int done;
    va_start (arg, format);
    done = vfprintf (stdout, format, arg);
    va_end (arg);

   return done;
}
