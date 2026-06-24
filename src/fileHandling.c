#include "fileHandling.h"

// SOURCE =  https://stackoverflow.com/questions/7623323/how-to-take-a-look-at-the-next-character-in-a-stream-in-c
int fpeek(FILE * const fp) {
    const int c = getc(fp);
    return c == EOF ? EOF : ungetc(c, fp);
}

SafeStr safeStrFromFileLine(FILE* fptr) {
    SafeStr safe_str = SafeStr_createEmpty(10);

    while (fpeek(fptr) != '\n') {
        SafeStr_appendChar(&safe_str, fgetc(fptr));
    }

    // skip past \n char
    fgetc(fptr);

    return safe_str;
}
