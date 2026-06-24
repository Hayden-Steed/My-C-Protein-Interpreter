#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <malloc.h>

typedef struct SafeStr {
    char* content;
    int length; // current length of string
    int allocated; // current allocated # of chars
} SafeStr;

SafeStr SafeStr_createFromExpr (char* init);

SafeStr SafeStr_createEmpty(int size);

// appends by one, reallocates by doubling if size exeeded, FALSE = failed
void SafeStr_appendChar(SafeStr* safe_str, char c);

// fits safe string allocated to exactly length + 1
void SafeString_fit(SafeStr* safe_str);

// appends by new string, exact size if reallocation needed
void SafeString_appendString(SafeStr* safeStr, char* string);


// searches string for keyword, returns index upon finding
int SafeString_findKey (SafeStr* safe_str, char* id);

