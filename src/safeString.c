#include "safeString.h"


SafeStr SafeStr_createEmpty(int size) {
    SafeStr safe_str;
    safe_str.length = 0;
    safe_str.allocated = size + 1;

    char* memory = (char*) malloc(safe_str.allocated * sizeof(char));
    memset(memory, '\0', size * sizeof(char));


    if (memory == NULL) {
        printf("ERROR (SafeStr_createEmpty): could not allocate space\n");
        exit(1);
    }

    safe_str.content = memory;

   
    return safe_str;
}


// NOTE: expensive (realloc called), used only when necessary
SafeStr SafeStr_createFromExpr (char* init) {
    SafeStr safe_str; 

    safe_str.length = strlen(init);
    safe_str.allocated = safe_str.length + 1;

    // initialze memory for string literal, copy to memory
    safe_str.content = (char*)malloc(safe_str.allocated * sizeof(char));    

     if (safe_str.content == NULL) {
        printf("ERROR (SafeStr_createFromExpr): could not allocate space\n");
        exit(1);
    }


    strcpy(safe_str.content, init);

    return safe_str;
}

void SafeStr_appendChar (SafeStr* safe_str, char c) {
    
    SafeStr new_safe_str = *safe_str;
    int new_length = new_safe_str.length + 1;


    // if we are to exeed bounds of string, reallocate to double the memory
    if (new_length + 1 >= new_safe_str.allocated) {
        // only double the non \0 terminal, then re-add it
        int new_alloc = (new_safe_str.allocated - 1) * 2 + 1;

        char* new_mem = (char*) realloc(new_safe_str.content, new_alloc * sizeof(char));
        
        if (new_mem == NULL) {
            printf("ERROR: SafeStr_appendChar could not reallocate %s\n", new_safe_str.content);
            exit(1);
        }

        new_safe_str.content = new_mem;
        new_safe_str.allocated = new_alloc;
    }

    // write character (over \0)
    new_safe_str.content[new_safe_str.length] = c;

    // re-add \0
    new_safe_str.content[new_safe_str.length + 1] = '\0';

    new_safe_str.length = new_length;
    
    *safe_str = new_safe_str;

}

void SafeString_fit(SafeStr* safe_str) {
    SafeStr new_safe_str = *safe_str;

    // if already exactly allocated we don't need to do anything
    if (new_safe_str.length + 1 == new_safe_str.allocated) {
        return;
    }


    int new_alloc = new_safe_str.length + 1;
    char* new_mem = realloc(new_safe_str.content, new_alloc);

    if (new_mem == NULL) {
        printf("ERROR: SafeStr_fit could not reallocate %s\n", new_safe_str.content);
        exit(1);
    }

    new_safe_str.content = new_mem;

    *safe_str = new_safe_str;
}

void SafeString_appendString(SafeStr* safe_str, char* string) {
    SafeStr new_safe_str = *safe_str;

    int append_length = strlen(string);
    int space_needed = new_safe_str.length + append_length + 1;
    if (space_needed >= new_safe_str.allocated) {

        char* new_mem = (char*) realloc(new_safe_str.content, space_needed * sizeof(char));
            
        if (new_mem == NULL) {
            printf("ERROR: SafeStr_appendString could not reallocate %s\n", new_safe_str.content);
            exit(1);
        }

        new_safe_str.content = new_mem;
        new_safe_str.allocated = space_needed;
    }
    
    
    strcat(new_safe_str.content, string);

    new_safe_str.length = new_safe_str.length + append_length;

    *safe_str = new_safe_str;
}

