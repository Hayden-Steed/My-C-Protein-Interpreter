#include "safeString.h"
#include "FASTA.h"
#include "fileHandling.h"


ProteinEntry ProteinEntry_ReadFromFile (FILE* fptr) {
    ProteinEntry pe;

    // first line in entry is description
    SafeStr ss_desc_line = safeStrFromFileLine(fptr);
    int desc_length = ss_desc_line.length;
    char* desc_line = ss_desc_line.content;

    // if it doesn't start with a >, we are not at a new protein entry
    if (desc_line[0] != '>') {
        printf("Error (ProteinEntry_ReadFromFile): fptr not at start of next protein\n");
        exit(1);
    }

    SafeStr ss_ref_seq = SafeStr_createEmpty(10);

    // line pointer, 1 to skip past '>'
    int line_pointer = 1;

    // copy line until space
    while (line_pointer < desc_length && desc_line[line_pointer] != ' ') {
        SafeStr_appendChar(&ss_ref_seq, desc_line[line_pointer]);
        line_pointer++;
    }

    pe.refSeq = ss_ref_seq.content;

    // next lets get the protein sequence
    SafeStr ss_aa_seq = SafeStr_createEmpty(1);

    while (fpeek(fptr) != '>' && fpeek(fptr) != EOF) {
        SafeStr next_line = safeStrFromFileLine(fptr);
        
        SafeString_appendString(&ss_aa_seq, next_line.content);
    }

    pe.AAs = ss_aa_seq.content;
    
    return pe;
}

FASTA FASTA_ReadFromFile(FILE* fptr) {
    // get approximate number of proteins so we can initialize memory once
    printf("getting protein count....\n\n");
    int p_count = 0;

    int c = '\0';
    while (c != EOF) {
        
        if (c == '>') {
            p_count ++;
            printf("protein count = %d\r", p_count);
        }
        
        c = fgetc(fptr);
    } 

    rewind(fptr);
    
    printf("protein count = %d\n", p_count);
    // reading 
    FASTA fasta;
    fasta.p_count = p_count;

    ProteinEntry* entrys = (ProteinEntry*)malloc(p_count * sizeof(ProteinEntry));

    printf("interpeting FASTA file...\n\n");


    // iterate through getting protein enetrys until we reach end of file
    int p_index = 0;
    while (fpeek(fptr) != EOF) {
        printf("(%d/%d) analyzed\r",p_index, p_count);

        ProteinEntry PE = ProteinEntry_ReadFromFile(fptr);

        entrys[p_index] = PE;

        p_index++;
    }

    printf("complete!               \n\n");

    fasta.entrys = entrys;
    return fasta;
}

ProteinEntry* FASTA_GetEntryFromReference (FASTA fasta, char* ref) {
    
    // search through every entry
    for (int i = 0; i < fasta.p_count; i++) {
        // compare id to given string
        if (strcmp(fasta.entrys[i].refSeq, ref) == 0) {

            // return if matched
            return &fasta.entrys[i];
        }
    }

    // return Null otherwise
    return NULL;
}
