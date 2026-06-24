

typedef struct ProteinEntry {
    char* refSeq;
    int existance;
    int sequenceVersion;
    char* AAs;
} ProteinEntry;

typedef struct FASTA {
    ProteinEntry* entrys;
    char* idenifiers; 
    int p_count;
} FASTA;


ProteinEntry ProteinEntry_ReadFromFile (FILE* fptr);

FASTA FASTA_ReadFromFile(FILE* fptr);

ProteinEntry* FASTA_GetEntryFromReference (FASTA fasta, char* uni_prot_kb);

