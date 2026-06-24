#include <stdio.h>
#include "safeString.h"
#include "fileHandling.h"
#include "FASTA.h"
#include "proteinAnalysis.h"


// example use
// exe -f Data/ecoli.faa -o output.csv -im


#define ARGS_INIT {.FN_ARG = "-f", \
    .OP_ARG = "-o", \
    .IMA_ARG = "-im",\
    .intermembrane_analysis = false,\
    .MW_ARG = "-mw", \
    .mw_analysis = false, \
    .NC_ARG = "-nc", \
    .net_charge_analysis = false,\
}

// only create using initializer above
typedef struct Args {
    char* file_name; 
    char* FN_ARG; // effectively constant

    char* output_name;
    char* OP_ARG;

    bool intermembrane_analysis;
    char* IMA_ARG;

    bool mw_analysis;
    char* MW_ARG;

    bool net_charge_analysis;
    char* NC_ARG;
    float net_charge_ph;

} Args;

Args interpretArgs (int argc, char** argv) {
    Args args = ARGS_INIT;

    // run through each argument
    for (int i = 0; i < argc; i++ ) {

        // check for argument in struct
        if (strcasecmp(argv[i], args.FN_ARG) == 0) {
            // in the case of arg with extra value 
            // get value after and add it
            i++;
            args.file_name = argv[i];
        }else if (strcasecmp(argv[i], args.OP_ARG) == 0) {
            i++;
            args.output_name = argv[i];
        }else if (strcasecmp(argv[i], args.IMA_ARG) == 0) {
            args.intermembrane_analysis = true;
        }else if (strcasecmp(argv[i], args.MW_ARG) == 0) {
            args.mw_analysis = true;
        } else if (strcasecmp(argv[i], args.NC_ARG) == 0) {
            args.net_charge_analysis = true;
            i++;
            args.net_charge_ph = atof(argv[i]);
        }else {
            printf("arg '%s' not recognized\n", argv[i]);
        }
    }
    
    return args;
}

int main (int argc, char** argv) {

    Args args = interpretArgs(argc, argv);

    char* file_path = args.file_name;
    FILE* i_fptr = fopen64(file_path, "r");

    //ERROR - file could not open
    if (i_fptr == NULL) {
        printf("ERROR (main): file could not open '%s'\n", file_path);
        exit(1);
    }    

    // interpret fasta file
    FASTA fasta = FASTA_ReadFromFile(i_fptr);    
    fclose(i_fptr);


    // if argument for integral protein analysis used
    // analyze proteins for if they are integral
    bool* is_integral;

    if (args.intermembrane_analysis) {
        // initialize data
        is_integral = (bool*)malloc(fasta.p_count * sizeof(bool));
        memset(is_integral, 0, fasta.p_count * sizeof(bool));

        // iterate through each protein
        for (int i = 0; i < fasta.p_count; i++) {
            ProteinEntry pe = fasta.entrys[i];

            // get if its integral
            bool integral = proteinIsIntegral(pe.AAs);
            // add to data
            is_integral[i] = integral;
        }
    }

    // if argument for molecular weights used
    // calculate mw for each protein
    float* mws;

    if (args.mw_analysis) {
        // initialize data
        mws = (float*)malloc(fasta.p_count * sizeof(float));
        memset(mws, 0, fasta.p_count * sizeof(float));

        // iterate through each protein
        for (int i = 0; i < fasta.p_count; i++) {
            ProteinEntry pe = fasta.entrys[i];

            // get mw 
            float mw = getProteinMolecularWeight(pe.AAs);

            // add to data
            mws[i] = mw;
        }
    }
    

    // if argument for net charges is used
    // calculate net charge for each protein
    float* net_charges;

    if (args.net_charge_analysis) {
        // initialize data
        net_charges = (float*)malloc(fasta.p_count * sizeof(float));
        memset(net_charges, 0, fasta.p_count * sizeof(float));

        // iterate through each protein
        for (int i = 0; i < fasta.p_count; i++) {
            ProteinEntry pe = fasta.entrys[i];

            // get net charge at given ph
            float net_charge = estimateProteinNetCharge(pe.AAs, args.net_charge_ph);

            // add to data
            net_charges[i] = net_charge;
        }
    }
    

    file_path = args.output_name;
    FILE* o_fptr = fopen64(file_path, "w+");
    if (o_fptr == NULL) {
        printf("ERROR (main): file could not be created '%s'\n", file_path);
    }

    for (int protein_index = 0; protein_index < fasta.p_count; protein_index++) {
        ProteinEntry pe = fasta.entrys[protein_index];

        // print out identifier and sequence
        fprintf(o_fptr, "%s", pe.refSeq);
        fprintf(o_fptr, ", %s", pe.AAs);
        

        // output any additional data requested
        if (args.mw_analysis) {
            float mw = mws[protein_index];
            fprintf(o_fptr, ", %f (g/mol)", mw);
        }
         if (args.net_charge_analysis) {
            float net_charge = net_charges[protein_index];
            fprintf(o_fptr, ", %f (charge)", net_charge);
        }

        if (args.intermembrane_analysis) {
            bool integ = is_integral[protein_index];
            fprintf(o_fptr, ", %s", integ ? "integral" : "non integral");
        }

        // new line
        fprintf(o_fptr, "\n");
    }
    
    return 0;
}