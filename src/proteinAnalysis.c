#include "proteinAnalysis.h"


// EVERY AA PROPERTY USES THIS KEY
// ALREADY SORTED FOR BINARY SEARCH
const char* AA_KEY = "ACDEFGHIKLMNPQRSTVWY";

// compare function for binary search
int compareAAs(const void* a, const void* b) {
    return (*(char*)a - *(char*)b);
}


// gets the index of the AA in key, applies that index to provided property list
// list must be same size, and matching AA to property in order to avoid errors
// use only macro below
void* getAAProp (const void* prop, size_t prop_size, char aa, void* error) {
      // find location (bsearch for speed)
    char* aa_location = (char*) bsearch(&aa, AA_KEY, 20, sizeof(char), compareAAs);

    if (aa_location == NULL) return error;

    // get index from pointer to location
    int aa_index = aa_location - AA_KEY;

    return (char*) prop + aa_index * prop_size;
}

int getAAIndex (char aa) {
    char* aa_location = (char*) bsearch(&aa, AA_KEY, 20, sizeof(char), compareAAs);

    // get index from pointer to location
    int aa_index = aa_location - AA_KEY;

    return aa_index;
}

/*
more usuable macro for this function
USE THIS ONLY, to avoid weird pointer arrithmatic errors

typeof(*list) = data type for property
1. calls function with the proper size for type of property, forwarding list and AA
2. converts null* to property type*
3. derefences property type pointer
*/



// list = property list
// aa = amino acid
// e = error code
#define GET_AA_PROP(list, aa, e) *(typeof(*list)*)getAAProp(list, sizeof(typeof(*list)), aa, &e);



bool searchInSeq(const char* seq, size_t start, size_t end, char aa) {
    for (size_t i = start; i <= end; i++) {
        if (seq[i] == aa) return true;
    }

    return false;
}

// returns [b]/[a] for henderson hasselbach equation
// short name, as I assume it would become essential eventually
float HH (float pH, float pKa) {
    return pow(10, pH - pKa);
}

// due to the tertiary structure of 
float estimateSequenceNetCharge (const char* sequence, size_t start, size_t end, float pH) {
    // source: https://www.vanderbilt.edu/AnS/Chemistry/Rizzo/stuff/AA/AminoAcids.html
    const float r_group_pKa[20] = {
        0, //a 
        8.37, //c
        3.90, //d
        4.07, //e
        0, 0, //f, g
        6.04, //h 
        0, //i
        10.54, //k
        0, 0, 0, 0, 0, // l, m, n, p, q
        12.48, //r
        0, 0, 0, 0,//s, t, v, w
        10.46 //y 
    }; 

    const char negative_ion_aa[4] = {
        'D', 'E', 'C', 'Y'
    };
    int negative_ion_aa_count[4] = {0, 0, 0, 0};

    const char positive_ion_aa[3] = {
        'H',
        'K',
        'R'
    };
    int positive_ion_aa_count[3] = {0, 0, 0};

    // start by counting number of aas
    for (size_t i = start; i <= end; i++) {
        char aa = sequence[i];

        // get if aa is negative
        for (size_t j = 0; j < 4; j++) {

            // if matching a negative aa
            if (aa == negative_ion_aa[j]) {

                // add one to ledger
                negative_ion_aa_count[j]++;

                // we are done searching
                break;
            }
        }

        // get if aa is positive
        for (size_t j = 0; j < 3; j++) {

            // if matching a positive aa
            if (aa == positive_ion_aa[j]) {

                // add one to ledger
                positive_ion_aa_count[j]++;

                // we are done searching
                break;
            }
        }
    }

    float net_charge = 0;

    for (size_t i = 0; i < 4; i++ ) {
        char aa = negative_ion_aa[i];
        
        // get R-group pka and run HH
        float pka_error = -100.0f;
        float aa_pKa = GET_AA_PROP(r_group_pKa, aa, pka_error);
        float b_over_a = HH(pH, aa_pKa);

        // if negative, base is negatively charged, get the truu proportion of base
        float proportion_base = b_over_a / (1 + b_over_a);
        float total_base = proportion_base * negative_ion_aa_count[i];
        net_charge -= total_base;

    }

    
    for (size_t i = 0; i < 3; i++ ) {
        char aa = positive_ion_aa[i];

        // get R-group pka and run HH
        float pka_error = -100.0f;
        float aa_pKa = GET_AA_PROP(r_group_pKa, aa, pka_error);
        float b_over_a = HH(pH, aa_pKa);
        
        // if positve, acid is positively charged
        float proportion_acid = 1 / (1 + b_over_a);
        float total_acid = proportion_acid * positive_ion_aa_count[i];
        net_charge += total_acid;

    }

    return net_charge;
}

float estimateProteinNetCharge(const char* sequence, float pH){
    int protein_length = strlen(sequence);
    return estimateSequenceNetCharge(sequence, 0, protein_length, pH);
}


// I actually dont know how to do an error code for this
// - 1 won't work
float getSequenceHydropathy (const char* sequence, size_t start, size_t end) {

    const float KDH[20] = {
        1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6,
        -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3
    };

    if (strlen(sequence) <= start) {
        printf("ERROR (getSequenceHydropathy): end index out of bounds");
        return -1;
    }

    float hp_sum = 0;

    for (size_t i = start; i <= end; i++) {

        // AA we are iterating through
        char aa = sequence[i];

        float error = -100.0;

        float aa_hp = GET_AA_PROP(KDH, aa, error);
        //sum 
        if (aa_hp != -100.0) {
            hp_sum +=  aa_hp;
        }
    }

    // return the average of hydropathy sum
    return hp_sum / ((end - start) + 1);
}


float getProteinMolecularWeight(const char* sequence) {

    int seq_length = strlen(sequence);

    float aa_molecular_weights[20] = {
        89.10, // a
        121.16, // c
        133.11, // d
        147.13, // e
        165.19, // f
        75.07, // g
        155.16, //h
        131.18, // i
        146.19, // k
        131.18, // l
        149.21, // m
        132.12, // n
        115.13, // p
        146.15, // q
        174.20, // r
        105.09, // s
        119.12, // t
        117.15, // v
        204.23, // w
        181.19, // y 
    };

    int aa_counts [20];
    memset(aa_counts, 0, 20 * sizeof(int));

    float molecular_weight = 0;
    for (int i = 0; i < seq_length; i++) {

        float mw_error = -1.0f;
        float aa_mw = GET_AA_PROP(aa_molecular_weights, sequence[i], mw_error);

        if (aa_mw != mw_error) {
            molecular_weight += GET_AA_PROP(aa_molecular_weights, sequence[i], mw_error);
        }

    }

    return molecular_weight;
}



bool proteinIsIntegral(const char* sequence) {

    int seq_length = strlen(sequence);

    // if less than this, protein too small to have both 
    // a signal and transmembrane sequence
    if (seq_length < 42) {
        return false;
    }

    //  -- find signal sequence -- 

    int start_frame = 0;
    // we assume it is 8 long
    int end_frame = start_frame + 8 - 1;

    bool signal_seq_found = false;

    for (;                                                             // we have already declared the iterators    
        end_frame <= 30 && !signal_seq_found;                          // until we reach the end, or when we have found a sequence
        start_frame++, end_frame++                                    
    ) {
        // get the hydropathy
        float frame_kb = getSequenceHydropathy(sequence, start_frame, end_frame);

        // if it is >= 2.5 we have found the sequence
        signal_seq_found = frame_kb >= 2.5;
    }

    // check sequence for lack of prolines
    if (signal_seq_found) {
        signal_seq_found = !searchInSeq(sequence, start_frame, end_frame, 'p');
    }

    // if signal sequence hasn't been found, no use searching for transmembrane seq
    if (!signal_seq_found) return false;


    // -- search for transmembrane seq --
    start_frame = 30;
    end_frame = start_frame + 11 - 1;
    bool transmembrane_seq_found = false;

    for (;                                                             // we have already declared the iterators 
        end_frame < seq_length && !transmembrane_seq_found;                   // unil we have reached the end, or when we have found a sequence
        start_frame++, end_frame++
    ) {

        // get the hydropathy
        float frame_kb = getSequenceHydropathy(sequence, start_frame, end_frame);

        // if it is >= 2 we have found the sequence
        transmembrane_seq_found = frame_kb >= 2.0;
    }

    // search sequence for lack of prolines
    if (transmembrane_seq_found) {
        transmembrane_seq_found = !searchInSeq(sequence, start_frame, end_frame, 'p');
    }

    return transmembrane_seq_found;
}