#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

float getSequenceHydropathy (const char* sequence, size_t start, size_t end);

bool proteinIsIntegral(const char* sequence);

float getProteinMolecularWeight(const char* sequence);

float estimateSequenceNetCharge (const char* sequence, size_t start, size_t end, float pH);
float estimateProteinNetCharge(const char* sequence, float pH);

typedef struct Range {
    float min;
    float max;
}Range;


// returns a rough min max range where pI will fall for sequence
// due to vast differences depending on folding, we are unconcerned
// with H-H calculated pI
Range estimateSequencePI (const char* sequence, size_t start, size_t end);