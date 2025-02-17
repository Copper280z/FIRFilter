#include <stdio.h>
#include <math.h>

/*
https://www.staff.ncl.ac.uk/oliver.hinton/eee305/Chapter4.pdf
*/

#define ASSERT_FINITE(val) if ((val) != (val)) { printf("NaN at line %d\n", __LINE__); return;}

typedef struct {
    int order;
    float* coeffs;
    float sample_freq;
} FIRFilter;

float sinc(float x) {
    if (fabs(x) < 1e-8)
        return 1.0;
    else
        return sin(x) / (x);
}

float HanningWindow(int i, float nTaps) {
    float val = 0.5 - 0.5 * cos(2 * M_PI * (i - nTaps/2.0));
    return val;
}

float HammingWindow(int i, float nTaps) {
    float val = 0.54 - 0.46 * cos(2 * M_PI * (i - nTaps/2.0));
    return val;
}

void calculateFIRLowPass(float cutoffFrequency, FIRFilter *filter) {
    int nTaps = filter->order;
    float angularCutoff = M_PI * cutoffFrequency / (filter->sample_freq / 2.0);

    float sum = 0.0;
    for (int i = 0; i < nTaps; ++i) {
        float sincVal = (angularCutoff / M_PI) * sinc(((float)i-nTaps/2.0) * angularCutoff);
        float windowFunction = HanningWindow(i, nTaps);
        filter->coeffs[i] = sincVal * windowFunction;
        sum += filter->coeffs[i];
    }

    ASSERT_FINITE(sum);
    
    for (int i = 0; i < nTaps; ++i) {        
        filter->coeffs[i] /= sum;
        ASSERT_FINITE(filter->coeffs[i]);
    }

}

void save_filter_to_file(const char *fname, FIRFilter *filter) {
    FILE *file = fopen(fname, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    fprintf(file, "%f\n", filter->sample_freq);
    for (int i=0; i<filter->order; i++) {
        fprintf(file, "%f\n",filter->coeffs[i]);
    }
    
    if (fclose(file) != 0) {
        perror("Error closing file");
        return;
    }
}

// Example usage
int main() {
    float cutoffFrequency = 1000.0; // Hz
    int nTaps = 31;
    float sampleFreq = 5000.0; // Hz
    FIRFilter fir;
    float *coefficients = (float *)malloc(nTaps * sizeof(float));
    fir.coeffs = coefficients;
    fir.order = nTaps;
    fir.sample_freq = sampleFreq;
    calculateFIRLowPass(cutoffFrequency, &fir);

    printf("FIR Low-Pass Filter Coefficients:\n");
    for (int i = 0; i < nTaps; ++i)
        printf("%f, ", coefficients[i]);

    save_filter_to_file("fir_filter.txt", &fir);

    free(coefficients);
    return 0;
}
