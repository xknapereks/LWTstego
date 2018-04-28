#ifndef EMBED_H_INCLUDED
#define EMBED_H_INCLUDED

const int starting_depth = 0;

const int comp_th_16bit = 30000;
const int comp_th_24bit = 8388000;
const float comp_rat_inv = 0.1;
const float int16_max = 32768.0;

//const int windowSize = 4096;
const int windowSize = 28672;

const bool kfiltering = true;
const bool fixEmbParams = false;
const int c_d1 = 5;
const int c_d2 = 3;
const int levelOffset = 42; // in dB

const unsigned write_mask[24] = {0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFFFFC, 0xFFFFFFF8,
                                 0xFFFFFFF0, 0xFFFFFFE0, 0xFFFFFFC0, 0xFFFFFF80,
                                 0xFFFFFF00, 0xFFFFFE00, 0xFFFFFC00, 0xFFFFF800,
                                 0xFFFFF000, 0xFFFFE000, 0xFFFFC000, 0xFFFF8000,
                                 0xFFFF0000, 0xFFFE0000, 0xFFFC0000, 0xFFF80000,
                                 0xFFF00000, 0xFFE00000, 0xFFC00000, 0xFF800000};

struct {
    const int thresholds[12] {-85, -76, -70, -63, -57, -51, -45, -39, -33, -27, -21, -15};
    const int d1d2[12][2] {{1, 0}, {3, 1}, {4, 2}, {5, 3}, {6, 4}, {7, 5}, {8, 6}, {9, 7}, {10, 8}, {11, 9}, {12, 10}, {13, 11}};
} embTable_16bit;

struct {
    const int thresholds[20] {-133, -124, -118, -112, -105, -99, -93, -87, -81, -75, -69, -63, -57, -51, -45, -39, -33, -27, -21, -15};
    const int d1d2[20][2] {{1, 0}, {3, 1}, {4, 2}, {5, 3}, {6, 4}, {7, 5}, {8, 6}, {9, 7}, {10, 8}, {11, 9},
                        {12, 10}, {13, 11}, {14, 12}, {15, 13}, {16, 14}, {17, 15}, {18, 16}, {19, 20}, {20, 21}, {21, 19}};
} embTable_24bit;


const int32_t header_mask1 = ~(((1 << 6) - 1) << (starting_depth)); //mask with length of detail1_scope bits
const int32_t header_mask2 = ~(((1 << 4) - 1) << (starting_depth)); //mask with length of detail2_scope bits

unsigned get_available_space(int32_t* samples, unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample);
unsigned LWT_LSBinsertion(int32_t* samples, unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample, unsigned char *messageBuffer, unsigned buffer_length);
int LWT_LSBretrieval(int32_t* samples, unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample, unsigned char *messageBuffer);


#endif // EMBED_H_INCLUDED
