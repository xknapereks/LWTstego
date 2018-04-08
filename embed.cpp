#include "DspFilters/Dsp.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "embed.h"
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <fstream>


using namespace std;

void write_text_to_log_file( const std::string &text ) {
    std::ofstream log_file(
        "log_file.txt", std::ios_base::out | std::ios_base::app );
    log_file << text << std::endl;
}

void pushToLogFile (int sample) {

    stringstream stream;
    //stream << fixed << setprecision(4) << (sample/32768.0);
    stream << sample;
    string s = stream.str();
    write_text_to_log_file(s);
}

int highestOneBit(int i) {
    //i = abs(i);
    if (i < 0)
        i = ~i;
    int j = 0;
    while (i > 0xFF) {
        i >>= 8;
        j += 8;
    }
    while (i > 0) {
        i >>= 1;
        j++;
    }
    return j;
    /*
    i |= (i >>  1);
    i |= (i >>  2);
    i |= (i >>  4);
    i |= (i >>  8);
    i |= (i >> 16);
    return i - ((unsigned)i >> 1);
    */
}

int16_t sign(int32_t x) {
    return (x > 0) - (x < 0);
}

int16_t staticCompressor (int32_t input) {
    //compresses input value if it's above comp_th. Compression ratio is 1/comp_rat_inv
    if (abs(input) > comp_th)
        return (sign(input) * round(pow(abs(input)/comp_th, comp_rat_inv) * comp_th));
    else
        return input;
}

bool checkOverflow(int32_t a, unsigned BitsPerSample) {
    if (BitsPerSample == 16) {
        if (abs(a) > (0x7FFF - 1)) {
            //printf("OVERFLOW: %i -> %i\n", a, b);
            return true; //overflow happened
        }
        else
            return false;
    }
    else if (BitsPerSample == 24) {
        if (abs(a) > (0x7FFFFF - 1)) {
            //printf("OVERFLOW: %i -> %i\n", a, b);
            return true; //overflow happened
        }
        else
            return false;
    }
}

int16_t fix (int32_t a) { //round towards zero
    if (a < 0) return (ceil(a));
    else return (floor(a));
}

void lifting_step (int32_t* xo, int32_t* xe) {
    //input xo, xe ... odd sample, even sample
    //output average, difference
    *xe -= *xo;
    *xo += fix(*xe/2);
}


void inv_lifting_step (int32_t* p, int32_t* d) {
    //input: p,d .... average, difference
    //output: odd sample, even sample
    *p -= fix(*d/2);
    *d += *p;
}

//float estimateLevel (int* samples, int cnt, int bitsPerSample) {
//    // computes RMS level of given integer samples of the given count cnt
//    const float normalizationConst = 1.0 * (1 << (bitsPerSample - 1));
//    float level = 0;
//    for (int i = 0; i < cnt; i++) {
//        float s = samples[i] / normalizationConst;
//        level += s*s;
////        level += s;
//    }
//    //level = log10(sqrt(level / cnt)); //RMS is defined with sqrt. However, when doing a logarithm afterwards, it's not efficient
//    level = (log10(level / cnt))*10;
//    return level;
//}

float estimateLevel (float** samples, int cnt, int bitsPerSample) {
    // computes RMS level of given float samples of the given count cnt
    //const float normalizationConst = 1.0 * (1 << (bitsPerSample - 1));
    float level = 0;
    for (int i = 0; i < cnt; i++) {
        float s = samples[0][i];
        level += s*s;
//        level += s;
    }
    //level = log10(sqrt(level / cnt)); //RMS is defined with sqrt. However, when doing a logarithm afterwards, it's not efficient
    level = (log10(level / cnt))*10;
    return level;
}

void iir_stage2 (float** samples, int cnt, float* xm_1, float* xm_2, float a1, float a2, float b0, float b1, float b2) {
    float xm = 0;
    for (int i = 0; i < cnt; i++) {
        xm = samples[0][i] - a1*(*xm_1) - a2*(*xm_2);
        samples[0][i] = b0 * xm + b1*(*xm_1) + b2*(*xm_2);
        *xm_2 = *xm_1;
        *xm_1 = xm;
    }
}

void pre_filter (float** samples, int cnt) {
    float xm_1 = 0;
    float xm_2 = 0;
    iir_stage2(samples, cnt, &xm_1, &xm_2, pre_a1, pre_a2, pre_b0, pre_b1, pre_b2);
}

void rlb_filter (float** samples, int cnt) {
    float xm_1 = 0;
    float xm_2 = 0;
    iir_stage2(samples, cnt, &xm_1, &xm_2, rlb_a1, rlb_a2, rlb_b0, rlb_b1, rlb_b2);
}

void setEmbeddingParams(int *samples, int NumChannels, int cnt, int sampleRate, int bitsPerSample, int index, int* d1, int* d2) {
    if (cnt > windowSize)
        cnt = windowSize;
    float level;
    int numSamples = cnt;
    const float normalizationConst = 1.0 * (1 << (bitsPerSample - 1));
    float* audioData[1];
//    float* audioDataFiltered[1];
    audioData[0] = new float[numSamples];
//    audioDataFiltered[0] = new float[numSamples];
    // float* audioData = new float[numSamples];
    for (int i = 0; i < numSamples; i++) {
        audioData[0][i] = samples[i] / normalizationConst;
//        audioDataFiltered[0][i] = 0;
    }

    if (fixEmbParams) {
        *d1 = c_d1;
        *d2 = c_d2;
    } else {
        if (kfiltering) {

//            pre_filter(audioData, numSamples);
//            rlb_filter(audioData, numSamples);

            {
                Dsp::Filter* f = new Dsp::FilterDesign
                  <Dsp::ChebyshevII::Design::HighShelf <5>, 1>;
                Dsp::Params params;
                params[0] = sampleRate; // sample rate
                params[1] = 5; // order
                params[2] = 1000; // corner frequency
                params[3] = 4; // shelf gain
                params[4] = 0.1; // passband ripple
                f->setParams (params);
                f->process (numSamples, audioData);
              }

                // set up a 2-channel RBJ High Pass with parameter smoothing,
                // but bypass virtual function overhead
              {


                // the difference here is that we don't go through a pointer.
                Dsp::SmoothedFilterDesign <Dsp::RBJ::Design::HighPass, 1> f (cnt);
                Dsp::Params params;
                params[0] = sampleRate; // sample rate
                params[1] = 100; // cutoff frequency
                params[2] = 0.25; // Q = 1.25
//                params[1] = 1000;
//                params[2] = 0.1;
                f.setParams (params);

                f.process (numSamples, audioData);
              }

//            Dsp::SimpleFilter <Dsp::RBJ::HighPass> f;
//            f.setup (44100, // sample rate Hz
//                     80, // cutoff frequency Hz
//                     1); // "Q" (resonance)


            level = estimateLevel(audioData, cnt, bitsPerSample);
//            level = estimateLevel(audioDataFiltered, cnt, bitsPerSample);
        } else {
            level = estimateLevel(audioData, cnt, bitsPerSample);
        }

        level -= audibleOffset;

        int i = 0;
        if (bitsPerSample == 16) {
            while (level > embTable_16bit.thresholds[i])
                i++;
            *d1 = embTable_16bit.d1d2[i][0];
            *d2 = embTable_16bit.d1d2[i][1];
        } else if (bitsPerSample == 24) {
            while (level > embTable_24bit.thresholds[i])
                i++;
            *d1 = embTable_24bit.d1d2[i][0];
            *d2 = embTable_24bit.d1d2[i][1];
        }
    }

    //printf("Output parameters: d1 = %i, d2 = %i\n", *d1, *d2);
//    printf("%i (%i): \t%.3f\t%i, %i\n", int (round(1.0*index/(sampleRate*NumChannels))), index, level + audibleOffset, *d1, *d2);
    //printf("%i, %i, %.3f, %i, %i\n", int (round(1.0*index/(sampleRate*NumChannels))), index, level + audibleOffset, *d1, *d2);
    //printf("%i ", int (round(1.0*index/(sampleRate*NumChannels))));
}


void forward_LWT (int32_t* LWTcoefs, long int pos, int NumChannels) {
    if (NumChannels == 2) { //interleaved stereo
        //left channel
        lifting_step(LWTcoefs + pos + 0, LWTcoefs + pos + 2);
        lifting_step(LWTcoefs + pos + 4, LWTcoefs + pos + 6);
        lifting_step(LWTcoefs + pos + 0, LWTcoefs + pos + 4);

        //right channel
        lifting_step(LWTcoefs + pos + 1, LWTcoefs + pos + 3);
        lifting_step(LWTcoefs + pos + 5, LWTcoefs + pos + 7);
        lifting_step(LWTcoefs + pos + 1, LWTcoefs + pos + 5);
    }
    else { //mono
        lifting_step(LWTcoefs + pos + 0, LWTcoefs + pos + 1);
        lifting_step(LWTcoefs + pos + 2, LWTcoefs + pos + 3);
        lifting_step(LWTcoefs + pos + 0, LWTcoefs + pos + 2);
    }
}

void backward_LWT (int32_t* LWTcoefs, long int pos, int NumChannels) {

    if (NumChannels == 2) { //interleaved stereo
        //left channel
        inv_lifting_step(LWTcoefs + pos + 0, LWTcoefs + pos + 4);
        inv_lifting_step(LWTcoefs + pos + 0, LWTcoefs + pos + 2);
        inv_lifting_step(LWTcoefs + pos + 4, LWTcoefs + pos + 6);

        //right channel
        inv_lifting_step(LWTcoefs + pos + 1, LWTcoefs + pos + 5);
        inv_lifting_step(LWTcoefs + pos + 1, LWTcoefs + pos + 3);
        inv_lifting_step(LWTcoefs + pos + 5, LWTcoefs + pos + 7);
    }
    else { //mono
        inv_lifting_step(LWTcoefs + pos + 0, LWTcoefs + pos + 2);
        inv_lifting_step(LWTcoefs + pos + 0, LWTcoefs + pos + 1);
        inv_lifting_step(LWTcoefs + pos + 2, LWTcoefs + pos + 3);
    }
}

unsigned get_bits_inc(unsigned char *buffer, unsigned *offset, int count) { //read sequence of bits from up to 3 adjacent bytes and output them shifted about starting_depth
    unsigned start_bit = *offset;
    unsigned end_bit = *offset + count - 1;
    unsigned start_byte = 0;
	unsigned end_byte = 0;

	while (start_bit >= 8) {
		end_bit -= 8;
		start_bit -= 8;
		start_byte++;
		end_byte++;
	}
	while (end_bit >= 8) {
		end_bit -= 8;
		end_byte++;
	}

	/* Get the bits. */
	unsigned result = (unsigned) (((unsigned char)(buffer[start_byte] << (start_bit))) >> (start_bit));
	if (start_byte != end_byte) {
		while (++start_byte != end_byte) {
			result <<= 8;
			result += buffer[start_byte];
		}
		result <<= end_bit + 1;
		result += buffer[end_byte] >> (7 - end_bit);
	} else if (end_bit != 7)
		result >>= (7 - end_bit);

    *offset += count;

	return result << starting_depth; //final placement
}

unsigned get_available_space(int32_t* samples, unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample) {
    int NumSamplesAbsolute = NumSamples * NumChannels;
    int detail1_scope;
    int detail2_scope;
    float availableSpace = 0;
    int i;
    for (i = 8; i < (NumSamplesAbsolute - windowSize); i += windowSize) {
        setEmbeddingParams(samples + i, NumChannels, NumSamplesAbsolute - i, SampleRate, BitsPerSample, i, &detail1_scope, &detail2_scope);
        availableSpace += ((4 * detail1_scope) + (2 * detail2_scope)) * (windowSize / 8);
    }
    //add residual space
    setEmbeddingParams(samples + i, NumChannels, NumSamplesAbsolute - i, SampleRate, BitsPerSample, i, &detail1_scope, &detail2_scope);
    availableSpace += ((4 * detail1_scope) + (2 * detail2_scope)) * ((NumSamplesAbsolute - i) / 8);
    return (unsigned) (floor(availableSpace/8)); //convert from bits to bytes, and return and unsigned integer
}

unsigned LWT_LSBinsertion(int32_t* samples, unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample, unsigned char *messageBuffer, unsigned buffer_length) {
    int NumSamplesAbsolute = NumSamples * NumChannels;
    int32_t* temp_samples = new int32_t[8];
    unsigned temp_mB_cursor = 0;
    bool overflow = false;
    int run = 0;
    float snr_numerator = 0;
    float snr_denominator = 0;
    int detail1_scope = 7;
    int detail2_scope = 5;
    int32_t write_mask1; //mask with length of detail1_scope bits
    int32_t write_mask2; //mask with length of detail2_scope bits
    bool isWindowBegin;
//    for (int i = 0; i < NumSamples; i++) {
//        pushToLogFile(samples[i]);
//    }
    //printf("%i\n%i\n%i\n%i\n", highestOneBit(2), highestOneBit(437), highestOneBit(24569), highestOneBit(26));

    //concealed message size - first 8 samples (4B integer). Fixed size. Dependent on starting depth only.
    forward_LWT(samples, 0, NumChannels);
    unsigned char* lengthBuffer = reinterpret_cast<unsigned char *>(&buffer_length);
    unsigned lB_cursor = 0;

    if (NumChannels == 2) { //interleaved stereo
        //apply write_mask
        samples[2] &= header_mask1;
        samples[3] &= header_mask1;
        samples[4] &= header_mask2;
        samples[5] &= header_mask2;
        samples[6] &= header_mask1;
        samples[7] &= header_mask1;

        //write length of the secret message into first block of data
        samples[2] = samples[2] | get_bits_inc(lengthBuffer, &lB_cursor, 6);
        samples[3] = samples[3] | get_bits_inc(lengthBuffer, &lB_cursor, 6);
        samples[4] = samples[4] | get_bits_inc(lengthBuffer, &lB_cursor, 4);
        samples[5] = samples[5] | get_bits_inc(lengthBuffer, &lB_cursor, 4);
        samples[6] = samples[6] | get_bits_inc(lengthBuffer, &lB_cursor, 6);
        samples[7] = samples[7] | get_bits_inc(lengthBuffer, &lB_cursor, 6);
    }
    else { //mono
        samples[1] &= header_mask1;
        samples[2] &= header_mask2;
        samples[3] &= header_mask1;

        samples[5] &= header_mask1;
        samples[6] &= header_mask2;
        samples[7] &= header_mask1;

        //write length of the secret message into first 2 blocks of data
        samples[1] = samples[1] | get_bits_inc(lengthBuffer, &lB_cursor, 6);
        samples[2] = samples[2] | get_bits_inc(lengthBuffer, &lB_cursor, 4);
        samples[3] = samples[3] | get_bits_inc(lengthBuffer, &lB_cursor, 6);
        samples[5] = samples[5] | get_bits_inc(lengthBuffer, &lB_cursor, 6);
        samples[6] = samples[6] | get_bits_inc(lengthBuffer, &lB_cursor, 4);
        samples[7] = samples[7] | get_bits_inc(lengthBuffer, &lB_cursor, 6);
    }

	backward_LWT(samples, 0, NumChannels);

    unsigned mB_cursor = 0;
    unsigned mbres = 0;
    int bitsLeftToWrite = buffer_length;


    int i = 8;
    printf("\nData insertion...\n");
    //for (i = 8; (mB_cursor / 8) < buffer_length; i += 8) {
    while(bitsLeftToWrite > 0) {
        //printf("%i/%i samples\n", i, NumSamples*2);
        isWindowBegin = ((i - 8) % windowSize) == 0;
        if (isWindowBegin) {
            setEmbeddingParams(samples + i, NumChannels, NumSamplesAbsolute - i, SampleRate, BitsPerSample, i, &detail1_scope, &detail2_scope);
        }
        memcpy(temp_samples, samples + i, 8 * sizeof(int32_t));
        temp_mB_cursor = mB_cursor;
        run = 0;
        while ((run == 0) || overflow) {

            if (overflow) {
                // apply compressor
                printf("Compressing 8 samples from sample %i\n", i);
                memcpy(samples + i, temp_samples, 8 * sizeof(int32_t)); //paste original samples (erase written data)
                mB_cursor = temp_mB_cursor; //shift message buffer cursor back
                if (run > 1) {
                    printf("ALERT: REPEATED COMPRESSION!\n"); //unwanted state
                }
                for (int j = 0; j < 8; j++)
                    samples[i + j] = staticCompressor(samples[i + j]); //compress all 8 samples in the block


                overflow = false;
            }
            run++; //we will check if compression is not being performed repeatedly above the same sample block

            forward_LWT(samples, i, NumChannels);


            if (NumChannels == 2) { //STEREO
                /*
                for (int j = 2; j < 8; j++) {
                    int order = highestOneBit(samples[i + j]);
                    if (order > orderThreshold) {
                        samples[i + j] &= write_mask[order - 2];
                        samples[i + j] = samples[i + j] | get_bits_inc(messageBuffer, &mB_cursor, order - orderThreshold);
                    }
                }
                */

                //apply write_mask
                samples[i + 2] &= write_mask[detail1_scope];
                samples[i + 3] &= write_mask[detail1_scope];
                samples[i + 4] &= write_mask[detail2_scope];
                samples[i + 5] &= write_mask[detail2_scope];
                samples[i + 6] &= write_mask[detail1_scope];
                samples[i + 7] &= write_mask[detail1_scope];

                //bits to write - change type to char
                samples[i + 2] = samples[i + 2] | get_bits_inc(messageBuffer, &mB_cursor, detail1_scope);
                samples[i + 3] = samples[i + 3] | get_bits_inc(messageBuffer, &mB_cursor, detail1_scope);
                samples[i + 4] = samples[i + 4] | get_bits_inc(messageBuffer, &mB_cursor, detail2_scope);
                samples[i + 5] = samples[i + 5] | get_bits_inc(messageBuffer, &mB_cursor, detail2_scope);
                samples[i + 6] = samples[i + 6] | get_bits_inc(messageBuffer, &mB_cursor, detail1_scope);
                samples[i + 7] = samples[i + 7] | get_bits_inc(messageBuffer, &mB_cursor, detail1_scope);


            }
            else { //MONO
                /*
                for (int j = 1; j < 8; j++) {
                    int order = highestOneBit(samples[i + j]);
                    if ((order > orderThreshold) && (j != 4)) {
                        samples[i + j] &= write_mask[order - 2];
                        samples[i + j] = samples[i + j] | get_bits_inc(messageBuffer, &mB_cursor, order - orderThreshold);
                    }
                }
                */

                samples[i + 1] &= write_mask[detail1_scope];
                samples[i + 2] &= write_mask[detail2_scope];
                samples[i + 3] &= write_mask[detail1_scope];
                samples[i + 5] &= write_mask[detail1_scope];
                samples[i + 6] &= write_mask[detail2_scope];
                samples[i + 7] &= write_mask[detail1_scope];

                //bits to write - change type to char
                samples[i + 1] = samples[i + 1] | get_bits_inc(messageBuffer, &mB_cursor, detail1_scope);
                samples[i + 2] = samples[i + 2] | get_bits_inc(messageBuffer, &mB_cursor, detail2_scope);
                samples[i + 3] = samples[i + 3] | get_bits_inc(messageBuffer, &mB_cursor, detail1_scope);
                samples[i + 5] = samples[i + 5] | get_bits_inc(messageBuffer, &mB_cursor, detail1_scope);
                samples[i + 6] = samples[i + 6] | get_bits_inc(messageBuffer, &mB_cursor, detail2_scope);
                samples[i + 7] = samples[i + 7] | get_bits_inc(messageBuffer, &mB_cursor, detail1_scope);

            }


            backward_LWT(samples, i, NumChannels);

            for (int j = 0; (j < 8 && !overflow); j++) {
                overflow = overflow || checkOverflow(samples[i + j], BitsPerSample);
            }

            //compute SNR
            if (!overflow) {
                for (int j = 0; j < 8; j++) {
                    snr_numerator += pow(samples[i + j] / int16_max, 2);
                    snr_denominator += pow((samples[i + j] - temp_samples[j]) / int16_max, 2);
                }
            }
        }
        bitsLeftToWrite = buffer_length - (mB_cursor / 8);
        i += 8;
    }
    //printf("Embedded to %i samples out of %i.\n", i, NumSamplesAbsolute);
    float snr = 10 * log10(snr_numerator / snr_denominator);
    printf("\nSNR: %f\n", snr);

    return i;
}

/*
Example of usage: if three bits, 011 are to be inserted to *buffer, then:
bits = 0110 0000 0000 0000 (shifted towards MSB) and
bits_length = 3
Important is, that all following bits are 0
*/
void push_bits_inc(unsigned char *buffer, unsigned bits, int bits_length, unsigned *bit_cursor) {
	unsigned byte_cursor = 0;
	unsigned start_bit = *bit_cursor;

	while (start_bit >= 8) {
		start_bit -= 8;
		byte_cursor++;
	}
	*bit_cursor += bits_length;

	// loop to write bits to a single byte. Parameters to hand over: byte_cursor, start_bit, bits, bits_length

	while (bits_length > 0) {
		if (start_bit == 0) //write all zeros to byte
			buffer[byte_cursor] = 0x00;

		//bits to write - change type to unsigned char. Important is, that bits is unsigned, otherwise arithmetic shift would add bits '1' from left hand side
		buffer[byte_cursor] = buffer[byte_cursor] | ((unsigned char) (bits >> (24 + start_bit)));
		bits_length -= 8 - start_bit;

		//following is necessary for the case of another run of the loop
		bits <<= 8 - start_bit; //erase written bits by shifting 'bits' variable
		start_bit = 0;
		byte_cursor++;
	}
}

int LWT_LSBretrieval(int32_t* samples, unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample, unsigned char *messageBuffer) {
    int NumSamplesAbsolute = NumSamples * NumChannels;
    int detail1_scope = 7;
    int detail2_scope = 5;

    forward_LWT(samples, 0, NumChannels);

    unsigned char* lengthBuffer = new unsigned char[4];
    unsigned lB_cursor = 0;

    if (NumChannels == 2) { //interleaved stereo
        push_bits_inc(lengthBuffer, (unsigned) ((samples[2] >> starting_depth) << 26), 6, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[3] >> starting_depth) << 26), 6, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[4] >> starting_depth) << 28), 4, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[5] >> starting_depth) << 28), 4, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[6] >> starting_depth) << 26), 6, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[7] >> starting_depth) << 26), 6, &lB_cursor);
    }
    else { //mono
        push_bits_inc(lengthBuffer, (unsigned) ((samples[1] >> starting_depth) << 26), 6, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[2] >> starting_depth) << 28), 4, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[3] >> starting_depth) << 26), 6, &lB_cursor);

        push_bits_inc(lengthBuffer, (unsigned) ((samples[5] >> starting_depth) << 26), 6, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[6] >> starting_depth) << 28), 4, &lB_cursor);
        push_bits_inc(lengthBuffer, (unsigned) ((samples[7] >> starting_depth) << 26), 6, &lB_cursor);

    }


    int buffer_length = *((int *) lengthBuffer);


    unsigned mB_cursor = 0;
    printf("\nData retrieval...\n");
    for (int i = 8; (mB_cursor / 8) < buffer_length; i += 8) {
        //printf("%i/%i samples\n", i, NumSamples);
        if (((i - 8) % windowSize) == 0) {
            setEmbeddingParams(samples + i, NumChannels, NumSamplesAbsolute - i, SampleRate, BitsPerSample, i, &detail1_scope, &detail2_scope);
        }

        forward_LWT(samples, i, NumChannels);
        if (NumChannels == 2) { //interleaved stereo
            /*
            for (int j = 2; j < 8; j++) {
                int order = highestOneBit(samples[i + j]);
                    if (order > orderThreshold) {
                        push_bits_inc(messageBuffer, (unsigned) ((samples[i + j]) << (32-(order-orderThreshold))), order-orderThreshold, &mB_cursor);
                    }
            }
            */

            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 2] >> starting_depth) << (32-detail1_scope)), detail1_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 3] >> starting_depth) << (32-detail1_scope)), detail1_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 4] >> starting_depth) << (32-detail2_scope)), detail2_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 5] >> starting_depth) << (32-detail2_scope)), detail2_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 6] >> starting_depth) << (32-detail1_scope)), detail1_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 7] >> starting_depth) << (32-detail1_scope)), detail1_scope, &mB_cursor);

        }
        else { //mono
            /*
            for (int j = 1; j < 8; j++) {
                int order = highestOneBit(samples[i + j]);
                    if ((order > orderThreshold) && (j != 4)) {
                        push_bits_inc(messageBuffer, (unsigned) ((samples[i + j]) << (32-(order-orderThreshold))), order-orderThreshold, &mB_cursor);
                    }
            }
            */

            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 1] >> starting_depth) << (32-detail1_scope)), detail1_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 2] >> starting_depth) << (32-detail2_scope)), detail2_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 3] >> starting_depth) << (32-detail1_scope)), detail1_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 5] >> starting_depth) << (32-detail1_scope)), detail1_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 6] >> starting_depth) << (32-detail2_scope)), detail2_scope, &mB_cursor);
            push_bits_inc(messageBuffer, (unsigned) ((samples[i + 7] >> starting_depth) << (32-detail1_scope)), detail1_scope, &mB_cursor);

        }
    }
    return buffer_length;
}
