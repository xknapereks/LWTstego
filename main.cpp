#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <string.h>
#include <embed.h>
#include <limits.h>
#include <strings.h>

using namespace std;

unsigned char *get_file_uchar(const char *dir, unsigned *buffer_length)
{
	std::ifstream file(dir, std::ios::in | std::ios::binary);
	file.seekg(0, file.end);
	*buffer_length = file.tellg();
	file.seekg(0, file.beg);
	char *buffer = new char[*buffer_length];
	file.read(buffer, *buffer_length);
	file.close();
	return reinterpret_cast<unsigned char *>(buffer);
}

int16_t *get_file_int16(const char *dir, unsigned *buffer_length)
{
	std::ifstream file(dir, std::ios::in | std::ios::binary);
	file.seekg(0, file.end);
	*buffer_length = file.tellg();
	file.seekg(0, file.beg);
	char *buffer = new char[*buffer_length];
	file.read(buffer, *buffer_length);
	file.close();
	return reinterpret_cast<int16_t *>(buffer);
}

int *get_file_int(const char *dir, unsigned *buffer_length)
{
	std::ifstream file(dir, std::ios::in | std::ios::binary);
	file.seekg(0, file.end);
	*buffer_length = file.tellg();
	file.seekg(0, file.beg);
	char *buffer = new char[*buffer_length];
	file.read(buffer, *buffer_length);
	file.close();
	return reinterpret_cast<int *>(buffer);
}

char* decimalNumToChar(int num) {
    char strarray[] = "00";
    if ((num > 99) || (num < 0))
        return strarray;
    strarray[0] = '0' + (num / 10);
    strarray[1] = '0' + (num % 10);
    return strarray;
}

void set_stego_file(char *dir, unsigned char *buffer, unsigned buffer_length)
{
//    dir[0] = '_';
//    char array1[strlen(dir)];

    char array3[] = "00";
    if ((audibleOffset < 99) && (audibleOffset > 0)) {
        array3[0] = '0' + (audibleOffset / 10);
        array3[1] = '0' + (audibleOffset % 10);
    }

    char array1[std::char_traits<char>::length(dir)];
    memset(array1, '\0', strlen(dir));
    strncpy(array1, dir, strlen(dir) - 4);
    char array2[] = "_stego";
    char array4[] = ".wav";
    char * newDir = new char[strlen(dir)+5+1 + 2];
    strcpy(newDir,array1);
    strcat(newDir,array2);
    strcat(newDir,array3);
    strcat(newDir,array4);
//    printf("\n%s", newDir);
    std::ofstream file(newDir, std::ios::out | std::ios::binary);
    delete [] newDir;
    file.write(reinterpret_cast<const char *>(buffer), buffer_length);
    file.close();
}

void set_retrieved_file(char *dir, unsigned char *buffer, unsigned buffer_length)
{
//    dir[0] = '_';
    char array1[strlen(dir)];
    memset(array1, '\0', strlen(dir));
    strncpy(array1, dir, strlen(dir) - 4);
    char array2[] = "_retrieved";
    char array3[8];
    memset(array3, '\0', 8);
    strncpy(array3, dir + strlen(dir) - 4, 4);
//    char* newDir = new char[strlen(dir)+9+1];
    char* newDir = new char[strlen(dir)+10];
    memset(newDir, '\0', strlen(dir) + 10);
    strcpy(newDir, array1);
    strcat(newDir, array2);
    strcat(newDir, array3);
//    printf("\n%s", newDir);
    std::ofstream file(newDir, std::ios::out | std::ios::binary);
    delete [] newDir;
    file.write(reinterpret_cast<const char *>(buffer), buffer_length);
    file.close();
}

int main(int argc, char **argv)
{
	if (argc > 3) {
		printf("Unexpected number of arguments.\n");
		return -1;
	} else if (argc == 1) {
		printf("No directory specified.\n");
		return -1;
	}

	try {
	    //load cover audio
	    unsigned wav_length = 0;
        unsigned char *wav = get_file_uchar(argv[1], &wav_length);

        int pos = 12;   // First Subchunk ID from 12 to 16

        // Keep iterating until we find the fmt chunk
        while(!(wav[pos]==0x66 && wav[pos+1]==0x6d && wav[pos+2]==0x74 && wav[pos+3]==0x20)) {
            pos += 4;
            int chunkSize = wav[pos] | (wav[pos + 1] << 8) | (wav[pos + 2] << 16) | (wav[pos + 3] << 24);
            pos += 4 + chunkSize;
        }

        /* READ "FMT " SUBCHUNK */
        int NumChannels = wav[pos + 10];
        int SampleRate = wav[pos + 12] | (wav[pos + 13] << 8);
        int BitsPerSample = wav[pos + 22];
        pos += 4;
        int Subchunk1Size = *(reinterpret_cast<int*>(wav+pos));
        pos += 4 + Subchunk1Size;

        printf("#############################################\n");
        printf("Cover audio properties:\n");
        printf("Name: %s\n", argv[1]);
        printf("Number of channels: %i.\nBits per sample: %i\nSample rate: %i\nFile size: %.2f kB\n", NumChannels, BitsPerSample, SampleRate, wav_length/1024.0);

        if ((NumChannels > 2)) {
            printf("File type not supported.\n");
            return 0;
        }

        // Keep iterating until we find the data chunk (i.e. 64 61 74 61 ...... (i.e. 100 97 116 97 in decimal))
        while(!(wav[pos]==100 && wav[pos+1]==97 && wav[pos+2]==116 && wav[pos+3]==97)) {
            pos += 4;
            int chunkSize = wav[pos] | (wav[pos + 1] << 8) | (wav[pos + 2] << 16) | (wav[pos + 3] << 24);
            pos += 4 + chunkSize;
        }
        pos += 4;
        int Subchunk2Size = *(reinterpret_cast<int*>(wav+pos));
        unsigned NumSamples = Subchunk2Size * 8 / BitsPerSample / NumChannels;
        printf("Audio data size: %.2f kB\nNumber of samples in each channel: %i\n", Subchunk2Size / 1024.0, NumSamples);
        printf("Audio duration: %.2f s\n", 1.0 * NumSamples/SampleRate);
        pos += 4;

        int16_t *wav_int16 = reinterpret_cast<int16_t *>(wav); //will work for 16-bit wav only!
        //all samples are now accessible. Note, that "left" and "right" samples are interleaved.
        //NumSamples is known

        int NumSamplesAbsolute = NumSamples * NumChannels;
        int32_t* wav_int = new int32_t[NumSamplesAbsolute];
        int j = 0;
        int i = pos;

        if (BitsPerSample == 24)
            while (j < NumSamplesAbsolute) {
                wav_int[j] = wav[i] | (wav[i+1] << 8) | (wav[i+2] << 16);

                if (wav_int[j] > 0x7FFFFF)
                    wav_int[j] |= 0xFF000000;

                j++;
                i += 3;
            }
        else if (BitsPerSample == 16)
            while (j < NumSamplesAbsolute) {
                wav_int[j] = wav[i] | (wav[i+1] << 8);
                if (wav_int[j] > 0x7FFF)
                    wav_int[j] |= 0xFFFF0000;
                j++;
                i += 2;
            }
        else
            printf("Invalid BitsPerSample!\n");

		//load secret data
		unsigned dataToInsert_length;
//        int16_t *dataToInsert = get_file_int16(argv[2], &dataToInsert_length);
        unsigned char *dataToInsert = get_file_uchar(argv[2], &dataToInsert_length);

        unsigned availableSpace = get_available_space(wav_int, NumSamples, NumChannels, SampleRate, BitsPerSample);
//        unsigned int availableSpace = (NumSamples * NumChannels * (2*detail1_scope + detail2_scope) / 32 - 4); //how many bits can be embedded into the cover audio
        printf("\nEmbedding properties:\n");
//        printf("starting_depth: %i; detail1_scope = %i; detail2_scope = %i\n", starting_depth, detail1_scope, detail2_scope);
        printf("Chosen audible offset: %i\n", audibleOffset);
        printf("%.2f kB available for secret data (%.2f %% of file size)\n", availableSpace / 1024.0, 100.0 * availableSpace / wav_length);
        printf("%.2f kB to insert\n", dataToInsert_length / 1024.0);



        if (dataToInsert_length > availableSpace) {
            printf("ERROR. Data to insert too long!\n");

            return 0;
        }

        //insert data to LWT coefficients
        unsigned changedSamplesCount = LWT_LSBinsertion(wav_int, NumSamples, NumChannels, SampleRate, BitsPerSample, dataToInsert, dataToInsert_length);
        /////////////////////////////////////////////////////////////////////////////////

        printf("Successfully embedded into first %.0f s of audio.\n", 1.0*changedSamplesCount/(NumChannels*SampleRate));
        //printf("%i\n%i\n%.2f\n\n", changedSamplesCount, (NumSamples * NumChannels), 1.0*changedSamplesCount/(SampleRate * NumChannels));

        // produce the final bitstream by overwriting the data subchunk
        int ii = 0;
        int jj = pos;
        if (BitsPerSample == 24) {
            while (jj < wav_length) {
                wav[jj] = wav_int[ii];
                wav[jj + 1] = wav_int[ii] >> 8;
                wav[jj + 2] = wav_int[ii] >> 16;
                jj += 3;
                ii++;
            }
        } else if (BitsPerSample == 16) {
            while (jj < wav_length) {
                /*if (wav_int[ii] < 0) { //no conversion is now necessary. The upper bits in the 32-bit variable are all ones and can be dropped
                    wav_int[ii] = ~(-wav_int[ii]) + 1;
                }*/

                wav[jj] = wav_int[ii]; //less valuable byte goes first (little endian)
                wav[jj + 1] = wav_int[ii] >> 8;
                jj += 2;
                ii++;
            }
        }

        set_stego_file(argv[1], wav, wav_length);
        printf("successfully written to file\n");


        //unsigned char* retrievedData = new unsigned char[availableSpace]; //availableSpace can be computed from stego audio as well
        unsigned char* retrievedData = new unsigned char[wav_length/4]; //rough definition of space for retrieved data
        int retrievedData_length = LWT_LSBretrieval(wav_int, NumSamples, NumChannels, SampleRate, BitsPerSample, retrievedData);
        printf("Number of retrieved kB: %.2f\n", retrievedData_length / 1024.0);

//        printf("\nSECRET\tRETRIEVED\n");
//        for (int i = 0; i < 100; i++) {
//            printf("%.2X %.2X\n", dataToInsert[i], retrievedData[i]);
//        }

        //retrieve hidden data
        set_retrieved_file(argv[2], retrievedData, retrievedData_length);
        /////////////////////////////////////////////////////////////////


    } catch (std::bad_alloc) {
        printf("File does not exist.\n");
        return -1;
	}

	return 0;
}
