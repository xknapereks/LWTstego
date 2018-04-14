# LWTstego
Dynamically Adjusting Insertion into LWT coefficients

This library contains the following files:

> main.cpp – includes:

>> file I/O operations

>>> get_file_uchar – reads any type of file into a returned unsigned char array. Used for loading both cover WAV files and files to be embedded

>>> set_stego_file – writes unsigned char buffer into a file on disk. Used for saving stego files. The following naming scheme is employed: if the original file name was e. g. song.wav, and the chosen level offset was 42, then the name of the stego file is song_stego42.wav

>>> set_retrieved_file – similar to the preceding function, but this one is for saving retrieved data files to disk. Adds “_retrieved” string before the filetype subfix

>> main – call this function to launch the whole algorithm (insertion and retrieval at once). It requires two files as parameters: cover audio (argv[1]) and secret data (argv[2]). This function manages file operations, parsing from / to arrays, console outputs and calls important get_available_space, LWT_LSBinsertion and LWT_LSBretrieval functions.

> embed.cpp – it includes the following functions:

>> unsigned get_available_space(int32_t* samples, unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample) – returns the accurate number of available bytes for insertion

>> unsigned LWT_LSBinsertion (int32_t* samples,
unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample, unsigned char *messageBuffer, unsigned buffer_length) – function performing the designed and explained dynamically adjusting insertion of *messageBuffer into *samples. Returns index of the last changed sample

>> int LWT_LSBretrieval(int32_t* samples, unsigned NumSamples, int NumChannels, int SampleRate, unsigned BitsPerSample, unsigned char *messageBuffer) – function retrieving dynamically inserted data from LWT coefficients. Returns number of retrieved bytes of secret data.  

>> void forward_LWT (int32_t* LWTcoefs, long int pos, int NumChannels) – computes asymetrical 2-level LWT (the first HP band is not being further decomposed). LWTcoefs is the vector of all coefficients, pos gives the index of first sample to be decomposed. It is designed to process 8 samples, i.e. two quadruples. These quadruples are ordered either subsequentially (for mono) or interleaved (for stereo). Decomposed coefficients are stored in-place, therefore no variable is returned.

>> void backward_LWT (int32_t* LWTcoefs, long int pos, int NumChannels) – see description of forward_LWT

>> void setEmbeddingParams(int *samples, int NumChannels, int cnt, int sampleRate, int bitsPerSample, int index, int* d1, int* d2) – important function performing K-filtering, computiation of mean-square and level mapping to the scopes of insertion (*d1 and *d2). cnt is the number of samples to be taken into account

>> int32_t staticCompressor (int32_t input, unsigned BitsPerSample) – function compressing peak of a single sample. Compression parameters are defined at embed.h

>> And supporting functions

> embed.h – adjustment of various parameters, such as compressor constants, level offset and mapping tables
> files of DspFilters library (https://github.com/vinniefalco/DSPFilters) used for K-Filtering
