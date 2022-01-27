/*
	Cole L - 5th January 2022 - https://github.com/cole8888/Gene2Pic
*/

#ifndef GENE2PIC_H
#define GENE2PIC_H

#define _GNU_SOURCE
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "LODEPNG/lodepng.h"
#include "NearestNeighbourUpscale.h"

#define DEFAULT_FILENAME "GenePic"
#define FILENAME_BUFFER_SIZE 255
#define PATH_SEPERATOR "/"
#define CHANNELS_PER_PIXEL_RGB 3

// Corresponding pixel values for each base colour. If you want, change these to the RGB values you want to use
#define CYTOSINE_COLOUR {6,   201, 150}
#define GUANINE_COLOUR  {17,  138, 178}
#define ADENINE_COLOUR  {239, 71,  111}
#define THYMINE_COLOUR  {255, 209, 102}

// Hard coded arguments. If you don't want to pass command line arguments for some reason.
// Cannot just specify one and collect the other from the commandline, must indicate all of them here.
#define USE_HARDCODED_ARGS false
#define INPUT_FILE_HARDCODED "inputSequence.txt"
#define SERPENTINE_HARDCODED false;
#define SCALE_HARDCODED 1

// Quickly find the length of the input file. This may not actually be the gene sequence length since
// characters like newlines or letters that are not a,t,c,g,u (upper and lower case) will be ignored.
long long int getFileLen(FILE *f);

// Finds the smallest dimensions for a square which can fit all the bases. 
// Image will have a section with black pixels at end if length is not a perfect square.
long long int findSquareSize(long long int len);

// Calculate the amount of time in seconds between the provided start and finish timespec structs.
double getElapsedTime(struct timespec start, struct timespec finish);

// Functions to determine whether or not a string starts with or ends with a particular string.
bool endsWith(char* str, char* toCheck);
bool startsWith(char* str, char* toCheck);

// Construct a filepath from two strings. sizeFile is the number of valid characters in the file string.
char *path_join(char* dir, char* file, int sizeFile);

// Determine how many digits are in a base 10 integer.
int getIntDigits(int num);

// Save the image, do not overwrite any previous images.
void saveImg(u_char* colours, long long int dim);

// Performs nearest neighbor upscaling of the original image where each pixel in the original image is expanded into an expanded pixel of size scale^2 in the upscaled image.
// Used for 24bit RBG images.
void upscaleNN_RGB(u_char *originalImg, u_char *upscaledImg, int dimX, int dimY, int scale);

// Assign each base in the sequence a coloured pixel in the image.
void base2colour(const char* gene, long long int dim, long long int len, int scale, bool serpentineLastRowFlip);

/* Flips every other row so that instead of:
	1->2->3
	<------
	4->5->6
	<------
	7->8->9

	It does:
	1->2->3
		  |
	6<-5<-4
	|
	7->8->9
*/
bool applySerpentine(char *geneSequence, long long int dim, long long int len);

// Read in the data from the sequence file and ignore any characters that are not ATCGU (upper or lowercase). Also convert lowercase to uppercase.
long long int readAndValidateInput(char *geneSequence, FILE *geneFile, long long int len);

// Main function, responsible for parsing the commandline arguments, opening the text file then coordinating other functions.
int main(int argc, char* argv[]);

#endif