/*
	Cole Lightfoot - 5th January 2022 - https://github.com/cole8888/Gene2Pic

	Creates an image representation of a genetic sequence.
	Capital letters are the same value as their lowercase counterparts.
	Thymine (T) and Uracil (U) are treated the same.

	"long long int" is necessary in some places because some sequences (such as the human genome) can overwhelm
	"long int" when we are using them to index elements in the colours array (since it is multiplied by 3 for RGB).
*/

#include "gene2pic.h"

// Quickly find the length of the input file. This may not actually be the gene sequence length since
// characters like newlines or letters that are not a,t,c,g,u (upper and lower case) will be ignored.
long long int getFileLen(FILE *f){
	fseek(f, 0, SEEK_END);	// Seek to the end of the file.
	long long int size = ftell(f);	// Get the current file pointer.
	fseek(f, 0, SEEK_SET);	// Seek back to the start.
	return size;
}

// Finds the smallest dimensions for a square which can fit all the bases. 
// Image will have a section with black pixels at end if length is not a perfect square.
long long int findSquareSize(long long int len){
	double dim = sqrt(len);
	return (long long int)ceil(dim);
}

// Calculate the amount of time in seconds between the provided start and finish timespec structs.
double getElapsedTime(struct timespec start, struct timespec finish){
	return (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
}

// Functions to determine whether or not a string starts with or ends with a particular string.
bool endsWith(char *str, char *toCheck){
	int n = strlen(str);
	int cl = strlen(toCheck);
	if(n < cl){
		return false;
	}
	char *subbuff = (char *)malloc(cl * sizeof(char));
	memcpy(subbuff, &str[n - cl], cl);
	return (strcmp(subbuff, toCheck) == 0);
}
bool startsWith(char *str, char* toCheck){
	int n = strlen(str);
	int cl = strlen(toCheck);
	if(n < cl){
		return false;
	}
	char *subbuff = (char *)malloc(cl * sizeof(char));
	memcpy(subbuff, &str[0 - cl], cl);
	return (strcmp(subbuff, toCheck) == 0);
}

// Construct a filepath from two strings. sizeFile is the number of valid characters in the file string.
char *path_join(char* dir, char* file, int sizeFile){
	int size = strlen(dir) + sizeFile + 2;
	char *buf = (char*)malloc(size*sizeof(char));
	if(NULL == buf){
		return NULL;
	}
	
	strcpy(buf, dir);

	// If the path does not end with the path seperator then we need to add it.
	if(!endsWith(dir, PATH_SEPERATOR)){
		strcat(buf, PATH_SEPERATOR);
	}
	// If the file starts with the path seperator, ignore it.
	if(startsWith(file, PATH_SEPERATOR)){
		char *filecopy = strdup(file);
		if(NULL == filecopy){
			free(buf);
			return NULL;
		}
		strcat(buf, filecopy + 1);
		free(filecopy);
	}
	else{
		strcat(buf, file);
	}
	return buf;
}

// Determine how many digits are in a base 10 integer.
int getIntDigits(int num){
	int digits = 0;
	while(num != 0){
		num /= 10;
		++digits;
	}
	return digits;
}

// Save the image, do not overwrite any previous images.
void saveImg(u_char *img, long long int dim){
	char *cwd = get_current_dir_name();	// Current working directory.
	char ext[] = ".png";	// File extention. Must be .PNG
	char buff[FILENAME_BUFFER_SIZE];	// Buffer to store the filename while we figure out what the file should be called.
	char *file;	// Once we settle on a filename, this will hold it.

	// See if an image with this name exists. If so, append a number and keep incrementing that number until we find one where that name is not in use.
	sprintf(buff, "%s%s", DEFAULT_FILENAME, ext);
	if(access(path_join(cwd, buff, strlen(DEFAULT_FILENAME) + strlen(ext)), F_OK) == 0){
		// File with this filename already exists. Let's try another one.
		int num = 2;	// Start by trying to append the number 2.
		int digits = 1;	// 2 has one digit.
		sprintf(buff, "%s%d%s", DEFAULT_FILENAME, num, ext);
		
		// Increment the appended number until we find a name that is not in use.
		while(access(path_join(cwd, buff, digits + strlen(DEFAULT_FILENAME) + strlen(ext)), F_OK) == 0){
			num++;
			digits = getIntDigits(num);
			sprintf(buff, "%s%d%s", DEFAULT_FILENAME, num, ext);
		}
		
		// Save the filename with the number at the end.
		sprintf(buff, "%s%d%s", DEFAULT_FILENAME, num, ext);
		file = path_join(cwd, buff, digits + strlen(DEFAULT_FILENAME) + strlen(ext));
	}
	else{
		// File with this name does not exist, lets use this name for the new file.
		sprintf(buff, "%s%s", DEFAULT_FILENAME, ext);
		file = path_join(cwd, buff, strlen(DEFAULT_FILENAME) + strlen(ext));
	}

	// If the filename is NULL for some reason, let the user know, but still try to save saving the image. (Since this might have taken a long time to compute and there is no need to exit for this error.)
	if(NULL == file){
		fprintf(stderr, "\nError when making filename, it is NULL. Will attempt saving anyway using filename gene2pic_backupName.png\n");
		file = "gene2pic_backupName.png";
	}

	// Start the timer and then save the image.
	struct timespec start, finish;
	clock_gettime(CLOCK_MONOTONIC, &start);
	unsigned error = lodepng_encode24_file(file, img, dim, dim);
	clock_gettime(CLOCK_MONOTONIC, &finish);

	// See if there was an issue when saving the image.
	if(error){
		fprintf(stderr, "\nUnable to save the image, lodepng returned an error.\nError %u: %s\n", error, lodepng_error_text(error));
	}
	else{
		printf("Saved to %s (%f secs)\n\n", file, getElapsedTime(start, finish));
	}
}

// Assign each base in the sequence a coloured pixel in the image.
void base2colour(const char *gene, long long int dim, long long int len, int scale, bool serpentineLastRowFlip){
	printf("\nStart assigning bases to colours...\n");
	// Corresponding pixel values for each base colour. (You can change these in the header file.)
	u_char baseColour[4][3] = {
		CYTOSINE_COLOUR,
		GUANINE_COLOUR,
		ADENINE_COLOUR,
		THYMINE_COLOUR
	};

	// Hold all the colour values which will then be turned into an image. No need to zero, will be completely overwritten.
	u_char *img = (u_char *)malloc(dim * dim * (long long int)3 * sizeof(u_char));
	if(NULL == img){
		fprintf(stderr, "Unable to allocate img array... May have run out of RAM.\n");
		exit(EXIT_FAILURE);
	}

	// Time how long it takes to go through all the bases.
	struct timespec start, finish;
	clock_gettime(CLOCK_MONOTONIC, &start);

	// long long int diff = len % dim;	// Number of bases on the incomplete row of the image.

	// Branchlessly assign the proper colours to each base. Done using multiprocessing to speed it up.
	#pragma omp parallel for
	for(long long int i = 0; i < len; i++){
		// Since we already removed all invalid characters at the start, we only actually need to test for 3 base pairs.
		u_int8_t baseId = 0;	// Default to cytosine.
		// baseId += ((gene[i] == 'C') ? 0 : 0);	// Cytosine	(Commented because we only need to check for 3 base pairs since we know all bases are valid, speeds up colour assignment.)
		baseId += ((gene[i] == 'G') ? 1 : 0);		// Guanine
		baseId += ((gene[i] == 'A') ? 2 : 0);		// Adenine
		baseId += ((gene[i] == 'T') ? 3 : 0);		// Thymine
		
		// Assign the colour to this pixel.
		img[(long long int)3 * i] = baseColour[baseId][0];						// R
		img[(long long int)3 * i + (long long int)1] = baseColour[baseId][1];	// G
		img[(long long int)3 * i + (long long int)2] = baseColour[baseId][2];	// B
	}
	
	// Number of bases does not fit perfectly into the image, will be blank spots we need to zero and fill out the incomplete row.
	if(len < dim * dim){
		// Make sure the pixels we did not overwrite at the end are zeroed (Prevents garbage from memory from sneaking into the image). (Not enough for multithreading to make sense).
		for(long long int i = 0; i < (dim * dim) - len; i++){
			img[(long long int)3 * (len + i)] = 0;						// R
			img[(long long int)3 * (len + i) + (long long int)1] = 0;	// G
			img[(long long int)3 * (len + i) + (long long int)2] = 0;	// B
		}
	}

	// If serpentine mode was selected and we previously detected that the last row is partially completed and it needs to be flipped, then handle the flipping of the last row here.
	// (Cannot do it earlier in program because the geneSequence array does not have characters at these blank positions.)
	if(serpentineLastRowFlip){
		long long int filledRows = len / dim;	// Number of rows which have been completely filled previously. (Tells us which row is the partially completed one and it is not necessarily the last row in the image.)
		// Swap the pixels from one half of this row with the ones in the other half making sure to place them in the right places so the row is flipped.
		for(long long int j = 0; j < dim / (long long int)2; j++){
			// Copy the data for the pixel we are about to overwrite so that we can swap it later.
			u_char tmpR = img[(filledRows * dim + j) * (long long int)3];						// R
			u_char tmpG = img[(filledRows * dim + j) * (long long int)3 + (long long int)1];	// G
			u_char tmpB = img[(filledRows * dim + j) * (long long int)3 + (long long int)2];	// B
			
			// Overwrite the pixel we just made a temporary copy of with the one in the corresponding position in the other half of the row.
			img[(filledRows * dim + j) * (long long int)3] = img[(filledRows * dim + dim - (long long int)1 - j) * 3];														// R
			img[(filledRows * dim + j) * (long long int)3 + (long long int)1] = img[(filledRows * dim + dim - (long long int)1 - j) * (long long int)3 + (long long int)1];	// G
			img[(filledRows * dim + j) * (long long int)3 + (long long int)2] = img[(filledRows * dim + dim - (long long int)1 - j) * (long long int)3 + (long long int)2];	// B
			
			// Overwrite the pixel we just copied from with the temporary copy of the pixel it overwrote.
			img[(filledRows * dim + dim - (long long int)1 - j) * (long long int)3] = tmpR;						// R
			img[(filledRows * dim + dim - (long long int)1 - j) * (long long int)3 + (long long int)1] = tmpG;	// G
			img[(filledRows * dim + dim - (long long int)1 - j) * (long long int)3 + (long long int)2] = tmpB;	// B
		}
	}

	// Stop the clock, we finished assigning colours to bases.
	clock_gettime(CLOCK_MONOTONIC, &finish);
	printf("Finished assigning colours to bases.\t(%f secs)\n", getElapsedTime(start, finish));

	// See if we should upscale the image.
	if(scale > 1){
		// We want to upscale the image.
		printf("\nStart upscaling the image...\n");
		
		// Allocate memory for the upscaled image.
		long long int scaledDim = dim * (long long int)scale;	// Dimmension of the upscaled image.
		u_char *upscaledImg = (u_char *)malloc(scaledDim * scaledDim * (long long int)3 * sizeof(u_char));
		if(NULL == upscaledImg){
			fprintf(stderr, "Unable to allocate upscaledImg array... May have run out of RAM.\n");
			exit(EXIT_FAILURE);
		}
		
		clock_gettime(CLOCK_MONOTONIC, &start);	// Start the timer.
		
		upscaleNN_RGB(img, upscaledImg, dim, dim, scale);	// Upscale the original image.
		
		clock_gettime(CLOCK_MONOTONIC, &finish);	// Stop the timer.
		printf("Finished upscaling the image.\t\t(%f secs)\n", getElapsedTime(start, finish));

		free(img);	// Free the original unscaled image.
		
		printf("\nStart saving the image...\n");
		saveImg(upscaledImg, scaledDim);	// Save the array as an image.
		free(upscaledImg);	// Free the upscaled image.
	}
	else{
		// We do not want to upscale the image. Save the 1:1 image.
		printf("\nStart saving the image...\n");
		saveImg(img, dim);	// Save the array as an image.
		free(img);	// Free the image.
	}
}

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

	This may be a more natural way to view the sequence since if you were to lay out the DNA strand using the
	non-serpentine method it would require you to snip the strand every time you want to go to a new row in the image.
	This way if you were to draw it out on paper, you could traverse the whole sequence without ever taking your pencil
	off the paper. 

	Unfortunately I need to wait until later in the program to flip the incomplete row, since the gene sequence char array
	does not have a character I can use to indicate it is blank. Adding one would get rid of the sneaky trick which lets
	me only check 3 bases instead of 4 in base2colour().

	Returns a boolean which indicates whether or not base2colour() needs to flip the incomplete row.
*/
bool applySerpentine(char *geneSequence, long long int dim, long long int len){
	printf("\nStart applying serpentine pattern to sequence...\n");
	long long int filledRows = len/dim;	// Number of rows which can be completely filled. Excludes the incomplete row near the end if it exists, that is handled in base2colour().

	// Start the timer.
	struct timespec start, finish;
	clock_gettime(CLOCK_MONOTONIC, &start);

	// For every second row, flip it so the bases at the start are the ones at the end and vice versa. Use multi-processing because why not.
	#pragma omp parallel for
	for(long long int i = 1; i < filledRows; i+=2){
		for(long long int j = 0; j < dim / (long long int)2; j++){
			char tmp = geneSequence[i * dim + j];
			geneSequence[i * dim + j] = geneSequence[i*dim + dim - (long long int)1 - j];
			geneSequence[i * dim + dim - (long long int)1 - j] = tmp;
		}
	}

	// Stop the timer.
	clock_gettime(CLOCK_MONOTONIC, &finish);
	printf("Finished applying serpentine.\t\t(%f secs)\n", getElapsedTime(start, finish));
	
	// Figure out if we need to ask base2colour() to flip the incomplete row (if it exists).
	if((filledRows)%2 != 0){
		// Incomplete row is on an uneven row index.
		if(len%dim != 0){
			// Incomplete row is actually incomplete and does have blank spots.
			return true;
		}
	}
	return false;
}

// Read in the data from the sequence file and ignore any characters that are not ATCGU (upper or lowercase). Also convert lowercase to uppercase.
long long int readAndValidateInput(char *geneSequence, FILE *geneFile, long long int len){
	printf("Start validation of input sequence...\n");
	// Length of the valid gene sequence.
	long long int validBaseCount = 0;

	// Initialize and start the timer.
	struct timespec start, finish;
	clock_gettime(CLOCK_MONOTONIC, &start);
	
	// Read the entire sequence all at once into the geneSequence array. (Slightly faster than reading one letter at a time.) 
	fread(geneSequence, sizeof(char), len, geneFile);

	// Remove any invalid characters from the sequence.
	// Cannot be multithreaded without using mutex when writing to the geneSequence array and incrementing validBaseCount, which slows it down even more.
	// So I'm leaving this single threaded since it is already pretty damn fast.
	for(long long int i = 0; i < len; i++){
		// Make a copy of the letter and change it to it's ordinal form so we can do math on it.
		u_int8_t normLetter = (u_int8_t)geneSequence[i];
		
		// Branchlessly convert the letter to uppercase.
		normLetter -= 32 * (normLetter >= 'a' && normLetter <= 'z');

		// Treat Uracil (U) as if it is Thymine (T) by subracting 1 from the ordinal value. (84 = T, 85 = U)
		normLetter -= normLetter == 'U' ? 1 : 0;

		// Determine if the current letter is a valid character. (This commented one is slower despite the fact it stops early when the letter is found. Intuitively I'd think this would have been faster.)
		// normLetter == 'C' ? geneSequence[validBaseCount++] = (char)normLetter : normLetter == 'G' ? geneSequence[validBaseCount++] = (char)normLetter :  normLetter == 'A' ? geneSequence[validBaseCount++] = (char)normLetter : normLetter == 'T' ? geneSequence[validBaseCount++] = (char)normLetter : 0;
		u_int8_t isValid = 0;
		isValid += normLetter == 'C' ? 1 : 0;
		isValid += normLetter == 'G' ? 1 : 0;
		isValid += normLetter == 'A' ? 1 : 0;
		isValid += normLetter == 'T' ? 1 : 0;
		isValid ? geneSequence[validBaseCount++] = (char)normLetter : 0; // If base is valid, then add it to the sequence and then increment the length counter.
	}

	// Stop the timer and figure out how long it took to read in and validate all the bases.
	clock_gettime(CLOCK_MONOTONIC, &finish);
	printf("Valid input sequence is %lld bases.\t(%f secs)\n", validBaseCount, getElapsedTime(start, finish));

	// Close the input file, we are done with it now.
	fclose(geneFile);
	
	return validBaseCount;
}

// Main function, responsible for parsing the commandline arguments, opening the text file then coordinating other functions.
int main(int argc, char *argv[]){
	int scale = SCALE_HARDCODED;
	bool serpentine = SERPENTINE_HARDCODED;
	char *inputFile = INPUT_FILE_HARDCODED;

	// See if we should check the commandline arguments or use hardcoded ones instead.
	if(!USE_HARDCODED_ARGS){
		if(argc < 2 || argc > 4){
			fprintf(stderr, "Incorrect number of arguments!\nAvailable usage modes:\n./gene2pic <INPUT_FILE>\n./gene2pic <INPUT_FILE> <SCALE>\n./gene2pic <INPUT_FILE> <SERPENTINE>\n./gene2pic <INPUT_FILE> <SERPENTINE> <SCALE>\n");
			return EXIT_FAILURE;
		}
		else if(argc == 3){
			// Only 2 arguments provided. Input file is first one and scale or serpentine is second. Find out which it is.

			char *temp;
			long value = strtol(argv[2], &temp, 10);
			if(temp != argv[2] && *temp == '\0'){
				// Argument is a number, so it is specifying a scale.
				scale = (int)value;
			}
			else if(strcmp(argv[2], "serpentine") == 0 || strcmp(argv[2], "SERPENTINE") == 0){
				// Argument matches the serpentine string, so enable serpentine mode.
				serpentine = true;
			}
			else{
				// Argument matches neither scale or serpentine, inform user.
				fprintf(stderr, "Invalid data in last argument.\nFor serpentine it must be \"serpentine\" or be blank, and for scale it must be a non-zero integer.\nUsage: ./gene2pic <INPUT_FILE> <SERPENTINE> <SCALE>\n");
				return EXIT_FAILURE;
			}
		}
		else if(argc == 4){
			// 3 arguments provided. Make sure they are in the right order and have acceptable values.

			if(strcmp(argv[2], "serpentine") == 0 || strcmp(argv[2], "SERPENTINE") == 0){
				// Argument matches the serpentine string, so enable serpentine mode.
				serpentine = true;
			}
			else{
				fprintf(stderr, "Invalid data for serpentine argument. Must be \"serpentine\" or left empty.\nUsage: ./gene2pic <INPUT_FILE> <SERPENTINE> <SCALE>\n");
				return EXIT_FAILURE;
			}

			char *temp;
			long value = strtol(argv[3], &temp, 10);
			if(temp != argv[3] && *temp == '\0'){
				// Argument is a number, so it is specifying a scale.
				scale = (int)value;
			}
			else{
				fprintf(stderr, "Invalid data in scale argument. Must be a non-zero integer.\nUsage: ./gene2pic <INPUT_FILE> <SERPENTINE> <SCALE>\n");
				return EXIT_FAILURE;
			}
		}

		// Make sure scale is not an invalid value.
		if(scale < 1){
			fprintf(stderr, "Invalid data in scale argument. Must be a non-zero integer.\nUsage: ./gene2pic <INPUT_FILE> <SERPENTINE> <SCALE>\n");
			return EXIT_FAILURE;
		}
		inputFile = argv[1];	// Read the filename from the user input.
	}

	// Start timer to see how long the whole program takes.
	struct timespec start, finish;
	clock_gettime(CLOCK_MONOTONIC, &start);

	FILE *geneFile = fopen(inputFile, "r");	// Get the file, open with read permissions.
	if(geneFile == (FILE *) NULL){
		// Error when opening the file or the file was not found.
		fprintf(stderr,"File %s not found!\n", inputFile);
		return EXIT_FAILURE;
	}

	// Find the length of the file (including invalid characters).
	long long int len = getFileLen(geneFile);
	printf("Input file is %lld characters.\n\n", len);

	// Place to hold the sequence in memory.
	char *geneSequence = (char *)malloc(len * sizeof(char));
	if(NULL == geneSequence){
		fprintf(stderr,"Unable to allocate geneSequence array. May have run out of RAM.");
		return EXIT_FAILURE;
	}

	// Remove any invalid characters from the input sequence and determine how many valid bases there are.
	long long int validBaseCount = readAndValidateInput(geneSequence, geneFile, len);
	if(validBaseCount < 1){
		// No valid bases.
		fprintf(stderr, "Input file has 0 valid characters... Exiting.\n");
		return EXIT_FAILURE;
	}

	// Find the optimal sized square dimmensions which can fit the sequence with the least amount of blank pixels as possible.
	long long int dim = findSquareSize(validBaseCount);

	// If we want to represent the sequence using a serpentine pattern.
	bool serpentineLastRowFlip = false;	// Flag to indicate if we need to flip the last row in base2colour.
	if(serpentine){
		serpentineLastRowFlip = applySerpentine(geneSequence, dim, validBaseCount);
	}

	// Start assigning colours to bases, upscales the image (if wanted), and then sends the finished array to saveImg().
	base2colour(geneSequence, dim, validBaseCount, scale, serpentineLastRowFlip);

	// Stop the timer.
	clock_gettime(CLOCK_MONOTONIC, &finish);
	printf("DONE. Took %f seconds.\n", getElapsedTime(start, finish));

	free(geneSequence);	// Free the sequence character array.
	return EXIT_SUCCESS;
}