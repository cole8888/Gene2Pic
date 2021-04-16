// Creates an image from a genetic sequence.
// Capital letters are the same value as their lowercase counterparts.
// Thymine (T) and Uracil (U) are treated the same.

// "long long int" is necesary in some places because some sequences (such as the human genome) can overwhelm
// "long int" when we are using them to index elements in the colours array (since it is multiplied by 3).

// Cole Lightfoot - 16th April 2021 - https://github.com/cole8888/Gene2Pic

#define _GNU_SOURCE
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "lodepng.h"

// Quickly find the length of the input file. This may not actually be the gene sequence length since
// characters like newlines or letters that are not a,t,c,g,u (upper and lower case) will be ignored.
long long int getFileLen(FILE* f){
    // Seek to the end of the file.
    fseek(f, 0, SEEK_END);
    // Get the current file pointer.
    long long int size = ftell(f);
    // Seek back to the end.
    fseek(f, 0, SEEK_SET);
    return size;
}

// Finds the smallest sized square which can fit all the bases. 
// Image will have blank spots at end if length is not a perfect square.
long long int findSquareSize(long long int len){
    double dim = sqrt(len);
    return (long long int)ceil(dim);
}

// Functions to determine whether or not a string starts with or ends with a particular string.
bool endsWith(char* str, char* toCheck){
    int n = strlen(str);
    int cl = strlen(toCheck);
    if(n < cl){
       return false;
    }
    char* subbuff = (char*)malloc(cl * sizeof(char));
    memcpy(subbuff, &str[n-cl], cl);
    return (strcmp(subbuff, toCheck) == 0);
}
bool startsWith(char* str, char* toCheck){
    int n = strlen(str);
    int cl = strlen(toCheck);
    if(n < cl){
       return false;
    }
    char* subbuff = (char*)malloc(cl * sizeof(char));
    memcpy(subbuff, &str[0-cl], cl);
    return (strcmp(subbuff, toCheck) == 0);
}

// Construct a filepath from two strings.
// sizeFile is the number of valid characters in the file string.
char * path_join(char* dir, char* file, int sizeFile){
    char* path_sep = "/";
    int size = strlen(dir) + sizeFile + 2;
    char* buf = (char*)malloc(size*sizeof(char));
    if(NULL == buf){
        return NULL;
    }
    strcpy(buf, dir);
    if(!endsWith(dir, path_sep)){
        strcat(buf, path_sep);
    }
    if(startsWith(file, path_sep)){
        char *filecopy = strdup(file);
        if(NULL == filecopy){
            free(buf);
            return NULL;
        }
        strcat(buf, ++filecopy);
        free(--filecopy);
    }
    else{
        strcat(buf, file);
    }
    return buf;
}

// Determine how many digits are in an integer.
int getIntDigits(int num){
    int digits = 0;
    while(num != 0){
        num /= 10;
        ++digits;
    }
    return digits;
}

// Save the image, do not overwrite any previous images.
void saveImg(u_char* colours, long long int dim, long long int len){
    char* cwd = get_current_dir_name();
    char fileName[] = "GenePic";
    char ext[] = ".png";
    char* file;
    char buff[100];
    // See if an image with this name exists. If so, append a number and keep incrementing that number until
    // we find one where that name is not in use.
    sprintf(buff, "%s%s", fileName, ext);
    if(access(path_join(cwd, buff, strlen(fileName)+strlen(ext)), F_OK) == 0){
        int num = 2;
        int digits = 1;
        sprintf(buff, "%s%d%s", fileName, num, ext);
        // Increment the appended number until we find a name that is not in use.
        while(access(path_join(cwd, buff, digits+strlen(fileName)+strlen(ext)), F_OK) == 0){
            num += 1;
            digits = getIntDigits(num);
            sprintf(buff, "%s%d%s", fileName, num, ext);
        }
        sprintf(buff, "%s%d%s", fileName, num, ext);
        file = path_join(cwd, buff, digits+strlen(fileName)+strlen(ext));
    }
    else{
        sprintf(buff, "%s%s", fileName, ext);
        file = path_join(cwd, buff, strlen(fileName)+strlen(ext));
    }
    // If the filename is NULL for some reason, let the user know and don't attempt saving the image.
    if(file == NULL){
        fprintf(stderr, "Unable to allocate memory when using path_join(). (Or some other problem with that function)");
    }
    else{
        // Save the image.
        clock_t begin = clock();
        unsigned error = lodepng_encode24_file(file, colours, dim, dim);
        clock_t end = clock();
        if(error){
            fprintf(stderr, "Error %u: %s\n", error, lodepng_error_text(error));
        }
        else{
            printf("Image with %lld bases saved to %s (%f secs)\n", len, file, (double)(end - begin) / CLOCKS_PER_SEC);
        }
    }
}

// Assign each pixel in the image a colour corresponding to the base it represents.
void colour2base(const char* gene, long long int dim, long long int len){
    // Corresponding pixel values for each base colour.
    // If you want, change these to the RGB values you want to use
    u_char baseColour[4][3] = {
    //   R    G    B
        {6,   201, 150}, // Cytosine
        {17,  138, 178}, // Guanine
        {239, 71,  111}, // Adenine
        {255, 209, 102}, // Thymine
    };

    // How many bases will be in the last row.
    long long int diff = len%dim;

    // Hold all the colour values which will then be turned into an image. No need to zero, will be overwritten.
    long long int colours_size = dim * dim * (long long int)3;
    u_char* colours = (u_char*)malloc(colours_size*sizeof(u_char));
    if(NULL == colours){
        fprintf(stderr,"Unable to allocate colours array... May have run out of RAM.");
        exit(1);
    }

    // Time how long it takes to go through all the bases.
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    // Branchlessly assign the proper colours to each base done using multiprocessing to speed it up.
    // This loop is for all rows except the last one. 
    #pragma omp parallel for
    for(long long int i = 0; i<(len-diff); i++){
        int baseId = 0;
        baseId += ((gene[i] == 'C') ? 1 : 0);     // Cytosine
        baseId += ((gene[i] == 'G') ? 2 : 0);     // Guanine
        baseId += ((gene[i] == 'A') ? 3 : 0);     // Adenine
        baseId += ((gene[i] == 'T') ? 4 : 0);     // Thymine
            
        colours[(long long int)3*i] = baseColour[baseId-1][0];                    // R
        colours[(long long int)3*i+(long long int)1] = baseColour[baseId-1][1];   // G
        colours[(long long int)3*i+(long long int)2] = baseColour[baseId-1][2];   // B
    }
    // For all the bases to be placed on the last row.
    for(long long int i = 0; i<diff; i++){
        int baseId = 0;
        baseId += (gene[len-diff+i] == 'C') ? 1 : 0;    // Cytosine
        baseId += (gene[len-diff+i] == 'G') ? 2 : 0;    // Guanine
        baseId += (gene[len-diff+i] == 'A') ? 3 : 0;    // Adenine
        baseId += (gene[len-diff+i] == 'T') ? 4 : 0;    // Thymine
        
        colours[(long long int)3*(len-diff) + (long long int)3*i] = baseColour[baseId-1][0];                     // R
        colours[(long long int)3*(len-diff) + (long long int)3*i + (long long int)1] = baseColour[baseId-1][1];  // G
        colours[(long long int)3*(len-diff) + (long long int)3*i + (long long int)2] = baseColour[baseId-1][2];  // B
    }
    // Make sure the elements we did not overwrite at the end of the last row are zeroed.
    for(long long int i = 0; i<(dim-diff); i++){
        colours[(long long int)3*(len-diff) + (long long int)3*(i+diff)] = 0;                    // R
        colours[(long long int)3*(len-diff) + (long long int)3*(i+diff) + (long long int)1] = 0; // G
        colours[(long long int)3*(len-diff) + (long long int)3*(i+diff) + (long long int)2] = 0; // B
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    double elapsed = (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("All bases assigned to colours, saving image... (%f secs)\n", elapsed);

    // Call the function to save the array as an image.
    saveImg(colours, dim, len);

    free(colours);
}

// Main function, responsible for opening the text file, retrieving the sequence and parsing it and storing it to prepare it for processing.
// Also handles the conversion to uppercase and disregarding invalid characters.
int main(int argc, char* argv[]){
    // Check arguments number, must be 2.
    if(argc != 2){
		fprintf(stderr, "Filename not provided in commandline arguments!\n");
		return EXIT_FAILURE;
	}
    // Get the file, open with read permissions
    FILE *geneFile = fopen(argv[1], "r");
    // Make sure the file was opened okay.
    if(geneFile == (FILE*) NULL){
		fprintf(stderr,"File %s not found\n", argv[1]);
		return EXIT_FAILURE;
	}
    // Find the length of the file (including invalid characters)
    long long int len = getFileLen(geneFile);
    printf("Input file is %lld characters.\n", len);

    // Holds the gene sequence in memory.
    char* geneSequence = (char*)malloc(len * sizeof(char));
    if(NULL == geneSequence){
        fprintf(stderr,"Unable to allocate geneSequence array... May have run out of RAM.");
        return EXIT_FAILURE;
    }

    // Length of the valid gene sequence.
    long long int len2 = 0;

    // Copy the valid characters into memory and determine the number of valid characters.
    // for all characters in the input file.
    // Unfortunately needs to be single threaded due to the nature of the operation.
    clock_t begin = clock();
    for(long long int i = 0; i < len; i++){
        // Take input one character at a time 
        char c = fgetc(geneFile);
        // Make a copy of the letter and change it to it's ordinal form so we can do math on it.
        int normLetter = (int)c;
        // Branchlessly convert the letter to uppercase.
        normLetter -= 32 * (normLetter >= 'a' && normLetter <= 'z');

        // Treat Uracil (U) as if it is Thymine (T) by subracting 1 from the ordinal value. (84 = T, 85 = U)
        normLetter -= normLetter == 'U' ? 1 : 0;

        // Determine if the current letter is a valid character.
        int isValid = 0;
        isValid += normLetter == 'C' ? 1 : 0;
        isValid += normLetter == 'G' ? 1 : 0;
        isValid += normLetter == 'A' ? 1 : 0;
        isValid += normLetter == 'T' ? 1 : 0;
        if(isValid == 1){
            // If so, then add it to the sequence and then increment the length counter.
            geneSequence[len2] = (char)normLetter;
            len2++;
        }
    }
    clock_t end = clock();

    // Close the input file, we are done with it now.
    fclose(geneFile);

    printf("Valid Gene sequence is %lld bases. (%f secs)\n", len2, (double)(end - begin) / CLOCKS_PER_SEC);

    // Check to make sure valid sequence is not less than 1 base.
    if(len2 < 1){
        fprintf(stderr, "Input file has 0 valid characters. Exiting.\n");
        return EXIT_FAILURE;
    }

    // Find the optimal sized square which can fit the gene sequence with the least amount of blank pixels as possible
    long long int dim = findSquareSize(len2);

    // Start assigning colours to bases and then send the finished array to the saveImg function.
    colour2base(geneSequence, dim, len2);
    return EXIT_SUCCESS;
}