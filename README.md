# Gene2Pic
Represents a genetic sequence as an image.

The Python and C program both do the same thing, except the Python version has more functionality and the C version is much faster and tends to produce smaller image file sizes. I'd reccomend using the C version because it is more efficient at encoding the images.



To run the Python script: `python3 gene_to_pic_multi_test.py`
Arguments:
- Read the sequence from a file: `--file X` (X is the name of the file)
- Scale up the image: `--scale X` (X is a positive integer)
- Manually set the number of threads: `--threads X` (X is a positive integer)
- Do not optimize the image: `-- no_optimize` (saves RAM and is faster, but larger images)
- Change the base colour(s): `--A FFFFFF --T FFFFFF --C FFFFFF --G FFFFFF` (FFFFFF is a valid colour hex code, only need to specify the ones you want to change)

To run the C program:
- Place the lodepng.c and lodepng.h files into the same directory as gene2pic.c
- Compile: `gcc gene2pic.c lodepng.c -Wall -Wextra -fopenmp -lm -o gene2pic`
- Run: `./gene2pic X` (X is the name of the file)
