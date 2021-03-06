# Gene2Pic
Represents a genetic sequence as an image.

Inspired by [this](https://www.reddit.com/r/dataisbeautiful/comments/mg1cxr/oc_entire_genome_of_covid_virus_sarscov2/) post. I didn't see the author provide any source code, so I made my own programs that could handle any genetic sequence.

The Python and C programs do almost the same thing. The C version is much faster, has serpentine mode, and tends to produce smaller image file sizes. I'd reccomend using the C version if you can because it is more efficient at encoding the images. Both programs automatically treat Uracil (U) as Thymine (T).

<hr>

To run the C program:
- Clone this repository: `git clone https://github.com/cole8888/Gene2Pic`
- Move into the folder: `cd Gene2Pic`
- Compile: `make`
- Run: `./gene2pic <INPUT_FILE>` (Where \<INPUT_FILE\> is the name of the file)

Optional arguments:
- Upscale the final image: `./gene2pic <INPUT_FILE> <SCALE>` (Where \<SCALE\> is a positive integer)
- Flip every second row: `./gene2pic <INPUT_FILE> <SERPENTINE>` (Where \<SERPENTINE\> is "serpentine" without the quotes)
- Upscale and flip every second row: `./gene2pic <INPUT_FILE> <SERPENTINE> <SCALE>`

On a Ryzen 3700X it is able to go through the entire Human genome in less than 24 seconds. Most of that time is spent reading from the disk and making sure that only valid characters are stored in memory. It also takes fairly long for lodepng to save such a huge image.

![Image](https://github.com/cole8888/Gene2Pic/blob/main/C_Example.png)

<hr>

To run the Python script: `python3 gene2pic.py`

Arguments:
- Read the sequence from a file: `--file X` (X is the name of the file)
- Scale up the image: `--scale X` (X is a positive integer)
- Manually set the number of threads: `--threads X` (X is a positive integer)
- Do not optimize the image: `--no_optimize` (saves RAM and is faster, but larger images)
- Change the base colour(s): `--A FFFFFF --T FFFFFF --C FFFFFF --G FFFFFF` (where FFFFFF is a valid colour hex code, only need to specify the ones you want to change)

![Image](https://github.com/cole8888/Gene2Pic/blob/main/Python_Example.png)

<hr>

I've provided several example genetic sequences as well as their expected outputs you can try out if you'd like. You can find additional genetic sequences at https://www.ncbi.nlm.nih.gov/genome/.
Before using a sequence, open it and make sure that there are no headers in the data. You can remove these using find and replace with regex.

Here's some examples:

Covid-19 (Scaled up 4X)

![Image](https://github.com/cole8888/Gene2Pic/blob/main/Example%20Images/Covid-19_scale4X.png)

Nanoarchaeum equitaans

![Image](https://github.com/cole8888/Gene2Pic/blob/main/Example%20Images/Nanoarchaeum%20equitaans.png)

Tremblaya

![Image](https://github.com/cole8888/Gene2Pic/blob/main/Example%20Images/Tremblaya.png)

If you want to scale up an image from the program after the fact in an image editor, you need to use "no interpolation" or "nearest neighbor". Also, you must save the images in a format which supports a lossless compression format such as with PNGs since any lossy compression formats like JPEG will actually make the image filesize larger than it would be otherwise and mess up the colours.

The C version uses lodepng which can be found here https://github.com/lvandeve/lodepng
