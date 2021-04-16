# Creates an image from a genetic sequence.
# Capital letters are the same value as their lowercase counterparts.
# Thymine (T) and Uracil (U) are treated the same.

# Cole Lightfoot - 30th March 2021

import sys
import math
import os.path
import argparse
import subprocess
import ctypes as c
import numpy as np
from PIL import Image
from PIL import ImageColor
import multiprocessing as mp

# Command line arguments.
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "--file",
    type=str,
    default=None,
    help='Directory of images you want to use. This should be a single level directory.',
)
parser.add_argument(
    "--scale",
    type=int,
    default=1,
    help="How much to multiply the scale of the image by. 1 = 1px per base, 2=4px, 4=8px, etc...",
)
parser.add_argument(
    "--threads",
    type=int,
    default=mp.cpu_count(),
    help="Number of threads to use. Automatically uses as many as your cpu has.",
)
parser.add_argument(
    "--A",
    type=str,
    default=None,
    help='Hex colour to use for Adenine',
)
parser.add_argument(
    "--T",
    type=str,
    default=None,
    help='Hex colour to use for Thymine',
)
parser.add_argument(
    "--C",
    type=str,
    default=None,
    help='Hex colour to use for Cytosine',
)
parser.add_argument(
    "--G",
    type=str,
    default=None,
    help='Hex colour to use for Guanine',
)
parser.add_argument('--no_optimize', action="store_true", help='Do not compress and optimize the final png, this will save time and RAM.')
args = parser.parse_args()

# Catch possible issues with arguments
if(args.scale <= 0):
    raise ValueError("Scale multiplier cannot be 0 or less.")
if(args.threads < 1):
    raise ValueError("Invalid number of threads. Must be greater than zero.")
if(args.threads > mp.cpu_count()):
    if(input("\nWarning, you are going to use " + str(args.threads) + " threads. Your CPU has " + str(mp.cpu_count()) + " threads.\nThis may cause lower performance than if you used " + str(mp.cpu_count()) + " threads.\nDo you want to continue? (y/N): ").lower() != "y"):
        exit()

# Finds the smallest sized square which can fit all the basepairs
def findSquareSize(gene):
    dim = math.sqrt(len(gene))
    # If length of sequence is not a perfect square, add 1. Image will have blank spots at end.
    if(dim%1 != 0):
        dim += 1
    return int(dim)

# Get the input sequence, we need to switch the terminal mode to allow
# inputs greater than 4095 bases. We get the sequence and switch it back.
def getGeneFromCLI():
    subprocess.check_call(["stty","-icanon"]) # Comment me if input errors happen.
    inputStr = input("Input the sequence: ")
    subprocess.check_call(["stty","icanon"]) # Comment me if input errors happen.
    return inputStr

# Get the contents of a file and remove newline characters, then return as a string.
def getGeneFromFile(f):
    try:
        with open(f, 'r') as file:
            return file.read().replace('\n', '')
    except:
        raise FileNotFoundError("Could not find / open the file.", f)

# Convert a hex value into a RGB value. Used for custom colours.
def hex2RGB(colour):
    # Check if user did not add "#" to start, add it if not.
    if(colour[0] != "#"):
        colour = "#" + colour
    if(len(colour) != 7):
        raise ValueError(colour + " Is not a valid hex code for a colour. Incorrect number of characters.")
    for i in colour:
        if(not(ord(i) > 47 and ord(i) < 58 or ord(i) > 64 and ord(i) < 71 or ord(i) > 97 and ord(i) < 103 or ord(i) == 35)):
            raise ValueError(colour + " Is not a valid hex code for a colour. Invalid character: " + str(i))
    return ImageColor.getcolor(colour, "RGB")

# Save the image and apply scaling
def saveImg(array, dim):
    # Generate the image.
    out = Image.fromarray(array, mode="RGB")

    # Scale up the image
    if(args.scale != 1):
        out = out.resize((dim*args.scale, dim*args.scale), resample=Image.NEAREST)

    # Save the image, do not overwrite any previous images.
    path = os.getcwd()
    fileName = "GenePic"
    ext = ".png"
    if(os.path.isfile(os.path.join(path, fileName + ext))):
        num = 2
        while(os.path.isfile(os.path.join(path, fileName + str(num) + ext))):
            num += 1
        if(not args.no_optimize):
            out.save(os.path.join(path, fileName + str(num) + ext), optimize = True, compress_level = 9)
        else:
            out.save(os.path.join(path, fileName + str(num) + ext))
        print("Image with " + str(len(sequence)) + " bases saved to " + os.path.join(path, fileName + str(num) + ext))
    else:
        if(not args.no_optimize):
            out.save(os.path.join(path, fileName + ext), optimize = True, compress_level = 9)
        else:
            out.save(os.path.join(path, fileName + ext))
        print("Image with " + str(len(sequence)) + " bases saved to " + os.path.join(path, fileName + ext))

# Associate the bases with their proper colours and places the rgb values into the numpy array.
def base2color(lock, mp_arr, mp_arr2, tdone, gene, splits, dim, id):
    # Turn the colours array into a numpy array, but still use the same shared memory between threads.
    arr = np.frombuffer(mp_arr.get_obj(), dtype=np.uint8)
    colours = arr.reshape((dim,dim,3))

    # Turn the tmp_colours array into a numpy array, but still use the same shared memory between threads.
    arr2 = np.frombuffer(mp_arr2.get_obj(), dtype=np.uint8)
    tmp_colours = arr2.reshape((args.threads,lastThreadSplits,dim,3))

    # Select the proper tmp_colours array for this thread.    
    curColours = tmp_colours[id]
    
    # Index of our position in the genetic sequence.
    index = 0
    # Vertical
    for i in range(splits):
        max = dim
        if(len(gene) - (i+1)*dim <= 0 and id+1 == args.threads):
            max = len(gene) - (i)*dim
        # Horizontal
        for j in range(max):
            if(ord(gene[index]) < 86):
                 # Cytosine
                if(ord(gene[index]) == 67):
                    curColours[i, j] = cytosineColour
                # Guanine
                elif(ord(gene[index]) == 71):
                    curColours[i, j] = guanineColour
                # Adenine
                elif(ord(gene[index]) == 65):
                    curColours[i, j] = adenineColour
                # Thymine or Uracil
                elif(ord(gene[index]) == 84 or ord(gene[index]) == 85):
                    curColours[i, j] = thymineColour
            #Lowercase
            elif(ord(gene[index]) > 85):
                # Cytosine
                if(ord(gene[index]) == 99):
                    curColours[i, j] = cytosineColour
                # Guanine
                elif(ord(gene[index]) == 103):
                    curColours[i, j] = guanineColour
                # Adenine
                elif(ord(gene[index]) == 97):
                    curColours[i, j] = adenineColour
                # Thymine or Uracil
                elif(ord(gene[index]) == 116 or ord(gene[index]) == 117):
                    curColours[i, j] = thymineColour

            # Keep track of where we are in the sequence and break once sequence is finished.
            index += 1
    # Clear RAM
    geneThread[id] = ""
    t_done.value+=1
    isLast(dim, colours, tmp_colours)

def isLast(dim, colours, tmp_colours):
    global t_done
    global splitsPerThread
    if(t_done.value >= args.threads):
        for i in range(args.threads):
            if(i == args.threads-1):
                colours[(args.threads - 1)*splitsPerThread:(args.threads-1)*splitsPerThread + lastThreadSplits, 0:dim] = tmp_colours[i]
            else:
                tmp = np.resize(tmp_colours[i], (splitsPerThread, dim, 3))
                colours[i*splitsPerThread:(i+1)*splitsPerThread, 0:dim] = tmp

        print("Done parsing sequence, saving image now...")

        # Save the image and apply scaling.
        saveImg(colours, dim)
    else:
        sys.exit()

def geneSplit(gene, dim):
    splitGene = []
    for i in range(dim):
        splitGene.append(gene[(dim*i):(dim*(i+1))])

    splitsPerThread = int(dim/(args.threads))
    lastThreadSplits = splitsPerThread

    if(dim%args.threads != 0):
        splitsPerThread = int(dim/(args.threads))
        lastThreadSplits = dim - (args.threads-1)*splitsPerThread

    geneThread = [gene[i:i+(splitsPerThread*dim)] for i in range(0, len(gene), (splitsPerThread*dim))]
    if(args.threads>1):
        geneThread[args.threads-1] += geneThread[args.threads]
        geneThread[args.threads] = ""
   
    return geneThread, lastThreadSplits, splitsPerThread

# Hold the genetic sequence as a string.
sequence = ""

# Default base colours.
adenineColour = [239, 71, 111]
thymineColour = [255, 209, 102]
cytosineColour = [6, 201, 150]
guanineColour = [17, 138, 178]
uracilColour = [204, 102, 255]

# See if custom colours are specified.
if(args.A is not None):
    adenineColour = hex2RGB(args.A)
if(args.T is not None):
    thymineColour = hex2RGB(args.T)
if(args.C is not None):
    cytosineColour = hex2RGB(args.C)
if(args.G is not None):
    guanineColour = hex2RGB(args.G)

# See how user wants to input the sequence.
if(args.file is not None):
    sequence = getGeneFromFile(args.file)
else:
    sequence = getGeneFromCLI()

sequenceLength = len(sequence)
print("Length: "+ str(sequenceLength) + " bases")

# Find optimal size square for the array.
dim = findSquareSize(sequence)

# See if there is any point to using multiple threads.
if(dim < args.threads):
    args.threads = dim

print("Array Dimmention: " + str(dim) + " bases")

geneThread, lastThreadSplits, splitsPerThread = geneSplit(sequence, dim)

# Values which need to be writable and viewable by all threads.
t_done = mp.Value('i', 0)
# Allocate memory for the arrays. We will then turn these into numpy arrays inside each thread.
mp_arr = mp.Array(c.c_uint8, dim*dim*3)
mp_arr2 = mp.Array(c.c_uint8, lastThreadSplits*dim*3*args.threads)

lock = mp.Lock()
j = splitsPerThread
# Start the threads
for i in range(args.threads):
    # If last thread, pass in the proper number of splits for the last thread.
    if(i == args.threads-1):
        j = lastThreadSplits
    mp.Process(target=base2color, args=(lock, mp_arr, mp_arr2, t_done, geneThread[i], j, dim, i,)).start()