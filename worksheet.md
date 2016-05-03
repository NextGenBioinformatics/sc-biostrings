# Biostrings Workshop

All of the questions below can be answered using the code we worked through last week. Have a go, and if you get stuck peek at the `biostrings-demonstration.R` file in this repository.

Before you begin, you'll need to load up R, point it at the directory where you have downloaded these files, and load the Biostrings library:

~~~.R
library(Biostrings)
~~~

If this doesn't work, you might need to install Biostrings, via the [Bioconductor installer](http://bioconductor.org/install/):

~~~.R
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
~~~

### 1: Create a DNAString object in R to hold the following DNA string:

~~~
TTGGGTAGGGGAGAAGAATTTTGGGGCGATGAAACTCTATGAAAAGTTTGGGTTTGTACCAGTTGGTAAG
~~~

###2: What is the letter frequency of the string?

You should be able to use a single command to get the frequency of A, T, G, and C in the string.

### 3: How many GA pairs are present in the string?
Use a very similar command to the answer to question 2.

### 4: What is the GC content of the string?
This is (again) very similar to the last two questions, but involves two steps.

### 5: Read in the first sequence in the "mystery.fasta" file.
Remember that we only want **one** sequence back, and that reading in a FASTA file will give us a set of sequences.

### 6: How long is the sequence?

### 7: Is the string you found above present in the sequence?  If so, where?

### 8: What is the overall GC content of the loaded sequence?
See above!

### 9: Plot the GC content of the sequence across a sliding window
This is a bit more tricky...

### 10: Take a guess as to which organism this sequence might be from.
