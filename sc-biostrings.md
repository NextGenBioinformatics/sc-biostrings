# Software Carpentry Biostrings Workshop

This set of examples demonstrates some of the functionality that the [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) package (part of the [Bioconductor](https://bioconductor.org) project) adds to R. We will work through some fundamental manipulations of nucleotide data.

## Install And Load the Biostrings Package

Before we begin, we might have to install the Biostrings library. We can check if it is installed by attempting to load it:

~~~R
library(Biostrings)
~~~

If the library is installed, you'll get lots of messages on the screen. However, if the library is **not** installed, you'll get an error message:

~~~R
Error in library(Biostrings) : there is no package called ‘Biostrings’
~~~

In this case, we need to install it from Bioconductor, using the .[Bioconductor installer](http://bioconductor.org/install/):

~~~R
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
~~~

Remember that installing the package is very different to loading it. So, once we've installed it, we still need to load it:

~~~R
library(Biostrings)
~~~

## Read In Nucleotide Data

There are several ways of getting nucleotide data into R. The most simple is to type (or copy & paste) it in. For example, the first 50 nucleotides from the human TP35 gene are:

~~~R
s <- DNAString("GAGACAGAGTCTCACTCTGTTGCACAGGCTGGAGTGCAGTGGCACAATCT")
~~~

Once we've loaded this, we can examine it. Running `print(s)`, shows us a representation of the object itself:

~~~R
> print(s)
  50-letter "DNAString" instance
seq: GAGACAGAGTCTCACTCTGTTGCACAGGCTGGAGTGCAGTGGCACAATCT
~~~

Although OK for short sequences, we often will need to load sequences from file. Reading a FASTA file is easy. For example, the file `rRNA.fasta` contains the sequences of all rRNA segences in the human genome. As this file contains 569 genes with a total length of 65,057 nucleotides, we don't want to type it in! We can load the FASTA directly from the GitHub page for this tutorial file with:

~~~R
rRNA <- readDNAStringSet("https://raw.githubusercontent.com/alastair-droop/sc-biostrings/master/data/rRNA.fasta")
~~~

When we examine the object we loaded, we see something a bit different to before:

~~~R
> print(rRNA)
  A DNAStringSet instance of length 569
      width seq                                                               names               
  [1]   120 AAACACGAACAGCCATGCCTGAACGGGCCCG...ACCAACTGGGAGTATCAGGTGCTCGAGGCTT ENSG00000200516|E...
  [2]   119 GTCTACGGCCATACCACCCTGAACGCGCCCG...GACCGCCTGGGAATACCGGGTGCTGTAGGCT ENSG00000199334|E...
  [3]   118 GTCTATGGCCGTAAGAGCCCGCAGGCACCCT...GACCATGTAGGAATGCCGGGTGCTCTAGGCT ENSG00000199806|E...
  [4]   115 GTCTAGGGCCATACCACACTGAAAGCGGCTG...AGAGATTTGGATGGGGACACAGAGCCAAACC ENSG00000199508|E...
  [5]   119 GTCTACGGCCATACCACCCTGAACGCGCCCG...GACCGCCTGGGAATACCGGGTGCTGTAGGCT ENSG00000201321|E...
  ...   ... ...
[565]   131 GTCTATGGTCATACTACCCTGAAAGTGCTTG...ATGATTTTGGCTTGTTTGAAGCATATAGGCC ENSG00000252231|E...
[566]   106 GTCTATGGCCATGTCACCCTGAACATTCTCC...ACTTGGATGGGAGAAAAGTGGGGCAGGGGCT ENSG00000252653|E...
[567]   108 TTCTACTGGCATACCACTCTGAACGTGTCTG...AGAGGGATGTTTGATCTGGCAATTGCTGAAG ENSG00000252553|E...
[568]   109 AGCTCCTGCCATAGCAACCTGAGGGCACCTC...GACTGCCTGAGCAGACCGGGGGCTGGAGGCT ENSG00000252182|E...
[569]   107 GTCTATGGCCGTAACACCCTGAATGCACCCG...ACTTGGATGGGAGAAAATGAAGAGCCAGGCT ENSG00000284736|E...
~~~

There are several things to note here:

1. We've loaded 569 DNA sequences, not one;
2. We have a `DNAStringSet` object (not a `DNAString`);
3. Each sequence has a different length; and
4. Each sequence has a name

We can get the *width* of each sequence:

~~~R
> width(rRNA)
  [1] 120 119 118 115 119 117 125 119 155 119 110 114 119 119 111 152 119  94 117 117 119 118 152
 [24] 118 120 117 119 109 113 119 115 119 119 126 119 111 121 100 113 110 116 111 115 119 119 120
 [47] 119 114 121 121 115 108 115 118  91 114 119 121 119 117 121 117 109 119 107 115 119 134 107
 ...
~~~

Similarly, we can get the sequence names:

~~~R
> names(rRNA)
  [1] "ENSG00000200516|ENST00000363646" "ENSG00000199334|ENST00000362464"
  [3] "ENSG00000199806|ENST00000362936" "ENSG00000199508|ENST00000362638"
  [5] "ENSG00000201321|ENST00000364451" "ENSG00000199843|ENST00000362973"
...
~~~

We can extract a single DNA sequence from the set using square brackets:

~~~R
> rRNA[1]
  A DNAStringSet instance of length 1
    width seq                                                                 names               
[1]   120 AAACACGAACAGCCATGCCTGAACGGGCCCGG...GACCAACTGGGAGTATCAGGTGCTCGAGGCTT ENSG00000200516|E...
~~~

What has this given us back? A `DNAStringSet` containing a single sequence. This might not be what we wanted; we might want a single `DNAString`. We can get this using the double-square-bracket syntax:

~~~R
> rRNA[[1]]
  120-letter "DNAString" instance
seq: AAACACGAACAGCCATGCCTGAACGGGCCCGGTTGCATCTGATTG...TACTTGGATGGGAGACCAACTGGGAGTATCAGGTGCTCGAGGCTT
~~~

Often, we know that out FASTA files contain only a single sequence, so we can read the file in and extract the first (and only) sequence from it in one line. We will now use this to download the human genome sequence for chromosome 17:

~~~R
> chr17 <- readDNAStringSet("https://raw.githubusercontent.com/alastair-droop/sc-biostrings/master/data/chr17.fasta.gz")[[1]]
trying URL 'https://raw.githubusercontent.com/alastair-droop/sc-biostrings/master/data/chr17.fasta.gz'
Content type 'application/octet-stream' length 23708649 bytes (22.6 MB)
==================================================
downloaded 22.6 MB

> print(chr17)
  83257441-letter "DNAString" instance
seq: NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
~~~

## Views and Matches

Once we have a sequence, we often want to identify multiple parts of it which might be of interest. For example, now that we have the genome sequence of chr17, we can extract all of the gene sequences. However, we need to know the regions we want to get back. We can get the start and end positions of all genes on chr17 from the GitHub project:

~~~R
chr17.genepos <- read.csv("https://raw.githubusercontent.com/alastair-droop/sc-biostrings/master/data/chr17-genes.csv", row.names=1)
~~~

This gives us a `data.frame` with three columns: name, start and end. We can use the start and end positions to pull out the gene sequences for these 3,014 genes.

~~~R
chr17.genes <- Views(chr17, start=chr17.genepos$start, end=chr17.genepos$end, names= chr17.genepos$name)
~~~

We now have a set of sequences from chr17 that correspond to the genes. We can now extract a specific gene by its name. For example, to pull out the sequence of TP53:

~~~R
> chr17.genes[['TP53']]
  25772-letter "DNAString" instance
seq: GAGACAGAGTCTCACTCTGTTGCACAGGCTGGAGTGCAGTGGCAC...AGCGCCAGTCTTGAGCACATGGGAGGGGAAAACCCCAATCCCATC
~~~

We can search for small sequences across our whole chromosome. For example, we can find the positions of all the `TATA` sequences:

~~~R
> matchPattern('TATA', chr17)
  Views on a 83257441-letter DNAString subject
subject: NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
views:
            start      end width
     [1]    61167    61170     4 [TATA]
     [2]    61483    61486     4 [TATA]
     [3]    61732    61735     4 [TATA]
     [4]    65802    65805     4 [TATA]
     [5]    65804    65807     4 [TATA]
     ...      ...      ...   ... ...
[336405] 83246832 83246835     4 [TATA]
[336406] 83246929 83246932     4 [TATA]
[336407] 83246975 83246978     4 [TATA]
[336408] 83247041 83247044     4 [TATA]
[336409] 83247062 83247065     4 [TATA]
~~~

As we can see, there are 336,409 instances of the exact sequence `TATA` on chr17.

## Alphabet Frequency

The alphabet frequency across DNA sequences is a common problem. Biostrings makes this easy. We can see the distribution of nucleotides very easily:

~~~R
> alphabetFrequency(chr17)
       A        C        G        T        M        R        W        S        Y        K 
22639499 18723944 18851500 22705261        0        0        0        0        0        0 
       V        H        D        B        N        -        +        . 
       0        0        0        0   337237        0        0        0 
~~~

We can see that there are many A, C, G, and T nucleotides, but also many N's.

We can look at the overall GC content similarly:

~~~R
> letterFrequency(chr17, letters='GC', as.prob=TRUE)
      G|C 
0.4513163 
~~~

Compare this to the genes:

~~~R
> mean(letterFrequency(chr17.genes, letters='GC', as.prob=TRUE)[,1])
[1] 0.497416
~~~
