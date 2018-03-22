# Software Carpentry Biostrings Workshop

This set of examples demonstrates some of the functionality that the [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) package (part of the [Bioconductor](https://bioconductor.org) project) adds to R. We will work through some fundamental manipulations of nucleotide data.

## Install And Load the Biostrings Package

Before we begin, we might have to install the Biostrings library. We can check if it is installed by attempting to load it:

~~~R
library(Biostrings)
~~~

If the library is installed, you'll get lots of messages on the screen. However, if the library is **not** installed, you'll get an error message:

~~~
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

There are several ways of getting nucleotide data into R. The most simple is to type (or copy & paste) it in. For example, the first 50 nucleotides from the human TP53 gene can be loaded into R as:

~~~R
s <- DNAString("GAGACAGAGTCTCACTCTGTTGCACAGGCTGGAGTGCAGTGGCACAATCT")
~~~

Once we've loaded this, we can examine it. Running `print(s)`, shows us a representation of the object itself:

~~~R
> print(s)
  50-letter "DNAString" instance
seq: GAGACAGAGTCTCACTCTGTTGCACAGGCTGGAGTGCAGTGGCACAATCT
~~~

Although OK for short sequences, we often will need to load sequences from file. Reading a FASTA file is easy. For example, the file `TP53-exons.fasta` contains the exons from the TP53 transcript [TP53-222](http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000141510;r=17:7661779-7687550;t=ENST00000617185) (ENST00000617185.4). As this file contains 12 exons with a total length of 2,724 nucleotides, we don't want to type it in! We can load the FASTA file with:

~~~R
TP53.222 <- readDNAStringSet("TP53-exons.fasta")
~~~

When we examine the object we loaded, we see something a bit different to before:

~~~R
> print(TP53.222)
  A DNAStringSet instance of length 12
     width seq                                                                              names               
 [1]   174 GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACT...CGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGG ENST00000617185 E...
 [2]   102 CAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCC...CCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACT ENST00000617185 E...
 [3]    22 ACTTCCTGAAAACAACGTTCTG                                                           ENST00000617185 E...
 [4]   279 TCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTG...TCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACG ENST00000617185 E...
 [5]   184 TACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCC...GCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATG ENST00000617185 E...
 ...   ... ...
 [8]   137 TGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGT...GCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAG ENST00000617185 E...
 [9]    74 CACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAG       ENST00000617185 E...
[10]   133 GACCAGACCAGCTTTCAAAAAGAAAATTGTTAAAGAGAG...GTTACTTCCTGATAAACTCGTCGTAAGTTGAAAATATT ENST00000617185 E...
[11]   107 ATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTG...GCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAG ENST00000617185 E...
[12]  1289 CCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCA...CAATAAAACTTTGCTGCCACCTGTGTGTCTGAGGGGTG ENST00000617185 E...
~~~

There are several things to note here:

1. We've loaded 12 DNA sequences, not one;
2. We have a `DNAStringSet` object (not a `DNAString`);
3. Each sequence has a different length; and
4. Each sequence has a name

We can get the *width* of each sequence:

~~~R
> width(TP53.222)
 [1]  174  102   22  279  184  113  110  137   74  133  107 1289
~~~

Similarly, we can get the sequence names:

~~~R
> names(TP53.222)
 [1] "ENST00000617185 ENSE00003753508 exon:protein_coding" "ENST00000617185 ENSE00002667911 exon:protein_coding"
 [3] "ENST00000617185 ENSE00002419584 exon:protein_coding" "ENST00000617185 ENSE00003625790 exon:protein_coding"
 [5] "ENST00000617185 ENSE00003518480 exon:protein_coding" "ENST00000617185 ENSE00003723991 exon:protein_coding"
 [7] "ENST00000617185 ENSE00003712342 exon:protein_coding" "ENST00000617185 ENSE00003725258 exon:protein_coding"
 [9] "ENST00000617185 ENSE00003786593 exon:protein_coding" "ENST00000617185 ENSE00003735852 exon:protein_coding"
[11] "ENST00000617185 ENSE00003634848 exon:protein_coding" "ENST00000617185 ENSE00003492844 exon:protein_coding"
~~~

We can extract a single DNA sequence from the set using square brackets:

~~~R
> TP53.222[1]
  A DNAStringSet instance of length 1
    width seq                                                                               names               
[1]   174 GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACT...GCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGG ENST00000617185 E...
~~~

What has this given us back? A `DNAStringSet` containing a single sequence. This might not be what we wanted; we might want a single `DNAString`. We can get this using the double-square-bracket syntax:

~~~R
> TP53.222[[1]]
  174-letter "DNAString" instance
seq: GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTT...GTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGG
~~~

We 

## Extract Subsequences


