library(Biostrings)

# 1: Create a DNA String:
s <- DNAString('ATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTG')

# Basic Statistics on the sequence:
length(s)
subseq(s, start=1, end=10)

# Frequencies:
alphabetFrequency(s)
letterFrequency(s, letters=c('C', 'G'))
dinucleotideFrequency(s)

# Conversions:
complement(s)
reverse(s)
reverseComplement(s)
codons(s)
translate(s, genetic.code=GENETIC_CODE)
GENETIC_CODE_TABLE
getGeneticCode('SGC1')

# Now, let's play with the whole Human genome:
library(BSgenome.Hsapiens.UCSC.hg38)
Hsapiens
Hsapiens$chr1
seqnames(Hsapiens)

# Extract a more substantial part:
TP53 <- reverseComplement(subseq(Hsapiens$chr17, start=7668402, end=7687550))
letterFrequencyInSlidingView(TP53, view.width=100, letters=DNA_BASES)

# Demonstrate PCR:
forward.primer <- DNAString('TCTCATGCTGGATCCCCACT')
reverse.primer <- DNAString('AGTCAGAGGACCAGGTCCTC')
res <- matchProbePair(forward.primer, reverse.primer, Hsapiens$chr17)

# Show PWM matching:
TP53.binding <- readDNAStringSet('TP53-binding.fasta', use.names=FALSE)
consensusString(TP53.binding)
pwm <- PWM(TP53.binding)
p.matches <- matchPWM(pwm, Hsapiens, min.score=0.95)
