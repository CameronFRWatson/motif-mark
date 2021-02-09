## Motif Marking Tool
Cameron Watson

### Description

*motifMark.py* is a tool for visualizing the location of short motifs in a fragment of gene or transcript. The program takes a FASTA file as input, with exons capitalized and introns in lower-case (see **Example Input**). Additionally, a file containing unique motifs to search for is also required, with one motif per line. Motifs can include IUPAC degenerate bases. This program will return all possible hits to each motif. The output is a graphic visualization of each input sequence, with the introns, exons, and various motifs all denoted (see **Example Output**). 

### Usage

*motifMark.py* requires the following Python modules:

 - argparse
 - re
 - pycario
 - numpy

The program can be called from the command-line with two flags for the input files:

```
./motifMark.py -f input.fasta -m motifs.txt
```

*current work in progress*

### Example Input

### Example Output