# Bioinformatics-Sequence-Alignment
Iımplementation of Pairwise sequence alignment of amino acid (protein) sequences via dynamic programming algorithms, local (Smith-Waterman) & global (Needleman-Wunch).


## Input


The input to this implementation should be a text file containing two sequences of amino acids, written on separate lines. Script accepts this file, along with the following additional arguments as a command line arguments:

-> alignment algorithm: local (Smith-Waterman) or global (Needleman-Wunsch)

-> scoring matrix (a square matrix in a specific format)

-> gap opening penalty (a negative integer)

-> gap extension penalty (a negative integer)



## Output

The script outputs the aligned sequences in the classical 3-line notation (first line for the first sequence, second line for the '|' characters for the positions where there is a match between 2 sequences and space characters when there is no match, and third line for the second sequence) including '–' character for gaps in the aligned sequences.


### How to run the script:

```console
python pairwise_alignment.py input_file.txt --algorithm local --scoring_matrix matrix.txt --gap_open -5 --gap_extend -2 
```


### Example:

```console
python pairwise_alignment.py input.txt global BLOSUM62.txt -5  -2 
```
