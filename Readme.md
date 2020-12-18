# Hirschberg Algorithm for Global Alignment
**Developed by:**

Jiayang Zheng: [zheng68@illinois.edu](zheng68@illinois.edu)

Kehang Chang: [kehangc2@illinois.edu](kehangc2@illinois.edu)

In this repository is an implementation of Hirschberg Algorithm proposed by [Dan Hirschberg](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm).
To provide a comparison of memory efficiency, we have also developed the Needleman-Wunsch algorithm.

In order to run the Hirschberg algorithm on custom input sequences, you would do
```bash
python run_hirschberg.py --path output_dir --file_name hirschberg_output.txt --delta score_matrix.txt 
--keys A,C,T,G,- --alignment1 sequence1.txt --alignment2 sequence2.txt --match 1 --mismatch -1 --gap -1
```
Where

**path:** This will be the output directory where the program store alignment score and corresponding alignments.
This arguments will be expected as a directory with "/" at the end. If no directory is specified, the program will
print the output in stdout.

**file_name:** This will be the filename where the alignment score and alignments will be stored. This argument has a 
default value of "hirschberg_output.txt" so the corresponding output will be in a text file name "hirschberg_output".

**delta:** This will be the file where a scoring function is stored. We provide BLOSUM62 scoring matrix in this repo.
If no delta file is specified, the program will make one using keys, match score, mismatch score, and gap penalty.

**keys:** This will be the list of symbols used by the sequences. When delta is given, this argument will be override by
the keys produced from delta.

**alignment1:** This will be the first sequence to be aligned. This argument can be either a string of sequence or a file that
contains the sequence.

**alignment2:** This will be the second sequence to be aligned. This argument can be either a string of sequence or a file that
contains the sequence.

**match:** This will be the match score for scoring function. This argument will only be effective when delta is not
specified. This argument has a default value of 1.

**mismatch:** This will be the mismatch score for scoring function. This argument will only be effective when delta is not
specified. This argument has a default value of -1.

**gap:** This will be the gap penalty for scoring function. This argument will only be effective when delta is not
specified. This argument has a default value of -1.

The command to run Needleman-Wunsch algorithm is similar:
```bash
python run_needleman_wunsh.py --path output_dir --file_name needleman_wunsch_output.txt --delta score_matrix.txt 
--keys A,C,T,G,- --alignment1 sequence1.txt --alignment2 sequence2.txt --match 1 --mismatch -1 --gap -1
```

This repository is meant to be a convenient tool for the usage of global alignment. As a sidenote, both `run_hirschberg.py`
and `run_needleman_wunsch.py` will give maximum memory usage and running time as you can see the resources used.
We provide two sequences `v.txt` and `w.txt` that use the scoring matrix stored in `BLOSUM62.txt` as a demonstration of
how Hirschberg can be much more space-efficient than Needleman-Wunsch.

#Scoring Matrix Format
The scoring matrix file should basically look like a matrix and have no other context besides the matrix.
We give an example scoring matrix BLOSUM62.txt as follows:
```
                           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                        A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
                        R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
                        N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
                        D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
                        C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
                        Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
                        E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
                        G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
                        H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
                        I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
                        L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
                        K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
                        M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
                        F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
                        P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
                        S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
                        T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
                        W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
                        Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
                        V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
                        B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
                        Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
                        X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
                        * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
``` 