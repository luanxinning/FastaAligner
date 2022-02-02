# FastaAligner  
## Description:
    This program will calculate the alignment scores of input aligned DNA sequnces, as well as their identities and gaps. If there are more than two sequencese, the result will be listed in pairs.  
    For parameters of calculating the score, this program allow the user to provide a file with different scores. If there is no parameter file, the script will use default parameters.  
    This program is able to check whether the input file is a FASTA file. If not, the program will detect it as an error, telling the user this isn't a fasta file and won't run the following part.  
## Procedure:
    1. Check whether the input file is a fasta file.  
    2. Check whether the user uses a parameter file. If so, use those parameters. If not, usedefault parameters.  
    3. Store the headers to be keys in a dictionary, and its sequence to be the key's value.  
    4. Pairwise the nucleotide sequences.  
    5. Calculate the score.  
    6. write the output file.   
