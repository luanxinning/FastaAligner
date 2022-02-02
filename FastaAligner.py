#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: DNA sequences Alignment
Date:2021-10-15
Author:Luan Xinning

Description:
    This program will calculate the alignment scores of input aligned DNA sequnces, as well as 
their identities and gaps. If there are more than two sequencese, the result will be listed in 
pairs.
    For parameters of calculating the score, this program allow the user to provide a file with 
different scores. If there is no parameter file, the script will use default parameters.
    This program is able to check whether the input file is a FASTA file. If not, the program will
detect it as an error, telling the user this isn't a fasta file and won't run the following part.

List of user defined functions:
    transition
        This function is used to test whether a substitution is a transition.    

Procedure:
    1. Check whether the input file is a fasta file.
    2. Check whether the user uses a parameter file. If so, use those parameters. If not, use
default parameters.
    3. Store the headers to be keys in a dictionary, and its sequence to be the key's value.
    4. Pairwise the nucleotide sequences
    5. Calculate the score
    6. write the output file   
    
Usage:
    python FastaAligner.py input_fasta.fna parameters.txt[optional] output_fasta.txt

"""

#######################################################################
# read the input file, which should be a fasta file
import sys  # the modle that we need to read files
f = open(sys.argv[1]) # open the fasta file
file = f.readlines()  # read it in lines. the lines will be store line by line in 'file'


#%%
# 1. check whether the input file is a fasta file
possible_chac = {'-', 'A', 'C', 'G', 'T'} # make a set that contains all of the possible 
# characters and symbols in DNA sequences
test = '' # define an empty string

for line in file:
    if not line.startswith('>'):  # for the lines in input file, if it doesn't start with
    # '>' which contained in headers, it is a DNA sequence.
        test += line.upper().rstrip() # with upper(), we can change all of the characters to
        # be uppercase. rstrip() is used to delete invisible ending symbols. Use '+=' to 
        # combine all of these lines.
test = set(test) # Change it into a set.

if len(test-possible_chac) != 0: # if the difference is not 0, that means besides the element
# in possible_chac, there are other elements in input file. If so, this file isn't a fasta file.
    raise Exception('The input file is not a fasta file. Please upload a fasta file')
    # Give an error and terminate running the program

# If the input file is a fasta file, the progrm continues


#%%
# 2. Check whether the user uses a parameter file. If so, use those parameters. If not, use
# default parameters.

if len(sys.argv) == 4: # in this case the user uses a parameters file 
    o = open(sys.argv[3],'w') # create the output file
    print('''Please upload the parameter file in this format:
    gap_value
    match
    transition_score
    tranversion
    ''') # Remind the user what the parameter file should look like
    p = open(sys.argv[2]) # open the parameter file
    parameters_file = p.readlines() # read it in lines. the lines will be store line by line
    # in 'parameters_file'
    parameters = [] # set an empty list
    for line in parameters_file:
        parameters.append(line)  # Store all of the lines from 'parameters_file' to parameters
    gap_value = int(parameters[0]) # get gap penalty
    match = int(parameters[1])     # get match score
    transition_score = int(parameters[2])  # get transition penalty
    tranversion = int(parameters[3])       # get transversion penalty
    p.close()   # close the parameter file
    
elif len(sys.argv) == 3:# in this case the user doesn't use a parameters file
    o = open(sys.argv[2],'w')  # create the output file
    # store all of the default parameters:
    gap_value = -1            # gap penalty
    match = 1                 # match score
    transition_score = -1     # transition penalty
    tranversion = -2          # transversion penalty


#%%
# 3. store the headers to be keys in a dictionary, and its sequence to be the key's value
fasta = {}  #Define an empty dictionary
key = ''    #Define key as an empty string
value = ''  #Define value as an empty string
for line in file:
    if line.startswith('>'): #Identify lines start with greater than sign,
        #this is the headers/IDs in the fasta file 
        key = line.strip()   # delete invisible ending symbols              
        value = ''           #Clean 'value' to be an empty string
        value2 = []          #clean 'value2' to be an empty list
    else:  # these lines should be sequences
        value1 = line.strip()  # delete invisible ending symbols
        value2.append(value1)  # store each line of sequences
        value = ''.join(value2) #combine all nucletide sequences from a single ID 
        #to be a string
        fasta[key] = value   #Identify this ID as a key in dictionary and set this string 
        # to be the header's value


#%%
# 4. Pairwise the nucleotide sequences
pairwise = []  # create an empty list to store all of the pairs
for i in range(len(fasta)):   # i will be used for the first sequence
        for j in range(i+1,len(fasta)):  # j will be used for the second sequence
            key1 = list(fasta.keys())[i] # get the ID of the first sequence
            key2 = list(fasta.keys())[j]  # get the ID of the second sequence
            value1 = list(fasta.values())[i] # get the nucletide sequence of the first sequence
            value2 = list(fasta.values())[j] ## get the nucletide sequence of the second sequence
            pairwise.append("{}\n{}\n{}\n{}".format(key1,value1,key2,value2))
            # store the 2 keys and 2 sequences into each element of 'pairwise' list


#%%
# 5. Calculate the score

# define an transition function
def transition(nucl1,nucl2):  # use nucl1 and nucl2 as input variables
    if nucl1 == 'A' or nucl1 == 'G': 
        if nucl2 == 'A' or nucl2 =='G': # Because we will exclude the possibilty of 'match',
        # we only need to figure out whether the two nucletide is A or G.
            return True  # If it is transition, return 'True' to the program.

Score = 0            # create 'Score' as 0 
scores = []          # create 'scores' as an empty list 
match_number = 0     # create 'match_number' as 0
matches = []         # create 'matchers' as an empty list
gap_number = 0       # create 'gap_number' as 0
gaps = []            # create 'gaps' as an empty list
for pairs in pairwise:  # for each pair of sequences:
    line1 = pairs.split('\n')[1].upper()  # extract the first sequence as line1
    line2 = pairs.split('\n')[3].upper()  # extract the second sequence as line2
    length = len(line1)         # get the length of an aligned sequence
    for i in range(length):     # enumerate all nucletides in the two sequences
        nucl1 = line1[i]        # get a nuclestide in the first sequence as nucl1
        nucl2 = line2[i]        # get a nuclestide in the second sequence as nucl2
        if nucl1 == '-' or nucl2 == '-':   # the condition of 'gap'
            Score += gap_value    # add the gap penalty to this pair's score
            gap_number += 1       # calculate how many gaps are there in this pair
        elif nucl1 == nucl2:      # the condition of 'match'
            Score += match        # add the match score to this pair's score
            match_number += 1     # calculate how many matched are there in this pair
        elif transition(nucl1, nucl2):    # the condition of transition
            Score += transition_score     # add the transition score to this pair's score
        else:                     # the condition of transversition
            Score += tranversion  # add the transvertion score to this pair's score
    scores.append(Score)          # after loop, store scores of each pair into 'scores'
    matches.append(match_number)  # store the number of matches of each pair into 'matches'
    gaps.append(gap_number)       # store the number of gaps of each pair into 'gaps'
    Score = 0                     # clear 'Score' into 0
    match_number = 0              # clear 'match_number' into 0
    gap_number = 0                # clear 'gap_number' into 0
    
#%%
# 6. write the output file
for i in range(len(pairwise)):    # enumerate within 'pairwise'
    identity_percentage = matches[i]/length     # calculate identity
    Identity = "Identity:{}/{}({}%)".format\
        (matches[i],length,round(identity_percentage * 100,2)) 
        # get the output format of identity
    gaps_percentage = gaps[i]/length            # calculate gaps
    Gaps = "Gaps:{}/{}({}%)".format\
        (gaps[i],length,round(gaps_percentage * 100,2))
        # get the output format of gaps
    a, b = pairwise[i].split('\n')[0::2]        # a,b will be the first and second sequence's header
    a = a[1:]   # deleate the greater than symble in 'a' header          
    b = b[1:]   # deleate the greater than symble in 'b' header
    o.write('{}-{}:{},{},Score={}\n'.format(a,b,Identity,Gaps,scores[i]))
    # write the output file
o.close()   # close o file
f.close()   # close f file


























    
    
   

        
