# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 16:20:41 2024

@author: cherr
"""

#%%
"""
Unit3 = functions on Motifs
"""
#%%

import sys 
sys.path.append('C:\\Users\\cherr\\OneDrive\Documents\\bioinformatics 23-24')


# Input:  A set of kmers Motifs where each Motif is length k
# Output: Count(Motifs)
#INVARIANT: Motifs is not empty
def Count(Motifs):
    k = len(Motifs[0])
    count = {} # initializing the count dictionary
    for nucleotide in "ATCG":
        row = [0] * k
        count[nucleotide] = row
    
    for i in range(k):
        for motif in Motifs:
            symbol = motif[i]
            count[symbol][i] = count[symbol][i] + 1
            
    return count

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
# INVARIANT: len(Motifs) > 0
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = Count(Motifs)
    for nucleotide in count.keys():
        row = [0] * k
        for i in range(len(count[nucleotide])):
            row[i] = count[nucleotide][i] / t
        profile[nucleotide] = row
    return profile


# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
# INVARIANT: len(Motifs) > 0
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for i in range(k):
        max_symbol = ''
        max_count = 0
        for symbol in "ATCG":
            if (count[symbol][i] > max_count):
                max_count = count[symbol][i]
                max_symbol = symbol
        consensus = consensus + max_symbol
    return consensus


# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = Consensus(Motifs)
    score = 0
    for i in range(k):
        for symbol in "ATCG":
            if (symbol != consensus[i]):
                score = score + count[symbol][i]
    return score


# Motif Finding Problem:
#  Given a collection of strings, find a set of k-mers, one from each string, that minimizes the score of the resulting motif.
#  Input: A collection of strings Dna and an integer k.
#  Output: A collection Motifs of k-mers, one from each string in Dna, minimizing Score(Motifs) among
#     all possible choices of k-mers.


#BruteForceMotifSearch
#considers every possible choice of k-mers Motifs from Dna (one k-mer from each string of n nucleotides)
#returns the collection Motifs having minimum score.
#ISSUES: runtime = (n-k+1)^t where n = length of DNA string and t is number of strings 
#(n-k+1)t)⋅k⋅t

# Input:  String Text and profile matrix Profile;
# INVARIANT: text must be same kmer length as profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    # insert your code here
    pr = 1
    for i in range(len(Text)):
        symbol = Text[i]
        pr = pr * Profile[symbol][i]
    return pr

#Input: text, k, profile
#Output: kmer in the text with the highest probability
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
def ProfileMostProbableKmer(text, k, profile):
    max_pr = -1.0
    max_pr_kmer = ''
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        pr = Pr(kmer, profile)
        if (pr > max_pr):
            max_pr = pr
            max_pr_kmer = kmer
    return max_pr_kmer
    

# Input:  Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: A list of best motifs
def GreedyMotifSearch(Dna, k, t):
    # type your GreedyMotifSearch code here.
    n = len(Dna[0])
    BestMotifs = []
    ## initalize "best" motifs
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs = Motifs
    
    return BestMotifs
    
#%%
"""Improved functions"""

# input: a list of strings Motifs as input 
# output: the count matrix of Motifs with pseudocounts as a dictionary of lists.
def CountWithPseudocounts(Motifs):
    k = len(Motifs[0])
    countMatrix = {}
    for nucleotide in "ATCG":
        row = [1] * k
        countMatrix[nucleotide] = row
    
    for i in range(k):
        for motif in Motifs:
            symbol = motif[i]
            countMatrix[symbol][i] = countMatrix[symbol][i] + 1
            
    return countMatrix

# input:  that takes a list of strings Motifs as input 
# output: the profile matrix of Motifs with pseudocounts as a dictionary of lists. Then add this function to Motifs.py. Make sure to use CountWithPseudocounts(Motifs) as a subroutine!
    
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)
    denominator = t + 4
    for nucleotide in count.keys():
        row = [0] * k
        for i in range(k):
            row[i] = count[nucleotide][i] / denominator 
        profile[nucleotide] = row
    return profile


# Input:  Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: A list of best motifs
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    # type your GreedyMotifSearch code here.
    n = len(Dna[0])
    BestMotifs = []
    ## initalize "best" motifs
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        #score is still original
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs = Motifs
    
    return BestMotifs