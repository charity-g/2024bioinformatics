# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 21:24:17 2024

@author: cherr
"""
import sys 
import os
import random
sys.path.append('C:\\Users\\cherr\\OneDrive\Documents\\bioinformatics 23-24')
from unit3 import ProfileMostProbableKmer, ProfileWithPseudocounts, Score, Pr, Profile


# input: profile matrix Profile corresponding to a list of strings Dna
# output: a list of the Profile-most probable k-mers in each string from Dna.
def Motifs(Profile, Dna):
    t = len(Dna)
    k = len(Profile['A'])
    mostProbableKmers = [''] * t
    for i in range(t):
        mostProbableKmers[i] = ProfileMostProbableKmer(Dna[i], k, Profile)
    return mostProbableKmers

# input: A list of strings Dna, and integers k and t
# output: returns a list of t strings, each string being k char in length.
# description: choose a random k-mer from each of t different strings Dna

def RandomMotifs(Dna, k, t):
    n = len(Dna[0])
    M = n - k
    randomMotifs = [''] * t
    for i in range(t):
        start = random.randint(0, M)
        randomMotifs[i] = Dna[i][start:start + k]
    return randomMotifs

# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    # insert your code here
    CurrMotifs = RandomMotifs(Dna, k, t)
    BestMotifs = CurrMotifs
    while True:
        Profile = ProfileWithPseudocounts(CurrMotifs)
        CurrMotifs = Motifs(Profile, Dna)
        if Score(CurrMotifs) < Score(BestMotifs):
            BestMotifs = CurrMotifs
        else:
            return BestMotifs 
        

#%%
"""Gibbs Sampler"""

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    # your code here
    kmer_sum = 0 
    normalizedPr = {}
    for kmer in Probabilities:
        kmer_sum += Probabilities[kmer]
    for kmer in Probabilities:
        normalizedPr[kmer] = (Probabilities[kmer] / kmer_sum)
    return normalizedPr

# Input:  A dictionary Probabilities whose keys are k-mers 
#         and whose values are the probabilities of these kmers and sum up to one
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    outputkmer = '' # output variable
    roll = random.uniform(0, 1)
    startrange = 0
    for kmer in Probabilities:
        endrange = startrange + Probabilities[kmer]
        if (startrange <= roll) and (roll <= endrange):
            outputkmer = kmer
            break
        else: 
            startrange = endrange
    return outputkmer

# Purpose: randomly chooses a k-mer from Text based on a profile matrix profile
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: kmer string
def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {}
    for i in range(n - k + 1):
        kmer = Text[i:i + k]
        probabilities[kmer] = Pr(kmer, profile)
        
    normalized_pr = Normalize(probabilities)
    return WeightedDie(normalized_pr)

# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: Collection of "best" Motifs
def GibbsSampler(Dna, k, t, N):
    CurrMotifs = RandomMotifs(Dna, k, t)
    BestMotifs = CurrMotifs
    for j in range(1, N):
        i = random.randint(0, t - 1)
        BestMotifs.pop(i)
        profile = ProfileWithPseudocounts(BestMotifs)
        motif_i = ProfileGeneratedString(Dna[i], profile, k)
        BestMotifs.insert(i, motif_i)
        if Score(CurrMotifs) < Score(BestMotifs):
            BestMotifs = CurrMotifs
    return BestMotifs

Dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
k = 8
t = len(Dna)
N = len(Dna[0])
res = GibbsSampler(Dna, k, t, N)
print(res)