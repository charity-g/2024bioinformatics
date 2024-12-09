# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:23:17 2023

@author: cherr
"""
import time
#%%
#UTIL FUNCTIONS

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 


#Input: A string Text and an integer k.
#Output: All possible k-mers in Text. 
def getKmers(Text, k):
    kmers = []
    for i in range(len(Text) - k + 1):
        kmers.append(Text[i:i+k])
    return kmers

#Input: A string Text and an integer k.
#Output: All most frequent k-mers in Text. 
def FrequencyMap(Text, k):
    #generate all k-length substr
    map_count = {}
    for i in range(len(Text) - k + 1):
        substr = Text[i:i+k]
        count = PatternCount(substr, Text)
        map_count[substr] = count
        
    return map_count

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for [key, val] in freq.items():
        if (val == m):
            words.append(key)
    return words


#Find the reverse complement of a DNA string.
#     Input: A DNA string Pattern.
#     Output: The reverse complement of Pattern
def ReverseComplement(Pattern):
    return Complement(Pattern[::-1])

def Complement(Pattern):
    complement_map = {'A': 'T',
                      'T': 'A',
                      'C': 'G',
                      'G': 'C'} 
    return ''.join(map(lambda char: complement_map[char] , Pattern))
    

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions


# Purpose: Find patterns forming clumps in a string.

# Output: All distinct k-mers forming (L, t)-clumps in Genome.
# input: string Genome, and integers k, L, and t
#        where k = length of k-mer
#              L = length of sliding window to chcek in
#              t = # of occurences k-mer must show up in
def ClumpFinder(Genome, k, L, t):
    clumping_patterns = set()
    for i in range(len(Genome) - L + 1):
        windowText = Genome[i:i + L]
        window_kmers = getKmers(windowText, k)
        for kmer in window_kmers:
            if (PatternCount(kmer, windowText) >= t):
                clumping_patterns.add(kmer)
    return clumping_patterns

