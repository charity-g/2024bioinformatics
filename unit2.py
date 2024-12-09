# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 16:09:51 2024

@author: cherr
"""
import sys 
import time
import os
sys.path.append('C:\\Users\\cherr\\OneDrive\Documents\\bioinformatics 23-24')
import unit1_parsingDNA as u1


# Input:  Strings Genome and symbol
# Output: FasterSymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = u1.PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def FasterSymbolArray_v2(Genome, symbol):
    t = time.time()
    n = len(Genome)
    array = {}
    ExtendedGenome = Genome + Genome[0:n//2]
    curr_symbol_count = 0
    for i in range(n + n//2 - 1):
        if (time.time() - t > 60):
            print('overtime')
            break
        if (ExtendedGenome[i] == symbol):
            curr_symbol_count = curr_symbol_count + 1
        if (i >= n//2 - 1):
            array[i - n//2 + 1] = curr_symbol_count
            print(f'array[{i - n//2 + 1}] = {curr_symbol_count}')
            if (ExtendedGenome[i - n//2 + 1] == symbol):
                curr_symbol_count = curr_symbol_count - 1
            
    return array

# DO NOT USE as this is Omega(n^2)
# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def SymbolArray(Genome, symbol):
    array = {}
        # n = len(Genome)
    # ExtendedGenome = Genome + Genome[0:n//2]
    # for i in range(n):
    #     count = PatternCount(symbol, ExtendedGenome[i:i+n//2])
    #     array[i] = count
    return array


# Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArrayList(Genome):
    # your code here
    skew = [0]
    for i in range(0, len(Genome)):
        if (Genome[i] == 'G'):
            skew.append(skew[i] + 1)
        elif (Genome[i] == 'C'):
            skew.append(skew[i] - 1)
        else:
            skew.append(skew[i])
    return skew

# Input:  A String Genome
# Output: The skew array of Genome as a dict.
def SkewArrayDict(Genome):
    # your code here
    skew = {0:0}
    for i in range(0, len(Genome)):
        if (Genome[i] == 'G'):
            skew[i + 1] = skew[i] + 1
        elif (Genome[i] == 'C'):
            skew[i+1] = skew[i] - 1
        else:
            skew[i+1] = skew[i]
    return skew

# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    positions = [] # output variable
    # your code here
    skew_array = SkewArrayDict(Genome)
    
    min_skew = len(Genome)
    for i in skew_array:
        if (min_skew > skew_array[i]):
            min_skew = skew_array[i]
            positions = [i]
        elif (min_skew == skew_array[i]):
            positions.append(i)
    
    return positions

# Input:  Two strings p and q of EQUAL LENGTH
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    # your code here
    hamming_dist = 0
    for i in range(len(p)):
        if (p[i] != q[i]):
            hamming_dist = hamming_dist + 1
    return hamming_dist

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    # your code here
    for i in range(len(Text) - len(Pattern) + 1):
        subText = Text[i:i+len(Pattern)]
        if (HammingDistance(subText, Pattern) <= d):
            positions.append(i)
    return positions


# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count = count+1
    return count
a=list(range(5))
b=a
a[2]=12