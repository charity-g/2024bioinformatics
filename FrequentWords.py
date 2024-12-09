# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 15:46:51 2023

@author: cherr
"""

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

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


#%%
FrequencyMap("Royalty", 3)

