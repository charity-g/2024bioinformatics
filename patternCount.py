# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:23:17 2023

@author: cherr
"""

#%% Sliding Window PatternCount

def PatternCount(Pattern, Text):
    pattern_len = len(Pattern)
    match_count = 0
    #sliding window
    for i in range(0, len(Text)):
        substr_to_match = Text[i:len(Text)]
        if (pattern_len <= len(substr_to_match)):
            #match pattern char to substr char
            match = True
            for j in range(0, pattern_len):
                if (substr_to_match[j] != Pattern[j]):
                    match = False
                    break
            if (match):
                match_count = 1 + match_count  
            
    return match_count
                
            
def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 
        
    
#%%
test1 = PatternCount("ACTAT", "ACAACTATGCATACTATCGGGAACTATCCT")
test2 = PatternCount("ATA", "CGATATATCCATAG")
test3 = PatternCount("GCG", "GCGCG")
#print(test3)

Text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
Pattern = "TGATCA"

# Finally, let's print the result of calling PatternCount on Text and Pattern.
print(PatternCount(Pattern, Text))

#%%
def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    # your code here
    return positions