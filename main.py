import sys 
import os
sys.path.append('C:\\Users\\cherr\\OneDrive\Documents\\bioinformatics 23-24')
import unit1_parsingDNA as u1

#%%
f = open(r"C:\Users\cherr\OneDrive\Documents\bioinformatics 23-24\ecoli.txt", "r")
ecoli_genome = f.readline()
f.close() 

a = u1.SymbolArray(ecoli_genome, 'C')
print(a)
