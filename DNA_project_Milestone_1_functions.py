# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 16:29:23 2021

@author: logan
"""

def s(dna):
    dna_0 = 'ACGT'
    dictionary = {'A':0,'C':0,'G':0,'T':0}
    for letters in dna:
        if letters in dna_0:
            dictionary[ letters ] += 1
    dna_count = dictionary
    return dna_count
print (s("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"))
if s("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC") == {'A': 20, 'C': 12, 'G': 17, 'T': 21}:
    print('True')
else:
    print('False')

def dna2rna(dna):
    rna = ''
    for character in dna:
        if character == 'T':
            rna = rna + 'U'
        elif character == 'A':
            rna = rna + 'A'
        elif character == 'G':
            rna = rna + 'G'
        elif character == 'C':
            rna = rna + 'C'
    return rna
print (dna2rna("GATGGAACTTGACTACGTAAATT"))
if dna2rna("GATGGAACTTGACTACGTAAATT") == "GAUGGAACUUGACUACGUAAAUU":
    print('True')
else:
    print('False')
    
def reverse_compliment(dna):
    compliment = ''
    for character in dna:
        if character == 'T':
            compliment = compliment + 'A'
        elif character == 'A':
            compliment = compliment + 'T'
        elif character == 'G':
            compliment = compliment + 'C'
        elif character == 'C':
            compliment = compliment + 'G'
    return compliment[::-1]
print (reverse_compliment("AAAACCCGGT"))
if reverse_compliment("AAAACCCGGT") == "ACCGGGTTTT":
    print('True')
else:
    print('False')