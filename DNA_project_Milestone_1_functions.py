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
"""
@author: Adesh Patel
"""
def rna2codon(rna):
  protein_string = ""
  string_length = len(rna)
  for i in range(int((string_length - string_length % 3) / 3)):
    index = i * 3
    if (rna[index] == 'U'):
      if (rna[index + 1] == 'U'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C'):
          protein_string += 'F'
        elif (rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'L'

      elif (rna[index + 1] == 'C'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'S'

      elif (rna[index + 1] == 'A'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C'):
          protein_string += 'Y'

      elif (rna[index + 1] == 'G'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C'):
          protein_string += 'C'

        elif (rna[index + 2] == 'G'):
          protein_string += 'W'

    elif (rna[index] == 'C'):
      if (rna[index + 1] == 'U'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'L'

      elif (rna[index + 1] == 'C'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'P'

      elif (rna[index + 1] == 'A'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C'):
          protein_string += 'H'

        elif (rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'Q'

      elif (rna[index + 1] == 'G'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'R'

    elif (rna[index] == 'A'):
      if (rna[index + 1] == 'U'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A'):
          protein_string += 'I'

        elif (rna[index + 2] == 'G'):
          protein_string += 'M'

      elif (rna[index + 1] == 'C'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'T'

      elif (rna[index + 1] == 'A'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C'):
          protein_string += 'N'

        elif (rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'K'

      elif (rna[index + 1] == 'G'):
        if (rna[index + 2] == 'U'or rna[index + 2] == 'C'):
          protein_string += 'S'

        elif (rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'R'

    elif (rna[index] == 'G'):
      if (rna[index + 1] == 'U'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'V'

      elif (rna[index + 1] == 'C'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'A'

      elif (rna[index + 1] == 'A'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C'):
          protein_string += 'D'

        elif (rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'E'

      elif (rna[index + 1] == 'G'):
        if (rna[index + 2] == 'U' or rna[index + 2] == 'C' or 
            rna[index + 2] == 'A' or rna[index + 2] == 'G'):
          protein_string += 'G'
    
  return protein_string;

def locate_substring(dna_snippet, dna):
    x = []
    for i in range(len(dna)):
        if dna.startswith(dna_snippet, i):
            x.append(i)
    return x

def hamming_dist(dna1, dna2):
    x = 0
    for i in range(len(dna1)):
        if dna1[i] != dna2[i]:
            x += 1
    return x
