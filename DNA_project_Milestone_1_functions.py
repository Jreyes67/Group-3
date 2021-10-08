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
    @author: Jose Reyes
    """
def mendels_law(hom, het, rec): #Given 
    
    ntotal = (hom + het + rec) 
    
    dominantHomo = (rec/ntotal) * ((rec - 1)/(ntotal - 1)) #rec = recessive
    
    heterozygous = (het/ntotal) * ((het - 1)/(ntotal - 1)) #het = heterozygous
    
    recessiveHomo = (hom/ntotal) * (hom/(ntotal - 1)) + (hom/ntotal) * (hom/(ntotal-1)) #hom = homozygous
    
    probabilityRecessive = dominantHomo + heterozygous * (1/4) + recessiveHomo * (1/2)
    
    print (1-probabilityRecessive)
   
def fibonacci_rabbits(n, k): #n = months and k = pairs of offspring
    baby, adult = 1,1
    for i in range(n - 1):
        baby,adult = adult,adult + (baby * k)
    print(baby)
    
def gc_content(dna_list):
    max_value = 0
    for i in range (len(dna_list)):
        nucleotides = dna_list[i] == ('G') + dna_list[j] == ('C')
        gc_con = ((nucleotides/len(dna_list[i]))*100) #Percentage of Nucleotides
        if max_value < nucleotides:
            index = i
            max_gc = gc_con
    print(index, max_gc)

    
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
"""
@author: Zyaire Howard
"""
def count_dom_phenotype(genotypes):
    prob = 0
    for index in range(len(genotypes)):
        if index == 0:
            prob = genotypes[0] * 1
        if index == 1:
            prob += genotypes[1] * 2
        if index == 2:
            prob += genotypes[2] * 3
        if index == 3:
            prob += genotypes[3] * 4
        if index == 4:
            prob += genotypes[4] * 5
        if index == 5:
            prob += 0
    return prob
print(count_dom_phenotype([9,1,0,2,42,7]))

rna2codon={
'UUU': 'BRAKE', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'BRAKE', 'CUG': 'L',
'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'BRAKE', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'BRAKE',
'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'BRAKE', 'CCC': 'BRAKE', 'CCA': 'P', 'CCG': 'P',
'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'BRAKE', 'GCG': 'BRAKE',
'UAU': 'BRAKE', 'UAC': 'Y', 'UAA': 'BRAKE', 'UAG': 'H', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'BRAKE', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'BRAKE',
'UGU': 'C', 'UGC': 'C', 'UGA': 'BRAKE', 'UGG': 'W', 'CGU': 'BRAKE', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
'AGU': 'BRAKE', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'BRAKE', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def source_rna(protein):
    count=0
    for key in rna2codon.keys():

        if rna2codon[key] in protein:
            count+=1

        elif rna2codon[key]=='BRAKE':
            count+=len(protein)

    return count%1000000

print(source_rna("GU"))

rna2codon = {
'UUU': 'BRAKE', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'BRAKE', 'CUG': 'L',
'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'BRAKE', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'BRAKE',
'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'BRAKE', 'CCC': 'BRAKE', 'CCA': 'P', 'CCG': 'P',
'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'BRAKE', 'GCG': 'BRAKE',
'UAU': 'BRAKE', 'UAC': 'Y', 'UAA': 'BRAKE', 'UAG': 'H', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'BRAKE', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'BRAKE',
'UGU': 'C', 'UGC': 'C', 'UGA': 'BRAKE', 'UGG': 'W', 'CGU': 'BRAKE', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
'AGU': 'BRAKE', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'BRAKE', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def replaceCodons(rnasequence):
    i = 0
    protein = ""
    rnasequence = rnasequence.replace("\n","")
    
    while (i < len(rnasequence)):
        
        if(rna2codon[rnasequence[i:i+3]]== 'BRAKE'):
            break
            
        protein += rna2codon[rnasequence[i:i + 3]]
        i += 3
    protein = protein.replace('BRAKE','')
    
    return protein

def rna_splyce(dna, intron_list):
    intron_list.sort(reverse = True)
    
    for i in intron_list:
        dna = dna.replace(i,"")
        
    rna = dna.replace("T","U")    
    output=replaceCodons(rna)
    
    return output

rna_splyce("GTATCCTAGTAGGTTATCGTATCG", ['TCAGTACATATG','GGATAGTGACACTAG'])
