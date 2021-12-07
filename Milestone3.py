#assert find_splice('GTA', 'ACGACATCACGTGACG') == [2, 6, 8] #added from milestone 3 relate page

import csv
f = "ms3-dna-mammuthus.txt"

def load_file(file):
    myfile = open( file )
    mydata = myfile.readlines()
    myfile.close()
    return mydata

def assemble_genome(dna_list):
    length = len(dna_list)
    combine = [[0 for _ in range(length)] for _ in range(length)]
    for a in range(length):
        for b in range(length):
            if a == b:
                continue
            c, d = dna_list[a], dna_list[b]
            letters = len(c)
            for z in range(1, letters):
                if d.startswith(c[z:]):
                    combine[a][b] = letters - z
                    break


def assemble_genome2(dna_list=load_file(f)):
  mammoth = dna_list[0]
  while len(mammoth) < 50:
    for i in range(0,len(dna_list)):
      string = dna_list[i]
      if string[0:8] == mammoth[-8:]:
        mammoth += string[8:]
        break
    for i in range(0,len(dna_list)):
      string = dna_list[i]
      if string[-8:] == mammoth[0:8]:
        mammoth = string + mammoth[8:]
        break
    print(len(mammoth))
  return mammoth

#print(assemble_genome2(load_file(f)))
