#f = open('ms3-dna-mammuthus.txt', 'r')
import csv
file = "ms3-dna-mammuthus.txt"

def load_file(file):
    myfile = open( file )
    mydata = myfile.readlines()
    myfile.close()
    return mydata

def assemble_genome2(dna_list=mydata):
  mammoth = dna_list[0]
  while len(mammoth) < 16770:
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
    print(mammoth)
    print(len(mammoth))
  return mammoth

    


