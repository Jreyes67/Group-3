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


def assemble_genome2(dna_list):
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
    

f = open('ms3-dna-100.txt', 'r')

print(assemble_genome2(f))
    


