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

ms3-dna-50 = [TCCCGATGTAGTCCCCACCG,
CCCCACCGCGACACAGTCAGT,
AGTAGACGCACTCGCCGTCCCGATG]

def assemble_genome2(dna_list):
    first_string = dna_list[i]
    beginning8 = first_string[-1:-8]
    
