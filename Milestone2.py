def find_splice(dna_motif, dna):
    test_value = 0
    pos_list = []
    for i in range(0,len(dna)):
        if dna[i] == dna_motif[0]:
            pos_list.append(i)
            test_value = i
            break
    for i in range(0,len(dna)):
        if dna[i] == dna_motif[1]:
            if i > test_value:
                pos_list.append(i)
                test_value = i
                break
    for i in range(0,len(dna)):
        if dna[i] == dna_motif[2]:
            if i > test_value:
                pos_list.append(i)
                test_value = i
                break
                
    return pos_list
print(find_splice("GTA","ACGACATCACGTGACG"))
if find_splice("GTA","ACGACATCACGTGACG") == [2, 6, 8]:
    print('True')
else:
    print('False')
    

def shared_motif(dna_list):   
    counter = 0  
    answer = ''         
    for order in range(len(dna_list)):
        for i in range(len(dna_list[order])):
            string = (dna_list[order][0:len(dna_list[order])-i])
            #print(string)
            for length1 in range(len(dna_list)):
                if string in dna_list[length1]:
                    #print('true')
                    counter += 1
                else:
                    #print('false')
                    counter = 0
                if counter == len(dna_list):
                    if len(string) > len(answer):
                        answer = string
    return answer

print(shared_motif(["GATTACA", "TAGACCA", "ATACA"]))
if shared_motif(["GATTACA", "TAGACCA", "ATACA"]) == "TA":
    print('True')
else:
    print('False')

print(shared_motif(["ATATACA", "ATACAGA", "GGTATACA"]))
if shared_motif(["ATATACA", "ATACAGA", "GGTATACA"]) == "ATACA":
    print('True')
else:
    print('False')
    

    
"Written by Zyaire Howard"
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

    def transform(a, cover):
        if cover == (1<<length) - 1:
            return dna_list[a]
        result = 'stop' * 320
        for b in range(length):
            if cover & (1<<b) == 0:
                z = combine[a][b]
                string = transform(b, cover | (1<<b))
                if len(dna_list[a] + string[z:]) < len(result):
                    result = dna_list[a] + string[z:]
        return result
    return min([transform(a, 1<<a) for a in range(length)], key=len)


print(assemble_genome(["ATTAGACCTG", "CCTGCCGGAA", "AGACCTGCCG", "GCCGGAATAC"]))
