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
    
