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

def get_edges(dict):
        key_list = dict.keys()
        list = []
        for i in key_list:
            list.append(i)
        list2 = []
        for i in range(0,len(list)):
            for j in range(i+1,len(list)):

                if(dict[list[i]][-3:] == dict[list[j]][:3] or dict[list[i]][:3] == dict[list[j]][-3:]):
                    list2.append((list[i],list[j]))          
        return list2
    

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

def reverse_complement(dna):
    complement = ''
    for character in dna:
        if character == 'T':
            complement = complement + 'A'
        elif character == 'A':
            complement = complement + 'T'
        elif character == 'G':
            complement = complement + 'C'
        elif character == 'C':
            complement = complement + 'G'
    return complement[::-1]

def rev_palindrome(dna):
    result = []
    for a in range(len(dna)-4):
        for b in range(a+3,min(len(dna),a+12)):
            string = dna[a:b+1]
            if reverse_complement(dna[a:b+1]) == string:
                result.append((a,b-a+1))

    return result
print(rev_palindrome("TCAATGCATGCGGGTCTATATGCAT"))

#Written by Adesh
import math 
def random_genome(dna, gc_content):
    dna = dna.upper()
    cg = len(dna.replace('A','').replace('T',''))
    at = len(dna.replace('C','').replace('G',''))
    result= []
    for i in range (0, len(gc_content)):
        prob = cg * math.log10(float(gc_content[i])/ 2) + at * math.log10((1- float(gc_content[i]))/ 2)
        result.append(round(prob, 3))
    return result
dna="ACGATACAA"
gc_content=[0.129,0.287,0.423,0.476,0.641,0.742,0.783]
a = random_genome(dna,gc_content)
print(a)

def fact(n):
    if n==1 or n==0:
        return 1
    return n*fact(n-1)


def perfect_match(rna):
    d= {'A':0,'C':0,'G':0,'U':0}
    for i in rna:
        d[i]+=1
    if d['A']==d['U'] and d['C']==d['G']:
        return fact(d['A'])*fact(d['C'])
    else:
        return "Perfect Match doesn't exist"

def main():
    rna = input("Enter RNA String : \n")
    total_match = perfect_match(rna)
    if type(total_match)==int:
        print("Total Number of perfect matching : ",total_match)
    else:
        print(total_match)
if __name__=="__main__":
    main()


