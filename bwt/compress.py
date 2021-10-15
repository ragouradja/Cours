import itertools
import random
import matplotlib.pyplot as plt
import time
import math
import copy

def transform(text):
    nb = len(text)
    perm = []
    for i in range(nb):
        pre = text[-i:]
        su = text[:-i]
        perm.append([pre+""+su])
    perm.sort()
    final = ""
    first_col = ""
    indexing = {}
    count_letter = []
    for i,liste in enumerate(perm):
        first_letter = "".join(liste)[0]
        last_letter = "".join(liste)[-1]
        if first_letter not in indexing:
            indexing[first_letter] = i

        count_letter += [final.count(last_letter)]
        
        final += last_letter
        first_col += first_letter
    return final, count_letter, indexing, first_col

def reverse_bwt(indexing, count_letter, final):
    original = []
    new_index = 0
    len_seq = len(final) - 1
    count = 0
    position = [0]*(len_seq+1)
    while 1:
        letter = final[new_index]
        if letter == "$":
            original.reverse()
            original += ["$"]
            break
        original.append(letter)
        position[new_index] = (len_seq - count)

        new_index = indexing[letter] + count_letter[new_index]
        count += 1
    return original,position

def seq_ran(n):
    nuc = ["A","T","C","G"]
    seq = ''
    for i in range(n):
        seq += random.choice(nuc)
    seq += "$"
    return seq

def FM_dico(final):
    dico_occ = {"$" : 0, "A" : 0, "C" : 0, "G": 0, "T": 0}
    list_dico = []
    list_dico.append(copy.deepcopy(dico_occ))
    for i in range(1,len(final) + 1):
        dico_occ[final[i-1]] += 1
        list_dico.append(copy.deepcopy(dico_occ))
    return list_dico
 
def search_read(read, final, list_dico, indexing, count_letter, first_col, position):
    d,f = (1,2)
    final_position = []
    for i in range(1, len(read)):
        first = read[-i-1]
        second = read[-i]
        sub_read = second +""+first
        d,f = match(indexing, list_dico, sub_read)
        for j in range(d, f + 1):
            if final[j] == first:
                if first_col[j] == second:
                    pos_first = position[j]
                    final_position += [pos_first, pos_first+1]

    print(final_position)

def match(indexing, list_dico, read):

    nt = read[-1]
    d = indexing[nt]
    f = indexing[nt] + list_dico[-1][nt] - 1
    for n in read[-2::-1]:
        d = indexing[n] + list_dico[d][n]
        f = indexing[n] + list_dico[f + 1][n] - 1
        if f < d:
            return (-1,-1)
    return d,f


def complexity():
    kb = 1000
    end_bwt = []
    end_reverse = []
    x = []
    for i in range(kb,kb*50,kb):
        x.append(i)
        sequence = seq_ran(i)
        start_bwt = time.time()
        final, count_letter, indexing = transform(sequence)
        end_bwt.append((time.time() - start_bwt))
        start_reverse = time.time()
        reverse_bwt(indexing, count_letter, final)
        end_reverse.append((time.time() - start_reverse))

    plt.figure(figsize=(20,10))
    plt.title("Time complexity")
    plt.plot(x, end_bwt, label = "bwt")
    plt.plot(x, end_reverse, label = "reverse")
    plt.legend()
    plt.xlabel("Taille")
    plt.ylabel("Time (s)")
    plt.show()
    
def pos(start, end,  pos_list):
    position = []
    for i in range(start, end+1):
        position.append(pos_list[i] + 1)
    return position

def read_fasta(fasta_file, reference, output):
    seq = ""
    to_write = []
    ref_genome = get_ref(reference)
    ref_genome += "$"
    read = ""
    zeros = 0
    one = 0
    two = 0
    with open(fasta_file, "r") as filin, open(output, "w") as filout:
        for line in filin:
            if line.startswith(">"):
                if read != "":
                    final, count_letter, indexing, first_col = transform(ref_genome)
                    origin, position = reverse_bwt(indexing, count_letter, final)

                    list_dico = FM_dico(final)
                    
                    d,f = match(indexing, list_dico, read)
                    if (d,f) == (-1,-1):
                        pos_in_seq = 0
                        zeros += 1
                    else:
                        pos_in_seq = pos(d,f,position)
                        if len(pos_in_seq) == 1:
                            one += 1
                        elif len(pos_in_seq) == 2:
                            two += 1
                    to_write = "{} \t {} \t {}\n".format(read, read_name, pos_in_seq)
                    filout.write(to_write)
                read_name = line[1:].strip()
                read = ""
            else:
                read += line.strip()

    print("0 : {} ; 1 : {}; 2 : {}".format(zeros, one, two))
def get_ref(fasta_file):
    seq = ""
    with open(fasta_file) as filin:
        for line in filin:
            if not line.startswith(">"):
                seq += line.strip()
    return seq
if __name__ == '__main__':
    path1 = "masterBI/reads_1000_10_patient_6.fna"
    path2 = "masterBI/reads_10000_30_patient_6.fna"
    
    ref = "masterBI/NC_045512-N.fna"
    read_fasta(path2,ref,"output2.txt")