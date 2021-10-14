import itertools
import random
import matplotlib.pyplot as plt
import time
import math

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
    print(text)
    reverse_bwt(indexing, count_letter, final)
    return final, count_letter, indexing

def reverse_bwt(indexing, count_letter, final):
    original = []
    new_index = 0
    len_seq = len(final) - 1
    count = 1
    position = [len_seq]*(len_seq+1)
    print(final)
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
    print(position)
    return original

def seq_ran(n):
    nuc = ["A","T","C","G"]
    seq = ''
    for i in range(n):
        seq += random.choice(nuc)
    seq += "$"
    return seq

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

if __name__ == '__main__':
    # complexity()
    text = "ATATCGT$"
    final, count_letter, indexing = transform(text)