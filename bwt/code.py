def gener_BWT_transform(text):
    """Burrows-Wheeler transformation
    """
    bwt_text = []
    for letter in text :
        text = text[-1:] + text[:-1]
        bwt_text.append(text)
    return "".join([i[-1] for i in sorted(bwt_text)])
 
 
def gener_P_N(transform):
    """N & P matrices for Burrows-Wheeler transformation
    """
    P = []
    dico = {"$" : 0, "A" : 0, "C" : 0, "G" : 0, "T" : 0}
    for lettre in transform :
        P.append(dico[lettre])
        dico[lettre] += 1
    N = {}
    count = 0
    for k in dico.keys() :
        N[k] = count
        count += dico[k]
    return P, N
 
def gener_retro_bwt(transform, P, N) :
    i = len(transform) - 1
    origin = "$"
    j = 0
    while i >= 1 :
        origin = transform[j] + origin
        j = N[transform[j]] + P[j]
        i -= 1
    return origin
 
def pos_BWT(transform, P, N):
    i = len(transform) - 2
    pos = [-1] * len(transform)
    j = 0
    while transform[j] != "$" :
        pos[j] = i
        j = N[transform[j]] + P[j]
        i -= 1
    return pos
 
def FM_index(transform):
    dict_list = [{"$" : 0, "A" : 0, "C" : 0, "G" : 0, "T" : 0}]
    for letter in transform :
        index_dict = dict_list[-1]
        index_dict[letter] += 1
        dict_list.append(index_dict.copy())
    return dict_list
 
def match(read, N, FM):
    nt = read[-1]
    d = N[nt] #debut
    f = N[nt] + FM[-1][nt] -1
    for nt in read[-2::-1]:
        d = N[nt] + FM[d][nt]
        f = N[nt] + FM[f+1][nt] -1
        print(d,f)
        if f < d :
            return (-1,-1)
    return d, f

if __name__ == '__main__':
    # complexity()
    text = "ATATCGT$"
    read = "CGT"
    final = gener_BWT_transform(text)
    P, N = gener_P_N(final)
    origin = gener_retro_bwt(final, P, N)
    dict_list = FM_index(final)
    d,f = match(read, N, dict_list)
    print(d,f)