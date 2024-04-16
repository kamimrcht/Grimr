import pickle

def del_inclusions(f:dict[int,list[set[int]]]):
    E_list = []
    keys = sorted(f.keys())
    for i in range(len(keys)):
        n = keys[i]
        while len(f[n]) > 0:
            E = f[n].pop()
            for j in range(i, len(keys)):
                m = keys[j]
                new_list = []
                while len(f[m]) > 0:
                    Ep = f[m].pop()
                    if not E.issubset(Ep):
                        new_list.append(Ep)
                f[m] = new_list
            E_list.append(list(E))
    return E_list

def simplify_query(labels_dict,All_list,Any_list,NotAll_list,NotAny_list):
    A = set()
    for label in All_list:
        A.update(labels_dict[label])
    D = set()
    for label in NotAny_list:
        D.update(labels_dict[label])

    if len(A.intersection(D))!=0:
        return [],[],[],[]
    else:
        f = dict()
        for label in Any_list:
            B = labels_dict[label].copy()
            if len(B.intersection(A))==0:
                B.difference_update(D)
                if len(B)>0:
                    if len(B) not in f.keys():
                        f[len(B)]=[]
                    f[len(B)].append(B)

        Bcal = del_inclusions(f)

        f = dict()
        for label in NotAll_list:
            C = labels_dict[label].copy()
            if len(C.intersection(D))==0:
                C.difference_update(A)
                if len(C)>0:
                    if len(C) not in f.keys():
                        f[len(C)]=[]
                    f[len(C)].append(C)

        Ccal = del_inclusions(f)

        return list(A),Bcal,Ccal,list(D)


with open("../test_files/metadata.pickle", "rb") as input_file:
    labels_dict = pickle.load(input_file)

# QUERY

# Each letter below corresponds to a label in `labels_dict`, i.e. a list of files

All_list = ['J','K']
Any_list = ['I','H']
NotAll_list = ['NOT-F']
NotAny_list = ['F']

A,B,C,D = simplify_query(labels_dict,All_list,Any_list,NotAll_list,NotAny_list)

with open("../test_files/query.txt", "a") as fichier:
    fichier.write("A"+"\t"+"ALL"+"\t"+str(A)+"\n")
    fichier.write("B"+"\t"+"ANY"+"\t"+str(B)+"\n")
    fichier.write("C"+"\t"+"NOT-ALL"+"\t"+str(C)+"\n")
    fichier.write("D"+"\t"+"NOT-ANY"+"\t"+str(D)+"\n")

# f = {'1':[set([1])], '2':[set([1,2]),set([2,3])], '3':[set([2,3,4]),set([2,4,5]),set([1,3,4])]}
#
# print(del_inclusions(f))