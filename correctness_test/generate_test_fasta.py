import sys
from copy import copy

def int_to_bin(num):
    binary = []
    while num != 0:
        bit = num % 2
        binary.append(bit)
        num = num // 2
    binary.reverse()
    return binary

def int_to_kmer(num,k):
    binary = int_to_bin(num)
    while len(binary)<2*k:
        binary.insert(0,0)

    dic = {(0,0): "A", (0,1):"C", (1,0):"G", (1,1): "T"}

    s=''
    for i in range(k):
        bits = tuple(binary[2*i:2*i+2])
        s+=dic[bits]

    assert len(s)==k
    return s

def kmer_set_to_fasta(F,filename):
    f = open(filename, "w")
    i=0
    for kmer in F:
        f.write("> kmer "+str(i+1)+"\n")
        f.write(kmer+"\n")
        i+=1

k=31
N = 100 # nb de k-mers par bloc
M = 22 # nb de blocs

kmer_list = []
nmax = 2**(2*k)-1

assert N*M <= nmax

for n in range(N*M):
    sys.stdout.write("\rprogression: %i" % n)
    kmer_list.append(int_to_kmer(n,k))
print("\n")

blocs = {}
for i in range(M):
    blocs[i]=[]
    for j in range(N):
        kmer = kmer_list.pop()
        blocs[i].append(kmer)

a = set()
for i in [1,2,3,4,5,6,7,8]:
    a.update(blocs[i])

b = set()
for i in [2,6,7,8,12,13,14,15]:
    b.update(blocs[i])

c = set()
for i in [4,5,6,7,10,11,12,13]:
    c.update(blocs[i])

d = set()
for i in [3,4,7,8,9,10,12,14]:
    d.update(blocs[i])

# F1 et F2

F1 = copy(a)
F1.update(blocs[16])

F2 = copy(a)
F2.update(blocs[17])

# F3 et F4

split_d = copy(d)
F3 =set()
F4 = set()
while len(split_d)>0:
    F3.add(split_d.pop())
    F4.add(split_d.pop())

# F5, F6, F7 et F8

F5 = set()
F6 = set()
F7 = set()
F8 = set()

split_b = copy(b)
while len(split_b)>0:
    elt = split_b.pop()
    F6.add(elt)
    F7.add(elt)
    elt = split_b.pop()
    F5.add(elt)
    F8.add(elt)

split_c = copy(c)
while len(split_c)>0:
    elt = split_c.pop()
    F5.add(elt)
    F6.add(elt)
    elt = split_c.pop()
    F7.add(elt)
    F8.add(elt)

F5.update(blocs[18])
F6.update(blocs[19])
F7.update(blocs[20])
F8.update(blocs[21])

print(len(F1),len(F2),len(F3),len(F4),len(F5),len(F6),len(F7),len(F8))

kmer_set_to_fasta(F1,"F1.fa")
kmer_set_to_fasta(F2,"F2.fa")
kmer_set_to_fasta(F3,"F3.fa")
kmer_set_to_fasta(F4,"F4.fa")
kmer_set_to_fasta(F5,"F5.fa")
kmer_set_to_fasta(F6,"F6.fa")
kmer_set_to_fasta(F7,"F7.fa")
kmer_set_to_fasta(F8,"F8.fa")

# A	ALL	[1,2]
# B	ANY	[[5,6],[7,8]]
# C	NOT-ALL	[[5,6],[7,8]]
# D	NOT-ANY	[3,4]

univers = F1 | F2 | F3 | F4 | F5 | F6 | F7 | F8

assert len(univers)==N*(M-1)

# QUERY 1 : A

query = F1 & F2
assert query == a
kmer_set_to_fasta(a,"expected_1.fa")

# QUERY 2 : B

query = (F5 | F6) & (F7 | F8)
assert query == b
kmer_set_to_fasta(b,"expected_2.fa")

# QUERY 3 : C

query = (F5 & F6) | (F7 & F8)
assert query == c
kmer_set_to_fasta(c,"expected_3.fa")

# QUERY 4 : D

query = F3 | F4
assert query == d
kmer_set_to_fasta(d,"expected_4.fa")

# QUERY 5 : (A & B)

query = a & b
sol = set(blocs[2]+blocs[6]+blocs[7]+blocs[8])
assert query == sol
kmer_set_to_fasta(sol,"expected_5.fa")

# QUERY 6 : U \ (C | D)

query = univers - (c | d)
sol = set(blocs[1]+blocs[2]+blocs[15]+blocs[16]+blocs[17]+blocs[18]+blocs[19]+blocs[20]+blocs[21])
assert query == sol
kmer_set_to_fasta(sol,"expected_6.fa")

# QUERY 7 : A \ C

query = a - c
sol = set(blocs[1]+blocs[2]+blocs[3]+blocs[8])
assert query == sol
kmer_set_to_fasta(sol,"expected_7.fa")

# QUERY 8 : A \ D

query = a - d
sol = set(blocs[1]+blocs[2]+blocs[5]+blocs[6])
assert query == sol
kmer_set_to_fasta(sol,"expected_8.fa")

# QUERY 9 : B \ C

query = b - c
sol = set(blocs[2]+blocs[8]+blocs[14]+blocs[15])
assert query == sol
kmer_set_to_fasta(sol,"expected_9.fa")

# QUERY 10 : B \ D

query = b - d
sol = set(blocs[2]+blocs[6]+blocs[13]+blocs[15])
assert query == sol
kmer_set_to_fasta(sol,"expected_10.fa")

# QUERY 11 : A \ (C | D)

query = a - (c | d)
sol = set(blocs[1]+blocs[2])
assert query == sol
kmer_set_to_fasta(sol,"expected_11.fa")

# QUERY 12 : B \ (C | D)

query = b - (c | d)
sol = set(blocs[2]+blocs[15])
assert query == sol
kmer_set_to_fasta(sol,"expected_12.fa")

# QUERY 13 : (A & B) \ C

query = (a & b) - c
sol = set(blocs[2]+blocs[8])
assert query == sol
kmer_set_to_fasta(sol,"expected_13.fa")

# QUERY 14 : (A & B) \ D

query = (a & b) - d
sol = set(blocs[2]+blocs[6])
assert query == sol
kmer_set_to_fasta(sol,"expected_14.fa")

# QUERY 15 : (A & B) \ (C | D)

query = (a & b) - (c | d)
sol = set(blocs[2])
assert query == sol
kmer_set_to_fasta(sol,"expected_15.fa")
