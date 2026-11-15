import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import sys, itertools, pickle, random
sys.path.insert(0,RIB)
from ribbon import transfer_system, count_with_matrices
from sympy import Matrix, Rational
W=3; SYS=transfer_system(W)
units=sys.argv[1].split(',')
def lets(s): return [frozenset(int(c) for c in t if c in '123') for t in s.split('.')]
U=[lets(u) for u in units]
def word(seq): 
    w=[frozenset()]
    for i in seq: w+=U[i]
    return [frozenset(x for x in a if x<=q) for q,a in enumerate(w)]
cache={}
def H(seq):
    seq=tuple(seq)
    if seq not in cache: cache[seq]=count_with_matrices(W,word(seq),SYS)
    return cache[seq]
K=int(sys.argv[2]) if len(sys.argv)>2 else 2
words=[()]
for k in range(1,K+1): words+=list(itertools.product(range(len(U)),repeat=k))
# modular rank to find size, then greedy exact selection
p=1000003
rows=[]; cols=[]
# greedy: build basis of suffixes (columns) and prefixes (rows) with exact Fractions via sympy on small matrices
def rank_of(Mat): return Matrix(Mat).rank()
J=[()]  # empty suffix forced
I=[()]
def Z(I,J): return [[H(a+b) for b in J] for a in I]
r=1
improved=True
while improved:
    improved=False
    for s in words:
        if s in J: continue
        for a in words:
            if a in I: continue
            M=Z(I+[a],J+[s])
            if rank_of(M)==r+1:
                I.append(a); J.append(s); r+=1; improved=True; break
        if improved: break
print('Hankel rank',r,'prefixes',I,'suffixes',J)
Zm=Matrix(Z(I,J)).T   # columns z(P_i) = (H(P_i S_j))_j
Zi=Zm.inv()
red=[]
for ui in range(len(U)):
    Zu=Matrix([[H(a+(ui,)+b) for b in J] for a in I]).T
    red.append(Zu*Zi)
z0=Matrix([[H(b)] for b in J]); fin=Matrix([[1 if b==() else 0 for b in J]])
# validate on random unit words
rng=random.Random(3); bad=0
for t in range(300):
    seq=tuple(rng.randrange(len(U)) for _ in range(rng.randint(0,6)))
    z=z0
    for i in seq: z=red[i]*z
    if (fin*z)[0]!=H(seq): bad+=1
print('validation mismatches',bad)
pickle.dump({'units':units,'red':red,'z0':z0,'fin':fin},open(sys.argv[3],'wb'))
