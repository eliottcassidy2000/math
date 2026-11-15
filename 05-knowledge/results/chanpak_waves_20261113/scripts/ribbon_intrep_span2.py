import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import numpy as np, itertools, sys
sys.path.insert(0, RIB)
from ribbon import letters, count_with_matrices, transfer_system
R=os.path.join(RIB,'span2_intrep.txt')
tok=open(R).read().split()
pos=0
def nxt():
    global pos; v=tok[pos]; pos+=1; return int(v)
dim=nxt(); nl=nxt()
init=np.array([nxt() for _ in range(dim)],dtype=object)
fin=np.array([nxt() for _ in range(dim)],dtype=object)
MATS={}
for _ in range(nl):
    li=nxt()
    M=np.array([[nxt() for _ in range(dim)] for _ in range(dim)],dtype=object)
    MATS[li]=M
L=letters(2)
print('letters order',[sorted(a) for a in L])
SYS=transfer_system(2)
def Hrep(word):
    y=init.copy()
    for a in word: y=MATS[L.index(a)].dot(y)
    return int(fin.dot(y))
# test (letters index: as in letters(2))
import random
rng=random.Random(1)
ok=True
for t in range(200):
    n=rng.randint(1,12)
    word=[frozenset()]+[L[rng.randrange(4)] for _ in range(n-1)]
    word=[frozenset(s for s in a if s<q+1) if q>0 else frozenset() for q,a in enumerate(word)]
    if Hrep(word)!=count_with_matrices(2,word,SYS): ok=False; print('mismatch',word); break
print('rep ok',ok)
