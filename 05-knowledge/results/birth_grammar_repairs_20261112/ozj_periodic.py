"""Finite periodic-word spectral probe in the inherited O/Z/J root grammar."""
from common import *
from itertools import permutations
from collections import Counter
import time

M={
'O':[[1,0,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]],
'Z':[[1,-1,1,0],[0,0,0,1],[0,0,2,-1],[0,0,0,1]],
'J':[[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]]}
PERMS=[(p,(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))) for p in permutations(range(4))]
def det_shift(A,x):
    from math import prod
    return sum(sign*prod((x if i==p[i] else 0)-A[i][p[i]] for i in range(4)) for p,sign in PERMS)
def step(A,label):
    return [[sum(M[label][i][k]*A[k][j] for k in range(4)) for j in range(4)] for i in range(4)]

def run(depth=9):
    start=time.time();stats=Counter();signs=Counter();ones=[];profiles=[]
    def visit(A,word):
        if word and 'J' in word:
            stats[len(word)]+=1;d=det_shift(A,1);signs[(d>0)-(d<0)]+=1
            if d==0:
                ones.append({'word':word,'matrix':A,'integer_eigenvalues':[b for b in [2,3,4,9] if det_shift(A,b)==0]})
        if len(word)==depth:return
        for label in M:visit(step(A,label),word+label)
    visit([[int(i==j) for j in range(4)] for i in range(4)],'')
    out={'scope':f'All O/Z/J words of lengths1..{depth} containing J, repeated as fixed periods; exact integer det(I-M). Matrices inherited from THM-4904. Finite spectral probe only, no all-length exclusion or coverage claim.','depth':depth,'word_counts_by_length':dict(stats),'determinant_signs':dict(signs),'eigenvalue_one_words':ones}
    save('ozj_periodic.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k!='eigenvalue_one_words'}),flush=True);print('eigenvalue1 words',len(ones),'seconds',time.time()-start,flush=True)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 9)
