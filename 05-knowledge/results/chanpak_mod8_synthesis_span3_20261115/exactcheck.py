import numpy as np
import sys as _s
tok=open(_s.argv[1] if len(_s.argv)>1 else '../../../04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt').read().split(); pos=0
def nx():
    global pos; v=int(tok[pos]); pos+=1; return v
dim=nx(); nl=nx(); init=[nx() for _ in range(dim)]; fin=[nx() for _ in range(dim)]
M={}
for _ in range(nl):
    li=nx(); M[li]=[[nx() for _ in range(dim)] for _ in range(dim)]
PR=[2147483629,2147483587]
def mm(A,B,p):
    A=np.array(A,dtype=object); B=np.array(B,dtype=object); return (A.dot(B))%p
def rank_mod(A,p):
    A=[[int(x)%p for x in row] for row in A]; r=0; rows=len(A); cols=len(A[0])
    for c in range(cols):
        piv=next((i for i in range(r,rows) if A[i][c]),None)
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]; inv=pow(A[r][c],-1,p); A[r]=[x*inv%p for x in A[r]]
        for i in range(rows):
            if i!=r and A[i][c]:
                f=A[i][c]; A[i]=[(x-f*y)%p for x,y in zip(A[i],A[r])]
        r+=1
    return r
def word_letters(w):
    out=[]
    for u in w.split('|'):
        out+= [0]*(len(u)-1)+[int(u[-1])]
    return out
def Hval(letters):
    v=np.array(init,dtype=object)
    for a in letters: v=np.array(M[a],dtype=object).dot(v)
    return int(np.array(fin,dtype=object).dot(v))
for w in ['05','2|7|05','2|05','1|2|4','1|3|4','1|04','3|04','1|1|7|4']:
    L=word_letters(w)
    ranks=[]
    for p in PR:
        P=np.identity(dim,dtype=object)
        for a in L: P=mm(M[a],P,p)
        Q=np.identity(dim,dtype=object)
        for _ in range(64): Q=mm(P,Q,p)
        ranks.append(rank_mod(Q.tolist(),p))
    hs=[Hval(L*k) for k in (1,2,3,4,5)]
    print('%-10s letters %-12s rank of P^64 mod two primes %s; H(word^k), k=1..5: %s'%(w,L,ranks,hs))
