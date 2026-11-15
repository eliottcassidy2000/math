import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import numpy as np, itertools, sys
sys.path.insert(0, RIB)
from ribbon import letters, tournament_from_word
w=int(sys.argv[1]); lmax=int(sys.argv[2]); fam=sys.argv[3] if len(sys.argv)>3 else 'all'
R=os.path.join(RIB,f'span{w}_intrep.txt')
tok=open(R).read().split(); pos=0
def nxt():
    global pos; v=tok[pos]; pos+=1; return int(v)
dim=nxt(); nl=nxt()
init=[nxt() for _ in range(dim)]; fin=[nxt() for _ in range(dim)]
MI={}
for _ in range(nl):
    li=nxt(); MI[li]=np.array([[nxt() for _ in range(dim)] for _ in range(dim)],dtype=np.int64)
MF={k:v.astype(float) for k,v in MI.items()}
L=letters(w); MASK=[sum(1<<(s-1) for s in a) for a in L]
def nm(a): return ''.join(map(str,sorted(a))) or '0'
NAMES=[nm(a) for a in L]
PR=[33554393, 33554383]
def rank_mod(A,p):
    A=A.copy()%p; r=0; rows,cols=A.shape
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]%p: piv=i;break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        inv=pow(int(A[r,c]),p-2,p)
        A[r]=(A[r]*inv)%p
        for i in range(rows):
            if i!=r and A[i,c]:
                A[i]=(A[i]-A[i,c]*A[r])%p
        r+=1
        if r==rows: break
    return r
def matmod(u,p,power):
    M=np.identity(dim,dtype=np.int64)
    for li in u: M=(MI[MASK[li]]@M)%p
    P=np.identity(dim,dtype=np.int64)
    for _ in range(power): P=(M@P)%p
    return P
def bulk_strong(c):
    l=len(c); k=max(6,30//l+4); word=[]
    for q in range(k*l):
        a=L[c[q%l]]; word.append(frozenset(s for s in a if s<=q))
    adj=tournament_from_word(w,word); n=len(adj)
    def reach(adjm,src):
        seen=1<<src; fr=1<<src
        while fr:
            v=(fr&-fr).bit_length()-1; fr&=fr-1
            new=adjm[v]&~seen; seen|=new; fr|=new
        return seen
    radj=[0]*n
    for p in range(n):
        o=adj[p]
        while o:
            qq=(o&-o).bit_length()-1; o&=o-1; radj[qq]|=1<<p
    mid=n//2; comp=reach(adj,mid)&reach(radj,mid)
    return all((comp>>t)&1 for t in range(n//3,2*n//3))
def canon(u): return min(u[i:]+u[:i] for i in range(len(u)))
cands=[]
if fam=='all':
    gen=(u for l in range(1,lmax+1) for u in itertools.product(range(nl),repeat=l))
else:
    floors=['0.1.3','0.3.1','0.13.2','0.1.123','0.3.12','0.12.3','1.2.3','1.12.3','0.1.23']
    def g():
        for B in floors:
            b=[NAMES.index(t) for t in B.split('.')]
            for j in (1,2):
                for sl in range(1,lmax+1):
                    for s in itertools.product(range(nl),repeat=sl): yield tuple(b*j+list(s))
    gen=g()
seen=set()
for u in gen:
    c=canon(u)
    if c in seen: continue
    seen.add(c); l=len(c)
    if not all(c!=c[k:]+c[:k] for k in range(1,l) if l%k==0): continue
    M=np.identity(dim)
    for li in c: M=MF[MASK[li]]@M
    rho=max(abs(np.linalg.eigvals(M)))
    if rho<1.5: continue
    g=rho**(1.0/l)
    if g>=1.7700: continue
    if not bulk_strong(list(c)): continue
    cands.append((g,c))
print('candidates below 1.7700 (bulk-strong):',len(cands))
nonb=[]
for g,c in sorted(cands):
    r=[rank_mod(matmod(c,p,dim),p) for p in PR]
    if max(r)!=1: nonb.append((g,'.'.join(NAMES[i] for i in c),r))
print('non-bottleneck among them (rank of M^dim mod two primes != 1):',len(nonb))
for t in nonb[:20]: print(t)
