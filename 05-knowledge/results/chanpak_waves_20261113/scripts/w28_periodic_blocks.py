import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import numpy as np, itertools, sys
sys.path.insert(0, RIB)
from ribbon import letters, tournament_from_word, is_strong
w=int(sys.argv[1]); lmax=int(sys.argv[2])
R=os.path.join(RIB,f'span{w}_intrep.txt')
tok=open(R).read().split(); pos=0
def nxt():
    global pos; v=tok[pos]; pos+=1; return int(v)
dim=nxt(); nl=nxt()
init=np.array([nxt() for _ in range(dim)],dtype=float)
fin=np.array([nxt() for _ in range(dim)],dtype=float)
MATS={}
for _ in range(nl):
    li=nxt(); MATS[li]=np.array([[nxt() for _ in range(dim)] for _ in range(dim)],dtype=float)
L=letters(w)
MASK=[sum(1<<(s-1) for s in a) for a in L]
def nm(a): return ''.join(map(str,sorted(a))) or '0'
def canon(u): return min(u[i:]+u[:i] for i in range(len(u)))
seen=set(); rows=[]
for l in range(1,lmax+1):
    for u in itertools.product(range(nl),repeat=l):
        c=canon(u)
        if c in seen: continue
        if not all(c!=c[k:]+c[:k] for k in range(1,l) if l%k==0): continue
        seen.add(c)
        M=np.identity(dim)
        for li in c: M=MATS[MASK[li]]@M
        ev=sorted(abs(np.linalg.eigvals(M)),reverse=True)
        rho=ev[0]
        if rho<=1e-9 or rho**(1.0/l)>1.80: continue
        # strongness of a long periodic word (truncate spans at start)
        k=max(6, (9*w)//l+6)
        word=[]
        for q in range(k*l):
            a=L[c[q%l]]; word.append(frozenset(s for s in a if s<=q))
        adj=tournament_from_word(w,word); n=len(adj)
        # middle third in one strong component?
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
        mid=n//2; a=reach(adj,mid); b=reach(radj,mid); comp=a&b
        lo=n//3; hi=2*n//3
        if any(not (comp>>t)&1 for t in range(lo,hi)): continue
        sub=ev[1]
        rows.append((rho**(1.0/l), sub/rho if rho>0 else 0, l, '.'.join(nm(L[i]) for i in c), round(rho,6)))
rows.sort()
print('strong periodic blocks, lowest growth: (growth/letter, sub/dominant ratio, len, block, rho)')
for r in rows[:30]: print(r)
print('lowest NON-bottleneck (sub/dominant > 1e-6):')
cnt=0
for r in rows:
    if r[1]>1e-6:
        print(r); cnt+=1
        if cnt>=25: break
