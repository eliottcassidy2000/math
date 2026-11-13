# exact strong/decomposable profile collisions within the C=1 slice, orders 6,8,10,12
import sys
from math import comb
def qvals(F,T):  # Q(t)=sum_j F_j binom(t,j), t=1..T
    return tuple(sum(F[j-1]*comb(t,j) for j in range(1,len(F)+1)) for t in range(1,T+1))
blocks={1:[('TT1',(1,))]}
for a in (3,5,7,9):
    blocks[a]=[]
    for line in open('prof_%d.txt'%a):
        p=line.split()
        if int(p[2])==1: blocks[a].append((p[1],tuple(map(int,p[7:]))))
print({a:len(v) for a,v in blocks.items()})
def c1list(n):
    L=[]
    for line in open('c1p_%d.txt'%n):
        left,right=line.split('|'); p=left.split(); L.append((tuple(map(int,right.split())),tuple(map(int,p[2:]))))
    return L
def wouldY(F):
    n=len(F); G=[]; prev=0
    for j in range(1,n+1):
        if F[j-1]%j: return None
        g=F[j-1]//j-prev
        if g<0: return None
        G.append(g); prev=g
    if G[-1]!=0: return None
    return tuple(G[:-1])
for n in (6,8,10,12):
    T=n
    vals={}
    for a in range(1,n,2):
        b=n-a
        if a>b: continue
        if a not in blocks or b not in blocks: continue
        for c1,F1 in blocks[a]:
            v1=qvals(F1,T)
            for c2,F2 in blocks[b]:
                v2=qvals(F2,T)
                key=tuple(x*y for x,y in zip(v1,v2))
                vals.setdefault(key,[]).append((a,c1,b,c2))
    L=c1list(n); hits=[]; wy=0
    for rows,F in L:
        key=qvals(F,T)
        if key in vals: hits.append((rows,F[0],vals[key]))
        if wouldY(F) is not None: wy+=1
    print('order',n,': C=1 classes',len(L),' collisions with a two-block (odd,odd) ordinal sum:',len(hits),'  would-be-Y integral:',wy)
    for rows,H,m in hits: print('   H',H,'rows',rows,'partners',m[:6])
