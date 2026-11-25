"""Reproduce the exact positive rational and signed-lattice certificates."""
import sys
from collections import defaultdict
from math import gcd

F=list(map(int,sys.stdin.read().split()))
n=F[0];p=1;rows=[]
for h in range(n):
    old,k=F[p:p+2];p+=2
    X=[tuple(F[p+4*j:p+4*j+4]) for j in range(k)]
    p+=4*k
    rows.append(X)
assert p==len(F) and n==53
assert sum(map(len,rows))==1848  # after the marker row is restricted, before removing its hits
assert len(rows[0])==1 and rows[0][0][0]==7
mark=rows[0][0][1:]
cols=[]
for j,R in enumerate(rows[1:]):
    for y,*t in R:
        if all(a!=b for a,b in zip(t,mark)):
            cols.append((j,*t))
assert len(cols)==1741
ids=[dict(zip(sorted({C[i] for C in cols}),range(i*52,(i+1)*52))) for i in range(4)]
C=[tuple(ids[i][x] for i,x in enumerate(col)) for col in cols]
L=[[] for i in range(208)]
for j,col in enumerate(C):
    for i in col:L[i].append(j)
S=1<<280;w=[S//35]*len(C)
for _ in range(200):
    for row in L:
        s=sum(w[j] for j in row)
        for j in row:w[j]=(w[j]*S+s//2)//s
err=max(abs(sum(w[j] for j in row)-S) for row in L)
assert err==15 and min(w)>4*10**82
assert min(w)>205*err*(1<<204)
print("Positive fractional cover: rounded error",err,"; min rounded weight",min(w))

# The equal-total lattice projects bijectively onto 205 free integer coordinates
# by deleting one row of each of the three central margin sets.
U=[{x:1 for x in row if x not in {103,155,207}} for row in C]
def eliminate(v, u, k, q):
    for a,b in u.items():
        m=v.get(a,0)-q*b
        if m:v[a]=m
        elif a in v:del v[a]

pivots=0;gcd_turns=0
while any(U):
    # Restrict each pivot selection to live columns. Use small-fill unit pivots.
    occ=defaultdict(list)
    for j,col in enumerate(U):
        for i in col:occ[i].append(j)
    P=min((((len(occ[i])-1)*(len(col)-1),len(col),len(occ[i]),i,j)
          for j,col in enumerate(U) for i,x in col.items() if abs(x)==1), default=None)
    if P is None:
        i=next(i for i,js in sorted(occ.items()) if gcd(*(U[j][i] for j in js))==1)
        J=list(occ[i])
        while True:
            j=min((j for j in J if i in U[j]),key=lambda j:(abs(U[j][i]),j))
            d=U[j][i]
            if abs(d)==1:break
            for a in J:
                if a!=j and i in U[a]:
                    eliminate(U[a],U[j],i,U[a][i]//d);gcd_turns+=1
    else:
        _,_,_,i,j=P
    u=U[j];d=u[i]
    for a in occ[i]:
        if a!=j and i in U[a]:
            eliminate(U[a],u,i,U[a][i]//d)
    U[j]={}
    pivots+=1
    if pivots==205:break
assert not any(U)
print("Signed integer span: unit pivots",pivots,"; Euclidean column reductions",gcd_turns)
