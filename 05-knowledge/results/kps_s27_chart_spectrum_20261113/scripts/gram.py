# H-spectrum of the grammar G generated from a vertex by ordinal joins X=>Y and cyclic joins C3[X,Y,Z],
# computed on full slot profiles F=(F_1..F_n) (the composite's profile depends only on the children's profiles).
import sys, itertools, time
from math import comb
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 9
def vals(F,T):  # Q(t), t=0..T
    return [sum(F[j-1]*comb(t,j) for j in range(1,len(F)+1)) for t in range(T+1)]
def fromvals(v,n):
    return tuple(sum((-1)**(j-i)*comb(j,i)*v[i] for i in range(j+1)) for j in range(1,n+1))
def ordinal(F,G):
    n=len(F)+len(G); a=vals(F,n); b=vals(G,n)
    return fromvals([x*y for x,y in zip(a,b)],n)
# itinerary monomials for one slot: start s in {0,1,2}, length L>=1 -> exponent vector
def mono(s,L):
    e=[0,0,0]
    for k in range(L): e[(s+k)%3]+=1
    return tuple(e)
# C[t][(a,b,c)] = coefficient of x^a y^b z^c in W^t, W = 1 + sum_{s,L} mono, truncated at a,b,c <= NMAX
def mul(A,B,cap):
    R={}
    for (a,b,c),u in A.items():
        for (d,e,f),w in B.items():
            k=(a+d,b+e,c+f)
            if k[0]<=cap[0] and k[1]<=cap[1] and k[2]<=cap[2]: R[k]=R.get(k,0)+u*w
    return R
Wcache={}
def Wpow(cap,t):
    key=(cap,t)
    if key in Wcache: return Wcache[key]
    if t==0: R={(0,0,0):1}
    else:
        W={(0,0,0):1}
        for s in range(3):
            L=1
            while True:
                e=mono(s,L)
                if e[0]>cap[0] or e[1]>cap[1] or e[2]>cap[2]: break
                W[e]=W.get(e,0)+1; L+=1
        R=mul(Wpow(cap,t-1),W,cap)
    Wcache[key]=R; return R
def c3join(F,G,K):
    a,b,c=len(F),len(G),len(K); n=a+b+c; cap=(a,b,c)
    v=[]
    for t in range(n+1):
        P=Wpow(cap,t); s=0
        for (i,j,k),coef in P.items():
            if i==0 or j==0 or k==0: continue
            s+=coef*F[i-1]*G[j-1]*K[k-1]
        v.append(s)
    return fromvals(v,n)
P={1:{(1,)}}
t0=time.time()
for n in range(2,NMAX+1):
    S=set()
    for a in range(1,n//2+1):
        for F in P[a]:
            for G in P[n-a]: S.add(ordinal(F,G))
    for a in range(1,n+1):
        for b in range(a,n+1):
            c=n-a-b
            if c<b: continue
            for F in P[a]:
                for G in P[b]:
                    for K in P[c]:
                        S.add(c3join(F,G,K))
    P[n]=S
    print('order',n,'profiles',len(S),'time %.1f'%(time.time()-t0),flush=True)
allH=set()
for n in P:
    for F in P[n]: allH.add(F[0])
B=max(allH)
odd=[h for h in range(1,1200,2) if h not in allH]
print('max H',B,'; odd values <1200 missing:',odd[:60],len(odd))
import pickle; pickle.dump(P,open('gram_%d.pkl'%NMAX,'wb'))
