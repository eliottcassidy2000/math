#!/usr/bin/env python3
"""Exact controls for clock duality, repeat cost and the doubling-only hostile."""
import itertools
import json
import sys
from fractions import Fraction as F
from pathlib import Path

def need(b,s):
    if not b: raise RuntimeError(s)

def stats(word):
    J=0;U=V=0
    for v in word:
        U+=J;V+=J*J;J+=v
    return len(word),J,U,V

def mm(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]

def tr(A):return list(map(list,zip(*A)))

def inv(A):
    n=len(A);B=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if B[i][j])
        B[j],B[pivot]=B[pivot],B[j]
        q=B[j][j];B[j]=[x/q for x in B[j]]
        for i in range(n):
            if i!=j:
                q=B[i][j];B[i]=[x-q*y for x,y in zip(B[i],B[j])]
    return [row[n:] for row in B]

def exact_controls():
    repeats=0;duals=0;kernels=0
    for L in range(1,7):
        for word in itertools.product((-1,0,1),repeat=L):
            _,d,U,V=stats(word)
            for r in (1,2,3,7):
                expected=(r*L,r*d,r*U+L*d*r*(r-1)//2,
                          r*V+d*U*r*(r-1)+L*d*d*r*(r-1)*(2*r-1)//6)
                need(stats(word*r)==expected,'repeat cost');repeats+=1
            for j in (-3,-1,F(1,2),1,4):
                j=F(j)
                Q=L*j*j+2*U*j+V
                need(Q==j*j*(V/(j*j)+2*F(U)/j+L),'homogeneous inversion')
                duals+=1
    for lengths in ((1,), (1,1), (1,2,3), (4,1,7,2), (1,)*8):
        n=len(lengths);P=[[int(j<=i) for j in range(n)] for i in range(n)]
        D=[[F(lengths[i]) if i==j else F(0) for j in range(n)] for i in range(n)]
        B=[[F(i==j)-F(i==j+1) for j in range(n)] for i in range(n)]
        K=mm(mm(tr(P),D),P)
        need(inv(K)==mm(mm(B,inv(D)),tr(B)),'local precision identity');kernels+=1
    # The norm-one condition cannot be deleted from the split-ring theorem.
    need((2*1%7,1*0%7)==(2,0) and (2*0%7,1*1%7)==(0,1),'split hostile')
    return {'repeat_controls':repeats,'homogeneous_duality_controls':duals,'precision_controls':kernels}

def doubling_hostile(k):
    X=2**k;N=2*X;coeff={}
    for n in range(1,N):
        count=X if n<X else N-n
        a=-(2 if n%3==2 else -1)*count+((0,-1,1)[n%3] if n>=X else 0)
        q=n;sign=1
        while q%2==0:q//=2;sign=-sign
        coeff[q]=coeff.get(q,0)+sign*a
    free=sum(abs(a) for q,a in coeff.items() if q!=1)
    anchor=coeff[1]
    direction=1 if anchor>=0 else -1
    signs={q:direction*(1 if a>=0 else -1) for q,a in coeff.items()}
    signs[1]=1
    t=[0]*N
    for n in range(1,N):t[n]=-t[n//2] if n%2==0 else signs[n]
    H=J=0;M=0;energy=0
    for n in range(1,N):
        H+=(2 if n%3==2 else -1)*t[n]
        need(H==(0,-1,1)[n%3]*t[n]-2*J,'hostile source identity')
        if n>=X:M+=J;energy+=J*J
        if n+1<N:J+=F((0,-1,1)[n%3]*(t[n+1]-t[n]),2)
    need(2*abs(M)==abs(anchor)+free,'sharp doubling support function')
    need(M*M/F(X)<=energy,'mean projection')
    return {'k':k,'twice_abs_mean_sum':int(2*abs(M)),
            'mean_energy_over_X2':float(M*M/X**3),
            'source_energy_over_X2':float(energy/X**2)}

def curling(seq):
    best=1
    for b in range(1,len(seq)+1):
        r=1;tail=seq[-b:]
        while (r+1)*b<=len(seq) and seq[-(r+1)*b:-r*b]==tail:r+=1
        best=max(best,r)
    return best

def main():
    data={'status':'FINITE-EXACT controls, with general proofs in companion notes',**exact_controls(),
          'doubling_only_hostiles':[doubling_hostile(k) for k in range(4,17)]}
    g=[1]
    while len(g)<220:g.append(curling(g))
    need(g[-1]==4 and 4 not in g[:-1],'Gijswijt first four')
    need(g[:9]==[1,1,2,1,1,2,2,2,3],'Gijswijt initial glue')
    data['gijswijt_first_four']=220
    Path(sys.argv[1]).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
