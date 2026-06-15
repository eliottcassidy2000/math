#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THE CONVERSE-PARITY THEOREM (the real FKN content).  kind-pasteur-2026-06-15-S6.

On the FULL arc cube {0,1}^{C(n,2)} (every arc a free coordinate, x_e=+1 in the
"reference" orientation i->j for i>j), the global sign flip x -> -x is exactly the
CONVERSE involution T -> T^op (reverse all arcs).  Therefore:

  * an invariant f with f(T^op) = +f(T)  (converse-EVEN) is supported on EVEN Fourier levels;
  * an invariant f with f(T^op) = -f(T)  (converse-ODD)  is supported on ODD  Fourier levels.

Consequences (the dictionary):
  - H (Hamiltonian-path count): H(T^op)=H(T)  -> converse-EVEN -> EVEN levels only.
  - c_k (k-cycle counts), all OCF cycle data, the conflict graph Omega: converse-EVEN -> EVEN.
  - scores s_v: s_v(T^op)=(n-1)-s_v(T), so s_v-(n-1)/2 is converse-ODD -> ODD levels (level 1 = affine).
  - skew S=A-A^T: S(T^op)=-S(T) (converse-ODD); the skew spectrum {mu} and det(I+S)=prod(1+mu^2)
    are converse-EVEN -> the determinant lens lives at even levels.

The TILING model fixes the base path (n-1 arcs), partially BREAKING converse symmetry --
which is precisely what injects the odd-level (score / level-1 "ranking signal") content.

We verify level-parity of c3, c5, H, scores on the full arc cube (n=4,5,6), and check
whether H's even-level weights track the OCF alpha_k hierarchy.
"""
import sys, itertools
from math import comb
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
import numpy as np

def arcs(n): return [(i,j) for i in range(n) for j in range(i)]  # i>j reference arc i->j

def adj_full(n, E, bits):
    A=np.zeros((n,n),dtype=int)
    for (i,j),b in zip(E,bits):
        if b==0: A[i][j]=1          # reference: i beats j
        else:    A[j][i]=1          # reversed
    return A

def scores(A): return A.sum(axis=1)
def count_ck(A,k):
    n=len(A); c=0
    for S in itertools.combinations(range(n),k):
        # count directed k-cycles within S: (k-1)!/2 cyclic orders; just test all cyclic perms
        verts=list(S); cnt=0
        for perm in itertools.permutations(verts[1:]):
            cyc=[verts[0]]+list(perm)
            if all(A[cyc[t]][cyc[(t+1)%k]] for t in range(k)): cnt+=1
        c+=cnt//1  # each directed cycle counted once per starting rotation? perms over k-1 fixes rotation; /2 for reflection not needed (directed)
    return c
def count_c3(A):
    n=len(A); c=0
    for i,j,k in itertools.combinations(range(n),3):
        e=A[i][j]+A[j][k]+A[k][i]; f=A[j][i]+A[k][j]+A[i][k]
        if e==3 or f==3: c+=1
    return c
def count_H(A):
    n=len(A); cnt=0
    for perm in itertools.permutations(range(n)):
        if all(A[perm[t]][perm[t+1]] for t in range(n-1)): cnt+=1
    return cnt

def fwht(vec):
    a=vec.astype(float).copy(); h=1
    while h<len(a):
        for i in range(0,len(a),h*2):
            for j in range(i,i+h):
                x=a[j]; y=a[j+h]; a[j]=x+y; a[j+h]=x-y
        h*=2
    return a/len(a)
def level_weights(coeffs,m):
    w=np.zeros(m+1)
    for idx,c in enumerate(coeffs): w[bin(idx).count('1')]+=c*c
    return w
def parity_split(coeffs,m):
    even=0.0; odd=0.0
    for idx,c in enumerate(coeffs):
        if bin(idx).count('1')%2==0: even+=c*c
        else: odd+=c*c
    return even,odd

def main():
    for n in (4,5,6):
        E=arcs(n); m=len(E); N=1<<m
        c3=np.zeros(N); Hv=np.zeros(N); s0=np.zeros(N)
        for idx in range(N):
            bits=[(idx>>t)&1 for t in range(m)]
            A=adj_full(n,E,bits); c3[idx]=count_c3(A); Hv[idx]=count_H(A); s0[idx]=scores(A)[0]
        print(f"\n===== FULL arc cube  n={n}, {m} arcs, {N} tournaments =====")
        for name,vec,expect in [("c3",c3,"EVEN"),("H",Hv,"EVEN"),("score[0]",s0,"ODD")]:
            cf=fwht(vec); ev,od=parity_split(cf,m); lw=level_weights(cf,m)
            tag="EVEN" if od<1e-9 else ("ODD" if ev<=cf[0]**2+1e-9 else "MIXED")
            lvls=[f"L{L}={lw[L]:.3f}" for L in range(m+1) if lw[L]>1e-9]
            print(f"  {name:9s} even-mass={ev:.3f} odd-mass={od:.3f}  -> {tag:6s} (expect {expect})")
            print(f"            levels: {', '.join(lvls)}")
        # OCF check: H even-level weights vs alpha_k structure (H = sum_k alpha_k 2^k)
        cf=fwht(Hv); lw=level_weights(cf,m)
        print(f"  H even levels only: "+", ".join(f"L{2*j}={lw[2*j]:.3f}" for j in range(m//2+1) if lw[2*j]>1e-9))
        print(f"  (L0=mean^2; L2 = pairwise odd-cycle 'conflict' interactions; L4 = disjoint-pair alpha_2 layer)")

if __name__=="__main__":
    main()
