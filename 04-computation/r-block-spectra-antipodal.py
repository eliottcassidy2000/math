#!/usr/bin/env python3
"""
r-block-spectra-antipodal.py   (klein-2026-06-29-S1)

HYP-B confirmed complement (R) commutes with the iso-class arc-flip adjacency A.
This refines WHY: the labeled arc-flip graph on 2^C(n,2) tournaments is the
HYPERCUBE Q_{C(n,2)} (one arc flip = one cube edge). Complement reverses every
arc = flips every bit = the ANTIPODAL MAP of the cube (the Borsuk-Ulam map).
A is the S_n-quotient of Q_d (d=C(n,2)); its eigenvalues are level eigenvalues
d-2k. The antipodal map acts as (-1)^k on level k, so:

  R-even (antipodal +1) block  <-> EVEN levels k -> eigenvalues  =  d (mod 4)
  R-odd  (antipodal -1) block  <-> ODD  levels k -> eigenvalues  =  d-2 (mod 4)

Test: build A, split into R-even / R-odd blocks (sym/antisym basis of the
complement involution), compute each block's spectrum, and check the mod-4
residue law. Confirms: complement = antipodal, R-split = cube level parity.
"""
import itertools
from math import comb
import numpy as np

def pairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def perm_tables(n):
    P=pairs(n); idx={p:k for k,p in enumerate(P)}; T=[]
    for perm in itertools.permutations(range(n)):
        row=[]
        for (i,j) in P:
            a,b=perm[i],perm[j]
            row.append((idx[(a,b)],False) if a<b else (idx[(b,a)],True))
        T.append(row)
    return T
def canonical(bits,T):
    best=None
    for row in T:
        v=0
        for q,(tgt,inv) in enumerate(row):
            bit=(bits>>tgt)&1
            if inv: bit^=1
            v|=bit<<q
        if best is None or v<best: best=v
    return best

def analyze(n):
    d=comb(n,2); T=perm_tables(n)
    # enumerate classes over ALL 2^d labeled tournaments (no fixed base path here:
    # we want the full arc-flip cube quotient). For n<=6 this is fine.
    class_of={}; classes=[]
    for bits in range(2**d):
        c=canonical(bits,T)
        class_of[bits]=c
        if c==bits: classes.append(c)   # canonical reps are their own min
    cidx={c:i for i,c in enumerate(classes)}; N=len(classes)
    # adjacency A[i][j] = #arcs whose flip takes rep_i to class j
    A=np.zeros((N,N))
    for c in classes:
        for a in range(d):
            A[cidx[c]][cidx[class_of[c^a if False else c^(1<<a)]]]+=1
    # complement permutation sigma (antipodal: flip all d bits)
    full=(1<<d)-1
    sigma={c:class_of[c^full] for c in classes}
    SC=sum(1 for c in classes if sigma[c]==c)
    # build sym/antisym basis of the involution P_sigma
    Eb=[]; Ob=[]
    done=set()
    for c in classes:
        if c in done: continue
        s=sigma[c]; i=cidx[c]; j=cidx[s]
        if s==c:
            v=np.zeros(N); v[i]=1.0; Eb.append(v); done.add(c)
        else:
            ve=np.zeros(N); ve[i]=ve[j]=1/np.sqrt(2); Eb.append(ve)
            vo=np.zeros(N); vo[i]=1/np.sqrt(2); vo[j]=-1/np.sqrt(2); Ob.append(vo)
            done.add(c); done.add(s)
    E=np.array(Eb).T; O=np.array(Ob).T if Ob else np.zeros((N,0))
    # A restricted to each block (A commutes with P_sigma so blocks are invariant)
    Aee=E.T@A@E
    Aoo=O.T@A@O if O.shape[1] else np.zeros((0,0))
    even_ev=sorted(round(float(x.real),2) for x in np.linalg.eigvals(Aee))
    odd_ev =sorted(round(float(x.real),2) for x in np.linalg.eigvals(Aoo)) if Aoo.size else []
    print(f"\n n={n}: d=C(n,2)={d},  classes N={N} (A000568), SC={SC}, "
          f"R-even dim={E.shape[1]}, R-odd dim={O.shape[1]}")
    print(f"   R-even block spectrum: {even_ev}")
    print(f"   R-odd  block spectrum: {odd_ev}")
    # mod-4 law: even-level eigenvalues = d mod 4 ; odd-level = d-2 mod 4
    em = set(int(round(x))%4 for x in even_ev)
    om = set(int(round(x))%4 for x in odd_ev)
    print(f"   predicted residues: R-even = d%4 = {d%4} ; R-odd = (d-2)%4 = {(d-2)%4}")
    print(f"   observed  residues: R-even = {em} ; R-odd = {om}   "
          f"=> antipodal/level-parity law {'HOLDS' if em<= {d%4} and (om<= {(d-2)%4} or not om) else 'FAILS'}")
    # where is the Perron and where is the most negative?
    allev=even_ev+odd_ev
    print(f"   Perron {d} in R-even: {d in [int(round(x)) for x in even_ev]} ; "
          f"global min eig = {min(allev)}  in {'R-even' if min(allev) in even_ev else 'R-odd'}")

if __name__=="__main__":
    for n in (3,4,5,6):
        analyze(n)
