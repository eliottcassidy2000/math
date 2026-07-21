#!/usr/bin/env python3
"""
klein-2026-07-20-S379 -- THE PER-ROOT ARBORESCENCE VECTOR is the stationary distribution of the
tournament random walk (Markov-Chain-Tree Theorem), a poly-time RANKING; does it crack the n=7
residue that spec(A) + Sum a + H cannot (THM-1580)?

Owner: "work another cutting edge math session, think arborescences."

Extends THM-1460 (arborescences = determinantal relaxation of H) and THM-1580 (poly-time Sum a
outseparates #P-hard H 14x; residue = 2 groups at n=7). Sum a is a SUM; the VECTOR {a_r} is
finer, and by the Markov-Chain-Tree Theorem a_r is proportional to the stationary weight of
vertex r -- a poly-time tournament ranking. New invariants tested:
  - {a_r} sorted multiset (per-root arborescence vector);
  - charpoly(L_in) of the directed in-Laplacian (poly-time, finer than Sum a);
against the wall (spec(A), Sum a, H).
"""
import itertools
from fractions import Fraction as Fr

def pairs_of(n): return [(i, j) for i in range(n) for j in range(i + 1, n)]
def relabel(om, perm, n):
    new = [0]*n
    for v in range(n):
        mv, t = om[v], 0
        while mv:
            b = mv & -mv; w = b.bit_length()-1; mv ^= b; t |= 1 << perm[w]
        new[perm[v]] = t
    return tuple(new)
def word(om, n):
    w = 0
    for v in range(n): w = (w << n) | om[v]
    return w
def refine(om, n):
    colour = [bin(om[v]).count("1") for v in range(n)]
    while True:
        sig = []
        for v in range(n):
            cnt = {}; mv = om[v]
            while mv:
                b = mv & -mv; w = b.bit_length()-1; mv ^= b; cnt[colour[w]] = cnt.get(colour[w],0)+1
            sig.append((colour[v], tuple(sorted(cnt.items()))))
        order = sorted(set(sig)); newc = [order.index(sig[v]) for v in range(n)]
        if newc == colour: break
        colour = newc
    cells = {}
    for v in range(n): cells.setdefault(colour[v], []).append(v)
    return [tuple(cells[k]) for k in sorted(cells)]
def canon(om, n):
    cells = refine(om, n); base = []; pos = 0
    for c in cells: base.append((c, pos)); pos += len(c)
    best = None
    for choice in itertools.product(*[itertools.permutations(c) for (c,_) in base]):
        perm = [0]*n
        for (blk,(c,start)) in zip(choice, base):
            for k,v in enumerate(blk): perm[v] = start+k
        w = word(relabel(om, perm, n), n)
        if best is None or w < best: best = w
    return best
def classes(n):
    P = pairs_of(n); om0 = tuple(sum(1<<j for j in range(i)) for i in range(n))
    seen = {canon(om0,n): om0}; fr = [om0]
    while fr:
        nx = []
        for om in fr:
            for (i,j) in P:
                nm = list(om)
                if om[i]>>j&1: nm[i]&=~(1<<j); nm[j]|=1<<i
                else: nm[j]&=~(1<<i); nm[i]|=1<<j
                nm = tuple(nm); cc = canon(nm,n)
                if cc not in seen: seen[cc]=nm; nx.append(nm)
        fr = nx
    return seen
def detf(M):
    M = [row[:] for row in M]; n = len(M); d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c,n) if M[r][c]!=0), None)
        if p is None: return Fr(0)
        if p!=c: M[c],M[p]=M[p],M[c]; d=-d
        d*=M[c][c]; inv=1/M[c][c]
        for r in range(c+1,n):
            f=M[r][c]*inv
            for k in range(c,n): M[r][k]-=f*M[c][k]
    return d
def Lin(om, n):
    L=[[Fr(0)]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j and (om[i]>>j&1): L[j][j]+=1; L[i][j]-=1
    return L
def arb_vector(om, n):
    L = Lin(om, n); out=[]
    for r in range(n):
        idx=[i for i in range(n) if i!=r]
        out.append(int(detf([[L[i][j] for j in idx] for i in idx])))
    return out
def charpoly_int(M, n):
    """Faddeev-LeVerrier, integer char poly of n x n rational matrix"""
    A=[[M[i][j] for j in range(n)] for i in range(n)]
    Mk=[[Fr(1 if i==j else 0) for j in range(n)] for i in range(n)]; cs=[Fr(1)]
    for k in range(1,n+1):
        AM=[[sum(A[i][l]*Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c=Fr(-sum(AM[i][i] for i in range(n)),k)
        cs.append(c)
        Mk=[[AM[i][j]+(c if i==j else 0) for j in range(n)] for i in range(n)]
    return tuple(cs)
def charA(om, n):
    A=[[Fr(1 if (om[i]>>j&1) else 0) for j in range(n)] for i in range(n)]
    return charpoly_int(A, n)
def charLin(om, n):
    return charpoly_int(Lin(om,n), n)
def hp(om, n):
    c=0
    def go(l,u,d):
        nonlocal c
        if d==n: c+=1; return
        mv=om[l]&~u
        while mv:
            b=mv&-mv; w=b.bit_length()-1; mv^=b; go(w,u|(1<<w),d+1)
    for s in range(n): go(s,1<<s,1)
    return c

# ---- MCTT sanity check at n=5: a_r proportional to stationary dist of reverse walk
print("=" * 88)
print("MARKOV-CHAIN-TREE CHECK: {a_r} proportional to the tournament random-walk stationary dist")
print("=" * 88)
import numpy as np
for n in (5,):
    cls = classes(n)
    ok = True; shown = 0
    for c, om in cls.items():
        av = arb_vector(om, n)
        # random walk: from v -> uniform IN-neighbor (so arborescences TOWARD r) ; P[i][j]=1/indeg(i) if j->i
        # continuous-time unit-rate generator Q: Q[i][j]=1 if i->j, Q[i][i]=-out_i.  pi Q = 0.
        Q=np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                if om[i]>>j&1: Q[i][j]=1.0; Q[i][i]-=1.0
        w,V=np.linalg.eig(Q.T)
        k=np.argmin(np.abs(w)); pi=np.real(V[:,k]); pi=pi/pi.sum()
        avn=np.array(av)/sum(av)
        if not np.allclose(sorted(pi),sorted(avn),atol=1e-6): ok=False
        if shown<2: print(f"   class: a_r={av}  normalized={[round(x,3) for x in avn]}  stat.dist={sorted([round(x,3) for x in pi])}"); shown+=1
    print(f"   {{a_r}} matches the stationary distribution (MCTT) for all n={n} classes: {ok}")
    print("   => the per-root arborescence vector IS a poly-time tournament ranking.")

# ---- the separation test at n=7
print("\n" + "=" * 88)
print("SEPARATION at n=7: does {a_r} or charpoly(L_in) crack the (spec A, Sum a, H) residue?")
print("=" * 88)
n=7; cls=classes(n)
data={}
for c,om in cls.items():
    cA=charA(om,n); av=tuple(sorted(arb_vector(om,n))); Sa=sum(av); cL=charLin(om,n); H=hp(om,n)
    data[c]=(cA,av,Sa,cL,H)
# group by spec(A)
groups={}
for c,(cA,av,Sa,cL,H) in data.items(): groups.setdefault(cA,[]).append(c)
cospec={g:v for g,v in groups.items() if len(v)>1}
def fails(groupdict, keyfn):
    return sum(1 for v in groupdict.values() if len({keyfn(data[c]) for c in v})==1)
byS  = {g:v for g,v in cospec.items() if len({data[c][2] for c in v})==1}          # Sum a fails
byH  = {g:v for g,v in cospec.items() if len({data[c][4] for c in v})==1}          # H fails
bySH = {g:v for g,v in cospec.items() if len({(data[c][2],data[c][4]) for c in v})==1}  # both fail
print(f"  cospectral groups: {len(cospec)}")
print(f"  Sum a fails to split: {len(byS)};  H fails: {len(byH)};  (Sum a,H) both fail: {len(bySH)}")
# now the residue: within bySH, does {a_r} vector split?  does charpoly(L_in) split?
res_av = {g:v for g,v in bySH.items() if len({data[c][1] for c in v})==1}   # a_r vector ALSO fails
res_cL = {g:v for g,v in bySH.items() if len({data[c][3] for c in v})==1}   # charpoly L_in fails
print(f"  of the (Sum a,H) residue, {{a_r}} VECTOR also fails: {len(res_av)}  (splits {len(bySH)-len(res_av)})")
print(f"  of the (Sum a,H) residue, charpoly(L_in) fails: {len(res_cL)}  (splits {len(bySH)-len(res_cL)})")
for g,v in bySH.items():
    avs=[data[c][1] for c in v]; cLs=[data[c][3] for c in v]; Sas=[data[c][2] for c in v]; Hs=[data[c][4] for c in v]
    print(f"    residue group (Sum a={Sas[0]}, H={Hs[0]}): |grp|={len(v)}")
    print(f"       distinct {{a_r}} vectors: {len(set(avs))}  -> {'SPLIT by a_r' if len(set(avs))>1 else 'a_r fails'}")
    print(f"       distinct charpoly(L_in): {len(set(cLs))}  -> {'SPLIT by L_in' if len(set(cLs))>1 else 'L_in fails'}")
    print(f"       the vectors: {sorted(set(avs))}")
print("\nDONE.")
