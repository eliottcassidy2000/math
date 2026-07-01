#!/usr/bin/env python3
"""
tournament_dihedral_atoms_klein.py  --  klein-2026-07-01-S70

RUNNERS-ON-THE-LOOP  x  TOURNAMENTS  x  DIHEDRAL, on the 7 odd points of (1/14)Z (= roots of -1).

Setup (owner's construction). The tight AP core {1..13} has lonely measure = 6 ATOMS at (Z/14)* =
{1,3,5,9,11,13}/14 (units, 6=phi(14)). ADD the point 7/14 = 1/2 (the iota-fixed point, 7 = odd non-unit).
The 7 points {1,3,5,7,9,11,13}/14 are the SEVENTH ROOTS OF -1 (z^7=-1), carrying a DIHEDRAL D_7 action:
rotation by 2/14=1/7 and reflection = antipode iota: t->1-t (the runner COMPLEMENT / tournament CONVERSE).
Label them by k=0..6 via point (2k+1)/14; differences live in Z/7 (the prime 7 = n/2 = the 7-VANISHING,
THM-503). Build tournaments on these 7 vertices (Paley/QR = the dihedral one) and apply the TILING model.

Computes: (A) verify the 6 atoms + their Verblunsky (|alpha_5|=1 = 6-atom termination); (B) the QR/Paley
and rotational tournaments on Z/7: scores, odd-cycle counts (OCF N), Redei Hamiltonian-path count H(T),
|Aut|; (C) the TILING (n=7) -- circulant => difference-striped => grid-symmetric (BLUE); (D) the invariant
dictionary runner<->measure<->tournament.
"""
import numpy as np
from itertools import combinations, permutations

# ---------- (A) the 6 atoms + Verblunsky ----------
def norms(v, t):
    x = (v*t) % 1.0
    return np.minimum(x, 1-x)

def verblunsky(c):
    K=len(c)-1; a=[1.0]; E=c[0]; al=[]
    for n in range(1,K+1):
        if E<=1e-13: break
        acc=c[n]+sum(a[j]*c[n-j] for j in range(1,n)); k=-acc/E
        k=max(-1.0,min(1.0,k)); a=[1.0]+[a[j]+k*a[n-j] for j in range(1,n)]+[k]; E*=1-k*k; al.append(k)
    return np.array(al)

def part_A(n=14):
    print("(A) THE 6 UNIT ATOMS of the tight AP core {1..13} and their Verblunsky sequence")
    core=list(range(1,n)); N=200004; t=np.arange(N)/N
    G=np.full(N,1.0)
    for v in core: G=np.minimum(G,norms(v,t))
    r=1.0/n; peak=G.max()
    # atoms = local maxima reaching ~1/14
    idx=np.where(G>peak-1e-4)[0]; atoms=sorted(set(round(14*t[i]) for i in idx))
    units=[u for u in range(1,n) if np.gcd(u,n)==1]
    print(f"    max min_v||vt|| = {peak:.5f} (=1/14={r:.5f}? tight); atom positions k/14 for k in {atoms}")
    print(f"    (Z/14)* units = {units}  (6=phi(14)); match => atoms ARE the units")
    # Verblunsky of the 6-atom measure (equal weights at units/14)
    K=12; c=np.array([np.mean([np.cos(2*np.pi*k*u/n) for u in units]) for k in range(K+1)])
    al=verblunsky(c)
    print(f"    6-atom measure moments c_1={c[1]:+.4f}; Verblunsky |alpha_n|: "+" ".join(f"{abs(x):.3f}" for x in al))
    print(f"    LADDER: |alpha_k| = 1/(6-k)? "+" ".join(f"{abs(al[k]):.3f}~{1/(6-k):.3f}" for k in range(len(al))))
    print(f"    terminates at n={len(al)} with |alpha_{{{len(al)-1}}}|={abs(al[-1]):.4f} (=1 => EXACTLY {len(al)} atoms = 6 units) [Verblunsky-Geronimus]")
    # moments ARE Ramanujan sums c_n(k)/phi(n)  (=> convergence with mac-mini/kps HYP-3793)
    def ramanujan(k,nn):
        g=np.gcd(k,nn); m=nn//g
        # c_n(k) = mu(m)*phi(n)/phi(m)
        def mu(x):
            x=int(x); r=1; d=2
            while d*d<=x:
                if x%d==0:
                    x//=d
                    if x%d==0: return 0
                    r=-r
                d+=1
            if x>1: r=-r
            return r
        def phi(x):
            x=int(x); res=x; d=2
            while d*d<=x:
                if x%d==0:
                    while x%d==0: x//=d
                    res-=res//d
                d+=1
            if x>1: res-=res//x
            return res
        return mu(m)*phi(nn)//phi(m)
    rk=[ramanujan(k,n)/len(units) for k in range(7)]
    print(f"    MOMENTS = RAMANUJAN SUMS c_14(k)/phi(14): c_k measured "+" ".join(f"{c[k]:+.3f}" for k in range(7)))
    print(f"                                              c_14(k)/6     "+" ".join(f"{v:+.3f}" for v in rk)+"  (=> HYP-3793 convergence)")

# ---------- (B) tournaments on Z/7 ----------
def circ_adj(S, p=7):
    return [[1 if (i-j)%p in S else 0 for j in range(p)] for i in range(p)]

def scores(A): return [sum(row) for row in A]
def count_cycles(A, L):
    # count DIRECTED L-cycles, each once (fix start = min vertex; direction set by the arcs)
    p=len(A); cnt=0
    for c in permutations(range(p), L):
        if c[0]!=min(c): continue                      # each directed cycle counted once
        if all(A[c[i]][c[(i+1)%L]] for i in range(L)): cnt+=1
    return cnt
def ham_paths(A):
    p=len(A); cnt=0
    for perm in permutations(range(p)):
        if all(A[perm[i]][perm[i+1]] for i in range(p-1)): cnt+=1
    return cnt
def aut_order(S, p=7):
    # affine maps x->m x + b, m in (Z/p)*, preserving connection set S (m*S=S)
    order=0
    for m in range(1,p):
        if set((m*s)%p for s in S)==set(S): order+=p   # each valid m gives p translations
    return order

def part_B():
    print("\n(B) TOURNAMENTS on Z/7 (differences = the prime 7 = n/2 = 7-vanishing locus)")
    tours={"Paley/QR {1,2,4}":{1,2,4}, "Rotational {1,2,3}":{1,2,3}}
    for name,S in tours.items():
        A=circ_adj(S); sc=scores(A)
        c3=count_cycles(A,3); c5=count_cycles(A,5); c7=count_cycles(A,7)
        N=c3+c5+c7; H=ham_paths(A); au=aut_order(S)
        conv = set((-s)%7 for s in S)
        selfconv = any(set((m*s)%7 for s in S)==conv for m in range(1,7))
        print(f"   {name}: scores={sc} (regular={len(set(sc))==1}); c3={c3} c5={c5} c7={c7}  N(OCF)={N}; H(Redei Ham-paths)={H} (odd={H%2==1}); |Aut|={au}; self-converse={selfconv}")

# ---------- (C) the tiling (n=7): circulant => difference-striped => grid-symmetric ----------
def part_C():
    print("\n(C) TILING model (n=7, staircase delta_5): circulant tournament => DIFFERENCE-STRIPED tiling")
    for name,S in {"Paley/QR":{1,2,4}, "Rotational":{1,2,3}}.items():
        # base path 6->5->...->0 (needs 1 in S). tile (x,y), x>y, x-y>=2 (0-indexed 0..6). set iff (x-y) in S
        tiles=[(x,y) for y in range(0,5) for x in range(6,y+1,-1) if x-y>=2]
        bits={(x,y):(1 if (x-y)%7 in S else 0) for (x,y) in tiles}
        set_diffs=sorted(set(x-y for (x,y),b in bits.items() if b))
        # grid symmetry: transform (x,y)->(6-y,6-x) preserves x-y => circulant tiling invariant
        gridsym=all(bits[(x,y)]==bits.get((6-y,6-x), bits[(x,y)]) for (x,y) in tiles if (6-y,6-x) in bits)
        nset=sum(bits.values())
        print(f"   {name}: 1 in S={1 in S} (base path OK); tiles set at differences {set_diffs} ({nset}/{len(tiles)} tiles); grid-symmetric(BLUE)={gridsym}")
    print("   => circulant/dihedral tournaments are difference-striped; grid-sym preserves x-y => all BLUE (self-converse).")

def part_D():
    print("\n(D) INVARIANT DICTIONARY  runner/measure  <->  Verblunsky/OPUC  <->  tournament/tiling")
    rows=[
     ("6 atoms = (Z/14)* units","|alpha_5|=1 (6-atom termination)","7 vertices (6 units + iota-fixed 7/14)"),
     ("antipode iota: t->1-t","conjugation z->1/z-bar (reflection)","tournament CONVERSE T^op = x->-x"),
     ("7-vanishing (prime n/2=7)","moment period tied to 7","tournament lives on Z/7 (QR/Paley)"),
     ("dihedral D_7 on roots of -1","rotational Verblunsky phases","circulant => Aut>=Z/7; Paley Aut=21"),
     ("runner speed M_v: t->vt","degree-v Blaschke product","arc/connection multiplier m in (Z/7)*"),
     ("binding depth M","alpha_0=-cos, termination level","OCF N / Ham-path count H (parity, Redei)"),
    ]
    w=max(len(a) for a,_,_ in rows)
    for a,b,c in rows: print(f"   {a:<{w}} | {b:<34} | {c}")
    print("   KEY BRIDGE: prime 2 = antipode/converse (iota, self-converse tournaments); prime 7 = n/2 = the")
    print("   7-vanishing = the vertex set Z/7 of the tournament; QR(7) = the Paley orientation = the dihedral one.")

if __name__=="__main__":
    part_A(); part_B(); part_C(); part_D()
