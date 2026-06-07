#!/usr/bin/env python3
"""
S642 / HYP-2320 — Commutator depth = solvability by radicals; the cube-root 3-cycle engine.

The user's picture: permuting the n roots acts on the n+1 coefficients; solvability by
radicals <=> the Galois group's DERIVED SERIES terminates <=> nested commutators die.
Ladder:  S2 (1 swap), S3 (1 commutator: der^2=1), S4 (2: der^3=1), S5 (never: A5 perfect).

This script:
  (A) computes the DERIVED LENGTH of S_n for n=2..6 (the exact commutator depth);
  (B) verifies A5 is PERFECT ([A5,A5]=A5) -- commutators never shrink -- and that the
      engine is two 3-cycles sharing ONE point (the cube-root atoms, S635);
  (C) the arc connections: 3-cycle = cube root of unity; solvable = tower of cyclic
      (root-of-unity) layers; and the constructibility corner (the 7-gon, 7=Phi_3(2)).
Pure-python permutation groups (no external libs).
"""
from itertools import permutations, product

# ---------- permutation helpers (tuples = images of 0..n-1) ----------
def comp(p, q):                      # (p*q)(i) = p(q(i))
    return tuple(p[q[i]] for i in range(len(p)))
def inv(p):
    r=[0]*len(p)
    for i,x in enumerate(p): r[x]=i
    return tuple(r)
def comm(a,b):                       # [a,b] = a b a^-1 b^-1
    return comp(comp(a,b), comp(inv(a),inv(b)))
def Sn(n):
    return [tuple(p) for p in permutations(range(n))]
def gen_subgroup(gens, n):
    ident=tuple(range(n)); G={ident}; frontier=[ident]
    while frontier:
        x=frontier.pop()
        for g in gens:
            y=comp(x,g)
            if y not in G: G.add(y); frontier.append(y)
    return G
def derived(G, n):                   # [G,G]
    gens={comm(a,b) for a in G for b in G}
    return gen_subgroup(gens, n)
def derived_length(G, n):
    H=set(G); k=0
    while len(H)>1:
        H2=derived(H,n)
        if H2==H:                    # perfect: never terminates
            return None
        H=H2; k+=1
    return k

print("="*66)
print("(A) derived length of S_n (= commutator depth = solvability)")
print("="*66)
print("  n | |S_n| | derived length | solvable? (radical formula?)")
for n in range(2,7):
    G=set(Sn(n)); dl=derived_length(G,n)
    sols = "YES (formula)" if dl is not None else "NO (A_%d perfect)"%n
    dlstr = str(dl) if dl is not None else "infinity"
    print(f"  {n} | {len(G):4d}  | {dlstr:>14} | {sols}")
print("  -> der^2(S3)=1 (cubic: 1 commutator), der^3(S4)=1 (quartic: 2),")
print("     S5 NEVER terminates (quintic: no formula) -- Abel-Ruffini.")

print("\n" + "="*66)
print("(B) the engine: A5 perfect, built from 3-cycles sharing ONE point")
print("="*66)
def threecycle(a,b,c,n):
    p=list(range(n)); p[a]=b; p[b]=c; p[c]=a; return tuple(p)
A5=gen_subgroup([threecycle(0,1,2,5), threecycle(2,3,4,5)],5)
print(f"  <(012),(234)> has order {len(A5)}  (= |A5| = 60? {len(A5)==60})")
dA5=derived(A5,5)
print(f"  [A5,A5] order = {len(dA5)}; A5 perfect ([A5,A5]=A5)? {dA5==A5}")
# nested commutators of the two overlapping 3-cycles never die:
s=threecycle(0,1,2,5); t=threecycle(2,3,4,5); ident=tuple(range(5))
chain=[("[s,t]",comm(s,t))]
x=comm(s,t)
for d in range(2,7):
    x = comm(x, s if d%2==0 else t)
    chain.append((f"depth-{d}", x))
for name,g in chain:
    print(f"  {name:10s}: {'NONTRIVIAL' if g!=ident else 'trivial'}  {g}")
print("  -> two 3-cycles overlapping in one point: commutators never vanish.")
# contrast: A3 (cubic) is abelian -> its commutators die at once
A3=gen_subgroup([threecycle(0,1,2,3)],3)
print(f"  A3 = <(012)> order {len(A3)}; [A3,A3] order {len(derived(A3,3))} (=1 => cube-root group abelian)")

print("\n" + "="*66)
print("(C) arc connections")
print("="*66)
print("  * 3-cycle = CUBE ROOT of unity (s^3=1, eigenvalues 1,w,w^2); the atom of A_n (S635).")
print("  * solvable by radicals = tower of CYCLIC (root-of-unity) layers; the derived")
print("    series = peeling abelian (cyclotomic) shells -- the same shell-tower as the arc.")
print("  * monodromy (opus S699p): a physical root-swap = a transposition picked up looping")
print("    the coefficients around the DISCRIMINANT (= the worry-set / branch locus); the")
print("    Galois group = the monodromy group of the FTA covering. Solvable = the monodromy")
print("    is built from abelian loops.")
print("  * CONSTRUCTIBILITY corner (Gauss-Wantzel): the regular 7-gon is NOT constructible")
print("    because [Q(zeta_7):Q] = phi(7) = 6 = 2*3 is divisible by 3 (a CUBE root), not a")
print("    power of 2. The same 7 = Phi_3(2) = N(3+w) the arc converges on (S637/8/40): the")
print("    cube-root obstruction to the heptagon. (7-gon IS solvable by radicals -- Galois")
print("    group Z/6 cyclic -- but needs a cube root, so ruler+compass cannot reach it.)")
from math import gcd
def phi(n): return sum(1 for k in range(1,n) if gcd(k,n)==1)
for m in [3,5,7,9,15,17,257]:
    pe = phi(m); pow2 = (pe & (pe-1))==0
    print(f"    {m}-gon: [Q(zeta_{m}):Q]=phi({m})={pe}  power-of-2? {pow2}  "
          f"constructible? {pow2}")
print("  3-gon,5-gon constructible (phi=2,4); 7-gon,9-gon NOT (phi=6 has a 3).")
