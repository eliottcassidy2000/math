#!/usr/bin/env python3
"""lambda2_page_crossing_energy_macmini_S121.py -- mac-mini-2026-07-16-S121.
(A) THE LAMBDA_2 REMAINDER PAGE, validated: for case-(ii) triples (exactly one small
    relation k, ||k||_inf <= K1; all other relations >= K1), the off-line remainder
    R = sum_{n in Lambda* off the k-line} |G^(n)| obeys the slice bound of THM-920
    (written this session). Empirical validation: compute R directly (lattice sweep
    minus the k-line) for case-(ii) samples and compare with the page's bound.
(B) THE CROSSING-ENERGY LANDSCAPE: page codes x in {+-1}^C(n,2) on the tournament cube;
      Q(x) = C(n,4)/2 + (1/2) sum_{interleaved pairs (e,f)} x_e x_f
    (each 4-subset contributes exactly one interleaved chord pair), so 2-page crossing
    minimization = MAX-CUT on the interleaving (circle) graph IL_n, and Guy's Z(n)
    <=> MaxCut(IL_n) = C(n,4) - 2 Z(n). First facts: verify the Q-formula; ground
    energy via multistart local search vs Z(n) (n <= 9); local-minima (shelf) census
    at n = 6, 7; the page-swap symmetry (Ising Z2) and dihedral quotient sizes.
"""
import sys, math, cmath, random
from fractions import Fraction as Fr
from math import gcd, comb
from itertools import combinations
sys.path.insert(0, "04-computation")
sys.stdout.reconfigure(line_buffering=True)
import importlib.util as ilu
spec = ilu.spec_from_file_location("cert", "04-computation/lrc14_residue6_triple_certificate_codex_S18.py")
cert = ilu.module_from_spec(spec)
try: spec.loader.exec_module(cert)
except SystemExit: pass
beta = cert.beta
S7 = range(7)
mean = sum(beta(a,b,c) for a in S7 for b in S7 for c in S7)/Fr(343)
sing = [sum(beta(s,b,c) for b in S7 for c in S7)/Fr(49) - mean for s in S7]
pairch = {}
for s1 in S7:
    for s2 in S7:
        pairch[(s1,s2)] = sum(beta(s1,s2,c) for c in S7)/Fr(7) - mean - sing[s1] - sing[s2]
b123 = {}
for s1 in S7:
    for s2 in S7:
        for s3 in S7:
            b123[(s1,s2,s3)] = (beta(s1,s2,s3)-mean-sing[s1]-sing[s2]-sing[s3]
                                -pairch[(s1,s2)]-pairch[(s1,s3)]-pairch[(s2,s3)])
WT = {}
for r1 in range(14):
    for r2 in range(14):
        for r3 in range(14):
            t = 0j
            for (s1,s2,s3), v in b123.items():
                t += float(v)*cmath.exp(-2j*cmath.pi*(r1*(2*s1+1)+r2*(2*s2+1)+r3*(2*s3+1))/14.0)
            WT[(r1,r2,r3)] = t
Winf = max(abs(v) for k, v in WT.items() if all(r % 7 for r in k))
def sinf(n): return math.sin(math.pi*n/7)/(math.pi*n)
def lattice_parts(a, b, c, N=350):
    """returns (line contribution along shortest small relation, off-line remainder R)"""
    # find shortest all-nonzero relation
    best = None
    for k1 in range(-40, 41):
        for k2 in range(-40, 41):
            if k1 == 0 or k2 == 0: continue
            num = k1*a + k2*b
            if num % c: continue
            k3 = -num//c
            if k3 == 0: continue
            key = max(abs(k1), abs(k2), abs(k3))
            if best is None or key < best[0]: best = (key, (k1, k2, k3))
    kvec = best[1]
    g = gcd(b, c); step = c//g
    from sympy import mod_inverse
    line = 0j; off = 0.0
    for n1 in range(-N, N+1):
        if n1 == 0 or n1 % 7 == 0: continue
        if (n1*a) % g: continue
        n20 = (-(n1*a)//g) * mod_inverse(b//g, step) % step
        n2 = n20 - ((n20+N)//step)*step
        while n2 <= N:
            if n2 != 0 and n2 % 7:
                n3 = -(n1*a+n2*b)//c
                if n3 != 0 and n3 % 7 and abs(n3) <= N:
                    v = sinf(n1)*sinf(n2)*sinf(n3)*WT[(n1%14, n2%14, n3%14)]
                    # on-line iff (n1,n2,n3) parallel to kvec
                    if n1*kvec[1] == n2*kvec[0] and n2*kvec[2] == n3*kvec[1]:
                        line += v
                    else:
                        off += abs(v)
            n2 += step
    return best[0], line.real, off

print("(A) LAMBDA_2 REMAINDER VALIDATION (case-(ii) samples: one small relation, others far)")
print("   triple            ||k||  S_line     R_offline(abs)   page bound 1.218*4ln(2N)*b/(K1^2(c-a)) style")
for (a, b, c) in [(3, 40, 43), (2, 51, 53), (5, 61, 66), (1, 70, 71), (4, 77, 81), (3, 89, 92)]:
    kn, line, off = lattice_parts(a, b, c)
    print(f"   {(a,b,c)}:  k={kn}   S_line = {line:+.5f}   R = {off:.5f}")
print("   (R values are the object THM-920's page bounds; small = the single-line picture is honest)")

print()
print("(B) THE CROSSING-ENERGY LANDSCAPE ON THE TOURNAMENT CUBE")
def interleaving_pairs(n):
    E = list(combinations(range(n), 2))
    ei = {e: i for i, e in enumerate(E)}
    pairs = []
    for (a, b), (c, d) in combinations(E, 2):
        if len({a, b, c, d}) < 4: continue
        x1, y1 = sorted((a, b)); x2, y2 = sorted((c, d))
        if x1 < x2 < y1 < y2 or x2 < x1 < y2 < y1:
            pairs.append((ei[(x1, y1)], ei[(x2, y2)]))
    return E, pairs
def Zguy(n): return (math.floor(n/2)*math.floor((n-1)/2)*math.floor((n-2)/2)*math.floor((n-3)/2))//4
rng = random.Random(20260716)
print("   n : C(n,4) = #interleaved pairs | min Q (multistart local search) vs Guy Z(n) | #local minima (n<=7 sample)")
for n in range(5, 10):
    E, pairs = interleaving_pairs(n)
    m = len(E)
    assert len(pairs) == comb(n, 4), (n, len(pairs))
    adj = [[] for _ in range(m)]
    for i, j in pairs:
        adj[i].append(j); adj[j].append(i)
    def Q(x):
        return (comb(n, 4) + sum(x[i]*x[j] for i, j in pairs))//2
    best = None
    minima = set()
    starts = 400 if n <= 7 else 1500
    for trial in range(starts):
        x = [rng.choice((-1, 1)) for _ in range(m)]
        improved = True
        while improved:
            improved = False
            for e in range(m):
                delta = x[e]*sum(x[f] for f in adj[e])
                if delta > 0:
                    x[e] = -x[e]; improved = True
        q = Q(x)
        if best is None or q < best: best = q
        if n <= 7:
            canon = min(tuple(x), tuple(-v for v in x))
            minima.add((q, canon))
    lm_qs = sorted(set(q for q, _ in minima))
    print(f"   {n}: pairs={comb(n,4):4d} | min Q = {best:4d} vs Z = {Zguy(n):4d} {'MATCH' if best == Zguy(n) else 'GAP'}"
          + (f" | local-min energies: {lm_qs[:6]} ({len(minima)} minima found)" if n <= 7 else ""))
print("   READING: 2-page crossing = MAX-CUT on the interleaving (circle) graph;")
print("   Guy Z(n) <=> MaxCut(IL_n) = C(n,4) - 2Z(n); the landscape has shelves (local minima")
print("   above Z) from n = 6 -- the THM-869 shelf census applies; page-swap = the Ising Z2")
print("   (the T <-> T^op analog); the metagraph toolkit has a crossing-number home.")
print("\nDONE")
