#!/usr/bin/env python3
"""
klein-2026-07-20-S339 -- THE PFAFFIAN IS THE ODD FUNCTION AND det IS THE EVEN ONE;
plus the GLOBAL Pfaffian spectrum (the hp-spectrum question for a second odd-valued invariant).

Owner: "each odd sized tournament corresponds to one of the natural numbers.  chase the high
leverage question, see the relation between odd valued functions and tournament adjacent ideas.
they both relate also to even concepts like even graphs and even functions."

WHAT IS ALREADY DONE BY OTHERS (not repeated here):
  - the canonical Cut (+) Cycle splitting at ODD n (bicycle space = 0 iff n odd) -- another
    agent's THM-1440 (seidel-spectra), 09:17 today.  That is "each odd-sized tournament IS a
    natural number": T <-> (tiling number in [0, 2^C(n-1,2)), score-parity vector).
  - kind-pasteur THM-1455: the LOCAL Pfaffian expansion [x^{n-2k}]char(S) = sum_{|A|=2k} Pf(S_A)^2,
    and a mod-16 law resting on k_4 = #{4-subsets with |Pf| = 3}.

WHAT IS NEW HERE:
  (1) Pf(S(T)) is ODD for every tournament of even order.  One line: S = J - I mod 2, so
      Pf(S) = Pf(J-I) = #{perfect matchings of K_n} = (n-1)!!, a product of odds.
      This EXPLAINS kind-pasteur's |Pf| in {1,3} at 4-subsets and generalises it: for a
      4-subset det = Pf^2 <= 9 and Pf odd, so Pf in {+-1,+-3} -- their dichotomy is forced.
  (2) Pf is an ODD (ALTERNATING) function and det = Pf^2 is an EVEN (invariant) one:
      Pf(P^T S P) = det(P) Pf(S) = sign(sigma) Pf(S).  That is the literal odd-function /
      even-function dichotomy the owner named, carried by the sign character of S_n.
  (3) COROLLARY: since Pf(S) is odd it is NONZERO, so no tournament of even order admits an
      automorphism of odd SIGN: Aut(T) <= A_n.  (Classical via Moon's |Aut| odd; this is an
      independent one-line route.)
  (4) THE GLOBAL PFAFFIAN SPECTRUM -- the hp-spectrum question for a second odd-valued
      invariant.  hp attains {odds} minus {7,21}.  Which odd values does Pf attain?  GAPS?
  (5) hp vs Pf at even n: two odd-valued invariants on the same object -- any congruence?
"""
import itertools, random
from math import comb
from fractions import Fraction

def pairs_of(n): return [(i, j) for i in range(n) for j in range(i + 1, n)]
def om_from_bits(n, bits, P):
    om = [0] * n
    for k, (i, j) in enumerate(P):
        if bits >> k & 1: om[i] |= 1 << j
        else:             om[j] |= 1 << i
    return tuple(om)
def skew(om, n):
    S = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j: S[i][j] = 1 if (om[i] >> j & 1) else -1
    return S
def pf(M):
    n = len(M)
    if n == 0: return 1
    if n % 2: return 0
    tot = 0
    for j in range(1, n):
        idx = [k for k in range(1, n) if k != j]
        sub = [[M[a][b] for b in idx] for a in idx]
        tot += ((-1) ** (j - 1)) * M[0][j] * pf(sub)
    return tot
def hp(om, n):
    cnt = 0
    def go(last, used, d):
        nonlocal cnt
        if d == n: cnt += 1; return
        mv = om[last] & ~used
        while mv:
            b = mv & -mv; w = b.bit_length() - 1; mv ^= b
            go(w, used | (1 << w), d + 1)
    for s in range(n): go(s, 1 << s, 1)
    return cnt
def relabel(om, perm, n):
    new = [0] * n
    for v in range(n):
        mv, t = om[v], 0
        while mv:
            b = mv & -mv; w = b.bit_length() - 1; mv ^= b
            t |= 1 << perm[w]
        new[perm[v]] = t
    return tuple(new)
def sgn(p):
    p = list(p); s = 1; seen = [False] * len(p)
    for i in range(len(p)):
        if seen[i]: continue
        j, L = i, 0
        while not seen[j]: seen[j] = True; j = p[j]; L += 1
        if L % 2 == 0: s = -s
    return s

print("=" * 78)
print("(2) Pf IS AN ODD (ALTERNATING) FUNCTION, det = Pf^2 IS AN EVEN (INVARIANT) ONE")
print("=" * 78)
random.seed(5)
for n in (4, 6):
    P = pairs_of(n); E = len(P); ok_pf = ok_det = True
    for _ in range(200):
        om = om_from_bits(n, random.getrandbits(E), P)
        S = skew(om, n); p0 = pf(S); d0 = p0 * p0
        perm = list(range(n)); random.shuffle(perm)
        om2 = relabel(om, perm, n); S2 = skew(om2, n)
        # perm as a permutation matrix acting: sigma sends v -> perm[v]
        if pf(S2) != sgn(perm) * p0: ok_pf = False
        if pf(S2) ** 2 != d0: ok_det = False
    print(f" n={n}: Pf(relabelled) = sign(sigma)*Pf : {ok_pf}    det invariant : {ok_det}")
print(" => Pf transforms by the SIGN CHARACTER (odd/alternating); det = Pf^2 is invariant (even).")

print("\n" + "=" * 78)
print("(3) COROLLARY: Pf odd => Pf != 0 => Aut(T) <= A_n for n even")
print("=" * 78)
for n in (4, 6):
    P = pairs_of(n); E = len(P); bad = 0; tested = 0
    for _ in range(300):
        om = om_from_bits(n, random.getrandbits(E), P)
        for perm in itertools.permutations(range(n)):
            if relabel(om, list(perm), n) == om:
                tested += 1
                if sgn(perm) == -1: bad += 1
    print(f" n={n}: automorphisms found {tested}, of ODD sign: {bad}   (must be 0)")

print("\n" + "=" * 78)
print("(4) THE GLOBAL PFAFFIAN SPECTRUM -- gaps?  (hp attains {odds} minus {7,21})")
print("=" * 78)
for n in (4, 6):
    P = pairs_of(n); E = len(P); vals = set()
    for bits in range(1 << E):
        vals.add(pf(skew(om_from_bits(n, bits, P), n)))
    pos = sorted(v for v in vals if v > 0)
    mx = max(pos); missing = [k for k in range(1, mx + 1, 2) if k not in vals]
    print(f" n={n} EXHAUSTIVE over {1<<E} tournaments:")
    print(f"    |Pf| values attained: {pos}   max = {mx}")
    print(f"    all odd? {all(v % 2 == 1 for v in pos)}   symmetric +-? {set(-v for v in pos) <= vals}")
    print(f"    MISSING odd values below the max: {missing if missing else 'NONE -- no gaps'}")
for n in (8,):
    P = pairs_of(n); E = len(P); vals = set()
    for _ in range(400000):
        vals.add(pf(skew(om_from_bits(n, random.getrandbits(E), P), n)))
    pos = sorted(v for v in vals if v > 0); mx = max(pos)
    missing = [k for k in range(1, mx + 1, 2) if k not in vals]
    print(f" n={n} SAMPLED (400k of 2^{E}):")
    print(f"    |Pf| attained: {pos}   max seen = {mx}")
    print(f"    MISSING odd below max: {missing if missing else 'NONE in sample'}")

print("\n" + "=" * 78)
print("(5) hp vs Pf: TWO ODD-VALUED INVARIANTS ON THE SAME TOURNAMENT -- any congruence?")
print("=" * 78)
for n in (4, 6):
    P = pairs_of(n); E = len(P)
    pairsHP = []
    for bits in range(1 << E):
        om = om_from_bits(n, bits, P)
        pairsHP.append((hp(om, n), abs(pf(skew(om, n)))))
    for m in (4, 8, 16):
        rel = {(h % m, p % m) for (h, p) in pairsHP}
        det = all(len({p for (hh, p) in rel if hh == a}) == 1 for a in {x for (x, _) in rel})
        print(f" n={n} mod {m}: distinct (hp,|Pf|) residue pairs = {len(rel)}"
              f"   hp residue DETERMINES |Pf| residue: {det}")
    print(f"    joint support mod 8: {sorted({(h%8,p%8) for (h,p) in pairsHP})}")
print("\nDONE.")
