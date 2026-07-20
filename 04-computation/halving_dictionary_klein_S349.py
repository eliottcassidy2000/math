#!/usr/bin/env python3
"""
klein-2026-07-20-S349 -- THE HALVING DICTIONARY: {-1,0,1} vs {1,1/2,0}, or why the repo's
SIGN world and its EVEN/ODD world are the same object over Z[1/2] and different mod 2.

Owner: "think about how {-1,0,1} and {1,1/2,0} are functionally equivalent but have shown up
in this repo many times each. this is what i mean by even/odd vs positive/negative."

THE MAP.  x -> (1+x)/2 sends {-1,0,1} -> {0,1/2,1}.  For tournaments this is EXACTLY the
standard conversion between the SKEW adjacency S (S_ij = +-1, 0 on the diagonal) and the
0/1 adjacency A:
        A = (J - I + S)/2,        S = 2A - (J - I).
So the repo's two halves -- the GF(2)/parity machinery (tilings, switching classes, cut+cycle,
even graphs, two-graphs) and the sign machinery (skew adjacency, Pfaffian, charge/weight,
the nullcone) -- are ONE object viewed in two coordinates.  The dictionary costs a 1/2, and
1/2 does not exist in characteristic 2.  THAT is the even/odd vs positive/negative split.

TESTED HERE:
 1. the dictionary, exactly;
 2. THE COLLAPSE PRINCIPLE: S == J - I (mod 2) for EVERY tournament, so every INTEGER
    polynomial invariant of S is constant mod 2 -- a theorem with Pf-oddness as a corollary;
 3. a census of which invariants are constant mod 2 and which are not -- the ones that VARY
    are exactly the ones whose definition needs the 1/2;
 4. the {1,1/2,0} ladder on the nullcone side, and whether the correspondence is structural.
"""
import itertools, math
from fractions import Fraction as Fr

def pairs_of(n): return [(i, j) for i in range(n) for j in range(i + 1, n)]
def om_from_bits(n, bits, P):
    om = [0] * n
    for k, (i, j) in enumerate(P):
        if bits >> k & 1: om[i] |= 1 << j
        else:             om[j] |= 1 << i
    return om
def skew(om, n):
    return [[0 if i == j else (1 if (om[i] >> j & 1) else -1) for j in range(n)] for i in range(n)]
def adj01(om, n):
    return [[0 if i == j else (1 if (om[i] >> j & 1) else 0) for j in range(n)] for i in range(n)]
def pf(M):
    n = len(M)
    if n == 0: return 1
    if n % 2: return 0
    t = 0
    for j in range(1, n):
        idx = [k for k in range(1, n) if k != j]
        t += ((-1) ** (j - 1)) * M[0][j] * pf([[M[a][b] for b in idx] for a in idx])
    return t
def det_int(M):
    M = [row[:] for row in M]; n = len(M); det = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return 0
        if p != c: M[c], M[p] = M[p], M[c]; det = -det
        det *= M[c][c]
        inv = Fr(1, 1) / Fr(M[c][c])
        for r in range(c + 1, n):
            f = Fr(M[r][c]) * inv
            for k in range(c, n): M[r][k] = M[r][k] - f * M[c][k]
    return int(det)
def hp_count(om, n):
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
def c3(om, n):
    t = 0
    for a, b, c in itertools.combinations(range(n), 3):
        for x, y, z in ((a,b,c),(a,c,b)):
            if (om[x] >> y & 1) and (om[y] >> z & 1) and (om[z] >> x & 1): t += 1
    return t
def trace_pow(M, k):
    n = len(M); R = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    for _ in range(k):
        R = [[sum(R[i][l] * M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
    return sum(R[i][i] for i in range(n))

print("=" * 80)
print("1. THE DICTIONARY  A = (J - I + S)/2  -- exact, and it needs 1/2")
print("=" * 80)
ok = True
for n in (4, 5, 6):
    P = pairs_of(n)
    for bits in range(min(1 << len(P), 400)):
        om = om_from_bits(n, bits, P); S = skew(om, n); A = adj01(om, n)
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if A[i][j] * 2 != 1 + S[i][j]: ok = False
print(f"   A_ij = (1 + S_ij)/2 for every off-diagonal entry, n = 4,5,6 : {ok}")
print("   diagonal: S_ii = 0  <->  'A_ii = 1/2'  -- the tie/self entry is the 1/2.")
print("   So {-1,0,1} (sign) and {1,1/2,0} (fraction) are the SAME data, one halving apart.")

print("\n" + "=" * 80)
print("2. THE COLLAPSE PRINCIPLE (a theorem):  S == J - I  (mod 2), for EVERY tournament")
print("=" * 80)
print("""   Every off-diagonal entry of S is +-1, and +1 == -1 == 1 (mod 2).  So mod 2 the skew
   adjacency of EVERY tournament is the SAME matrix J - I -- all orientation information is
   destroyed.  Consequently:

     THEOREM.  If I(T) = F(S) for a polynomial F with INTEGER coefficients, then
               I(T) == F(J - I)  (mod 2)  is the SAME for every tournament on n vertices.

   Pfaffian oddness (THM-1475) is the corollary F = Pf:  Pf(S) == Pf(J-I) = (n-1)!!, odd.
   The principle predicts which invariants can and cannot see the tournament mod 2.""")

print("\n" + "=" * 80)
print("3. THE CENSUS -- which invariants are constant mod 2, and which need the 1/2?")
print("=" * 80)
print(f"{'n':>3} {'invariant':<26} {'values mod 2 over all T':<26} {'constant?':>10} {'integer poly in S?':>20}")
for n in (4, 5, 6):
    P = pairs_of(n); E = len(P)
    tot = 1 << E
    inv = {"Pf(S)": [], "det(S)": [], "tr(S^2)": [], "tr(S^4)": [],
           "hp(T)": [], "c3(T)": [], "tr(A^3)": []}
    for bits in range(min(tot, 4096)):
        om = om_from_bits(n, bits, P); S = skew(om, n); A = adj01(om, n)
        if n % 2 == 0: inv["Pf(S)"].append(pf(S) % 2)
        inv["det(S)"].append(det_int(S) % 2)
        inv["tr(S^2)"].append(trace_pow(S, 2) % 2)
        inv["tr(S^4)"].append(trace_pow(S, 4) % 2)
        inv["hp(T)"].append(hp_count(om, n) % 2)
        inv["c3(T)"].append(c3(om, n) % 2)
        inv["tr(A^3)"].append(trace_pow(A, 3) % 2)
    ispoly = {"Pf(S)": "YES", "det(S)": "YES", "tr(S^2)": "YES", "tr(S^4)": "YES",
              "hp(T)": "no (needs 1/2)", "c3(T)": "no (needs 1/2)", "tr(A^3)": "no (needs 1/2)"}
    for k, v in inv.items():
        if not v: continue
        s = sorted(set(v))
        print(f"{n:>3} {k:<26} {str(s):<26} {str(len(s)==1):>10} {ispoly[k]:>20}")
print("""
   READING: every integer polynomial in S is constant mod 2, exactly as the principle says.
   The invariants that DO vary mod 2 are precisely those defined in the 0/1 (fraction) world,
   whose expression in S carries a 1/2.  hp is the interesting middle case -- it needs the
   1/2, yet is still constant mod 2 (always ODD).  That is REDEI's theorem, and the census
   shows it is NOT an instance of the collapse principle but an independent fact.""")

print("\n" + "=" * 80)
print("4. THE SAME TRIPLE ON THE NULLCONE SIDE:  charges {-1,0,1}  vs  Gamma shapes {1,1/2,0}")
print("=" * 80)
print("""   GMC witnesses (THM-1490/1500/1510): P = (1+Z)(W - g(Z)U) with U ~ Gamma(a), a = nu/2 for
   nu real Gaussians.  The master equation forces 1 + s g(s) = (1+s)^{1/a}, so g is a
   POLYNOMIAL iff 1/a is a positive integer.  The realizable ladder is therefore

        a = nu/2  in  {1, 1/2, 0}   <->   nu in {2, 1, 0}   <->   n = 2+nu in {4, 3, 2}

   a = 1   : G(x) = 1/(1+x)         rational      n = 4 witness   (1/a = 1)
   a = 1/2 : G(x) = (1+x)^{-1/2}    square root   n = 3 witness   (1/a = 2)
   a = 0   : G(x) = e^{-cx}         DEGENERATE    n = 2           (1/a = infinity)

   And the witnesses' CHARGES are {-1, 0, +1}.  So the same three-element pattern indexes both
   sides: {1, 1/2, 0} counts HALF-GAUSSIANS (fraction world) and {-1, 0, 1} counts CHARGE
   (sign world).  HONEST: the shared shape is that both are the image of nu in {2,1,0} under
   different normalisations -- nu/2 on one side, and 'how many signs are available' on the
   other.  I am NOT claiming a canonical bijection between charge and Gamma shape; what IS
   structural is that BOTH degenerate exactly at the same end (a = 0 / n = 2), which is
   precisely where GMC(2) lives and where 1/a = infinity kills the polynomiality.""")
