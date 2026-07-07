#!/usr/bin/env python3
r"""
lrc_paley_bridge_kps_S64.py   (kind-pasteur-2026-07-07-S64, HYP-4927 part 2)

THE PALEY BRIDGE + the composite-14 tournament + the witness-semicircle tournament.

S64 part 1 found: the good-set-leader rule T_good is always the TRANSITIVE consensus
order (leader-on-average is a total order; intertwining trivially true).  The MOD-p rule
is the rich one.  Here:

THEOREM CANDIDATE (the concrete realization of THM-126 Paley-flatness):
   T_modp(AP_{p-1}) = the PALEY TOURNAMENT T_p   (p = 3 mod 4 prime).
   i.e. arrow i->j iff (v_j - v_i) mod p in QR_p, on {1,...,p-1}, IS the circulant
   Paley tournament.  So the LRC density-floor minimizer (the AP, at t=1/p = the p-th
   roots of unity) maps under this pair-statistic cutoff to the tournament H-MAXIMIZER
   (Paley).  The SAME root-of-unity arithmetic that MINIMIZES M (spectral flatness)
   MAXIMIZES H -- THM-126 made a literal tournament-construction, not an analogy.

THE COMPOSITE-14 RULE (S63 catalog #1, CRT): Z/14* ~ Z/6 cyclic (gen 3): the QR
   subgroup is {1,5,9,11? } -- compute; -1 = 13 = 3^3 is a NON-residue, so a Paley-like
   antisymmetry EXISTS on the units.  Non-unit differences (even, or =7) are the
   COMPOSITE OBSTRUCTION (klein-S151's 14=2*7) -- they need a tiebreak = the tournament
   image of the sieve's hard residues.

THE WITNESS-SEMICIRCLE RULE (genuinely cyclic, reads the floor geometry): at the
   maxgap-witness time t* (center of the widest gap over Good), the phases {frac(v_i t*)}
   are points on the circle; i->j iff frac(v_j t*) lies in the clockwise open semicircle
   ahead of frac(v_i t*).  A rotational/local tournament -- HAS 3-cycles, and its
   structure is the density-floor's fingerprint at the optimal time.
"""
from fractions import Fraction as F
import random

# ------------------------------------------------------------------ QR sets
def QRset(p):
    return {(a * a) % p for a in range(1, p)}

# ------------------------------------------------------------------ tournaments
def T_modp(E, p):
    QR = QRset(p); k = len(E); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            d = (E[j] - E[i]) % p
            if d in QR: A[i][j] = 1
    # sanity: antisymmetric off p|d
    return A

def paley_tournament(p):
    QR = QRset(p); A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i == j: continue
            if (j - i) % p in QR: A[i][j] = 1
    return A

# ------------------------------------------------------------------ invariants
def scores(A):
    return tuple(sorted(sum(r) for r in A))

def c3(A):
    k = len(A); n = 0
    for i in range(k):
        for j in range(k):
            if A[i][j]:
                for l in range(k):
                    if A[j][l] and A[l][i]: n += 1
    return n // 3

def ham_paths(A):
    k = len(A); dp = [[0]*k for _ in range(1 << k)]
    for i in range(k): dp[1 << i][i] = 1
    for mask in range(1 << k):
        row = dp[mask]
        for last in range(k):
            c = row[last]
            if not c: continue
            Al = A[last]
            for nxt in range(k):
                if not (mask >> nxt) & 1 and Al[nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    full = (1 << k) - 1
    return sum(dp[full])

def complement(A):
    k = len(A); return [[A[j][i] for j in range(k)] for i in range(k)]

def relabel(A, s):
    k = len(A); return [[A[s[i]][s[j]] for j in range(k)] for i in range(k)]

def canon_key(A):
    """cheap iso invariant: sorted multiset of (out-nbr score profile)."""
    k = len(A); sc = [sum(A[i]) for i in range(k)]
    prof = []
    for i in range(k):
        outp = tuple(sorted(sc[j] for j in range(k) if A[i][j]))
        inp = tuple(sorted(sc[j] for j in range(k) if A[j][i]))
        prof.append((sc[i], outp, inp))
    return tuple(sorted(prof))

def steps(E):
    E = sorted(E); return tuple(E[i+1]-E[i] for i in range(len(E)-1))

# ------------------------------------------------------------------ PART 1: the Paley theorem
print("=" * 98)
print("PART 1 -- THE PALEY BRIDGE:  T_modp(AP_{p-1}) == Paley tournament T_p ?  (p = 3 mod 4)")
print("=" * 98)
for p in (3, 7, 11, 19, 23):
    AP = list(range(1, p))          # {1,...,p-1}, the LRC density-floor minimizer base
    Amp = T_modp(AP, p)
    # the AP {1..p-1} has vertex i = value i+1; differences (j-i) exactly, mod p in 1..p-1
    Tp = paley_tournament(p)         # on p vertices 0..p-1
    # T_modp(AP) is on p-1 vertices; the Paley tournament restricted to nonzero = the
    # SUBTOURNAMENT of T_p induced by {1,...,p-1} (drop vertex 0).  Compare iso invariants.
    Tp_sub = [[Tp[a][b] for b in range(1, p)] for a in range(1, p)]
    same = (canon_key(Amp) == canon_key(Tp_sub))
    exact = (Amp == Tp_sub)
    print(f"  p={p:2d}: T_modp(AP_{p-1}) scores={scores(Amp)}  c3={c3(Amp)}  H={ham_paths(Amp) if p<=13 else 'skip'}"
          f"   == T_p minus vertex0 (exact:{exact}, iso:{same})")
print("  (T_p itself, full p vertices: p=7 -> H=189 c3=14 regular; the LRC AP builds it.)")
for p in (7, 11):
    Tp = paley_tournament(p)
    print(f"    Paley T_{p}: scores={scores(Tp)} c3={c3(Tp)} H={ham_paths(Tp) if p<=13 else 'big'} (canonical H-maximizer)")

# ------------------------------------------------------------------ PART 2: composite-14
print()
print("=" * 98)
print("PART 2 -- THE COMPOSITE-14 CRT RULE: Z/14* structure and the AP {1..13} image")
print("=" * 98)
units14 = [u for u in range(1, 14) if __import__('math').gcd(u, 14) == 1]
# Z/14* cyclic gen 3
g = 3; powers = {}
x = 1
for e in range(6):
    powers[x] = e; x = (x * g) % 14
QR14 = {u for u in units14 if powers[u] % 2 == 0}   # even powers = squares
print(f"  Z/14* = {units14}, cyclic <3>; QR (even powers) = {sorted(QR14)}; -1=13 power {powers[13]} ({'NR' if powers[13]%2 else 'R'})")
def T_crt14(E):
    k = len(E); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            d = (E[j] - E[i]) % 14
            if d in QR14: A[i][j] = 1
            elif d in units14: pass                       # antisymmetric on units (13=NR handles reverse)
            else:                                          # non-unit d: the COMPOSITE obstruction
                # CRT tiebreak: (d mod 2, d mod 7); use mod-7 QR on the 7-part, parity on 2-part
                d7 = d % 7
                if d7 in QRset(7): A[i][j] = 1
                elif d7 != 0: pass
                else:                                      # d == 0 or 7 mod 14: pure 7-obstruction
                    A[i][j] = 1 if (E[j] > E[i]) else 0    # last-resort orientation by value
    # fix antisymmetry for unit non-QR (13-type): if neither set, orient by NR
    for i in range(k):
        for j in range(i+1, k):
            if A[i][j] == 0 and A[j][i] == 0:
                d = (E[j]-E[i]) % 14
                A[i][j] = 1 if d in units14 else (1 if E[j] > E[i] else 0)
                if A[i][j]==0: A[j][i]=1
    return A
AP13 = list(range(1, 14))
Ac = T_crt14(AP13)
# how many pairs fall in the "hard" (non-unit-difference) class?
hard = sum(1 for i in range(13) for j in range(i+1,13) if __import__('math').gcd((AP13[j]-AP13[i])%14 or 14,14)>1)
print(f"  AP {{1..13}} under CRT-14: scores={scores(Ac)} c3={c3(Ac)} H={ham_paths(Ac)}")
print(f"    of C(13,2)={13*12//2} pairs, {hard} have NON-UNIT difference mod 14 (the 2- or 7-divisible")
print(f"    = the composite obstruction, klein-S151's 14=2*7, now VISIBLE as the tournament's irregular arcs)")

# ------------------------------------------------------------------ PART 3: witness-semicircle
print()
print("=" * 98)
print("PART 3 -- THE WITNESS-SEMICIRCLE tournament (genuinely cyclic; reads the floor geometry)")
print("=" * 98)
def maxgap_witness(E, res=200000):
    """t* = a time achieving (near) the widest gap; return phases at t*."""
    best = (-1, None)
    for r in range(res):
        x = (r + .5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        g = max(ph[i+1]-ph[i] for i in range(len(ph)-1))
        g = max(g, ph[0]+1-ph[-1])
        if g > best[0]: best = (g, x)
    return best[1]

def T_witness(E, tstar):
    k = len(E); ph = [(E[i]*tstar) % 1.0 for i in range(k)]
    A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            d = (ph[j] - ph[i]) % 1.0        # clockwise distance i->j
            if 0 < d < 0.5: A[i][j] = 1       # j in clockwise semicircle ahead of i
            elif d > 0.5: A[j][i] = 1
    # antisymmetrize any exact-0.5 ties by value
    for i in range(k):
        for j in range(i+1, k):
            if A[i][j]==A[j][i]:
                A[i][j], A[j][i] = (1,0) if E[i]<E[j] else (0,1)
    return A
zoo = {
    "AP {1..13}": list(range(1,14)),
    "record": [2,4,6,8,10,11,12,13,14,16,18,20,22],
    "GW {1..11,13,24}": list(range(1,12))+[13,24],
    "primes13": [2,3,5,7,11,13,17,19,23,29,31,37,41],
}
for nm, E in zoo.items():
    ts = maxgap_witness(E, 60000)
    A = T_witness(E, ts)
    palin = steps(E) == steps(E)[::-1]
    print(f"  {nm:22s} t*={ts:.4f}  scores={scores(A)} c3={c3(A)} H={ham_paths(A)}  {'[palindromic]' if palin else ''}")

print()
print("=" * 98)
print("PART 4 -- reversal->complement intertwining for the mod-p and CRT-14 rules (all families)")
print("=" * 98)
for label, rule in [("mod7", lambda E: T_modp(E, 7)), ("crt14", T_crt14)]:
    allok = True
    for nm, E in list(zoo.items()) + [("AP8", list(range(1,9)))]:
        Estar = sorted(max(E)+min(E)-e for e in E)
        A = rule(E); As = rule(Estar)
        k = len(E); sigma = list(range(k-1,-1,-1))
        ok = (canon_key(As) == canon_key(relabel(complement(A), sigma)))
        allok &= ok
    print(f"  {label}: reversal(E) -> complement(T) [iso] for all tested families: {allok}")
print()
print("DONE.")
