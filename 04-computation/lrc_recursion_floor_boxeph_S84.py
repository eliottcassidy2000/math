#!/usr/bin/env python3
"""
The sharp descent recursion and the sub-factor-2 floor  (boxeph-2026-07-18-S84)
===============================================================================

THM-1008 used the TIGHT base M(V')>=1/13, losing a factor 2 to 1/26 at rho=1.
Using the ACTUAL M of the kept set gives the sharp recursion (elementary, same
one-runner perturbation):
    M(V) >= rho * M(V \ {r}) / (rho_r + 1),   rho_r = v_r / max(V\{r}),
removing runner r and fixing it back.  Take the BEST r.  This beats 1/26 unless
the best kept-12-subset is TIGHT (M=1/13), i.e. a dilated AP {d,..,12d} -- which
forces V to be a dilated tight family (non-primitive, dilation-reducible).

This script computes, for compact covering families (rho<13, covers 2..14):
  - M(V) exactly (pair-sum bound), and M(V\{r}) for every r,
  - the best-removal recursion floor  max_r rho_r M(V\{r})/(rho_r+1),
  - whether it reaches 1/14  (=> the recursion alone closes the compact residual?),
  - the equality/near-tight structure.
"""
from fractions import Fraction as F
from math import gcd
import itertools, random

def norm(x):
    r = x % 1
    return min(r, 1 - r)

def exact_M(V):
    """M(V)=max_t min_v||v t||, over t=a/q with q|v_i+v_j, q<=2max (THM-999)."""
    if len(V) == 1:
        return F(1, 2)  # single runner: t=1/(2v)
    best = F(0)
    qs = set([14])
    for i in range(len(V)):
        for j in range(i, len(V)):
            s = V[i] + V[j]
            for d in range(1, s + 1):
                if s % d == 0:
                    qs.add(d)
    for q in qs:
        for a in range(1, q):
            if gcd(a, q) == 1:
                m = min(min((v * a) % q, q - (v * a) % q) for v in V)
                c = F(m, q)
                if c > best:
                    best = c
    return best

def is_covering(V, n=14):
    return all(any(v % b == 0 for v in V) for b in range(2, n + 1))

def recursion_floor_best(V):
    """max_r rho_r * M(V\{r}) / (rho_r+1), rho_r = v_r/max(V\{r})."""
    best = F(0); bestr = None
    for idx, r in enumerate(V):
        kept = V[:idx] + V[idx+1:]
        Bk = max(kept)
        rho = F(r, Bk)
        Mk = exact_M(kept)
        fl = rho * Mk / (rho + 1)
        if fl > best:
            best = fl; bestr = (r, Mk, rho)
    return best, bestr

def base_floor(V):
    vs = sorted(V); return F(vs[-1], 13*(vs[-1]+vs[-2]))  # THM-1008

def is_dilated_ap(W):
    """Is W a dilated AP {d,2d,...,kd}?"""
    Ws = sorted(W); d = Ws[0]
    return all(Ws[i] == (i+1)*d for i in range(len(Ws)))

FAMILIES = [
    ("kps floor-min M=1/9", [3,4,11,12,13,15,18,20,24,42,55,64,67]),
    ("{2..14}",              list(range(2,15))),
    ("{1,3,4..14}",          [1]+list(range(3,15))),
    ("{1,2,4,5..14}",        [1,2]+list(range(4,15))),
    ("2*{1..13} dilated tight", [2*i for i in range(1,14)]),
    ("{1..10,13,22,84} multi", list(range(1,11))+[13,22,84]),
    ("residue {1..11,13,84}", list(range(1,12))+[13,84]),
]

print("="*100)
print("SHARP DESCENT RECURSION  M(V) >= rho * M(V\\{r})/(rho_r+1)  — best removal r")
print("="*100)
allbeat = True; allreach = True
for name, V in FAMILIES:
    V = sorted(V)
    M = exact_M(V)
    bf = base_floor(V)
    rf, info = recursion_floor_best(V)
    cov = is_covering(V)
    rho = F(V[-1], V[-2])
    beats26 = rf > bf
    reach14 = rf >= F(1,14)
    allbeat &= (beats26 or rf==bf)  # dilated tight may equal
    allreach &= reach14
    r_, Mk_, rho_ = info
    print(f"\n{name}")
    print(f"  cover={cov} rho={float(rho):.2f}  M={M} ({float(M):.4f})")
    print(f"  base floor (1008)={bf} ({float(bf):.4f})   RECURSION floor={rf} ({float(rf):.4f})")
    print(f"     best removal r={r_} (M(kept)={Mk_}={float(Mk_):.4f}, rho_r={float(rho_):.2f})")
    print(f"  beats 1/26-scale base: {beats26}   REACHES 1/14: {reach14}")

# stress: random compact covering families
print("\n" + "="*100)
print("RANDOM COMPACT COVERING families (rho<13): does best-removal recursion reach 1/14?")
print("="*100)
random.seed(7)
found = 0; reached = 0; fails = []
attempts = 0
while found < 20 and attempts < 20000:
    attempts += 1
    V = sorted(random.sample(range(1, 45), 13))
    if not is_covering(V):
        continue
    if F(V[-1], V[-2]) >= 13:
        continue  # skip dominant
    found += 1
    rf, info = recursion_floor_best(V)
    r14 = rf >= F(1,14)
    reached += r14
    if not r14:
        fails.append((V, rf, exact_M(V)))
print(f"compact covering sampled: {found}; recursion reached 1/14: {reached}/{found}")
for V, rf, M in fails[:6]:
    print(f"  FAIL rf={float(rf):.4f} M={float(M):.4f}  V={V}")
print("\nREADING: if recursion reaches 1/14 on all compact covering families, the")
print("sub-factor-2 floor = the recursion + a floor on M(V\\{r}); equality only at")
print("dilated tight APs (M(kept)=1/13), which are non-primitive/dilation-reducible.")
