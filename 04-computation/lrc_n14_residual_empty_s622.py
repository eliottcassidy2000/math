#!/usr/bin/env python3
"""
S622 — the decisive test:  TWISTED-SHELL DODGE (m<=27) + LEMMA B (dominance) cover
every multiple-of-14 config?  A config in neither's reach is the TRUE residual.

  Lemma B (THM-398 Crit B'): S loose if G(S\{v}) has a component longer than v's arc 2/(n v).
  Twisted-shell dodge: S loose if some a/m, 2<=m<=27, has all a*s_i off the band (14*dist>m).

Report any config evading BOTH (the genuine open core), and its exact M.
"""
from fractions import Fraction as Fr
from math import gcd
import random, sys
sys.path.insert(0, '04-computation')
from lrc_tight_enum_s621 import norm
from lrc_n14_flowshells_s622 import gap_and_argmax, safe_components

N = 14; LVL = Fr(1, N)

def twisted_escape(S, mmax=27):
    for m in range(2, mmax+1):
        for a in range(1, m):
            if gcd(a, m) != 1: continue
            if all(14*min((a*s) % m, m-((a*s) % m)) > m for s in S):
                return (m, a)
    return None

def lemmaB(S):
    """does some multiple v=14w have a safe component of S\\{v} longer than 2/(n v)?"""
    for v in S:
        if v % N: continue
        Sp = [x for x in S if x != v]
        comps = safe_components(Sp, LVL)
        longest = max((b-a for a, b in comps), default=Fr(0))
        if longest > Fr(2, N*v):
            return True
    return False

def rand_cfg(rng):
    while True:
        S = {14*rng.randint(1,2)}
        for m in rng.sample([2,3,4,5,6,7,8,9,10,11,12,13], rng.randint(7,12)):
            S.add(m*rng.randint(1,3))
        while len(S) < 13: S.add(rng.randint(1,45))
        S = tuple(sorted(S))
        if len(S)==13:
            g=0
            for x in S: g=gcd(g,x)
            if g==1: return S

if __name__ == "__main__":
    rng = random.Random(11)
    N_CFG = 6000
    residual = []; only_lemmaB = 0; only_twist = 0; both = 0
    for _ in range(N_CFG):
        S = rand_cfg(rng)
        t = twisted_escape(S) is not None
        b = lemmaB(S)
        if t and b: both += 1
        elif t: only_twist += 1
        elif b: only_lemmaB += 1
        else:
            M, arg = gap_and_argmax(S)
            residual.append((S, str(M), float(M), [str(x) for x in arg]))
    print(f"tested {N_CFG} multiple-of-14 configs:")
    print(f"  covered by BOTH criteria: {both}")
    print(f"  twisted-shell(m<=27) only: {only_twist}")
    print(f"  Lemma B only:              {only_lemmaB}")
    print(f"  TRUE RESIDUAL (neither):   {len(residual)}")
    for r in residual[:8]:
        print(f"    {r[0]}  M={r[1]}={r[2]:.5f}  loose={r[2]>1/14}  t*={r[3]}")
    cov = both + only_twist + only_lemmaB
    print(f"\n=> twisted-shell(m<=27) ∪ Lemma B covers {cov}/{N_CFG} = {100*cov/N_CFG:.2f}% of multiple-of-14 configs")
