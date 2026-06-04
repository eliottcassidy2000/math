#!/usr/bin/env python3
"""
S623 — the CONSTRUCTIVE route to C'(n) at small n, to illuminate n=14.
C'(n): a primitive speed set S of n-1 distinct positive ints with n|v for some v is LOOSE (M>1/n).
Constructive criteria:
  (1-clock)  t=1/m, m in {2..n-1}: works iff no s_i ≡ 0 (mod m)            [M >= 1/m > 1/n]
  (twist)    t=a/m, a in (Z/m)*, n+1<=m<=2n-1: works iff all a*s_i mod m off band (n*dist>m)
  (B')       some multiple v=nw: G(S\{v}) has a component longer than 2/(n v)   [dominance/measure]
We EXHAUSTIVELY enumerate primitive multiple-of-n configs in a box and record, for each, the
MINIMAL m that escapes (1-clock or twist), and fall back to B' otherwise.  Focus: n=5 (2n-1=9=3^2),
the smallest prime-power shell case = baby version of n=14 (3^3).
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_n14_flowshells_s622 import gap_and_argmax, safe_components

def dodge_min_m(S, n):
    """minimal m (2..2n-1) with a twisted/1-clock escape; returns (m, a, kind) or None."""
    for m in range(2, 2*n):
        for a in range(1, m):
            if gcd(a, m) != 1: continue
            if all(n*min((a*s) % m, m-((a*s) % m)) > m for s in S):
                kind = '1clock' if m <= n-1 else 'twist'
                return (m, a, kind)
    return None

def bprime(S, n):
    lvl = Fr(1, n)
    for v in S:
        if v % n: continue
        Sp = [x for x in S if x != v]
        comps = safe_components(Sp, lvl)
        if max((b-a for a, b in comps), default=Fr(0)) > Fr(2, n*v):
            return True
    return False

def study(n, M):
    lvl = Fr(1, n); k = n-1
    total = covered = twist_used = bprime_only = 0
    mdist = {}; resid = []; shell_examples = {}
    for S in itertools.combinations(range(1, M+1), k):
        if not any(v % n == 0 for v in S): continue
        if reduce(gcd, S) != 1: continue
        total += 1
        d = dodge_min_m(S, n)
        if d:
            m, a, kind = d
            mdist[m] = mdist.get(m, 0) + 1
            if kind == 'twist':
                twist_used += 1
                shell_examples.setdefault(m, (S, a))
            covered += 1
        elif bprime(S, n):
            bprime_only += 1; covered += 1
        else:
            resid.append(S)
    print(f"\n n={n} (2n-1={2*n-1}), box<= {M}, {k} speeds:  {total} primitive multiple-of-{n} configs")
    print(f"   covered by dodge(m<=2n-1) OR B': {covered}/{total}   TRUE residual: {len(resid)} {resid[:3]}")
    print(f"   needed a deep-shell twist: {twist_used};  B'-only (no dodge<=2n-1): {bprime_only}")
    print(f"   minimal-escaping-m distribution: {dict(sorted(mdist.items()))}")
    for m in sorted(x for x in mdist if x > n-1):
        S, a = shell_examples[m]
        print(f"     twist shell m={m}: e.g. a={a}, S={S}, residues mod {m}={tuple(s%m for s in S)}")
    return resid

if __name__ == "__main__":
    for n, M in [(4, 24), (5, 26), (6, 22), (7, 20), (8, 19)]:
        study(n, M)
