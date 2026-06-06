#!/usr/bin/env python3
"""
S630 part 2 — quantitative C'(19): is min M over multiple-of-19 configs = 2/37 (THM-415 pattern,
unramified 37), and does shell-dodge(m<=37) ∪ B' cover ALL multiple-of-19 configs (residual empty,
like n=14 in S622)?  Strong evidence for LRC(19).
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, random, sys
sys.path.insert(0, '04-computation')
from lrc_fast_s628 import gap_shells, good_mult_fast

N = 19; LVL = Fr(1, N)

def twisted_dodge(S, Mmax=37):
    """loose if some a/m (m<=Mmax) puts all runners off the central band (n*dist>m)."""
    for m in range(2, Mmax+1):
        for a in range(1, m):
            if gcd(a, m) != 1: continue
            if all(N*min((a*v) % m, m-((a*v) % m)) > m for v in S):
                return (m, a)
    return None

def dominant(S):
    S = sorted(S); return S[-1] > (N-1)*S[-2]

if __name__ == "__main__":
    print("[1] min M over multiple-of-19 configs — is it 2/37? (gap_shells, exact for near-tight)")
    print(f"    2/(2n-1)=2/37={float(Fr(2,37)):.5f}   1/19={float(LVL):.5f}")
    best = (Fr(1), None)
    # targeted near-AP search: AP {1..18} with subsets shifted by +37 (lift) / one made mult of 19
    base = list(range(1, 19))
    cands = []
    # AP variants: replace element j by j+37k or by a multiple of 19
    for j in range(18):
        for repl in [19, 38, base[j]+37, base[j]+2*37]:
            S = base[:]; S[j] = repl
            if len(set(S)) == 18: cands.append(tuple(sorted(S)))
    # plus random multiple-of-19 near-tight
    rng = random.Random(1)
    for _ in range(20000):
        S = set(base); S.discard(rng.choice(base)); S.add(19*rng.randint(1,2))
        while len(S) < 18: S.add(rng.randint(1, 40))
        if len(S) == 18: cands.append(tuple(sorted(S)))
    for S in cands:
        if not any(v % 19 == 0 for v in S): continue
        if reduce(gcd, S) != 1: continue
        g = gap_shells(list(S), N)              # exact for near-tight; lower bound else
        if g < best[0]: best = (g, S)
    print(f"    min gap_shells found = {best[0]}={float(best[0]):.5f}  ==2/37? {best[0]==Fr(2,37)}  at {best[1]}")

    print("\n[2] residual: shell-dodge(m<=37) ∪ B'(dominant) over multiple-of-19 configs")
    tested = covered = dodge_only = dom_only = 0; resid = []
    for _ in range(20000):
        S = {19*rng.randint(1, 2)}
        while len(S) < 18: S.add(rng.randint(1, 45))
        S = tuple(sorted(S))
        if len(S) != 18 or reduce(gcd, S) != 1: continue
        tested += 1
        d = twisted_dodge(S, 37) is not None
        b = dominant(S)
        if d or b: covered += 1
        else: resid.append(S)
        if d and not b: dodge_only += 1
        if b and not d: dom_only += 1
    print(f"    tested {tested}: covered {covered}  residual {len(resid)}")
    print(f"    dodge(m<=37)-only {dodge_only}  B'-only {dom_only}")
    for S in resid[:5]:
        g = gap_shells(list(S), N, Mmax=60)
        print(f"      RESID {S}  gap_shells(m<=60)={g}={float(g):.5f} (>1/19? {g>LVL})")
    print(f"\n    => shell-dodge(m<=2n-1=37) ∪ B' covers {100*covered/max(tested,1):.2f}% of multiple-of-19 configs")
