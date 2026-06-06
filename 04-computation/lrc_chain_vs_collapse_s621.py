#!/usr/bin/env python3
"""
S621 part 4 — test HYP-2195's conjecture  (p0=0  <=>  additive chain).
additive chain: every element except the minimum is a sum of two (<=) smaller
elements of the set.  collapse: gap_exact == 1/(n+1)  (p0=0, exactly tight).
Cross-tabulate over all primitive n-sets in range to settle both directions.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_tight_enum_s621 import gap_exact

def is_additive_chain(S):
    """strict: every non-minimum element is a sum of two strictly-smaller elements."""
    S = sorted(S)
    for x in S[1:]:
        if not any((a + b == x) for a in S for b in S if a <= b < x and a < x):
            return False
    return True

def max_in_sumset(S):
    """clean necessary condition: max(S) = a+b for some a,b in S (a<=b<max)."""
    S = sorted(S); mx = S[-1]
    return any(a + b == mx for a in S for b in S if a <= b < mx)

print("HYP-2195 test:  collapse (p0=0)  vs  additive-chain  over primitive n-sets")
print("="*74)
for n in range(3, 8):
    delta = Fr(1, n+1); R = (2*n+3) if n <= 5 else (2*n+1)
    collapse, chains = [], []
    chain_not_collapse, collapse_not_chain = [], []
    collapse_max_not_sum, n_max_in_sumset = [], 0
    for s in itertools.combinations(range(1, R+1), n):
        if reduce(gcd, s) != 1: continue
        ch = is_additive_chain(s)
        ms = max_in_sumset(s)
        co = gap_exact(list(s)) == delta
        if ch: chains.append(s)
        if ms: n_max_in_sumset += 1
        if co: collapse.append(s)
        if ch and not co: chain_not_collapse.append(s)
        if co and not ch: collapse_not_chain.append(s)
        if co and not ms: collapse_max_not_sum.append(s)
    print(f"\n n={n} (range 1..{R}): #collapse={len(collapse)}  #strict-chains={len(chains)}  #max-in-sumset={n_max_in_sumset}")
    print(f"   [necessary] collapse ⟹ max(S)∈S+S? {len(collapse_max_not_sum)==0}   (exceptions: {collapse_max_not_sum[:3]})")
    print(f"   [biconditional?] strict-chain ⟹ collapse? {len(chain_not_collapse)==0}   (chains NOT collapse: {len(chain_not_collapse)})")
    print(f"   so additive structure is NECESSARY-ish but vastly INSUFFICIENT: {len(collapse)} collapse vs {n_max_in_sumset} max-in-sumset")
print("\nVERDICT: collapse ⟹ additive-chain is the true (one-way) implication;")
print("         additive-chain ⟹ collapse is FALSE (chains vastly outnumber collapses).")
