#!/usr/bin/env python3
"""
lrc14_w24_census_verify_klein.py  --  klein-2026-07-02-S113

SAFE forward step toward unconditional LRC(14): verify (Python, no Lean emission) that the W=24 band WOULD
close -- i.e., every covering 13-family inside [1..24] is lonely (finds a G2/kernel-gate witness). This
quantifies the next band step (CoveringFarLonely 22 -> CoveringFarLonely 24) for whoever builds the Lean data
file, and confirms there is no obstruction in the 22<max<=24 shell.

Reuses the S112 W=22 emitter's covering() + gate_ok() + find_witness() (the kernel-gate mirror). Does NOT emit
Lean; only counts + verifies. If any covering family fails to find a witness, that would be a counterexample
signal (printed loudly).
"""
from itertools import combinations
from math import gcd
import time

def covering(S):
    return all(any(s % q == 0 for s in S) for q in range(2, 15))

def gate_ok(S, a, q):
    for s in S:
        m = (s * a) % q
        if 14 * min(m, q - m) < q:
            return False
    return True

def find_witness(S):
    qs = set()
    for i, u in enumerate(S):
        qs.add(2 * u)
        for v in S[i+1:]:
            qs.add(u + v); qs.add(v - u)
    qs.discard(0)
    for q in sorted(qs):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            if gate_ok(S, a, q):
                return (a, q)
    return None

def census(W):
    t0 = time.time()
    total_cov = 0; with_wit = 0; fails = []; shell = 0; shell_wit = 0
    for S in combinations(range(1, W + 1), 13):
        if not covering(S):
            continue
        total_cov += 1
        in_shell = (S[-1] > 22)          # max entry in {23,24} = retired from CoveringFarLonely 22
        if in_shell: shell += 1
        w = find_witness(S)
        if w is None:
            fails.append(S)
        else:
            with_wit += 1
            if in_shell: shell_wit += 1
    return dict(W=W, total_cov=total_cov, with_wit=with_wit, fails=fails,
               shell=shell, shell_wit=shell_wit, secs=time.time()-t0)

if __name__ == "__main__":
    print("="*72)
    print("W=24 BAND CENSUS VERIFICATION (does the next band close?)")
    print("="*72)
    r = census(24)
    print(f"W={r['W']}: covering families in [1..24] = {r['total_cov']}; "
          f"with G2/kernel-gate witness = {r['with_wit']}; FAILURES = {len(r['fails'])}  ({r['secs']:.0f}s)")
    print(f"  22<max<=24 SHELL (retired from CoveringFarLonely 22 by moving to W=24): {r['shell']} families, "
          f"all lonely = {r['shell']==r['shell_wit']}")
    if r['fails']:
        print("  !!! COVERING FAMILIES WITHOUT A WITNESS (would need a different route) !!!")
        for S in r['fails'][:20]:
            print("   ", S)
    else:
        print("  => EVERY covering family in [1..24] is lonely (G2 gate). The W=24 band is CLOSEABLE.")
        print(f"  => building LRCWindowData24.lean (native_decide over C(24,13)=2496144) would shrink the")
        print(f"     remaining hypothesis from CoveringFarLonely 22 to CoveringFarLonely 24 ({r['shell']} families retired).")
    print("\nHONEST: this verifies the FINITE window shell only; the infinite far tail (max>24) remains")
    print("CoveringFarLonely -- the band shrinks but never closes it. The closure route is the peel/rate")
    print("descent (kind-pasteur peel20 + the proved rate_core), not band extension.")
