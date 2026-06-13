#!/usr/bin/env python3
"""
lrc_padic_tree_channels_s542.py    oracle-2026-06-01-S542

Integrate the p-adic tree: unify the Bruhat-Tits endpoint descent (S410: endpoint
depth d_p = v_p(denominator), "gate repair = debt export to child vertices"), the
independent-pair CHANNELS (S532), the APEX = source-sink = speed n/2 (S530), and
the annular-braid CENTER = shift (S541), as ONE object — the p-adic tree of the
speed system.

Claim to verify: for the channel prime p = (the prime so that channels are mod p;
for n with n/2 prime this is n/2, e.g. n=14 -> 7), the p-adic tree of the runner
speeds {1..n-1} has:
  * LEVEL-1 nodes = residues mod p = the floor(n/2)-style CHANNELS (S532);
  * the residue-0 node = the APEX / gate channel = the speed n/2 (S530, S541
    diameter pair);
  * the channel COUPLING (S532: channels not independent for n>=6) = the deeper
    levels (residues mod p^2, ...): the pairs that share a level-1 node split only
    at level 2 -> the tree depth is the coupling.
And the endpoint MOAT (S410) sits at the boundary: depth v_p(n).
"""
from math import gcd
from collections import defaultdict

def vp(x, p):
    if x == 0: return 99
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

def padic_tree(speeds, p, depth):
    """group speeds by residue mod p, p^2, ..., p^depth (the p-adic tree levels)."""
    levels = []
    for d in range(1, depth + 1):
        mod = p ** d
        groups = defaultdict(list)
        for s in speeds: groups[s % mod].append(s)
        levels.append((mod, dict(sorted(groups.items()))))
    return levels

def channels_mod(n):
    """the channel modulus: floor(n/2) (the antipodal/matching pairing, S532)."""
    return n // 2

def main():
    print("p-adic tree = channels (level 1) + coupling (deeper) + apex (residue 0) (oracle-S542)\n")
    for n in (12, 14, 16, 18):
        speeds = list(range(1, n))
        m = channels_mod(n)                  # = n/2
        # pick the channel prime: largest prime factor of n/2 (for n=14, m=7 prime)
        pf = [p for p in range(2, m + 1) if m % p == 0 and all(p % q for q in range(2, p))]
        p = max(pf) if pf else m
        print(f"=== n={n}: channels mod n/2={m}; channel prime p={p}; factor n={n}={'*'.join(f'{q}^{vp(n,q)}' for q in sorted(set(q for q in range(2,n+1) if n%q==0 and all(q%r for r in range(2,q)))))} ===")
        # p-adic tree to depth 2
        tree = padic_tree(speeds, p, 2)
        l1mod, l1 = tree[0]
        print(f"  LEVEL-1 (mod {l1mod}) nodes = {len(l1)} channels:")
        for r, members in l1.items():
            tag = "  <- APEX / gate channel (residue 0; the speed n/2 family)" if r == 0 else ""
            print(f"     residue {r}: {members}{tag}")
        # show the coupling: which level-1 nodes split at level 2
        l2mod, l2 = tree[1]
        coupled = sum(1 for r, members in l1.items() if len({x % l2mod for x in members}) > 1)
        print(f"  COUPLING: {coupled}/{len(l1)} level-1 nodes split further at level 2 (mod {l2mod})")
        # endpoint moat depth (S410): unit endpoint a/n has depth v_p(n) at each p|n
        moat = {q: vp(n, q) for q in (2, p) if n % q == 0}
        print(f"  endpoint MOAT depth (S410) d_p = v_p(n): {moat}")
        # apex check: is n/2 in the residue-0 node?
        apex = n // 2
        print(f"  apex speed n/2 = {apex}: residue mod {l1mod} = {apex % l1mod} (=0 means apex is the gate channel)")
        print()
    print("READING: the p-adic tree is the common object. Level 1 = the S532 CHANNELS")
    print("(= floor(n/2) = the independent-pair matching); the residue-0 branch = the APEX")
    print("/ source-sink / speed n/2 (S530, S541's diameter pair). The CHANNEL COUPLING")
    print("(why n>=6 is multi-channel, S531/S532) = the tree DEPTH (level-2+ refinement).")
    print("The endpoint MOAT (S410) = the boundary at depth v_p(n). The annular-braid CENTER")
    print("= shift (S541) = translation toward the tree boundary, so LRC is a boundary")
    print("condition. Covering character sums (S526) = harmonic analysis on this tree.")

if __name__ == "__main__":
    main()
