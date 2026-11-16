#!/usr/bin/env python3
"""address_moves.py n [I Z g]: address-level small moves for the synchronized reader of variable_reader_coverage /
live_exterior_reader.  Exact count for an n-bit address b (then one final ordinary sink):
    H(b) = I + Z (2 D - rev_n(b)) + g (2 b - J(b)),   D = 2^n,   J(b) = sum_{i<j set} (2^(j-i-1) + [j = i+1]).
Defaults (I, Z, g) = (-11, 82, 9) (the first reader after the common synchronizing append).  We enumerate all 2^n
addresses, and for each gap in {2, 4, 6, 8} count the addresses b having a partner b' with H(b') - H(b) = gap, the
distribution of Hamming distances of the closest such partner, and whether the edit avoids the top k and bottom k bits
(k = n // 5), i.e. keeps both exterior phases.  Also reports the mod-9 residues reachable among addresses with fixed
top/bottom phases (the 'prescribe the count modulo 9' question)."""
import sys
from collections import Counter, defaultdict
n = int(sys.argv[1]); I, Z, g = (int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (-11, 82, 9)
D = 1 << n
def rev(b):
    r = 0
    for i in range(n):
        if b >> i & 1: r |= 1 << (n - 1 - i)
    return r
def J(b):
    bits = [i for i in range(n) if b >> i & 1]; s = 0
    for a in range(len(bits)):
        for c in range(a + 1, len(bits)):
            i, j = bits[a], bits[c]; s += (1 << (j - i - 1)) + (1 if j == i + 1 else 0)
    return s
H = [I + Z * (2 * D - rev(b)) + g * (2 * b - J(b)) for b in range(D)]
byval = defaultdict(list)
for b, h in enumerate(H): byval[h].append(b)
k = n // 5
print(f"n = {n}, (I, Z, g) = ({I}, {Z}, {g}); addresses {D}; distinct counts {len(byval)}; all counts odd: {all(h % 2 for h in H)}; range [{min(H)}, {max(H)}]")
for gap in (2, 4, 6, 8):
    have = 0; dist = Counter(); interior = 0; examples = []
    for b, h in enumerate(H):
        partners = byval.get(h + gap, [])
        if not partners: continue
        have += 1
        best = min(partners, key=lambda c: bin(b ^ c).count('1')); hd = bin(b ^ best).count('1'); dist[hd] += 1
        x = b ^ best; lowmask = (1 << k) - 1; highmask = ((1 << k) - 1) << (n - k)
        if not (x & lowmask) and not (x & highmask): interior += 1
        if len(examples) < 3 and hd <= 3: examples.append((b, best, bin(x)))
    print(f"gap +{gap}: {have} of {D} addresses have a partner ({100*have/D:.1f} %); Hamming distance of the closest partner: {sorted(dist.items())[:8]}; edits avoiding the top and bottom {k} bits: {interior}; examples (b, b', xor): {examples}")
# mod 9 within fixed phases: group addresses by (top k bits, bottom k bits) and report how many groups reach all 9 residues
groups = defaultdict(set)
for b, h in enumerate(H):
    groups[(b >> (n - k), b & ((1 << k) - 1))].add(h % 9)
full = sum(1 for s in groups.values() if len(s) == 9)
print(f"phase cells (top {k} bits, bottom {k} bits): {len(groups)}; cells whose addresses reach every residue mod 9: {full}; residue-count distribution: {sorted(Counter(len(s) for s in groups.values()).items())}")
