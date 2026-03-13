#!/usr/bin/env python3
"""
tradeoff_mechanism.py -- kind-pasteur-2026-03-13-S60

Investigate the N vs alpha_2 trade-off mechanism.

Key observation: Paley maximizes N but minimizes alpha_2.
Since H = 1 + 2N + 4*alpha_2 + 8*alpha_3, maximizing N HELPS
but minimizing alpha_2 HURTS. Why does N win?

Hypothesis: The cycle count N and disjoint pair count alpha_2 are
ANTI-correlated because more cycles means more overlaps.

Also: why exactly 4 distinct H values at p=11? What symmetry
explains the 2:10:10:10 distribution?
"""

from itertools import combinations
from collections import defaultdict

def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


p = 11
m = (p - 1) // 2

print(f"={'='*70}")
print(f"  TRADE-OFF MECHANISM AT p={p}")
print(f"={'='*70}")

# The 4 classes at p=11 have specific S-types
# Let's characterize by the number of QR elements in S
S_qr_set = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
print(f"  QR = {sorted(S_qr_set)}")
print(f"  NQR = {sorted(set(range(1,p)) - S_qr_set)}")

# For each pair {j, p-j}, one is QR and one is NQR (since -1 is NQR at p=11)
# The orientation choice is: which of {j, p-j} goes in S
# Paley picks all QR: S = QR
# Count how many QR elements each orientation picks

print(f"\n  Pairs: ", end="")
for j in range(1, m+1):
    qr_j = j if j in S_qr_set else p-j
    nqr_j = p-j if j in S_qr_set else j
    print(f"({qr_j},{nqr_j})", end=" ")
print()

class_map = {}
for bits in range(1 << m):
    S = []
    n_qr = 0
    for j in range(m):
        if bits & (1 << j):
            s = j + 1
        else:
            s = p - (j + 1)
        S.append(s)
        if s in S_qr_set:
            n_qr += 1
    class_map[bits] = n_qr

# Group by n_qr
qr_groups = defaultdict(list)
for bits, n_qr in class_map.items():
    qr_groups[n_qr].append(bits)

print(f"\n  Orientations grouped by #QR elements in S:")
for n_qr in sorted(qr_groups):
    bits_list = qr_groups[n_qr]
    print(f"    #QR={n_qr}: {len(bits_list)} orientations, bits={bits_list[:5]}{'...' if len(bits_list)>5 else ''}")

# From the data: H classes are {2, 29}, {4,5,7,11,15,16,20,24,26,27},
#                               {0,3,6,10,13,18,21,25,28,31}, {1,8,9,12,14,17,19,22,23,30}
# What's the #QR for each?
print(f"\n  H class membership by #QR:")
H_classes = {
    95095: [2, 29],
    93467: [4, 5, 7, 11, 15, 16, 20, 24, 26, 27],
    93027: [0, 3, 6, 10, 13, 18, 21, 25, 28, 31],
    92411: [1, 8, 9, 12, 14, 17, 19, 22, 23, 30]
}

for H_val, bits_list in sorted(H_classes.items(), reverse=True):
    qr_counts = [class_map[b] for b in bits_list]
    print(f"  H={H_val}: #QR values = {sorted(set(qr_counts))}, distribution = {sorted(qr_counts)}")

# KEY: It's not just #QR that matters. Let's look at which specific pairs are flipped.
# The 5 pairs at p=11 are {1,10}, {2,9}, {3,8}, {4,7}, {5,6}
# QR = {1,3,4,5,9}, NQR = {2,6,7,8,10}
# So pairs (QR, NQR) are: (1,10), (9,2), (3,8), (4,7), (5,6)
# bit j selects j+1 (vs p-(j+1)):
#   bit 0: 1 vs 10 → QR=1 is bit=1
#   bit 1: 2 vs 9 → QR=9 is bit=0
#   bit 2: 3 vs 8 → QR=3 is bit=1
#   bit 3: 4 vs 7 → QR=4 is bit=1
#   bit 4: 5 vs 6 → QR=5 is bit=1
# Paley: picks all QR = bits where we get QR from each pair
# bit0=1(pick 1=QR), bit1=0(pick 9=QR), bit2=1(pick 3=QR), bit3=1(pick 4=QR), bit4=1(pick 5=QR)
# = 0b11101 = 29 ✓

# The complement (all NQR): bits = 0b00010 = 2
# This should also give H=95095 since the complement tournament has same H (reversal)

print(f"\n  Paley = bits {29} = 0b{29:05b}")
print(f"  Anti-Paley = bits {2} = 0b{2:05b}")
print(f"  These are complements: bit pattern XOR = 0b{29^2:05b} = {29^2}")

# So the 4 H-classes correspond to some orbit structure of the bit patterns
# under the Z_5 symmetry? Or the QR/NQR structure?

# Let's see: what symmetry maps between orientations within the same H-class?
# Two tournaments T_S and T_S' are isomorphic iff S' = a*S mod p for some a coprime to p.
# For circulant tournaments on Z_p, the isomorphism group is generated by:
#   multiplication by elements of Z_p*
# Two orientation sets S and S' give isomorphic tournaments iff
# S' = a*S mod p for some a in Z_p*.

# BUT: we're restricting to |S| = m = (p-1)/2 (tournaments).
# The multiplier group acts on S.
# At p=11, Z_11* = {1,2,...,10}.
# QR = {1,3,4,5,9} = subgroup of index 2.
# Multiplying S_QR by a QR element gives S_QR back.
# Multiplying S_QR by a NQR element gives S_NQR = complement.

# For general S: aS mod p gives another valid orientation set.
# The orbit of S under Z_p* is the isomorphism class.

print(f"\n{'='*70}")
print(f"  ISOMORPHISM CLASSES via Z_p* action")
print(f"{'='*70}")

def orientation_to_set(bits, p, m):
    S = set()
    for j in range(m):
        if bits & (1 << j):
            S.add((j + 1) % p)
        else:
            S.add(p - (j + 1))
    return frozenset(S)

# Group orientations by isomorphism class
iso_classes = {}
for bits in range(1 << m):
    S = orientation_to_set(bits, p, m)
    # Find canonical form: min of all a*S mod p for a in Z_p*
    orbits = []
    for a in range(1, p):
        aS = frozenset((a * s) % p for s in S)
        orbits.append(aS)
    canon = min(orbits)
    if canon not in iso_classes:
        iso_classes[canon] = []
    iso_classes[canon].append(bits)

print(f"  Number of isomorphism classes: {len(iso_classes)}")
for canon, bits_list in sorted(iso_classes.items(), key=lambda x: min(x[1])):
    print(f"    {sorted(canon)}: bits = {bits_list}")

# Now deeper analysis: per-length cycle contribution to alpha_2
print(f"\n{'='*70}")
print(f"  PER-LENGTH DISJOINT PAIR ANALYSIS")
print(f"{'='*70}")

# For each H-class representative, break down alpha_2 by (k1, k2) pairs
for H_val, bits_list in sorted(H_classes.items(), reverse=True):
    bits = bits_list[0]
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A = build_adj(p, S)

    active_vsets = []
    for k in range(3, p+1, 2):
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            nc = count_ham_cycles(A, list(subset))
            if nc > 0:
                active_vsets.append((fs, k, nc))

    # alpha_2 by (k1, k2) type
    disj_by_type = defaultdict(int)
    for i in range(len(active_vsets)):
        for j in range(i+1, len(active_vsets)):
            V1, k1, n1 = active_vsets[i]
            V2, k2, n2 = active_vsets[j]
            if not (V1 & V2):
                key = tuple(sorted([k1, k2]))
                disj_by_type[key] += n1 * n2

    print(f"\n  H={H_val} (bits={bits}):")
    total_a2 = sum(disj_by_type.values())
    for kpair in sorted(disj_by_type):
        val = disj_by_type[kpair]
        print(f"    disj({kpair[0]},{kpair[1]}) = {val:>8} ({100*val/total_a2:>5.1f}%)")
    print(f"    TOTAL alpha_2 = {total_a2}")

# What about the c_k patterns?
print(f"\n{'='*70}")
print(f"  CYCLE COUNT PER LENGTH")
print(f"{'='*70}")

c_k_data = {}
for H_val, bits_list in sorted(H_classes.items(), reverse=True):
    bits = bits_list[0]
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A = build_adj(p, S)
    c_k = defaultdict(int)
    for k in range(3, p+1, 2):
        for subset in combinations(range(p), k):
            nc = count_ham_cycles(A, list(subset))
            c_k[k] += nc

    c_k_data[H_val] = dict(sorted(c_k.items()))

print(f"  {'H':>8} {'c3':>6} {'c5':>6} {'c7':>6} {'c9':>6} {'c11':>6} {'N':>8}")
for H_val in sorted(c_k_data, reverse=True):
    ck = c_k_data[H_val]
    N = sum(ck.values())
    print(f"  {H_val:>8} {ck.get(3,0):>6} {ck.get(5,0):>6} {ck.get(7,0):>6} {ck.get(9,0):>6} {ck.get(11,0):>6} {N:>8}")

# Differences from Paley
paley_ck = c_k_data[95095]
print(f"\n  Differences from Paley:")
print(f"  {'H':>8} {'dc3':>6} {'dc5':>6} {'dc7':>6} {'dc9':>6} {'dc11':>6} {'dN':>8}")
for H_val in sorted(c_k_data, reverse=True):
    if H_val == 95095:
        continue
    ck = c_k_data[H_val]
    diffs = {k: ck.get(k,0) - paley_ck.get(k,0) for k in [3,5,7,9,11]}
    dN = sum(diffs.values())
    print(f"  {H_val:>8} {diffs[3]:>+6} {diffs[5]:>+6} {diffs[7]:>+6} {diffs[9]:>+6} {diffs[11]:>+6} {dN:>+8}")

# The crucial question: are the differences all proportional?
# Or do different lengths vary independently?
print(f"\n  Ratios of differences (normalized by c5 change):")
for H_val in sorted(c_k_data, reverse=True):
    if H_val == 95095:
        continue
    ck = c_k_data[H_val]
    dc5 = ck.get(5,0) - paley_ck.get(5,0)
    if dc5 != 0:
        ratios = {}
        for k in [3,5,7,9,11]:
            dk = ck.get(k,0) - paley_ck.get(k,0)
            ratios[k] = dk / dc5
        print(f"  H={H_val}: dc/dc5 = {', '.join(f'c{k}:{ratios[k]:.3f}' for k in [3,5,7,9,11])}")

print("\nDONE.")
