#!/usr/bin/env python3
"""Audit A (independent), THM-4640 (5a)-(5d): two-position relations and carry derivability.

(i)  Brute force: all pairs of equal-length exponent sequences (length 2..5, exponents 0..EX) with equal
     val(B) = sum_j 3^(|B|-j) 2^(B_j) that differ in exactly one or two positions; classify the two-position ones.
(ii) 3^D (2^x - 1) = 2^y - 1: exhaustive in a box, by exact integer test.
(iii) Carry BFS without a state cap, for every canonical primitive pair with exponents <= E: from the padded short
     side, with j = 0..JMAX EXTRA leading 1s on both sides (the anchor relation (1).s ~ s allows any padding), and
     with exponent window [FLOOR, CAP]; report component sizes and whether the long side is reached.
Usage: python3 a5_carries.py [E=14] [CAP=30] [JMAX=3] [EX=12]
"""
import sys, time, itertools
from collections import deque, Counter, defaultdict

E = int(sys.argv[1]) if len(sys.argv) > 1 else 14
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 30
JMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 3
EX = int(sys.argv[4]) if len(sys.argv) > 4 else 12
t0 = time.time()


def val(B):
    L = len(B)
    return sum(3 ** (L - 1 - j) * 2 ** b for j, b in enumerate(B))


def K_of(s):
    K = -1
    for e in s:
        K = 3 * K + (1 << e)
    return K


# (i) brute force two-position classification
types = Counter(); one_pos = 0
for L in range(2, 6):
    byv = defaultdict(list)
    for B in itertools.product(range(0, EX + 1), repeat=L):
        byv[val(B)].append(B)
    for v, Bs in byv.items():
        for B1, B2 in itertools.combinations(Bs, 2):
            diff = [j for j in range(L) if B1[j] != B2[j]]
            if len(diff) == 1:
                one_pos += 1
            if len(diff) == 2:
                j, k = diff
                # normalise: the sequence with the larger entry at position j first
                if B1[j] < B2[j]:
                    B1, B2 = B2, B1
                a, b, c, d = B1[j], B2[j], B1[k], B2[k]
                Delta = k - j
                shape = (Delta, a - b, c - b, d - b)  # relative to b
                types[shape] += 1
print(f"(i) equal-length sequences, lengths 2..5, exponents 0..{EX}: one-position relations {one_pos}; two-position "
      f"relation shapes (Delta, a-b, c-b, d-b) with counts: {dict(types)}  [{time.time()-t0:.1f}s]", flush=True)

# (ii) 3^D (2^x - 1) = 2^y - 1
sols = []
for D in range(1, 201):
    p = 3 ** D
    for x in range(1, 1201):
        v = p * ((1 << x) - 1) + 1
        if v & (v - 1) == 0:
            sols.append((D, x, v.bit_length() - 1))
print(f"(ii) 3^D (2^x - 1) = 2^y - 1, D <= 200, x <= 1200: {sols}", flush=True)


# (iii) carry BFS
def nbrs(B, floor):
    out = []
    m = len(B)
    for j in range(m - 1):
        x, y = B[j], B[j + 1]
        if x == y + 1 and y + 2 <= CAP:                 # (b+1, b) -> (b, b+2)
            out.append(B[:j] + (y, y + 2) + B[j + 2:])
        if y == x + 2 and x >= floor:                   # (b, b+2) -> (b+1, b)
            out.append(B[:j] + (x + 1, x) + B[j + 2:])
    for j in range(m - 2):
        x, y = B[j], B[j + 2]
        if x == y + 3 and y + 6 <= CAP:                 # (b+3, *, b) -> (b, *, b+6)
            out.append(B[:j] + (y,) + B[j + 1:j + 2] + (y + 6,) + B[j + 3:])
        if y == x + 6 and x >= floor:                   # (b, *, b+6) -> (b+3, *, b)
            out.append(B[:j] + (x + 3,) + B[j + 1:j + 2] + (x,) + B[j + 3:])
    return out


print("    neighbours of (1,1,7):", nbrs((1, 1, 7), 0), " val check:", val((1, 1, 7)), val((4, 1, 1)), flush=True)

reps = defaultdict(list)
for n in range(E + 1):
    for s in itertools.combinations(range(2, E + 1), n):
        reps[K_of(s)].append(s)
pairs = []
for k, rs in reps.items():
    if k == -1:
        continue
    for s, t in itertools.combinations(rs, 2):
        if s[-1] != t[-1]:
            if len(s) > len(t):
                s, t = t, s
            pairs.append((k, s, t))
print(f"(iii) {len(pairs)} canonical primitive pairs with exponents <= {E}", flush=True)
for floor in (0, -6):
    for j in range(0, JMAX + 1):
        found = 0; sizes = []
        for k, s, t in pairs:
            L = len(t) + j
            src = (1,) * (L - len(s)) + s
            tgt = (1,) * j + t
            seen = {src}; dq = deque([src]); hit = False
            while dq:
                B = dq.popleft()
                for C in nbrs(B, floor):
                    if C not in seen:
                        seen.add(C); dq.append(C)
            hit = tgt in seen
            found += hit; sizes.append(len(seen))
        print(f"    floor {floor:2d}, extra padding j = {j}: derivable {found}/{len(pairs)}; component size max "
              f"{max(sizes)}, mean {sum(sizes)/len(sizes):.2f}, size-1 components {sizes.count(1)}  "
              f"[{time.time()-t0:.1f}s]", flush=True)
