#!/usr/bin/env python3
"""seams.py: seams around an even centre from one base.  A contact seam is a pair (T + v_X, T + v_{X u x}) whose counts
differ by 2 (centre c = their mean); an arc-flip seam is a pair (T, T with one arc reversed) with counts differing by 2.
(1) The arc-flip lemma H(T^{a<->b}) - H(T) = M(b, a) - M(a, b), checked by direct DP on orders <= 7, and the identity
sum over arcs a->b of M(a, b) = (n - 1) h.  (2) Contact seams from bases of orders <= 8 (values by THM-5266 (A)): the even
centres covered in [4, 2000], the misses; the digits |Delta|/2 of contact flips (which multipliers occur); the +2
window (3H+2, 3H+4) around 3(H+1) for odd counts H <= 300; how many seams change the quotient (number of strong
components).  (3) Arc-flip seams from orders <= 8: centres covered in [4, 2000]; flip digits."""
import numpy as np, itertools
from collections import Counter
exec(open('nearjoin_check.py').read().split('# (1)')[0])
P = np.load('profiles8.npz')
# (1) arc-flip lemma on orders 3-7
bad = tot = 0; ident_bad = 0
for n in range(3, 8):
    A, h, M = P[f'A{n}'], P[f'h{n}'], P[f'M{n}']
    for i in range(len(h)):
        T = A[i].tolist(); s = 0
        for a in range(n):
            for b in range(n):
                if T[a][b]:
                    s += M[i][a][b]; T2 = [r[:] for r in T]; T2[a][b] = 0; T2[b][a] = 1; tot += 1
                    if ham_all(T2)[1][-1] != h[i] + M[i][b][a] - M[i][a][b]: bad += 1
        if s != (n - 1) * h[i]: ident_bad += 1
print(f"(1) arc-flip lemma: {tot} flips on orders 3-7, {bad} mismatches; sum_(a->b) M(a,b) = (n-1)h fails {ident_bad} times")
# (2) contact seams, orders 1-8
LIM = 2000; centers = {}; digits = Counter(); seam_pairs = 0; quot_change = 0; window_hit = set(); counts_seen = set()
def ncomp(adj):
    n = len(adj); reach = []
    for v in range(n):
        r = 1 << v; f = [v]
        while f:
            a = f.pop()
            for b in range(n):
                if adj[a][b] and not (r >> b) & 1: r |= 1 << b; f.append(b)
        reach.append(r)
    comps = set()
    for v in range(n):
        comps.add(reach[v] & sum(1 << u for u in range(n) if (reach[u] >> v) & 1))
    return len(comps)
for n in range(1, 9):
    A, h, e, B, M = P[f'A{n}'], P[f'h{n}'], P[f'e{n}'], P[f'B{n}'], P[f'M{n}']
    k = len(h); N = 1 << n
    X = np.array([[(m >> x) & 1 for x in range(n)] for m in range(N)], dtype=np.int64)   # N x n
    lin = X @ (B - e).T                               # N x k
    quad = np.einsum('mi,kij,mj->mk', X, M, X)        # includes y == x terms: M diag is 0
    val = h[None, :] + lin - quad                     # N x k  values H(T + v_X)
    counts_seen.update(val.ravel().tolist())
    for i in range(k):
        v = val[:, i]
        for m in range(N):
            for x in range(n):
                if (m >> x) & 1: continue
                m2 = m | (1 << x); d = int(v[m2] - v[m]); digits[abs(d) // 2] += 1
                if abs(d) == 2:
                    seam_pairs += 1; c = int((v[m] + v[m2]) // 2)
                    if c <= LIM and c not in centers:
                        centers[c] = (n, i, m, x)
                        T = A[i].tolist()
                        U1 = [r[:] + [0 if (m >> j) & 1 else 1] for j, r in enumerate(T)] + [[1 if (m >> j) & 1 else 0 for j in range(n)] + [0]]
                        U2 = [r[:] + [0 if (m2 >> j) & 1 else 1] for j, r in enumerate(T)] + [[1 if (m2 >> j) & 1 else 0 for j in range(n)] + [0]]
                        if ncomp(U1) != ncomp(U2): quot_change += 1
    print(f"    order {n}: classes {k}, contact values up to {int(val.max())}, centres so far {len(centers)}")
missing = [c for c in range(4, LIM + 1, 2) if c not in centers]
print(f"(2) contact seams from bases of orders <= 8: {seam_pairs} seam pairs; even centres in [4, {LIM}] covered: {LIM//2 - 1 - len(missing)} of {LIM//2 - 1}; missing: {missing}")
print(f"    first-witness seams whose two tournaments have different numbers of strong components: {quot_change} of {len(centers)}")
print(f"    contact-flip digits |Delta|/2 occurring (orders <= 8): every value in [0, {max(d for d in range(0, 400) if all(digits[j] > 0 for j in range(d + 1)))}] occurs; largest {max(digits)}; share of flips with digit 1: {digits[1] / sum(digits.values()):.4f}, digit 0: {digits[0] / sum(digits.values()):.4f}")
odd_counts = sorted(c for c in counts_seen if c % 2 and c <= 300)
win = [H for H in odd_counts if (3 * H + 3) in centers]
print(f"(2w) the +2 window around 3(H+1): of the {len(odd_counts)} odd counts H <= 300 realised here, {len(win)} have (3H+2, 3H+4) as a contact seam from a base of order <= 8; failures: {[H for H in odd_counts if H not in win]}")
# (3) arc-flip seams
fcent = {}; fdig = Counter()
for n in range(3, 9):
    A, h, M = P[f'A{n}'], P[f'h{n}'], P[f'M{n}']
    for i in range(len(h)):
        T = A[i]
        for a in range(n):
            for b in range(n):
                if T[a][b]:
                    d = int(M[i][b][a] - M[i][a][b]); fdig[d // 2] += 1
                    if abs(d) == 2:
                        c = int(h[i]) + d // 2
                        if c <= LIM and c not in fcent: fcent[c] = (n, i, a, b)
fmissing = [c for c in range(4, LIM + 1, 2) if c not in fcent]
print(f"(3) arc-flip seams from orders <= 8: even centres in [4, {LIM}] covered: {LIM//2 - 1 - len(fmissing)}; missing: {fmissing[:60]}{' ...' if len(fmissing) > 60 else ''}")
print(f"    flip digits Delta/2: all values in [-{max(d for d in range(0, 400) if all(fdig[-j] > 0 for j in range(d + 1)))}, {max(d for d in range(0, 400) if all(fdig[j] > 0 for j in range(d + 1)))}] occur; share of count-preserving flips: {fdig[0] / sum(fdig.values()):.4f}; share with |digit| = 1: {(fdig[1] + fdig[-1]) / sum(fdig.values()):.4f}")
