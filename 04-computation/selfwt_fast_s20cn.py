#!/usr/bin/env python3
"""
selfwt_fast_s20cn.py -- kind-pasteur-2026-03-22-S20cn

FAST verification of self_wt_7 via EXHAUSTIVE enumeration.

Key optimization: for each of 2^21 tournaments, compute (score, c3).
Then for each arc flip, check if the flipped tournament has the SAME
(score, c3). If YES, they're likely isomorphic (self-flip).
If NO, definitely not isomorphic.

At n<=6: (score, c3) is a COMPLETE invariant (determines iso class).
At n=7: (score, c3) might have rare collisions, giving slight OVERestimates.

We also compute EXACT twin counts (much faster than iso check).

Predicted:
  total_twins = m * 2^{m-n+2} = 21 * 2^16 = 1,376,256
  NTE(7) = 151,200
  self_wt_7 = 1,527,456

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys
from math import comb
import time
sys.stdout.reconfigure(line_buffering=True)

n = 7
m = comb(n, 2)  # 21
total = 1 << m   # 2097152

PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

# Precompute: for each arc k and vertex a, how does flipping k affect score of a?
# Flipping arc k = {i,j} changes:
#   score[i] by +1 or -1 (if i was beating j, now losing, so -1; vice versa)
#   score[j] by -1 or +1 (opposite of i)
# All other scores unchanged.

print("=" * 70)
print("  FAST self_wt COMPUTATION AT n=7")
print("=" * 70)

t0 = time.time()

# EXACT twin pair enumeration
# A twin pair {i,j} means T[a][i] = T[a][j] for all a != i,j.
# We check this for ALL 2^21 tournaments.

total_twins = 0
total_self_fp = 0  # self-flips detected by (score,c3) fingerprint
total_self_fp_nontwin = 0

# For efficiency, process in batches
# For each tournament (represented as bits):
#   - Compute score sequence
#   - Compute c3
#   - For each arc k, compute score/c3 of flipped tournament
#   - If match: possible self-flip
#   - Also check twin condition for arc k

def compute_scores(bits):
    """Compute sorted score sequence."""
    s = [0] * n
    for k, (i,j) in enumerate(PAIRS):
        if bits & (1 << k):
            s[i] += 1
        else:
            s[j] += 1
    return s  # unsorted, to track individual scores

def compute_c3(bits):
    """Count directed 3-cycles."""
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations of the 3-cycle
                ij = 1 if bits & (1 << PAIRS.index((i,j))) else 0
                jk = 1 if bits & (1 << PAIRS.index((j,k))) else 0
                ik = 1 if bits & (1 << PAIRS.index((i,k))) else 0
                # i->j->k->i: ij=1, jk=1, ki=1 (ik=0)
                if ij and jk and not ik:
                    c3 += 1
                # i->k->j->i: ik=1, kj=1 (jk=0), ji=1 (ij=0)
                if not ij and not jk and ik:
                    c3 += 1
    return c3

# Precompute pair index lookup
pair_idx = {}
for k, (i,j) in enumerate(PAIRS):
    pair_idx[(i,j)] = k

# For c3 delta: when we flip arc k={p,q}, the change in c3 depends on
# the triangles containing {p,q}.
# Triangle {p,q,r}: contributes to c3 based on the 3 arc directions.
# Flipping {p,q} changes only the p-q direction, affecting all triangles through {p,q}.
# Number of such triangles: n-2 = 5.

# For each arc k, precompute the triangles it participates in
arc_triangles = {}
for k, (p,q) in enumerate(PAIRS):
    triangles = []
    for r in range(n):
        if r == p or r == q: continue
        # Triangle {p,q,r}
        # Arcs: (p,q), (p,r), (q,r)
        # We need the bit positions of (p,r) and (q,r)
        pr = pair_idx[(min(p,r), max(p,r))]
        qr = pair_idx[(min(q,r), max(q,r))]
        triangles.append((r, pr, qr))
    arc_triangles[k] = triangles

def c3_delta(bits, k):
    """Compute change in c3 when flipping arc k."""
    p, q = PAIRS[k]
    pq_bit = 1 if bits & (1 << k) else 0  # 1 if p->q, 0 if q->p

    delta = 0
    for r, pr_idx, qr_idx in arc_triangles[k]:
        # pr_idx stores the pair index for {p,r}. The bit is 1 if min(p,r)->max(p,r).
        # We need: pr = 1 iff p->r.
        raw_pr = 1 if bits & (1 << pr_idx) else 0
        pr_bit = raw_pr if p < r else 1 - raw_pr

        raw_qr = 1 if bits & (1 << qr_idx) else 0
        qr_bit = raw_qr if q < r else 1 - raw_qr

        # Before flip: count 3-cycle contribution of {p,q,r}
        # Directions: p->q if pq_bit, p->r if pr_bit, q->r if qr_bit
        # After flip: pq_bit flips
        old_c3 = is_3cycle(pq_bit, pr_bit, qr_bit, p, q, r)
        new_c3 = is_3cycle(1 - pq_bit, pr_bit, qr_bit, p, q, r)
        delta += new_c3 - old_c3

    return delta

def is_3cycle(pq, pr, qr, p, q, r):
    """Check if {p,q,r} forms a directed 3-cycle given arc directions."""
    # p->q: pq=1, q->p: pq=0
    # p->r: pr=1, r->p: pr=0
    # q->r: qr=1, r->q: qr=0
    # 3-cycle p->q->r->p: pq=1, qr=1, pr=0 (r->p)
    # 3-cycle p->r->q->p: pq=0, qr=0, pr=1
    if pq and qr and not pr:
        return 1
    if not pq and not qr and pr:
        return 1
    return 0

print(f"  Enumerating all {total} tournaments...")

for bits in range(total):
    scores = compute_scores(bits)
    sorted_score = tuple(sorted(scores))
    c3 = compute_c3(bits)

    for k in range(m):
        i, j = PAIRS[k]

        # ALWAYS check twin (independent of fingerprint)
        is_twin = all(
            ((bits >> pair_idx[(min(a,i), max(a,i))]) & 1 if a < i else
             1 - ((bits >> pair_idx[(min(i,a), max(i,a))]) & 1))
            ==
            ((bits >> pair_idx[(min(a,j), max(a,j))]) & 1 if a < j else
             1 - ((bits >> pair_idx[(min(j,a), max(j,a))]) & 1))
            for a in range(n) if a != i and a != j
        )

        if is_twin:
            total_twins += 1

        # Score delta: flipping {i,j}
        if bits & (1 << k):
            new_scores_i = scores[i] - 1
            new_scores_j = scores[j] + 1
        else:
            new_scores_i = scores[i] + 1
            new_scores_j = scores[j] - 1

        new_scores = list(scores)
        new_scores[i] = new_scores_i
        new_scores[j] = new_scores_j
        new_sorted = tuple(sorted(new_scores))

        if new_sorted != sorted_score:
            continue

        dc3 = c3_delta(bits, k)
        if dc3 != 0:
            continue

        # Same (score, c3) => likely self-flip
        total_self_fp += 1

        if not is_twin:
            total_self_fp_nontwin += 1

    if (bits + 1) % 200000 == 0:
        elapsed = time.time() - t0
        rate = (bits + 1) / elapsed
        eta = (total - bits - 1) / rate
        print(f"    {bits+1}/{total} ({rate:.0f}/s, ETA {eta:.0f}s) twins={total_twins} fp_self={total_self_fp}")

elapsed = time.time() - t0
print(f"\n  Done in {elapsed:.1f}s")

print(f"\n  RESULTS:")
print(f"  Total twin-pair self-flips: {total_twins}")
print(f"  Predicted twins: {m * (1 << (m-n+2))}")
print(f"  Match: {total_twins == m * (1 << (m-n+2))}")

print(f"\n  Total (score,c3)-matched self-flips: {total_self_fp}")
print(f"  Predicted self_wt: {m * (1 << (m-n+2)) + 21 * 7200}")
print(f"  Non-twin (score,c3) matches: {total_self_fp_nontwin}")
print(f"  Predicted NTE(7): 151200")

if total_self_fp > 0:
    nte_measured = total_self_fp - total_twins
    print(f"\n  NTE measured (approx): {nte_measured}")
    print(f"  NTE predicted: 151200")
    print(f"  NTE match: {nte_measured == 151200}")
    if nte_measured != 151200:
        print(f"  Difference: {nte_measured - 151200}")
        print(f"  (Positive = false positives from score/c3 collisions)")
