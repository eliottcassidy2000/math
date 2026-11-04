"""Covering intervals of the root-response numerations from the triple census (data/triples_n.txt: h e B strong count).

Families over a rooted base (T, t) with triple (h, e, B):
  two-step  : h + 2B            (the C2 child of the lineage; the sole-exit tower at height 2; the near-join at k = 1)
  near-join : h + 2^k B, k <= 8 (HYP-9453, T -> M_k -> z with z -> t reversed)
  tower     : h_k = h A_k - e C_k + B D_k, 3 <= k <= 10 (PROOF.md)
Output: for each family the largest X with every odd number in [69, X] attained, the next gaps, and the odd numbers
below 69 that the family misses (7 and 21 are never attained, THM-1370).
usage: python3 coverage.py <max order> [data dir]   (orders 4..max; order 9 is regenerated with dedupe_triples.sh, 29 MB)"""
import sys, os
maxn = int(sys.argv[1]); D = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
def trib(a0, a1, a2, K):
    s = [a0, a1, a2]
    while len(s) <= K: s.append(s[-1] + s[-2] + s[-3])
    return s
K = 10; A, C, Dd = trib(1, 1, 1, K), trib(0, 1, 0, K), trib(0, 1, 2, K)
two = set(); two_s = set(); nj = set(); tw = set(); pairs = 0; triples = 0
for n in range(4, maxn + 1):
    for line in open(os.path.join(D, f'triples_{n}.txt')):
        h, e, B, st, cnt = map(int, line.split()); triples += 1; pairs += cnt
        v2 = h + 2 * B; two.add(v2)
        if st: two_s.add(v2)
        for k in range(1, 9): nj.add(h + (1 << k) * B)
        for k in range(3, K + 1): tw.add(h * A[k] - e * C[k] + B * Dd[k])
def interval(S, start=69):
    X = start
    while X + 2 in S or X + 2 in (7, 21): X += 2
    gaps = [v for v in range(X + 2, X + 402, 2) if v not in S]
    missed = [v for v in range(1, 69, 2) if v not in S and v not in (7, 21)]
    return X, gaps[:8], missed
print(f"bases of orders 4..{maxn}: {pairs} rooted tournaments, {triples} distinct (h, e, B, strong) rows")
for name, S in [("two-step h+2B, all bases", two), ("two-step, strong bases only", two_s), ("near-join h+2^k B (k<=8)", nj),
                ("tower h_k (3<=k<=10)", tw), ("two-step + near-join", two | nj), ("two-step + tower", two | tw), ("all three", two | nj | tw)]:
    X, gaps, missed = interval(S)
    print(f"{name:30s}: every odd in [69, {X}] attained; next gaps {gaps}; missed below 69: {missed}; |values| = {len(S)}")
