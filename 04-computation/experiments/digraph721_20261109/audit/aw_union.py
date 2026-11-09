"""AUDIT D1: top-band values realised by disjoint unions (refutes 'they need prime graphs', PROOF.md remark 7.2 /
README results 4).  For k = 2^(n-1) - j in the top band (M(n-1), 2^(n-1)) with 5 | k: P_3 u H, H the session's
(independently verified) even witness for k/5 on <= n - 3 vertices.  Writes aw_union_topband.txt for aw_verify.
usage: python3 aw_union.py ../witnesses.txt"""
import sys
def M(n): return 2**n - n if n % 2 else 2**n - n - 1
W = {}
for line in open(sys.argv[1]):
    left, _, right = line.partition('|'); t, s, n, k = left.split()[:4]
    if s == 'E': W[int(k)] = (int(n), [tuple(map(int, e.split('-'))) for e in right.split()])
P3 = [(0, 1), (1, 2)]
with open('aw_union_topband.txt', 'w') as out:
    for n, k in ((10, 505), (11, 1015), (12, 2045), (13, 4085)):
        assert M(n - 1) < k < 2**(n - 1) and k % 5 == 0
        m, H = W[k // 5]; assert m <= n - 3
        out.write(f"union-P3x{k // 5} E {n} {k} | " + " ".join(f"{a}-{c}" for a, c in P3 + [(a + 3, c + 3) for a, c in H]) + "\n")
        print(f"order {n}: top-band value {k} = 2^{n - 1} - {2**(n - 1) - k} = 5 * {k // 5} (realiser of {k // 5} on {m} vertices)")
