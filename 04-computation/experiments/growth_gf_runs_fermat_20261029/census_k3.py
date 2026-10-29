"""Census of the ordered three-run observer over all tournaments of order n (gentourng classes)."""
import sys, itertools, json
from runs import *

def scores(A):
    return tuple(sorted(sum(r) for r in A))

def is_strong(A):
    n = len(A)
    reach = [set([i]) for i in range(n)]
    for _ in range(n):
        for i in range(n):
            for j in range(n):
                if A[i][j]:
                    reach[i] |= reach[j]
    return all(len(r) == n for r in reach)

def lifted_boundaries(n, tuples):
    """Run-symmetric lifts of the boundaries of all tetrahedra on the n vertices (as column vectors)."""
    index = {t: i for i, t in enumerate(tuples)}
    vecs = []
    for Q in itertools.combinations(range(n), 4):
        v = [0] * len(tuples)
        for pos, y in enumerate(Q):
            tri = tuple(x for x in Q if x != y)
            sign = -1 if pos % 2 else 1
            for perm in itertools.permutations(tri):
                v[index[perm]] += sign
        vecs.append(v)
    return vecs

def pdiv(summary, p):
    return sum(m for d, m in summary.items() if d and d % p == 0)

n = int(sys.argv[1])
rows = []
classes = tournaments(n)
cmds = []
data = []
for bits in classes:
    A = tour_from_bits(n, bits)
    tuples, M = run_matrix(A, 3)
    nb = len(tuples) // 6
    N, P, S, L0, L1 = unglued_k3(M, nb)
    B = lifted_boundaries(n, tuples)
    MB = [[sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))] for v in B]
    kills = all(all(x == 0 for x in col) for col in MB)
    data.append((bits, A, M, kills, len(B)))
    for X in (M, N, P, S, L0, L1):
        cmds.append('matsnf(%s)' % gp_matrix(X))
    # rank of the boundary span, and rank of [M ; plus-projection] etc.
    cmds.append('matrank(%s)' % gp_matrix(B) if B else '0')
    # dimension of the kernel of M intersected with run-symmetric vectors: rank of P gives it
out = gp_batch(cmds)
per = 7
print('n=%d classes=%d' % (n, len(classes)))
print('%-12s %-14s %-6s %-5s | %-26s | %-26s | %-10s %-10s %-22s | kills dB  rank dB  spur2 spur3' % ('bits','scores','strong','H','M3','unglued','plus','minus','std col'))
for c, (bits, A, M, kills, nB) in enumerate(data):
    snf = [snf_summary(parse_vec(out[per*c+i])) for i in range(6)]
    rB = out[per*c+6]
    H = sum(1 for w in itertools.permutations(range(n)) if all(A[w[i]][w[i+1]] for i in range(n-1)))
    sp2 = pdiv(snf[1], 2) - pdiv(snf[0], 2)
    sp3 = pdiv(snf[1], 3) - pdiv(snf[0], 3)
    print('%-12s %-14s %-6s %-5d | %-26s | %-26s | %-10s %-10s %-22s | %-5s %-7s %-5d %-5d' % (
        bits, ''.join(map(str, scores(A))), is_strong(A), H, fmt_snf(snf[0]), fmt_snf(snf[1]), fmt_snf(snf[2]),
        fmt_snf(snf[3]), fmt_snf(snf[4]) + (' =' if snf[4] == snf[5] else ' !=' + fmt_snf(snf[5])), kills, rB, sp2, sp3))
