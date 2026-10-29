import itertools, time
from runlib import *
# brute-force check of M_k against direct enumeration of path covers as sets of paths (n<=5)
def brute(A, k):
    n = len(A)
    tuples = ordered_tuples(n, k)
    idx = {t: i for i, t in enumerate(tuples)}
    M = [[0]*len(tuples) for _ in tuples]
    # enumerate ordered k-tuples of vertex-disjoint directed paths covering V
    paths = []
    for L in range(1, n+1):
        for seq in itertools.permutations(range(n), L):
            if all(A[seq[i]][seq[i+1]] for i in range(L-1)):
                paths.append(seq)
    def rec(i, used, chosen):
        if i == k:
            if used == (1<<n)-1:
                s = tuple(p[0] for p in chosen); t = tuple(p[-1] for p in chosen)
                M[idx[s]][idx[t]] += 1
            return
        for p in paths:
            m = 0
            for v in p: m |= 1<<v
            if m & used: continue
            rec(i+1, used|m, chosen+[p])
    rec(0, 0, [])
    return M
for n in (3,4,5):
    for bits in tournaments(n):
        A = adj_from_bits(n, bits)
        for k in range(1, n+1):
            if k == 1 and n > 5: continue
            _, M = run_matrix(A, k)
            assert M == brute(A, k), (n, bits, k)
print('subset-DP observer == brute-force path-cover enumeration for all classes n=3,4,5, all k')
A = circulant(5, {1,2})
t0=time.time()
tup, M = run_matrix(A, 3)
print('T5 M3 SNF', fmt(snf_summary(snf_diag(M))), 'mass', sum(map(sum,M)), 'rank', rank_Q(M), time.time()-t0)
N = unglued(M, 3)
print('T5 unglued N SNF', fmt(snf_summary(snf_diag(N))))
for name, co in (('plus',[0]),('minus',[1]),('std(2,4)',[2,4]),('std(3,5)',[3,5])):
    print(name, fmt(snf_summary(snf_diag(sector(N,3,co)))))
# block diagonal check
s=6; r=len(N)//s
blocks=[[0],[1],[2,4],[3,5]]
ok=True
for a in range(4):
    for b in range(4):
        if a==b: continue
        for o in range(r):
            for o2 in range(r):
                for c in blocks[a]:
                    for c2 in blocks[b]:
                        if N[o*s+c][o2*s+c2]: ok=False
print('N block diagonal in sectors (u),(v),(a,c),(b,d):', ok)
