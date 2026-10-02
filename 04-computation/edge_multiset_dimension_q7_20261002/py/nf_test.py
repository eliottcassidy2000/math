"""nf_test.py -- empirical test of the normal-form reduction (the WLOG step of the search).

For random k-subsets S of Q_d: apply the construction of the proof (coordinate j of minimal |beta_j| ->
transposed to d-1; flip d-1 if beta_{d-1} < 0; split S = A' x {0} u B' x {1}; h in Aut(Q_{d-1}) maps A'
to its lex-min sorted tuple, computed here by brute force over all 2^(d-1) (d-1)! elements, B = h(B')).
Checks: A is in the representative file, |beta_i(A x {0} u B x {1})| >= a - b for i < d-1, a >= b, and
the normal form lies in the same Aut(Q_d)-orbit as S (canonical forms via src/canon).
usage: python3 nf_test.py d k nsets seed
"""
import sys, os, random, itertools, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..'))

def main():
    d, k, nsets, seed = map(int, sys.argv[1:5])
    n = d - 1; N1 = 1 << n
    rng = random.Random(seed)
    group = []
    for pi in itertools.permutations(range(n)):
        for t in range(N1):
            group.append([sum(1 << pi[i] for i in range(n) if ((v ^ t) >> i) & 1) for v in range(N1)])
    reps_cache = {}
    def reps(a):
        if a not in reps_cache:
            R = np.fromfile(os.path.join(ROOT, 'data', 'q%d' % n, 'reps_n%d_a%d.bin' % (n, a)), dtype='<u8'); R.sort()
            reps_cache[a] = R
        return reps_cache[a]
    bad = 0; pairs = []; stats = {}
    for _ in range(nsets):
        S = rng.sample(range(1 << d), k)
        beta = [sum(1 if not (s >> i) & 1 else -1 for s in S) for i in range(d)]
        j = min(range(d), key=lambda i: abs(beta[i]))
        def sw(x):
            bj, bl = (x >> j) & 1, (x >> (d - 1)) & 1
            x &= ~((1 << j) | (1 << (d - 1))); return x | (bj << (d - 1)) | (bl << j)
        T = [sw(s) for s in S]
        if sum(1 if not (s >> (d - 1)) & 1 else -1 for s in T) < 0: T = [s ^ (1 << (d - 1)) for s in T]
        Ap = [s for s in T if s < N1]; Bp = [s - N1 for s in T if s >= N1]
        best = None
        for g in group:
            img = tuple(sorted(g[x] for x in Ap))
            if best is None or img < best[0]: best = (img, g)
        A = list(best[0]); g = best[1]; B = sorted(g[x] for x in Bp)
        a, b = len(A), len(B)
        mask = sum(1 << x for x in A)
        R = reps(a); pos = np.searchsorted(R, np.uint64(mask))
        inreps = pos < len(R) and int(R[pos]) == mask
        NF = A + [x + N1 for x in B]
        bnf = [sum(1 if not (s >> i) & 1 else -1 for s in NF) for i in range(d)]
        good = inreps and a >= b and bnf[d - 1] == a - b and all(abs(bnf[i]) >= a - b for i in range(d - 1)) and a - b == min(abs(x) for x in beta)
        bad += not good
        stats[a] = stats.get(a, 0) + 1
        pairs.append((S, NF))
    inp = '\n'.join(' '.join(map(str, X)) for p in pairs for X in p) + '\n'
    co = subprocess.run([os.path.join(ROOT, 'src', 'canon'), str(d)], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
    co = [l.split()[0] for l in co if l.strip()]
    same = all(co[2 * i] == co[2 * i + 1] for i in range(len(pairs)))
    print('NFTEST d=%d k=%d sets=%d domain_failures=%d same_orbit=%s a-distribution=%s' % (d, k, nsets, bad, same, dict(sorted(stats.items()))))

if __name__ == '__main__':
    main()
