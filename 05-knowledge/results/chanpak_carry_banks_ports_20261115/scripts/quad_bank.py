# Abstract complete banks: quadratic functions f(y) = sum D_i y_i + sum_{i<j} s_ij y_i y_j on {0,1}^k that are
# bijections onto 2^k consecutive integers.  In a tournament every pair of head contacts interacts: s_ij != 0
# (|s_ij| = W(a,b) + W(b,a) >= H(B - a - b) >= 1).  Which k admit a bijection with ALL s_ij nonzero?
import itertools, sys
def search(k, require_nonzero=True, limit=5):
    N = 1 << k; sols = []; count = 0
    # f is determined by its values on vertices of weight <= 2; enumerate those values by backtracking,
    # then the rest is forced; values must be a permutation of {c..c+N-1} with f(0) = 0 (c in [-(N-1), 0]).
    low = [m for m in range(N) if bin(m).count('1') <= 2 and m]
    singles = [1 << i for i in range(k)]; pairs = [(1 << i) | (1 << j) for i in range(k) for j in range(i + 1, k)]
    def full(vals):
        D = [vals[1 << i] for i in range(k)]
        s = {(i, j): vals[(1 << i) | (1 << j)] - D[i] - D[j] for i in range(k) for j in range(i + 1, k)}
        f = []
        for m in range(N):
            v = sum(D[i] for i in range(k) if m >> i & 1) + sum(s[i, j] for i in range(k) for j in range(i + 1, k) if (m >> i & 1) and (m >> j & 1))
            f.append(v)
        return f, D, s
    for c in range(-(N - 1), 1):
        rng = [v for v in range(c, c + N) if v != 0]
        for sv in itertools.permutations(rng, k):
            if any(sv[i] > sv[i + 1] for i in range(k - 1)) and False: pass
            used = set(sv) | {0}
            rest = [v for v in rng if v not in used]
            for pv in itertools.permutations(rest, len(pairs)):
                vals = {0: 0}
                for m, v in zip(singles, sv): vals[m] = v
                for m, v in zip(pairs, pv): vals[m] = v
                f, D, s = full(vals)
                if require_nonzero and any(x == 0 for x in s.values()): continue
                if sorted(f) == list(range(c, c + N)):
                    count += 1
                    if len(sols) < limit: sols.append((D, s, f))
    return count, sols
for k in (2, 3):
    cnt, sols = search(k)
    print("k = %d: quadratic bijections onto 2^k consecutive integers with all pair interactions nonzero: %d" % (k, cnt))
    for D, s, f in sols[:4]: print("   D =", D, " s =", dict(s), " values", f)
cnt0, _ = search(3, require_nonzero=False)
print("k = 3 without the nonzero condition:", cnt0)
