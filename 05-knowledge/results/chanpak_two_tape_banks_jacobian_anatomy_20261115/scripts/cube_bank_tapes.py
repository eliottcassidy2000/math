# Complete banks of head contacts inside a block C, with arbitrary tapes above and below (exact arithmetic).
# Digits S in C, block context X_C (contacts of the head in C - S), tape data (alpha, u1, beta, z1).
# In cut units the marginal of digit s (other digits off) is
#     Delta_s = m_s + kap_odd[s](beta + z1) + kap_even[s] beta - lam_odd[s](alpha + u1) - lam_even[s] alpha,
# and the pair interaction is omega_st = W(s,t) + W(t,s) (block-internal, independent of the tapes).
# A complete k-bank needs g(Y) = sum Delta_s y_s - sum omega_st y_s y_t to be a bijection of the k-cube onto 2^k
# consecutive integers.  For each admissible Delta the tape data solve a linear system; a rank-one system gives a
# whole line of tape data (the leading coefficient varies while the bank stays exact).
import sys, itertools
from fractions import Fraction as Fr
from master_formula_check import block_stats, brute
def read_classes(path, m):
    out = []
    for line in open(path):
        t = line.split()
        if len(t) < m: continue
        masks = [int(v) for v in t[:m]]
        out.append([[1 if masks[a] >> b & 1 else 0 for b in range(m)] for a in range(m)])
    return out
_cache = {}
def admissible(om):
    """all Delta tuples for which sum Delta y - sum om y y is a bijection onto 2^k consecutive integers."""
    key = tuple(om)
    if key in _cache: return _cache[key]
    k = {1: 2, 3: 3}[len(om)]; N = 1 << k; res = []
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for D in itertools.product(range(-(N - 1), N), repeat=k):
        if 0 in D: continue
        vals = []
        for Y in range(N):
            v = sum(D[i] for i in range(k) if Y >> i & 1) - sum(om[p] for p, (i, j) in enumerate(pairs) if (Y >> i & 1) and (Y >> j & 1))
            vals.append(v)
        if max(vals) - min(vals) == N - 1 and len(set(vals)) == N: res.append(D)
    _cache[key] = res; return res
def solve(rows, rhs):
    """integer solutions (alpha, beta) >= 0 of  -lam_s alpha + kap_s beta = rhs_s.  Returns ('point', (a, b)) / ('line', ...) / None."""
    k = len(rows)
    # find two independent rows
    for i in range(k):
        for j in range(i + 1, k):
            det = rows[i][0] * rows[j][1] - rows[i][1] * rows[j][0]
            if det != 0:
                a = Fr(rhs[i] * rows[j][1] - rows[i][1] * rhs[j], det); b = Fr(rows[i][0] * rhs[j] - rhs[i] * rows[j][0], det)
                if a.denominator != 1 or b.denominator != 1 or a < 0 or b < 0: return None
                if all(rows[t][0] * a + rows[t][1] * b == rhs[t] for t in range(k)): return ('point', (int(a), int(b)))
                return None
    # rank one: all rows proportional to rows[0]
    r0 = rows[0]
    for t in range(1, k):
        if rhs[t] * r0[0] != rhs[0] * rows[t][0]: return None
    # -lam alpha + kap beta = rhs0 : list the first solutions
    lam, kap = -r0[0], r0[1]; sols = []
    for a in range(0, 4 * (lam + kap) + abs(rhs[0]) + 2):
        num = rhs[0] + lam * a
        if num >= 0 and num % kap == 0: sols.append((a, num // kap))
        if len(sols) >= 3: break
    return ('line', sols) if sols else None
def search(m, path, k):
    found = []; nblocks = 0; pts = 0; lines = 0
    for C in read_classes(path, m):
        nblocks += 1; h, W, G, lam, kap = block_stats(C, m); full = (1 << m) - 1
        for S in itertools.combinations(range(m), k):
            prs = [(S[i], S[j]) for i in range(k) for j in range(i + 1, k)]
            om = [W[s][t][0] + W[t][s][0] for s, t in prs]
            if max(om) > 2 * ((1 << k) - 1): continue                       # bank rigidity
            adm = admissible(om)
            if not adm: continue
            rest = [v for v in range(m) if v not in S]
            for bits in range(1 << len(rest)):
                XC = sum(1 << rest[i] for i in range(len(rest)) if bits >> i & 1)
                ms = [sum(W[s][b][0] for b in range(m) if b != s and not XC >> b & 1) - sum(W[a][s][0] for a in range(m) if XC >> a & 1) for s in S]
                for u1 in (0, 1):
                    for z1 in (0, 1):
                        rows = [(-(lam[s][0] + lam[s][1]), kap[s][0] + kap[s][1]) for s in S]
                        for D in adm:
                            rhs = [D[i] - ms[i] - kap[S[i]][1] * z1 + lam[S[i]][1] * u1 for i in range(k)]
                            sol = solve(rows, rhs)
                            if sol:
                                if sol[0] == 'point': pts += 1
                                else: lines += 1
                                found.append((C, S, XC, u1, z1, D, tuple(om), sol, h[full]))
    return nblocks, pts, lines, found
if __name__ == '__main__':
    k = int(sys.argv[1]); Q = sys.argv[2]
    for m in [int(v) for v in sys.argv[3:]]:
        nb, pts, lines, found = search(m, "%s/classes_order%d.txt" % (Q, m), k)
        print("k=%d, blocks of order %d: %d classes; solutions: %d isolated tape data, %d lines" % (k, m, nb, pts, lines), flush=True)
        seen = set()
        for C, S, XC, u1, z1, D, om, sol, hC in found:
            key = (sol[0], hC, om, D)
            if key in seen: continue
            seen.add(key)
            if len(seen) <= 14 or sol[0] == 'line':
                print("   %s: H(C)=%d digits %s X_C=%d u1=%d z1=%d Delta=%s omega=%s tape data %s block %s" % (sol[0], hC, S, XC, u1, z1, D, om, sol[1], [sum(C[a][b] << b for b in range(m)) for a in range(m)]))
