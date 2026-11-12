# v_cells.py - independent recount of the collision cells and NET counts (verifier, lane energy)
# eta from q's recursion (different from the lane's direct eta recursion); cells from xor; exact int64.
import sys, time, math
import numpy as np
from numba import njit

@njit(cache=True)
def eta_words(n):
    N = 1 << n
    q = np.zeros(N, dtype=np.int64)
    for h in range(1, N):
        p = h >> 1
        if h & 1:
            q[h] = q[p] + (p + 1) // 2
        else:
            q[h] = q[p]
    lo = 1 << (n - 1)
    e = np.empty(N - lo, dtype=np.int64)
    for h in range(lo, N):
        e[h - lo] = 2 * h + 1 - 2 * q[h]
    return e

@njit(cache=True)
def rev_bits(h, n):
    r = 0
    for i in range(n):
        if (h >> i) & 1:
            r |= 1 << (n - 1 - i)
    return r

@njit(cache=True)
def cells(vals, idx, n, cell, R, stats):
    # stats[0] = pairs, [1] = violations of ratio restriction, [2] = run-free pairs, [3] = run-free pairs violating |s-t'|<=R+3
    # stats[4] = max |s-t'| over run-free pairs, [5] = pairs with X or Xi even (must be 0), [6] = pairs with |gamma|>2^w (0)
    m = vals.shape[0]
    i = 0
    while i < m:
        j = i + 1
        while j < m and vals[j] == vals[i]:
            j += 1
        for a in range(i, j):
            for b in range(a + 1, j):
                h = idx[a]; hp = idx[b]
                d = h ^ hp
                bl = 0; dd = d
                while dd > 0:
                    bl += 1; dd >>= 1
                tp = n - bl
                s = 0; dd = d
                while (dd & 1) == 0:
                    s += 1; dd >>= 1
                w = n - tp - s
                cell[tp, s] += 1
                stats[0] += 1
                x = (h >> s) & ((1 << w) - 1); xp = (hp >> s) & ((1 << w) - 1)
                X = x - xp
                Xi = rev_bits(x, w) - rev_bits(xp, w)
                if (X & 1) == 0 or (Xi & 1) == 0:
                    stats[5] += 1
                P = h >> (s + w)
                Z = h & ((1 << s) - 1)
                c = (1 << (s + 1)) - rev_bits(Z, s)
                g = c * X - P * Xi
                if abs(g) > (1 << w):
                    stats[6] += 1
                dcell = tp - s
                if dcell >= 4:
                    if not (abs(Xi) < (1 << (w - dcell + 3))):
                        stats[1] += 1
                # run-free test for both words
                rf = True
                for word in (h, hp):
                    run = 1; mx = 1
                    for k in range(1, n):
                        if ((word >> k) & 1) == ((word >> (k - 1)) & 1):
                            run += 1
                            if run > mx: mx = run
                        else:
                            run = 1
                    if mx > R: rf = False
                if rf:
                    stats[2] += 1
                    if abs(s - tp) > R + 3: stats[3] += 1
                    if abs(s - tp) > stats[4]: stats[4] = abs(s - tp)
        i = j

@njit(cache=True)
def net_count(prod_sorted, thr):
    # unordered pairs with |a-b| < thr (thr integer, products exact)
    m = prod_sorted.shape[0]
    cnt = 0
    lo = 0
    for hi in range(m):
        while prod_sorted[hi] - prod_sorted[lo] >= thr:
            lo += 1
        cnt += hi - lo
    return cnt

@njit(cache=True)
def products(n):
    lo = 1 << (n - 1)
    out = np.empty(lo, dtype=np.int64)
    for h in range(lo, 1 << n):
        out[h - lo] = h * ((1 << (n + 1)) - rev_bits(h, n))
    return out

def parse_lane(path):
    # returns {n: (distinct, pairs, {(tp,s): (pairs, dwp)})}
    res = {}
    cur = None
    mode = None
    for line in open(path):
        if line.startswith("n="):
            parts = dict(p.split("=") for p in line.split() if "=" in p)
            cur = int(parts["n"]); res[cur] = [int(parts["distinct"]), int(parts["unordered_pairs"]), {}]
            mode = None
        elif line.startswith("cells:"):
            mode = "cells"
        elif line.startswith("by "):
            mode = None
        elif mode == "cells" and cur is not None:
            f = line.split()
            if len(f) >= 5:
                res[cur][2][(int(f[0]), int(f[1]))] = (int(f[3]), int(f[4]))
    return res

def main():
    n0 = int(sys.argv[1]); n1 = int(sys.argv[2])
    lane = {}
    for p in ("../cells_16_20.out", "../cells_21_26.out"):
        try:
            lane.update(parse_lane(p))
        except Exception as ex:
            print("could not parse", p, ex)
    for n in range(n0, n1 + 1):
        t0 = time.time()
        e = eta_words(n)
        order = np.argsort(e, kind="stable")
        vals = e[order]
        idx = (order + (1 << (n - 1))).astype(np.int64)
        distinct = int(np.count_nonzero(np.diff(vals)) + 1)
        cell = np.zeros((n + 1, n + 1), dtype=np.int64)
        stats = np.zeros(8, dtype=np.int64)
        R = math.ceil(math.log2(n)) + 1
        cells(vals, idx, n, cell, R, stats)
        words = 1 << (n - 1)
        pairs = int(stats[0])
        # symmetry
        sym = True
        for tp in range(1, n):
            for s in range(0, n):
                if s + 1 <= n and tp - 1 >= 0 and cell[tp, s] != cell[s + 1, tp - 1]:
                    sym = False
        msg = ""
        if n in lane:
            ld, lp, lc = lane[n]
            same_cells = all(int(cell[tp, s]) == v[0] for (tp, s), v in lc.items()) and int(cell.sum()) == sum(v[0] for v in lc.values())
            msg = " | lane: distinct %s, pairs %s, all cell counts %s" % (ld == distinct, lp == pairs, same_cells)
        print("n=%d words=%d distinct=%d ratio=%.5f unordered_pairs=%d pairs/word=%.5f E_n/2^n=%.4f symmetry=%s%s" % (
            n, words, distinct, distinct / words, pairs, pairs / words, (words + 2 * pairs) / 2.0 ** n, sym, msg))
        print("   X,Xi odd violations=%d  |gamma|>2^w violations=%d  ratio-restriction violations (cells d>=4)=%d  run-free pairs=%d, violating |s-t'|<=R+3: %d (max |s-t'|=%d, R=%d)" % (
            stats[5], stats[6], stats[1], stats[2], stats[3], stats[4], R))
        # constants C_d: Coll(t', s)/2^(n-t') for t' = 2..6
        for d in (-3, -2, -1, 0, 1, 2, 3):
            row = []
            for tp in range(1, 8):
                s = tp + d
                if 0 <= s < n:
                    row.append("%.4f" % (cell[tp, s] / 2.0 ** (n - tp)))
                else:
                    row.append("  -   ")
            print("   d=%2d  Coll(t',t'+d)/2^(n-t') for t'=1..7: %s" % (d, " ".join(row)))
        # prefix sums
        cum = 0; rows = []
        for tp in range(n - 1, 0, -1):
            cum += int(cell[tp].sum()); rows.append((tp, cum))
        rows.reverse()
        print("   ordered pairs with common prefix >= t per word times 2^t, t=1..8: %s" % " ".join("%.3f" % (2 * c / words * 2.0 ** tp) for tp, c in rows[:8]))
        # NET
        pr = products(n); pr.sort()
        outs = []
        for T8 in (8, 32, 128, 3 * (n + 1), 5 * n, 6 * n):   # T = T8/8 : 1, 4, 16, 3(n+1)/8, 5n/8, 3n/4
            thr = T8 << (n - 3)                              # T * 2^n
            c = net_count(pr, thr)
            T = T8 / 8.0
            outs.append("T=%.3f: %.4f" % (T, c / (T * words)))
        print("   NET unordered pairs / (T * words): " + " | ".join(outs) + "   (%.1fs)" % (time.time() - t0))
        sys.stdout.flush()

if __name__ == "__main__":
    main()
