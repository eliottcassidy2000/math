# census.py - the complete rooted census of orders 1..9.
#
# A rooted tournament (T, v) of order m is (B, S): B = T - v of order m - 1 and S = N+(v).  So the rooted
# tournaments of order m are obtained from the isomorphism classes of order m - 1 and all 2^(m-1) cuts, and by
# F1, F2-type deletions and F4 the rooted state at the new vertex is
#     x = H(T) = sum_{b in S} s_b + sum_{a notin S} e_a + sum_{a notin S, b in S} Q_ab   (THM-4097 (4))
#     s(T, v) = sum_{b in S} s_b,            e(T, v) = sum_{a notin S} e_a,
#     Yplus(T, v) = H(B) + sum_{a notin S} Yplus(B, a),   Yminus(T, v) = H(B) + sum_{b in S} Yminus(B, b),
#     d = H(B).
# Output: rooted_le9.npz with the distinct rows (x, s, e, Yplus, Yminus, d, strong flag, out-degree of the root)
# of each order, reps_le8.npz with the
# class representatives, codes9.npy with a list of labelled tournaments of order 9 containing every class.
# Checks: class counts 1, 1, 2, 4, 12, 56, 456, 6880; the rows of every order m <= 8 agree as a set with the rows
# computed straight from the definitions at every root of every class of order m (and, m <= 6, of every labelled
# tournament); the same at order 9 on the list codes9; the H-values of orders 9 agree with the class census
# out/hist_9.txt of the mac-mini packet (1,520 values).
import os
import sys
import time
import numpy as np
from numba import njit
import core

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + ((" : " + detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def uniq(rows):
    return np.unique(rows, axis=0)


@njit(cache=True)
def filtered_cuts(n, code, keep):
    # cuts S of the new vertex w for which w has the least out-degree of B + w
    out = np.zeros(n, np.int64)
    core.out_masks(code, n, out)
    deg = np.zeros(n, np.int64)
    for v in range(n):
        m = out[v]
        c = 0
        while m:
            c += 1
            m &= m - 1
        deg[v] = c
    nk = 0
    for S in range(1 << n):
        dw = 0
        m = S
        while m:
            dw += 1
            m &= m - 1
        ok = True
        for v in range(n):
            dv = deg[v] + (0 if (S >> v) & 1 else 1)
            if dv < dw:
                ok = False
                break
        if ok:
            keep[nk] = S
            nk += 1
    return nk


def main():
    t0 = time.time()
    reps = {}
    cache = os.path.join(HERE, "reps_le8.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        reps = {n: z["n%d" % n] for n in range(1, 9)}
        print("class representatives loaded from reps_le8.npz")
    else:
        for n in range(1, 9):
            reps[n] = core.class_reps(n, core.perm_table(n))
            print("order %d: %d classes (%.1f s)" % (n, len(reps[n]), time.time() - t0), flush=True)
        np.savez_compressed(cache, **{"n%d" % n: reps[n] for n in range(1, 9)})
    check("class counts 1, 1, 2, 4, 12, 56, 456, 6880 (A000568)",
          [len(reps[n]) for n in range(1, 9)] == [1, 1, 2, 4, 12, 56, 456, 6880])

    rows = {1: np.array([[1, 1, 1, 1, 1, 1, 1, 0]], dtype=np.int64)}
    for n in range(1, 9):
        buf = np.zeros((len(reps[n]) << n, 8), np.int64)
        r = core.rooted_rows_from_bases(n, reps[n], buf)
        assert r == buf.shape[0]
        rows[n + 1] = uniq(buf)
        print("rooted order %d: %d (class, cut) pairs, %d distinct rows (x, s, e, Y+, Y-, d, strong, outdeg)  (%.1f s)"
              % (n + 1, r, len(rows[n + 1]), time.time() - t0), flush=True)

    # independent recomputation straight from the definitions, at every root of every class
    for m in range(1, 9):
        buf = np.zeros((len(reps[m]) * m, 10), np.int64)
        r = core.rooted_rows_direct(m, reps[m], buf)
        direct = uniq(buf[:r, :8])
        check("order %d: cut construction = direct definition at every root of every class (%d rows)"
              % (m, len(direct)), direct.shape == rows[m].shape and np.array_equal(direct, rows[m]))
    for m in range(1, 7):
        codes = np.arange(1 << (m * (m - 1) // 2), dtype=np.int64)
        buf = np.zeros((len(codes) * m, 10), np.int64)
        r = core.rooted_rows_direct(m, codes, buf)
        direct = uniq(buf[:r, :8])
        check("order %d: the same on all %d labelled tournaments" % (m, len(codes)),
              direct.shape == rows[m].shape and np.array_equal(direct, rows[m]))

    # a list of labelled tournaments of order 9 containing every class: B + w with w of least out-degree
    keep = np.zeros(256, np.int64)
    codes9 = []
    for code in reps[8]:
        nk = filtered_cuts(8, int(code), keep)
        for S in keep[:nk]:
            codes9.append(core.extend_code(int(code), 8, int(S)))
    codes9 = np.array(sorted(set(codes9)), dtype=np.int64)
    print("order 9: %d labelled tournaments in the list (new vertex of least out-degree)" % len(codes9), flush=True)
    np.save(os.path.join(HERE, "codes9.npy"), codes9)
    parts = []
    CH = 40000
    buf = np.zeros((CH * 9, 10), np.int64)
    for i in range(0, len(codes9), CH):
        r = core.rooted_rows_direct(9, codes9[i:i + CH], buf)
        parts.append(uniq(buf[:r, :8]))
    direct9 = uniq(np.concatenate(parts))
    check("order 9: cut construction over the 6,880 classes of order 8 = direct definition at every root of the "
          "list codes9 (%d rows)" % len(direct9),
          direct9.shape == rows[9].shape and np.array_equal(direct9, rows[9]))
    print("(%.1f s)" % (time.time() - t0), flush=True)

    # H-values of order 9 against the class census of the mac-mini packet
    hist = os.path.join(HERE, "..", "..", "..", "..", "05-knowledge", "results",
                        "tilted_ensemble_local_limit_20261110", "out", "hist_9.txt")
    if os.path.exists(hist):
        hv = sorted(int(l.split()[0]) for l in open(hist) if l.strip())
        mine = sorted(set(int(v) for v in rows[9][:, 0]))
        check("H-values of order 9: %d values, equal to the support of out/hist_9.txt" % len(mine), hv == mine)
    else:
        check("hist_9.txt of the mac-mini packet found", False)

    # known maxima of H: 1, 1, 3, 5, 15, 45, 189, 661, 3357
    mx = [int(rows[m][:, 0].max()) for m in range(1, 10)]
    check("max H at orders 1..9 = 1, 1, 3, 5, 15, 45, 189, 661, 3357", mx == [1, 1, 3, 5, 15, 45, 189, 661, 3357],
          str(mx))
    # converse symmetry: (x, s, e, Y+, Y-, d) -> (x, e, s, Y-, Y+, d) preserves the set
    for m in range(1, 10):
        sw = rows[m][:, [0, 2, 1, 4, 3, 5, 6, 7]].copy()
        sw[:, 7] = (m - 1) - sw[:, 7]
        sw = uniq(sw)
        check("order %d: the row set is invariant under converse (s <-> e, Y+ <-> Y-)" % m,
              np.array_equal(sw, rows[m]))

    np.savez_compressed(os.path.join(HERE, "rooted_le9.npz"), **{"m%d" % m: rows[m] for m in range(1, 10)})
    print("saved rooted_le9.npz (%.1f s)" % (time.time() - t0))
    if FAIL:
        print("FAILED:", FAIL)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
