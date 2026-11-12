"""d2_independent_direct_closure.py - INDEPENDENT closure of the direct joins of one-apex blocks (addendum D;
audit of scripts/direct_joins.py).

No normal form: every binary word w of length 0..NMAX is a block with its own data
    g (closed formula of b3, checked there against the subset DP),  n = length,
    P = value of the letters 1 (first spine vertex = top bit),  Q = value of the letters 0 (first = unit bit),
    O = 1 + Q//2, E = (Q+1)//2,   odd = 1 + P//2, evn = (P+1)//2,   t_e = t_o = 2^(n-1) (n >= 1), (1, 0) (n = 0).
Law D1 (checked by brute force in b2):
    A' = A g + 2 (o E + e O),   o' = A odd + o t_o + e t_e,   e' = A evn + o t_e + e t_o.
First block: (g, odd, evn) for every word.

Word length.  A successor with A' <= X needs O <= X/2 (if e >= 1) or E <= X/2 (if o >= 1; o >= 1 always), hence
Q <= X < 2^(K+1): all letters 0 are at positions <= K, i.e. w = w' 1^m with |w'| <= K + 1, and A' depends on w'
only.  NMAX = K + 2 is used.

Modes
  literal K   every state with A <= X is stored AND expanded (no pruning at all).  Run for K = 8 and K = 11
              (242,785 states at K = 11).
  numpy K     a state is expanded iff A + 2 e <= X.  This is exact: every successor has
              A'' = A g + 2 (o E + e O) >= A + 2 e  (g >= 1, O >= 1), with equality for the empty word.
  numba K     same rule; words sorted by O and scanned while A + 2 e O <= X (exact cut); own hash set (numba set).
  blocks K L1 L2 ...   bonus (addendum E, not in the assigned scope): the numba closure with every block restricted
              to words of length <= L.  Output: e1_finite_blocks_bonus.out.
The three modes are compared with one another where they overlap, with the stored runs, and with the author's
closure called directly (full list of misses).  Output: d2_independent_direct_closure.out.
"""
import os
import sys
import time
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from b3_independent_closure import g_closed          # noqa: E402


def word_arrays(nmax):
    g, O, E, od, ev, te, to = [], [], [], [], [], [], []
    for n in range(0, nmax + 1):
        for w in product((0, 1), repeat=n):
            P = 0
            for a in w:
                P = 2 * P + a
            Q = sum(1 << j for j, a in enumerate(w) if not a)
            g.append(g_closed(w))
            O.append(1 + Q // 2)
            E.append((Q + 1) // 2)
            od.append(1 + P // 2)
            ev.append((P + 1) // 2)
            if n >= 1:
                te.append(1 << (n - 1))
                to.append(1 << (n - 1))
            else:
                te.append(1)
                to.append(0)
    arr = lambda v: np.array(v, dtype=np.int64)
    return arr(g), arr(O), arr(E), arr(od), arr(ev), arr(te), arr(to)


def closure_numpy(K, literal, nmax=None):
    X = 1 << K
    nmax = K + 2 if nmax is None else nmax
    g, O, E, od, ev, te, to = word_arrays(nmax)
    hit = np.zeros(X + 1, dtype=bool)
    seen = set()
    stack = []
    sel = g <= X
    hit[g[sel]] = True
    for st in set(zip(g[sel].tolist(), od[sel].tolist(), ev[sel].tolist())):
        if literal or st[0] + 2 * st[2] <= X:
            seen.add(st)
            stack.append(st)
    expanded = 0
    while stack:
        A, o, e = stack.pop()
        expanded += 1
        vals = g * A + 2 * (o * E + e * O)
        ok = vals <= X
        if not ok.any():
            continue
        A2 = vals[ok]
        o2 = A * od[ok] + o * to[ok] + e * te[ok]
        e2 = A * ev[ok] + o * te[ok] + e * to[ok]
        hit[A2] = True
        if literal:
            cand = set(zip(A2.tolist(), o2.tolist(), e2.tolist()))
        else:
            live = (A2 + 2 * e2) <= X
            cand = set(zip(A2[live].tolist(), o2[live].tolist(), e2[live].tolist()))
        for st in cand:
            if st not in seen:
                seen.add(st)
                stack.append(st)
    odd = np.arange(1, X + 1, 2)
    return odd[~hit[odd]].tolist(), expanded, len(seen)


def closure_numba(K, nmax=None):
    """nmax: optional cap on the word length of every block (addendum E: finitely many blocks)"""
    from numba import njit

    @njit(cache=False)
    def run(X, g, O, E, od, ev, te, to, hit):
        seen = set()
        stack = []
        nw = g.shape[0]
        for c in range(nw):
            if g[c] <= X:
                hit[g[c]] = 1
                if g[c] + 2 * ev[c] <= X:
                    key = (g[c] << 44) | (od[c] << 22) | ev[c]
                    if key not in seen:
                        seen.add(key)
                        stack.append(key)
        expanded = 0
        mask = (1 << 22) - 1
        while len(stack) > 0:
            key = stack.pop()
            A = key >> 44
            o = (key >> 22) & mask
            e = key & mask
            expanded += 1
            for c in range(nw):
                if e >= 1 and A + 2 * e * O[c] > X:
                    break                               # words are sorted by O
                A2 = A * g[c] + 2 * (o * E[c] + e * O[c])
                if A2 > X:
                    continue
                hit[A2] = 1
                e2 = A * ev[c] + o * te[c] + e * to[c]
                if A2 + 2 * e2 > X:
                    continue
                o2 = A * od[c] + o * to[c] + e * te[c]
                if o2 > mask:
                    continue                            # cannot happen for a live state: o2 <= e2 + A2 (see note)
                k2 = (A2 << 44) | (o2 << 22) | e2
                if k2 not in seen:
                    seen.add(k2)
                    stack.append(k2)
        return expanded, len(seen)

    X = 1 << K
    g, O, E, od, ev, te, to = word_arrays(K + 2 if nmax is None else nmax)
    order = np.argsort(O, kind="stable")
    g, O, E, od, ev, te, to = g[order], O[order], E[order], od[order], ev[order], te[order], to[order]
    hit = np.zeros(X + 1, dtype=np.uint8)
    expanded, nseen = run(X, g, O, E, od, ev, te, to, hit)
    odd = np.arange(1, X + 1, 2)
    return odd[hit[odd] == 0].tolist(), int(expanded), int(nseen)


def author_closure(K):
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "scripts")))
    import direct_joins as dj
    from chain_semigroup import g_list, L_of, R_of
    X = 1 << K
    half = X // 2
    G = g_list(4 * X + 8)
    cores = [0] + [p for p in range(2, 2 * X + 6, 2)]
    coreP = np.array(cores, dtype=np.int64)
    coreG = np.array([G[p] for p in cores], dtype=np.int64)
    coreQ = np.array([R_of(p) for p in cores], dtype=np.int64)
    coreN = np.array([L_of(p) for p in cores], dtype=np.int64)
    hit = np.zeros(X + 1, dtype=np.uint8)
    tbits = 25
    table = np.full(1 << tbits, -1, dtype=np.int64)
    cap = 1 << (tbits - 1)
    stK = np.zeros(cap, dtype=np.int64)
    ns = 0
    for P in range(0, 4 * X + 8):
        A, o, e = G[P], 1 + P // 2, (P + 1) // 2
        if A <= X:
            hit[A] = 1
            if o <= half and e <= half:
                key = (A << 40) | (o << 20) | e
                if dj.insert(table, key):
                    stK[ns] = key
                    ns += 1
    states = dj.closure(X, coreP, coreG, coreQ, coreN, hit, table, stK, ns)
    odd = np.arange(1, X + 1, 2)
    return odd[hit[odd] == 0].tolist(), int(states)


def stored(K):
    path = os.path.join(HERE, "..", "runs", "direct_closure_%d.txt" % K)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        lines = f.read().splitlines()
    count = int(lines[0].split("odd numbers missed:")[1])
    lst = [int(s) for s in lines[1].split("[")[1].split("]")[0].split(",")]
    return count, lst


def main():
    mode = sys.argv[1]
    K = int(sys.argv[2])
    t0 = time.time()
    if mode == "blocks":
        # bonus check of addendum E (not in the audit's assigned scope): words of length <= L only
        for L in [int(a) for a in sys.argv[3:]]:
            t0 = time.time()
            miss, expanded, nseen = closure_numba(K, L)
            X = 1 << K
            hist = [sum(1 for m in miss if (1 << k) <= m < (1 << (k + 1))) for k in range(8, K)]
            print("words of length <= %d (%d blocks), closure below 2^%d: %d states expanded; odd numbers missed: %d; largest miss %d; misses in [2^k, 2^(k+1)), k = 8..%d: %s  [%.1f s]"
                  % (L, (1 << (L + 1)) - 1, K, expanded, len(miss), miss[-1] if miss else 0, K - 1, hist, time.time() - t0))
        return
    if mode == "literal":
        miss, expanded, nseen = closure_numpy(K, True)
        label = "literal (every state with A <= X stored and expanded)"
    elif mode == "numpy":
        miss, expanded, nseen = closure_numpy(K, False)
        label = "numpy (expanded iff A + 2e <= X)"
    else:
        miss, expanded, nseen = closure_numba(K)
        label = "numba (expanded iff A + 2e <= X; words sorted by O, exact cut)"
    print("independent closure of direct joins below 2^%d, mode %s: %d states expanded, %d stored; odd numbers missed: %d  [%.1f s]"
          % (K, label, expanded, nseen, len(miss), time.time() - t0))
    print("   misses:", miss)
    if len(sys.argv) > 3 and sys.argv[3] == "noauthor":
        ref = None
    else:
        ref, states = author_closure(K)
        print("author's closure below 2^%d: %d live states, %d misses; full lists identical: %s" % (K, states, len(ref), ref == miss))
    st = stored(K)
    if st is not None:
        print("stored runs/direct_closure_%d.txt: count %d, list identical to mine: %s" % (K, st[0], st[1] == miss and st[0] == len(miss)))


if __name__ == "__main__":
    main()
