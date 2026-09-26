#!/usr/bin/env python3
"""Orchestrator audit of lane `drift` (entropy/merge law for sign strategies),
written from the note's statements; the lane's scripts were not read.

For random and structured sign strategies of the q n +- 1 cube (q = 3, 5, 7, 9, 23;
k = 3..9) compute every closed class of the uniform-lift chain, its stationary
law pi, and check
  Theorem 1:  pi(odd) = 1/2 - sum_{s in R} pi(s) g(s)   (g = odd steps removed
              from the next k-1 forced steps by the flip at s)
  Theorem 2:  1 - h(pi(odd)) <= sum_merge (pi(s)+pi(s*)) h(pi(s)/(pi(s)+pi(s*)))
              <= pi(R u R*) <= pi(odd),   s* = s - 2 q^-1 mod 2^k
  Cor 2.1:    p0 = root of p + h(p) = 1 is 0.2270922, and log_23 2 < p0 < log_21 2;
              exhaustively no class-(i) strategy for q = 23 at k <= 4 (all 2^(2^(k-1))).
"""
import math, random
import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def h(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def parity_word(x, q, m):
    w = []
    for _ in range(m):
        w.append(x & 1)
        x = (q * x + 1) // 2 if x & 1 else x // 2
    return w


def targets(k, q, sig):
    N, H = 1 << k, 1 << (k - 1)
    t = []
    for s in range(N):
        if s % 2 == 0:
            tt = (s // 2) % H
        else:
            tt = ((q * s + sig[s // 2]) // 2) % H
        t.append(tt)
    return t


def closed_class_laws(k, q, sig):
    N, H = 1 << k, 1 << (k - 1)
    t = targets(k, q, sig)
    rows, cols = [], []
    for s in range(N):
        rows += [s, s]
        cols += [t[s], t[s] + H]
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N))
    ncomp, lab = connected_components(A, directed=True, connection="strong")
    laws = []
    for c in range(ncomp):
        nodes = np.nonzero(lab == c)[0]
        S = set(nodes.tolist())
        if any(t[s] not in S or t[s] + H not in S for s in nodes):
            continue  # not closed
        idx = {s: i for i, s in enumerate(nodes)}
        m = len(nodes)
        P = np.zeros((m, m))
        for s in nodes:
            P[idx[s], idx[t[s]]] += 0.5
            P[idx[s], idx[t[s] + H]] += 0.5
        # solve pi (P - I) = 0, sum pi = 1
        M = np.vstack([(P - np.eye(m)).T, np.ones(m)])
        b = np.zeros(m + 1); b[-1] = 1
        pi, *_ = np.linalg.lstsq(M, b, rcond=None)
        assert np.abs(pi @ P - pi).max() < 1e-10 and pi.min() > -1e-12
        full = np.zeros(N)
        full[nodes] = pi
        laws.append(full)
    return laws


def audit(k, q, sig):
    N, H = 1 << k, 1 << (k - 1)
    R = [2 * i + 1 for i in range(H) if sig[i] == -1]
    qi = pow(q, -1, N)
    g = {}
    for s in R:
        wplus = parity_word(s, q, k)[1:k]
        ym = ((q * s - 1) // 2) % H
        u = parity_word(ym, q, k - 1)
        g[s] = sum(wplus) - sum(u)
    worst = 0.0
    for pi in closed_class_laws(k, q, sig):
        p = pi[1::2].sum()
        # Theorem 1
        rhs = 0.5 - sum(pi[s] * g[s] for s in R)
        assert abs(p - rhs) < 1e-9, (k, q, p, rhs)
        # Theorem 2
        Rset = set(R)
        merge = 0.0
        RuRs = set(R)
        for s in R:
            ss = (s - 2 * qi) % N
            assert ss % 2 == 1
            if ss not in Rset:
                RuRs.add(ss)
                a, b = pi[s], pi[ss]
                if a + b > 0:
                    merge += (a + b) * h(a / (a + b))
        mass = sum(pi[x] for x in RuRs)
        assert 1 - h(p) <= merge + 1e-9, (k, q, p, merge)
        assert merge <= mass + 1e-12 and mass <= p + 1e-12
        worst = max(worst, (1 - h(p)) - merge)
    return worst


random.seed(20260926)
count = 0
for q in (3, 5, 7, 9, 23):
    for k in range(3, 10):
        H = 1 << (k - 1)
        strategies = [[1] * H, [-1] * H]
        for dens in (0.05, 0.2, 0.5, 0.8):
            for _ in range(6 if k <= 7 else 3):
                strategies.append([-1 if random.random() < dens else 1 for _ in range(H)])
        for sig in strategies:
            audit(k, q, sig)
            count += 1
check(True, f"Theorems 1 and 2 hold on every closed class of {count} strategies (q = 3,5,7,9,23; k = 3..9)")

lo, hi = 0.0, 0.5
for _ in range(200):
    mid = (lo + hi) / 2
    if mid + h(mid) < 1:
        lo = mid
    else:
        hi = mid
p0 = lo
check(abs(p0 - 0.2270922) < 1e-6 and math.log(2) / math.log(23) < p0 < math.log(2) / math.log(21),
      f"p0 = {p0:.7f}; log_23 2 = {math.log(2)/math.log(23):.4f} < p0 < log_21 2 = {math.log(2)/math.log(21):.4f}")
for q, c in ((5, 0.01391), (7, 0.06051), (9, 0.10062)):
    assert abs((1 - h(math.log(2) / math.log(q))) - c) < 5e-5
check(True, "the constants 1 - h(log_q 2) = 0.01391, 0.06051, 0.10062 for q = 5, 7, 9")

# exhaustive: q = 23 has no class-(i) strategy at k <= 4 (floating Bellman-Ford on log weights)
def class_i(k, q, sig):
    N, H = 1 << k, 1 << (k - 1)
    t = targets(k, q, sig)
    w = [math.log(q / 2) if s % 2 else -math.log(2) for s in range(N)]
    D = [0.0] * N
    for it in range(N + 1):
        ch = False
        for s in range(N):
            v = D[s] + w[s]
            for tg in (t[s], t[s] + H):
                if v > D[tg] + 1e-12:
                    D[tg] = v; ch = True
        if not ch:
            return True
    return False

import itertools
for k in (2, 3, 4):
    H = 1 << (k - 1)
    for sig in itertools.product((1, -1), repeat=H):
        assert not class_i(k, 23, list(sig))
        assert not class_i(k, 7, list(sig))
check(True, "exhaustively no class-(i) strategy for q = 7 and q = 23 at k = 2, 3, 4")
