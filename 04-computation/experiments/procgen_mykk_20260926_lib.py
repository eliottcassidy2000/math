#!/usr/bin/env python3
"""
procgen_mykk_20260926_lib.py -- exact helpers for the expanding-cycle feedback lane ("mykk")
(session collatz-procgen-20260922, 2026-09-26).

Objects
-------
* de Bruijn graph B(2,k).  A node is a binary word x = x_0 x_1 ... x_{k-1}, encoded as the integer
  X = sum_j x_j 2^j (x_0 = least significant bit).  Edges: x -> x_1 ... x_{k-1} b, i.e.
  X -> (X >> 1) | (b << (k-1)).  A node is ODD iff x_0 = 1 (iff X is odd).
  (By Terras, this is the parity graph G_0 of T_q(n) = n/2, (q n + 1)/2 in parity-word coordinates, any odd q.)
* A closed walk of length p with a odd nodes has density a/p.  A threshold object decides "expanding":
      Thr('gt', a, b):   density >  a/b
      Thr('ge', a, b):   density >= a/b
      Thr('log', q=q):   q^a > 2^p   (density > c_q = log_q 2; equality impossible)
* FVS_c(k)   = min |R|, R a node set meeting every expanding cycle of B(2,k);
  FVS^odd_c(k): the same with R inside the odd nodes;
  N_c(k)     = number of necklaces (rotation classes of k-words) of density > c;
  nu_c(k)    = max number of pairwise node-disjoint expanding cycles (packing number).
  Always N_c <= nu_c <= FVS_c <= FVS^odd_c.

Everything used as a proof step is exact: integer weights, integer/Fraction comparisons, exact
cyclotomic zero tests.  Floating point is only used for the SIGN of a sine sum that has been shown to be
nonzero, with an explicit separation check (|s| > 1e-7 >> rounding error).
"""
import os
import sys
import time
import json
import gzip
import base64
import hashlib
from fractions import Fraction
from math import comb, gcd, log, sin, pi

import numpy as np
# NOTE: ortools must be imported before highspy (both ship a HiGHS library; loading highspy's first breaks
# ortools' cp_model_helper symbol resolution).
from ortools.sat.python import cp_model  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SCR = os.path.join(ROOT, 'scratch', 'procgen_mykk')


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


def peak_rss_mb():
    import resource
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024.0 * 1024.0) if sys.platform == 'darwin' else r / 1024.0


# ============================================================================ thresholds
class Thr:
    """Expanding-cycle threshold; all decisions exact."""

    def __init__(self, kind, a=None, b=None, q=None):
        assert kind in ('gt', 'ge', 'log')
        self.kind = kind
        if kind == 'log':
            assert q is not None and q % 2 == 1 and q >= 3
            self.q = q
            self.frac = None
        else:
            f = Fraction(a, b)
            self.frac = f
            self.q = None
        self._wcache = {}

    def label(self):
        if self.kind == 'log':
            return "log_%d 2" % self.q
        return ("> %s" if self.kind == 'gt' else ">= %s") % self.frac

    def value(self):
        return log(2) / log(self.q) if self.kind == 'log' else float(self.frac)

    def expanding(self, a, p):
        """exact: is a closed walk with a odd nodes and length p >= 1 expanding?"""
        if self.kind == 'log':
            return self.q ** a > 2 ** p
        if self.kind == 'gt':
            return a * self.frac.denominator > p * self.frac.numerator
        return a * self.frac.denominator >= p * self.frac.numerator

    def weights(self, nmax):
        """integers (u, v), v >= 1, such that for every 1 <= p <= nmax and 0 <= a <= p:
        expanding(a, p)  <=>  a*v - p*u > 0.   (Node weight v*[odd] - u; a closed walk of length <= nmax
        is expanding iff its weight is positive.  A positive closed walk contains a positive SIMPLE cycle.)"""
        if nmax in self._wcache:
            return self._wcache[nmax]
        if self.kind == 'gt':
            res = (self.frac.numerator, self.frac.denominator)
        else:
            # largest fraction u/v with v <= nmax that is NOT expanding-level: u/v < c (log) or u/v < a/b (ge)
            best = None
            for v in range(1, nmax + 1):
                if self.kind == 'ge':
                    u = (v * self.frac.numerator - 1) // self.frac.denominator
                else:
                    u = int(v * log(2) / log(self.q))
                    while self.q ** (u + 1) < 2 ** v:
                        u += 1
                    while u >= 0 and not self.q ** u < 2 ** v:
                        u -= 1
                if u < 0:
                    continue
                f = Fraction(u, v)
                if best is None or f > best:
                    best = f
            res = (best.numerator, best.denominator)
        # exact self-test on all (a,p), p <= min(nmax, 64); the general statement is by maximality of u/v
        u, v = res
        for p in range(1, min(nmax, 64) + 1):
            for a in range(p + 1):
                check(self.expanding(a, p) == (a * v - p * u > 0), "weights self-test %s" % self.label())
        self._wcache[nmax] = res
        return res

    def min_ones(self, k):
        """least a with expanding(a, k) (necklaces of length k with >= this many ones are expanding)"""
        for a in range(k + 1):
            if self.expanding(a, k):
                return a
        return k + 1


# ============================================================================ de Bruijn basics
def mask(k):
    return (1 << k) - 1


def succ(X, b, k):
    return (X >> 1) | (b << (k - 1))


def rotL(X, k):
    """x_0 x_1 .. x_{k-1} -> x_1 .. x_{k-1} x_0  (the necklace successor, = succ(X, x_0))"""
    return (X >> 1) | ((X & 1) << (k - 1))


def rotR(X, k):
    """x_0 .. x_{k-1} -> x_{k-1} x_0 .. x_{k-2}  (inverse of rotL)"""
    return ((X << 1) & mask(k)) | (X >> (k - 1))


def word_str(X, k):
    return ''.join(str((X >> j) & 1) for j in range(k))


def popcount_arr(A):
    A = A.astype(np.uint64)
    c = np.zeros(A.shape, dtype=np.int64)
    while True:
        nz = A != 0
        if not nz.any():
            break
        c += (A & np.uint64(1)).astype(np.int64)
        A = A >> np.uint64(1)
    return c


_DB = {}


def db_arrays(k):
    """numpy arrays: pred0, pred1 (the two predecessors of each node), odd indicator"""
    if k not in _DB:
        n = 1 << k
        X = np.arange(n, dtype=np.int64)
        p0 = (X << 1) & (n - 1)
        p1 = p0 | 1
        _DB[k] = (p0, p1, (X & 1).astype(np.int64))
    return _DB[k]


def is_edge(x, y, k):
    return (x >> 1) == (y & (mask(k) >> 1))


def cycle_ok(cyc, k):
    """cyc = list of nodes in walk order: a closed walk of B(2,k) with distinct nodes (simple cycle)"""
    if len(cyc) == 0 or len(set(cyc)) != len(cyc):
        return False
    return all(is_edge(cyc[i], cyc[(i + 1) % len(cyc)], k) for i in range(len(cyc)))


def cycle_from_word(u, l, k):
    """nodes (k-windows) of the periodic sequence (u_0 .. u_{l-1})^infinity, positions 0..l-1"""
    bits = [(u >> t) & 1 for t in range(l)]
    out = []
    for t in range(l):
        X = 0
        for j in range(k):
            X |= bits[(t + j) % l] << j
        out.append(X)
    return out


def word_of_cycle(cyc):
    """the periodic sequence read along a closed walk: z_t = first letter of node t"""
    return [x & 1 for x in cyc]


# ============================================================================ necklaces
_NECK = {}


def necklaces(k):
    """(neck_id array, list of necklaces as node lists in rotL order starting at the least integer)"""
    if k in _NECK:
        return _NECK[k]
    n = 1 << k
    nid = np.full(n, -1, dtype=np.int64)
    neck = []
    for X in range(n):
        if nid[X] >= 0:
            continue
        orb = [X]
        Y = rotL(X, k)
        while Y != X:
            orb.append(Y)
            Y = rotL(Y, k)
        for Y in orb:
            nid[Y] = len(neck)
        neck.append(orb)
    _NECK[k] = (nid, neck)
    return _NECK[k]


def Z_formula(k):
    """number of binary necklaces of length k: (1/k) sum_{d|k} phi(d) 2^(k/d)"""
    def phi(m):
        return sum(1 for i in range(1, m + 1) if gcd(i, m) == 1)
    s = sum(phi(d) * 2 ** (k // d) for d in range(1, k + 1) if k % d == 0)
    assert s % k == 0
    return s // k


def necklaces_with_ones(k, j):
    """number of binary necklaces of length k with exactly j ones (Burnside)"""
    def phi(m):
        return sum(1 for i in range(1, m + 1) if gcd(i, m) == 1)
    # (1/k) sum_{d | gcd(k,j)} phi(d) C(k/d, j/d)
    g = gcd(k, j) if j > 0 else k
    s = sum(phi(d) * comb(k // d, j // d) for d in range(1, g + 1) if g % d == 0)
    assert s % k == 0
    return s // k


def N_count(k, thr):
    """N_c(k): necklaces of length k (a closed walk of length k: the rotation orbit, walked k/period times)
    whose density is > c"""
    return sum(necklaces_with_ones(k, j) for j in range(k + 1) if thr.expanding(j, k))


def necklace_constraints(k, thr):
    """the expanding necklaces as node lists (these are node-disjoint expanding cycles)"""
    nid, neck = necklaces(k)
    out = []
    for orb in neck:
        a = sum(x & 1 for x in orb)
        if thr.expanding(a, len(orb)):
            out.append(list(orb))
    return out


# ============================================================================ positive cycles
NEG = -(1 << 60)


def positive_cycle(k, alive, w):
    """Longest-path Bellman-Ford on the subgraph of B(2,k) induced by `alive` (bool array), with node
    weights w (int64; the weight of a closed walk is the sum over its nodes).  Returns a positive-weight
    simple cycle (list of nodes in walk order) or None; if None, also returns the potential d (max weight of
    a walk ending at each node, i.e. d(y) >= d(x) + w(x) on every alive edge x -> y).
    Early detection: every 16 rounds the parent pointers of a few just-improved nodes are followed; a cycle in
    the parent graph is returned if its weight is positive (always the case for Bellman-Ford parent cycles)."""
    n = 1 << k
    p0, p1, _ = db_arrays(k)
    d = np.where(alive, 0, NEG).astype(np.int64)
    par = np.full(n, -1, dtype=np.int64)
    a0 = alive[p0]
    a1 = alive[p1]
    na = int(alive.sum())

    def parent_cycle(v):
        seen = {}
        path = []
        while v >= 0 and v not in seen:
            seen[v] = len(path)
            path.append(v)
            v = int(par[v])
        if v < 0:
            return None
        cyc = path[seen[v]:]
        cyc.reverse()
        return cyc

    for it in range(na + 2):
        c0 = np.where(a0, d[p0] + w[p0], NEG)
        c1 = np.where(a1, d[p1] + w[p1], NEG)
        best = np.maximum(c0, c1)
        arg = np.where(c1 > c0, p1, p0)
        imp = alive & (best > d)
        if not imp.any():
            return None, d
        d[imp] = best[imp]
        par[imp] = arg[imp]
        if it % 16 == 15 or it >= na:
            cand = np.nonzero(imp)[0]
            for v in cand[:3].tolist():
                cyc = parent_cycle(v)
                if cyc is not None and sum(int(w[x]) for x in cyc) > 0 and len(set(cyc)) == len(cyc):
                    return cyc, None
        if it >= na:
            v = int(np.nonzero(imp)[0][0])
            for _ in range(n + 1):
                v = int(par[v])
            cyc = [v]
            u = int(par[v])
            while u != v:
                cyc.append(u)
                u = int(par[u])
            cyc.reverse()
            return cyc, None
    raise RuntimeError("Bellman-Ford did not terminate")


def node_weights(k, thr, nmax=None):
    u, v = thr.weights((1 << k) if nmax is None else nmax)
    _, _, odd = db_arrays(k)
    return (v * odd - u).astype(np.int64), (u, v)


def has_expanding_cycle(k, thr, removed):
    """exact: does B(2,k) - removed contain an expanding cycle?  Returns (cycle or None, potential)."""
    n = 1 << k
    alive = np.ones(n, dtype=bool)
    if len(removed):
        alive[np.asarray(list(removed), dtype=np.int64)] = False
    w, _ = node_weights(k, thr)
    cyc, d = positive_cycle(k, alive, w)
    if cyc is not None:
        check(cycle_ok(cyc, k), "BF returned a non-cycle")
        a = sum(x & 1 for x in cyc)
        check(thr.expanding(a, len(cyc)), "BF returned a non-expanding cycle")
        check(all(alive[x] for x in cyc), "BF cycle meets removed set")
    return cyc, d


def verify_potential(k, thr, removed, d):
    """exact re-check of a potential certificate: d(y) >= d(x) + w(x) on every edge x->y of B(2,k)-removed.
    Summing around any closed walk of B(2,k)-removed gives weight <= 0, i.e. no expanding cycle."""
    n = 1 << k
    alive = np.ones(n, dtype=bool)
    if len(removed):
        alive[np.asarray(list(removed), dtype=np.int64)] = False
    w, _ = node_weights(k, thr)
    p0, p1, _ = db_arrays(k)
    ok = True
    for p in (p0, p1):
        m = alive & alive[p]
        ok &= bool(np.all(d[m] >= d[p][m] + w[p][m]))
    return ok


def disjoint_positive_cycles(k, thr, removed, limit=10 ** 9):
    """repeatedly extract an expanding cycle and delete its nodes; returns node-disjoint expanding cycles
    of B(2,k) - removed"""
    n = 1 << k
    alive = np.ones(n, dtype=bool)
    if len(removed):
        alive[np.asarray(list(removed), dtype=np.int64)] = False
    w, _ = node_weights(k, thr)
    out = []
    while len(out) < limit:
        cyc, _ = positive_cycle(k, alive, w)
        if cyc is None:
            break
        out.append(cyc)
        alive[np.asarray(cyc, dtype=np.int64)] = False
    return out


def shortest_positive_cycles(k, alive, w, lmin=1, lmax=64, want=200, extra_lengths=2):
    """Expanding (positive-weight) simple cycles of minimal length in the subgraph induced by `alive`:
    all-pairs DP over walk lengths (D_l[s, y] = max weight of a walk of exactly l steps from s to y; the
    weight of a step is w(source node)).  Scans l = 1..lmax; from the first length l0 with some
    D_l[s, s] > 0 it collects cycles for lengths l0 .. l0+extra_lengths (at most `want`, distinct start
    nodes), reconstructing each by a single-source DP with predecessors.  A positive closed walk may repeat
    nodes; it is decomposed and its densest simple sub-cycle (again positive) is returned.
    int32 arithmetic (|weights| * lmax < 2^30 is asserted); memory ~ 12 m^2 bytes, m = #alive."""
    n = 1 << k
    p0, p1, _ = db_arrays(k)
    idx = np.nonzero(alive)[0]
    m = len(idx)
    if m == 0:
        return []
    check(int(np.abs(w).max()) * (lmax + 1) < (1 << 29), "int32 overflow risk in shortest_positive_cycles")
    pos = np.full(n, -1, dtype=np.int64)
    pos[idx] = np.arange(m)
    q0 = np.where(alive[p0[idx]], pos[p0[idx]], -1)
    q1 = np.where(alive[p1[idx]], pos[p1[idx]], -1)
    wc = w[idx].astype(np.int32)
    NEGI = np.int32(-(1 << 30))
    g0 = np.maximum(q0, 0)
    g1 = np.maximum(q1, 0)
    wq0 = np.where(q0 >= 0, wc[g0], 0).astype(np.int32)
    wq1 = np.where(q1 >= 0, wc[g1], 0).astype(np.int32)
    ok0 = q0 >= 0
    ok1 = q1 >= 0
    D = np.full((m, m), NEGI, dtype=np.int32)
    D[np.arange(m), np.arange(m)] = 0
    found = []
    l0 = None
    starts_done = set()
    for l in range(1, lmax + 1):
        A = D[:, g0]
        A += wq0
        A[:, ~ok0] = NEGI
        B = D[:, g1]
        B += wq1
        B[:, ~ok1] = NEGI
        np.maximum(A, B, out=A)
        del B
        A[A < NEGI // 2] = NEGI
        D = A
        diag = np.diagonal(D)
        pos_s = np.nonzero(diag > 0)[0]
        if len(pos_s) and l >= lmin:
            if l0 is None:
                l0 = l
            for s0 in pos_s:
                if len(found) >= want:
                    break
                s0 = int(s0)
                if s0 in starts_done:
                    continue
                cyc = _reconstruct_positive(idx, q0, q1, wc, s0, l)
                if cyc is None:
                    continue
                starts_done.update(int(pos[x]) for x in cyc)
                found.append(cyc)
        if l0 is not None and (l >= l0 + extra_lengths or len(found) >= want):
            break
    return found


def _reconstruct_positive(idx, q0, q1, wc, s, l):
    """single-source DP from compressed node s for exactly l steps with predecessors; returns the densest
    simple cycle (original node ids) in the best closed walk s -> s of length l (None if not positive)"""
    m = len(idx)
    NEGI = -(1 << 40)
    d = np.full(m, NEGI, dtype=np.int64)
    d[s] = 0
    g0 = np.maximum(q0, 0)
    g1 = np.maximum(q1, 0)
    w0 = np.where(q0 >= 0, wc[g0], 0).astype(np.int64)
    w1 = np.where(q1 >= 0, wc[g1], 0).astype(np.int64)
    pred = np.full((l + 1, m), -1, dtype=np.int64)
    for t in range(1, l + 1):
        c0 = np.where((q0 >= 0) & (d[g0] > NEGI // 2), d[g0] + w0, NEGI)
        c1 = np.where((q1 >= 0) & (d[g1] > NEGI // 2), d[g1] + w1, NEGI)
        nd = np.maximum(c0, c1)
        pr = np.where(c1 > c0, q1, q0)
        pr[nd <= NEGI // 2] = -1
        pred[t] = pr
        d = nd
    if d[s] <= 0:
        return None
    walk = [s]
    y = s
    for t in range(l, 0, -1):
        y = int(pred[t][y])
        walk.append(y)
    walk.reverse()
    nodes = [int(idx[x]) for x in walk[:-1]]
    return _densest_simple_subcycle(nodes)


def _densest_simple_subcycle(nodes):
    """decompose a closed walk (nodes in walk order, last -> first) into simple cycles and return the one of
    maximum odd density (its density is >= the walk's density: mediant property)"""
    stack = []
    where = {}
    cycles = []
    for x in nodes + nodes[:1]:
        if x in where:
            i = where[x]
            cyc = stack[i:]
            cycles.append(cyc)
            for y in cyc:
                del where[y]
            stack = stack[:i]
        where[x] = len(stack)
        stack.append(x)
    best = None
    for c in cycles:
        dens = Fraction(sum(v & 1 for v in c), len(c))
        if best is None or dens > best[0]:
            best = (dens, c)
    return best[1] if best else None


def max_density_karp(k, alive):
    """exact maximum cycle density (Fraction) of the subgraph induced by alive (Karp, integer DP)"""
    n = 1 << k
    p0, p1, odd = db_arrays(k)
    idx = np.nonzero(alive)[0]
    m = len(idx)
    if m == 0:
        return None
    NEGI = -(1 << 40)
    D = np.full((m + 1, n), NEGI, dtype=np.int64)
    D[0, alive] = 0
    a0 = alive[p0]
    a1 = alive[p1]
    for t in range(1, m + 1):
        prev = D[t - 1]
        c0 = np.where(a0 & (prev[p0] > NEGI // 2), prev[p0] + odd[p0], NEGI)
        c1 = np.where(a1 & (prev[p1] > NEGI // 2), prev[p1] + odd[p1], NEGI)
        D[t] = np.where(alive, np.maximum(c0, c1), NEGI)
    best = None
    for v in idx:
        if D[m, v] <= NEGI // 2:
            continue
        worst = None
        for t in range(m):
            if D[t, v] <= NEGI // 2:
                continue
            f = Fraction(int(D[m, v] - D[t, v]), m - t)
            if worst is None or f < worst:
                worst = f
        if worst is not None and (best is None or worst > best):
            best = worst
    return best


# ============================================================================ short-cycle pools
def primitive_necklace_reps(l):
    """all binary words of length l (LSB = first letter) that are the unique least rotation of a primitive
    necklace, as a numpy int64 array (vectorized; memory ~ 40 * 2^l bytes)"""
    full = (1 << l) - 1
    out = []
    CH = 1 << 20
    for s in range(0, 1 << l, CH):
        W = np.arange(s, min(1 << l, s + CH), dtype=np.int64)
        ok = np.ones(len(W), dtype=bool)
        for r in range(1, l):
            rot = ((W >> r) | (W << (l - r))) & full
            ok &= W < rot
        out.append(W[ok])
    return np.concatenate(out) if out else np.zeros(0, dtype=np.int64)


def pool_cycles(k, thr, lmax, simple_only=True):
    """all simple expanding cycles of B(2,k) of length l <= lmax, as a dict l -> int64 array (rows = cycles,
    columns = nodes in walk order).  Built from primitive necklaces of length l (period-l sequences)."""
    pools = {}
    for l in range(1, lmax + 1):
        reps = primitive_necklace_reps(l)
        if len(reps) == 0:
            continue
        ones = popcount_arr(reps)
        keep = np.array([thr.expanding(int(a), l) for a in range(l + 1)])[ones]
        reps = reps[keep]
        if len(reps) == 0:
            continue
        full = (1 << l) - 1
        cols = []
        for t in range(l):
            R = ((reps >> t) | (reps << (l - t))) & full if t else reps.copy()
            # periodic extension to >= k bits
            E = R.copy()
            span = l
            while span < k:
                E = E | (R << span)
                span += l
            cols.append(E & ((1 << k) - 1))
        M = np.stack(cols, axis=1)
        if simple_only and l > 1:
            S = np.sort(M, axis=1)
            dup = (S[:, 1:] == S[:, :-1]).any(axis=1)
            M = M[~dup]
        if len(M):
            pools[l] = M
    return pools


def pool_violated(pools, removed_mask, cap=None):
    """cycles of the pools that avoid the removed set (bool array); shortest first"""
    out = []
    for l in sorted(pools):
        M = pools[l]
        hit = removed_mask[M].any(axis=1)
        sel = M[~hit]
        for row in sel:
            out.append([int(x) for x in row])
            if cap is not None and len(out) >= cap:
                return out
    return out


# ============================================================================ Mykkeltveit construction
def cyclotomic_poly(k):
    """integer coefficients (low degree first) of the k-th cyclotomic polynomial"""
    polys = {}

    def div(a, b):
        a = list(a)
        out = [0] * (len(a) - len(b) + 1)
        for i in range(len(a) - len(b), -1, -1):
            c = a[i + len(b) - 1] // b[-1]
            out[i] = c
            for j in range(len(b)):
                a[i + j] -= c * b[j]
        assert all(x == 0 for x in a[:len(b) - 1])
        return out
    for d in range(1, k + 1):
        if k % d:
            continue
        p = [-1] + [0] * (d - 1) + [1]
        for e in range(1, d):
            if d % e == 0:
                p = div(p, polys[e])
        polys[d] = p
    return polys[k]


def poly_rem(a, b):
    a = list(a)
    while len(a) >= len(b):
        c = a[-1]
        if c != 0:
            s = len(a) - len(b)
            for j in range(len(b)):
                a[s + j] -= c * b[j]   # b monic
        a.pop()
    return a


def zero_in_cyclotomic(coeffs, k, phi=None):
    """exact test: sum_j coeffs[j] * omega^j == 0 for omega = exp(2 pi i / k)"""
    if phi is None:
        phi = cyclotomic_poly(k)
    return all(c == 0 for c in poly_rem(coeffs, phi))


def sine_data(k):
    """for every node X: exact flags w(X)==0, s(X)==0, and float s(X) = sum_j x_j sin(2 pi j/k)"""
    n = 1 << k
    phi = cyclotomic_poly(k)
    sins = [sin(2 * pi * j / k) for j in range(k)]
    s = np.zeros(n)
    wzero = np.zeros(n, dtype=bool)
    szero = np.zeros(n, dtype=bool)
    for X in range(n):
        bits = [(X >> j) & 1 for j in range(k)]
        s[X] = sum(b * sv for b, sv in zip(bits, sins))
        wzero[X] = zero_in_cyclotomic(bits, k, phi)
        szero[X] = zero_in_cyclotomic([bits[j] - bits[(-j) % k] for j in range(k)], k, phi)
    # separation: exact zeros are float-small, exact non-zeros are float-large
    check(np.all(np.abs(s[szero]) < 1e-9), "exact sine zero but float not small (k=%d)" % k)
    if (~szero).any():
        check(np.min(np.abs(s[~szero])) > 1e-7, "a nonzero sine sum is too small to sign safely (k=%d)" % k)
    return s, wzero, szero


def sign_s(s, szero):
    sg = np.sign(s).astype(np.int64)
    sg[szero] = 0
    return sg


def mykkeltveit_set(k, zero_choice='odd_min', thr=None):
    """Mykkeltveit's set M_k (k >= 3):
       S* = { y : s(y) <= 0 < s(rotR(y)) }   (the node where the necklace's weight vector, rotating clockwise by
                                              2 pi/k per step, has just crossed the positive real axis)
       plus one node of every zero-weight necklace (w == 0 on the whole necklace).
       If thr is given, zero-weight necklaces that are not expanding are omitted (they never need hitting).
    Returns (sorted list of nodes, S* list, list of zero-weight necklace ids)."""
    n = 1 << k
    s, wzero, szero = sine_data(k)
    sg = sign_s(s, szero)
    nid, neck = necklaces(k)
    X = np.arange(n, dtype=np.int64)
    RR = ((X << 1) & (n - 1)) | (X >> (k - 1))
    Sstar = np.nonzero((sg <= 0) & (sg[RR] > 0))[0]
    zero_necks = [i for i, orb in enumerate(neck) if all(wzero[x] for x in orb)]
    # consistency: a necklace is zero-weight iff it has no S* node; otherwise exactly one S* node
    cnt = np.bincount(nid[Sstar], minlength=len(neck))
    for i, orb in enumerate(neck):
        if i in set(zero_necks):
            check(cnt[i] == 0, "zero-weight necklace with an S* node")
        else:
            check(cnt[i] == 1, "nonzero necklace without exactly one S* node (k=%d)" % k)
    chosen = []
    for i in zero_necks:
        orb = neck[i]
        if thr is not None:
            a = sum(x & 1 for x in orb)
            if not thr.expanding(a, len(orb)):
                continue
        if zero_choice == 'odd_min':
            odds = [x for x in orb if x & 1]
            chosen.append(min(odds) if odds else min(orb))
        else:
            chosen.append(min(orb))
    M = sorted(set(int(x) for x in Sstar) | set(chosen))
    return M, [int(x) for x in Sstar], zero_necks


def is_acyclic(k, removed):
    """exact: B(2,k) - removed has no cycle at all (Kahn's algorithm)"""
    n = 1 << k
    alive = np.ones(n, dtype=bool)
    if len(removed):
        alive[np.asarray(list(removed), dtype=np.int64)] = False
    p0, p1, _ = db_arrays(k)
    X = np.arange(n, dtype=np.int64)
    s0 = X >> 1
    s1 = s0 | (1 << (k - 1))
    indeg = (alive[p0] & alive).astype(np.int64) + (alive[p1] & alive).astype(np.int64)
    indeg[~alive] = 0
    # the loops 0..0 and 1..1: p0(0)=0 counts itself; handled uniformly
    stack = [int(x) for x in np.nonzero(alive & (indeg == 0))[0]]
    seen = 0
    removed_flag = np.zeros(n, dtype=bool)
    while stack:
        x = stack.pop()
        removed_flag[x] = True
        seen += 1
        for y in (int(s0[x]), int(s1[x])):
            if alive[y] and not removed_flag[y]:
                indeg[y] -= 1
                if indeg[y] == 0:
                    stack.append(y)
    return seen == int(alive.sum())


def longest_path_acyclic(k, removed):
    """number of nodes on a longest path of the DAG B(2,k) - removed (None if it has a cycle)"""
    n = 1 << k
    alive = np.ones(n, dtype=bool)
    if len(removed):
        alive[np.asarray(list(removed), dtype=np.int64)] = False
    p0, p1, _ = db_arrays(k)
    X = np.arange(n, dtype=np.int64)
    s0 = X >> 1
    s1 = s0 | (1 << (k - 1))
    indeg = (alive[p0] & alive).astype(np.int64) + (alive[p1] & alive).astype(np.int64)
    indeg[~alive] = 0
    depth = np.where(alive, 1, 0).astype(np.int64)
    stack = [int(x) for x in np.nonzero(alive & (indeg == 0))[0]]
    seen = 0
    while stack:
        x = stack.pop()
        seen += 1
        for y in (int(s0[x]), int(s1[x])):
            if alive[y]:
                if depth[x] + 1 > depth[y]:
                    depth[y] = depth[x] + 1
                indeg[y] -= 1
                if indeg[y] == 0:
                    stack.append(y)
    if seen != int(alive.sum()):
        return None
    return int(depth.max())


# ============================================================================ hitting-set master problems
def solve_cover_highs(cycles, allowed, time_limit=None, threads=1, incumbent=None, full=False):
    """min |R|, R subset of `allowed` (sorted node list), meeting every cycle (HiGHS MIP).
    Returns (value, R) or (None, None) if infeasible; with full=True returns
    (value, R, optimal, lower_bound) where on a time-out R is the best incumbent and lower_bound = the
    rounded-up MIP dual bound (a valid lower bound for the minimum hitting set of `cycles`)."""
    import highspy
    pos = {v: i for i, v in enumerate(allowed)}
    h = highspy.Highs()
    h.silent()
    h.setOptionValue('threads', threads)
    h.setOptionValue('mip_rel_gap', 0.0)
    h.setOptionValue('mip_abs_gap', 0.5)
    if time_limit is not None:
        h.setOptionValue('time_limit', float(time_limit))
    nv = len(allowed)
    h.addVars(nv, np.zeros(nv), np.ones(nv))
    h.changeColsCost(nv, np.arange(nv, dtype=np.int32), np.ones(nv))
    h.changeColsIntegrality(nv, np.arange(nv, dtype=np.int32),
                            np.array([highspy.HighsVarType.kInteger] * nv))
    starts, idxs = [], []
    for c in cycles:
        cols = sorted(set(pos[v] for v in c if v in pos))
        if not cols:
            return (None, None, True, None) if full else (None, None)
        starts.append(len(idxs))
        idxs.extend(cols)
    nr = len(starts)
    if nr:
        h.addRows(nr, np.ones(nr), np.full(nr, highspy.kHighsInf), len(idxs),
                  np.array(starts, dtype=np.int32), np.array(idxs, dtype=np.int32), np.ones(len(idxs)))
    h.setMinimize()
    if incumbent is not None:
        vals = np.zeros(nv)
        for v in incumbent:
            if v in pos:
                vals[pos[v]] = 1.0
        h.setSolution(nv, np.arange(nv, dtype=np.int32), vals)
    h.run()
    st = h.getModelStatus()
    if st == highspy.HighsModelStatus.kInfeasible:
        return (None, None, True, None) if full else (None, None)
    optimal = (st == highspy.HighsModelStatus.kOptimal)
    if not optimal and not full:
        raise RuntimeError("HiGHS status %s" % h.modelStatusToString(st))
    sol = h.getSolution().col_value
    R = [allowed[i] for i in range(nv) if sol[i] > 0.5]
    Rs = set(R)
    for c in cycles:
        assert any(v in Rs for v in c)
    if not full:
        return len(R), R
    if optimal:
        return len(R), R, True, len(R)
    db = h.getInfo().mip_dual_bound
    return len(R), R, False, int(np.ceil(db - 1e-6))


def minimalize(k, thr, R, pools=None, order=None):
    """drop nodes of a feasible set R (no expanding cycle in B(2,k) - R) while it stays feasible;
    returns an inclusion-minimal feasible subset"""
    n = 1 << k
    R = list(R)
    if order is None:
        order = sorted(R, key=lambda x: (bin(x).count('1'), x))
    cur = set(R)
    for r in order:
        trial = cur - {r}
        rem = np.zeros(n, dtype=bool)
        rem[np.asarray(sorted(trial), dtype=np.int64)] = True
        if pools is not None and pool_violated(pools, rem, cap=1):
            continue
        cyc, _ = has_expanding_cycle(k, thr, sorted(trial))
        if cyc is None:
            cur = trial
    return sorted(cur)


def cover_infeasible_below_cpsat(cycles, allowed, bound, workers=2, time_limit=None):
    """independent re-proof with CP-SAT: there is no R subset allowed, |R| <= bound, meeting all cycles.
    Returns True iff CP-SAT proves infeasibility."""
    from ortools.sat.python import cp_model
    m = cp_model.CpModel()
    x = {v: m.NewBoolVar("x%d" % v) for v in allowed}
    for c in cycles:
        lits = [x[v] for v in set(c) if v in x]
        if not lits:
            return True
        m.AddBoolOr(lits)
    m.Add(sum(x.values()) <= bound)
    s = cp_model.CpSolver()
    s.parameters.num_search_workers = workers
    if time_limit:
        s.parameters.max_time_in_seconds = float(time_limit)
    st = s.Solve(m)
    return st == cp_model.INFEASIBLE


def max_packing_cpsat(cycles, workers=2, time_limit=60, lb_hint=None):
    """max number of pairwise node-disjoint cycles among `cycles` (CP-SAT). Returns (value, chosen idx, optimal?)"""
    from ortools.sat.python import cp_model
    m = cp_model.CpModel()
    y = [m.NewBoolVar("y%d" % i) for i in range(len(cycles))]
    bynode = {}
    for i, c in enumerate(cycles):
        for v in set(c):
            bynode.setdefault(v, []).append(i)
    for v, lst in bynode.items():
        if len(lst) > 1:
            m.AddAtMostOne([y[i] for i in lst])
    m.Maximize(sum(y))
    s = cp_model.CpSolver()
    s.parameters.num_search_workers = workers
    s.parameters.max_time_in_seconds = float(time_limit)
    st = s.Solve(m)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return 0, [], False
    ch = [i for i in range(len(cycles)) if s.Value(y[i])]
    return len(ch), ch, st == cp_model.OPTIMAL


# ============================================================================ implicit hitting set
def _separate(k, thr, pools, w, removed, per_round, dp_lmax, use_dp=True):
    """violated expanding cycles for the removed set (bool mask): pool first, then DP, then Bellman-Ford"""
    new = pool_violated(pools, removed, cap=per_round)
    src = 'pool'
    if len(new) < per_round // 4 and use_dp:
        more = shortest_positive_cycles(k, ~removed, w, lmax=dp_lmax, want=per_round)
        if more:
            new.extend(more)
            src = 'dp'
    if not new:
        R = np.nonzero(removed)[0].tolist()
        new = disjoint_positive_cycles(k, thr, R, limit=per_round)
        src = 'bf'
    return new, src


def _greedy_repair(R, cycles_new, allowed_set):
    """add nodes (greedily, most-hitting first) so that every cycle in cycles_new is hit"""
    R = set(R)
    todo = [c for c in cycles_new if not any(x in R for x in c)]
    while todo:
        cnt = {}
        for c in todo:
            for x in c:
                if x in allowed_set:
                    cnt[x] = cnt.get(x, 0) + 1
        best = max(cnt, key=lambda x: (cnt[x], x & 1, -x))
        R.add(best)
        todo = [c for c in todo if best not in c]
    return R


def ihs_fvs(k, thr, odd_only=False, lmax=None, per_round=400, log=None, time_limit=None,
            seed_cycles=None, master_threads=1, verbose=False, ckpt=None, dp_lmax=96, repair=True,
            reverse=True, mip_time=None, incumbent=None):
    """Exact FVS_c(k) (or FVS^odd_c(k)) by an implicit hitting set:
       master = min hitting set of the expanding cycles found so far (HiGHS MIP; a lower bound);
       separation = exact expanding-cycle search in B(2,k) - R: (1) the pool of all simple expanding cycles
       of length <= lmax, (2) minimal-length expanding cycles by an all-pairs length DP (lengths <= dp_lmax),
       (3) Bellman-Ford (any length), extracting node-disjoint cycles.
       After each master solve a greedy REPAIR loop (MaxHS-style non-optimal hitting sets) extends the master
       solution until it is feasible, collecting the violated cycles on the way and giving an upper bound.
       Terminates when the master optimum is feasible (then optimal) or equals the best feasible size.
    ckpt: optional path; the cycle list is saved there (json) every round and reloaded if present.
    Returns dict(value, R, cycles (= LB certificate), rounds, optimal, lb, ub, time)."""
    t0 = time.time()
    n = 1 << k
    allowed = [v for v in range(n) if (v & 1 or not odd_only)]
    allowed_set = set(allowed)
    if lmax is None:
        lmax = min(2 * k + 2, 20)
    pools = pool_cycles(k, thr, lmax)
    w, _ = node_weights(k, thr)
    cycles = [list(c) for c in necklace_constraints(k, thr)]
    if seed_cycles:
        cycles.extend([list(c) for c in seed_cycles])
    ck_best, ck_lb = None, 0
    if ckpt and os.path.exists(ckpt):
        with open(ckpt) as f:
            ckd = json.load(f)
        cycles.extend(ckd['cycles'])
        ck_best, ck_lb = ckd.get('best_R'), ckd.get('lb', 0)
    uniq = []
    seen = set()
    for c in cycles:
        key = tuple(sorted(c))
        if key not in seen:
            seen.add(key)
            uniq.append(list(c))
    cycles = uniq
    rounds = 0
    lb = 0
    best_R = None
    for cand in (ck_best, incumbent):
        if cand is None:
            continue
        cyc0, _ = has_expanding_cycle(k, thr, cand)
        if cyc0 is None and (not odd_only or all(x & 1 for x in cand)):
            if best_R is None or len(cand) < len(best_R):
                best_R = sorted(cand)

    rev_tab = None
    if reverse:
        rev_tab = [int(format(x, '0%db' % k)[::-1], 2) for x in range(n)]

    def add(cs):
        a = 0
        for c in cs:
            variants = [list(c)]
            if rev_tab is not None:
                variants.append([rev_tab[x] for x in reversed(c)])
            for cc in variants:
                key = tuple(sorted(cc))
                if key not in seen:
                    seen.add(key)
                    cycles.append(cc)
                    a += 1
        return a

    while True:
        rounds += 1
        if best_R is not None:
            inc = [x for x in best_R]
        else:
            inc = None
        val, R, opt, bnd = solve_cover_highs(cycles, allowed, threads=master_threads, incumbent=inc,
                                             time_limit=mip_time, full=True)
        check(val is not None, "master infeasible")
        lb = max(lb, bnd)
        removed = np.zeros(n, dtype=bool)
        removed[np.asarray(R, dtype=np.int64)] = True
        new, src = _separate(k, thr, pools, w, removed, per_round, dp_lmax)
        added = add(new)
        if not new:
            if best_R is None or len(R) < len(best_R):
                best_R = sorted(R)
            if opt:
                if verbose and log:
                    log("    round %d: LB=%d optimal  |cycles|=%d  (%.1fs)" % (rounds, val, len(cycles),
                                                                          time.time() - t0))
                return dict(value=val, R=sorted(R), cycles=cycles, rounds=rounds, optimal=True,
                            lb=val, ub=val, time=time.time() - t0)
        else:
            check(added > 0, "separation returned only known cycles (should be impossible)")
        # repair loop: extend R greedily until feasible, collecting cycles; then make it inclusion-minimal
        rep_rounds = 0
        if repair and new:
            Rr = _greedy_repair(R, new, allowed_set)
            while True:
                rep_rounds += 1
                rem2 = np.zeros(n, dtype=bool)
                rem2[np.asarray(sorted(Rr), dtype=np.int64)] = True
                new2, _ = _separate(k, thr, pools, w, rem2, per_round, dp_lmax, use_dp=False)
                if not new2:
                    break
                add(new2)
                Rr = _greedy_repair(Rr, new2, allowed_set)
            Rr = minimalize(k, thr, sorted(Rr), pools=pools)
            if best_R is None or len(Rr) < len(best_R):
                best_R = sorted(Rr)
        if ckpt:
            with open(ckpt + '.tmp', 'w') as f:
                json.dump({'cycles': cycles, 'best_R': best_R, 'lb': lb}, f)
            os.replace(ckpt + '.tmp', ckpt)
        if verbose and log:
            log("    round %d: LB=%d%s UB=%s |cycles|=%d new=%d (%s) repair %d  (%.1fs)" % (
                rounds, lb, '' if opt else '(bound)', len(best_R) if best_R else '-', len(cycles), added, src,
                rep_rounds, time.time() - t0))
        if best_R is not None and len(best_R) == lb:
            cyc, _ = has_expanding_cycle(k, thr, best_R)
            check(cyc is None, "best set not feasible")
            return dict(value=lb, R=best_R, cycles=cycles, rounds=rounds, optimal=True,
                        lb=lb, ub=lb, time=time.time() - t0)
        if time_limit is not None and time.time() - t0 > time_limit:
            return dict(value=None, R=best_R, cycles=cycles, rounds=rounds, optimal=False, lb=lb,
                        ub=len(best_R) if best_R else None, time=time.time() - t0)


def rho_max_removed(k, R):
    """exact maximum cycle density of B(2,k) - R (None if acyclic).  A float Karp computation proposes a/p
    (denominator <= number of alive nodes); it is accepted only after two exact Bellman-Ford checks
    (no cycle of density > a/p, some cycle of density >= a/p); otherwise the exact integer Karp is used."""
    n = 1 << k
    alive = np.ones(n, dtype=bool)
    if len(R):
        alive[np.asarray(list(R), dtype=np.int64)] = False
    m = int(alive.sum())
    if m == 0 or is_acyclic(k, R):
        return None
    p0, p1, odd = db_arrays(k)
    NEGF = -1e18
    D = np.full((m + 1, n), NEGF)
    D[0, alive] = 0.0
    a0 = alive[p0]
    a1 = alive[p1]
    for t in range(1, m + 1):
        prev = D[t - 1]
        c0 = np.where(a0, prev[p0] + odd[p0], NEGF)
        c1 = np.where(a1, prev[p1] + odd[p1], NEGF)
        D[t] = np.where(alive, np.maximum(c0, c1), NEGF)
    valid = D[m] > NEGF / 2
    with np.errstate(invalid='ignore', divide='ignore'):
        ts = np.arange(m)[:, None]
        ratio = (D[m][None, :] - D[:m]) / (m - ts)
        ratio[D[:m] <= NEGF / 2] = np.inf
        per_v = ratio.min(axis=0)
    per_v[~valid] = -np.inf
    est = float(per_v.max())
    cand = Fraction(est).limit_denominator(m)
    ok1 = has_expanding_cycle(k, Thr('gt', cand.numerator, cand.denominator), R)[0] is None
    ok2 = cand == 0 or has_expanding_cycle(k, Thr('ge', cand.numerator, cand.denominator), R)[0] is not None
    if ok1 and ok2:
        return cand
    return max_density_karp(k, alive)


def step_function(k, log=None, odd_only=False):
    """The exact step function c -> FVS_c(k) on [0, 1), as a list of segments from the top:
         dict(top=b_i, bottom=b_{i+1}, V=V_i, R=R_i, cert=cycles_i)
       meaning FVS_c(k) = V_i for c in [b_{i+1}, b_i), certified by
         * R_i: |R_i| = V_i and rho_max(B(2,k) - R_i) = b_{i+1}   (so FVS_c <= V_i for c >= b_{i+1});
         * cert_i: cycles of density >= b_i whose minimum hitting set is V_i  (so FVS_c >= V_i for c < b_i).
       Method: from the current value V, c*(V+1) = min over |R| = V of rho_max(B - R) is found by descent:
       rho = rho_max(B - R) (Karp); FVS(>= rho) is computed by IHS; if it is still V the new optimal R has
       smaller rho_max, otherwise rho is the next breakpoint."""
    allowed = [v for v in range(1 << k) if (v & 1 or not odd_only)]
    top = Fraction(1, 1)
    V = 1
    R = [(1 << k) - 1]
    cert = [[(1 << k) - 1]]
    segs = []
    while True:
        rho = rho_max_removed(k, R)
        r = None
        while rho is not None and rho > 0:
            thr = Thr('ge', rho.numerator, rho.denominator)
            r = ihs_fvs(k, thr, odd_only=odd_only)
            if r['value'] > V:
                break
            R = r['R']
            rho2 = rho_max_removed(k, R)
            check(rho2 is None or rho2 < rho, "descent did not decrease rho_max")
            rho = rho2
            r = None
        bottom = Fraction(0) if (rho is None or rho == 0) else rho
        segs.append(dict(top=top, bottom=bottom, V=V, R=list(R), cert=cert))
        if log:
            log("   k=%d: FVS_c = %d on [%s, %s)" % (k, V, bottom, top))
        if bottom == 0:
            break
        top = bottom
        V = r['value']
        R = r['R']
        cert = shrink_certificate(r['cycles'], allowed, V)
    return segs


def cycle_to_word(cyc):
    """compact encoding of a cycle of B(2,k): (length p, integer sum_t z_t 2^t), z_t = first letter of node t"""
    return [len(cyc), sum((x & 1) << t for t, x in enumerate(cyc))]


def word_to_cycle(enc, k):
    p, u = enc
    return cycle_from_word(u, p, k)


def rc2_min_cover(cycles, allowed):
    """independent exact minimum hitting set (pysat RC2 MaxSAT, core-guided): returns (value, set)"""
    from pysat.examples.rc2 import RC2
    from pysat.formula import WCNF
    idx = {v: i + 1 for i, v in enumerate(allowed)}
    wc = WCNF()
    for c in cycles:
        cl = sorted(set(idx[v] for v in c if v in idx))
        if not cl:
            return None, None
        wc.append(cl)
    used = sorted(set(v for c in cycles for v in c if v in idx))
    for v in used:
        wc.append([-idx[v]], weight=1)
    with RC2(wc) as rc:
        m = rc.compute()
        if m is None:
            return None, None
        R = [allowed[i - 1] for i in m if i > 0 and allowed[i - 1] in set(used)]
        return rc.cost, R


def _rc2_worker(cycles, allowed, q):
    v, _ = rc2_min_cover(cycles, allowed)
    q.put(v)


def rc2_min_cover_timed(cycles, allowed, timeout=60):
    """rc2_min_cover in a forked subprocess with a wall-clock limit; returns the optimum or None on time-out"""
    import multiprocessing as mp
    ctx = mp.get_context('fork')
    q = ctx.Queue()
    pr = ctx.Process(target=_rc2_worker, args=(cycles, allowed, q))
    pr.start()
    pr.join(timeout)
    if pr.is_alive():
        pr.terminate()
        pr.join()
        return None
    return q.get() if not q.empty() else None


def shrink_certificate(cycles, allowed, value, batch=50, time_limit=600):
    """a (usually much) smaller subset S of `cycles` with the same minimum hitting-set value: IHS restricted to
    the finite list (start from the shortest cycles; add violated ones in batches).  If this takes longer than
    time_limit seconds, the full list is returned (it is a valid certificate by construction)."""
    t0 = time.time()
    order = sorted(range(len(cycles)), key=lambda i: len(cycles[i]))
    S = [cycles[i] for i in order[:batch]]
    inS = set(order[:batch])
    while True:
        if time.time() - t0 > time_limit:
            return list(cycles)
        v, R = solve_cover_highs(S, allowed)
        Rs = set(R)
        viol = [i for i in order if i not in inS and not any(x in Rs for x in cycles[i])]
        if not viol:
            check(v == value, "shrunk certificate lost the bound (%s vs %s)" % (v, value))
            return S
        for i in viol[:batch]:
            S.append(cycles[i])
            inS.add(i)


def verify_fvs_value(k, thr, R, lb_cycles, value, odd_only=False, packing=None, rc2_timeout=60):
    """exact re-verification of FVS = value:
       (1) |R| = value, R inside the allowed set, B(2,k) - R has no expanding cycle (Bellman-Ford) and the
           returned potential is re-checked edge by edge;
       (2) lower bound: either a PACKING of `value` pairwise node-disjoint expanding cycles, or a list of
           expanding cycles whose minimum hitting set is re-proved = value by RC2 (independent of the HiGHS search)."""
    n = 1 << k
    check(len(R) == value and len(set(R)) == value, "UB set size")
    if odd_only:
        check(all(v & 1 for v in R), "odd-only set has an even node")
    cyc, d = has_expanding_cycle(k, thr, R)
    check(cyc is None, "UB set leaves an expanding cycle")
    check(verify_potential(k, thr, R, d), "potential certificate fails")
    allowed = [v for v in range(n) if (v & 1 or not odd_only)]
    if packing is not None:
        check(len(packing) == value, "packing size")
        used = set()
        for c in packing:
            check(cycle_ok(c, k), "packing cycle not simple")
            check(thr.expanding(sum(x & 1 for x in c), len(c)), "packing cycle not expanding")
            check(not (used & set(c)), "packing cycles not disjoint")
            used |= set(c)
        return 'packing'
    for c in lb_cycles:
        check(cycle_ok(c, k), "certificate cycle is not a simple cycle")
        check(thr.expanding(sum(x & 1 for x in c), len(c)), "certificate cycle not expanding")
    v = rc2_min_cover_timed(lb_cycles, allowed, timeout=rc2_timeout)
    if v is not None:
        check(v == value, "RC2 lower bound %s != %s" % (v, value))
        return 'rc2'
    v, _ = solve_cover_highs(lb_cycles, allowed)
    check(v == value, "HiGHS lower bound %s != %s" % (v, value))
    return 'highs'


# ============================================================================ fractional packings (LP duals)
def lp_dual_packing(cycles, allowed, denom=10 ** 6):
    """LP relaxation of the hitting set over `cycles` restricted to `allowed` nodes; returns an exact rational
    fractional packing y (list of Fractions, one per cycle) with sum_{C contains v} y_C <= 1 for every allowed
    node v, obtained by rounding the HiGHS row duals down and rescaling; and its value sum y."""
    import highspy
    pos = {v: i for i, v in enumerate(allowed)}
    h = highspy.Highs()
    h.silent()
    nv = len(allowed)
    h.addVars(nv, np.zeros(nv), np.full(nv, highspy.kHighsInf))
    h.changeColsCost(nv, np.arange(nv, dtype=np.int32), np.ones(nv))
    starts, idxs = [], []
    rows = []
    for c in cycles:
        cols = sorted(set(pos[v] for v in c if v in pos))
        rows.append(cols)
        starts.append(len(idxs))
        idxs.extend(cols)
    nr = len(starts)
    h.addRows(nr, np.ones(nr), np.full(nr, highspy.kHighsInf), len(idxs),
              np.array(starts, dtype=np.int32), np.array(idxs, dtype=np.int32), np.ones(len(idxs)))
    h.setMinimize()
    h.run()
    duals = h.getSolution().row_dual
    y = [Fraction(max(0, int(np.floor(max(0.0, d) * denom))), denom) for d in duals]
    load = {}
    for yc, cols in zip(y, rows):
        if yc:
            for cidx in cols:
                load[cidx] = load.get(cidx, 0) + yc
    mx = max(load.values()) if load else Fraction(1)
    if mx > 1:
        y = [yc / mx for yc in y]
    return y, sum(y), h.getInfo().objective_function_value


def verify_frac_packing(cycles, y, allowed):
    """exact check: y >= 0 and every allowed node carries total weight <= 1; returns sum y (then every
    allowed hitting set of the cycles has size >= sum y, by double counting)"""
    aset = set(allowed)
    load = {}
    for c, yc in zip(cycles, y):
        check(yc >= 0, "negative packing weight")
        for v in set(c):
            if v in aset:
                load[v] = load.get(v, 0) + yc
    check(all(val <= 1 for val in load.values()), "fractional packing overloads a node")
    return sum(y)


# ============================================================================ fractional covers (tau*)
def _bf_longest_from(k, alive_mask, w, src_nodes_w, targets_ok):
    """longest walks from a virtual source: d[y] = max over (p, walk p -> z_1..z_m -> y, z_i in alive_mask)
    handled by caller; this helper runs Bellman-Ford from one start node p through the zero-cost subgraph."""
    raise NotImplementedError


def cheap_cycle_below(k, thr, x, D):
    """Exact separation for a fractional cover x (dict node -> nonnegative integer UNITS, x_v = units/D):
    find an expanding cycle C of B(2,k) with sum_{v in C} units(v) <= D-1 (i.e. x(C) < 1), or return None.
    Method: P = support of x, Zs = the other nodes.  (0) a positive cycle inside B[Zs] is returned directly;
    otherwise longest walks through Zs between P-nodes (Bellman-Ford, finite since B[Zs] has no positive cycle)
    give a matrix A[p][p'], and a max-plus DP over (P-node, units used <= D-1) finds a positive closed P-walk.
    The corresponding closed walk of B(2,k) has positive weight and x-weight < 1; its densest simple sub-cycle is
    returned (it is expanding and has x-weight < 1)."""
    n = 1 << k
    w, _ = node_weights(k, thr)
    units = np.zeros(n, dtype=np.int64)
    for v, u in x.items():
        units[v] = u
    P = [int(v) for v in np.nonzero(units > 0)[0]]
    Zmask = units == 0
    cyc, _ = positive_cycle(k, Zmask.copy(), w)
    if cyc is not None:
        return cyc
    p0, p1, _ = db_arrays(k)
    X = np.arange(n, dtype=np.int64)
    s0 = X >> 1
    s1 = s0 | (1 << (k - 1))
    NEGI = -(1 << 60)
    Pidx = {p: i for i, p in enumerate(P)}
    m = len(P)
    A = np.full((m, m), NEGI, dtype=np.int64)
    parents = {}
    for p in P:
        # d[y]: best weight of walk p -> (Zs)* -> y, weight counts nodes left (p and the Zs nodes)
        d = np.full(n, NEGI, dtype=np.int64)
        par = np.full(n, -1, dtype=np.int64)
        for y in (int(s0[p]), int(s1[p])):
            if w[p] > d[y]:
                d[y] = int(w[p])
                par[y] = p
        for it in range(n + 2):
            # relax edges out of Zs nodes
            src = np.nonzero(Zmask & (d > NEGI // 2))[0]
            changed = False
            for b_arr in (s0, s1):
                tgt = b_arr[src]
                cand = d[src] + w[src]
                # process maxima per target
                order = np.argsort(-cand, kind='stable')
                tgt_o = tgt[order]
                cand_o = cand[order]
                src_o = src[order]
                seen_t = set()
                for t_, c_, s_ in zip(tgt_o.tolist(), cand_o.tolist(), src_o.tolist()):
                    if t_ in seen_t:
                        continue
                    seen_t.add(t_)
                    if c_ > d[t_]:
                        d[t_] = c_
                        par[t_] = s_
                        changed = True
            if not changed:
                break
            check(it < n + 1, "positive cycle in zero-cost subgraph (should have been caught)")
        parents[p] = par.copy()
        for q2 in P:
            A[Pidx[p], Pidx[q2]] = d[q2]
    D1 = D - 1
    for s_i in range(m):
        best = np.full((m, D1 + 1), NEGI, dtype=np.int64)
        pj = np.full((m, D1 + 1), -1, dtype=np.int64)
        pb = np.full((m, D1 + 1), -1, dtype=np.int64)
        best[s_i, 0] = 0
        for b in range(D1 + 1):
            for j in np.nonzero(best[:, b] > NEGI // 2)[0].tolist():
                nb = b + int(units[P[j]])
                if nb > D1:
                    continue
                vals = np.where(A[j] > NEGI // 2, best[j, b] + A[j], NEGI)
                if vals[s_i] > 0:
                    seq = [s_i, j]
                    cj, cb = j, b
                    while pj[cj, cb] >= 0:
                        cj, cb = int(pj[cj, cb]), int(pb[cj, cb])
                        seq.append(cj)
                    seq.reverse()
                    walk = []
                    for a_, b_ in zip(seq[:-1], seq[1:]):
                        pa, pbn = P[a_], P[b_]
                        par = parents[pa]
                        path = []
                        y = pbn
                        while True:
                            y = int(par[y])
                            path.append(y)
                            if y == pa:
                                break
                        walk.extend(reversed(path))
                    c = _densest_simple_subcycle(walk)
                    check(cycle_ok(c, k), "separation produced a non-cycle")
                    return c
                imp = vals > best[:, nb]
                if imp.any():
                    best[imp, nb] = vals[imp]
                    pj[imp, nb] = j
                    pb[imp, nb] = b
    return None


def tau_star(k, thr, cycles, max_rounds=500, log=None, denom=720):
    """cutting-plane computation of a fractional cover: LP over the known cycles, rounded UP to multiples of
    1/denom (so it stays feasible for the known rows), then exact separation (cheap_cycle_below).
    Returns (value as Fraction, units dict, denom, cycles).  The returned x is a certified fractional cover of
    ALL expanding cycles (separation found nothing), so nu_c(k) <= nu*_c(k) = tau*_c(k) <= value."""
    import highspy
    n = 1 << k
    cycles = [list(c) for c in cycles]
    for rnd in range(max_rounds):
        h = highspy.Highs()
        h.silent()
        h.addVars(n, np.zeros(n), np.ones(n))
        h.changeColsCost(n, np.arange(n, dtype=np.int32), np.ones(n))
        starts, idxs = [], []
        for c in cycles:
            starts.append(len(idxs))
            idxs.extend(sorted(set(c)))
        nr = len(starts)
        h.addRows(nr, np.ones(nr), np.full(nr, highspy.kHighsInf), len(idxs),
                  np.array(starts, dtype=np.int32), np.array(idxs, dtype=np.int32), np.ones(len(idxs)))
        h.setMinimize()
        h.run()
        xv = h.getSolution().col_value
        lpval = h.getInfo().objective_function_value
        best_choice = None
        for D in list(range(1, 61)) + [denom]:
            units = {}
            for v in range(n):
                if xv[v] > 1e-9:
                    units[v] = int(np.ceil(xv[v] * D - 1e-6))
            for c in cycles:
                tot = sum(units.get(v, 0) for v in set(c))
                if tot < D:
                    v0 = max(set(c), key=lambda v: units.get(v, 0))
                    units[v0] = units.get(v0, 0) + D - tot
            val = Fraction(sum(units.values()), D)
            if best_choice is None or val < best_choice[0] - Fraction(1, 10 ** 9):
                best_choice = (val, units, D)
            if float(val) <= lpval + 1e-7:
                break
        val, units, D = best_choice
        new = cheap_cycle_below(k, thr, units, D)
        if log:
            log("   tau* round %d: LP %.4f  rounded cover %s = %.4f  %s" % (
                rnd, h.getInfo().objective_function_value, val, float(val), 'violated' if new else 'CERTIFIED'))
        if new is None:
            return val, units, D, cycles
        cycles.append(new)
    raise RuntimeError("tau_star did not converge")


# ============================================================================ Terras coordinates (q-maps)
def Tq(n, q):
    return n // 2 if n % 2 == 0 else (q * n + 1) // 2


def parity_word_int(r, k, q):
    """first k parities of T_q starting from r, encoded LSB-first (the de Bruijn node of residue r)"""
    X = 0
    x = r
    for j in range(k):
        X |= (x & 1) << j
        x = Tq(x, q)
    return X


def terras_tables(k, q):
    """res2word[r] = de Bruijn node of residue r mod 2^k; word2res inverse (bijection by Terras)"""
    n = 1 << k
    res2word = np.array([parity_word_int(r, k, q) for r in range(n)], dtype=np.int64)
    check(len(set(res2word.tolist())) == n, "Terras map not a bijection (k=%d, q=%d)" % (k, q))
    word2res = np.zeros(n, dtype=np.int64)
    word2res[res2word] = np.arange(n)
    return res2word, word2res


def edit_map(q, k, R_words):
    """the periodic edit G of T_q at level k: G(n) = 1 if the parity word of n mod 2^k lies in R_words,
    else T_q(n).  Returns (G, set of edited residues)."""
    res2word, word2res = terras_tables(k, q)
    Rres = set(int(word2res[x]) for x in R_words)
    M = (1 << k) - 1

    def G(n):
        return 1 if (n & M) in Rres else Tq(n, q)
    return G, Rres


def ballot_walk_length(q, k, R_words):
    """Lookahead certificate for the periodic edit G (= 1 on the classes of R_words, T_q elsewhere).
    A BALLOT WALK of length j is a sequence of nodes x_0 .. x_{j-1} of H = B(2,k) - R (consecutive nodes joined
    by edges) such that the prefix multipliers q^(a_i)/2^i exceed 1 for i = 1..j (a_i = number of odd nodes among
    x_0..x_{i-1}).  Returns (L, n0, counts):
      L  = maximal length of a ballot walk (finite iff H has no expanding cycle; +infinity otherwise);
      n0 = max over walks x_0..x_j (x_0..x_{j-1} ballot, x_j in H, multiplier after j+1 steps < 1) of
           floor(c/(2^(j+1) - q^(a_(j+1)))), where T_q^(j+1)(n) = (q^a n + c)/2^(j+1) on the class of the walk;
      counts[j-1] = number of ballot walks of length j.
    Consequence (note, Theorem 3): every n > max(n0, 1) satisfies G^i(n) < n for some 1 <= i <= L + 1."""
    n = 1 << k
    alive = np.ones(n, dtype=bool)
    if len(R_words):
        alive[np.asarray(list(R_words), dtype=np.int64)] = False
    X = np.arange(n)
    succs = [(int(x) >> 1, (int(x) >> 1) | (1 << (k - 1))) for x in range(n)]
    odd = [int(x) & 1 for x in range(n)]
    # state: dict (x, a) -> (count, max c) for walks of length j whose NEXT node is x (x in H)
    cur = {}
    for x in range(n):
        if alive[x]:
            cur[(x, 0)] = (1, 0)
    j = 0
    counts = []
    n0 = 0
    while True:
        nxt = {}
        for (x, a), (cnt, cmax) in cur.items():
            par = odd[x]
            a2 = a + par
            c2 = q * cmax + (1 << j) if par else cmax
            if q ** a2 > 2 ** (j + 1):
                for y in succs[x]:
                    old = nxt.get((y, a2))
                    if old is None:
                        nxt[(y, a2)] = (cnt, c2)
                    else:
                        nxt[(y, a2)] = (old[0] + cnt, max(old[1], c2))
            else:
                # first descent at step j+1 for n large: threshold c2 / (2^(j+1) - q^a2)
                n0 = max(n0, c2 // (2 ** (j + 1) - q ** a2))
        j += 1
        tot = sum(v[0] for v in nxt.values())
        if tot == 0:
            return j - 1, n0, counts
        counts.append(tot)
        cur = {key: v for key, v in nxt.items() if alive[key[0]]}
        if j > 5000:
            raise RuntimeError("ballot walks too long (expanding cycle left?)")


def orbit_census(G, nmax, cap=10 ** 6):
    """iterate G from every 1 <= n <= nmax; returns (dict cycle(min elt) -> cycle list, max steps to cycle,
    number of starts entering each cycle).  Uses memoized cycle membership."""
    fate = {}
    cycles = {}
    maxsteps = 0
    for n0 in range(1, nmax + 1):
        path = []
        seen = {}
        x = n0
        steps = 0
        while x not in fate and x not in seen:
            seen[x] = len(path)
            path.append(x)
            x = G(x)
            steps += 1
            if steps > cap:
                raise RuntimeError("orbit of %d too long" % n0)
        if x in fate:
            f = fate[x]
        else:
            cyc = path[seen[x]:]
            f = min(cyc)
            cycles[f] = cyc
            for y in cyc:
                fate[y] = f
        for y in path:
            fate[y] = f
        maxsteps = max(maxsteps, steps)
    hits = {}
    for n0 in range(1, nmax + 1):
        hits[fate[n0]] = hits.get(fate[n0], 0) + 1
    return cycles, maxsteps, hits


def first_descent_time(G, n, L):
    """least 1 <= j <= L with G^j(n) < n, or None"""
    x = n
    for j in range(1, L + 1):
        x = G(x)
        if x < n:
            return j
    return None


def rising_witness(q, k, R_words, cyc, L):
    """for an expanding cycle `cyc` of B(2,k) avoiding R_words: the least positive integer n whose first
    k+L parities follow the cycle from its ballot rotation; then G^j(n) = T_q^j(n) > n for 1 <= j <= L.
    Returns (n, start index)."""
    z = [x & 1 for x in cyc]
    p = len(z)
    # ballot start: after the last minimum of the height walk (weights a*log q - j*log 2, exact by integers)
    best = None
    for t0 in range(p):
        ok = True
        a = 0
        for j in range(1, p + 1):
            a += z[(t0 + j - 1) % p]
            if not q ** a > 2 ** j:
                ok = False
                break
        if ok:
            best = t0
            break
    check(best is not None, "no ballot rotation of an expanding cycle (impossible by the cycle lemma)")
    t0 = best
    bits = [z[(t0 + i) % p] for i in range(k + L)]
    # find n mod 2^(k+L) with these parities (Terras inverse, digit by digit)
    n = 0
    for m in range(1, k + L + 1):
        # choose bit m-1 of n so that the first m parities match
        for cand in (n, n + (1 << (m - 1))):
            x = cand
            ok = True
            for i in range(m):
                if (x & 1) != bits[i]:
                    ok = False
                    break
                x = Tq(x, q)
            if ok:
                n = cand
                break
        else:
            raise RuntimeError("Terras inverse failed")
    if n == 0:
        n = 1 << (k + L)
    return n, t0


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


def pack_json(obj):
    return base64.b64encode(gzip.compress(json.dumps(obj, separators=(',', ':')).encode())).decode()


def unpack_json(s):
    return json.loads(gzip.decompress(base64.b64decode(s)).decode())
