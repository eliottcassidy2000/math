#!/usr/bin/env python3
"""
procgen_seven3_20260926_lib.py -- helpers of lane "seven3" (session collatz-procgen-20260922, 2026-09-26):
itinerary-coded sign strategies of 7n +- 1 and their exact rho_max via Markov refinement.

Setting (THM-4486, THM-4508).  q = 7 (the code keeps q as a parameter).  A sign rule sigma gives every odd 2-adic
integer a sign; T(x) = x/2 (x even), (q x + sigma(x))/2 (x odd).  rho_max(sigma) = the largest odd density of a
periodic orbit of T (Lemma C of THM-4508: = the densest cycle of the parity graph at any level >= the depth of sigma).

MH (max-halving): the sign s with q x + s = 0 mod 4.  MH itinerary of odd x: (s_1, v_1), (s_2, v_2), ... with
x_i = (q x_(i-1) + s_i) / 2^(v_i); a prefix of n symbols is one residue class mod 2^(1 + sum v_i) (Lemma MH(ii)).

A RULE is a finite partition of the odd 2-adic integers into residue classes (c mod 2^d) with one sign each; it is
stored as a dict {(c, d): sign}.  Rules are generated from MH-itinerary cylinders by word_class().

MARKOV REFINEMENT (the exact evaluator of this lane).  Leaves of a binary trie (classes (c, d), both parities) such
that sigma is constant on each odd leaf and the T-image of every leaf (a class of depth d - 1) is a union of leaves.
Then (Lemma M of the note) the closed walks of the leaf graph L -> L' (L' inside T(L)) are in bijection with the
periodic orbits of T, with the same odd density; rho_max = the maximum mean cycle of the leaf graph (weight 1 on odd
leaves).  Exactness: a witness cycle is turned into its rational periodic point and re-walked with exact
arithmetic; the upper bound is an integer potential checked on every edge (Lemma G1 on the leaf graph).
"""
import os
import sys
import hashlib
import subprocess
from fractions import Fraction as Fr
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SCR = os.path.join(ROOT, 'scratch', 'procgen_seven3')
Q = 7


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ MH arithmetic
def mh(x, q=Q):
    """MH sign of odd x (int): the s in {+1,-1} with q x + s = 0 mod 4"""
    return 1 if (q * x + 1) % 4 == 0 else -1


def v2(n):
    return (n & -n).bit_length() - 1


def itin_of_residue(c, d, q=Q, maxsym=10 ** 9):
    """the MH itinerary symbols of the class c mod 2^d that are determined by the class, plus the undetermined tail:
    returns (symbols [(s, v)], (s_last, vmin) or None).  The last entry (s, vmin) means: next sign s, valuation >= vmin
    (not determined).  None means the next sign itself is undetermined (fewer than 2 bits left)."""
    out = []
    x, m = c % (1 << d), d          # x known mod 2^m
    while len(out) < maxsym:
        if m < 2 or x % 2 == 0:
            return out, None
        s = mh(x, q)
        t = (q * x + s) % (1 << m)
        if t == 0:
            return out, (s, m)
        v = v2(t)
        if v + 1 > m:
            return out, (s, m)
        out.append((s, v))
        x, m = (t >> v), m - v
    return out, None


def word_class(word, q=Q):
    """the residue class (c, d) of the MH-itinerary cylinder of `word`: a list of (s, v) with v an int >= 2 (exact
    valuation) or ('ge', V) (valuation >= V, allowed only as the last symbol).  Lemma MH(ii)."""
    qi = None
    if word and isinstance(word[-1][1], tuple):
        s, (_, V) = word[-1]
        d = V
        c = (-s * pow(q, -1, 1 << d)) % (1 << d)
        rest = word[:-1]
    else:
        c, d = 1, 1
        rest = word
    for (s, v) in reversed(rest):
        assert isinstance(v, int) and v >= 2
        d2 = d + v
        c = ((c << v) - s) * pow(q, -1, 1 << d2) % (1 << d2)
        d = d2
    # sanity: the class has the requested itinerary
    return c, d


def check_word_class(word, q=Q):
    c, d = word_class(word, q)
    sy, tail = itin_of_residue(c, d, q)
    exact = [w for w in word if isinstance(w[1], int)]
    assert sy[:len(exact)] == exact, (word, sy)
    if word and isinstance(word[-1][1], tuple):
        s, (_, V) = word[-1]
        assert len(sy) == len(exact) and tail is not None and tail[0] == s and tail[1] >= V, (word, sy, tail)
    return c, d


# ------------------------------------------------------------------------------------------------ rules
def mh_rule(depth=2, q=Q):
    """max-halving as a rule on classes mod 2^depth"""
    return {(c, depth): mh(c, q) for c in range(1, 1 << depth, 2)}


def rule_sign(rule, x, dmax):
    """sign of the odd integer x under the rule (search the leaf containing x)"""
    for d in range(1, dmax + 1):
        key = (x % (1 << d), d)
        if key in rule:
            return rule[key]
    raise KeyError(x)


def check_partition(rule):
    """the odd leaves partition the odd 2-adic integers (measure sums to 1/2, no leaf inside another)"""
    tot = Fr(0)
    ds = sorted({d for (_, d) in rule})
    for (c, d) in rule:
        assert c % 2 == 1 and 0 <= c < (1 << d)
        for e in ds:
            if e < d and (c % (1 << e), e) in rule:
                raise AssertionError('nested leaves %s' % ((c, d),))
        tot += Fr(1, 1 << d)
    assert tot == Fr(1, 2), tot
    return True


def flip_rule(flips, q=Q, base_depth=2):
    """MH except on the given odd classes [(c, d)] (pairwise disjoint), where the sign is flipped.  Returns the
    rule (a partition with signs): the flipped classes plus the complement cut into classes."""
    flips = sorted(set((c % (1 << d), d) for (c, d) in flips), key=lambda t: t[1])
    # complement: descend the binary trie from the odd class (1,1)
    fl = set(flips)
    depthset = {}
    for (c, d) in flips:
        for e in range(1, d):
            depthset[(c % (1 << e), e)] = True    # internal nodes on the path
    rule = {}
    stack = [(1, 1)]
    while stack:
        c, d = stack.pop()
        if (c, d) in fl:
            rule[(c, d)] = -mh(c, q) if d >= 2 else None
            if d < 2:
                raise ValueError('flip class too shallow for a sign')
            continue
        if (c, d) in depthset or d < base_depth:
            stack.append((c, d + 1))
            stack.append((c + (1 << d), d + 1))
            continue
        rule[(c, d)] = mh(c, q)
    return rule


# ------------------------------------------------------------------------------------------------ Markov refinement
class Markov:
    """leaf trie refined until images are unions of leaves (see module docstring)"""

    def __init__(self, rule, q=Q):
        self.q = q
        self.leaf = {}        # (c, d) -> sign (odd) or 0 (even)
        self.internal = set()
        for (c, d), s in rule.items():
            self.leaf[(c, d)] = s
        # internal nodes on the paths to the odd leaves
        for (c, d) in rule:
            for e in range(0, d):
                self.internal.add((c % (1 << e), e))
        # even part: leaf (0, 1) unless something is under it (never, rules are odd)
        self.leaf[(0, 1)] = 0
        self.internal.add((0, 0))
        self.dmax = max(d for (_, d) in self.leaf)
        self._refine()

    def image(self, c, d):
        s = self.leaf[(c, d)]
        if c % 2 == 0:
            t = c >> 1
        else:
            t = (self.q * c + s) >> 1
        return t % (1 << (d - 1)), d - 1

    def _split(self, c, d):
        s = self.leaf.pop((c, d))
        self.internal.add((c, d))
        a, b = (c, d + 1), (c + (1 << d), d + 1)
        self.leaf[a] = s if a[0] % 2 else 0
        self.leaf[b] = s if b[0] % 2 else 0
        return a, b

    def _ensure_union(self, c, d, work):
        """make the class (c, d) a union of leaves (split the leaf containing it, along the path)"""
        for e in range(0, d + 1):
            key = (c % (1 << e), e)
            if key in self.leaf:
                if e == d:
                    return
                # split down to depth d along the path
                cur = key
                while cur[1] < d:
                    a, b = self._split(*cur)
                    work.extend([a, b])
                    nxt = (c % (1 << (cur[1] + 1)), cur[1] + 1)
                    cur = nxt
                return
            if key not in self.internal:
                raise AssertionError('broken trie at %s' % (key,))
        return   # (c, d) is an internal node: a union of leaves

    def _refine(self):
        work = list(self.leaf.keys())
        while work:
            L = work.pop()
            if L not in self.leaf:
                continue
            c, d = L
            if d == 0:
                continue
            ic, idd = self.image(c, d)
            self._ensure_union(ic, idd, work)

    def leaves_under(self, c, d):
        out, st = [], [(c, d)]
        while st:
            k = st.pop()
            if k in self.leaf:
                out.append(k)
            else:
                cc, dd = k
                st.append((cc, dd + 1))
                st.append((cc + (1 << dd), dd + 1))
        return out

    def graph(self):
        nodes = sorted(self.leaf.keys(), key=lambda t: (t[1], t[0]))
        idx = {n: i for i, n in enumerate(nodes)}
        succ = []
        for (c, d) in nodes:
            ic, idd = self.image(c, d)
            succ.append([idx[k] for k in self.leaves_under(ic, idd)])
        w = [c & 1 for (c, d) in nodes]
        return nodes, w, succ


# ------------------------------------------------------------------------------------------------ max mean cycle
def howard(w, succ, iters=10000):
    """Howard policy iteration for the maximum mean cycle (weights w on nodes).  Floating point search; the caller
    certifies the result exactly (certify_upper + cycle check).  Returns (best cycle as node list)."""
    n = len(w)
    pol = [s[0] for s in succ]
    for s_ in range(n):   # prefer an odd successor initially
        for t in succ[s_]:
            if w[t]:
                pol[s_] = t
                break
    for it in range(iters):
        # value determination on the functional graph
        eta = [None] * n
        h = [0.0] * n
        state = [0] * n     # 0 new, 1 on stack, 2 done
        for s0 in range(n):
            if state[s0]:
                continue
            path = []
            u = s0
            while state[u] == 0:
                state[u] = 1
                path.append(u)
                u = pol[u]
            if state[u] == 1:
                # new cycle starting at u
                j = path.index(u)
                cyc = path[j:]
                m = sum(w[x] for x in cyc) / len(cyc)
                # h on the cycle: root u, h(u)=0, h(v) = w(v) - m + h(pol v)
                h[u] = 0.0
                eta[u] = m
                for x in reversed(cyc[1:]):
                    h[x] = w[x] - m + h[pol[x]]
                    eta[x] = m
                for x in cyc:
                    state[x] = 2
                path = path[:j]
            for x in reversed(path):
                eta[x] = eta[pol[x]]
                h[x] = w[x] - eta[x] + h[pol[x]]
                state[x] = 2
        # improvement
        changed = False
        eps = 1e-9
        for s_ in range(n):
            best_eta = eta[pol[s_]]
            bt = pol[s_]
            for t in succ[s_]:
                if eta[t] > best_eta + eps:
                    best_eta, bt = eta[t], t
            if bt != pol[s_]:
                pol[s_] = bt
                changed = True
        if changed:
            continue
        for s_ in range(n):
            cur = h[pol[s_]]
            bt = pol[s_]
            for t in succ[s_]:
                if abs(eta[t] - eta[s_]) < eps and h[t] > cur + eps:
                    cur, bt = h[t], t
            if bt != pol[s_]:
                pol[s_] = bt
                changed = True
        if not changed:
            break
    # extract the best cycle of the final policy
    best, bestcyc = -1.0, None
    seen = [False] * n
    for s0 in range(n):
        if seen[s0]:
            continue
        path, u, onp = [], s0, {}
        while not seen[u] and u not in onp:
            onp[u] = len(path)
            path.append(u)
            u = pol[u]
        if u in onp:
            cyc = path[onp[u]:]
            m = sum(w[x] for x in cyc) / len(cyc)
            if m > best:
                best, bestcyc = m, cyc
        for x in path:
            seen[x] = True
    return bestcyc


def certify_upper(w, succ, a, p, maxsweep=10 ** 7):
    """exact integer potential psi with psi(t) + e(s) <= psi(s) on every edge, e = p - a (odd), -a (even), by the
    least fixed point iteration (Lemma G1).  Returns psi (list) or None if the iteration exceeds the bound (then a
    denser cycle exists)."""
    n = len(w)
    e = [(p - a) if w[i] else -a for i in range(n)]
    bound = n * max(p - a, 1) + 1
    psi = [0] * n
    # reverse adjacency for a worklist
    pred = [[] for _ in range(n)]
    for s_ in range(n):
        for t in succ[s_]:
            pred[t].append(s_)
    from collections import deque
    dq = deque(range(n))
    inq = [True] * n
    while dq:
        s_ = dq.popleft()
        inq[s_] = False
        val = e[s_] + max(psi[t] for t in succ[s_])
        if val > psi[s_]:
            psi[s_] = val
            if val > bound:
                return None
            for r in pred[s_]:
                if not inq[r]:
                    inq[r] = True
                    dq.append(r)
    # exact check
    for s_ in range(n):
        for t in succ[s_]:
            assert psi[t] + e[s_] <= psi[s_]
    return psi


# ------------------------------------------------------------------------------------------------ cycles -> rationals
def cycle_rational(nodes_cyc, signs, q=Q):
    """the periodic point of the walk through classes nodes_cyc = [(c, d)] with signs (0 for even): the fixed point of
    the composed affine maps; then an exact re-walk: x_i in class i, parity and sign as given.  Returns (x0, odd, len)."""
    n = len(nodes_cyc)
    a = 0
    num = Fr(1)
    const = Fr(0)
    for (c, d), s in zip(nodes_cyc, signs):
        if c % 2 == 0:
            num, const = num / 2, const / 2
        else:
            num, const = q * num / 2, (q * const + s) / 2
            a += 1
    # x = num x + const  => x = const / (1 - num)
    x0 = const / (1 - num)
    x = x0
    for (c, d), s in zip(nodes_cyc, signs):
        D = x.denominator
        assert D % 2 == 1
        r = x.numerator * pow(D, -1, 1 << d) % (1 << d)
        assert r == c, ('class mismatch', x, c, d)
        if c % 2 == 0:
            x = x / 2
        else:
            x = (q * x + s) / 2
    assert x == x0
    return x0, a, n


# ------------------------------------------------------------------------------------------------ uniform-level cross-check
def build_engines():
    """compile the read-only seven2 engine rhomax.c (and this lane's C engines) into SCR/bin"""
    import shutil
    b = os.path.join(SCR, 'bin')
    os.makedirs(b, exist_ok=True)
    out = {}
    for name, src in [('rhomax', 'procgen_seven2_20260926_rhomax.c'), ('game', 'procgen_seven3_20260926_game.c'),
                      ('markov', 'procgen_seven3_20260926_markov.c'), ('adv', 'procgen_seven3_20260926_adv.c')]:
        p = os.path.join(HERE, src)
        exe = os.path.join(b, name + '_' + sha(p)[:12])
        if not os.path.exists(exe):
            subprocess.run(['cc', '-O2', '-o', exe, p], check=True)
        out[name] = exe
    return out


def rule_minus_array(rule, k):
    """the level-k strategy of a rule of depth <= k, as the uint8 array over odd nodes 2i+1 (1 = sign '-')"""
    import numpy as np
    H = 1 << (k - 1)
    arr = np.full(H, 255, dtype=np.uint8)
    for (c, d), s in rule.items():
        if d > k:
            raise ValueError('rule deeper than k')
        arr[(c - 1) // 2::(1 << (d - 1))] = 1 if s < 0 else 0
    assert not (arr == 255).any()
    return arr


def rhomax_uniform(rule, k, tag='u', q=Q):
    """exact rho_max of the rule at the uniform level k by the seven2 engine (independent of the Markov code):
    returns (Fraction, witness cycle of G_sigma as residues mod 2^k)"""
    import numpy as np
    eng = build_engines()['rhomax']
    arr = rule_minus_array(rule, k)
    bits = np.packbits(arr.astype(np.uint8), bitorder='little')
    path = os.path.join(SCR, 'sig_%s_%d.bin' % (tag, k))
    bits.tofile(path)
    r = subprocess.run([eng, str(q), str(k), path], capture_output=True, text=True)
    os.remove(path)
    val, cyc = None, None
    for line in r.stdout.splitlines():
        if line.startswith('RHOMAX '):
            a, p = map(int, line.split()[1:3])
            val = Fr(a, p)
        if line.startswith('CYCLE'):
            cyc = list(map(int, line.split()[1:]))
    if val is None:
        raise RuntimeError(r.stdout + r.stderr)
    return val, cyc


def check_uniform_cycle(rule, k, cyc, q=Q):
    """exact check of a closed walk of G_sigma at level k and of its rational periodic point; returns (x0, a, n)"""
    N, H = 1 << k, 1 << (k - 1)
    dmax = max(d for (_, d) in rule)
    signs = []
    for i, x in enumerate(cyc):
        y = cyc[(i + 1) % len(cyc)]
        if x % 2 == 0:
            t, s = x >> 1, 0
        else:
            s = rule_sign(rule, x, dmax)
            t = (q * x + s) >> 1
        assert (t - y) % H == 0
        signs.append(s)
    return cycle_rational([(x, k) for x in cyc], signs, q)


def rho_markov(rule, q=Q, want_cert=True):
    """exact rho_max of a rule by Markov refinement: (Fraction, number of leaves, periodic point x0, cycle length)"""
    M = Markov(rule, q)
    nodes, w, succ = M.graph()
    cyc = howard(w, succ)
    a = sum(w[i] for i in cyc)
    p = len(cyc)
    g = gcd(a, p)
    x0, aa, nn = cycle_rational([nodes[i] for i in cyc], [M.leaf[nodes[i]] for i in cyc], q)
    assert (aa, nn) == (a, p)
    if want_cert:
        psi = certify_upper(w, succ, a // g, p // g)
        if psi is None:
            raise RuntimeError('Howard missed a denser cycle')
    return Fr(a, p), len(nodes), x0, p


# ------------------------------------------------------------------------------------------------ the level-k game (task D)
def game_lfp(side, k, F, q=Q, cap=None, tag='g', mhtie=False):
    """least fixed point of Min's (side='min') or Max's (side='max') operator at level k and threshold F, by this
    lane's C engine; returns (pot int64 array, sig-or-lift uint8 array) or None if the run diverged past cap"""
    import numpy as np
    eng = build_engines()['game']
    N = 1 << k
    cap = cap or N * F.denominator + 1
    pre = os.path.join(SCR, '%s_%s_%d' % (tag, side, k))
    r = subprocess.run([eng, side, str(q), str(k), str(F.numerator), str(F.denominator), str(cap), pre] +
                       (['mhtie'] if mhtie else []), capture_output=True, text=True, check=True)
    if not r.stdout.startswith('FINITE'):
        return None
    pot = np.fromfile(pre + '.pot', dtype=np.int32)
    aux = np.fromfile(pre + ('.sig' if side == 'min' else '.lift'), dtype=np.uint8)
    os.remove(pre + '.pot')
    os.remove(pre + ('.sig' if side == 'min' else '.lift'))
    return pot, aux


def check_min_cert(k, F, pot, sig, q=Q, chunk=1 << 20):
    """exact Lemma G1 check on G_sigma at level k: pot(t) + e(s) <= pot(s) for both lifts t of the sigma-option
    (chunked, int64 arithmetic)"""
    import numpy as np
    N, H = 1 << k, 1 << (k - 1)
    for lo in range(0, N, chunk):
        s = np.arange(lo, min(N, lo + chunk), dtype=np.int64)
        odd = (s & 1) == 1
        sign = np.ones(len(s), dtype=np.int64)
        sign[odd] = np.where(sig[(s[odd] - 1) >> 1] == 1, -1, 1)
        P = np.where(odd, ((q * s + sign) >> 1) % H, (s >> 1) % H)
        e = np.where(odd, F.denominator - F.numerator, -F.numerator)
        ps = pot[s].astype(np.int64)
        ok = (pot[P].astype(np.int64) + e <= ps) & (pot[P + H].astype(np.int64) + e <= ps) & (ps >= 0)
        if not ok.all():
            return False
    return True


def check_max_cert(k, F, pot, lift, q=Q, chunk=1 << 20):
    """exact Lemma G2 check with W = Z/N: for every node and every option P, t = P + lift(P) H has pot(t) <= pot(s) + e(s)"""
    import numpy as np
    N, H = 1 << k, 1 << (k - 1)
    for lo in range(0, N, chunk):
        s = np.arange(lo, min(N, lo + chunk), dtype=np.int64)
        odd = (s & 1) == 1
        e = np.where(odd, F.denominator - F.numerator, -F.numerator)
        ps = pot[s].astype(np.int64)
        for sg in (1, -1):
            P = np.where(odd, ((q * s + sg) >> 1) % H, (s >> 1) % H)
            t = P + lift[P].astype(np.int64) * H
            if not (pot[t].astype(np.int64) <= ps + e).all():
                return False
    return True


# ------------------------------------------------------------------------------------------------ restricted game (search only)
_RLIB = None


def rlib():
    """the seven2 restricted-game engine (read-only source) compiled as a shared library into SCR/bin"""
    global _RLIB
    if _RLIB is None:
        import ctypes
        src = os.path.join(HERE, 'procgen_seven2_20260926_restrict.c')
        so = os.path.join(SCR, 'bin', 'restrict_%s.so' % sha(src)[:12])
        os.makedirs(os.path.dirname(so), exist_ok=True)
        if not os.path.exists(so):
            subprocess.run(['cc', '-O2', '-shared', '-fPIC', '-o', so, src], check=True)
        L = ctypes.CDLL(so)
        L.r_init.argtypes = [ctypes.c_int, ctypes.c_int64]
        L.r_mask_class.argtypes = [ctypes.c_int64, ctypes.c_int, ctypes.c_int]
        L.r_solve.argtypes = [ctypes.c_int64, ctypes.c_int64, ctypes.c_int64]
        L.r_solve.restype = ctypes.c_int64
        L.r_get_mask.argtypes = [ctypes.c_void_p]
        L.r_set_mask.argtypes = [ctypes.c_void_p]
        L.r_strategy.argtypes = [ctypes.c_void_p]
        L.r_get_values.argtypes = [ctypes.c_void_p]
        _RLIB = L
    return _RLIB


def greedy_tree(k, F, order='bfs', q=Q, cap=None, prefer=None):
    """greedy decision tree at level k and threshold F (Min's restricted energy game must stay finite): classes are
    processed from a queue (breadth first by depth); a class becomes a leaf with the preferred sign (MH first unless
    `prefer` says otherwise) if the restricted game stays finite, else it is split by the next bit.
    Returns the rule {(c, d): sign}.  Exploratory: the rule is re-evaluated exactly afterwards."""
    import numpy as np
    cap = cap or 20 * F.denominator + 200
    L = rlib()
    L.r_init(k, q)
    L.r_reset_values()
    H = 1 << (k - 1)
    if L.r_solve(F.numerator, F.denominator, cap) != 0:
        raise RuntimeError('infeasible at the start')
    mask = np.full(H, 3, dtype=np.uint8)
    rule = {}
    from collections import deque
    dq = deque([(1, 1)])
    while dq:
        c, d = dq.popleft()
        if d >= 2:
            s0 = mh(c, q) if prefer is None else prefer(c, d)
            tries = [s0, -s0]
        else:
            tries = []
        done = False
        for s in tries:
            m = 1 if s > 0 else 2
            L.r_save()
            L.r_get_mask(mask.ctypes.data)
            L.r_mask_class(c, d, m)
            if L.r_solve(F.numerator, F.denominator, cap) == 0:
                rule[(c, d)] = s
                done = True
                break
            L.r_restore()
            L.r_set_mask(mask.ctypes.data)
        if not done:
            if d >= k:
                raise RuntimeError('a residue admits no sign')
            dq.append((c, d + 1))
            dq.append((c + (1 << d), d + 1))
    return rule


# ------------------------------------------------------------------------------------------------ itinerary patterns
def parse_pattern(pat):
    """a negation-normalized itinerary pattern, e.g. '2 =2 !3 =>=5' (first sign +, then '=' same sign as the previous
    symbol, '!' opposite sign; valuation an int or '>=V' (last only), or '*' for 'any valuation' (last only, = '>=2')).
    Returns the two signed words (for first sign + and -)."""
    toks = pat.split()
    out = []
    for s0 in (1, -1):
        w, s = [], s0
        for i, t in enumerate(toks):
            if i > 0:
                rel, t = t[0], t[1:]
                s = s if rel == '=' else -s
            if t == '*':
                t = '>=2'
            if t.startswith('>='):
                w.append((s, ('ge', int(t[2:]))))
            else:
                w.append((s, int(t)))
        out.append(w)
    return out


def pattern_rule(flip_patterns, q=Q):
    """MH flipped on the union of the cylinders of the given normalized patterns (both signs); the cylinders must be
    pairwise disjoint"""
    fl = []
    for pat in flip_patterns:
        for w in parse_pattern(pat):
            fl.append(check_word_class(w, q))
    # disjointness
    fl = sorted(set(fl), key=lambda t: t[1])
    for i, (c, d) in enumerate(fl):
        for (c2, d2) in fl[:i]:
            if c % (1 << d2) == c2:
                raise ValueError('overlapping flip cylinders %s %s' % ((c, d), (c2, d2)))
    return flip_rule(fl, q)


def fmt_class(c, d, q=Q):
    sy, tail = itin_of_residue(c, d, q)
    s = ' '.join(('+' if a > 0 else '-') + str(v) for a, v in sy)
    if tail:
        s += ' ' + ('+' if tail[0] > 0 else '-') + '>=' + str(tail[1])
    return s.strip()


def fmt_rel(c, d, q=Q):
    """normalized (relative-sign) itinerary of a class"""
    sy, tail = itin_of_residue(c, d, q)
    toks, prev = [], None
    for a, v in sy + ([(tail[0], '>=' + str(tail[1]))] if tail else []):
        toks.append(str(v) if prev is None else ('=' if a == prev else '!') + str(v))
        prev = a
    return ' '.join(toks)


def rho_markov_c(rule, q=Q, check=True, want_crit=False):
    """exact rho_max of a rule by the C Markov engine (Howard + exact potential checked on every edge); the witness
    cycle is re-checked here as a rational periodic point (Lemma M).  Returns (Fraction, leaves, x0, cycle length)."""
    eng = build_engines()['markov']
    inp = '%d\n' % q + ''.join('%d %d %d\n' % (c, d, s) for (c, d), s in rule.items())
    r = subprocess.run([eng], input=inp, capture_output=True, text=True)
    val, nl, cyc, cert, crit = None, None, None, False, None
    for line in r.stdout.splitlines():
        t = line.split()
        if t[0] == 'LEAVES':
            nl = int(t[1])
        elif t[0] == 'RHO':
            val = Fr(int(t[1]), int(t[2]))
        elif t[0] == 'CERT':
            cert = (t[1] == 'ok')
        elif t[0] == 'CRIT':
            crit = (int(t[1]), float(t[2]))
        elif t[0] == 'CYCLE':
            a = list(map(int, t[1:]))
            cyc = [((a[i], a[i + 1]), a[i + 2]) for i in range(0, len(a), 3)]
    if val is None or not cert:
        raise RuntimeError('markov engine failed: ' + r.stdout[-500:] + r.stderr[-500:])
    x0, aa, nn = cycle_rational([c for c, _ in cyc], [s for _, s in cyc], q)
    if check:
        # the witness signs must be the rule's signs at the leaves (odd leaves lie inside rule classes)
        dmax = max(d for (_, d) in rule)
        for (c, d), s in cyc:
            if c % 2:
                assert rule_sign(rule, c, dmax) == s if d >= dmax else True
        assert Fr(aa, nn) == val
    if want_crit:
        return val, nl, x0, nn, crit
    return val, nl, x0, nn


def adversary_value(k, lift, q=Q, tag='a'):
    """exact value of one Max lift strategy at level k (C engine): the certified lower bound F (Lemma G2 on W, checked
    exactly inside the engine) and 'maximal' = from every node Min reaches a cycle of density <= F (exact lfp).
    Returns (Fraction, |W|, maximal ok?, witness cycle nodes)."""
    import numpy as np
    eng = build_engines()['adv']
    p = os.path.join(SCR, 'lift_%s_%d.bin' % (tag, k))
    np.asarray(lift, dtype=np.uint8).tofile(p)
    r = subprocess.run([eng, str(q), str(k), p], capture_output=True, text=True)
    os.remove(p)
    out = {}
    for line in r.stdout.splitlines():
        t = line.split()
        out[t[0]] = t[1:]
    if 'VALUE' not in out or out.get('CERT', ['x'])[0] != 'ok':
        raise RuntimeError(r.stdout + r.stderr)
    cyc = list(map(int, out['CYCLE']))
    return Fr(int(out['VALUE'][0]), int(out['VALUE'][1])), int(out['W'][0]), out['MAXIMAL'][0] == 'ok', cyc


def check_adv_cycle(k, lift, cyc, q=Q):
    """the witness is a closed walk of Min's graph against the lift strategy (some sign at each odd node)"""
    N, H = 1 << k, 1 << (k - 1)
    for i, s in enumerate(cyc):
        t = cyc[(i + 1) % len(cyc)]
        if s % 2 == 0:
            P = (s >> 1) % H
            assert t == P + int(lift[P]) * H
        else:
            ok = False
            for sg in (1, -1):
                P = ((q * s + sg) >> 1) % H
                ok |= (t == P + int(lift[P]) * H)
            assert ok
    return sum(s & 1 for s in cyc), len(cyc)


def compress_strategy(minus, k, q=Q):
    """the minimal decision tree (bit trie) of a level-k strategy (uint8 over odd nodes, 1 = '-'): a class is a leaf iff
    all its odd residues carry the same sign.  Returns the rule {(c, d): sign}."""
    import numpy as np
    cur = np.asarray(minus, dtype=np.int8).copy()        # depth k, index i <-> residue 2i+1
    levels = {k: cur}
    d = k
    while d > 1:
        half = 1 << (d - 2)
        a, b = cur[:half], cur[half:]
        par = np.where(a == b, a, -1).astype(np.int8)
        d -= 1
        levels[d] = par
        cur = par
    rule = {}
    # leaves: pure at depth d with a mixed parent (or d = 1)
    for d in range(1, k + 1):
        arr = levels[d]
        if d == 1:
            if arr[0] >= 0:
                rule[(1, 1)] = -1 if arr[0] == 1 else 1
            continue
        par = levels[d - 1]
        half = 1 << (d - 2)
        idx = np.flatnonzero(arr >= 0)
        pidx = idx % half
        sel = idx[par[pidx] < 0]
        for i in sel.tolist():
            rule[(2 * i + 1, d)] = -1 if arr[i] == 1 else 1
    return rule
