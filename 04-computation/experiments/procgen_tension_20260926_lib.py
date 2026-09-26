#!/usr/bin/env python3
"""
procgen_tension_20260926_lib.py -- shared exact helpers for the tension lane
(session collatz-procgen-20260922, 2026-09-26: ranks = certificates, Christoffel maximizers,
Christoffel rigidity, the negation duality, typed dictionary).

Setting (THM-4474).  A level-k sign strategy sigma: odd residues mod 2^k -> {+1,-1};
    T_sigma(n) = n/2 (n even),   (3n + sigma(n mod 2^k))/2 (n odd).
Mask encoding (as in the cube lanes): bit i <-> residue 2i+1, bit value 1 = minus.
Mask 0 is Collatz (3n+1); the all-ones mask is 3n-1.
Parity graph G_sigma on Z/2^k: s -> the two lifts mod 2^k of T(s) mod 2^(k-1).
Node weight w(s) = log(3/2) (s odd), -log 2 (s even).  A cycle with a odd nodes and length p
is expanding iff 3^a > 2^p (equality never occurs); class (i) <=> every cycle contracting.

Everything used as a proof step or printed as a claim is exact (integers / Fractions);
floating point is used only for display and for explicitly labelled numerical sanity checks.
"""
from fractions import Fraction
from math import comb, gcd, log
import numpy as np

LOG2 = log(2.0)
LOG3 = log(3.0)
CRIT = LOG2 / LOG3          # c = log_3 2 (display only)


def check(cond, msg):
    """every printed claim goes through here; a failure aborts the run"""
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


# ----------------------------------------------------------------------------- rationals vs c
def below_c(f):
    """exact test f < log_3 2  (f = a/p >= 0):  3^a < 2^p"""
    f = Fraction(f)
    return 3 ** f.numerator < 2 ** f.denominator


def best_lower(D):
    """F_D: the largest a/d < log_3 2 with 1 <= d <= D (exact)"""
    best = Fraction(0)
    for d in range(1, D + 1):
        a = 0
        while 3 ** (a + 1) < 2 ** d:
            a += 1
        best = max(best, Fraction(a, d))
    return best


# ----------------------------------------------------------------------------- strategies
def sign(mask, r):
    return -1 if (mask >> ((r - 1) // 2)) & 1 else 1


def mask_from_minus(minus_residues):
    m = 0
    for r in minus_residues:
        m |= 1 << ((r - 1) // 2)
    return m


def minus_residues(mask, k):
    return [2 * i + 1 for i in range(1 << (k - 1)) if (mask >> i) & 1]


def all_minus(k):
    return (1 << (1 << (k - 1))) - 1


def chi4_mask(k):
    """chi_{-4}: + on 1 mod 4, - on 3 mod 4 (the all-d strategy)"""
    return mask_from_minus(r for r in range(3, 1 << k, 4))


def nu(mask, k):
    """(nu sigma)(m) = -sigma(-m mod 2^k)"""
    M = 1 << k
    out = 0
    for r in range(1, M, 2):
        if -sign(mask, M - r) == -1:
            out |= 1 << ((r - 1) // 2)
    return out


def uset(mask, k):
    """u-residues: sigma(r) != chi_{-4}(r) (exactly one halving follows the odd step)"""
    return [r for r in range(1, 1 << k, 2) if sign(mask, r) != (1 if r % 4 == 1 else -1)]


def lift(mask, k, K):
    """the same map viewed as a level-K strategy (K >= k)"""
    out = 0
    for r in range(1, 1 << K, 2):
        if sign(mask, r % (1 << k)) == -1:
            out |= 1 << ((r - 1) // 2)
    return out


def step(n, k, mask):
    """T_sigma on an integer (any sign) or on a residue"""
    if n % 2 == 0:
        return n // 2
    return (3 * n + sign(mask, n % (1 << k))) // 2


# ----------------------------------------------------------------------------- the parity graph
def succ(k, mask, s):
    M = 1 << k
    H = M >> 1
    if s % 2 == 0:
        b = (s // 2) % H
    else:
        b = ((3 * s + sign(mask, s)) // 2) % H
    return (b, b + H)


def edges(k, mask):
    return [(s, t) for s in range(1 << k) for t in succ(k, mask, s)]


def adjacency(k, mask):
    return {s: succ(k, mask, s) for s in range(1 << k)}


def karp(k, mask):
    """exact maximum odd density over the cycles of G_sigma (Karp, integer DP)"""
    M = 1 << k
    n = M
    NEG = None
    pred = {t: [] for t in range(M)}
    for s in range(M):
        for t in succ(k, mask, s):
            pred[t].append(s)
    D = [[0] * M]
    for _ in range(n):
        prev = D[-1]
        cur = []
        for t in range(M):
            vals = [prev[s] + (s & 1) for s in pred[t] if prev[s] is not NEG]
            cur.append(max(vals) if vals else NEG)
        D.append(cur)
    best = None
    for v in range(M):
        if D[n][v] is NEG:
            continue
        cands = [Fraction(D[n][v] - D[i][v], n - i) for i in range(n) if D[i][v] is not NEG]
        mn = min(cands)
        if best is None or mn > best:
            best = mn
    return best


def _preds_vec(k):
    M = 1 << k
    H = M >> 1
    sp = np.zeros(H, dtype=np.int64)
    sm = np.zeros(H, dtype=np.int64)
    for s in range(1, M, 2):
        sp[((3 * s + 1) // 2) % H] = s
        sm[((3 * s - 1) // 2) % H] = s
    return sp, sm


def karp_batch(k, masks):
    """vectorized exact Karp over many strategies of one level.
    D_j depends on a node only through its class mod 2^(k-1) (both lifts have the same in-neighbours),
    so the DP runs on Z/2^(k-1).  Returns exact (numerator, denominator) arrays of rho_max."""
    M = 1 << k
    H = M >> 1
    n = M
    NEG = -10 ** 6
    sp, sm = _preds_vec(k)
    S = len(masks)
    if k <= 6:
        masks = np.asarray(masks, dtype=np.int64)
        okp = ((masks[:, None] >> ((sp - 1) // 2)[None, :]) & 1) == 0   # s_+(b) has sign +
        okm = ((masks[:, None] >> ((sm - 1) // 2)[None, :]) & 1) == 1   # s_-(b) has sign -
    else:   # masks wider than 63 bits: build the sign tables in Python
        okp = np.array([[((int(m) >> int((s - 1) // 2)) & 1) == 0 for s in sp] for m in masks], dtype=bool)
        okm = np.array([[((int(m) >> int((s - 1) // 2)) & 1) == 1 for s in sm] for m in masks], dtype=bool)
    ev = np.array([(2 * b) % H for b in range(H)])
    E = np.zeros((n + 1, S, H), dtype=np.int32)
    for j in range(1, n + 1):
        prev = E[j - 1]
        c0 = prev[:, ev]
        c1 = np.where(okp, prev[:, sp % H] + 1, NEG)
        c2 = np.where(okm, prev[:, sm % H] + 1, NEG)
        cur = np.maximum(np.maximum(c0, c1), c2)
        cur[cur < NEG // 2] = NEG
        E[j] = cur
    En = E[n].astype(np.int64)
    bn = np.full((S, H), 10 ** 9, dtype=np.int64)
    bd = np.ones((S, H), dtype=np.int64)
    for i in range(n):
        Ei = E[i].astype(np.int64)
        num = np.where(Ei > NEG // 2, En - Ei, 10 ** 9)
        den = n - i
        less = num * bd < bn * den
        bn = np.where(less, num, bn)
        bd = np.where(less, den, bd)
    valid = En > NEG // 2
    rn = np.full(S, -1, dtype=np.int64)
    rd = np.ones(S, dtype=np.int64)
    for b in range(H):
        gt = valid[:, b] & (bn[:, b] * rd > rn * bd[:, b])
        rn = np.where(gt, bn[:, b], rn)
        rd = np.where(gt, bd[:, b], rd)
    return rn, rd


def rho_all(k, batch=8192):
    """exact rho_max of every level-k strategy (k <= 5), as a list of Fractions indexed by mask"""
    NM = 1 << (1 << (k - 1))
    out = []
    for st in range(0, NM, batch):
        rn, rd = karp_batch(k, np.arange(st, min(NM, st + batch)))
        out.extend(Fraction(int(a), int(b)) for a, b in zip(rn, rd))
    return out


def potential(k, mask, F):
    """Lemma P: least integer psi >= 0 with psi(t) <= psi(s) - w_F(s) on every edge, where
    F = q/r, w_F = r-q (odd), -q (even); psi(s) = sup over walks from s of partial sums of w_F.
    Returns None iff some cycle has density > F."""
    F = Fraction(F)
    q, r = F.numerator, F.denominator
    M = 1 << k
    adj = adjacency(k, mask)
    w = [(r - q) if s & 1 else -q for s in range(M)]
    psi = [0] * M
    for _ in range(M + 2):
        new = [max(0, w[s] + max(psi[t] for t in adj[s])) for s in range(M)]
        if new == psi:
            return psi
        psi = new
    return None


def potential_ok(k, mask, F, psi):
    F = Fraction(F)
    q, r = F.numerator, F.denominator
    if min(psi) < 0:
        return False
    for s in range(1 << k):
        ws = (r - q) if s & 1 else -q
        for t in succ(k, mask, s):
            if not psi[t] <= psi[s] - ws:
                return False
    return True


def tight_edges(k, mask, F, psi):
    F = Fraction(F)
    q, r = F.numerator, F.denominator
    out = []
    for s in range(1 << k):
        ws = (r - q) if s & 1 else -q
        for t in succ(k, mask, s):
            if psi[t] == psi[s] - ws:
                out.append((s, t))
    return out


def simple_cycles_adj(adj):
    """all simple cycles of a small digraph given as {v: iterable of successors}; each cycle once,
    listed from its least node (includes self-loops)"""
    out = []
    for st in sorted(adj):
        stack = [(st, [st], {st})]
        while stack:
            v, path, seen = stack.pop()
            for t in adj[v]:
                if t == st:
                    out.append(list(path))
                elif t > st and t not in seen:
                    stack.append((t, path + [t], seen | {t}))
    return out


def density(cyc):
    return Fraction(sum(1 for s in cyc if s & 1), len(cyc))


def cycle_word(cyc):
    return ''.join('1' if s & 1 else '0' for s in cyc)


def canon(word):
    """lexicographically largest rotation"""
    return max(word[i:] + word[:i] for i in range(len(word)))


# ----------------------------------------------------------------------------- Collatz words
def collatz_word(r, L):
    """Collatz (all +) parity word of length L of the residue / integer r (first L letters
    depend only on r mod 2^L)"""
    w = []
    n = r
    for _ in range(L):
        w.append(n & 1)
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
    return w


def ballot(word):
    """all prefixes j have 3^(a_j) > 2^j"""
    a = 0
    for j, b in enumerate(word, 1):
        a += b
        if not 3 ** a > 2 ** j:
            return False
    return True


def bad_set(k):
    return [r for r in range(1, 1 << k, 2) if ballot(collatz_word(r, k))]


def sigma_k_mask(k):
    """THM-4479 / cube-distance Theorem 1: flip exactly Bad_k"""
    return mask_from_minus(bad_set(k))


def first_descent(word):
    """least j with 3^(a_j) < 2^j, or None"""
    a = 0
    for j, b in enumerate(word, 1):
        a += b
        if 3 ** a < 2 ** j:
            return j
    return None


def is_fd_word(word):
    """first-descent word: descends exactly at its last letter"""
    return first_descent(word) == len(word)


def upper_christoffel(a, d):
    """upper Christoffel word of slope a/d (prefix counts ceil(a j / d)), as a 0/1 list"""
    w = []
    for j in range(1, d + 1):
        w.append(-(-a * j // d) - (-(-a * (j - 1) // d)))
    return w


def lower_christoffel(a, d):
    w = []
    for j in range(1, d + 1):
        w.append(a * j // d - a * (j - 1) // d)
    return w


def balanced(word):
    """cyclically balanced: all cyclic factors of equal length differ by <= 1 in their number of ones"""
    L = len(word)
    ww = word + word
    for m in range(1, L):
        cnt = [sum(ww[i:i + m]) for i in range(L)]
        if max(cnt) - min(cnt) > 1:
            return False
    return True


def strictly_above_words(a, d, m=1):
    """binary words of length m d with m a ones whose prefix counts exceed (a/d) j for 0 < j < m d
    (exact enumeration by DP over prefixes)"""
    L = m * d
    out = []

    def rec(j, cnt, w):
        if j == L:
            if cnt == m * a:
                out.append(list(w))
            return
        for b in (1, 0):
            c2 = cnt + b
            j2 = j + 1
            if c2 > m * a or (m * a - c2) > (L - j2):
                continue
            if j2 < L and not c2 * d > a * j2:
                continue
            w.append(b)
            rec(j2, c2, w)
            w.pop()
    rec(0, 0, [])
    return out


def necklace_count(n, j):
    """binary necklaces of length n with exactly j ones: (1/n) sum_{e | gcd(n,j)} phi(e) C(n/e, j/e)"""
    g = gcd(n, j)
    tot = 0
    for e in range(1, g + 1):
        if g % e == 0:
            tot += _phi(e) * comb(n // e, j // e)
    return tot // n


def _phi(x):
    res, p, y = x, 2, x
    while p * p <= y:
        if y % p == 0:
            while y % p == 0:
                y //= p
            res -= res // p
        p += 1
    if y > 1:
        res -= res // y
    return res


def v2(x):
    """2-adic valuation of a nonzero rational"""
    x = Fraction(x)
    if x == 0:
        raise ValueError
    n, d = x.numerator, x.denominator
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    while d % 2 == 0:
        d //= 2
        v -= 1
    return v


def periodic_point(word_signs):
    """the 2-adic periodic point (a rational) of the signed parity word:
    letters 0 (even: x -> x/2) or +1/-1 (odd: x -> (3x + s)/2); fixed point of the composite"""
    A, B = Fraction(1), Fraction(0)     # composite map x -> A x + B
    for s in word_signs:
        if s == 0:
            A, B = A / 2, B / 2
        else:
            A, B = 3 * A / 2, (3 * B + s) / 2
    return B / (1 - A)


def claim(cond, msg):
    """check and print one claim of the output"""
    check(cond, msg)
    print("  ok: " + msg, flush=True)


def cycle_in_edges(edge_list):
    """some directed cycle (node list) inside an edge list, or None"""
    adj = {}
    for s, t in edge_list:
        adj.setdefault(s, []).append(t)
    color = {}
    for root in list(adj):
        if root in color:
            continue
        stack = [(root, iter(adj.get(root, [])))]
        path = [root]
        color[root] = 1
        while stack:
            v, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                color[v] = 2
                stack.pop()
                path.pop()
                continue
            if color.get(nxt) == 1:
                return path[path.index(nxt):]
            if nxt not in color:
                color[nxt] = 1
                stack.append((nxt, iter(adj.get(nxt, []))))
                path.append(nxt)
    return None


def signed_word(k, mask, cyc):
    return [0 if s % 2 == 0 else sign(mask, s) for s in cyc]


def residue_of(x, m):
    """a rational with odd denominator, reduced mod 2^m"""
    x = Fraction(x)
    M = 1 << m
    return (x.numerator * pow(x.denominator, -1, M)) % M


def karp_single(k, mask):
    """exact rho_max of one strategy with O(2^k) memory (two DP passes, numpy); for k up to ~14"""
    M = 1 << k
    H = M >> 1
    n = M
    NEG = -10 ** 9
    sp, sm = _preds_vec(k)
    okp = np.array([sign(mask, int(s)) == 1 for s in sp])
    okm = np.array([sign(mask, int(s)) == -1 for s in sm])
    ev = np.array([(2 * b) % H for b in range(H)])

    def stepE(prev):
        c0 = prev[ev]
        c1 = np.where(okp, prev[sp % H] + 1, NEG)
        c2 = np.where(okm, prev[sm % H] + 1, NEG)
        cur = np.maximum(np.maximum(c0, c1), c2)
        cur[cur < NEG // 2] = NEG
        return cur
    E = np.zeros(H, dtype=np.int64)
    for _ in range(n):
        E = stepE(E)
    En = E.copy()
    bn = np.full(H, 10 ** 12, dtype=np.int64)
    bd = np.ones(H, dtype=np.int64)
    E = np.zeros(H, dtype=np.int64)
    for i in range(n):
        num = np.where(E > NEG // 2, En - E, 10 ** 12)
        less = num * bd < bn * (n - i)
        bn = np.where(less, num, bn)
        bd = np.where(less, n - i, bd)
        E = stepE(E)
    best = None
    for b in range(H):
        if En[b] > NEG // 2:
            f = Fraction(int(bn[b]), int(bd[b]))
            if best is None or f > best:
                best = f
    return best


def block_decomposition(k, mask_sigma_k, cyc):
    """blocks of the sigma_k-orbit around a cycle of G_(sigma_k): at a node in Bad_k an F-block of length 2,
    elsewhere a C-block of length d = first descent of the Collatz word of the node.  Returns the list of
    (start index, length, kind) over one period of the (eventually periodic) block sequence, or None."""
    p = len(cyc)
    bad = set(bad_set(k))
    i = 0
    seen = {}
    blocks = []
    while (i % p) not in seen:
        seen[i % p] = len(blocks)
        s = cyc[i % p]
        if s in bad:
            L, kind = 2, 'F'
        else:
            L, kind = first_descent(collatz_word(s, k)), 'C'
            if L is None:
                return None
        blocks.append((i % p, L, kind))
        i += L
    j = seen[i % p]
    per = blocks[j:]
    if sum(L for _, L, _ in per) % p != 0:
        return None
    return per


def as_cycle(k, mask, cyc):
    """verify that an ordered node list is a simple cycle of G_sigma and return it"""
    check(len(set(cyc)) == len(cyc) and len(cyc) >= 1, "as_cycle: repeated node")
    for i in range(len(cyc)):
        check(cyc[(i + 1) % len(cyc)] in succ(k, mask, cyc[i]), "as_cycle: missing edge")
    return list(cyc)


def potential_np(k, mask, F, max_iter=None):
    """numpy version of potential() for large k (same least integer potential, or None)"""
    F = Fraction(F)
    q, r = F.numerator, F.denominator
    M = 1 << k
    H = M >> 1
    s = np.arange(M)
    sg = np.array([0 if x % 2 == 0 else sign(mask, x) for x in range(M)], dtype=np.int64)
    b = np.where(s % 2 == 0, (s // 2) % H, ((3 * s + sg) // 2) % H)
    t0, t1 = b, b + H
    w = np.where(s % 2 == 1, r - q, -q).astype(np.int64)
    psi = np.zeros(M, dtype=np.int64)
    for _ in range((M + 2) if max_iter is None else max_iter):
        new = np.maximum(0, w + np.maximum(psi[t0], psi[t1]))
        if np.array_equal(new, psi):
            ok = bool(np.all(psi[t0] <= psi - w) and np.all(psi[t1] <= psi - w) and psi.min() >= 0)
            return [int(v) for v in psi] if ok else None
        psi = new
    return None
