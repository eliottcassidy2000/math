#!/usr/bin/env python3
"""
procgen_rank_20260926_lib.py -- shared exact tools for the two-place Lyapunov (rank) lane
(session collatz-procgen-20260922, 2026-09-26).

Conventions (as in THM-4474/THM-4482):
  * T is the shortcut Collatz map T(x) = x/2 (x even), (3x+1)/2 (x odd), on Z and on Z_(2)
    (rationals with odd denominator; parity = parity of the numerator).
  * T_sigma is a level-k sign strategy: (3x + sigma(x mod 2^k))/2 on odd x.
  * v2 of a nonzero rational = v2(numerator) - v2(denominator); H(u/D) = max(|u|, |D|).
  * A "bank" is a list of (center, weight) pairs; Phi(n) = sum weight * v2(n - center).
Every printed claim of the lane goes through check(...), which raises on failure.
"""
import math
import resource
import sys
from fractions import Fraction

LN2 = math.log(2.0)
LN3 = math.log(3.0)
KAPPA = math.log(1.5)             # chi(-1) = log(3/2), the Lyapunov exponent of the fixed point -1
C_CRIT = LN2 / LN3                # log_3 2
H_BITS = -(C_CRIT * math.log2(C_CRIT) + (1 - C_CRIT) * math.log2(1 - C_CRIT))   # h(log_3 2) = 0.94996


def check(cond, msg):
    """Every claim printed by the lane is a check: raise on failure, print 'ok: msg' on success."""
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)
    print("ok: " + msg)
    sys.stdout.flush()


def say(msg=""):
    print(msg)
    sys.stdout.flush()


def peak_rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024.0 * 1024.0) if sys.platform == "darwin" else r / 1024.0


# ----------------------------------------------------------------------------------------------
# valuations, heights, the map
# ----------------------------------------------------------------------------------------------
def v2int(n):
    if n == 0:
        raise ValueError("v2(0)")
    n = abs(n)
    return (n & -n).bit_length() - 1


def v2(q):
    q = Fraction(q)
    if q == 0:
        raise ValueError("v2(0)")
    return v2int(q.numerator) - v2int(q.denominator)


def height(q):
    q = Fraction(q)
    return max(abs(q.numerator), abs(q.denominator))


def is_odd(x):
    """parity of an element of Z_(2) (odd denominator assumed)"""
    x = Fraction(x)
    assert x.denominator % 2 == 1
    return x.numerator % 2 != 0


def T(x):
    """shortcut Collatz map on ints (fast path) or on Fractions with odd denominator"""
    if isinstance(x, int):
        return (3 * x + 1) >> 1 if x & 1 else x >> 1
    x = Fraction(x)
    return (3 * x + 1) / 2 if is_odd(x) else x / 2


def T_sigma(x, sig, k):
    """level-k sign strategy on ints; sig maps odd residues mod 2^k to +-1"""
    if x & 1:
        return (3 * x + sig[x % (1 << k)]) >> 1
    return x >> 1


def orbit_to_one(n, cap=10 ** 7):
    """list n, T n, ..., 1 (shortcut map); raises if cap exceeded"""
    out = [n]
    while n != 1:
        n = T(n)
        out.append(n)
        if len(out) > cap:
            raise RuntimeError("orbit too long")
    return out


def parities_of_orbit(orb, tail=0):
    """parity word along an orbit list; optionally append `tail` letters of the (1 0)^inf tail after 1"""
    w = [x & 1 for x in orb]
    if tail and orb[-1] == 1:
        # after 1: 1 -> 2 -> 1 -> 2 ...; the parity of 1 is already in w, continue 0,1,0,1,...
        w.extend([0, 1] * (tail // 2 + 1))
    return w


def parity_word(x, L):
    """first L parities of the orbit of x (int or Fraction)"""
    out = []
    for _ in range(L):
        if isinstance(x, int):
            out.append(x & 1)
        else:
            out.append(1 if is_odd(x) else 0)
        x = T(x)
    return out


# ----------------------------------------------------------------------------------------------
# periodic and preperiodic points
# ----------------------------------------------------------------------------------------------
def affine_of_word(word):
    """T^p along a parity word as x -> (A x + C)/B, A = 3^a, B = 2^p"""
    A, B, C = 1, 1, 0
    for b in word:
        if b:
            A, C = 3 * A, 3 * C + B
        B *= 2
    return A, B, C


def periodic_point(word):
    """the unique fixed point in Z_2 of T^p along `word` (Lagarias; Banach): C/(B - A)"""
    A, B, C = affine_of_word(word)
    return Fraction(C, B - A)


def is_expanding_word(word):
    a = sum(word)
    return 3 ** a > 2 ** len(word)


def chi_of_word(word):
    """Lyapunov exponent log(3^a/2^p)/p of the cycle of `word`"""
    a, p = sum(word), len(word)
    return (a * LN3 - p * LN2) / p


def cycle_points(word):
    """the p points of the orbit of the periodic point of `word` (rotations)"""
    x = periodic_point(word)
    pts = [x]
    for _ in range(len(word) - 1):
        pts.append(T(pts[-1]))
    assert T(pts[-1]) == x
    return pts


def preimages(x):
    """the T-preimages of x in Z_(2): 2x always, (2x-1)/3 if that is odd"""
    x = Fraction(x)
    out = [2 * x]
    y = (2 * x - 1) / 3
    if y.denominator % 2 == 1 and is_odd(y):
        out.append(y)
    return out


def common_prefix(a, b):
    m = min(len(a), len(b))
    for i in range(m):
        if a[i] != b[i]:
            return i
    return m


# ----------------------------------------------------------------------------------------------
# banks
# ----------------------------------------------------------------------------------------------
def Phi(n, bank):
    """sum c * v2(n - beta) (n int or Fraction not equal to any center)"""
    s = 0.0
    for beta, c in bank:
        s += c * v2(Fraction(n) - beta)
    return s


def Phi_star(x, bank):
    """the potential at x with the atom at x removed"""
    s = 0.0
    for beta, c in bank:
        if beta != x:
            s += c * v2(Fraction(x) - beta)
    return s


def atom(x, bank):
    return sum(c for beta, c in bank if beta == x)


def good_integer(z, D, lo_exp=1):
    """the unique integer y in [2^(D+lo_exp), 2^(D+lo_exp+1)) with v2(y - z) = D exactly (z in Z_(2))"""
    z = Fraction(z)
    M = 1 << (D + 1)
    r = (z.numerator * pow(z.denominator, -1, M)) % M          # z mod 2^(D+1)
    r = (r + (1 << D)) % M                                      # flip bit D: exact depth D
    base = 1 << (D + lo_exp)
    y = base + ((r - base) % M)
    assert v2(Fraction(y) - z) == D and base <= y < base + M
    return y


# ----------------------------------------------------------------------------------------------
# necklaces and the forced mass
# ----------------------------------------------------------------------------------------------
def mobius(n):
    res, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    if m > 1:
        res = -res
    return res


def n_primitive_necklaces(p, a):
    g = math.gcd(p, a) if a > 0 else p
    tot = 0
    for d in range(1, g + 1):
        if g % d == 0 and a % d == 0:
            tot += mobius(d) * math.comb(p // d, a // d)
    assert tot % p == 0
    return tot // p


# ----------------------------------------------------------------------------------------------
# strategy-cube graphs
# ----------------------------------------------------------------------------------------------
def cube_graph(k, sig):
    """G_sigma at level k: node s -> the two lifts of T_sigma(s) mod 2^(k-1)"""
    M = 1 << k
    H = M >> 1
    adj = []
    for s in range(M):
        t = s // 2 if s % 2 == 0 else (3 * s + sig[s]) // 2
        t %= H
        adj.append((t, t + H))
    return adj


def sccs(adj):
    """iterative Tarjan; returns the list of strongly connected components"""
    n = len(adj)
    index = [-1] * n
    low = [0] * n
    onst = [False] * n
    st = []
    comps = []
    counter = 0
    for root in range(n):
        if index[root] != -1:
            continue
        work = [(root, 0)]
        while work:
            v, i = work[-1]
            if i == 0:
                index[v] = low[v] = counter
                counter += 1
                st.append(v)
                onst[v] = True
            if i < len(adj[v]):
                work[-1] = (v, i + 1)
                w = adj[v][i]
                if index[w] == -1:
                    work.append((w, 0))
                elif onst[w]:
                    low[v] = min(low[v], index[w])
            else:
                work.pop()
                if work:
                    u = work[-1][0]
                    low[u] = min(low[u], low[v])
                if low[v] == index[v]:
                    comp = []
                    while True:
                        w = st.pop()
                        onst[w] = False
                        comp.append(w)
                        if w == v:
                            break
                    comps.append(comp)
    return comps


def karp_max_mean(nodes, adj):
    """maximum cycle mean of w = log(3/2) (odd), -log 2 (even) on the subgraph induced by `nodes`"""
    S = set(nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    w = [KAPPA if v % 2 else -LN2 for v in nodes]
    NEG = -1e300
    D = [[NEG] * n for _ in range(n + 1)]
    D[0][0] = 0.0
    for m in range(1, n + 1):
        prev, cur = D[m - 1], D[m]
        for i, v in enumerate(nodes):
            if prev[i] == NEG:
                continue
            val = prev[i] + w[i]
            for u in adj[v]:
                if u in S:
                    j = idx[u]
                    if val > cur[j]:
                        cur[j] = val
    best = NEG
    for j in range(n):
        if D[n][j] == NEG:
            continue
        worst = min((D[n][j] - D[m][j]) / (n - m) for m in range(n) if D[m][j] != NEG)
        best = max(best, worst)
    return best


def classify_strategy(k, sig):
    """'i' (no expanding cycle), "i'" (every SCC containing an expanding cycle is that single cycle),
    or 'other'; also returns the isolated expanding SCCs"""
    adj = cube_graph(k, sig)
    iso, bad, any_exp = [], False, False
    for comp in sccs(adj):
        S = set(comp)
        internal = [sum(1 for u in adj[v] if u in S) for v in comp]
        if len(comp) == 1 and internal[0] == 0:
            continue
        if all(d == 1 for d in internal):
            odd = sum(1 for v in comp if v % 2)
            if 3 ** odd > 2 ** len(comp):
                iso.append(sorted(comp))
                any_exp = True
        else:
            if karp_max_mean(comp, adj) > -1e-12:
                bad = True
                any_exp = True
    if not any_exp:
        return "i", iso
    return ("i'" if not bad else "other"), iso
