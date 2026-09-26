#!/usr/bin/env python3
"""
procgen_cube_20260925_core.py -- exact pure-Python reference implementation for the
strategy cube of Althofer's 3n+-1 game (session collatz-procgen-20260922, cube lane).

A sign strategy of level k is sigma : {odd residues mod 2^k} -> {+1,-1}; the map is
    T_sigma(n) = n/2 (n even),  (3n + sigma(n mod 2^k))/2 (n odd).
Encoding: mask bit i  <->  residue r = 2i+1;  bit = 1 means sigma(r) = -1.
So mask 0 is Collatz (all +) and the all-ones mask is 3n-1.

Everything here is exact (integers / Fractions) and deliberately simple; the C engine
(procgen_cube_20260925_engine.c) is the fast path and is cross-checked against this file.

Objects (all PROVED facts used are stated in the note, section 1):
  * the parity-transition graph G_sigma on Z/2^k: s -> the two lifts of T(s) mod 2^(k-1);
    paths of length L <-> classes mod 2^(k+L) (bijection), step weight
    w(s) = log(3/2) (s odd), -log 2 (s even);
  * a cycle is EXPANDING iff 3^a > 2^p (a odd steps, length p), CONTRACTING iff 3^a < 2^p;
  * Bad_L = classes mod 2^(k+L) whose first L prefixes all have 3^(a_j) > 2^j.
"""
from fractions import Fraction
import math
import sys

LOG2 = math.log(2.0)
LOG3 = math.log(3.0)


# ----------------------------------------------------------------------------- basics
def sig_table(k, mask):
    """list indexed by residue mod 2^k; entries at odd residues are +-1, 0 at even ones."""
    K = 1 << k
    t = [0] * K
    for i in range(K >> 1):
        t[2 * i + 1] = -1 if (mask >> i) & 1 else 1
    return t


def T(n, k, sig):
    if n % 2 == 0:
        return n // 2
    return (3 * n + sig[n % (1 << k)]) // 2


def negate_mask(k, mask):
    """nu: sigma'(m) = -sigma(-m mod 2^k).  In bits: complement of the bit-reversed mask."""
    H = 1 << (k - 1)
    out = 0
    for i in range(H):
        j = H - 1 - i            # index of -m when m = 2i+1
        bit = (mask >> j) & 1
        out |= (1 - bit) << i
    return out


def lift_mask(k, mask):
    """the level-(k+1) strategy with the same map."""
    H = 1 << (k - 1)
    return mask | (mask << H)


def is_primitive(k, mask):
    """not a lift of a level-(k-1) strategy."""
    if k == 1:
        return True
    Hh = 1 << (k - 2)
    lo = mask & ((1 << Hh) - 1)
    hi = mask >> Hh
    return lo != hi


def ud_word(k, mask):
    """u/d reading: at odd r, 'd' iff 4 | 3r + sigma(r) (>= 2 halvings), 'u' iff one halving."""
    sig = sig_table(k, mask)
    w = []
    for r in range(1, 1 << k, 2):
        w.append('d' if (3 * r + sig[r]) % 4 == 0 else 'u')
    return ''.join(w)


def sigma_string(k, mask):
    sig = sig_table(k, mask)
    return ''.join('+' if sig[r] > 0 else '-' for r in range(1, 1 << k, 2))


# ----------------------------------------------------------------------------- graph
def graph(k, mask):
    K = 1 << k
    H = K >> 1
    sig = sig_table(k, mask)
    succ = []
    for s in range(K):
        t = s // 2 if s % 2 == 0 else (3 * s + sig[s]) // 2
        t0 = t % H if H else 0
        succ.append((t0, t0 + H))
    return succ


def tarjan_scc(succ):
    """iterative Tarjan; returns list of SCCs (lists of nodes)."""
    n = len(succ)
    index = [None] * n
    low = [0] * n
    onstack = [False] * n
    stack = []
    comps = []
    counter = 0
    for root in range(n):
        if index[root] is not None:
            continue
        work = [(root, 0)]
        index[root] = low[root] = counter
        counter += 1
        stack.append(root)
        onstack[root] = True
        while work:
            v, ei = work[-1]
            if ei < len(succ[v]):
                work[-1] = (v, ei + 1)
                w = succ[v][ei]
                if index[w] is None:
                    index[w] = low[w] = counter
                    counter += 1
                    stack.append(w)
                    onstack[w] = True
                    work.append((w, 0))
                elif onstack[w]:
                    low[v] = min(low[v], index[w])
            else:
                work.pop()
                if work:
                    u = work[-1][0]
                    low[u] = min(low[u], low[v])
                if low[v] == index[v]:
                    comp = []
                    while True:
                        w = stack.pop()
                        onstack[w] = False
                        comp.append(w)
                        if w == v:
                            break
                    comps.append(sorted(comp))
    return comps


def bottom_sccs(succ, comps):
    out = []
    for c in comps:
        cs = set(c)
        if all(t in cs for v in c for t in succ[v]):
            # a singleton without self-loop cannot be bottom (out-degree 2), fine
            out.append(c)
    return out


def simple_cycles(succ, nodes=None):
    """Johnson-style enumeration of simple cycles (as node lists) inside `nodes`.
    Used only for small graphs as an independent exact check."""
    if nodes is None:
        nodes = list(range(len(succ)))
    nodes = sorted(nodes)
    allowed = set(nodes)
    cycles = []
    # simple DFS from each start s, only through nodes > s (canonical: s is the minimum)
    for s in nodes:
        stack = [(s, iter([t for t in succ[s] if t in allowed and t >= s]))]
        path = [s]
        onpath = {s}
        while stack:
            v, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                stack.pop()
                onpath.discard(path.pop())
                continue
            if nxt == s:
                cycles.append(list(path))
            elif nxt not in onpath:
                path.append(nxt)
                onpath.add(nxt)
                stack.append((nxt, iter([t for t in succ[nxt] if t in allowed and t >= s])))
    # remove duplicates caused by the two parallel edges? succ entries are distinct nodes
    uniq = {}
    for c in cycles:
        uniq[tuple(c)] = c
    return list(uniq.values())


def cycle_ap(cyc):
    a = sum(1 for v in cyc if v % 2 == 1)
    return a, len(cyc)


def expanding(a, p):
    return 3 ** a > 2 ** p


# ----------------------------------------------------------------------------- exact cycle densities
def karp_density(succ, nodes, maximize=True):
    """Exact max (or min) odd density a/p over the cycles inside `nodes`, returned as a pair
    (a, p) of integers realising it (a Karp witness ratio, reduced), or None if acyclic.
    Every walk of length j has weight a*log3 - j*log2, so maximising the weight of j-walks is
    maximising a: Karp's theorem runs on integer odd counts and is exact."""
    nodes = sorted(nodes)
    inset = set(nodes)
    n = len(nodes)
    D = [{v: 0 for v in nodes}]
    for j in range(1, n + 1):
        cur = {}
        for u, du in D[-1].items():
            val = du + (u & 1)
            for t in succ[u]:
                if t in inset:
                    if t not in cur or (val > cur[t] if maximize else val < cur[t]):
                        cur[t] = val
        D.append(cur)
    best = None
    for v in nodes:
        if v not in D[n]:
            continue
        worst = None
        for j in range(n):
            if v not in D[j]:
                continue
            num, den = D[n][v] - D[j][v], n - j
            if worst is None:
                worst = (num, den)
            else:
                # maximize: take min over j ; minimize: take max over j
                if (num * worst[1] < worst[0] * den) if maximize else (num * worst[1] > worst[0] * den):
                    worst = (num, den)
        if worst is None:
            continue
        if best is None or ((worst[0] * best[1] > best[0] * worst[1]) if maximize else (worst[0] * best[1] < best[0] * worst[1])):
            best = worst
    if best is None:
        return None
    g = math.gcd(best[0], best[1])
    return (best[0] // g, best[1] // g)


def density_vs_critical(ap):
    """sign of a/p - log_3 2 (exact): +1 expanding (3^a > 2^p), -1 contracting."""
    a, p = ap
    return 1 if 3 ** a > 2 ** p else -1


# ----------------------------------------------------------------------------- exact stationary law
def stationary(succ, comp):
    """exact stationary distribution (Fractions) of the uniform-edge chain on a bottom SCC."""
    idx = {v: i for i, v in enumerate(comp)}
    n = len(comp)
    # equations: pi_j = sum_i pi_i P_ij  ->  sum_i pi_i (P_ij - delta_ij) = 0 ; replace last by sum = 1
    A = [[Fraction(0)] * n for _ in range(n)]
    for i, v in enumerate(comp):
        for t in succ[v]:
            A[idx[t]][i] += Fraction(1, 2)   # row j = equation for pi_j, column i
        A[i][i] -= 1
    b = [Fraction(0)] * n
    A[n - 1] = [Fraction(1)] * n
    b[n - 1] = Fraction(1)
    # Gaussian elimination
    M = [row[:] + [bb] for row, bb in zip(A, b)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [x - f * y for x, y in zip(M[r], M[col])]
    pi = {v: M[idx[v]][n] for v in comp}
    return pi


def drift_sign_exact(pi_odd):
    """sign of pi_odd*log3 - log2, exactly: compare 3^p with 2^q for pi_odd = p/q."""
    p, q = pi_odd.numerator, pi_odd.denominator
    lhs, rhs = 3 ** p, 2 ** q
    return (lhs > rhs) - (lhs < rhs)


# ----------------------------------------------------------------------------- exceptional counts
def amin_table(Lmax):
    """amin[j] = least a with 3^a > 2^j (exact)."""
    out = []
    for j in range(Lmax + 1):
        a = 0
        while 3 ** a <= 2 ** j:
            a += 1
        out.append(a)
    return out


def bad_counts_dp(k, mask, Lmax):
    """|Bad_L| for L = 0..Lmax, by the Terras-style DP over (node, a)."""
    succ = graph(k, mask)
    K = 1 << k
    amin = amin_table(Lmax + 1)
    cur = {(s, 0): 1 for s in range(K)}
    counts = [K]
    for j in range(Lmax):
        nxt = {}
        for (s, a), c in cur.items():
            a2 = a + (s & 1)
            if a2 >= amin[j + 1]:
                for t in succ[s]:
                    nxt[(t, a2)] = nxt.get((t, a2), 0) + c
        cur = nxt
        counts.append(sum(cur.values()))
    return counts


def bad_counts_brute(k, mask, Lmax):
    """|Bad_L| by iterating T on every representative of Z/2^(k+L) (independent path)."""
    sig = sig_table(k, mask)
    out = []
    for L in range(Lmax + 1):
        M = 1 << (k + L)
        cnt = 0
        for x in range(M):
            n = x
            a = 0
            ok = True
            for j in range(1, L + 1):
                if n % 2:
                    a += 1
                n = T(n, k, sig)
                if 3 ** a < 2 ** j:
                    ok = False
                    break
            if ok:
                cnt += 1
        out.append(cnt)
    return out


# ----------------------------------------------------------------------------- (i) certificates
def descent_certificate(k, mask, Lmax=200):
    """If Bad_L = 0 for some L <= Lmax: return (L_min, n0, nclasses) where every n > n0 has
    T^j(n) < n at the first j <= L_min with 3^(a_j) < 2^j (class-determined), and nclasses is the
    number of stopping classes.  Uses a max-plus DP on the carry c (T^j(n) = (3^a n + c)/2^j)."""
    succ = graph(k, mask)
    sig = sig_table(k, mask)
    K = 1 << k
    # state (node, a) at depth j, all prefixes 1..j bad ; value = max carry c over such paths
    cur = {(s, 0): 0 for s in range(K)}
    n0 = 0
    nstop = 0
    ncount = {(s, 0): 1 for s in range(K)}
    for j in range(Lmax):
        nxt = {}
        nxtc = {}
        for (s, a), c in cur.items():
            if s & 1:
                a2 = a + 1
                c2 = 3 * c + sig[s] * (1 << j)
            else:
                a2 = a
                c2 = c
            den = (1 << (j + 1)) - 3 ** a2
            mult = ncount[(s, a)]
            if den > 0:
                # stopping: every class through this state stops here; threshold c2/den
                # (max over the paths -> max c2 gives the max threshold)
                thr = c2 // den if c2 >= 0 else -1
                n0 = max(n0, thr)
                nstop += mult
            else:
                for t in succ[s]:
                    key = (t, a2)
                    if key not in nxt or c2 > nxt[key]:
                        nxt[key] = c2
                    nxtc[key] = nxtc.get(key, 0) + mult
        cur = nxt
        ncount = nxtc
        if not cur:
            return j + 1, max(n0, 0), nstop
    return None


def cycles_below(k, mask, n0, extra=64):
    """all cycles met by orbits of 1..max(n0, extra); valid as the complete cycle list when the
    descent certificate with threshold n0 holds (every orbit then enters [1, n0])."""
    sig = sig_table(k, mask)
    top = max(n0, extra)
    found = {}
    known = {}
    for n in range(1, top + 1):
        seen = {}
        x = n
        i = 0
        while x not in seen and x not in known:
            seen[x] = i
            x = T(x, k, sig)
            i += 1
            if i > 10 ** 6:
                raise RuntimeError("orbit too long in finite check")
        if x in known:
            cid = known[x]
        else:
            # new cycle through x
            cyc = [x]
            y = T(x, k, sig)
            while y != x:
                cyc.append(y)
                y = T(y, k, sig)
            cid = min(cyc)
            found[cid] = cyc
            for y in cyc:
                known[y] = cid
        for y in seen:
            known[y] = cid
    return found


# ----------------------------------------------------------------------------- (ii) certificates
def divergence_certificate(k, mask, comp, Mmax=400):
    """For a bottom SCC comp: find the least M such that every path of length M inside comp has
    3^(a_M) >= 2^(M+1).  Then with W = max over such paths of max(0, -Off_M), every n in a class of
    comp with n > 2W satisfies T^(jM)(n) - 2W >= 2^j (n - 2W) -> infinity.
    Off_M = sum over odd steps i < M of sigma_i 2^i / 3^(a_(i+1)).  Returns (M, W) or None."""
    succ = graph(k, mask)
    sig = sig_table(k, mask)
    cs = set(comp)
    # min a over paths of length M from each node (within comp; comp is closed)
    amin = {s: 0 for s in comp}
    M_found = None
    for M in range(1, Mmax + 1):
        amin = {s: (s & 1) + min(amin[t] for t in succ[s]) for s in comp}
        if all(3 ** amin[s] >= 2 ** (M + 1) for s in comp):
            M_found = M
            break
    if M_found is None:
        return None
    M = M_found
    # W: max over paths of length M of max_{j<=M}(-Off_j)  (take the max over prefixes too)
    # DP over (node, a) at depth j: value = max of -Off_j ; -Off_{j+1} = -Off_j - sigma 2^j / 3^(a+1)
    cur = {(s, 0): Fraction(0) for s in comp}
    W = Fraction(0)
    for j in range(M):
        nxt = {}
        for (s, a), v in cur.items():
            if s & 1:
                a2 = a + 1
                v2 = v - Fraction(sig[s] * (1 << j), 3 ** a2)
            else:
                a2 = a
                v2 = v
            for t in succ[s]:
                key = (t, a2)
                if key not in nxt or v2 > nxt[key]:
                    nxt[key] = v2
        cur = nxt
        W = max(W, max(cur.values()))
    return M, W


# ----------------------------------------------------------------------------- cycle search
def cycle_search(k, mask, N, cap=10 ** 30, stepmax=100000):
    """positive starts 1..N.  Returns (cycles{min: (len, odd)}, escapes, unresolved)."""
    sig = sig_table(k, mask)
    fate = [0] * (N + 1)   # 0 unknown, >0 cycle id (min element), -1 escape, -2 unresolved
    cycles = {}
    esc = unres = 0
    for n in range(1, N + 1):
        if fate[n]:
            continue
        path = []
        seen = {}
        x = n
        res = None
        steps = 0
        while True:
            if x <= N and fate[x]:
                res = fate[x]
                break
            if x > cap:
                res = -1
                break
            if x in seen:
                # cycle among the path from seen[x]
                cyc = path[seen[x]:]
                m = min(cyc)
                a = sum(1 for y in cyc if y % 2)
                cycles[m] = (len(cyc), a)
                res = m
                break
            seen[x] = len(path)
            path.append(x)
            x = T(x, k, sig)
            steps += 1
            if steps > stepmax:
                res = -2
                break
        for y in path:
            if y <= N:
                fate[y] = res
    # per-start counts (every start whose fate is escape / unresolved)
    esc = sum(1 for n in range(1, N + 1) if fate[n] == -1)
    unres = sum(1 for n in range(1, N + 1) if fate[n] == -2)
    return cycles, esc, unres, fate


if __name__ == "__main__":
    # tiny smoke test
    for k in (1, 2):
        for mask in range(1 << (1 << (k - 1))):
            print(k, sigma_string(k, mask), ud_word(k, mask), graph(k, mask))
