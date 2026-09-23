#!/usr/bin/env python3
"""collatz_procgen_20260922_mirror_escape.py -- the Q1 escape near -1 (mirror of loops_escapes.py).

Setting.  n = 2^m u - 1, u odd, m >= 2 (the -1 thread at precision exactly m).
(A) A*(m) = least number of multiplications A over E-paths from -1 with exactly m halvings that end
    (after trailing multiplications) at an ODD value v; exit price P(m) = 3^{A*}/2^{m+1}.
    Exact: layered search from -1 with the safe pruning  a <= Amax - (m-b) log_3 2 + log_3|w|  (a state w
    after b halvings needs a'' >= (m-b) log_3 2 - log_3|w| more multiplications), Amax = an explicit
    upper bound (the path (MH)^m has A = m; an a0-loop has A = a0(m)).
(B) Exact integer checks of the exit lemma on random n: every state of every E-path with <= m halvings
    (enumerated with a <= a0(m)+2) exceeds n; the least value at the (m+1)-th halving is
    (3^{A*} u - |v|)/2 > n; the canonical a0-loop exit lands exactly at (3^{A*} u - 1)/2.
(C) Landing classes: y = (3^A u - 1)/2 is a bijection u (odd mod 2^{r+1}) <-> y mod 2^r; landing precision
    m' = v_2(3^A u + 1) - 1: for u = 1,3 mod 8, m' <= 1; for u = 5,7 mod 8, m' = v_2(A - lambda(u)) + 1 if
    A = lambda(u) mod 2 (lambda(u) = 2-adic log_3(-1/u)), else m' = 0.
(D) Chain recursion u_{t+1} = (3^{A*(m_t)} u_t + 1)/2^{m_{t+1}+1} and chain statistics on integers.
Usage: python3 ..._mirror_escape.py [MMAX]
"""
import sys, random, math, collections
from fractions import Fraction as F

TH = math.log(2) / math.log(3)

def a0(K):
    a = 0
    while 3 ** a <= 2 ** (K + 1): a += 1
    return a

def log3abs(w): return math.log(abs(w)) / math.log(3)

def endpoints(m, Amax):
    """all (v odd, A) reachable from -1 by paths with exactly m halvings (then trailing multiplications
    until odd), with A <= Amax; exact thanks to the safe pruning."""
    layer = {-1: 0}
    for b in range(m):
        nxt = {}
        for v, a in layer.items():
            x, aa = v, a
            if x % 2: x, aa = 3 * x + 1, aa + 1
            while aa <= Amax:
                w = x // 2
                if aa <= Amax - (m - b - 1) * TH + log3abs(w) + 1e-9 and aa < nxt.get(w, 10 ** 9): nxt[w] = aa
                x, aa = 9 * x + 4, aa + 2
        layer = nxt
    out = {}
    for w, a in layer.items():
        v, A = (w, a) if w % 2 else (3 * w + 1, a + 1)
        if A <= Amax and A < out.get(v, 10 ** 9): out[v] = A
    return out

def Astar(m):
    Amax = max(m, a0(m)) if m <= 9 else a0(m)
    ep = endpoints(m, Amax)
    A = min(ep.values())
    return A, sorted(v for v, a in ep.items() if a == A), ep

def partA(MMAX):
    print(f"(A) exit prices at the (m+1)-th halving, m = 2..{MMAX}  [A* exact; P(m) = 3^A*/2^(m+1)]")
    rows = []
    for m in range(2, MMAX + 1):
        A, vs, ep = Astar(m)
        # the bound 3^A >= 2^(m+1)|v| + 1 on every endpoint
        assert all(3 ** a >= 2 ** (m + 1) * abs(v) + 1 for v, a in ep.items())
        P = F(3 ** A, 2 ** (m + 1)); eta = a0(m) - (m + 1) * TH
        rows.append((m, A, P))
        print(f"   m={m:3d} a0(m)={a0(m):3d} A*={A:3d} {'=a0' if A == a0(m) else '>a0'} P={float(P):9.5f}"
              f"  3^eta(m+1)={3 ** eta:8.5f}  lnP/(m+1)={math.log(P)/(m+1):.4f}  endpoints v at A*: {vs[:4]}")
    return rows

def partB(rows, ntests=400, seed=1):
    print("(B) exit lemma on random integers n = 2^m u - 1 (u odd, 2 <= m <= 14):")
    rng = random.Random(seed); ok1 = ok2 = ok3 = True
    Ainfo = {m: Astar(m) for m in range(2, 15)}
    for t in range(ntests):
        m = rng.randint(2, 14); u = 2 * rng.randint(0, 10 ** 6) + 1; n = 2 ** m * u - 1
        lim = a0(m) + 2 if m >= 10 else m + 2
        stack = [(n, 0, 0)]
        while stack:                       # all E-path states with <= m halvings and a <= lim
            x, a, b = stack.pop()
            if (a, b) != (0, 0) and x <= n: ok1 = False
            if a < lim: stack.append((3 * x + 1, a + 1, b))
            if x % 2 == 0 and b < m: stack.append((x // 2, a, b + 1))
        A, vs, ep = Ainfo[m]
        least = min((vv + 3 ** aa * u) // 2 for vv, aa in ep.items())   # value at the (m+1)-th halving
        if least <= n: ok2 = False
        vmax = min(vs)                      # most negative A*-endpoint
        if least != (3 ** A * u + vmax) // 2: ok3 = False
        if vs == [-1] and least != (3 ** A * u - 1) // 2: ok3 = False
    print(f"   all states with <= m halvings exceed n: {ok1}; least value at the (m+1)-th halving exceeds n: {ok2};"
          f" least value = (3^A* u - |v_min|)/2 over the A*-endpoints, = (3^A* u - 1)/2 when the only one is -1: {ok3}")

def v2(x):
    x = abs(x); c = 0
    while x % 2 == 0: x //= 2; c += 1
    return c

def partC(r=14):
    print("(C) landing classes y = (3^A u - 1)/2:")
    for A in (7, 12, 53):
        ys = set(((3 ** A * u - 1) // 2) % 2 ** r for u in range(1, 2 ** (r + 1), 2))
        print(f"   A={A}: u odd mod 2^{r+1} -> y mod 2^{r}: {len(ys)} classes (bijection: {len(ys) == 2 ** r})")
    k = 16; mod = 2 ** (k + 2); bad = cnt = skipped = 0
    for u in range(1, 3001, 2):
        target = (-pow(u, -1, mod)) % mod
        lam = None
        if u % 8 in (5, 7):
            # 2-adic log: 3^lam = -1/u mod 2^(k+2), lam mod 2^k
            x3 = 1
            for x in range(2 ** k):
                if x3 == target: lam = x; break
                x3 = x3 * 3 % mod
        for A in range(1, 80):
            mp = v2(3 ** A * u + 1) - 1
            if u % 8 in (1, 3):
                pred_ok = mp <= 1
            else:
                d = A - lam
                if d % 2 == 0 and v2(d) >= k - 2: skipped += 1; continue
                pred = v2(d) + 1 if d % 2 == 0 else 0
                pred_ok = (mp == pred)
            cnt += 1; bad += (not pred_ok)
    print(f"   landing precision vs the 2-adic-log formula on {cnt} pairs (odd u < 3000, A < 80): mismatches {bad}"
          f" (skipped {skipped} beyond the {k}-bit resolution of lambda)")

def chain_stats(rows, NBITS=40, samples=200000, EXH=22, seed=7):
    """(D) canonical chains: n = 2^m u - 1 -> y = (3^{A*(m)} u - 1)/2; repeat while y = 2^{m'} u' - 1, m' >= 2."""
    Acache = {m: A for m, A, P in rows}
    def Aof(m): return Acache[m] if m in Acache else a0(m)   # A*(m) = a0(m) for m >= 10 (loop DP)
    def run(n):
        m = v2(n + 1); u = (n + 1) >> m; links = 0; logP = 0.0; ms = [m]
        while m >= 2:
            A = Aof(m); y = (3 ** A * u - 1) // 2
            logP += A * math.log(3) - (m + 1) * math.log(2); links += 1
            m = v2(y + 1); u = (y + 1) >> m; ms.append(m)
        return links, logP, ms
    print(f"(D) canonical exit chains (repeated landings on the -1 thread with precision >= 2):")
    rng = random.Random(seed); hist = collections.Counter(); maxl = (0,); maxlp = (0,)
    for _ in range(samples):
        n = (rng.randrange(1 << (NBITS - 3)) << 3) | 7
        l, lp, ms = run(n); hist[l] += 1
        if l > maxl[0]: maxl = (l, n, ms)
        if lp > maxlp[0]: maxlp = (lp, n, ms)
    tot = sum(hist.values())
    print(f"   random n = 7 mod 8 below 2^{NBITS} ({samples}): length histogram {dict(sorted(hist.items()))}")
    print("   P(length >= k): " + ", ".join(f"{k}:{sum(v for l,v in hist.items() if l>=k)/tot:.4f}" for k in range(1, 9)))
    print(f"   longest {maxl[0]} (n={maxl[1]}, precisions {maxl[2]}); largest cumulative price {math.exp(maxlp[0]):.2f}"
          f" (n={maxlp[1]}, precisions {maxlp[2]})")
    hist2 = collections.Counter(); best = (0,); bestl = (0,)
    for n in range(7, 1 << EXH, 8):
        l, lp, ms = run(n); hist2[l] += 1
        if lp > best[0]: best = (lp, n, ms)
        if l > bestl[0]: bestl = (l, n, ms)
    print(f"   all n = 7 mod 8 below 2^{EXH}: {dict(sorted(hist2.items()))}; longest {bestl[0]} at n={bestl[1]} {bestl[2]};"
          f" max cumulative price {math.exp(best[0]):.2f} at n={best[1]} ({best[2]}), i.e. n^{best[0]/math.log(best[1]):.3f}")

if __name__ == '__main__':
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rows = partA(MMAX)
    print("   seesaw argument list:", " ".join(f"{m}:{A}" for m, A, P in rows))
    partB(rows)
    partC()
    chain_stats(rows)
