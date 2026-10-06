#!/usr/bin/env python3
"""Independent audit (written from scratch) of Section 3 of
05-knowledge/results/chessboard_weave_20261006.md:
Prop 3.1 (integer cycles of T(n) = n/2, (3n+1)/2 and their shapes (p, A); Gersonides; Fibonacci leapers;
continued fraction of log_2 3; no cycle at shape (5,8)), Syracuse leap densities,
Prop 3.2 (3-adic valuation shells are single h-cycles; f maps into class 1 mod 3), plus a numerical look
at the spectrum of the walk P = (f* + h*)/2 cited from the five-papers note.
Exact integer arithmetic except the (clearly labelled) floating-point spectrum look.
"""
from math import gcd, comb
from collections import Counter
import itertools

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


# ---------- cycles met from |x0| <= 10^5 ----------
print("=== Prop 3.1: cycles of T on Z from |x0| <= 10^5 ===")
cycles = {}
known = {}  # element -> cycle id (min element)
for x0 in range(-10 ** 5, 10 ** 5 + 1):
    path = []
    seen = {}
    x = x0
    while x not in seen and x not in known:
        seen[x] = len(path)
        path.append(x)
        x = T(x)
    if x in known:
        cid = known[x]
    else:
        cyc = path[seen[x]:]
        cid = min(cyc)
        cycles[cid] = cyc
        for y in cyc:
            known[y] = cid
    for y in path:
        known.setdefault(y, cid)
for cid in sorted(cycles):
    cyc = cycles[cid]
    p = sum(1 for y in cyc if y % 2)
    A = len(cyc)  # each T step divides by 2 exactly once
    print(f"cycle min {cid}: length {len(cyc)}, elements {sorted(cyc) if len(cyc) < 12 else sorted(cyc)}, shape (p,A)=({p},{A}), |2^A-3^p|={abs(2**A-3**p)}")
shapes = sorted((sum(1 for y in c if y % 2), len(c)) for c in cycles.values())
check(shapes == [(0, 1), (1, 1), (1, 2), (2, 3), (7, 11)], "five shapes")

# ---------- Gersonides ----------
sols = [(p, A) for p in range(0, 400) for A in range(1, 700) if abs(2 ** A - 3 ** p) == 1]
print("solutions of |2^A - 3^p| = 1, A>=1, p<400, A<700:", sols)
check(sols == [(0, 1), (1, 1), (1, 2), (2, 3)], "Gersonides list")

# ---------- Fibonacci leapers ----------
F = [0, 1]
for _ in range(60):
    F.append(F[-1] + F[-2])
L = [2, 1]
for _ in range(60):
    L.append(L[-1] + L[-2])
v = (0, 1)
walk = [v]
for n in range(6):
    v = (v[1], v[0] + v[1])  # Q = [[0,1],[1,1]]
    walk.append(v)
print("Q^n (0,1)^T, n=0..6:", walk)
check(walk[:4] == [(0, 1), (1, 1), (1, 2), (2, 3)], "first four Fibonacci leapers")
print("(L_4, L_5) =", (L[4], L[5]), " (2,3)+(5,8) =", (7, 11), " (2,3)+(3,5) =", (5, 8))
check((L[4], L[5]) == (F[3] + F[5], F[4] + F[6]), "Lucas = F_{n-1}+F_{n+1}")


# ---------- continued fraction of log_2 3, exactly (Stern-Brocot with 2^A vs 3^p) ----------
def cf_log2_3(nterms):
    # x = log_2 3 ; compare A/p < x  <=>  2^A < 3^p
    def less(A, p):  # A/p < x ?
        return 2 ** A < 3 ** p
    terms = []
    # Stern-Brocot descent: left=(0,1), right=(1,0) in (num,den)=(A,p)
    ln, ld, rn, rd = 0, 1, 1, 0
    direction, run = None, 0
    while len(terms) < nterms:
        mn, md = ln + rn, ld + rd
        go_right = less(mn, md)  # mediant < x : move right
        d = 'R' if go_right else 'L'
        if d == direction:
            run += 1
        else:
            if direction is not None:
                terms.append(run)
            direction, run = d, 1
        if go_right:
            ln, ld = mn, md
        else:
            rn, rd = mn, md
    return terms


cf = cf_log2_3(9)
print("log_2 3 continued fraction (exact):", cf)
check(cf[:6] == [1, 1, 1, 2, 2, 3], "CF of log2 3")


def convergents(cf):
    h0, h1, k0, k1 = 1, cf[0], 0, 1
    out = [(h1, k1)]
    for a in cf[1:]:
        h0, h1 = h1, a * h1 + h0
        k0, k1 = k1, a * k1 + k0
        out.append((h1, k1))
    return out


print("convergents A/p of log_2 3 as (p,A):", [(k, h) for h, k in convergents(cf)[:6]])
print("convergents of phi as (p,A):", [(k, h) for h, k in convergents([1] * 8)[:7]])


# ---------- integer cycles of given shape via the cycle equation ----------
def integer_cycles_of_shape(p, A):
    D = 2 ** A - 3 ** p
    res = set()
    for ones in itertools.combinations(range(A), p):
        # parity vector eps_i, i=0..A-1; T^A x = (3^p x + r)/2^A with r = sum_{i: eps_i=1} 3^{#ones after i} 2^i
        r = 0
        cnt_after = p
        for i in range(A):
            if i in ones:
                cnt_after -= 1
                r += 3 ** cnt_after * 2 ** i
        if D != 0 and r % D == 0:
            x = r // D
            # verify the orbit realises the word
            y, ok = x, True
            for i in range(A):
                if (y % 2 == 1) != (i in ones):
                    ok = False
                    break
                y = T(y)
            if ok and y == x:
                res.add(x)
    return res


print("integer points on cycles of shape (5,8):", integer_cycles_of_shape(5, 8), " (2^8-3^5 =", 2 ** 8 - 3 ** 5, ")")
check(integer_cycles_of_shape(5, 8) == set(), "no (5,8) cycle")
print("integer points on cycles of shape (3,5):", integer_cycles_of_shape(3, 5), " (2^5-3^3 =", 2 ** 5 - 3 ** 3, ")")
allc = set()
for A in range(1, 19):
    for p in range(0, A + 1):
        allc |= integer_cycles_of_shape(p, A)
print("all integers lying on T-cycles with A <= 18 (any word):", sorted(allc))

# ---------- Syracuse leap densities ----------
print("=== Syracuse leaps ===")
K = 20
cnt = Counter()
for n in range(1, 2 ** K, 2):
    m = 3 * n + 1
    v = (m & -m).bit_length() - 1
    cnt[v] += 1
tot = 2 ** (K - 1)
print("v_2(3n+1) census over odd n < 2^20 (v: count, expected 2^(19-v)):", [(v, cnt[v], 2 ** (K - 1 - v)) for v in range(1, 8)])
check(all(cnt[v] == 2 ** (K - 1 - v) for v in range(1, K - 1)), "Geom(1/2) exact below 2^20 for v<=18")
# density of residues whose first k leaps are all v in {1,2}: exact over odd residues mod 2^(2k+1)
for k in range(1, 9):
    M = 2 ** (2 * k + 1)
    good = 0
    for n in range(1, M, 2):
        x, ok = n, True
        for _ in range(k):
            m = 3 * x + 1
            v = (m & -m).bit_length() - 1
            if v > 2:
                ok = False
                break
            x = m >> v
        good += ok
    check(good * 4 ** k == (M // 2) * 3 ** k, f"(3/4)^k density k={k}")
print("density of 'first k leaps in {1,2}' = (3/4)^k exactly over odd residues mod 2^(2k+1): checked k=1..8")

# ---------- Prop 3.2 ----------
print("=== Prop 3.2 ===")


def v3(x, k):
    if x == 0:
        return k
    j = 0
    while x % 3 == 0:
        x //= 3
        j += 1
    return j


for r in range(1, 16):
    M = 3 ** r
    order = 1
    y = 2 % M
    while y != 1:
        y = y * 2 % M
        order += 1
    check(order == 2 * 3 ** (r - 1), f"2 primitive root mod 3^{r}")
for k in range(1, 9):
    M = 3 ** k
    inv2 = pow(2, -1, M)
    h = lambda x: x * inv2 % M
    for j in range(0, k + 1):
        Vj = [x for x in range(M) if v3(x, k) == j]
        size_ok = (len(Vj) == (2 * 3 ** (k - 1 - j) if j < k else 1))
        check(size_ok, f"|V_{j}| k={k}")
        check(all(v3(h(x), k) == j for x in Vj), f"h preserves V_{j} k={k}")
        # single cycle
        x0 = Vj[0]
        x, L = h(x0), 1
        while x != x0:
            x = h(x)
            L += 1
        check(L == len(Vj), f"V_{j} single h-cycle k={k}")
    check(all((3 * x + 1) % M % 3 == 1 for x in range(M)), f"f into class 1 mod 3 k={k}")
print("Prop 3.2 checked k=1..8 (and 2 primitive root mod 3^r, r<=15)")

# ---------- look at the cited spectrum (floating point, informational only) ----------
try:
    import numpy as np
    print("--- informational: spectrum of P g(x) = (g(3x+1) + g(x/2))/2 ---")
    for k in range(1, 6):
        M = 3 ** k
        inv2 = pow(2, -1, M)
        for label, dom in (("Z/3^k", list(range(M))), ("units", [x for x in range(M) if x % 3])):
            idx = {x: i for i, x in enumerate(dom)}
            P = np.zeros((len(dom), len(dom)))
            for x in dom:
                P[idx[x], idx[(3 * x + 1) % M]] += 0.5
                P[idx[x], idx[x * inv2 % M]] += 0.5
            ev = np.linalg.eigvals(P)
            mods = Counter(round(abs(e), 6) for e in ev)
            print(f"k={k} {label}: |lambda| multiset {dict(sorted(mods.items()))}; real eigenvalues {sorted(set(round(e.real, 6) for e in ev if abs(e.imag) < 1e-9))}")
except Exception as e:
    print("numpy spectrum skipped:", e)

print()
print("TOTAL FAILURES:", len(fails))
