#!/usr/bin/env python3
"""Auditor 'citations' (S91): independent re-derivation of the counts quoted from J. Allikvere, arXiv:2607.28813,
Part III (king strip P_3 x P_n), and of the arithmetic the S91 README attaches to them.

Written from the PAPER'S definitions only (Lemmas 14.1, 14.2, Section 15, Theorem 16.1, Remark 16.3, Appendix A);
no code of the audited session is imported or read by this script.

Paper's data used as targets:
  counts of weight-n admissible sequences: 188 (n = 5), 556 (n = 6), 1660 (n = 7), 5004 (n = 8), 44,972 (n = 10);
  184 tight transitions; 320 minimum-weight 4-prefixes; 12 accepting states; 19-state recurrent core L;
  106 minimum 4-prefixes ending in L; eigenvector w (A w = 3 w) with the entries of Appendix A.
README claims checked:
  charpoly(core) = x^14 (x-3)(x-1)(x+1)(x^2+x+2);
  a(n) = 2a(n-1) + a(n-2) + 6a(n-3);
  42 a(n) = 32*3^n - 9 V_n - 21 U_n, V_n^2 + 7 U_n^2 = 2^(n+2);
  |U_n| = 1 exactly at n = 1, 2, 3, 5, 13 (Ramanujan-Nagell).
"""
import itertools
import sys

import sympy

OK = True


def check(name, cond, detail=""):
    global OK
    OK = OK and bool(cond)
    print(("PASS  " if cond else "FAIL  ") + name + ("   [" + str(detail) + "]" if detail != "" else ""))


COLS = list(itertools.product((0, 1), repeat=3))  # (p, m, q): rows 0, 1, 2
Z = (0, 0, 0)
INF = float("inf")


def ok_mid(u, v, w):
    """(A), (B), (C) of the paper at the middle column v of (u, v, w)."""
    pu, mu, qu = u
    pv, mv, qv = v
    pw, mw, qw = w
    a = (pv, pu + pw) != (qv, qu + qw)
    b = (pv != mv) or (qu + qv + qw >= 1)
    c = (qv != mv) or (pu + pv + pw >= 1)
    return a and b and c


wt = {t: sum(t) for t in COLS}


# ---------------------------------------------------------------- 1. counts by a weight-resolved dynamic programme
def count_weight_n(n):
    """number of admissible sequences t_1..t_n (t_0 = t_(n+1) = 000) of total weight exactly n."""
    # state: (t_(i-1), t_i) -> dict weight -> count
    cur = {}
    for t1 in COLS:
        if wt[t1] <= n:
            cur.setdefault((Z, t1), {}).setdefault(wt[t1], 0)
            cur[(Z, t1)][wt[t1]] += 1
    for i in range(1, n):  # choose t_(i+1), check the conditions at column i
        nxt = {}
        for (u, v), d in cur.items():
            for w in COLS:
                if not ok_mid(u, v, w):
                    continue
                for wgt, c in d.items():
                    nw = wgt + wt[w]
                    if nw <= n:
                        e = nxt.setdefault((v, w), {})
                        e[nw] = e.get(nw, 0) + c
        cur = nxt
    tot = 0
    below = 0
    for (u, v), d in cur.items():
        if ok_mid(u, v, Z):
            tot += d.get(n, 0)
            below += sum(c for wgt, c in d.items() if wgt < n)
    return tot, below


a = {}
low = {}
for n in range(1, 41):
    a[n], low[n] = count_weight_n(n)
print("a(1..12) =", [a[n] for n in range(1, 13)])
check("paper: a(5) = 188, a(6) = 556, a(7) = 1660, a(8) = 5004, a(10) = 44972",
      (a[5], a[6], a[7], a[8], a[10]) == (188, 556, 1660, 5004, 44972), (a[5], a[6], a[7], a[8], a[10]))
check("paper, Prop. 15.1: no admissible sequence of length n >= 3 has weight < n (n <= 40)",
      all(low[n] == 0 for n in range(3, 41)), [n for n in range(1, 41) if low[n]])

# brute force for small n (independent of the DP)
def brute(n):
    c = 0
    for seq in itertools.product(COLS, repeat=n):
        if sum(wt[t] for t in seq) != n:
            continue
        ext = (Z,) + seq + (Z,)
        if all(ok_mid(ext[i - 1], ext[i], ext[i + 1]) for i in range(1, n + 1)):
            c += 1
    return c


check("brute force agrees with the DP for n = 1..6", all(brute(n) == a[n] for n in range(1, 7)),
      [brute(n) for n in range(1, 7)])

# ---------------------------------------------------------------- 2. README recurrence and closed form
bad = [n for n in range(4, 41) if a[n] != 2 * a[n - 1] + a[n - 2] + 6 * a[n - 3]]
check("README: a(n) = 2a(n-1) + a(n-2) + 6a(n-3) for 4 <= n <= 40", not bad, bad)
U = {0: 0, 1: 1}
V = {0: 2, 1: -1}
for n in range(2, 20001):
    U[n] = -U[n - 1] - 2 * U[n - 2]
    V[n] = -V[n - 1] - 2 * V[n - 2]
check("Lucas pair of x^2 + x + 2: V_n^2 + 7 U_n^2 = 2^(n+2) for n <= 3000",
      all(V[n] ** 2 + 7 * U[n] ** 2 == 2 ** (n + 2) for n in range(0, 3001)))
badc = [n for n in range(1, 41) if 42 * a[n] != 32 * 3 ** n - 9 * V[n] - 21 * U[n]]
check("README: 42 a(n) = 32*3^n - 9 V_n - 21 U_n for 1 <= n <= 40", not badc, badc)
ones = [n for n in range(1, 20001) if abs(U[n]) == 1]
check("README: |U_n| = 1 exactly at n = 1, 2, 3, 5, 13 for n <= 20000", ones == [1, 2, 3, 5, 13], ones)
check("Ramanujan-Nagell list: |V_n| at those n is 1, 3, 5, 11, 181 and x^2 + 7 = 2^m at m = 3, 4, 5, 7, 15",
      [abs(V[n]) for n in ones] == [1, 3, 5, 11, 181] and all(abs(V[n]) ** 2 + 7 == 2 ** (n + 2) for n in ones)
      and [n + 2 for n in ones] == [3, 4, 5, 7, 15])
# direct search of x^2 + 7 = 2^m for m <= 4000 (exact integer square roots)
from math import isqrt
sol = [m for m in range(3, 4001) if isqrt(2 ** m - 7) ** 2 == 2 ** m - 7]
check("direct: x^2 + 7 = 2^m has solutions exactly at m = 3, 4, 5, 7, 15 for m <= 4000", sol == [3, 4, 5, 7, 15], sol)
x = sympy.symbols("x")
check("x^3 - 2x^2 - x - 6 = (x - 3)(x^2 + x + 2)", sympy.expand((x - 3) * (x ** 2 + x + 2) - (x ** 3 - 2 * x ** 2 - x - 6)) == 0)

# ---------------------------------------------------------------- 3. the automaton of Theorem 16.1 and its core
STATES = [(u, v) for u in COLS for v in COLS]
V1 = {(u, v): (wt[v] if u == Z else INF) for (u, v) in STATES}


def T(Vv):
    out = {}
    for (v, w) in STATES:
        best = INF
        for u in COLS:
            if Vv[(u, v)] < INF and ok_mid(u, v, w):
                best = min(best, Vv[(u, v)] + wt[w])
        out[(v, w)] = best
    return out


V2 = T(V1)
V3 = T(V2)
V4 = T(V3)
V5 = T(V4)
check("paper: V_5 = V_4 + 1 on all 64 states", all((V5[s] == V4[s] + 1) or (V5[s] == INF and V4[s] == INF) for s in STATES))
fin = [s for s in STATES if V4[s] < INF]
check("paper: V_4 finite at 62 states; infinite at (000,000) and (000,010)",
      len(fin) == 62 and sorted(s for s in STATES if V4[s] == INF) == [(Z, Z), (Z, (0, 1, 0))])
# Appendix A table (u v V4), transcribed from the paper
APP = """000 001 4 010 010 4 100 010 4 110 010 5
000 011 5 010 011 5 100 011 5 110 011 6
000 100 4 010 100 4 100 100 4 110 100 5
000 101 4 010 101 5 100 101 5 110 101 6
000 110 5 010 110 5 100 110 5 110 110 6
000 111 5 010 111 6 100 111 6 110 111 7
001 000 3 011 000 4 101 000 4 111 000 5
001 001 4 011 001 5 101 001 4 111 001 5
001 010 4 011 010 5 101 010 5 111 010 6
001 011 5 011 011 6 101 011 5 111 011 6
001 100 4 011 100 5 101 100 4 111 100 5
001 101 5 011 101 6 101 101 6 111 101 7
001 110 5 011 110 6 101 110 5 111 110 6
001 111 6 011 111 7 101 111 7 111 111 8
010 000 3 100 000 3 110 000 4
010 001 4 100 001 4 110 001 5""".split()
tup = lambda s: tuple(int(ch) for ch in s)
app = {(tup(APP[i]), tup(APP[i + 1])): int(APP[i + 2]) for i in range(0, len(APP), 3)}
check("paper, Appendix A: the 62 printed values of V_4 equal the recomputed ones",
      len(app) == 62 and all(V4[s] == app[s] for s in app))

tight = [((u, v), (v, w)) for (u, v) in fin for w in COLS
         if ok_mid(u, v, w) and V4[(v, w)] < INF and wt[w] + V4[(u, v)] - V4[(v, w)] == 1]
check("paper: 184 tight transitions", len(tight) == 184, len(tight))

# minimum-weight 4-prefixes: admissible prefixes t1..t4 (conditions at columns 1..3) of weight V4(t3, t4)
pref = []
for seq in itertools.product(COLS, repeat=4):
    ext = (Z,) + seq
    if all(ok_mid(ext[i - 1], ext[i], ext[i + 1]) for i in range(1, 4)):
        if sum(wt[t] for t in seq) == V4[(seq[2], seq[3])]:
            pref.append(seq)
check("paper: 320 minimum-weight 4-prefixes", len(pref) == 320, len(pref))
acc = [s for s in fin if V4[s] == 4 and ok_mid(s[0], s[1], Z)]
check("paper: 12 accepting states (V_4 = 4, legal last column)", len(acc) == 12, len(acc))

succ = {}
pred = {}
for s, t in tight:
    succ.setdefault(s, set()).add(t)
    pred.setdefault(t, set()).add(s)


def closure(start, nb):
    seen = set(start)
    stack = list(start)
    while stack:
        s = stack.pop()
        for t in nb.get(s, ()):
            if t not in seen:
                seen.add(t)
                stack.append(t)
    return seen


entry = {(p[2], p[3]) for p in pref}
reach = closure(entry, succ)
coreach = closure(acc, pred)
L = sorted(reach & coreach)
check("paper: the core L has 19 states", len(L) == 19, len(L))
check("paper: second components of L are 000, 001, 010, 100, 101",
      sorted({s[1] for s in L}) == sorted([Z, (0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 0, 1)]))
check("paper: all 12 accepting states lie in L", all(s in L for s in acc))
check("paper: 106 minimum 4-prefixes end in L", sum(1 for p in pref if (p[2], p[3]) in L) == 106,
      sum(1 for p in pref if (p[2], p[3]) in L))
idx = {s: i for i, s in enumerate(L)}
A = sympy.zeros(19, 19)
for s, t in tight:
    if s in idx and t in idx:
        A[idx[s], idx[t]] += 1
# strong connectivity inside L
succL = {s: {t for t in succ.get(s, ()) if t in idx} for s in L}
predL = {s: {t for t in pred.get(s, ()) if t in idx} for s in L}
check("paper: L is strongly connected and has loops",
      closure([L[0]], succL) == set(L) and closure([L[0]], predL) == set(L) and any(A[i, i] for i in range(19)))
WAPP = """000 001 9 001 001 9 010 010 6 100 101 2
000 100 9 001 010 5 010 100 9 101 000 6
000 101 6 001 100 9 100 000 2 101 001 9
001 000 2 001 101 2 100 001 9 101 100 9
010 001 9 100 010 5 100 100 9""".split()
wapp = {(tup(WAPP[i]), tup(WAPP[i + 1])): int(WAPP[i + 2]) for i in range(0, len(WAPP), 3)}
check("paper, Appendix A: the 19 printed states are L", sorted(wapp) == L)
wv = sympy.Matrix([wapp[s] for s in L])
check("paper, Appendix A: A w = 3 w for the printed vector w", A * wv == 3 * wv)
cp = sympy.factor(A.charpoly(x).as_expr())
target = x ** 14 * (x - 3) * (x - 1) * (x + 1) * (x ** 2 + x + 2)
check("README: charpoly of the core = x^14 (x-3)(x-1)(x+1)(x^2+x+2)", sympy.expand(cp - target) == 0, cp)

# automaton path counts = a(n) (Theorem 16.1), n = 5..40
ent = sympy.zeros(1, 19)
for p in pref:
    s = (p[2], p[3])
    if s in idx:
        ent[0, idx[s]] += 1
fin_vec = sympy.Matrix([1 if s in acc else 0 for s in L])
vec = ent
paths = {4: (vec * fin_vec)[0, 0]}
for n in range(5, 41):
    vec = vec * A
    paths[n] = (vec * fin_vec)[0, 0]
check("paper, Thm 16.1 / Remark 16.3: automaton path counts equal a(n) for 5 <= n <= 40",
      all(paths[n] == a[n] for n in range(5, 41)), "n = 4: automaton %s, direct %s" % (paths[4], a[4]))
# does the +-1 part of the spectrum contribute?  minimal recurrence of a(n)
seq = [a[n] for n in range(1, 41)]
r3 = all(seq[i] == 2 * seq[i - 1] + seq[i - 2] + 6 * seq[i - 3] for i in range(3, 40))
print("order-3 recurrence holds from n = 4 on:", r3, "| first values a(1..4) =", seq[:4])
print("a(n)/3^n at n = 40:", float(sympy.Rational(a[40], 3 ** 40)), " (16/21 = %.6f)" % (16 / 21))
print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
sys.exit(0 if OK else 1)
