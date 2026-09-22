#!/usr/bin/env python3
"""Independent audit of lane minus_sheet_positive_control (session collatz-mod6-20260922).

Recomputes every key number of collatz_mod6_20260922_minus_sheet_positive_control.out
with separately written code (different loop structure, pure-Python cross-checks
where affordable), and adds the audit refinements:
  A0  entropy of 8/pi^2
  A1  T_- basin census to 10^7 (numpy, sequential root resolution) and a pure-Python
      memoised census to 10^6 as a second method; plus-sheet control
  A2  Terras prefix-descent counts (DP on words) versus direct residue counts, J <= 20
  A3  sigma vs sigma_c to 2*10^5 (pure Python), residue-level strict-early-descent
      check at K = 20 (independent numpy code)
  A4  greedy G_- to 10^6: statuses, histogram head, peak; the one-step law is EXACT
      (m certified in one greedy step iff m mod 9 in {1,2,5,7}); BFS rescue of m = 4;
      arrow count of the SCC note's compound rescue; forward BFS distances from 1
  A5  Berggren table readings and the corrected coupling breakdown; the k=1 <-> k=2
      B3 bijection with the mod-8 characterisation of k=2 sources; the y = x edge case
  A6  budget sums and cycle gates
All checks raise explicitly.
"""
import math
import sys
from collections import deque
from fractions import Fraction

import numpy as np


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def v2(x):
    return (x & -x).bit_length() - 1


print("=== A0: entropy ===")
p = 8 / math.pi ** 2
H = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
print(f"8/pi^2 = {p:.6f}  H2 = {H:.5f} bits; |H-0.704| = {abs(H - 0.704):.5f}")
check(abs(H - 0.70028) < 1e-5 and abs(H - 0.704) > 3e-3, "entropy")

# ----------------------------------------------------------------------------
print()
print("=== A1: basin census (independent code) ===")
XMAX = 10 ** 7


def first_below(X, sign, cap=6000, tform=True):
    """fb[n] = first orbit value < n (0 if none within cap); also the max first-descent
    time and max value seen, and the list of non-descenders.  tform=True: odd n ->
    (3n+b)/2; tform=False: the C-form odd n -> 3n+b (the explorer's S1 iterates this)."""
    n = np.arange(X + 1, dtype=np.int64)
    fb = np.zeros(X + 1, dtype=np.int64)
    idx = n[1:]
    cur = idx.copy()
    tmax = 0
    vmax = 0
    for t in range(1, cap + 1):
        odd = (cur & 1) == 1
        cur = np.where(odd, (3 * cur + sign) >> (1 if tform else 0), cur >> 1)
        vmax = max(vmax, int(cur.max()))
        check(vmax < (1 << 61), "overflow")
        d = cur < idx
        if d.any():
            fb[idx[d]] = cur[d]
            tmax = t
        idx = idx[~d]
        cur = cur[~d]
        if idx.size == 0:
            break
    return fb, tmax, vmax, idx.tolist()


_, tmax_mT, vmax_mT, never_mT = first_below(XMAX, -1, tform=True)
print(f"minus T-form (odd n -> (3n-1)/2): max first-descent time {tmax_mT}, max value {vmax_mT}, never-descending {never_mT}")
fb, tmax_m, vmax_m, never_m = first_below(XMAX, -1, tform=False)
print(f"minus C-form (odd n -> 3n-1): max first-descent time {tmax_m}, max value {vmax_m}, never-descending {never_m}")
print("DEFINITION MISMATCH: the explorer's S1 census iterates the C-form; its 445 and 30541433029400 are C-form numbers, "
      f"the T-form gives {tmax_mT} and {vmax_mT} = {vmax_m}/2 (basins identical)")
check(tmax_m == 445 and vmax_m == 30541433029400 and never_m == [1, 5, 17], "minus census header (C-form)")
check(tmax_mT == 273 and 2 * vmax_mT == vmax_m and never_mT == [1, 5, 17], "minus census header (T-form)")
root = np.zeros(XMAX + 1, dtype=np.int64)
fbl = fb.tolist()
rl = [0] * (XMAX + 1)
for m in range(1, XMAX + 1):
    f = fbl[m]
    rl[m] = m if f == 0 else rl[f]
root = np.array(rl, dtype=np.int64)
del rl, fbl
expected = {10 ** 4: ((3244, 1605), (3213, 1623), (3543, 1772)),
            10 ** 5: ((33030, 16553), (32104, 16026), (34866, 17421)),
            10 ** 6: ((327679, 163486), (323351, 162122), (348970, 174392)),
            10 ** 7: ((3273791, 1636054), (3244985, 1623149), (3481224, 1740797))}
for X in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    sub = root[1:X + 1]
    subo = root[1:X + 1:2]
    got = tuple((int((sub == r).sum()), int((subo == r).sum())) for r in (1, 5, 17))
    print(f"X={X}: basins (all, odd) for 1,5,17: {got}")
    check(got == expected[X], f"basin counts at X={X}")
del root, fb
# second method, pure Python memoised, n <= 10^6
X2 = 10 ** 6
rt = [0] * (X2 + 1)
for m in range(1, X2 + 1):
    if m in (1, 5, 17):
        rt[m] = m
        continue
    x = m
    while x >= m:
        x = (3 * x - 1) >> 1 if x & 1 else x >> 1
    rt[m] = rt[x]
cnt = tuple((sum(1 for m in range(1, X2 + 1) if rt[m] == r),
             sum(1 for m in range(1, X2 + 1, 2) if rt[m] == r)) for r in (1, 5, 17))
print(f"pure-Python memoised census to 10^6: {cnt}")
check(cnt == expected[10 ** 6], "second-method census at 10^6")
del rt
_, tmax_pT, vmax_pT, never_pT = first_below(XMAX, +1, tform=True)
print(f"plus T-form: max first-descent time {tmax_pT}, max value {vmax_pT}, never-descending {never_pT}")
fbp, tmax_p, vmax_p, never_p = first_below(XMAX, +1, tform=False)
print(f"plus C-form: max first-descent time {tmax_p}, max value {vmax_p}, never-descending {never_p}")
check(tmax_p == 401 and vmax_p == 60342610919632 and never_p == [1] and 2 * vmax_pT == vmax_p and never_pT == [1], "plus census header")
del fbp

# ----------------------------------------------------------------------------
print()
print("=== A2: Terras prefix-descent counts, DP on words vs direct residues ===")
JMAX = 20


def no_descent_words(J):
    dp = {0: 1}
    for i in range(1, J + 1):
        nd = {}
        for a, c in dp.items():
            for s in (0, 1):
                if 3 ** (a + s) >= 2 ** i:
                    nd[a + s] = nd.get(a + s, 0) + c
        dp = nd
    return sum(dp.values())


def residue_descend_count(J, sign):
    r = np.arange(1 << J, dtype=np.int64)
    cur = r.copy()
    a = np.zeros(1 << J, dtype=np.int64)
    ok = np.zeros(1 << J, dtype=bool)
    for i in range(1, J + 1):
        odd = (cur & 1) == 1
        a += odd
        cur = np.where(odd, (3 * cur + sign) >> 1, cur >> 1)
        ok |= (3.0 ** a) < float(2 ** i)
    return int(ok.sum())


exp_desc = [1, 3, 6, 13, 28, 56, 115, 237, 474, 960, 1920, 3870, 7825, 15650, 31473, 63422, 126844, 254649, 509298, 1021248]
exp_nd = [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328]
for J in range(1, JMAX + 1):
    nd = no_descent_words(J)
    cp = residue_descend_count(J, +1)
    cm = residue_descend_count(J, -1)
    check(cp == cm == (1 << J) - nd == exp_desc[J - 1] and nd == exp_nd[J - 1], f"J={J}")
print(f"J=1..20: DP no-descent words {exp_nd}")
print(f"J=1..20: descend counts (both sheets equal, = 2^J - DP) {exp_desc}; fraction at J=20 = {exp_desc[-1] / 2 ** 20:.6f}")

# ----------------------------------------------------------------------------
print()
print("=== A3: sigma vs sigma_c (pure Python to 2*10^5), residue check K=20 ===")
NS = 2 * 10 ** 5
for sign, name in ((+1, "plus"), (-1, "minus")):
    inf = []
    lt = gt = 0
    for n in range(1, NS + 1):
        x = n
        a = 0
        sa = sc = None
        for k in range(1, 500):
            if x & 1:
                x = (3 * x + sign) >> 1
                a += 1
            else:
                x >>= 1
            if sa is None and x < n:
                sa = k
            if sc is None and 3 ** a < 2 ** k:
                sc = k
            if sa is not None and sc is not None:
                break
        if sa is None:
            inf.append(n)
        elif sa < sc:
            lt += 1
        elif sa > sc:
            gt += 1
    print(f"{name}: infinite sigma {inf}; sigma<sigma_c: {lt}; sigma>sigma_c: {gt}")
    check(lt == 0 and gt == 0 and inf == ([1] if sign == 1 else [1, 5, 17]), name + " stopping times")
K = 20
r = np.arange(1, 1 << K, dtype=np.int64)
cur = r.copy()
a = np.zeros(r.size, dtype=np.int64)
ta = np.full(r.size, 99, dtype=np.int64)
tc = np.full(r.size, 99, dtype=np.int64)
for i in range(1, K + 1):
    odd = (cur & 1) == 1
    a += odd
    cur = np.where(odd, (3 * cur - 1) >> 1, cur >> 1)
    ta = np.where((cur < r) & (ta == 99), i, ta)
    tc = np.where(((3.0 ** a) < float(2 ** i)) & (tc == 99), i, tc)
n_act = int((ta < 99).sum())
n_coef = int((tc < 99).sum())
strict = int(((ta < 99) & (ta < tc)).sum())
late = int(((tc < 99) & (ta > tc)).sum())
print(f"K=20 minus residues in [1,2^20): descend within 20: {n_act}; coefficient time <= 20: {n_coef}; "
      f"strict early (ta<tc): {strict}; actual later than coefficient (ta>tc): {late}; 2^20-1 = {(1 << K) - 1}")
check(n_act == 1021247 and n_coef == 1021247 and strict == 0 and late == 0, "residue check")
check(n_coef == exp_desc[-1] - 1, "coefficient-descending residues = word count minus the r=0 word")
del r, cur, a, ta, tc

# ----------------------------------------------------------------------------
print()
print("=== A4: greedy G_- to 10^6, exact one-step law, BFS rescues ===")
KM = {1: 1, 2: 0, 4: 3, 5: 0, 7: 1, 8: 2}
for rr, k in KM.items():
    check((rr * 2 ** k) % 9 in (2, 5) and all((rr * 2 ** j) % 9 not in (2, 5) for j in range(k)), "k table")
    check(((rr * 2 ** k) + 1) % 3 == 0 and (((rr * 2 ** k) + 1) // 3) % 3 != 0, "G_- lands off 3Z")


def Gm(m):
    k = KM[m % 9]
    return (m * 2 ** k + 1) // 3


below = one = cyc = cap = 0
hist = {}
peak = (0, 0)
onestep_res = set()
for m in range(2, 10 ** 6 + 1):
    if m % 3 == 0:
        continue
    x = m
    pk = m
    s = 0
    while True:
        x = Gm(x)
        s += 1
        pk = max(pk, x)
        if x < m:
            below += 1
            break
        if x == 1:
            one += 1
            break
        if x == m:
            cyc += 1
            break
        if s > 2000:
            cap += 1
            break
    hist[s] = hist.get(s, 0) + 1
    if s == 1 and x < m:
        onestep_res.add(m % 9)
    if pk > peak[0]:
        peak = (pk, m)
tot = below + one + cyc + cap
print(f"starts {tot}: below {below}, one {one}, cycle {cyc}, cap {cap}; histogram head {[(k, hist[k]) for k in sorted(hist)[:4]]}; max steps {max(hist)}")
print(f"peak {peak[0]} at m={peak[1]}; 3*2^24-1 = {3 * 2 ** 24 - 1}; (3^13-1)/2 = {(3 ** 13 - 1) // 2}; 3*2^24+1 = {3 * 2 ** 24 + 1}")
check(tot == 666666 and below == 666665 and one == 0 and cyc == 1 and cap == 0, "greedy statuses")
check(hist[1] == 444444 and hist[2] == 148148 and hist[3] == 24692 and hist[4] == 24691 and max(hist) == 35 and hist[35] == 1, "histogram")
check(peak == (3 * 2 ** 24 - 1, (3 ** 13 - 1) // 2) and peak[1] == 797161, "peak")
# exact one-step law: G_-(m) < m iff k_-(m mod 9) <= 1 (m >= 2): k=0 gives (m+1)/3 < m; k=1 gives (2m+1)/3 < m iff m > 1;
# k=2,3 give (4m+1)/3, (8m+1)/3 > m.
one_res = {rr for rr, k in KM.items() if k <= 1}
cnt_res = sum(1 for m in range(2, 10 ** 6 + 1) if m % 9 in one_res)
print(f"one-step residues mod 9 observed {sorted(onestep_res)} = predicted {sorted(one_res)}; #m in [2,10^6] with m mod 9 in that set = {cnt_res}; 444444/666666 = {Fraction(444444, 666666)}")
check(onestep_res == one_res and cnt_res == 444444 and Fraction(444444, 666666) == Fraction(2, 3), "exact one-step law")
check(one == 0, "status 'one' is empty by construction: 1 < m triggers 'below' first")


def bfs_inv(m, capv):
    par = {m: None}
    dq = deque([m])
    while dq:
        x = dq.popleft()
        ys = [2 * x]
        if (x + 1) % 3 == 0 and ((x + 1) // 3) % 3:
            ys.append((x + 1) // 3)
        for y in ys:
            if y in par or y > capv:
                continue
            par[y] = x
            if y < m:
                path = [y]
                while par[path[-1]] is not None:
                    path.append(par[path[-1]])
                return path[::-1]
            dq.append(y)


p4 = bfs_inv(4, 10 ** 6)
print(f"BFS inverse rescue of 4: {len(p4) - 1} arrows: {p4}")
check(len(p4) - 1 == 14 and p4[-1] == 2, "m=4 rescue length")
# the SCC note's compound path 4 -> 11 -> 59 -> 20 -> 7 -> 5 -> 2 with k-word (3,4,0,0,1,0)
comp = [4, 11, 59, 20, 7, 5, 2]
kw = [3, 4, 0, 0, 1, 0]
for (u, w), k in zip(zip(comp, comp[1:]), kw):
    check((u * 2 ** k + 1) // 3 == w and (u * 2 ** k + 1) % 3 == 0, "compound path step")
print(f"SCC-note compound rescue {comp} k-word {kw}: total arrows = sum(k+1) = {sum(k + 1 for k in kw)}")
check(sum(k + 1 for k in kw) == 14, "compound rescue arrow count")
# forward BFS distances from 1 in E_- (n -> 3n-1 all n, n -> n/2 even)
targets = [5, 7, 17, 25, 37, 55, 41, 61, 91]
dist = {1: 0}
dq = deque([1])
while dq and not all(t in dist for t in targets):
    x = dq.popleft()
    ys = [3 * x - 1] + ([x // 2] if x % 2 == 0 else [])
    for y in ys:
        if y in dist or y > 5 * 10 ** 6:
            continue
        dist[y] = dist[x] + 1
        dq.append(y)
got = {t: dist[t] for t in targets}
print(f"forward BFS distances from 1: {got}")
check(got == {5: 2, 7: 4, 17: 13, 25: 15, 37: 17, 55: 19, 41: 4, 61: 6, 91: 8}, "forward distances")

# ----------------------------------------------------------------------------
print()
print("=== A5: Berggren readings, corrected coupling breakdown, k=1<->k=2 bijection ===")


def is_edge(s, t, b):
    val = 3 * s + b
    return val > 0 and val % t == 0 and (val // t) >= 2 and (val // t) & (val // t - 1) == 0


def readings(p, q):
    out = set()
    for b in (1, -1):
        for s, t in ((p, q), (q, p)):
            if is_edge(s, t, b):
                out.add((s, t, b, v2(3 * s + b)))
    return out


def odd_cycle(n0, sign):
    c = [n0]
    x = n0
    while True:
        y = 3 * x + sign
        y >>= v2(y)
        if y == n0:
            return c
        c.append(y)
        x = y


edges = []
for sign, mins in ((1, [1]), (-1, [1, 5, 17])):
    for n0 in mins:
        c = odd_cycle(n0, sign)
        for i, x in enumerate(c):
            y = c[(i + 1) % len(c)]
            edges.append((sign, x, y, v2(3 * x + sign)))
check(len(edges) == 11, "11 cycle edges")
child_other = []
parent_other = []
for sign, x, y, k in edges:
    s, t = max(x, y), min(x, y)
    kids = [(s + 2 * t, t), (2 * s + t, s), (2 * s - t, s)]
    if any(any(b == -sign for _, _, b, _ in readings(p, q)) for p, q in kids):
        child_other.append((sign, x, y))
    if s == t:
        par = (s, t)
    elif s > 3 * t:
        par = (s - 2 * t, t)
    elif s > 2 * t:
        par = (t, s - 2 * t)
    else:
        par = (t, 2 * t - s)
    if any(b == -sign for _, _, b, _ in readings(*par)):
        parent_other.append((sign, x, y, par))
print(f"edges with an other-sheet child ({len(child_other)}): {child_other}")
print(f"edges that are B-children of an other-sheet edge ({len(parent_other)}): {parent_other}")
check(len(child_other) == 9 and len(parent_other) == 5, "coupling counts 9 and 5")
k1_minus = [(x, y) for sign, x, y, k in edges if sign == -1 and k == 1]
check(len(k1_minus) == 7 and (1, 1) in k1_minus, "seven k=1 minus edges incl. 1->1")
check(set(child_other) == {(-1, x, y) for x, y in k1_minus} | {(1, 1, 1), (-1, 7, 5)},
      "corrected breakdown: 7 k=1 minus edges (incl. 1->1) + plus 1->1 + 7->5 (shares the pair of 5->7)")
check((-1, 91, 17) not in child_other and all(e[1:3] != (91, 17) for e in parent_other), "91->17 uncoupled")
print("corrected breakdown of the 9: the seven k=1 minus edges (1->1 included) + the plus edge 1->1 + the k=2 edge 7->5 (pair {5,7})")
# bijection with the y = x edge case and the mod-8 law
XB = 10 ** 5
for b in (1, -1):
    k1 = [(x, (3 * x + b) // 2) for x in range(1, XB + 1, 2) if v2(3 * x + b) == 1]
    check(all(y >= x for x, y in k1) and [x for x, y in k1 if y == x] == ([1] if b == -1 else []),
          "k=1 edges have y >= x, equality only at the minus edge 1->1")
    check(all(x % 4 == (3 if b == 1 else 1) for x, _ in k1), "k=1 sources mod 4")
    img = {(2 * x + b, y) for x, y in k1}
    umax = 2 * XB + 1
    k2 = {(u, (3 * u - b) // 4) for u in range(1, umax + 1, 2) if v2(3 * u - b) == 2}
    check(all(u % 8 == (7 if b == 1 else 1) for u, _ in k2), "k=2 sources of sheet -b are u = 7 mod 8 (minus) / 1 mod 8 (plus)")
    check(img <= k2 and {e for e in k2 if e[0] <= 2 * XB - 1} <= img, "bijection onto k=2 edges of the other sheet")
    check(all(is_edge(u, y, -b) for u, y in img), "images are (3,-b,2) edges")
    print(f"sheet b={b:+d}: {len(k1)} k=1 edges x<=10^5 (x = {3 if b == 1 else 1} mod 4) -> {len(img)} k=2 edges of sheet {-b:+d} "
          f"(sources 2x{b:+d} = {7 if b == 1 else 1} mod 8); 3-adic aside: the mod-4 law u = {3 if b == 1 else 1} mod 4 is only necessary")
    check(len(k1) == 25000 and len(img) == 25000, "25000")

# ----------------------------------------------------------------------------
print()
print("=== A6: budget sums and gates ===")


def budget(n0, sign):
    x, K, i, seen, qs = n0, 0, 0, {}, []
    while x not in seen:
        seen[x] = i
        qs.append(Fraction(2 ** K, 3 ** i))
        k = v2(3 * x + sign)
        x = (3 * x + sign) >> k
        K += k
        i += 1
    i0 = seen[x]
    ratio = Fraction(2 ** K, 3 ** i) / qs[i0]
    if ratio >= 1:
        return None, ratio
    return sum(qs[:i0], Fraction(0)) + sum(qs[i0:], Fraction(0)) / (1 - ratio), ratio


starts = (1, 3, 5, 7, 9, 11, 17, 25, 37, 41, 55, 61, 91, 27, 1001)
sums = []
for n0 in starts:
    val, ratio = budget(n0, -1)
    check(val == 3 * n0 and ratio in (Fraction(2, 3), Fraction(8, 9), Fraction(2048, 2187)), f"budget {n0}")
    sums.append(int(val))
print(f"minus budgets 3 n0 for {starts}: {sums}")
check(budget(1, +1)[1] == Fraction(4, 3), "plus ratio 4/3")
gates = []
for n0 in (1, 5, 17):
    c = odd_cycle(n0, -1)
    L = len(c)
    ks = [v2(3 * x - 1) for x in c]
    K = sum(ks)
    B = sum(3 ** (L - 1 - i) * 2 ** sum(ks[:i]) for i in range(L))
    check(n0 * (3 ** L - 2 ** K) == B, "gate")
    gates.append((L, K, B, 3 ** L - 2 ** K))
print(f"gates (L, K, B, 3^L-2^K) for 1, 5, 17: {gates}")
check(gates == [(1, 1, 1, 1), (2, 3, 5, 1), (7, 11, 2363, 139)], "gate values")
print()
print("AUDIT ALL CHECKS PASSED")
