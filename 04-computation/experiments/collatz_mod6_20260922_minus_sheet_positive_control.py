#!/usr/bin/env python3
"""collatz_mod6_20260922_minus_sheet_positive_control.py

Lane minus_sheet_positive_control (session collatz-mod6-20260922).

The 3n-1 sheet T_-(n) = (3n-1)/2^v on the positive integers has three known
cycles.  Every "descent" heuristic for 3n+1 is run here against that sheet as a
positive control: a statement that also holds on the 3n-1 sheet cannot, by
itself, imply convergence to 1.

Sections (each prints its own header; S0 is the session-lead entropy probe):
  S0  binary entropy of 8/pi^2 (session lead probe)
  S1  basin census of T_- on [1, 10^7] (all n and odd n), plus-sheet control
  S2  Terras-type prefix-descent densities mod 2^J, J <= 20, both sheets
      (sign-blind parity-vector bijection), actual vs coefficient stopping time
  S3  E_- reachability from 1: greedy G_- to 10^6, BFS fallback, explicit paths
  S4  Berggren B3 sheet coupling of the cycle edges (children and parents)
  S5  additive budget sum q_i = 3 n_0 on the minus-sheet cycles (exact)
All checks use explicit `raise`; timing goes to stderr only.
"""
import math
import sys
import time
from collections import deque
from fractions import Fraction

import numpy as np

T0 = time.time()


def tlog(msg):
    print(f"[timing] {msg} t={time.time() - T0:.1f}s", file=sys.stderr)


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


# ----------------------------------------------------------------------------
print("=== S0: binary entropy of 8/pi^2 (session lead probe) ===")
p = 8 / math.pi ** 2
H = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
print(f"8/pi^2 = {p:.6f}   H2(8/pi^2) = {H:.5f} bits   (paste says 0.704; 'zero entropy' is false: H>0)")
check(abs(H - 0.70028) < 1e-5, "entropy value")

# ----------------------------------------------------------------------------
print()
print("=== S1: basin census of T_- (n -> n/2 even, n -> 3n-1 odd) on [1, 10^7] ===")
XMAX = 10 ** 7
STEP_CAP = 5000
SAFE = 1 << 61


def basin_census(X, sign, roots):
    n = np.arange(1, X + 1, dtype=np.int64)
    target = np.zeros(X + 1, dtype=np.int64)
    root_mask = np.isin(n, np.array(sorted(roots), dtype=np.int64))
    active = n[~root_mask]
    cur = active.copy()
    steps = 0
    maxval = 0
    unresolved = None
    while active.size:
        odd = (cur & 1).astype(bool)
        nxt = cur >> 1
        np.multiply(cur, 3, out=cur)
        cur += sign
        np.copyto(cur, nxt, where=~odd)
        del nxt
        m = int(cur.max())
        if m > maxval:
            maxval = m
        check(m < SAFE, "int64 overflow guard")
        steps += 1
        done = cur < active
        target[active[done]] = cur[done]
        keep = ~done
        active = active[keep]
        cur = cur[keep]
        if steps >= STEP_CAP:
            unresolved = active.copy()
            break
    for r in roots:
        target[r] = r
    never = n[target[1:] == n].tolist()
    del active, cur
    # pointer jumping to the cycle minimum (in place; target[n] < n or a root)
    root = target
    del target
    root[0] = 0
    rounds = 0
    while True:
        new = root[root]
        rounds += 1
        if np.array_equal(new, root):
            del new
            break
        root = new
        del new
        check(rounds < 200, "pointer jumping did not converge")
    return root, never, steps, maxval, unresolved


roots_minus = [1, 5, 17]
root_m, never, steps_m, maxval_m, unres_m = basin_census(XMAX, -1, roots_minus)
tlog("minus census")
print(f"minus sheet: max stopping time (all n <= 10^7, first value < n) = {steps_m}; "
      f"max orbit value seen = {maxval_m}; unresolved after {STEP_CAP} steps = "
      f"{0 if unres_m is None else unres_m.size}")
check(unres_m is None, "minus sheet: some start did not descend within the step cap (escape/new cycle candidate)")
check(bool(np.all(np.isin(root_m[1:], roots_minus))), "every root is a known cycle minimum")
print(f"minus sheet: n <= 10^7 whose orbit never drops below n: {never}")
check(never == roots_minus, "exactly the three cycle minima never descend")
print()
print("X          | basin{1}  all n   odd n | basin{5,7} all n  odd n | basin{17..91} all n  odd n")
for X in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    sub = root_m[1:X + 1]
    odds = np.arange(1, X + 1, 2, dtype=np.int64)
    subo = root_m[odds]
    row = []
    for r in roots_minus:
        ca = int(np.count_nonzero(sub == r))
        co = int(np.count_nonzero(subo == r))
        row.append((ca, co))
    check(sum(c for c, _ in row) == X, "basin counts partition [1,X]")
    check(sum(c for _, c in row) == (X + 1) // 2, "odd basin counts partition odd [1,X]")
    print(f"{X:<10d} | " + " | ".join(f"{ca:>9d} {ca / X:.5f} {co:>8d} {co / ((X + 1) // 2):.5f}" for ca, co in row))
# plus-sheet control
del root_m
root_p, never_p, steps_p, maxval_p, unres_p = basin_census(XMAX, +1, [1])
tlog("plus census")
print(f"plus sheet control: every n <= 10^7 in basin{{1}}: {bool(np.all(root_p[1:] == 1))}; max stopping time = {steps_p}; "
      f"max orbit value seen = {maxval_p}; unresolved = {0 if unres_p is None else unres_p.size}")
check(unres_p is None and bool(np.all(root_p[1:] == 1)), "plus sheet control")
print(f"plus sheet: n <= 10^7 whose orbit never drops below n: {never_p}")
check(never_p == [1], "plus: only 1 never descends")
del root_p

# ----------------------------------------------------------------------------
print()
print("=== S2: Terras-type prefix-descent density mod 2^J, both sheets (J <= 20) ===")
JMAX = 20
thr = {}
for i in range(1, JMAX + 1):
    a = 0
    while 3 ** (a + 1) < 2 ** i:
        a += 1
    thr[i] = a  # largest a with 3^a < 2^i


def prefix_descent_counts(J, sign):
    r = np.arange(1 << J, dtype=np.int64)
    cur = r.copy()
    a = np.zeros(1 << J, dtype=np.int64)
    flag = np.zeros(1 << J, dtype=bool)
    words = np.zeros(1 << J, dtype=np.int64)
    for i in range(1, J + 1):
        odd = (cur & 1).astype(bool)
        words |= odd.astype(np.int64) << (i - 1)
        cur = np.where(odd, (3 * cur + sign) >> 1, cur >> 1)
        a += odd
        flag |= a <= thr[i]
    return int(flag.sum()), int(np.unique(words).size)


def word_count_no_descent(J):
    # number of parity words of length J with 3^{a_i} >= 2^i for every i <= J
    # dynamic programming over (i, a)
    dp = {0: 1}
    for i in range(1, J + 1):
        nd = {}
        for a, c in dp.items():
            for step in (0, 1):
                a2 = a + step
                if 3 ** a2 >= 2 ** i:
                    nd[a2] = nd.get(a2, 0) + c
        dp = nd
    return sum(dp.values())


print("J  | #descend(+) #descend(-) | words(+) words(-) = 2^J | no-descent words (DP) | fraction")
for J in range(1, JMAX + 1):
    cp, wp = prefix_descent_counts(J, +1)
    cm, wm = prefix_descent_counts(J, -1)
    nd = word_count_no_descent(J)
    check(cp == cm, f"prefix descent counts differ at J={J}")
    check(wp == (1 << J) and wm == (1 << J), f"parity map not a bijection at J={J}")
    check(cp + nd == (1 << J), f"DP word count disagrees at J={J}")
    print(f"{J:<2d} | {cp:>11d} {cm:>11d} | {wp:>8d} {wm:>8d} {1 << J:>8d} | {nd:>10d} | {cp / (1 << J):.6f}")
tlog("terras")

# actual vs coefficient stopping time on integers (T map on all n), both sheets
print()
print("actual sigma(n) = min k: T^k(n) < n; coefficient sigma_c(n) = min k: 3^{a_k} < 2^k  (n <= 10^6)")
NS = 10 ** 6
KCAP = 400
for sign, name in ((+1, "plus"), (-1, "minus")):
    strict = []
    inf_act = []
    for n in range(1, NS + 1):
        cur = n
        a = 0
        sa = None
        sc = None
        for k in range(1, KCAP + 1):
            if cur & 1:
                cur = (3 * cur + sign) >> 1
                a += 1
            else:
                cur >>= 1
            if sa is None and cur < n:
                sa = k
            if sc is None and 3 ** a < 2 ** k:
                sc = k
            if sa is not None and sc is not None:
                break
        if sa is None:
            inf_act.append(n)
        elif sc is None or sa < sc:
            strict.append(n)
        if sc is not None and sa is not None and sa > sc:
            raise RuntimeError(f"actual > coefficient stopping time at n={n} sheet {name}")
    print(f"{name} sheet: n <= 10^6 with sigma actual infinite (within {KCAP} steps): {inf_act}; "
          f"count with sigma_actual < sigma_coefficient: {len(strict)}; smallest: {strict[:12]}")
    if sign == +1:
        check(strict == [] and inf_act == [1], "plus sheet: actual = coefficient stopping time")
    else:
        check(inf_act == [1, 5, 17], "minus sheet: only the cycle minima never descend")
        check(strict == [], "minus sheet: no strict early descent below 10^6")
tlog("stopping times")

# Reduction lemma (proved in the note): a strict early descent at first-descent time k for some n
# forces one for the residue r = n mod 2^k in [1, 2^k) at the same time k.  Exact check for k <= 20:
print()
print("residue-level check of the reduction lemma: r in [1, 2^K), K = 20, minus sheet")
K = 20
r = np.arange(1, 1 << K, dtype=np.int64)
cur = r.copy()
a = np.zeros(r.size, dtype=np.int64)
act = np.zeros(r.size, dtype=np.int64)   # first descent time within K steps (0 = none)
coef = np.zeros(r.size, dtype=np.int64)  # first coefficient time within K steps (0 = none)
for i in range(1, K + 1):
    odd = (cur & 1).astype(bool)
    cur = np.where(odd, (3 * cur - 1) >> 1, cur >> 1)
    a += odd
    d = (cur < r) & (act == 0)
    act[d] = i
    c = (a <= thr[i]) & (coef == 0)
    coef[c] = i
bad = r[(act > 0) & ((coef == 0) | (act < coef))]
n_act = int(np.count_nonzero(act > 0))
n_eq = int(np.count_nonzero((act > 0) & (act == coef)))
print(f"residues with first descent within {K} steps: {n_act}; of these with actual = coefficient time: {n_eq}; strict early descents: {bad.size}")
check(bad.size == 0, "reduction-lemma check: no residue has a strict early descent within 20 steps")
check(bool(np.all((coef == 0) | (act == coef))), "coefficient descent implies actual descent (B<0)")
del r, cur, a, act, coef

# ----------------------------------------------------------------------------
print()
print("=== S3: E_- reachability from 1 (arrows n->n/2 even, n->3n-1 all n): greedy G_- and BFS fallback ===")
KMIN = {1: 1, 2: 0, 4: 3, 5: 0, 7: 1, 8: 2}  # k_-(r): minimal k with 2^k r in {2,5} mod 9
for r, k in KMIN.items():
    check((r << k) % 9 in (2, 5), "k_- table")
    check(all((r << j) % 9 not in (2, 5) for j in range(k)), "k_- minimality")


def G_minus(m):
    k = KMIN[m % 9]
    return ((m << k) + 1) // 3


MQ = 10 ** 6
GCAP = 2000
status = {"below": 0, "one": 0, "cycle": 0, "cap": 0}
hard = []
steps_hist = {}
peak_max = (0, 0)
for m in range(2, MQ + 1):
    if m % 3 == 0:
        continue
    cur = m
    seen_cycle = False
    s = 0
    pk = m
    while True:
        cur = G_minus(cur)
        s += 1
        if cur > pk:
            pk = cur
        if cur < m:
            status["below"] += 1
            break
        if cur == 1:
            status["one"] += 1
            break
        if cur == m or s > GCAP:
            status["cycle" if cur == m else "cap"] += 1
            hard.append((m, "cycle" if cur == m else "cap"))
            break
    steps_hist[s] = steps_hist.get(s, 0) + 1
    if pk > peak_max[0]:
        peak_max = (pk, m)
tot = sum(status.values())
print(f"greedy G_-(m) = (2^k m + 1)/3 on 2 <= m <= 10^6, 3 !| m: {tot} starts; statuses {status} (step cap {GCAP})")
print(f"greedily certified (below or one): {status['below'] + status['one']} / {tot} = "
      f"{(status['below'] + status['one']) / tot:.6f}; hard cases: {hard}")
print(f"greedy steps to certify, histogram (steps: count): " + ", ".join(f"{k}: {steps_hist[k]}" for k in sorted(steps_hist)))
print(f"largest greedy peak before certification: {peak_max[0]} at m = {peak_max[1]}")
check(hard == [(4, "cycle")], "the only greedy-hard start is m = 4 (2-cycle {4,11})")
check(G_minus(4) == 11 and G_minus(11) == 4, "G_- 2-cycle {4,11}")


def bfs_inverse_below(m, cap):
    """BFS over inverse moves x->2x and x->(x+1)/3 (integral, not 0 mod 3) to a value < m."""
    parent = {m: None}
    dq = deque([m])
    while dq:
        x = dq.popleft()
        nxt = [2 * x]
        if (x + 1) % 3 == 0 and ((x + 1) // 3) % 3 != 0:
            nxt.append((x + 1) // 3)
        for y in nxt:
            if y in parent or y > cap:
                continue
            parent[y] = x
            if y < m:
                path = [y]
                while parent[path[-1]] is not None:
                    path.append(parent[path[-1]])
                return path[::-1]
            dq.append(y)
    return None


path4 = bfs_inverse_below(4, 10 ** 6)
print(f"BFS fallback for m = 4 (inverse moves): {' -> '.join(map(str, path4))}")
check(path4 is not None and path4[0] == 4 and path4[-1] < 4, "BFS rescue of m = 4")

# explicit forward E_- paths from 1 to the odd members of the three cycles
targets = [5, 7, 17, 25, 37, 55, 41, 61, 91]
CAPV = 5 * 10 ** 6
parent = {1: None}
dq = deque([1])
found = {}
while dq and len(found) < len(targets):
    x = dq.popleft()
    nxt = [3 * x - 1]
    if x % 2 == 0:
        nxt.append(x // 2)
    for y in nxt:
        if y in parent or y > CAPV:
            continue
        parent[y] = x
        if y in targets:
            found[y] = True
        dq.append(y)
check(len(found) == len(targets), "BFS from 1 reaches every odd cycle member")
print("BFS-shortest forward E_- paths from 1 (values capped at 5*10^6):")
for t in targets:
    path = [t]
    while parent[path[-1]] is not None:
        path.append(parent[path[-1]])
    path = path[::-1]
    print(f"  1 -> {t}: {len(path) - 1} arrows: {' -> '.join(map(str, path))}")
check(parent[5] == 2 and parent[2] == 1, "1 -> 2 -> 5")
tlog("E_- reachability")

# ----------------------------------------------------------------------------
print()
print("=== S4: Berggren B3 sheet coupling of the cycle edges ===")


def v2(x):
    return (x & -x).bit_length() - 1


def odd_cycle(n0, sign):
    cyc = [n0]
    x = n0
    while True:
        y = 3 * x + sign
        y >>= v2(y)
        if y == n0:
            break
        cyc.append(y)
        x = y
    return cyc


def edges_of(cyc, sign):
    out = []
    for i, x in enumerate(cyc):
        y = cyc[(i + 1) % len(cyc)]
        k = v2(3 * x + sign)
        check((3 * x + sign) == (y << k), "edge identity")
        out.append((x, y, k))
    return out


def readings(p, q):
    """All (3,b',j) readings of the unordered pair {p,q}, p>q, as an edge in either orientation."""
    out = []
    for b in (+1, -1):
        for (s, t) in ((p, q), (q, p)):
            val = 3 * s + b
            if val > 0 and val % t == 0 and (val // t) & (val // t - 1) == 0 and val // t >= 2:
                rd = (s, t, b, (val // t).bit_length() - 1)
                if rd not in out:
                    out.append(rd)
    return out


def children(x, y):
    if y > x:
        s, t = y, x
    else:
        s, t = x, y
    return {"B1": (s + 2 * t, t), "B2": (2 * s + t, s), "B3": (2 * s - t, s)}


def berggren_parent(s, t):
    if s == t:
        return ("B3-fixed", (s, t))
    r = s / t
    if r > 3:
        return ("B1", (s - 2 * t, t))
    if r > 2:
        return ("B2", (t, s - 2 * t))
    return ("B3", (t, 2 * t - s))


sheet_name = {+1: "plus", -1: "minus"}
all_edges = []
for sign, mins in ((+1, [1]), (-1, [1, 5, 17])):
    for n0 in mins:
        cyc = odd_cycle(n0, sign)
        for e in edges_of(cyc, sign):
            all_edges.append((sign, n0, e))
print("edge (sheet, cycle-min): x->y k | children B1,B2,B3 with (3,+-1) readings | Berggren parent pair and its readings")
count_child_other = 0
count_parent_other = 0
kcount = {}
for sign, n0, (x, y, k) in all_edges:
    ch = children(x, y)
    ch_str = []
    other_child = False
    for name, (p, q) in ch.items():
        rd = readings(p, q)
        rd_s = ",".join(f"{s}->{t}({sheet_name[b]},k={j})" for s, t, b, j in rd) or "none"
        ch_str.append(f"{name}=({p},{q}):{rd_s}")
        if any(b == -sign for _, _, b, _ in rd):
            other_child = True
    s, t = max(x, y), min(x, y)
    pname, (ps, pt) = berggren_parent(s, t)
    prd = readings(ps, pt) if ps != pt else readings(1, 1)
    if ps == pt:
        prd = [(1, 1, +1, 2), (1, 1, -1, 1)]
    prd_s = ",".join(f"{a}->{c}({sheet_name[b]},k={j})" for a, c, b, j in prd) or "none"
    other_parent = any(b == -sign for _, _, b, _ in prd)
    count_child_other += other_child
    count_parent_other += other_parent
    kcount[(sign, k)] = kcount.get((sign, k), 0) + 1
    print(f"({sheet_name[sign]},{n0}): {x}->{y} k={k} | " + " ".join(ch_str) +
          f" | parent via {pname}: ({ps},{pt}):{prd_s}")
print(f"minus-sheet cycle edges: {sum(v for (s, _), v in kcount.items() if s == -1)} "
      f"(k=1: {kcount.get((-1, 1), 0)}, k=2: {kcount.get((-1, 2), 0)}, k=4: {kcount.get((-1, 4), 0)}); "
      f"plus-sheet cycle edges: {kcount.get((+1, 2), 0)} (k=2)")
print(f"edges with a child on the other sheet: {count_child_other}; edges that are a Berggren child of an other-sheet edge: {count_parent_other} (of {len(all_edges)})")
check(kcount == {(+1, 2): 1, (-1, 1): 7, (-1, 2): 2, (-1, 4): 1}, "k census of cycle edges")
pairs = sorted(set((max(x, y), min(x, y)) for sign, _, (x, y, _) in all_edges if sign == -1))
print(f"distinct unordered pairs among the minus-sheet cycle edges: {len(pairs)} ({{5,7}} carries both 5->7 k=1 and 7->5 k=2)")
check(count_child_other == 9 and count_parent_other == 5 and len(pairs) == 9, "coupling counts")

# the proved bijection: B3 on k=1 edges of one sheet = k=2 edges of the other sheet, x <= 10^5
XB = 10 ** 5
for sign in (+1, -1):
    k1 = [(x, (3 * x + sign) >> 1) for x in range(1, XB + 1, 2) if v2(3 * x + sign) == 1]
    k2_other = set((u, (3 * u - sign) >> 2) for u in range(1, 2 * XB + 2, 2) if v2(3 * u - sign) == 2)
    img = set((2 * y - x, y) for x, y in k1)  # B3 child in case A (y > x): (2y - x, y) = (2x + sign, y)
    check(all(2 * y - x == 2 * x + sign for x, y in k1), "B3 child source is 2x + b")
    check(img <= k2_other, "B3 image of k=1 edges lies in the k=2 edges of the other sheet")
    missing = [e for e in k2_other if e not in img and e[0] <= 2 * XB - 1]
    check(not missing, "every k=2 edge of the other sheet is a B3 child")
    print(f"{sheet_name[sign]} sheet: {len(k1)} k=1 edges with x <= 10^5 (x = {3 if sign == 1 else 1} mod 4) map by B3 onto "
          f"{len(img)} k=2 edges of the {sheet_name[-sign]} sheet (source 2x{'+' if sign == 1 else '-'}1, i.e. u = {3 if sign == 1 else 1} mod 4); bijection verified")
tlog("berggren")

# ----------------------------------------------------------------------------
print()
print("=== S5: additive budget sum_{i>=0} q_i = 3 n_0 on the minus sheet (exact rational) ===")


def budget_sum(n0, sign):
    """sum_{i>=0} 2^{K_i}/3^i along the odd-skeleton orbit, exact once the orbit is periodic."""
    orbit = []
    x = n0
    K = 0
    i = 0
    seen = {}
    qs = []
    while x not in seen:
        seen[x] = i
        qs.append(Fraction(2 ** K, 3 ** i))
        k = v2(3 * x + sign)
        x = (3 * x + sign) >> k
        K += k
        i += 1
    i0 = seen[x]
    pre = sum(qs[:i0], Fraction(0))
    per = sum(qs[i0:], Fraction(0))
    ratio = Fraction(2 ** K, 3 ** i) / qs[i0]  # 2^{K_period}/3^{L_period}
    if ratio >= 1:
        return None, ratio
    return pre + per / (1 - ratio), ratio


for n0 in (1, 3, 5, 7, 9, 11, 17, 25, 37, 41, 55, 61, 91, 27, 1001):
    val, ratio = budget_sum(n0, -1)
    check(val == 3 * n0, f"budget identity at n0={n0}")
    print(f"minus sheet n0={n0:<5d}: sum q_i = {val} = 3 n0; period ratio 2^K/3^L = {ratio} < 1")
val, ratio = budget_sum(1, +1)
print(f"plus sheet n0=1: period ratio 2^K/3^L = {ratio} > 1, sum q_i diverges (D_L -> +infinity on the cycle)")
check(val is None and ratio == Fraction(4, 3), "plus-sheet cycle ratio")
# the three minus cycles: 2^K < 3^L
for n0 in (1, 5, 17):
    cyc = odd_cycle(n0, -1)
    L = len(cyc)
    K = sum(v2(3 * x - 1) for x in cyc)
    B = sum(3 ** (L - 1 - i) * 2 ** sum(v2(3 * c - 1) for c in cyc[:i]) for i in range(L))
    check(n0 * (3 ** L - 2 ** K) == B, "cycle gate n0 = B/(3^L - 2^K)")
    print(f"minus cycle min {n0}: L={L} odd steps, K={K} halvings, 2^K={2 ** K} < 3^L={3 ** L}, gate B={B}, n0 = B/(3^L-2^K) = {B}/{3 ** L - 2 ** K}")
tlog("budget")
import hashlib
with open(__file__, "rb") as fh:
    print(f"source sha256 {hashlib.sha256(fh.read()).hexdigest()}")
print()
print("ALL CHECKS PASSED")
