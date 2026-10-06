#!/usr/bin/env python3
"""Task B proof checks (chessboard weave, LRC reading), 2026-10-06.  Exact arithmetic.

P1  n = 2: closed form of tau(a,b) and the theorem  max tau = 4/9, only at (1,3)  (the camel);
    tau > 1/3 iff v = (1,3k), tau(1,3k) = 1/3 + 1/(9k); max over vmin >= 2 is 4/15 at (2,5).
P2  n = 3, speed 1 present: max tau(1,b,w) = 7/16, only at (1,3,12)  (hand proof; the
    printed safe sets and the finite residue (1,6,w), 7 <= w <= 11, are the inputs).
P3  the family F_N = {1,3,4,...,N} (n = N-1 speeds, threshold 1/N):
    Safe(F_N) cap (0,1/2] = [tau_N, 1/2 - 1/(4N)]  (N >= 4),  tau_N = 1/2 - (N-2)/(2 N k_o),
    k_o = largest odd element of F_N  (hand proof; exact check N = 3..40).
P4  kill-and-replace lemma: T tight with Safe(T) in (1/N)Z, s in T, m = kN not in T  =>
    Safe((T - {s}) + {m}) lies in the danger set {t : ||s t|| < 1/N} of the removed speed.
    (proof one line; exact check on all tight sets of the A-script, m <= 60).
P5  the other families met in the census: tau(1,3,12k) = 5/12 + 1/(48k) (n=3),
    tau(1,3,4,5k) = 2/5 + 1/(25k) (n=4)  (exact checks k <= 200).
Reproduce: python3 chessboard_weave_20261006_lrc_B_proofs.py > chessboard_weave_20261006_lrc_B_proofs.out
"""
import sys
import time
from fractions import Fraction
from math import gcd

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from chessboard_weave_20261006_lrc_core import safe_components, tau, check


def safe_half(v, N=None):
    comps, D = safe_components(v, N)
    return [(Fraction(a, D), Fraction(b, D)) for a, b in comps if 2 * a <= D]


def fmt(iv):
    return ", ".join(f"{{{a}}}" if a == b else f"[{a},{b}]" for a, b in iv)


def dist(x):
    x = x - (x.numerator // x.denominator)
    return min(x, 1 - x)


T0 = time.time()
# ------------------------------------------------------------------ P1
print("=" * 78)
print("P1. n = 2 (threshold 1/3).  THEOREM: for coprime 1 <= a < b,")
print("      tau(a,b) = min{ t >= 1/(3a) : ||t b|| >= 1/3 }  and  1/(3a) <= tau(a,b) <= 2/(3a);")
print("    tau(1,b) = 1/3 if 3 does not divide b, tau(1,3k) = 1/3 + 1/(9k).  Hence max tau = 4/9,")
print("    attained ONLY at (1,3), and tau > 1/3 iff v = (1,3k).")
print("""    PROOF.  Safe(a) = U_m [(m+1/3)/a, (m+2/3)/a]; its first interval is
    I0 = [1/(3a), 2/(3a)], so tau >= 1/(3a).  The danger set of b is a union of OPEN
    intervals of length 2/(3b) separated by safe closed intervals.  If b <= 2a then
    1/(3a) <= 2/(3b), i.e. 1/(3a) lies in b's first safe interval [1/(3b), 2/(3b)],
    so tau = 1/(3a).  If b > 2a then |I0| = 1/(3a) > 2/(3b), so the closed interval I0
    cannot sit inside one open danger interval of b: I0 contains a safe point of b and
    tau is the first one, <= 2/(3a).  For a >= 2 this gives tau <= 1/3.  For a = 1,
    I0 = [1/3, 2/3]; t = 1/3 is safe for b iff 3 does not divide b (||b/3|| = 1/3);
    if b = 3k the next safe point of b is (k + 1/3)/(3k) = 1/3 + 1/(9k) <= 4/9 < 2/3,
    with equality iff k = 1.  QED.
    For vmin = a >= 3: tau <= 2/(3a) <= 2/9.  For a = 2 (b odd >= 5): tau <= 1/6 + 2/(3b),
    which is <= 11/42 < 4/15 for b >= 7, and tau(2,5) = 4/15.  So max over vmin >= 2 is
    4/15, only at (2,5).""")
cnt = 0
above = []
best2 = (Fraction(0), None)
for b in range(2, 401):
    for a in range(1, b):
        if gcd(a, b) != 1:
            continue
        cnt += 1
        t = tau((a, b))
        # closed form: first t >= 1/(3a) with ||t b|| >= 1/3
        x = Fraction(1, 3 * a)
        y = x * b
        fl = y.numerator // y.denominator
        fr = y - fl
        if fr < Fraction(1, 3):
            cf = (fl + Fraction(1, 3)) / b
        elif fr > Fraction(2, 3):
            cf = (fl + 1 + Fraction(1, 3)) / b
        else:
            cf = x
        check(t == cf, (a, b, t, cf))
        check(Fraction(1, 3 * a) <= t <= Fraction(2, 3 * a), (a, b))
        if t > Fraction(1, 3):
            above.append(((a, b), t))
        if a >= 2 and t > best2[0]:
            best2 = (t, (a, b))
check(all(v[0] == 1 and v[1] % 3 == 0 and t == Fraction(1, 3) + Fraction(1, 3 * v[1])
          for v, t in above), "tau > 1/3 classification")
check(len(above) == 400 // 3, len(above))
check(best2 == (Fraction(4, 15), (2, 5)), best2)
print(f"    exact check over all {cnt} coprime pairs with b <= 400: closed form OK, "
      f"tau in [1/(3a), 2/(3a)] OK,")
print(f"    the {len(above)} pairs with tau > 1/3 are exactly (1,3k), tau = 1/3 + 1/(9k) OK; "
      f"max over a >= 2 = {best2[0]} at {best2[1]} OK.")

# ------------------------------------------------------------------ P2
print()
print("=" * 78)
print("P2. n = 3 (threshold 1/4), triples containing 1:  max tau(1,b,w) = 7/16, ONLY at (1,3,12).")
print("    Inputs (exact): Safe_{1/4}(1,b) cap [0,1/2] for b = 2..7:")
for b in range(2, 8):
    print(f"      b={b}: {fmt(safe_half((1, b), 4))}")
print("""    PROOF (w > b throughout; danger set of w = open intervals of length 1/(2w)).
    A closed interval J inside Safe(1,b) with |J| >= 1/(2w) contains a safe point of w.
    b = 2: J = [1/4,3/8], |J| = 1/8 >= 1/(2w) for w >= 4, so tau <= 3/8; (1,2,3): tau = 1/4.
    b = 3: Safe(1,3) = {1/4} u [5/12,1/2].  4 does not divide w: tau = 1/4.  w = 4m,
           3 does not divide m: ||5w/12|| = 1/3, tau = 5/12.  w = 12k: 5/12 is killed and
           the next safe point of w is (5k+1/4)/(12k) = 5/12 + 1/(48k) <= 7/16, equality iff k=1.
    b = 4: J = [5/16,7/16]; if [5/16,7/16) had no safe point of w it would lie in one
           danger interval, forcing 1/8 <= 1/(2w), w <= 4: impossible.  So tau < 7/16.
    b = 5: J = [1/4,7/20], |J| = 1/10 >= 1/(2w) (w >= 6): tau <= 7/20.
    b = 6: J = [1/4,7/24], |J| = 1/24 >= 1/(2w) for w >= 12: tau <= 7/24; w = 7..11 below.
    b = 7: J = [9/28,11/28], |J| = 1/14 > 1/(2w) (w >= 8): tau <= 11/28.
    b >= 8: [1/4, 7/16) has length 3/16 >= 3/(2b), so it contains a full safe interval of b,
           of length 1/(2b) > 1/(2w), hence a safe point of w: tau < 7/16.   QED (with the
           finite residue below).""")
res = {w: tau((1, 6, w)) for w in range(7, 12)}
print(f"    finite residue tau(1,6,w), w = 7..11: {', '.join(f'{w}: {t}' for w, t in res.items())}")
check(all(t < Fraction(7, 16) for t in res.values()), "residue")
for w in range(4, 3001):
    t = tau((1, 3, w))
    if w % 4:
        exp = Fraction(1, 4)
    elif w % 12:
        exp = Fraction(5, 12)
    else:
        exp = Fraction(5, 12) + Fraction(1, 48 * (w // 12))
    check(t == exp, (w, t, exp))
print("    family check tau(1,3,w) in {1/4, 5/12, 5/12 + 1/(48k)} exactly as stated, w = 4..3000: OK")

# ------------------------------------------------------------------ P3
print()
print("=" * 78)
print("P3. Family F_N = {1,3,4,...,N} = AP{1..N} minus {2} (n = N-1 speeds, threshold 1/N).")
print("""    THEOREM.  For N >= 4, Safe(F_N) cap (0,1/2] = [tau_N, 1/2 - 1/(4N)] with
       tau_N = 1/2 - (N-2)/(2N k_o),  k_o = N (N odd), N-1 (N even);  for N = 3 the component is [4/9, 5/9].
    i.e. tau_N = 1/2 - 1/(2N) + 1/N^2 (N odd),  1/2 - 1/(2N) + 1/(2N(N-1)) (N even).
    PROOF.  (i) The AP {1..N-1} is tight with Safe = {a/N}: if ||kt|| >= 1/N for k <= N-1,
    the N points 0,t,..,(N-1)t are pairwise >= 1/N apart on a circle of length 1, hence
    equally spaced, hence t in (1/N)Z.  (ii) So if t is safe for F_N and ||2t|| >= 1/N then
    t = a/N, but then ||N t|| = 0: contradiction.  Hence every safe t has ||2t|| < 1/N, and
    with t >= 1/N (speed 1) this leaves t = 1/2 - s, 0 <= s < 1/(2N).  (iii) There, for k <= N,
    ks < 1/2, so odd k gives ||kt|| = 1/2 - ks >= 1/N iff s <= (N-2)/(2Nk) (binding at k_o),
    and even k gives ||kt|| = ks >= 1/N iff s >= 1/(Nk) (binding at k = 4).  QED
    COROLLARY.  For every n >= 2:  sup over n-sets of tau >= tau_{n+1} = 1/2 - 1/(2(n+1)) + O(1/n^2),
    so the first lonely time can be pushed to 1/2 - O(1/n).""")
for N in range(3, 41):
    v = tuple([1] + list(range(3, N + 1)))
    ko = N if N % 2 else N - 1
    tN = Fraction(1, 2) - Fraction(N - 2, 2 * N * ko)
    S = safe_half(v)
    if N == 3:
        exp = [(Fraction(4, 9), Fraction(5, 9))]   # the component straddles 1/2
    else:
        exp = [(tN, Fraction(1, 2) - Fraction(1, 4 * N))]
    check(S == exp, (N, S, exp))
    check(tau(v) == tN, N)
print("    exact check N = 3..40 (Safe cap (0,1/2] equals the stated interval): OK")
print("    tau_N for N = 3..12: " + ", ".join(
    f"{N}:{Fraction(1, 2) - Fraction(N - 2, 2 * N * (N if N % 2 else N - 1))}" for N in range(3, 13)))

# ------------------------------------------------------------------ P4
print()
print("=" * 78)
print("P4. Kill-and-replace lemma.  If T is tight with Safe(T) in (1/N)Z (N = |T|+1), s in T and")
print("    m = kN is not in T, then every safe t of v = (T - {s}) + {m} has ||s t|| < 1/N.")
print("    PROOF: otherwise t is safe for T, so t in (1/N)Z and ||m t|| = 0.  QED")
TIGHT = [(1, 2), (1, 2, 3), (1, 2, 3, 4), (1, 3, 4, 7), (1, 2, 3, 4, 5), (1, 3, 4, 5, 9),
         (1, 2, 3, 4, 5, 6), (1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12),
         (1, 4, 5, 6, 7, 11, 13), tuple(range(1, 9))]
ncase = 0
for T in TIGHT:
    N = len(T) + 1
    ST, D = safe_components(T)
    check(all(lo == hi and (N * lo) % D == 0 for lo, hi in ST), ("not on grid", T))
    for s in T:
        rest = [x for x in T if x != s]
        for k in range(1, 60 // N + 1):
            m = k * N
            if m in T:
                continue
            v = tuple(sorted(rest + [m]))
            comps, D2 = safe_components(v)
            check(comps, ("LRC", v))
            # Safe(v) cap {||s t|| >= 1/N} = Safe_{1/N}(v + {s}) must be EMPTY
            both, _ = safe_components(tuple(sorted(set(v) | {s})), N)
            check(not both, (T, s, m))
            ncase += 1
print(f"    exact check: {ncase} triples (T, s, m) from the {len(TIGHT)} tight sets, m = kN <= 60: OK")
print("    (checked as: Safe_{1/N}(v + {s}) is empty, while Safe(v) is nonempty)")

# ------------------------------------------------------------------ P5
print()
print("=" * 78)
print("P5. Families met in the census (exact checks):")
for k in range(1, 201):
    check(tau((1, 3, 12 * k)) == Fraction(5, 12) + Fraction(1, 48 * k), k)
print("    tau(1,3,12k) = 5/12 + 1/(48k), k = 1..200: OK (proved in P2)")
print(f"    Safe_(1/5)(1,3,4) cap [0,1/2] = {fmt(safe_half((1, 3, 4), 5))}")
for k in range(1, 201):
    check(tau((1, 3, 4, 5 * k)) == Fraction(2, 5) + Fraction(1, 25 * k), k)
print("    tau(1,3,4,5k) = 2/5 + 1/(25k), k = 1..200: OK  [5k kills 1/5 and 2/5 (multiples of 1/5);")
print("     the next safe point of 5k after 2/5 is (2k+1/5)/(5k) = 2/5 + 1/(25k), inside [2/5,9/20]]")

# ------------------------------------------------------------------ P6
print()
print("=" * 78)
print("P6. Good period q(v): two easy facts.")
print("""    (a) q is UNBOUNDED for every n >= 2: if L = lcm(1..Q) is a speed, every t = k/q with
        q <= Q has ||L t|| = 0, so q(v) > Q (e.g. v = (1, L) is primitive).
    (b) q(v) <= 2 vmax whenever LRC holds for v: delta is attained at some t* = m/(v_i+v_j)
        (A1), f(t*) = delta >= 1/(n+1), so t* is a good time with denominator <= v_i+v_j <= 2vmax.
    (c) den(tau) >= q(v) always (tau is itself a good time).""")
from chessboard_weave_20261006_lrc_core import good_period, primitive_sets, lcm_list
row = []
for Q in range(2, 13):
    L = lcm_list(range(1, Q + 1))
    q, t = good_period((1, L))
    check(q > Q, (Q, q))
    row.append(f"Q={Q}: q(1,{L})={q}")
print("    (a) exact: " + "; ".join(row))
worst = {}
for n, B in {2: 200, 3: 60, 4: 30, 5: 18}.items():
    best = (Fraction(0), None)
    for v in primitive_sets(n, B):
        q, t = good_period(v)
        check(q <= 2 * v[-1], v)
        r = Fraction(q, v[-1])
        if r > best[0]:
            best = (r, v, q)
    worst[n] = best
print("    (b) exact check q(v) <= 2 vmax on the required universes: OK; max q/vmax: " +
      "; ".join(f"n={n}: {b[0]} at {b[1]} (q={b[2]})" for n, b in worst.items()))

# ------------------------------------------------------------------ P7
print()
print("=" * 78)
print("P7. Independent re-computation of the headline tau values by two other methods")
print("    (brute grid 1/((n+1) lcm v) and the fixed-point chase):")
from chessboard_weave_20261006_lrc_core import tau_brute, tau_chase
for v in [(1, 3), (1, 3, 12), (1, 3, 4, 5), (1, 3, 5, 8), (1, 5, 6, 7, 8), (1, 3, 4, 5, 30),
          (1, 3, 4, 5, 7, 18), (1, 3, 4, 5, 6, 7), (1, 5, 6, 7, 8, 11, 13), (2, 5), (2, 5, 8)]:
    t1, t2, t3 = tau(v), tau_brute(v), tau_chase(v)
    check(t1 == t2 == t3, (v, t1, t2, t3))
    print(f"    tau{v} = {t1}  (sweep = brute = chase)")

print(f"\nTotal time {time.time() - T0:.1f}s")
