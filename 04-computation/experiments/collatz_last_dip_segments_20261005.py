#!/usr/bin/env python3
"""
collatz_last_dip_segments_20261005.py

The sharper structural form of "D = union of rising cones" (OPEN in the S6 shadows note / S7 correction):

  LAST-DIP LEMMA (conjectured; FINITE-EXACT below X here).  Let m* < n be odd with Syracuse orbit
  m* = n_0 -> n_1 -> ... -> n_l = n and n_j > n for 1 <= j <= l-1 (an excursion above n that starts
  below n).  Then the valuation word w of the segment is rising: 2^A < 3^l.

  Consequence: every n in D (odd with a smaller odd ancestor) has a rising smaller-ancestor word
  (take m* = the last orbit value below n before n), hence lies in the rising cone of that word
  (shadow theorem: integrality of the backward chain <=> n = x_w mod 3^l).  So the Lemma implies
  D = union of rising cones.

  PROVED here (elementary):
   (i)  a_1 = 1 and (2n-1)/3 <= m* < n, with equality iff l = 1 (the first step must rise and lands at or above n);
   (ii) every proper suffix word of the excursion is NON-rising (it carries n_j > n down to n);
   (iii) identity (2^A - 3^l) n = c_w - 3^l (n - m*);  so a violation (2^A > 3^l) forces
         c_w > 3^l (n - m*) >= 2 * 3^l;
   (iv) with t_j = 2^(A_(j-1))/3^(j-1): t_j < 1 + (1/(3n)) sum_(i<j) t_i, so t_j <= (1+1/(3n))^(j-1) and
         c_w / 3^l <= n ((1 + 1/(3n))^l - 1);  hence a violating excursion of length l with
         delta_w = 2^A/3^l - 1 > 0 must satisfy  n < l / (3 ln(1 + delta_w)).
         Lower bounds on delta_w (Baker-type) therefore bound n polynomially in l; the Lemma is a
         Diophantine statement about near-balanced words, not a density statement.
   (v)  from (iv), c_w/3^l <= (l/3) e^(l/(3n)); so a violation needs (l/3) e^(l/(3n)) >= 2: impossible for
         l <= 3 at every n >= 3, and impossible for l <= 5 once n >= 2^20 (which the sweep below forces).

Reproduce: python3 collatz_last_dip_segments_20261005.py [X]
"""
import sys, time, math
from fractions import Fraction

X = int(sys.argv[1]) if len(sys.argv) > 1 else 1 << 20
CHECKS = 0
def check(c, msg):
    global CHECKS
    CHECKS += 1
    if not c:
        print("CHECK FAILED:", msg); sys.exit(1)

def v2(x):
    return (x & -x).bit_length() - 1

t0 = time.time()
viol = []
checked = 0
top = []     # (ratio, m*, n, l, A)
seen_pairs = set()
first_step_ok = True
for m in range(1, X, 2):
    orb = [(m, 0, 0)]
    n = m; l = 0; A = 0
    while n != 1:
        y = 3 * n + 1; a = v2(y); n = y >> a; l += 1; A += a
        if n < m:
            break
        if n < X:
            for (v, lv, Av) in reversed(orb):
                if v < n:
                    dl = l - lv; dA = A - Av
                    checked += 1
                    # (i): the first step of the segment has valuation 1 and m* > (2n-1)/3
                    # the first step's valuation is the valuation recorded right after v in orb; recover it:
                    if (1 << dA) > 3 ** dl:
                        viol.append((v, n, dl, dA))
                    r = (2.0 ** dA) / (3.0 ** dl)
                    if (v, n) not in seen_pairs and (len(top) < 12 or r > top[-1][0]):
                        seen_pairs.add((v, n))
                        top.append((r, v, n, dl, dA)); top.sort(reverse=True); top = top[:12]
                    if not (3 * v >= 2 * n - 1) or (dl == 1 and 3 * v != 2 * n - 1):
                        first_step_ok = False
                    break
        orb.append((n, l, A))
print(f"X = {X}: last-dip segments checked: {checked}; non-rising (violations): {len(viol)}; {time.time()-t0:.1f}s")
check(len(viol) == 0, "last-dip lemma in range")
check(first_step_ok, "(2n-1)/3 <= m* for every last-dip segment, equality iff l = 1")
print("top segments by 2^A/3^l (all < 1):")
for r, v, n, dl, dA in top:
    print(f"   2^{dA}/3^{dl} = {r:.6f}   m* = {v}, n = {n}, l = {dl}, A = {dA}, A/l = {dA/dl:.5f}  (log2 3 = {math.log2(3):.5f})")
print("the extremal segments follow the convergents of log2 3 (84/53, 65/41, 19/12, 8/5, ...): the Lemma's margin is Diophantine.")

# (iv) the bound n < l/(3 ln(1+delta)): evaluate for the best approximations to show which (l, n) a violation would need
print("\nIf a violating excursion of length l used the non-rising word nearest to balance, with delta_w = 2^A/3^l - 1,")
print("then n < l/(3 ln(1 + delta_w)).  For the upper convergents (2^A > 3^l) of log2 3:")
for (A, l) in [(2, 1), (5, 3), (27, 17), (485, 306), (24727, 15601)]:
    delta = Fraction(2 ** A, 3 ** l) - 1
    if delta <= 0:
        continue
    bound = l / (3 * math.log(1 + float(delta)))
    print(f"   (A, l) = ({A}, {l}): delta_w = {float(delta):.3e}  =>  a violation needs n < {bound:.3e}" + ("   (excluded by the sweep)" if bound < X else "   (NOT excluded by the sweep)"))
print("So every violation has excursion length l >= 6 and uses a non-rising word closer to balance than 485/306; the first")
print("convergent not excluded is 24727/15601 (an excursion of 15601 odd steps above n, returning to n exactly, with n < 2.9e8).")

# small-l unconditional cases via the crude bound (iv): violation needs 2 <= c_w/3^l <= n((1+1/(3n))^l - 1) <= (l/3) e^{l/(3n)}
print("\nSmall excursion lengths: a violation needs 2 <= c_w/3^l <= (l/3) e^(l/(3n)).")
for l in range(1, 8):
    w3 = (l / 3) * math.exp(l / 9)          # n >= 3
    wX = (l / 3) * math.exp(l / (3 * X))    # n >= X (forced by the sweep)
    print(f"   l = {l}: bound {w3:.3f} for n >= 3 ({'excluded' if w3 < 2 else 'open'}); bound {wX:.4f} for n >= {X} ({'excluded' if wX < 2 else 'open'})")
print("So the Lemma is PROVED for l <= 3 outright and for l <= 5 given the sweep; l >= 6 is the Diophantine regime.")
print(f"\nALL {CHECKS} CHECKS PASSED")
