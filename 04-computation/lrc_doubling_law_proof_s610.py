#!/usr/bin/env python3
"""lrc_doubling_law_proof_s610.py -- toward a PROOF of the doubling-sporadic
mod-3 law (HYP-2169), with each step rigorous or explicitly verified.

THEOREM (target). For even n >= 6, with V* = AP[(n-2)->2(n-2)] =
{1,...,n-3, n-1, 2n-4}, the loneliness radius satisfies
        M(V*) = 1/n  (tight)   <=>   3 | (2n-1).

PROOF STRUCTURE (labels: [PROVED] rigorous; [VERIFIED] checked n<=40, not yet
general):

  STEP 1 [PROVED]  M(V*) >= 1/n  for all even n >= 6.
    Witness t = 1/n. Then ||v/n|| = dist(v mod n)/n; the antipodal pair {1, n-1}
    binds at exactly 1/n, every AP runner v has dist(v mod n) >= 1, and the new
    runner satisfies ||(2n-4)/n|| = dist((n-4) mod n)/n = min(4, n-4)/n >= 1/n
    (n >= 6). The removed runner n-2 is irrelevant. So min_{v in V*} ||v/n|| =
    1/n, hence M(V*) >= 1/n.  (n=4 is degenerate: 2(n-2)=4=n.)

  STEP 2 [VERIFIED]  M(V*) = max(1/n, P),  P := max_m min_{v in V*} ||v m/(2n-1)||.
    The ONLY pinch family that can exceed 1/n is the pair-sum lattice t=m/(2n-1)
    of the mirror pair {3, 2n-4}. (Verified for all even n in [4,40].)

  STEP 3 [PROVED]  2n-4 = (2n-1) - 3 == -3 (mod 2n-1).
    So on the (2n-1)-lattice, runner 2n-4 mirrors runner 3:
    ||(2n-4) m/(2n-1)|| = ||3 m/(2n-1)||. The pair pinch distance is ||3m/(2n-1)||.

  STEP 4 [PROVED, number theory]  The pair {3, 2n-4} attains distance 2/(2n-1)
    for some m  <=>  2 in 3*(Z/(2n-1))  <=>  gcd(3, 2n-1)=1  <=>  3 does NOT
    divide 2n-1.  (3 reaches residue 2 mod (2n-1) iff gcd(3,2n-1) | 2 iff gcd=1,
    since gcd in {1,3} and 3 does not divide 2.)

  STEP 5 [VERIFIED]  If 3 does not divide 2n-1:  P = 2/(2n-1) > 1/n, so V* is
    LOOSE with M(V*) = 2/(2n-1). (At the m of Step 4 the other runners stay
    >= 2/(2n-1); verified n<=40. 2/(2n-1) > 1/n since 2n > 2n-1.)

  STEP 6 [VERIFIED]  If 3 | (2n-1):  P <= 1/n, so M(V*) = 1/n (TIGHT). On the
    (2n-1)-lattice, whenever the pair {3,2n-4} is far (>= 3/(2n-1)), some other
    runner is within 1/n (the multiple-of-3 shell, S592/S593, supplies a blocker;
    e.g. the runner (2n-1)/3 in the AP). Verified n<=40.

  Combining 5,6 with 1,2:  M(V*)=1/n  <=>  P<=1/n  <=>  3|(2n-1).  QED modulo
  the two verified covering lemmas (Steps 2,5,6).

This file runs each step's verification.
Session: claude-2026-06-03-S610 (lrc-doubling-law-proof).
"""
import sys
sys.stdout.reconfigure(line_buffering=True)
from fractions import Fraction as F
from itertools import combinations
from math import gcd

def dist(x): x = x % 1; return min(x, 1-x)
def Vstar(n):
    AP = list(range(1, n)); a = n-2
    return tuple(sorted([x for x in AP if x != a]+[2*a]))
def M_full(V):
    cand = {F(0)}
    for a, b in combinations(V, 2):
        for m in range(a+b+1): cand.add(F(m, a+b))
    return max(min(dist(v*t) for v in V) for t in cand)
def P_lattice(V, n):
    return max(min(dist(v*F(m, 2*n-1)) for v in V) for m in range(2*n))

NS = list(range(6, 41, 2))
print("\n  TOWARD A PROOF OF THE DOUBLING-SPORADIC MOD-3 LAW\n" + "=" * 70)

# STEP 1
ok1 = all(min(dist(v*F(1, n)) for v in Vstar(n)) == F(1, n) for n in NS)
print(f"\n  STEP 1 [PROVED]  M(V*) >= 1/n via t=1/n witness (n=6..40): {'OK' if ok1 else 'FAIL'}")

# STEP 2
ok2 = all(M_full(Vstar(n)) == max(F(1, n), P_lattice(Vstar(n), n)) for n in NS)
print(f"  STEP 2 [VERIFIED]  M(V*) = max(1/n, P_2n-1) (n=6..40): {'OK' if ok2 else 'FAIL'}")

# STEP 3
ok3 = all((2*n-4) % (2*n-1) == (2*n-1-3) for n in NS)
print(f"  STEP 3 [PROVED]  2n-4 == -3 (mod 2n-1): {'OK' if ok3 else 'FAIL'}")

# STEP 4
def reaches2(n):
    M = 2*n-1
    return any((3*m) % M == 2 or (3*m) % M == M-2 for m in range(M))
ok4 = all(reaches2(n) == (gcd(3, 2*n-1) == 1) for n in NS)
print(f"  STEP 4 [PROVED]  pair reaches dist 2/(2n-1) <=> gcd(3,2n-1)=1 <=> 3 nmid 2n-1: {'OK' if ok4 else 'FAIL'}")

# STEP 5 & 6 and the law
print(f"\n  STEPS 5-6 + LAW  [VERIFIED]")
print(f"  {'n':>3} {'3|2n-1':>7} {'M(V*)':>7} {'1/n':>6} {'2/(2n-1)':>9} {'tight':>6} {'law ok':>7}")
ok_law = True
for n in NS:
    M = M_full(Vstar(n)); div = (2*n-1) % 3 == 0; tight = (M == F(1, n))
    law = (tight == div); ok_law &= law
    pred = F(1, n) if div else F(2, 2*n-1)
    flag = '' if M == pred else '  <-VALUE?'
    print(f"  {n:>3} {str(div):>7} {str(M):>7} {str(F(1,n)):>6} {str(F(2,2*n-1)):>9} "
          f"{str(tight):>6} {str(law):>7}{flag}")
print(f"\n  => M(V*) = 1/n  <=>  3 | (2n-1)   [law holds n=6..40: {ok_law}]")
print("""
  STATUS: Steps 1,3,4 are rigorous; Steps 2,5,6 are verified covering lemmas
  (n<=40). The law is thus reduced to two explicit, finite-flavored covering
  statements about the (2n-1) pinch lattice of V* -- a concrete proof target,
  not an open-ended search. Loose value M(V*)=2/(2n-1) is a clean by-product.
""")
