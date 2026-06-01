---
id: HYP-2045
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S547
related:
  - HYP-2044
  - HYP-2008
  - HYP-2039
---

# HYP-2045: at n=2q the (q,q) cycle type is the loneliest config's maximal-odd automorphism; the apex is the marked reflection's co-observer

**Answer to HYP-2044's open question (vigorously, verified q=3,5,7,11):** YES, with
precision.

**VERIFIED (`lrc_n2q_qq_cycle_control_s547.py`):**
- The loneliest config at n=2q is the regular 2q-gon (roots of unity), t*=1/(2q),
  collar=1/n; symmetry group D_{2q}.
- (A) MARKED reflection v<->2q-v fixes exactly {0, q} = OBSERVER and APEX. So the apex
  (speed q=n/2, order-2 element, position 1/2) is the unique runner fixed by the same
  Z_2 that fixes the observer -- the antipodal co-observer (= why the apex is special).
- (B) rotation-by-2 of the 2q-gon has cycle type exactly (q,q) (even vs odd positions)
  AND preserves the half-turn relation = a genuine AUTOMORPHISM of the loneliest
  tournament. It is the MAXIMAL odd-order rotation (rotation-by-1 = a 2q-cycle = even =
  a Burnside 0). Observer in the even q-cycle, apex in the odd q-cycle, bracketing
  runners +-1 (odd) in the apex's cycle.
- BURNSIDE TIE (= S546): the equal pair (q,q) gives the maximal between-cycle term
  gcd(q,q)=q; that very (q,q) is the loneliest config's automorphism. Doubled-prime
  Burnside dominance and the LRC extremiser's symmetry are ONE fact. All carry q
  because apex = n/2 = order-2 element = (q,q) cycle length = gcd = q.

**HONEST NEGATIVE:** the apex is NOT the cascade-trap runner. The zero conditional
clearance (S545) falls on the LAST runner 2q-1, never the apex q (verified). The apex's
role is geometric (reflection co-fixed point), not the cascade bottleneck.

**OPEN:** (A) does (q,q)-automorphism + the marked reflection uniquely single out the
regular 2q-gon as the tight observer-source config at n=2q (symmetry-forces-extremiser
=> pins AP-tightness at the doubled primes)? (B) the trap = last runner 2q-1 = -1 mod 2q
= reflection-image of 1; relate to the apex.

**Files:** `04-computation/lrc_n2q_qq_cycle_control_s547.py` (+.out). Reflection:
`07-reflections/lrc-n2q-the-qq-cycle-is-the-loneliest-configs-symmetry-s547.md`.
