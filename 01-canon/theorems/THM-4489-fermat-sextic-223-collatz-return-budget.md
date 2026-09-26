---
id: THM-4489
title: "223 is the largest projectively empty Fermat-sextic prime; Collatz 223-returns descend through halving budget eight and fail at nine"
status: "PROVED elementary transfers; CITED genus/Hasse-Weil plus FINITE-EXACT complete census; independently audited"
source: "crossroads223-20260926 geometry lane and independent root projective-chart audit"
depends_on:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md
proof: 05-knowledge/results/crossroads223_20260926_geometry.md
scripts:
  - 04-computation/experiments/crossroads223_20260926_geometry.py
  - 04-computation/experiments/crossroads223_20260926_sextic_audit.py
---

# THM-4489: Fermat sextics and 223 returns

**PROVED in the scopes below.** Reserved in `5a0338e38`, promoted after
the root independently checked every projective chart at all 84 primes
through 433. The curve theorem uses cited input and a complete finite
certificate. No claim of literature priority or Collatz convergence.

## Complete local classification

For C_p: x^6+y^6+z^6=0 in P^2(F_p),

    C_p(F_p) is empty iff p in {7,31,67,79,139,223}.
    C_p(F_p) has no point with xyz nonzero iff
      p in {2,5,7,13,31,61,67,79,97,139,157,223,277}.

Thus 223 is largest for projective emptiness, and 277 for torus emptiness.
For p>3 this is a smooth genus-ten curve. The CITED Hasse-Weil inequality
N_p>=p+1-20sqrt(p) is positive for p>=401, and exceeds the at most 18
coordinate-zero points for p>=439. Exact integer square comparisons
prove these cutoffs; the finite prime census completes both classifications.
At the singular characteristics the exact counts (projective,torus) are
(3,0) for p=2 and (4,4) for p=3; no smooth-curve bound is used there.

The primary source for genus and Hasse-Weil is
[Achinger, Lecture 21, Example 21.2.4 and Theorem 21.3.1/Remark 21.3.2](https://achinger.impan.pl/ag2026/lec21.pdf).
The [complete proof note](../../05-knowledge/results/crossroads223_20260926_geometry.md)
supplies the exact multiplicative-subgroup count. A second implementation
directly enumerates 4,633,895 affine-chart cases plus all points at infinity,
without importing the geometry script. Both lists agree.

## Exact Collatz transfer and sharp return budget

For U(n)=(3n+1)/2^k on positive odd integers, k=v_2(3n+1),
ord_223(2)=37 and ord_223(3)=222. Write H=<2> and C_a=3^aH, a mod6.
The unique forbidden nonzero coset transition is C_5 to C_3: it would
give u+v+1=0 for u,v in H, exactly a torus point on C_223.
The quotient forgets k's size and the source's integer height.

For a positive valuation word w=(k_1,...,k_r), S_i=sum_(j<=i)k_j,

    U_w(n)=(3^r n+C_w)/2^S_r,
    C_w=sum_(i=0)^(r-1)3^(r-1-i)2^S_i.

A 223-to-223 return requires 223|C_w. Exhausting all 255 positive
compositions with total cost at most eight leaves only (2,3,1,1) and
(2,3,1,2), with maps (81n+223)/128 and (81n+223)/256. Both strictly
descend for EVERY positive 223-multiple, at all integer heights.

At cost nine, (1,1,1,1,2,1,1,1) gives the expanding FIRST return

    14495 -> 21743 -> 32615 -> 48923 -> 73385 -> 55039
          -> 82559 -> 123839 -> 185759.

Hence the budget-eight result cannot be replaced by a universal return
drift claim. The missing coordinate is odd-step count versus total
halving count, not merely the divisibility of a carry.

The minimal carry 223=27+36+96+64 also gives the rational cycle
223/47,179/47,73/47,133/47, because 128-81=47. Scaling the additive
shift to 47 makes it an integer cycle of 3n+47. Carry 223 itself occurs
exactly for words (2,3,1,k), k>=1; none has its gap 2^(6+k)-81 dividing
223, so no integer 3n+1 cycle has this exact unreduced carry. This says
nothing about larger carries divisible by 223.

Reproduce the main note's controls with
`python 04-computation/experiments/crossroads223_20260926_geometry.py --rank`
and the independent census with
`python 04-computation/experiments/crossroads223_20260926_sextic_audit.py`.
Both pass normally and with `-O`; finite universes and hostile controls
are recorded alongside the full proofs.
