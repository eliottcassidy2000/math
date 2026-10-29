# Exact high rank certificates and the meaning of twenty one

**FINITE-EXACT rank lower bounds / CITED classical character lemma and
Mordell--Weil theorem / CONDITIONAL rank-distribution model.** Data retrieved
2026-10-29 from the ICARM Elliptic Curve Rank Leaderboard. This work
acknowledges ICARM and NSF Grant DMS 2425401. The published curves and points
are not discoveries of this packet.

## E1 Three rank certificates from small finite fields

The following independently recomputed certificates use exact rational
substitution and finite-field quadratic characters. No analytic rank,
floating regulator, GRH or BSD is used for these lower bounds.

| Public source | Verified point count and character rank | Largest retained prime | Prime excluding rational two torsion |
|---|---:|---:|---:|
| [ICARM 273](https://elliptic-rank.icarm.cloud/curve/273) | 30 | 257 | 59 |
| [ICARM 302](https://elliptic-rank.icarm.cloud/curve/302) | 31 | 311 | 59 |
| [ICARM 4239](https://elliptic-rank.icarm.cloud/curve/4239) | 32 | 199 | 7 |

The three local `icarm_*.json` snapshots contain the supplied data;
[rank_certificate.json](rank_certificate.json) pins them by SHA256 and lists
every retained (prime, root, binary row). These are rank **lower bounds**.
The commentary on curves 273 and 302 gives conditional exact ranks under
GRH and BSD. This packet neither establishes upper bounds nor saturation.

For a generalized Weierstrass equation let

\[
Y=2y+a_1x+a_3,\quad Y^2=g(x)=4x^3+b_2x^2+2b_4x+b_6,
\]

where b2=a1^2+4a2, b4=a1*a3+2a4, b6=a3^2+4a6. At a good odd prime p,
each root r of g gives the classical homomorphism E(Fp)->F2 whose value is
the quadratic character bit of x-r; at (r,0) use g'(r), and at infinity
use zero. Multiplication by 4 does not affect a quadratic character.
The script skips primes dividing a point's x denominator.

The homomorphism follows from the chord identity: for the monic cubic
f=g/4 and a line through P,Q,-(P+Q), evaluating f-line^2 at r shows that
the product of the three x-r values is a square. The derivative gives the
limiting value at the two-torsion point. This is the standard 2-descent
character used by the [ICARM verification method](https://elliptic-rank.icarm.cloud/).

A prime at which g has no root excludes rational two-torsion: good
reduction at an odd prime injects rational two-torsion. By Mordell--Weil,
dim_F2 E(Q)/2E(Q)=rank E(Q) in this case. The square binary matrix with one
column per supplied point is nonsingular, proving their independence
modulo torsion and the stated lower bound.

The verifier also computes complete addition tables at one single-root
prime and one three-root prime per curve, checking every character on
every pair of finite-field points (23500 pair checks in total). Deliberately
duplicating the last point column drops the rank by exactly one. A second
row-list elimination agrees with the incremental elimination.

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/rank_certificate.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/rank_certificate.py
```

## E2 The free bit is an explicit relation

At a completely split good prime, g=4(x-r1)(x-r2)(x-r3). On every curve
point the product of its three descent characters is one, including at
the roots when the derivative convention is used. Therefore

\[
\boxed{\ell_{r1}+\ell_{r2}+\ell_{r3}=0\quad\text{over F}_2.}
\]

The third character bit is free once the first two are known. Retaining
prime and root labels identifies which relation justifies the omission.
This is the same linear dependence mechanism used in the repo's
square-class/Hamming banks. It does not identify their underlying numbers,
curve points, or operations. A bank certifying r independent columns needs
at least r binary rows; our retained banks attain that row count. This
does not minimize the number of primes or the computation cost.

The two square-root choices for y have the same character response.
That loss is harmless for linear independence modulo two, but can destroy
the positivity of a fruit triple under an elliptic parametrization. Thus
even a minimal exact certificate for one predicate need not preserve a
different target.

## E3 Twenty one is a heuristic threshold rather than a maximum

[Park--Poonen--Voight--Wood](https://arxiv.org/html/1602.01431v3) propose a
random alternating-matrix model predicting only finitely many rational
elliptic curves of rank greater than 21. This is not a theorem about all
curves or a claim that 21 is the typical rank. Its exponent balance is

\[
H^{5/6}\,H^{-(r-1)/24}=H^{(21-r)/24}.
\]

The first factor counts curves by height; the second is the model's rare
rank tail. A finite collection of rank-30,31,32 curves does not logically
contradict that prediction. The paper explicitly discusses limitations
for selected families. Turning a curated search record into a frequency
estimate would require its selection rule and height denominators.

The supplied ICARM curve 273 and curve 302 have published rational-parameter
specializations in their commentary. Those formulas preserve useful family
provenance; their generic rank and the frequency of specialization jumps
are separate questions from the verified ranks of these two fibers.
A promising test is to certify common sections over Q(T), then record
which additional independent character columns appear at each fiber.

For the Paley tournament on seven points, the number 21 is the order of
its affine automorphism group, 7*3. The arithmetic model's 21 comes from
the displayed height exponents. No map equating these invariants has been
constructed. The substantive bridge is the use of alternating forms and
retained realization data, not equality of the two integers.

## Recovered fruit equation boundary

The September [fruit-rank packet](../collatz_mod6_20260921_fruit_rank_positive_multiples.md)
provides the exact projective map

\[
[a:b:c]=[56Z-X+Y:56Z-X-Y:-56Z-12X]
\]

from y^2=x^3+109x^2+224x to the cubic underlying
a/(b+c)+b/(a+c)+c/(a+b)=4. It retains the nonzero pair-sum guards and
the real positivity window. Its rank-one calculation, six torsion cosets,
and positive-multiple census are inherited results, not new high-rank
computations here. The canonical hostile is positive 9G whose three
ternary children 26G,27G,28G are all nonpositive. Any proposed synthesis
must preserve that real chamber in addition to its finite-field bits.
