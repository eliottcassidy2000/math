---
id: THM-4490
title: "Uniform growing-prefix Collatz multiplicative independence from incident carry primes and S-unit bounds"
status: "PROVED application of a CITED S-unit theorem; independently audited; FINITE-EXACT controls"
source: "crossroads223-20260926 geometry lane, root refinement and independent automata audit"
depends_on:
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md
proofs:
  - 05-knowledge/results/crossroads223_20260926_geometry.md
  - 05-knowledge/results/crossroads223_20260926_geometry_growth.md
scripts:
  - 04-computation/experiments/crossroads223_20260926_geometry.py
  - 04-computation/experiments/crossroads223_20260926_geometry_growth.py
  - 04-computation/experiments/crossroads223_20260926_geometry_growth_audit.py
---

# THM-4490: affine-word S-unit specialization

**PROVED with cited input, independently audited.** Reserved in `5a0338e38`;
the promoted statement includes the subsequent incident-bank and factorial
refinements. The dependency on THM-4476 is only its exact Terras bijection.
This is a theorem about the distribution of starting integers. Collatz
and G2's independence on every injective chronological orbit remain OPEN.

## Exact bound and growing horizon

For odd n>0 define m_0=n, m_(j+1)=(3m_j+1)/2^k_(j+1), with the exact
positive valuation k_(j+1)=v_2(3m_j+1). Let B_N(I) count odd n in I for
which 3m_0,...,3m_(N-1) are multiplicatively dependent. For N>=2 put

    P_N=96^(N(N-1)/2),
    a_N=max{a>=0:a!<P_N},
    E_N=binom(3N,N)binom(N,2)2^(16a_N+40).

For EVERY interval I of H consecutive positive integers,

    B_N(I) <= E_N+(H/2)(27/32)^N+(27/4)^N.              (1)

In particular, with Q=log_2 H and

    N(H)=floor(sqrt(Q log_2 Q)/8),

the fraction of odd sources in I having full first-slot multiplicative
rank N(H) tends to one, uniformly in the location of I. More precisely,

    B_N(I) <= (H/2)(27/32)^N+H^(beta+o(1))=o(H),
    beta=(5+log_2 3)/8=0.82312031259...<1.              (2)

Full first-slot independence implies independence of the actual pairs
(-3m_j,3m_j+1). First-slot signs impose only two-torsion; rational
valuation rank detects independence because any integer kernel can be
doubled. No injectivity assumption on the selected words is used.

## Mechanism and proof

Fix a length-N positive valuation word, total S_N. Its sources form one
exact class a mod2^(S_N+1), and its nodes are affine forms

    L_i(t)=A_i t+B_i,
    A_i=3^i 2^(S_N+1-S_i), B_i=(3^i a+C_i)/2^S_i.

The normalized carries C_i/3^i strictly increase, so the forms have
distinct roots. Their nonzero determinants are

    A_j B_i-A_i B_j=-3^i 2^(S_N+1-S_j) C_(i,j),

where C_(i,j) is the subword carry. Let S_i^prime contain 2,3 and the
prime divisors of all carries incident to i. A prime outside this set
that divides L_i divides no other node. Therefore every supported node
of a multiplicative relation must be an S_i^prime-unit. There are at
least two supported nodes, because 3L_i>=3.

For such a pair i,j the identity

    A_j L_i/D_ij + (-A_i L_j/D_ij)=1

lies in the product of two rational unit groups of rank s_i+s_j. The
map t to this solution is injective. **CITED:** Beukers--Schlickewei,
[Theorem 1.1, author-hosted paper](https://webspace.science.uu.nl/~beuke106/s-units.pdf),
bounds solutions in the rational closure of a rank-r subgroup of (C*)^2
by 2^(8r+8). Thus each fixed word has finitely many exceptions.

For S_N<=3N, a subword of length r and total valuation K has carry
C<2^K(3/2)^r, coprime to 6. The product Q_i of all incident carries
is less than 96^(N(N-1)/2). If it has t_i distinct prime factors,
t_i!<=Q_i, so s_i<=2+a_N. Summing the S-unit bound over node pairs
and all binom(3N,N) words gives E_N, a GLOBAL count independent of
source heights.

The discarded event S_N>3N, with time-zero parity fixed odd, is exactly
sum_(j<N)binom(3N,j) classes modulo2^(3N+1). Each class meets I at most
H/2^(3N+1)+1 times. The binomial bound by (27/4)^N gives the last two
terms of (1). This is exact residue counting, not a temporal randomness
assumption. Finally a_N~[(5+log_2 3)/4]N^2/log_2 N gives (2).

The [complete proof and conditioned refinements](../../05-knowledge/results/crossroads223_20260926_geometry_growth.md)
include all strict inequalities, prime-bank types, sign constraints,
source collisions and the odd-source normalization. The
[fixed-word proof and high-height probes](../../05-knowledge/results/crossroads223_20260926_geometry.md)
give additional background without a novelty claim for basic S-unit
specialization.

## Conditioned no-descent candidates

Let D(X) contain n<=X with no descent through floor(log_2 n) half-steps.
For N=N(X) above, the number in D(X) with deficient first-slot rank is

    at most X^((3log_2 3-1)/8+o(1)) = X^(0.46936094...+o(1)).

This follows without a lower bound on #D(X): except for exp(O(N)) small
sources, no descent forces total valuation S_N<=ceil(N log_2 3)+1.
Otherwise that many half-steps contain at most N odd letters, have
multiplier <=1/2 and additive carry smaller than n/2, a contradiction.
Repeating the incident-bank bound at this smaller valuation budget gives
the exponent. The complete note states its exact integer finite bound.

Combining with the repaired THM-4487 no-descent count X^(0.949956...+o(1))
shows that almost all these candidates already have full growing-prefix
rank. This is a limitation of rank as a standalone filter. The absolute
failure bound is uniform on translated intervals of length X; the
relative statement uses the no-descent count on [1,X] and is not asserted
in every translated interval.

## Limits and validation

The total exception bound in (2) has leading power H^1 because of its
valuation tail; it is not a fixed power saving. Neither independence
nor density one forces descent on a particular integer orbit. Indeed,
n=(4^a-1)/3 goes directly to 1, so repeated terminal nodes give rank
failure at arbitrarily large starting values when words are allowed to vary.
This is outside G2's injective-orbit scope and explains why its tail
cannot simply be discarded.

The independent audit checked the local unit-group rank s_i+s_j, the
factorial inversion, and the L+1-bit Terras normalization. Exact controls
exhaust 3,597 words with N=2..5 and total<=3N, plus independent direct
residue censuses; high-height controls certify 60 prescribed cylinders
through 222 odd transitions and 617-bit nodes. The cited unit theorem,
rather than these finite controls, supplies finiteness at all heights.
