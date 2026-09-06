# Independent audit of concurrent synthesis reports

**Date: 2026-09-05. Status: VERIFIED proof review of incoming work.**
Inspected incoming commit `09c725138` and its integration through `7330df37a`.
This review accepts the arguments within the dependencies stated below; it
does not promote the authors' reserved canon files or claim new priority.
The source reports retain their own `PROVED CANDIDATE` or `CONDITIONAL`
labels. No exhaustive finite computation was rerun for this follow-on audit.

## 1. Cumulative signed-cycle gaps

Source: [D5/D6 closure](even-graph-d5-d6-closure-synthesis-sep05.md).
Independent reviewer: `wildcard_portfolio`.

The all-order induction and equality classification are accepted relative
to the report's exact finite bases and local comparisons L6/L7. The review
read the primary script/output and the inherited
[THM-4083 / balanced deletion](../../01-canon/theorems/THM-4083-even-graph-cumulative-d3-d4-spectral-gap.md)
and [THM-4078 / Fourier normalization](../../01-canon/theorems/THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation.md).
It independently checked the polynomial induction gaps, positivity after
`n=9+t`, the balanced-deletion trichotomy, and spectral multiplicities.

The critical transfer is exact restriction: each negative `k`-cycle belongs
to `binom(n-k,r-k)` induced `r`-vertex subgraphs. L6 and L7 therefore supply
`c6>=c4` and `c7>=2c4+c5` in the induction range. Keeping the cycle lengths
separate repairs the contribution missed by scalar deletion. Cases with
at least two balanced deletions give exactly a one-edge switching class;
one balanced deletion gives strictness; no balanced deletions give strictness
in the induction range. Equality has labelled multiplicity `binom(n,2)` and
relabelling quotient multiplicity one.

Consequence: the universal D5/D6 statements are accepted with their finite
exact inputs. THM-4416 was still `RESERVED / UNPROVED EMPTY STUB` at the
inspected revision and is not cited here as proved canon.

## 2. Width-two moment returns and parabolic dynamics

Source: [width-two parabolic seed](nc2-width-two-parabolic-seed-synthesis-sep05.md).
Independent reviewer: `lrc_anchor`.

The review accepts the algebraic argument and its stated classical input.
For `T=-R(0)z/R(z)`, the root-swap involution satisfies
`iota=T-b*z^(2m*+1)+...`; composing gives
`T^2=z-2b*z^(2m*+1)+...`. The two leading defects add, with the stated sign.
Thus the fixed-point multiplicity supplies exactly `m*` petal cycles.

The primary source was independently checked:
[Milnor, Dynamics in One Complex Variable, Section 10.3](https://abel.math.harvard.edu/archive/118r_spring_05/docs/milnor.pdf),
printed page 10-1 (PDF page 66), and its following sentence explicitly give
at least one distinct critical point for each petal cycle. Sections 7.4--7.5
exclude exact preimages of the fixed point. The same critical-point claim
appears in [Epstein, page 2](https://arxiv.org/pdf/math/9902158).
This is a `CITED` input, not a new result of either research session.

Poles map to infinity and then zero, so they are unavailable to the required
basins. At a nonzero root of `R` of multiplicity `e`, `R-zR'` has exactly
multiplicity `e-1`. Division by `gcd(R,R')` removes precisely the pole roots,
leaving a polynomial of degree equal to the number of distinct roots of `R`.
The resulting critical-point budget is sound. No unresolved mathematical
gap was found in this audit; the source's `CONDITIONAL` status is retained
until its author explicitly completes promotion.

The intersection with this session's trinomial result is precise. When the
sole negative charge is `-2`, collision forces `AB<=2`; hence `(A,B)=(1,2)`
and `a=AB`. These are `(-2,b,2b+2)` with odd `b`, already in the proved free
two-ray case. Reflected `(-c,1,2)` has `g=1` and singleton first return.
Dynamics adds arbitrary many terms at width two and a stronger root-count
bound; it does not resolve the all-trinomial two-rung coprimality question.

## 3. Precision, tensors, and the first dyadic mixed cluster

Source: [prime-wall synthesis, Sections 3 and 4A](confluent-twojet-prime-wall-synthesis-sep05.md).
Independent reviewer: `niche_moments`.

Both sections pass. DVR Smith normal form with largest exponent `L` gives
sharp observation precision `N+L` for coefficient recovery modulo `p^N`.
The finite kernel size is `p^(sum min(N,alpha_i))`. Invertible Smith changes
of basis tensor, so independent rectangular grids have all exponent sums.
The full source degree box and independent grid are essential hypotheses.

The dyadic formula has a direct audit independent of the 69-minor table.
After removing two unit pivots, use rows `(u,r)=(1,0),(1,1),(2,0),(2,1)`
and columns `q=2,3,4,5`. The residual integer entries are
`q^r*u^(q-r)` with scale weight `q-r`.

- One-minor valuation is `min(e+1,2e)`, attained by the column-two
  derivative and value entries.
- The unique weight-three two-minor uses the two derivative rows and
  columns two and three; its residual determinant is 12. Every other
  two-minor has weight at least four, with a unit weight-four witness.
  This gives `min(3e+2,4e)`.
- Both rows at `u=2` are divisible by four, so every three-minor has
  residual valuation at least two. Minimum weight is seven, attained
  with determinant 12 by columns two, three, four and rows value-one,
  derivative-one, derivative-two. Hence the third valuation is `7e+2`.
- The full determinant is `16*d^12`, where `d=2^e`, giving `12e+4`.

Successive differences prove all stated depths, not just the sampled
values of `e`. CRT for consecutive nodes `0,...,4` gives clusters `(0,2,4)`
and `(1,3)` and the full partition
`(0,0,0,0,2,2,2,2,5,7)`. Thus `p=2,n=p^2+1` is handled by the accepted
incoming proof. The next quadratic-boundary problem should specify odd
primes, or later dyadic mixed trees.

The pair `(0,8)` has fourth exponent eight, whereas adjoining sixteen
changes that exponent to seven. Same-residue adjoining therefore cannot
preserve an old Smith prefix; retain the full weighted minor filtration.

## Integration boundary

These reviews concern specific arguments and dependencies. They do not
turn reserved theorem IDs into proved entries, identify an external
novelty claim, solve LRC(14) or planar JC2, or replace the four primary
reports' independent proof and computational audits.
