# Audit, recovered mechanisms and exact boundaries

## Proof review

This is a self-audit with independent computational paths, **not** an
independent-agent audit. No theorem ID is reserved or promoted into canon.

1. **Types and ordering.** Integer carry `C`, normalized carry `B`, size
   `J`, valuation clock `D` and branch count `m` are distinct. Standard
   Collatz clock is `D+J`. Children in the grammar are ordered; a
   depth-first drawing is not the previous breadth-first drawing.
2. **Source guards.** Formula (2) specifies an odd cylinder modulo
   `2^{D+1}`, hence normalized odd-Haar mass `2^{-D}`. Finite moments
   are summaries, whereas `(J,D,C)` retains the source word. No independent
   sampling is introduced when conditioning to an original cylinder.
3. **Finite versus infinite expectations.** Slope tilting is explicitly
   a change of probability law. Subcriticality of the tilted tree proves
   total tilted mass one. Nonnegative finite-tree sums justify every
   moment limit, including divergence. Signed return atoms are never
   treated termwise positive.
4. **Cubic branches.** Every generating function is the minimal
   nonnegative branch obtained by iteration. Other positive cubic roots
   are not substituted for it. At the slope-moment double root the
   expectation is finite; the threshold is inclusive.
5. **Exponential budget.** Rational supersolutions at `t=4/3` establish
   convergence before the carry resolvent is divided. Conditioning costs
   exactly the factor permitted by the mass of the cylinder. The retained
   clock is not silently replaced by a tree depth or branch budget.
6. **Target.** Bounds on a source measure and on an auxiliary stopping
   rule do not prove descent, universal entry or exclusion of a null set.
   Both necessary/sufficient mistakes have explicit controls in the notes.

## Recovered niche work and the connection contracts

The search started from filename families (continuants, renewal, cocycles,
moment flows, reciprocal jets) and followed dependencies. It did not treat
an exhaustive chronological log as current truth. An October reference
search found THM-3186 in the 105 atlas, THM-3051 in the recent cursor/CM
packet, and THM-1770 in two October Collatz notes. Thus these are recovered
mechanisms, not claims of never-rediscovered results.

| Source | Map into this packet | Preserved predicate | Lost information / needed sidecar | Cheapest decisive test |
|---|---|---|---|---|
| [THM-3186, full exterior continuant](../../../01-canon/theorems/THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall.md) | ordered hidden-state exit sum becomes ordered affine carry sum | exact composition of an extension coordinate | its exterior geometry is not identified with Collatz; keep `C` | equal clocks, distinct carries |
| [THM-3192, reciprocal coefficient jets](../../../01-canon/theorems/THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return.md) | stable coefficient filtration inspires the moment-order filtration in (7) | exact recurrence after a declared truncation | initial boundary history survives; deeper jets can restore a hidden mode | compare full distributions with finite jets |
| [THM-3051, Stieltjes multiplier/Gamma flow](../../../01-canon/theorems/THM-3051-stieltjes-multiplier-gamma-flow-and-moving-lower-hankel-boundary.md) | retain one measure when taking moments; name its slope tilt explicitly | positivity under an actual positive pushforward | a family of good marginal laws is not one common source law | original mean clock infinite, tilted mean finite |
| [THM-1770, corrected first-return renewal](../../../01-canon/theorems/THM-1770-the-localisation-lemma-as-first-return-renewal-and-the-pair-only-closure.md) | word-language factorization, with all contributions retained | exact free-word grammar | scalar cancellations do not isolate atoms | use signed two-atom cancellation as hostile |
| [September Catalan/Ramsey note](../collatz_catalan_ramsey_20260927.md) | binary first passage generalized to ordered `d`-ary first extinction | unique first-return parsing | its orbit tournament was an invalid carrier; keep two orders as inversion graph | `(1,3,2)` |
| [THM-2966, July critical-run fair extractor](../../../01-canon/theorems/THM-2966-spine-normal-form-for-critical-run-fair-extractors.md) | explicit conditional prefix shares and retained boundary terms | exact source probabilities under refinement | its unknown-bias fairness and pathwise deadline are not our fixed Haar law or averaged budget | a null infinite ray survives decreasing cylinder mass |

THM-3192 is stronger than the bare continuant recurrence for its own
coefficient-degree problem: it already supplies a moving reciprocal gauge,
Smith layers and the reappearance of the hidden wedge. We inherit that
warning rather than claim its homogeneous shallow-jet recurrence as new.
Its source term vanishes only after the stated jet/index boundary; its
initial jets still contain earlier source history. Our clock truncation
is likewise legitimate for finite coefficients, not for deleting all
high-clock exceptional sources.

## Independent computations

`check_recurrences.py` uses standard-library integer/Fraction arithmetic.
It enumerates **all** first-extinction valuation words with `D<=16`,
with no source-descent or carry filter. There are 626 words for `(3,2)`.
It compares their full coefficient dictionaries with (7) for carry
moment orders 0 through 4; checks the compact first-jet formula and
shifted-product associativity; checks each word's two positive integer
source representatives directly; and checks every chronological cut in
the concatenation law.

A separate enumeration of ordered tree shapes integrates every infinite
letter-gap sum exactly. This checks the mass, slope **and normalized
carry** coefficients of (10)–(13), through branch count 5 for `d<=3` and
through 3 for larger arities. Tested rules are `(3,2),(3,3),(7,3),(2,2),
(1,2),(15,4)`. These checks are independent of the bivariate jet algebra.

The incoming extension checks all odd arities 3 through 101: exact least
residue clock, original probability `1/d`, base-four tilted probability,
and the strict bounds (27). Ordered-shape integration checks the new
five- and nine-child carry recurrences through branch count 3; direct
valuation-word enumeration through `D=16` checks their clock grammars
and native integer guards. The plus-one arity identity is checked for
`r=1,...,8`. The five-child second-moment supersolution and the coarse
two-branch budget obstruction are exact rational/integer checks.

Exact rational supersolutions check the refinement constants. Native
source-cylinder intersections through `K=6`, with all words `D<=16` for
`(3,3)`, check conditional carry generating-function bounds. The proof,
not that finite census, supplies the unbounded tails.

Decimal bisection with 70-digit precision probes (18) for arities 3,7,15
at `1-t=10^{-4},10^{-8},10^{-12}`. These numerical limit probes are
**VERIFIED**, not exact asymptotic certificates; (18) has its own proof.
The receipt records the exact check count and all parameters. Checks
remain active under `python -O`.

## Historical corrections

The Catalan note's two false implications were checked independently:
all carry chains for doubling `n<=1000` in base three; and the proposed
value-order orientation on all 5913 permutations through order 7. Its
smallest hostile sequence `(1,3,2)` already decides the claimed tournament
equivalence. The old script computed LIS/LDS, not tournament transitivity;
those numerical subsequence values survive and are relabelled accurately.

The source note, producer, one synthesis sentence and mistakes ledger are
repaired. The script's non-certified-fraction comparison also had a
factor-two normalization error. `historical_checks.txt` is a fresh replay
of only the affected producer sections 1,2,5; the unmodified original full
output is historical provenance. No new literature claim or revalidation
of the unrelated Ramsey computations is implied.

## Reproduction

`python3 -B reproduce.py` replays the checker in normal and optimized modes,
compares both JSON objects with `checks.json`, and reproduces the affected
historical sections. No symbolic package, network or random seed is needed.
All-size proofs and open directions are in RECURRENCES.md, MOMENTS.md and
REFINEMENT.md, with the incoming-work extension in ODD_ARITIES.md. A finite
receipt is never used to establish an all-size claim.
