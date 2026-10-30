# Same-agent audit and finite exact receipts

**PROVED displayed arguments / FINITE-EXACT native checks.**
No independent-agent audit or formal theorem build was performed.
The all-size claims rely on the proofs, not on extrapolating a finite census.

## Universes and independent paths

The verifier uses only the Python standard library and a pinned inherited
81/27 helper. Its local guards raise exceptions and remain active under
Python -O. Normal and optimized receipts are byte-identical.

| Experiment | Exact universe | Outcome |
|---|---|---|
| Signed-affine classification | All simple labelled supports n=1..4; every unordered pair absent or s=+/-1,b=0/1 | 15,756 systems:4504 integral,2128 half-integral-only,9124 real-infeasible |
| Independent real feasibility | Tree propagation versus Fraction Gaussian elimination on every system | Agreement throughout; mod-two assignments separately enumerated |
| Canonical double covers | Every real-compatible system in that census | Floor/ceiling assignment verified on actual lifted edges |
| Matching covers | All presentations on n=3 with the same edge alphabet, every edge matching in S_d | d=2:729; d=3:15625; all odd-cover equivalences hold |
| Even-cover hostile | Same d=2 census |64 presentations repair a base integer obstruction;224 remain infeasible;441 already base-feasible |
| Arbitrary-cover compression | Every integral cover solution produced in the cover census | Fibre extrema satisfy the canonical signed double cover |
| Exact all-modulus counts | All131 systems n=1..3, every q=2..12 |1441 formula-versus-full-assignment comparisons agree |
| Local marginals | Boolean disequality triangle, denominators1..13; a denominator5 signed-triangle example with independently reordered local tuple lists | Half-denominator hostile; odd local medians synchronize shared variables |
| Positive-face criterion |625 simplex moment systems:three coefficients and target in[-2,2]; odd denominators through31 |377 odd-feasible,48 real-only,200 real-infeasible |
| Sharp denominator exponent | Same625 systems; all denominators1..32 enumerated | Among425 real-feasible systems, exponents0/1/2 occur377/36/12 times, agreeing with the exact content formula |
| Convex median bound |18165 cost evaluations on multisets of size1,3,5,7,9 from[-2,2], shifted absolute/squared/hinge costs | Bound and sharp examples verified |
| Marked points | Strict-majority samples over the same finite alphabet; one marked zero among balanced +/-1 fibres | Majority preservation and unbounded factor-d hostile |
| Involution functions | f,r=0,1,2 and q=2,3,4,5 |36 complete function spaces agree with gcd(q,2)^f q^r |
| Integral parity gluing | All pairs(a,b) in[-5,5]^2 | Same-parity sum/difference reconstruction |
| Native arithmetic hostiles |620 fixed-clock triples of positive odd integers below64; mixed-clock(1,3,5);5x+1 cycle; wrapped phases; mixed squares; injectivity | Fixed-clock median works; each unguarded extension fails as stated |
| Higher clock obstruction | Cyclic q-state triangle q=2..12, denominator divisibility candidates through32 and explicit q-sheet lift | Exact minimum-degree mechanism retained; this is algebraic candidate checking, not enumeration of all covers |
| Recovered odd cycles | Four explicitly pinned five/six-vertex tournaments | Independent subset-path DP and literal cycle enumeration replay H=11/13 and the H9 coefficient collision |
| Labelled tournament majority | Three transitive named tournaments on three vertices | Majority cyclic; unlabelled class does not determine the result |
| Hamming affine fibres | All eight syndrome classes of16 masks; all6528 unordered triples with repetitions | XOR preserves every class; masks0,7,25 give a median hostile |
| Late incoming ordered carry | Two seven-step words, D=-139 and N0=2363/2123 | Rational equations verified; one integral and one nonintegral, both with odd denominators |
| Mixed cyclic clock | One permutation with cycles of lengths2 and3 | No fixed point; covers of degrees2 and3, so every primewise minimum zero does not imply a base solution |
|81/27 control | Inherited native helper, separately reported | Common labelled cyclic triples67; group orders81/27; regular mask counts1620/1728 |

The local counter records **196,773 checks**; this is a count of guard
evaluations, not a number of independent mathematical examples. Internal
checks of the inherited81/27 routine are not included in that count.

The positive-face experiment uses a specialized two-row gcd formula.
There is no generic Smith-normal-form implementation in verify.py.
For arbitrary matrices, P4–P4b are proved symbolically. Their correctness
must not be presented as an exhaustive computational audit of generic
Smith reduction or of all finite constraint languages.

## Scope and adversarial audit

- P1 fixes a function on each edge and copies unary restrictions to every
  sheet. It does not silently reorder or identify named ports.
- P1's minima/maxima compression is existential for a fixed network.
  Later support queries, joint distributions, counts and group actions
  need more information.
- P2's half-integrality is a consequence of slopes+/-1. General rational
  affine slopes, arbitrary modular rotations and unfixed square labels
  are outside that theorem.
- The interval extension is proved in P2 but was not separately exhausted
  in the native checker. Arbitrary unary holes/global distinctness are
  not promoted into the real-plus-mod-two criterion.
- P3 assumes actual allowed tuples and exact common marginals, not only
  pairwise feasibility or positive masses at selected variables.
- P4 uses the maximal support face. Its explicit negative-coordinate
  hostile blocks dropping nonnegativity before the Smith calculation.
- P4b classifies the least two-adic denominator exponent, not the least
  total denominator or smallest connected cover.
- P5 applies to additive nonnegative convex unary costs. Negative weights,
  joint objectives, marked-point indicators and unproved native potential
  convexity are excluded.
- P6 specifies the involution or restriction/transfer identity. Odd
  vertex count and odd automorphism-group order are not substituted
  for either.
- The signed-affine proof follows tree expressions through both positive
  zero rows and negative pin rows. Testing parity alone without real
  compatibility would be false.
- Modular all-order counts retain Delta, the negative pin and whether a
  component is balanced; a Gram matrix or a bare rank would erase data.

## Reproduction

From the repository root, run:

```text
python -X utf8 -B 05-knowledge/results/odd_descent_local_relaxations_20261030/verify.py --output 05-knowledge/results/odd_descent_local_relaxations_20261030/certificate.json
python -O -X utf8 -B 05-knowledge/results/odd_descent_local_relaxations_20261030/verify.py --output 05-knowledge/results/odd_descent_local_relaxations_20261030/certificate.optimized.json
python -X utf8 -B 05-knowledge/results/odd_descent_local_relaxations_20261030/seal.py
```

The manifest hashes exact LF packet bytes and exact Git blobs of inherited
sources. It records the base revision and the separately read incoming
Session118 frontier. No timing claim, formal axiom audit, or general
open-conjecture consequence is included.
