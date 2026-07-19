# LRC(14): the decorated pair-obligation complex

## The common object behind the newest routes

Several apparently different LRC(14) arguments now factor through the same
kind of carrier.  A bare pair of runners is not the vertex.  The useful vertex
is a **pair obligation** together with the data needed to interpret it:

```text
(labelled pair, signed difference or beat, phase cell,
 denominator/gcd, endpoint owner, metric weight).
```

Relations among these vertices form small proof-bearing 2-cells.

1. In THM-1210, continuum four-comb BAD forces all six pair differences of
   four phase points into the three-band set
   `A=[1/7,2/7] union [3/7,4/7] union [5/7,6/7]`.  Deleting to a non-AP
   triple leaves the additive circuit `(p,q,p+q)`.  Its six torus alcoves,
   not the runner order, carry the measure bound.
2. In THM-1156/1178/1179, an exact two-tooth seam has a `chi_7` colour
   switch, but an open cover cannot use a zero seam without a third owner.
   The proof cell is therefore a labelled pair plus a third-support
   hyperedge and an endpoint quantum.
3. In THM-1192, a sum/difference beat of a defining pair produces a safe
   rational point for both runners.  The remaining four or five combs must
   cover that point.  The proof cell is
   `defining pair -> beat numerator -> complementary coverers`.
4. THM-1191 (the thirteen-adic pair-floor refutation) shows why scalarizing
   these cells is unsafe: the exponent stalk preserves an alternating sign
   and an exponentially decaying weight.  Its tournament sees the sign but
   not the magnitude.
5. THM-1215 identifies one terminal beat stalk.  When a defining difference
   is `q=14d` and all six fast speeds have gcd `d` with `q`, every reduced
   danger mask is the same singleton `{0}` in `Z/14Z`.  A slow-gap numerator
   block has at least six consecutive points as soon as `q>=7a`, so it must
   contain a point outside that singleton.  Here the decorated complex
   genuinely collapses, but only because the phase block is retained and
   the masks coincide exactly.

The shared lesson is that pair information becomes useful only after its
compatibility cells and sidecars are retained.

## Static and dynamic quotients

The continuum quotient is static.  At one phase `u`, BAD supplies a labelled
`K4` of three-band obligations.  One additive triangle is enough for an upper
bound.  After common-gcd reduction and the shear

```text
(x,y) -> (x,z=x+y),       z={(p+q)u},
```

the section is empty off `A` and has at most two intervals on `A`.  Hence

```text
J(p,q) <= 3/49 + 6/[7(p+q)].
```

The analytic tail begins at `p+q=26`; exactly 99 reduced pairs remain.  The
unique equality pair is `(1,2)`.  Retaining all three missing `K4` edges then
forces the fourth frequency to fill the missing point of a four-term AP.

The beat-puncture quotient is dynamic.  It does not ask for all pair
obligations at one generic phase.  Each pair chooses its own rational beat
point `p/(u+v)` or `p/|u-v|`, and the cover must transport to that point.
Replacing the truncated numerator block by its average kill density destroys
the very phase which the theorem exposes.  Thus the continuum shear and the
beat quotient are complementary:

```text
continuum triangle: one phase, preserve additive compatibility;
beat puncture: many rational phases, preserve complementary-cover transport.
```

Boxeph S120's straddle formulation of the earlier pinch lemma
(THM-401/THM-1170) makes the sum-beat half structural, not merely
opportunistic: every maximum of `min_v ||vt||` has two active
runners straddling adjacent integer sides and occurs at
`t=(a_i+a_j)/(v_i+v_j)`.  Hence the global witness search is supported on at
most `binom(n,2)` sum-beat moduli.  For the open twelve-speed rigidity branch,
the missing statement is now exactly that a non-AP interior packing cannot
keep every one of those straddles at or below `1/13`.  The active pair is a
decorated obligation; the other ten phases are its feasibility sidecar.
This locates the hard object without proving Tao's `n=12` uniqueness.

THM-1215 is the first all-scale example where this transport can be completed
without averaging.  Its block sidecar supplies six consecutive residues and
its gcd sidecar identifies every fast mask with `{0}`.  The resulting witness
is literal, not merely a deficit in a union bound.  This suggests treating
mask coincidence as a terminal reduction rule, alongside additive-triangle
and seam-debt discharge.

## What tournament fingerprints do and do not retain

The runner-order tournament is transitive in the clustered and continuum
arguments.  Its score histogram, SCCs, cycles, and unique Hamiltonian path
contain no metric obstruction.  A sign tournament on exponent levels can be
strongly connected and diagnostically useful, but THM-1191 shows that it
still loses the weights needed for the pair sum.

THM-1200 adds an important typing guardrail.  The seven in the analytic
danger/safe partition is `14/2`, a parity feature of the LRC threshold; the
seven in the Paley/Fano tournament is the prime field size.  Their numerical
coincidence does not identify their circuits.  A Fano triple may organize
which pair obligations to sum, but its validity still has to be witnessed by
the `chi_7`, gcd, phase, and endpoint data of the analytic wall.  Conversely,
an analytic seven-wall inequality does not inherit quadratic-residue
symmetry merely because both objects have seven labels.

The meaningful binary switch is therefore local to a decorated obligation:

```text
three-band membership, strict seam sign, chi_7 colour,
or survival of a beat numerator.
```

Completing those switches to a tournament is optional telemetry.  The proof
object is the resulting weighted 2-complex.  A contraction is legitimate
only if it reconstructs the phase cell, complement relation, and endpoint
owner used by every incident 2-cell.

Klein's `+/-`-fibration gives a precise type for this distinction.  A Cayley
tournament chooses a section of
`(Z/qZ)^* -> (Z/qZ)^*/{+/-1}`, while the LRC divisor deck lives in the base.
The symmetric three-band predicate in THM-1210 also factors through that
base.  Seam direction, owner, and beat transport are not a global section;
they are local lift data attached only where a 2-cell needs them.  This
explains why the decorated complex can borrow tournament language without
becoming a tournament, and why self-complementary/spine arguments do not
transfer: negation acts freely on the LRC unit classes for `q>=3`.

Candidate vertex sets explicitly challenged in this synthesis are runners,
gaps, fixed circle sections, section boundaries, wall events, residues,
cover arcs, Fourier modes, matroid circuits, valuation levels, beat points,
and proof obligations.  Pair obligations are currently the smallest common
choice, but only with their sidecars.

## A recursive proof program

The combined routes suggest a more disciplined recursion for the remaining
global problem.

1. Extract a minimal obstruction as a decorated obligation complex, not as a
   speed subset alone.
2. Delete or contract a vertex only when every incident additive, seam, and
   complement-cover cell has an explicit reconstruction map.
3. Discharge rank-two additive cells by the THM-1210 triangle bound.
4. Discharge zero-seam cells by third-support/owner debt and retain their gcd
   quantum instead of a uniform scalar floor.
5. Discharge rational beat cells by a positive puncture deficit; if every
   deficit vanishes, retain the full truncated-block phase as the residual
   state rather than averaging it away.
6. If all reduced masks on a beat cell coincide, test the retained numerator
   block directly against that common mask; THM-1215 is the singleton-mask
   model.
7. Route valuation-stalk residuals before applying any finite-height or
   phase-free pair estimate.

This reframes the next hard theorem.  It is not merely a stronger pair-overlap
constant.  It is a **phase-preserving circuit-elimination theorem**: every
minimal decorated obstruction should expose either an additive triangle, a
positive seam/owner debt, a beat puncture, or a named valuation stalk.  A
successful statement would explain why the current exact finite banks have
no survivors without pretending that those banks are uniform proofs.

## Honest boundary

THM-1210 closes the sharp continuum four-comb ceiling and its equality locus.
THM-1214 closes the complete clustered five-killer stratum, and THM-1169 plus
THM-1212 close the clustered six-killer stratum.  The decorated-complex
viewpoint does not by itself prove global LRC(14), CrownCollapse14, or the
uniform non-AP twelve-speed equality classification.  The historical uniform
`q<=25` good-period claim is false (MISTAKE-143); the new number 25 in
THM-1210 is the unrelated and rigorously derived finite condition
`p+q<=25` for additive-triangle gaps.
