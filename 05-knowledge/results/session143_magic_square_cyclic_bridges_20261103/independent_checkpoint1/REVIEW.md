# Session 143 checkpoint 1 — independent proof review

**Disposition: PASS for the scoped tournament and cyclic-operator results.**
Reviewed against source revision `8286d710498a3b4e9b9f1ac683c6eb77716f3f4b`.
No repository edits or publication were performed by this audit lane.

## Independently established

- Direct ordered-arc orbit enumeration gives 12 tournament classes, eight
  self-converse classes, ten classes modulo converse, and nine adjacency
  characteristic polynomials. Exact Newton identities agree with independent
  directed-cycle counts. Spectrum plus score sequence separates these 12
  classes; the spectrum alone does not.
- The unique regular class has polynomial
  `(x-2)(x^4+2x^3+4x^2+3x+1)`. Its automorphism group has order five and its
  five anti-automorphisms are involutions. The four labeled circulant
  orientations have precisely the affine order-20 setwise stabilizer, with
  translation kernel of order five on the family. An affine family action is
  not a 5:4 automorphism action on a single tournament.
- Reviewed Burnside, trace, regular-uniqueness, cyclotomic inverse, and
  normalizer arguments in `tournaments/RESEARCH.md`; no gap found at the
  stated order-five scope. These are elementary proofs or finite-exact
  classifications, not extrapolations to larger orders.
- Exact arithmetic in Q[z]/Phi5 independently verifies all four roots of
  `q=x^4+2x^3+4x^2+3x+1` and Fourier multiplicities for the specified
  group-algebra operator `rho(g)+rho(g)^2`. Imported 5A/5B character values
  give `(891,870,870,870,870)` and `(871,875,875,875)` respectively. Their
  traces 42 and -8 exclude loopless adjacency operators in any basis.
- The field intersection argument is valid: coprime conductors give degree
  64 for Q(zeta85), so Q(zeta5) and Q(zeta17) intersect only in Q. The quartic
  roots generate Q(zeta5). The function field Q(X0(13))=Q(t) is not the
  number field Q(zeta17); specialization is not a whole-function-field map.

## Independent primary-source checks

- [McKay catalogue](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html):
  the listed order-five counts are 12 total, eight self-converse, one regular.
- [ATLAS Baby Monster](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B):
  order, two order-five classes, their centralizers, and a maximal subgroup
  `5:4 x HS:2` are present. This is cited subgroup occurrence, not execution
  of its straight-line program by this session.
- [Breuer–Magaard–Wilson protocol](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/htm/ctblbm.htm):
  sections 6–7 supply rational power-map data and the character in ATLAS
  class order. Local snapshot transcription verifies the 184-position
  alignment and values 21, -4. The underlying character-table reconstruction
  is the authors' result; this audit did not rerun GAP or construct matrices.
- [Steel thesis](https://magma.maths.usyd.edu.au/~allan/reps/AllanSteelPhD.pdf),
  section 5.11, printed page 119: minimal faithful ordinary degree 4371,
  realizable over Q. Since B is simple, this excludes nontrivial
  characteristic-zero 24-dimensional representations. It says nothing here
  about modular or nonlinear actions.
- [Daniels et al.](https://math.mit.edu/~drew/mcom3213.pdf), proposition 6.2
  proof, printed page 446: X0(13) has genus zero and the stated rational
  j-parametrization. Therefore its weight-two cuspidal space is zero; this
  is not a denial of all Hecke correspondences or other weights.
- [Conway–Elkies–Martin](https://arxiv.org/html/math/0508630v1), introduction
  and section 2.2: the puzzle has 13 points, 12 counters and a hole, while
  M13 denotes a pseudogroup rather than a five-point Mathieu group.

## Audit correction and controls

The original sporadic checker used bare assertions. This audit requested
explicit checks; v2 now contains no `assert` nodes and rejects deliberately
wrong eigenvalue, character and transcription inputs under both Python modes.
Its GAP script remains **NOT EXECUTED**. The independent checker here uses
no producer imports, has explicit exception checks, and rejects a forced
wrong spectrum count under ordinary Python and `-O`.

Source checks pin exact named repository blobs; they are not a re-audit of
all theorems or historical experiments in the source checkout. No octonion,
magic-square, universal-conjecture, or priority claim is made in this first
checkpoint audit.
