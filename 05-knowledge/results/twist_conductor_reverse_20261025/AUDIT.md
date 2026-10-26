# Internal proof and exact-control audit

2026-10-25. **PASS within the stated scopes.** This is a single-session
internal review with alternative exact computations, not an independently
assigned proof audit. Incoming Session64 has its own separately attributed
audits; those do not automatically audit this packet's new algebraic claims.

## Load-bearing mathematical checks

* I1's conductor uses the same fraction field. Every target character is
  distinct after multiplication by one homogeneous component, excluding
  hidden cancellation. The y=0 term already forces membership in B. The
  valuation version explicitly needs reflexive target ideals.
* I2 is limited to **two complementary binary packets** with disjoint
  height-one radicand supports. It does not absorb the earlier three-packet
  counterexample. Normalization is a separate hypothesis.
* I3's normality comes from the sign invariant ring and remains true after
  extending constants to Q(sqrt(2)). The conductor is nonprincipal even
  after that extension. The maximal ideal hostile has height two, so all
  height-one tests genuinely miss it.
* I4's normalized DVR has the integral ratio w with w^2=2. It is the same
  field as the product order. The original radicands share pi; this is
  outside the disjoint-support theorem. The repaired basis test retains
  the ratio, its relation, and both original ports.
* F1--F2 retain the specified twist and component. The U,D restriction is
  recomputed, not guessed from degree. The rational parametrization has an
  explicit inverse on its stated open chart. Omitted points are handled in
  F3's direct root equations, not ignored by parametrization.
* F3's rational exclusion uses the written infinite descent; a rational scan
  is not its proof. The E=0 chart is separately excluded. Point counts refer
  to the full original boundary versus its chosen component explicitly.
* F4 concerns the marked deck group, not every automorphism of the curve.
  K5 and K3,3 are excluded only as subgraphs of the specific deck graph;
  no minor or chromatic implication is claimed. The star/triangle pair
  verifies the necessity of higher-arity dependencies.
* F5's elliptic factor has a character outside the three-dimensional Y
  character space. Field degrees, reconstruction formulas, common t and
  normalization prove the fibre-product claim. The degree-zero half-divisor
  is nonprincipal geometrically because the double cover is connected.
  An independent point on the elliptic curve is not a compatible point.
* F6 recovers the existing stronger five-factor theorem, extends its exact
  Q-twists, and checks the h-invariant differential characters. The Prym
  has dimension two; it is not identified with the elliptic auxiliary curve.
* Port-code counts divide by GL(2,2) only after verifying all assignments
  span. Largest dimension is distinguished from inclusion-maximality.
  The b=6 and b=7 controls prevent extending the special Petersen rule.
  The nine-extension codec is exhaustively inverted through b=9. Six cases
  change one old quotient label, necessarily; their reversibility is not
  falsely promoted to preserving frozen input palettes or projective lines.
* D1 is a conditional probability law using the incoming harmonic function;
  it is not a new deterministic orbit map. All infinite-alphabet tails are
  retained. The uniform-integrability obstruction follows from absolute
  continuity of integrals, with Holder's quantitative bound and a sharp
  cylinder-density example. The conclusion concerns this proposed observer,
  not the truth or falsity of Collatz.

## Exact finite universes

| Checker | Universe / independently checked predicate |
|---|---|
| `check_geometry.py` | 1,377 rational parameter substitutions; exact polynomial gcds; all 15 nonidentity deck elements and all 45 pairs of free involutions; original U/D characters; explicit elliptic squareclass |
| `check_ideals.py` | monomial exponent box 0..12 in the sign-invariant ring; both ideals P and m; normalized DVR exponent 0..12 and both w parities |
| `check_port_codes.py` | all 3^b column assignments for b=5,7,9; complete GL(2,2) orbit counts; 921 recursive/projective round trips; b=5/6/7 disjoint-port hostiles |
| `check_golden_tilt.py` | independent Q and Q(sqrt(5)) horizon recursions, horizons 1..32 and all initial counters <= horizon: 528 comparisons |
| `check_repairs.py` | complete blue loci of six failed affine cells: exactly 20 literal sources; directedness, blackness, XOR, exact labelled endpoint and converse transport |
| `check_reverse.py` | N=0..60, K=1..8: 488 independent odd-partition comparisons; 73,557 structural edges; previous absolute engine for N=0..40 and N=50; order50 and invalid-type hostiles |
| `check_incoming.py` | exact shared-valuation hostile; all projective F3 conic and central-line points; ideal containment direction |

All load-bearing gates raise exceptions and remain active under `python -O`.
Every normal and optimized receipt pair matches byte for byte. Optional SAT
searches are witness producers only. A timeout and a failed small word bank
are recorded as unknown/finite failures, never as nonconnectivity.

## Corrections and provenance

The initial golden test represented q^6 using the coordinates `(9,-4)` of
the basis `(1,sqrt(5))` inside arithmetic using `(1,q)`. The exact check
rejected it. The repaired coordinates are `(5,-8)`, with both bases labelled
in the receipt. The mathematical identity q^6=9-4sqrt(5) was unchanged.

During internal review I2's two-packet hypothesis was made explicit to avoid
an unintended generalization to multiple packets. The earlier proved
multi-packet hostile remains live. No previously proved theorem is retracted.

The subsequent incoming THM-4772/4773 review did require corrections to
published implications: reversed ideal containment, an invalid unreduced
single-root shortcut, and the infinite-field requirement for avoiding
finitely many conic points. These are recorded in INCOMING.md and the
mistakes ledger, with exact witnesses in `check_incoming.py`; the repaired
source checksums are refreshed. No broader audit of those papers is claimed.

The source search recovered THM-261 (A4 roots), THM-4701 (Clebsch/Petersen
finite-circle graph), THM-2685 (divisor parity and class group), THM-3397
(torsor versus effective valuation), and Session56's LCA report and ten-card
niche atlas. The incoming Session64 proofs and audits were read as listed
in INCOMING.md. Classical Humbert claims are attributed to the primary
Hidalgo paper, particularly Theorem 1 and section6; the present fixed-twist
equations and their rational compatibility are shown directly.
The late Humbert search completed with four hits; its curve-factor sidecar
was read in full and integrated in F6. This matters because the five-factor
isogeny was already stronger than merely finding an elliptic quotient.

No universal Collatz, Antihydra halting, unrestricted black-core connectivity,
Hadwiger--Nelson, planar Jacobian or LRC(14) conclusion is claimed.
