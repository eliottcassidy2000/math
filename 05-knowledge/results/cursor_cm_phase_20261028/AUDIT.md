# Proof and computation audit

**Same-agent alternate-path audit.** No independent referee or formal Lean
verification is claimed. Initial source59cb22d0e93b; first checkpoint547025d2ff
integrated concurrent work through2f84e34bd6. Mathematical statements use
local IDs plus their filenames; no scarce global theorem ID was reserved.

## Executed universes and controls

| Check | Universe | Outcome and limit |
|---|---|---|
| Direct cursor | All path presentations through8 | No collision; at8 all2097152 inputs |
| Prime-order cursor cover | Every source with nontrivial endpoint group through9; rigid sources handled by proof | No collision; n9 bank386054, dual nontrivial580; all268435456 possible presentations covered logically |
| Degree endpoint decoder | All tournament classes through5, all directed spanning forests with k>=2 modulo source action, all faithful component orders | 663 forest orbits,5346 ordered-cover orbits decoded; no endpoint merger through5 |
| Phase-loss families | Class-only m=3,5,7; length-aware m=3,5,7; equal-length nine-vertex control | Partner endpoint image exceeds G; source endpoint image is trivial in every recorded case |
| Hostile81/27 | Two literal twelve-vertex tournaments; every singleton port decoded | Same67 labelled cyclic triples and H=462915; automorphism orders81/27 remain distinct |
| Cheap short C3 context | Same small universe, one unmarked context | 26 forest failures:2 at3,2 at4,22 at5; not universal response collisions |
| Eta/CM arithmetic | Exact series through q^32; eight discriminant pairs; five Paley character sums; quoted unit traces | 93 checks; exact data separated from certified/numerical checks |
| Quadratic trace carriers | All states for a=1,2,5,12; m=1..30; seven power menus, plus explicit moduli | 763611 checks; measured kernels/images equal Smith predictions |

The endpoint checker makes265920 explicit checks, using RuntimeError gates
that remain active with Python optimization. It relabels each constructed
context by a varying cyclic shift and recovers the old relation and every
ordered pair intrinsically. Its proof establishes invariance under every
permutation; the finite tests do not enumerate all permutations of each
large context. Strongness is explicitly checked for every hostile source.

`reproduce.py` creates a fresh temporary repository-shaped mirror, copies
only the required scripts and the two inherited helpers, and replays
endpoint_kernel, cm_series, unit_clocks, run_probe and cursor_symmetry under
both normal Python and `-O`. All ten outputs agree exactly with the saved
JSON receipts. Python3.13.14 was used. [audit_receipt.json](audit_receipt.json)
pins the LF-normalized source and receipt hashes.

The full direct order-eight census was executed once in the producer.
The default replay uses the complete prime-order cover through nine rather
than repeating that much larger rigid-input scan. `--direct-eight` adds
normal and optimized replays of the direct census if desired.

## Alternative paths and their dependence

- Cursor coverage has two paths at7/8: literal bit enumeration and Cauchy's
  prime-order cover. Both share the inherited exact GI oracle. Neither is
  called an independent GI implementation.
- The endpoint theorem's output decoder does not use GI on its large
  constructed tournaments. It reconstructs marked roles by unique degrees
  and arcs, then compares the recovered relation literally. Small action
  sizes use the inherited oracle separately.
- Eta coefficients are computed by the odd product and independently by
  product_(j>=1)(1+q^j)^(-24). Formal inversion then computes the Fricke trace.
  Algebraic unit traces and class numbers use integer recurrences and reduced
  quadratic forms. Five Paley ranking sums are a second arithmetic observable.
- Exponential bounds use rational Machin arctangent intervals and90-digit
  outward rounding. `sqrt` and `exp` round correctly; one adjacent Decimal
  is taken on each side. Radical/product comparisons are numerical checks,
  explicitly not proofs of the cited class-invariant evaluations.
- The quadratic-clock theorem has a symbolic Smith proof and exhaustive
  direct enumeration of image/kernel sizes. Mixed-clock closure is proved
  by the ring Z+gZu; the finite census alone would not establish it.

## Proof-critical boundaries checked

E2 fixes an aligned pair and compares orderings of its supplied intact runs.
It observes ordinary source/partner classes; it does not retain the final
Hamilton path or a blue flag. Endpoint groups are counted through their
PERMUTATION IMAGES: a hidden interior kernel must not inflate the count.
Length colours are optional and change the smallest hostile from7 to8.
For k=1 there is no order ambiguity; the degree construction is for k>=2.

E6 permits general child tournaments on both native sheets. In the eight-
vertex length-aware example the partner is not strong, so its pre-complement
repair is not available if a different menu permits strong children only.
The auxiliary marker children in E2 are always strong. Whether the two-sheet
intersection always equals the faithful simultaneous group remains OPEN.

Q1 works on all m² states of the specified quadratic ORDER. The37 unit's real
order has discriminant148, not37; all-state fibers are not norm-one fibers.
The real units, imaginary CM fields and Fricke function field are distinct
sources. Arbitrary selected power clocks, inverse clocks and conjugation
are covered; arbitrary nonlinear arithmetic observations are not.

Inherited gap topology is static and target-sensitive: all-runner failure
gives Ind(C14)~S4; a designated runner gives a contractible complex. The
star rank obstruction concerns exact real linear energy interfaces only.
It rules out a rank bound from homotopy type, not every nonlinear encoding.
Lucas sign-evenness, odd tournament symmetry, temporal reflection and
Frobenius are distinct actions. No phase is declared free from group order.

## Counter and attribution repairs

The first checkpoint's cursor table called20/160/2256 distinct candidate
presentations. They are tests with multiplicity. For rows
(62,100,88,82,98,76,1), ten tested permutations yield only two different
competitors. The label is repaired in CURSOR.md and the mechanism logged
in [MISTAKES.md](../../../01-canon/MISTAKES.md). Coverage and all receipts
are unchanged; no performance claim uses these counters as distinct states.

The class-number table uses primitive discriminants of the displayed CM
points; for the selected7 mod8 values the conductor-two counts coincide,
but the orders are not identified. 'Quadratic surds imply a near integer'
is explicitly refuted at130. No raw minimal-polynomial degree is inferred
from the excerpt's cubic/quintic/septic shorthand.

The Fricke series, Monster character constituents, cycle independence
homotopy types and cited CM values are credited as classical or inherited.
Four repository search waves are reproducible as a historical scan at the
stated source; substring hits are not proof dependencies. No universal
modular hub, major-open-problem solution or new all-order enumeration
algorithm is asserted.
