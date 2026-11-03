# Pairing symmetries, four binary digits, and growth recurrences

**PROVED scoped results; INHERITED PROVED count formulas; FINITE-EXACT
controls. The unrestricted Hamilton-path count spectrum remains OPEN.**
Research date: 2026-11-02 America/Denver; directory date is UTC.

The latest request joins three different recursive objects: coordinate
regions, tournament orientations, and isomorphism classes. The useful
unifying principle is to specify **which action and which observable a
recurrence preserves**. A power of two in a count permits an abstract
binary action; it does not identify that action with converse, matching
edge flips, or an operation preserving Hamilton paths.

This packet adds:

- [Compatible pairings](PAIRINGS.md): the maximum rank of a free elementary
  binary action commuting with converse is exactly
  min(v2(a_n),v2(c_n)), where c_n counts self-converse classes. For odd
  n>=3 this is (n-1)/2; for n=2 mod4 it is n/2-1. These are abstract
  class actions with an explicit construction, not native metagraph moves.
- [Four-bit formulas](FOUR_BITS.md): exact formulas modulo16 for
  a_n/2^floor(n/2), using coprime pairs and triples and explicit small-order
  corrections. These extend the inherited modulo4 formulas; the
  fixed-precision algorithm itself is inherited.
- A weighted CRT evaluator for the even-order four-bit formula, replacing
  the generic cubic partition scan by a linear gcd scan or a sum over
  squarefree divisors after factoring. Factoring cost is included; this
  accelerates the specified low digits, not the full integer count.
- [Recovery and limitations](RECOVERY.md): typed +2 versus doubling and
  seven-term rules, the toothpick boundary, and the joint strong-component
  recurrence retaining Hamilton count. The latter is recovered prior work.
- A corrected **minimum four-frame cover** of the twelve order-five classes.
  The old June reflection said five. All 256 frames are checked; there
  are 136 full frames, 15 distinct full class blocks, no disjoint three-block
  tiling, and a displayed four-block cover. The failed product decomposition
  survives the correction.
- A minimum fixed-point formula for any H-preserving involution commuting
  with converse. At orders3,4,5,6 the unavoidable counts are2,2,4,6.
  The inherited transitive-class obstruction already rules out a free
  H-preserving pairing at every order.

## Inheritance and concept board

Start revision: e30c89e9fcdc. The closest proved mechanisms are
[THM-305, odd-order valuation](../../../01-canon/theorems/THM-305-two-adic-valuation-tournament-count.md),
[October25 precision P1–P11](../port_residual_integral_20261025/precision/README.md),
[October21 converse channels](../converse_tournament_count_channels_20261021/session22-converse-channels/README.md),
and [THM-553, tile addresses](../../../01-canon/theorems/THM-553-tile-address-recurrence.md).
The canonical hostile remains the actual 81/27 pair. The corrected near
miss is a count factorization interpreted as a native class action.
The least-used sidecars are the self-converse fixed locus and the odd
part of each centralizer denominator.

Anchor: compatibility of class-count symmetries with the H spectrum.
Niche: next binary digits beyond the existing totient checksum.
Wildcard: exact matching-frame covers versus toothpick reflection.

| Live concept | Retained data | Consequence / boundary |
|---|---|---|
| +2 and doubling | order, H, charge, actual source | different recurrences cannot share an untyped scalar |
| Seven-term geometry | named coordinate sets and intersections | exponentiates before quotienting; does not average multiplicatively |
| Converse | ordinary class and fixed locus | exact commuting-action rank |
| H spectrum | strong atoms and ordered SCC word | multiplicative support; source-free growth still fails |
| Binary precision | cycle lengths, gcds, odd denominators | exact low digits; not full-count reconstruction |
| Toothpicks | physical exposed ports and birth time | real dyadic recurrence; no literal converse bijection |

The incoming commits 090fcaa1e7 through a28303347b were read before
close-out. [Session149](../session149_tournament_spectrum_20261103/README.md)
already recovered the parity proof, folding and a sorted complete pairing.
The [joint parity packet](../tournament_joint_parity_20261102/README.md)
adds an alternating-trail pairing, joint spectrum recurrences and a toothpick
port audit. Credit belongs there for those constructions. Our contribution
is the *commuting rank and observer-specific defect*, the next low digits,
and the exact four-frame correction. No historical or literature priority
is asserted for the elementary finite-action lemma.

The later incoming947b6c3852 [third-bit packet](../fibonacci_parity_precision_20261103/PRECISION.md)
was published before our first push and receives credit for the explicit
modulo8 formula and order40 hostile. Our fourth-bit formula supplies its
requested next layer, including the shifted totient term. Its factorized
modulo8 simplification agrees with our weighted CRT expression after reduction.

The existing “Controlled forgetting and unlabeled quotients require a
sidecar” method card is used. No new meta-pattern is promoted: the evidence
refines its fixed-locus requirement.

## Verification and continuation

Run from the repository root:

    python -X utf8 05-knowledge/results/paired_recursions_20261103/verify.py
    python -O -X utf8 05-knowledge/results/paired_recursions_20261103/verify.py --output receipt-optimized.json

The two receipts agree. Full rational Burnside sums through40 independently
check the new modular formulas; the inherited bounded-cycle evaluator
checks them through200. A separate signed pair-orbit implementation checks
self-converse counts through12. All labelled orientations are partitioned
into classes through6; explicit free actions, H-preserving pairings and both
SCC recurrences are checked on those actual classes. See [audit](AUDIT.md)
and [handoff](HANDOFF.md).

None of these results proves that7 and21 are the only missing positive odd
H values, a native free binary action of the arithmetic rank, or a compact
complete recurrence for all individual metagraph adjacencies.

