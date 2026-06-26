# Roth-Minkowski Diophantine Lattice Fence

codex-2026-06-26-S226

User request: merge Roth's theorem and Minkowski's theorem into the LRC14
proof search.

The useful interpretation is Roth's theorem in Diophantine approximation,
paired with Minkowski's geometry-of-numbers theorem.  The additive-combinatorics
Roth theorem may become relevant for density/AP packets, but the active bridge
to support-six, Farey, height, and low-wall work is the Diophantine one.

Minkowski says: once a lattice, covolume, and symmetric convex body are named,
volume pressure can force a short nonzero lattice vector.  In this repo that is
not a theorem by itself, because the quotient may have forgotten which relation
lattice, which convex body, and which packet wall supplied the volume.

Roth says: once an algebraic irrational target and height scale are named, the
inequality `|alpha - p/q| < q^(-2-epsilon)` has only finitely many rational
solutions.  In this repo that is also not a theorem by itself, because the
finite exceptional approximants and low-height walls are exactly the kind of
payload controlled forgetting tends to erase.

The combined carrier is a fence:

```text
finite low-height wall deletion
  -> Minkowski relation-lattice tail
  -> Roth algebraic-near-miss fence
  -> named packet exit or residual debt
```

This clarifies the older "execute Minkowski count" language around the
support-six/HYP-2612-HYP-2614 frontier.  Before the count can be promoted, a
packet row needs a sidecar with:

```text
relation_lattice
covolume
successive_minima_profile
convex_body_id
algebraic_target
field_degree
height_bound
approximation_exponent
epsilon_margin
exceptional_approximants
low_height_wall_class
deleted_anti_cosets
residue_signed_tail
diophantine_exit
```

LRC consequence: volume gates and approximation exponents should be treated
like observer-cut payloads.  They are legal to forget only when reconstructed,
annihilated, descended to a packet family, stopped at AP/GW, or named as
residual debt.  Raw volume and raw exponent scalars are scouts.

Next work: add these fields to a small HYP-2963/HYP-2614 support-six sample,
then prove the finite low-height deletion before invoking Minkowski and the
finite exceptional-approximant ledger before invoking Roth.
