---
id: THM-2643
title: "Degree-five/six Keller stabilizer gate and regular block-quotient census"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a transitive
  action Omega=G/H and N the normal closure of H,
  N is proper exactly when Omega has a nontrivial quotient G-set on which
  the induced quotient group acts regularly; the canonical such quotient is
  G/N.  THM-2633 therefore excludes precisely one of the five transitive
  quintic groups (regular C5) and nine of the sixteen transitive sextic
  groups.  In degree six every two-block C2 quotient and every three-block
  C3 quotient is fatal.  An S3 block action is not fatal by itself: the
  three nonregular imprimitive actions 6T7, 6T8, and 6T11 survive because
  their reflection stabilizer normally generates S3, while 6T3 is killed by
  a separate regular C2 quotient.  This is the exact
  C2/C3 free-factor content: an isolated cyclic tower is forbidden, not an
  S3 block action solely as such.  The census is action-sensitive; natural quartic
  A4 survives although its sextic C2-coset action is excluded.  No degree
  bound, A4/S4 quartic exclusion, JC(2), general JC, or DC(2) follows.
source: root-2026-07-28-keller-stabilizer-census
depends_on:
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
script: 04-computation/keller_degree_five_six_stabilizer_census_thm2643.py
output: 05-knowledge/results/keller_degree_five_six_stabilizer_census_thm2643.out
script_sha256: de99329cf123aac62376de63a33eab0d444f6119c6d478c435b540ba2c2177ae
output_sha256: f2ca355ef681eb4fc5d6a0a5d3af0fae46e959834e41166615536ae6ffd3ab6a
hash_basis: LF-normalized bytes
---

# THM-2643 -- regular sheet quotients are forbidden Keller shadows

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  THM-2633 says
that a point stabilizer in the geometric sheet action
of a polynomial Keller map must normally generate the full monodromy group.
The present theorem identifies the exact permutation-theoretic failure and
classifies every transitive action through degree six.

## 1. The stabilizer closure is a regular-quotient detector

Let a finite group `G` act transitively on

```text
Omega=G/H,                         N=normal_closure_G(H).  (1)
```

Since `H<=N`, there is a canonical `G`-map

```text
pi:G/H -> G/N.                                             (2)
```

The induced faithful group `G/N` acts regularly on the target.  Therefore

```text
N<G  ==>  Omega has a nontrivial regular quotient.         (3)
```

Conversely, suppose `Omega` maps equivariantly onto a set `Y` such that the
induced quotient `G/K` acts regularly and `|Y|>1`.  The stabilizer in `G` of
the image of `H` is exactly `K`; since `H` fixes that image,

```text
H<=K,                       N<=K<G.                        (4)
```

Thus

```text
N<G  iff  Omega has a nontrivial regular quotient G-set.   (5)
```

Here “nontrivial” means that the target has more than one point.  The target
may equal `Omega` when the original action is regular; in block language the
singleton partition must be retained.  Excluding singleton partitions would
incorrectly miss regular `C_5`, `C_6`, and `S_3`.

There is another immediate boundary.  If `g` fixes the sheet `xH`, then

```text
g in xHx^(-1) <= N.                                       (6)
```

Hence every element of `G\N` is a derangement.  THM-2633 makes (5) a
necessary Keller gate:

```text
Keller monodromy  ==>  N=G  ==>  no nontrivial regular sheet quotient. (7)
```

## 2. Prime degree and the complete quintic census

For a faithful transitive action of prime degree `p`, the `N`-orbits are
blocks.  If `N` is transitive, then `G=NH=N`.  If `N<G`, prime degree forces
the `N`-orbits to be singletons; faithfulness gives `N=1`, hence `H=1` and
`G=C_p`.  Conversely regular `C_p` has trivial point stabilizer.  Therefore

```text
in faithful prime degree, N<G iff G=C_p acts regularly.    (8)
```

For degree five this gives the whole table:

| label | group | `|H|` | `|N|` | `G/N` | Keller gate |
|---|---:|---:|---:|---:|---|
| `5T1` | `C5` | 1 | 1 | `C5` | **excluded** |
| `5T2` | `D10` | 2 | 10 | 1 | passes |
| `5T3` | `C5:C4` | 4 | 20 | 1 | passes |
| `5T4` | `A5` | 12 | 60 | 1 | passes |
| `5T5` | `S5` | 24 | 120 | 1 | passes |

“Passes” means only that THM-2633 does not exclude the action.

## 3. The complete sextic census

For the standard transitive-group labels, direct closure of every point
stabilizer gives:

| label | group/action | `|H|` | `|N|` | `G/N` | Keller gate |
|---|---|---:|---:|---:|---|
| `6T1` | regular `C6` | 1 | 1 | `C6` | **excluded** |
| `6T2` | regular `S3` | 1 | 1 | `S3` | **excluded** |
| `6T3` | `D12`, hexagon | 2 | 6 | `C2` | **excluded** |
| `6T4` | `A4` on `C2` cosets | 2 | 4 | `C3` | **excluded** |
| `6T5` | `S3 x C3` | 3 | 9 | `C2` | **excluded** |
| `6T6` | `A4 x C2` | 4 | 8 | `C3` | **excluded** |
| `6T7` | `S4`, edge action | 4 | 24 | 1 | passes |
| `6T8` | `S4` on `C4` cosets | 4 | 24 | 1 | passes |
| `6T9` | `S3 x S3` | 6 | 18 | `C2` | **excluded** |
| `6T10` | `C3^2:C4` | 6 | 18 | `C2` | **excluded** |
| `6T11` | `S4 x C2` | 8 | 48 | 1 | passes |
| `6T12` | `A5`, degree six | 10 | 60 | 1 | passes |
| `6T13` | `S3 wr C2` | 12 | 36 | `C2` | **excluded** |
| `6T14` | `S5=PGL2(5)`, exotic | 20 | 120 | 1 | passes |
| `6T15` | `A6` | 60 | 360 | 1 | passes |
| `6T16` | `S6` | 120 | 720 | 1 | passes |

Thus nine of sixteen sextic actions are impossible for a polynomial Keller
map in every dimension.  In every excluded row the full complement `G\N`
is fixed-point-free.  Each row also has a prime-cyclic derangement character:
use `C2` or `C3` on the displayed cyclic quotients, a prime quotient of
`C6`, and sign on the regular `S3`.  Hence degree at most six does not yet
need the strictly nonabelian strength of (5), although the normal-generation
form organizes the census losslessly.

## 4. Why two and three are the exact sextic fault line

There are ten set partitions of six sheets into two triples and fifteen
partitions into three pairs.  The exact invariant-block census has the
following consequence:

```text
two blocks:  induced C2 is regular and is always fatal;
three blocks: induced C3 is regular and fatal;
three blocks: induced S3 is not regular and can survive.   (9)
```

The imprimitive survivors `6T7`, `6T8`, and `6T11` have an induced `S3`
action on three blocks.  A block stabilizer maps to a reflection `C2`; its
conjugates normally generate `S3`.  By contrast, the isolated `C3` block
actions in `6T4` and `6T6` have trivial point stabilizer in the induced
`C3` action, so the ternary tower itself is a forbidden regular quotient.

This is the rigorous `C2*C3` interpretation: the problem is not that a
ternary tree exists, but that its reflection sidecar has been forgotten.
The cyclic `C3` face is killed by the Keller gate; the coupled dihedral `S3`
face survives it.

## 5. The quartic resolvent comparison is action-sensitive

The same abstract group behaves differently on different sheet sets.

- Natural quartic `A4` has `H=C3` and `N=A4`, so it passes.  Sextic `6T4`
  has `H=C2`, `N=V4`, and regular quotient `C3`, so it is excluded.
- Natural quartic `S4` has `H=S3` and passes.  Its two sextic actions `6T7`
  and `6T8` have respectively a nonnormal `V4` and a `C4` stabilizer; both
  normally generate `S4` and pass.

Therefore the three-matchings resolvent action is not automatically a block
quotient of the four-root action.  There is no canonical equivariant map
from one root to one matching.  The identity

```text
disc(quartic)=disc(resolvent cubic)                        (10)
```

says that both actions see the same sign character; it does not identify
their point stabilizers, affine inverse branches, or Jelonek components.
THM-2633 applies to a resolvent quotient only when the original root
stabilizer is killed by an actual regular sheet quotient.  This is the exact
action sidecar missing from a purely discriminantal transfer.

## 6. Scope and exact reproduction

This theorem excludes monodromy **types**, conditional on a Keller map having
the stated finite generic degree.  It gives no degree bound and does not
exclude the remaining quintic/sextic types or quartic `A4,S4`.  It proves no
Jacobian or Dixmier conjecture.

Run

```bash
python 04-computation/keller_degree_five_six_stabilizer_census_thm2643.py
python -O 04-computation/keller_degree_five_six_stabilizer_census_thm2643.py
```

The companion starts from explicit standard permutation generators.  It
independently closes `G`, the point stabilizer, and its normal closure;
checks every outside cycle type; enumerates all `10+15` sextic block systems
plus the singleton quotient; verifies (5) in all twenty-one degree-five/six
actions; and checks natural quartic `A4,S4` as action-sensitive controls.
Normal and optimized runs must agree after LF normalization and end in
`ALL CHECKS PASSED`.
