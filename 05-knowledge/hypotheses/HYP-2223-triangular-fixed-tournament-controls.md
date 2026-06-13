---
id: HYP-2223
name: triangular-fixed-tournament-controls
status: OPEN
session: codex-2026-06-04-S647
parents: [HYP-2220, HYP-2200, HYP-2202, HYP-2203, HYP-2217, THM-115]
---

# HYP-2223: Perfect Tournament Controls Turn 7/21 into Spine/Fiber Dimensions

## Claim

Triangular fixed controls from HYP-2220 become especially sharp in tournament
language.  Every tournament on `m` vertices has

```text
A_m = C(m,2)
```

arcs.  If one Hamiltonian path is chosen as a fixed section, then `m-1` arcs
are forced and the remaining off-path deformation fiber has

```text
C(m,2) - (m-1) = C(m-1,2)
```

free arcs.  Thus a fixed path decomposes the triangular carrier as:

```text
C(m,2) = (m-1) + C(m-1,2).
```

At the second even-perfect tournament control,

```text
C(8,2) = 28 = 7 + 21.
```

The two permanent Hamiltonian-path gap values `7` and `21` appear here as
dimensions: `7` is the forced path spine length, and `21` is the off-path
triangular fiber size.  The exhaustive `n=8` H-spectrum still says neither
`7` nor `21` is realized as a Hamiltonian-path count.

So the proposed guardrail is:

> Before transferring a scalar across H-gaps, LRC, or unit-distance carriers,
> decide whether the scalar is a count, section length, deformation dimension,
> shell modulus, or side-channel mass.

## S647 Evidence

`04-computation/triangular_fixed_tournament_controls_s647.py` stores
`05-knowledge/results/triangular_fixed_tournament_controls_s647.out`.

The script proves and computes:

1. Tournament arc counts are triangular:

```text
A_m = C(m,2).
```

2. Fixing a Hamiltonian path leaves a smaller triangular deformation budget:

```text
fixed spine = m-1
off-path fiber = C(m-1,2)
fixed-path fiber size = 2^C(m-1,2).
```

3. Euclid-Euler perfect numbers are tournament arc-count controls:

```text
m=4:   C(4,2)=6    = 3  + 3
m=8:   C(8,2)=28   = 7  + 21
m=32:  C(32,2)=496 = 31 + 465
m=128: C(128,2)=8128 = 127 + 8001
```

where the perfect rows occur when `m=2^p` and `2^p-1` is Mersenne prime.

4. The `m=8` fixed-path fiber is exact without a new `2^21` enumeration.  By
double-counting pairs `(T,P)` where `P` is a Hamiltonian path of `T`, if
`L_h` is the labelled count of `8`-vertex tournaments with `H(T)=h`, then one
fixed labelled base path sees

```text
fixed_count_h = h * L_h / 8!.
```

Using the full labelled `n=8` spectrum from monad S1, S647 reconstructs:

```text
total fixed-path tournaments = 2097152 = 2^21
distinct H values = 320
H range = [1, 661]
mean H in fixed-path fiber = 388.6875
fixed_count(H=7) = 0
fixed_count(H=21) = 0
```

Selected fixed-fiber counts:

```text
H=1   count=1
H=3   count=6
H=5   count=25
H=35  count=140
H=39  count=182
H=49  count=735
H=63  count=126
H=189 count=3798
H=661 count=3966
```

Thus `7` and `21` are not absent because the perfect `28` control cannot see
them.  It sees them exactly, but as path/fiber dimensions rather than as
values of `H(T)`.

## Interpretation

This reframes the persistent gap pair `{7,21}`.

The naive scalar story says:

```text
7 and 21 are special numbers.
```

The carrier story says:

```text
7  = section length at the m=8 perfect tournament control
21 = off-section triangular deformation dimension at the same control
28 = total perfect arc count
```

That is a role mismatch, not a role equality.  A number can be load-bearing
inside a carrier without being realizable in a different output coordinate.

This aligns with the unit-distance and LRC lessons:

- Unit-distance `n=21` can use `21` as vertex count and `57=20+37` as
  spine/bulk edge carrier without realizing literal tournament `H=21`.
- LRC `n=14` uses `C=27` as a shell clock and `C(14,2)=91` as a triangular
  pair-count carrier whose aliquot shadow is `21`, but the proof burden
  remains lift/carry/owner data.
- Tournament `H=21` is forbidden because path-count output obeys strong
  component and OCF side channels, not because `21` never appears elsewhere.

## Connection to HYP-2220

HYP-2220 identified even perfect numbers as aliquot fixed controls:

```text
s(C(2^p,2)) = C(2^p,2)
```

when `2^p-1` is Mersenne prime.  HYP-2223 adds the Hamiltonian-path section:

```text
C(2^p,2) = (2^p-1) + C(2^p-1,2).
```

For `p=3`, this becomes:

```text
28 = 7 + 21.
```

This makes the `7/21` pair a section/fiber decomposition of a perfect
tournament control, while THM-115 says `{7,21}` are permanent `H`-gaps.

## Use Toward Main Proofs

The practical move is not to treat `28=7+21` as a proof of H-gap
impossibility.  The practical move is to build proof ledgers with typed scalar
roles:

```text
raw count
fixed section length
off-section deformation dimension
shell modulus
aliquot shadow
side-channel mass
```

For LRC `n=14`, this suggests rewriting candidate quotient checks in a
section/fiber style:

1. Choose or prove the existence of a section that plays the role of a fixed
   Hamiltonian path: a witness shell, spine, orbit representative, or owner
   route.
2. Put the remaining uncertainty into an explicit deformation fiber.
3. Track which visible scalars belong to the section, which to the fiber, and
   which to the final predicate.

This is close to HYP-2217's retained-section/bulk split and HYP-2203's
unit-spine gauge separation.  The new value is the exact tournament control:
at `m=8`, the known H-gap values are dimensions of the perfect arc carrier.

## Tournament Analysis

S647 uses proof lenses as Tournament Analysis vertices.  The lens tournament is
transitive:

```text
vertices=7
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
hamiltonian_paths=1
```

Tie Hamiltonian path:

```text
off_path_triangular_fiber
n8_perfect_control_28_equals_7_plus_21
fixed_hamiltonian_section
h_gap_role_mismatch_guardrail
weighted_fixed_path_fiber_from_h_spectrum
euclid_euler_arc_count_controls
raw_7_21_numerology
```

Pairwise observable: which lens best preserves the distinction between
Hamiltonian-path count, fixed section, off-path fiber, and scalar control.
Switch/gauge: proof transfer value minus scalar-numerology risk.

## Assumption Challenge

Alternate vertices considered:

```text
tournament vertices
arcs
fixed path edges
off-path cells
Hamiltonian paths
conflict-graph cycles
LRC shell residues
unit-distance spine/bulk obligations
proof lenses
```

Chosen computation: off-path cells are the finite tiling variables, while
Tournament Analysis uses proof lenses.

Preserved information: path-section size, deformation-fiber size, exact
fixed-path H-spectrum through the perfect `m=8` control, and the `7/21` gap
guardrail.

Destroyed information: tournament isomorphism class, endpoint distribution,
SCC decomposition, Omega side channels, and LRC carry/owner labels.

Challenged assumption: a scalar appearing in several problems has the same
role in each.  Here `7` and `21` are forbidden H-counts but natural dimensions.

## Next Steps

1. Add a typed scalar ledger to cross-problem scripts: `count`, `section`,
   `fiber`, `modulus`, `aliquot`, `bulk`, `owner/carry`.
2. Re-run the H-gap carrier packets with a role label attached to every
   occurrence of `7`, `21`, `27`, and `28`.
3. For LRC `n=14`, compare the `27` shell clock with the punctured perfect
   control `28-1`, but keep this as a typed analogy until a carry/owner
   statement is proved.
4. For unit-distance, mark `n-1` unit-spine edges as a section and
   `u(n)-(n-1)` as fiber/bulk before mapping to tournament observables.

## See Also

`04-computation/triangular_fixed_tournament_controls_s647.py`;
`05-knowledge/results/triangular_fixed_tournament_controls_s647.out`;
`07-reflections/triangular-fixed-tournament-controls-s647.md`;
HYP-2220; HYP-2200; HYP-2202; HYP-2203; HYP-2217; THM-115.
