---
id: THM-1242
title: THE SEVEN-WALL PAIR-BEAT SUNFLOWER LAW — six common-zero masks force a hole on every nontrivial common clock, while the first mixed full clock is the sharp q=15 sunflower
status: PROVED (all-Q common-clock escape; exact seven/eight-wall criticality; minimal mixed full clock and q=15 sunflower; dependency-free referee; sorry-free Lean finite core)
source: codex-2026-07-19-S78 continuation with path-inequality agent
depends_on: [THM-1216, THM-1217]
related: [THM-1219, THM-1237, THM-1238]
script: 04-computation/lrc14_seven_wall_pair_beat_sunflower_thm1242.py
output: 05-knowledge/results/lrc14_seven_wall_pair_beat_sunflower_thm1242.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCSevenWallBeatSunflower.lean
script_sha256: 56f1b5ea957a4170c9fcf9914776022083ded3b494f0658000d943385659cb77
output_sha256: a20eca7fc135d7e8a9a84ebfea8d3b0b6ee70d2f7cb9e1f2e14ea0ce24ccb8fd
formalization_sha256: 762f4cbace3388e9703983500f5597d2e35ebb5264a06b3f1f95bfd34a1ea9b1
---

# THM-1242 — seven-wall pair-beat sunflower law

## 1. Seven-wall common-clock theorem

Let `s1,...,s7` be seven positive integer wall speeds, select a defining pair,
and let

```text
q=|si-sj|                    or                    q=si+sj. (1)
```

At every beat `t=p/q`, the two defining strict danger masks coincide by sign
symmetry.  Thus there are six distinct mask obligations: the pair mask and
the other five wall masks.

Assume a common gcd

```text
gcd(s_l,q)=d                 for l=1,...,7,             (2)
Q=q/d>=2.
```

Put

```text
h=ceil(Q/14),
A(Q)=2h-1,
B7(Q)=6A(Q)-4=12h-10.                                  (3)
```

Then

```text
B7(Q)<=Q                                                    (4)
```

for every `Q>=2`, and every `B7(Q)` consecutive beat numerators contain a
residue safe for all seven walls.

Consequently, if a closed interval `I` of length `L` contains a consecutive
`q`-beat block and

```text
qL>=B7(Q),                                               (5)
```

then `I` contains a seven-wall-safe beat.  On a protected needle
`L>=1/(7m)`, the clean supplier

```text
q/m>=7B7(Q)                                              (6)
```

suffices.

## 2. Proof and exact cardinality criticality

Each common-clock unit mask has exactly `A(Q)` residues and contains zero.
Adjoining six such masks loses at least one unit at each of five insertions,
so their union `U` obeys

```text
|U|<=1+6(A-1)=6A-5=B7(Q)-1.                            (7)
```

For `2<=Q<=14`, `h=1` and `B7=2<=Q`.  For `h>=2`, the first modulus in the
ceiling shell is `14h-13`, hence

```text
Q-B7 >=(14h-13)-(12h-10)=2h-3>=1.                      (8)
```

This proves (4).  A block of at most `Q` consecutive integers has distinct
residues; by (7), one of its first `B7` residues lies outside `U`.  This proves
the escape statement and (5)--(6).

The count is exactly critical in the number of walls.  For `R` walls a
defining pair leaves `R-1` obligations, so the common-zero threshold is

```text
B_R(Q)=2+2(R-1)(h-1).                                   (9)
```

At `R=7`, this is (3).  At `R=8`,

```text
B_8(Q)=14h-12.                                         (10)
```

For every noninitial shell `h>=2`, its first modulus is

```text
Q=14h-13,
B_8(Q)=Q+1.                                            (11)
```

Thus seven walls are the last cardinality at radius `1/14` for which common
zero incidence alone forces a hole on every nontrivial common reduced clock.
The first eight-wall failures are `Q=15,29,43,57,71,...`.  The common clock
`Q=1` remains the unique seven-wall obstruction.

## 3. Mixed clocks: the sharp q=15 sunflower

Nontrivial individual quotient periods do not imply a mixed-clock hole.  Take

```text
q=15,
(s1,...,s7)=(1,16,5,3,2,4,7),                         (12)
```

with defining difference `16-1=15`.  On `Z/15Z`, the six distinct masks are

```text
M_pair ={0,1,14},
M_5    ={0,3,6,9,12},
M_3    ={0,5,10},
M_2    ={0,7,8},
M_4    ={0,4,11},
M_7    ={0,2,13}.                                      (13)
```

Their quotient periods are

```text
(15,15,3,5,15,15,15),                                  (14)
```

so no wall has period one.  Nevertheless

```text
union M_i=Z/15Z.                                        (15)
```

More precisely, every two distinct masks intersect exactly in `{0}`, and
their nonzero petals partition the fourteen nonzero residues with sizes

```text
2,4,2,2,2,2.                                           (16)
```

The common-zero cardinality cap and every Hunter-tree cap are both saturated:

```text
sum_i |M_i|-5=3+5+3+3+3+3-5=15.                       (17)
```

This full clock is minimal among mixed clocks with no quotient period one.
If the master modulus `L<=14`, every divisor period `Q>1` also satisfies
`Q<=14`, so its strict reduced danger window is `{0}`.  Its lift is the proper
subgroup

```text
{p in Z/LZ:p=0 mod Q},                                 (18)
```

which misses the generator residue one.  No union of such lifts can fill the
clock.  Equation (13) shows that `L=15` is the exact first failure.

## 4. Structural and tournament audit

The common-clock proof uses mask obligations plus one distinguished incidence
point.  The mixed obstruction is a six-petal sunflower embedded in a cyclic
word.  Mask size alone sees the equality in (17) but not which residue is the
hole—or that no hole remains.

For Tournament Analysis, take pairwise observable

```text
w_ij=|Mi intersect Mj|.                                (19)
```

In the sunflower every `w_ij=1`; any tournament orientation is pure tie
gauge.  Ordering by period or speed gives a transitive tournament with one
Hamiltonian path, but all orientations destroy the petal placement.  The
faithful object is the common-zero sunflower together with cyclic residue
addresses and the consecutive beat block.

We challenged runners, pair beats, quotient clocks, master residues, masks,
intersection-tree edges, petals, and proof obligations as vertices.  The
smallest faithful carrier is

```text
(six lifted masks, common zero, cyclic petal word, block phase).           (20)
```

THM-1237's positioned forest and THM-1242's cyclic sunflower are orthogonal
projections: one retains interval discrepancy and loses residue placement;
the other retains exact beat placement and loses off-grid coverage.

## 5. Verification and scope

The dependency-free referee checks (4) for all `Q=2,...,10000`, `4,023,512`
common-clock mask unions through `Q=70`, the all-wall critical formula, every
small-master lifted subgroup, and the complete sunflower identities
(13)--(17).  Normal and optimized outputs are byte-identical.

The Lean module kernel-checks the six-mask common-zero cap, seven-wall
threshold, exact eight-wall shell failure, escape arithmetic, and the full
`Fin 15` sunflower union/intersection/cardinality ledger using ordinary kernel
reduction, not `native_decide`.  THM-1217 supplies the speed-to-mask bridge;
there are no proof placeholders.

Frozen hashes are

```text
source         56f1b5ea957a4170c9fcf9914776022083ded3b494f0658000d943385659cb77
output         a20eca7fc135d7e8a9a84ebfea8d3b0b6ee70d2f7cb9e1f2e14ea0ce24ccb8fd
formalization  762f4cbace3388e9703983500f5597d2e35ebb5264a06b3f1f95bfd34a1ea9b1
```

THM-1242 proves the common-clock theorem and identifies the exact mixed-clock
failure.  It does not force a useful pair beat on every protected needle,
exclude the `q=15` sunflower inside a global packet, or prove LRC(14).
