---
id: LEM-043
title: THE CONE FLOOR AND THE FULL-PERIOD WALL THEOREMS. (A) DEFECT BOUNDS: 1/49 − g/(7b) ≤ μ(D_a ∩ D_b) ≤ 1/49 + g/(7b). (B) EXACT CONE FLOORS: c₀(7/3) = 1/63 at (4,9), requiring only the reduced finite range b′ ≤ 31 because b′ ≥ 32 is closed by the defect tail; c₀(13) = 1/91 at (1,13), requiring only b′ ≤ 15 because b′ ≥ 16 is closed by the tail. (C) THE c = 7 DICHOTOMY: a sorted seven-block is lacunary (all adjacent ratios ≥ 7/3) or one adjacent pair supplies a 1/63 Hunter credit. (D) THE CONSECUTIVE c = 8 THEOREM: for every v ≥ 1, the seven adjacent pairs in {v,…,v+7} have reduced starts covering all residues mod 7; exactly the r = 0 and r = 6 edges equal 1/49 and the other five are strict, so their sum is > 1/7 and Hunter gives positive full-circle good measure. (E) DIFFERENCE COMB: D_a ∩ D_b ⊆ {t : ‖(b−a)t‖ < 1/7}. (F) SCOPE: none of these full-period credits may be inherited by an arbitrary child window without a separate positioning/aggregation lemma.
status: PROVED ((A) unimodal sampling bounds; (B) exact finite checks on b′ ≤ 31 and b′ ≤ 15 plus sharp tail inequalities; (C) Hunter + cone floor; (D) LEM-042 residue formula + seven-residue census + Hunter; (E) triangle inequality) + REFEREED EXACT (500 random defect pairs; exhaustive supersets of both finite cone scans; five seven-blocks; consecutive eight-block; beat censuses; clustered beat-miss 3/3)
source: boxeph-2026-07-17-S70 (owner directive: work the remaining proof surface; integrate incoming — opus lacunary branch + THM-956, klein hunter ledger, kps clustered closer, LEM-042)
depends_on: [LEM-042 (the exact pair formula), klein path_hunter_add_le (kernel-pure), opus S336 lacunary branch, kps cite_cluster7_lonely (the clustered face)]
script: 04-computation/lrc14_cone_floor_wall_boxeph_S70.py -> 05-knowledge/results/lrc14_cone_floor_wall_boxeph_S70.out
---

# LEM-043 — the cone floor; the wall dichotomy

At c = 7 the full-circle union budget is exactly zero, so one positive pair
credit crosses the wall.  The credit exists uniformly on the
lacunary-complement cone: `c0(7/3)=1/63`.  At c = 8, the coarse cone floor is
too small, but the exact consecutive formula closes every consecutive
eight-block.  These are full-period statements, not inherited-window
statements.

## Exact finite tails for the cone floors

After reducing `(a,b)` by `g=gcd(a,b)`, write `(a',b')=(a/g,b/g)`.  The
defect bound becomes

```text
mu(D_a cap D_b) >= 1/49 - 1/(7b').
```

For the `7/3` cone, exact integer evaluation on the finite set

```text
gcd(a',b')=1,  a'<b'<=31,  b'/a'<=7/3
```

has minimum `1/63`, attained at `(4,9)`.  Every omitted pair has `b'>=32`,
and hence

```text
mu >= 1/49 - 1/(7*32) = 25/1568 > 1/63.
```

Thus `b'<=31`, not `b'<=700`, is the complete finite obligation.  Likewise,
for the ratio-`13` cone the finite range `b'<=15` has minimum `1/91` at
`(1,13)`, while every `b'>=16` obeys

```text
mu >= 1/49 - 1/(7*16) = 9/784 > 1/91.
```

## Every consecutive eight-block crosses the full-period wall

Fix `v>=1` and let `A_i=D_{v+i}` for `0<=i<=7`.  Hunter's path inequality
and `mu(A_i)=1/7` give

```text
mu((Union_i A_i)^c)
  >= 1 - 8/7 + sum_{i=0}^6 mu(A_i cap A_{i+1}).
```

The starts `v,v+1,...,v+6` represent every residue modulo `7` exactly once.
By LEM-042, the edges with start residue `0` and `6` equal `1/49`, while the
other five are strictly greater than `1/49`.  Therefore

```text
sum_{i=0}^6 mu(A_i cap A_{i+1}) > 7/49 = 1/7,
```

and the complement has positive measure.  This is uniform as a theorem in
the integer start `v`; it does not assert a positive margin independent of
`v` (the excess tends to zero).  The earlier assertion that *every* edge is
strict was false: exactly five of the seven are strict.

## Inherited-window guardrail

All measures in the cone floors, Hunter calculation, and consecutive theorem
are Haar measures on the whole circle.  On an arbitrary inherited interval
`I`, neither `mu(D_v cap I)/mu(I)=1/7` nor a positive normalized pair credit
is available; LEM-042 gives explicit zero-credit windows.  Accordingly these
lemmas do not compose arbitrary blocks in an LRC14 tower without a separate
positioned-window or aggregation argument.

## Evidence log
- [x] defect bounds ±g/(7b) (500 exact pairs)
- [x] c₀(7/3) = 1/63 from b′ ≤ 31 + b′ ≥ 32 tail
- [x] c₀(13) = 1/91 from b′ ≤ 15 + b′ ≥ 16 tail
- [x] dichotomy on five blocks
- [x] every consecutive eight-block crosses: two equality edges, five strict
- [x] difference-comb subset law; beat censuses; clustered beat-miss 3/3
- [ ] named: Lean rendering — (A)/(B) are integer decides; (C) composes with
      path_hunter_add_le (already kernel-pure); the all-`v` consecutive
      residue formula still needs a kernel proof beyond the existing bounded check
