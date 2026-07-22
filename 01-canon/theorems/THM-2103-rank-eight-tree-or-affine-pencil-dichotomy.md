---
id: THM-2103
title: "The rank-eight tree-or-affine-pencil dichotomy fails on six exact dyadic rows"
status: >
  PROVED NO-GO plus exact bounded classification. The conjecture that every
  eight-terminal transverse plane has either maximum restricted-overlap tree
  weight greater than 5/49 or a nondegenerate signed affine pencil is false.
  In the complete dyadic coefficient cube c_i=(x_i,2^i), x_i in {0,1,2},
  exactly seven of 3^8 rows have tau<=5/49: one is an affine pencil and six
  are not. All six nonpencil no-go rows nevertheless have an explicit safe
  point of denominator six (five at (1/2,1/6), one at (1/6,1/6)). In the
  complete consecutive cube c_i=(x_i,i+1), x_i in {-1,0,1}, no row has
  tau<=5/49. Thus affine rank is still insufficient; the next carrier must
  retain a small-clock residue alternative. This does not prove rank eight
  or LRC(14).
source: codex-2026-07-22-LRC-rank-eight-tree-pencil
depends_on:
  - THM-2098
  - THM-2099
related:
  - THM-1065
  - THM-2069
  - THM-2096
  - THM-2104
  - THM-2126
  - THM-2114
script: 04-computation/lrc14_rank8_tree_pencil_dichotomy_scout_codex_20260722.py
output: 05-knowledge/results/lrc14_rank8_tree_pencil_dichotomy_scout_codex_20260722.out
script_sha256: 99afc601febd8c258dc3b524ecc560f2542a5eda4dff641ed2ffad518e51a8e3
output_sha256: f6d7fd654c06c1a4a46c0e807b4977bceb28bd5f24546d3b15e4997d5d571ad2
hash_basis: working-tree files with LF line endings
---

# THM-2103 -- the tree-or-pencil conjecture is false

Let the guard be

```text
g=(1,0),
```

and for eight transverse terminal columns put

```text
w_ij=measure({||c_i.X||<1/14,
              ||c_j.X||<1/14,
              ||g.X||>1/7}),
tau=max_(spanning trees T) sum_(ij in T)w_ij.           (1)
```

THM-2098 says an eight-band transverse cover must have `tau<=5/49`.
THM-2099 says every nondegenerate signed affine pencil has a mixed-safe
point. It was therefore natural to propose

```text
tau>5/49, or some sign gauge makes the c_i a
nondegenerate affine pencil.                            (2)
```

Equation (2) is false.

## 1. An explicit nonpencil counterexample to the certificate

Take

```text
(c_0,...,c_7)=
((0,1),(2,2),(1,4),(0,8),(2,16),(1,32),(0,64),(1,128)). (3)
```

The exact THM-2099 primitive-relation formula and Kruskal's algorithm give

```text
tau=2101896757/20742005760
   =5/49-14634443/20742005760
   <5/49.                                               (4)
```

No choice among the `2^7` relative terminal signs makes the eight signed
columns collinear on a nondegenerate affine line. Thus (3) satisfies neither
alternative of (2).

It is not a mixed cover. At

```text
X=(1/2,1/6),                                           (5)
```

the guard has distance `1/2`, and the eight terminal distances are

```text
(1/6,1/3,1/6,1/3,1/3,1/6,1/3,1/6).                   (6)
```

Hence all nine characters are strictly safe.

## 2. Complete two-cube classification

The companion exhausts two exact coefficient cubes. Distinct second
coordinates ensure that a sufficiently large parameter direction `(1,D)`
gives an odd positive guard and distinct positive terminal specializations.

### Consecutive cube

For

```text
c_i=(x_i,i+1),                 x_i in {-1,0,1},         (7)
```

all `3^8=6561` rows satisfy `tau>5/49`. The exact minimum is

```text
tau_min=3167/30870=5/49+17/30870.                      (8)
```

### Dyadic cube

For

```text
c_i=(x_i,2^i),                 x_i in {0,1,2},          (9)
```

the exact classification is

```text
tau<=5/49 rows                         7,
nondegenerate affine-pencil rows       1,
nonpencil rows                         6.              (10)
```

The six nonpencil rows and their maximum trees are printed in the frozen
transcript. Every row is independently checked against all `2^7` sign gauges.
Five have the safe point `(1/2,1/6)`; the remaining row has `(1/6,1/6)`.
Their terminal distances lie in `{1/6,1/3,1/2}`.

The small clock is not accidental. At denominator six, `2^i` enters the
eventual two-cycle `2,4 mod 6`, while the first coefficients `x_i mod 2`
choose whether the half-turn is added. The phase sees the residue-labelled
orbit that the unordered edge weights and affine rank both forget.

## 3. Correct frontier effect

The no-go chain is now exact:

```text
binary zero graph
  -> weighted maximum tree                 insufficient by THM-2099,
  -> weighted tree + signed affine rank    insufficient by THM-2103,
  -> weighted tree + affine rank + clock residue       live.     (11)
```

This does not say a denominator-six clock always exists when `tau<=5/49`.
The two cubes are finite diagnostics, and the primitive relation height is
unbounded. The lawful next target is a tree-or-bounded-clock inverse theorem,
with an analytic tail or a proof that the clock modulus is forced by the
low-edge relation graph.

The transfer from THM-2069 is typed: a rational torus point is an evaluation
residue for the character columns, but codeword support alone does not retain
the two unequal safety radii. The residue vector must be decorated by the
guard/terminal thresholds.

## 4. Assumption challenge and Tournament Analysis

The challenged assumption was that the missing sidecar after pair weights is
only affine rank. Row (3) preserves nonpencil rank and still defeats the tree
gate. Its denominator-six phase is the hidden coordinate.

The pair observable is (1). Kruskal ordering is the correct graphic-matroid
consumer. Turning the complete graph into a tournament by orienting toward
larger incident weight, with label tie-breaks, changes scores, cycles, SCCs,
edge flips, and Hamiltonian paths without changing (4). Conversely, those
fingerprints do not recover the residue vector at (5). Candidate vertices
challenged here were terminal columns, positive edges, relation types,
affine-pencil classes, clock residues, and proof obligations. The faithful
provisional carrier is the weighted relation graph plus signed affine rank
and a threshold-labelled clock evaluation. QED.

## 5. Exact referee and scope

The dependency-free companion imports THM-2099's exact edge provider, caches
every pair type, exhausts `2*3^8=13122` rows, runs Kruskal exactly with
`Fraction`, checks every relative sign gauge, and searches only denominator
at most fourteen for proof-readable safe points after a no-go row is found.
The six displayed safe phases are then checked exactly, not inferred from a
grid approximation.

This proves the failure of (2) and the finite classifications (8)--(10). It
does not prove a uniform clock bound, exclude a rank-eight cover, discharge
the vertical branch, or prove LRC(14).
