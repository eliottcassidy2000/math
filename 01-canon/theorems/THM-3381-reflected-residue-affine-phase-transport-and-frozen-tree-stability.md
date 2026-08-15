---
id: THM-3381
title: "Reflected-residue affine phase transport and frozen-tree stability"
status: >
  PROVED analytic theorem + VERIFIED-EXACT controls. Exact affine conjugacy
  retains the centered phase defect under an arbitrary integer drift change;
  a quantitative symmetric-difference bound yields a sufficient frozen-tree
  Hunter stability gate. A unit residue change exactly refutes the
  residue-blind pair-floor transfer, while a strict noncanonical packet
  verifies that the repaired gate is nonvacuous. No arbitrary-residue,
  projected-wall, physical-entry, or LRC(14) closure follows.
source: root/lrc-math-2026-08-14
audit: >
  independent line audit of the conjugacy, boundary-neighbourhood pullback,
  and tree coefficients; exact-rational hostile and positive controls;
  byte-identical ordinary/optimized replay
depends_on: []
related:
  - THM-3355-disconnected-low-affine-tail-and-reflected-branch-closure
  - THM-3360-uniform-reflected-high-pair-floor-by-admissible-affine-tails
  - THM-3382-independent-superunit-affine-tail-and-reflected-residue-closure
script: 04-computation/lrc14_reflected_residue_phase_stability_20260814.py
output: 05-knowledge/results/lrc14_reflected_residue_phase_stability_20260814.out
script_sha256: fedaa46451ef1582e9d6a472e4996380ad08060810c1cf2642c4b9121ce3f383
output_sha256: b41ef513589888543c3899ab0e11b684bbd0f3075f802c88adbdab09d223e4ef
hash_basis: LF-normalized bytes
---

# THM-3381 -- affine phase transport and frozen-tree stability

**PROVED analytic theorem + VERIFIED-EXACT controls.**

This theorem isolates the coordinate lost when the matched-residue floor of
THM-3360/3376 is transported to another residue.  The stable object is a
located phase together with a tree margin, not a reduced level ratio.

## 1. Exact affine conjugacy

Let `L>0`, `0<=j<L`, `z>0`, and let the integer `h` satisfy `z+h>0`.  On the
real line define

```text
B_z(j)={u in R: ||z(j+u)/L||<1/14},
A_z(j)=B_z(j) intersect [0,1].                         (1)
```

For any integer `m`, put

```text
Delta=hj-mL,                 phi_h(u)=(zu-Delta)/(z+h). (2)
```

Then

```text
(z+h)(j+phi_h(u))/L=z(j+u)/L+m,
B_(z+h)(j)=phi_h(B_z(j)).                              (3)
```

Thus `(3)` is an exact orientation-preserving affine conjugacy of the
untruncated periodic combs.  Intersecting with `[0,1]` remains literal; the
endpoint clipping is not silently identified under `phi_h`.

## 2. Symmetric-difference transport

Choose `m` so that `Delta` is a centered residue modulo `L`, and set

```text
eta=(|Delta|+|h|)/L,             alpha=z/L.             (4)
```

At the same `u in [0,1]`, the circular displacement of the old and new
phases is at most `eta`.  Their danger indicators can differ only when the
old phase lies within circular distance `eta` of one of the two boundary
points of the radius-`1/14` arc.  This boundary neighbourhood has circle
measure at most `min(1,4eta)`.

The old phase traverses an interval of length `alpha`.  For a periodic set of
circle measure `s`, its intersection with that interval has length at most
`(floor(alpha)+1)s`; pulling back by slope `alpha` and clipping at one gives

```text
sigma(z,h):=mu(A_z(j) triangle A_(z+h)(j))
 <= min(1,4eta(1+1/alpha))
 =  min(1,4eta(1+L/z)).                               (5)
```

No reflected-residue or high-pair hypothesis is used in `(3)--(5)`.

## 3. Frozen-tree Hunter stability

Let `A_1,...,A_6` be measurable clauses, let `T` be a spanning tree on their
six vertices, and define its Hunter margin

```text
M_T(A)=sum_(ij in T) mu(A_i intersect A_j)
       -[sum_i mu(A_i)-6/7].                           (6)
```

Suppose `mu(A_i triangle A_i')<=sigma_i`.  For every edge,

```text
mu(A_i' intersect A_j')
 >=mu(A_i intersect A_j)-sigma_i-sigma_j,              (7)
```

and `mu(A_i')<=mu(A_i)+sigma_i`.  Summing `(7)` over the tree and the
singleton debt gives the exact gate

```text
M_T(A') >= M_T(A)-sum_i (deg_T(i)+1)sigma_i.           (8)
```

Consequently a previously certified tree survives any simultaneous integer
drift perturbation for which its old margin strictly exceeds the right-hand
loss budget obtained from `(5)`.  The degree term records edge loss; the
additional `+1` is the changed singleton debt.

## 4. The residue-blind pair floor is false

On the upper-median cell

```text
E=(1,2,3,4,6,12),             (L,j)=(168,90),
```

the canonical high `3:5` pair on labels `(12,4)` is

```text
z=3L-12=492,                  w=5L-4=836.
```

Exact interval arithmetic gives

```text
mu(A_z intersect A_w)=6/209>Dmax/5.                   (9)
```

Changing only `z` to `z-1=491` gives

```text
mu(A_(z-1) intersect A_w)=0.                          (10)
```

For `h=-1,m=-1`, the centered defect is `Delta=78`, hence the cell-phase
jump is `13/28`.  This is a minimal one-coordinate, unit-displacement
counterexample.  It refutes any transfer that retains only the old high
ratio, `|h|`, or canonical edge floor while discarding `Delta`.  It is not an
LRC packet counterexample.

## 5. The repaired gate is nonvacuous

For

```text
E=(1,2,3,8,9,14),  (L,j)=(7056,3780),
q=(3,3,3,5,5,5),
```

the displayed cross-`K_(3,3)` tree in the companion has margin

```text
M=71440713312252560278/527394888495258135905.          (11)
```

Perturb the first drift by `h=28`.  Since `28j=15L`, `Delta=0` and
`eta=1/252`.  The exact symmetric difference and the bound `(5)` are

```text
sigma=67424/16616095
 <=28223/1333521.                                     (12)
```

That vertex has tree degree two, and even the coarse gate is strict:

```text
M-3(28223/1333521)
 =113864784228062404699/1582184665485774407715>0.      (13)
```

The literal perturbed margin is also positive.  Hence `(8)` certifies a
genuinely noncanonical, distinct-drift packet; it is not merely a formal
repair of `(10)`.

## 6. Audit and scope

The exact-rational companion reconstructs every interval in `(9)--(13)`,
checks the canonical target `Dmax/5`, the hostile phase defect, the selected
tree and its degree, the exact/coarse stability floors, and the literal
perturbed margin.  It contains no floating-point arithmetic or
assert-dependent gate.  Ordinary and optimized outputs are byte-identical.
Reproduce with

```bash
python3 04-computation/lrc14_reflected_residue_phase_stability_20260814.py
python3 -O 04-computation/lrc14_reflected_residue_phase_stability_20260814.py
```

The source and output hashes are pinned in the frontmatter.  The larger
coherent-shift censuses in the related reflection are finite applications and
boundary diagnostics, not part of the universal analytic statement here.

This theorem supplies a sufficient perturbation chamber.  It does not prove
that every noncanonical packet lies in that chamber, does not preserve a
cross-only tree under large shifts, and does not close arbitrary `k=1`, a
projected wall, the physical entry, or LRC(14).

**QED.**
