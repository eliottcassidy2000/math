---
id: THM-3381
title: "Reflected-residue affine phase transport and frozen-tree stability"
status: >
  PROVED analytic lemma + VERIFIED-EXACT controls + independently audited.
  Changing a reflected LRC drift z to z+h transports its untruncated danger
  comb by an exact affine map carrying the centered phase defect
  Delta=hj-mL.  The clipped clauses satisfy an explicit L1 bound, and every
  frozen Hunter tree has a degree-weighted stability inequality.  A one-unit
  residue perturbation makes a canonical positive pair overlap vanish,
  refuting residue-blind extension of THM-3360; a distinct noncanonical packet
  has a strictly positive coarse stability margin, so the repaired chamber is
  nonvacuous.  Arbitrary residues, changing trees/cells, arbitrary k=1, entry,
  the rung, and LRC(14) remain open.
source: root/repository-archaeology-second-pass-2026-08-14 selective recovery
recovered_from: origin/codex/lrc-math-20260812@73860320139774bca3313c986794d1a4ebf2db37
depends_on:
  - THM-3355-disconnected-low-affine-tail-and-reflected-branch-closure
  - THM-3360-uniform-reflected-high-pair-floor-by-admissible-affine-tails
related:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-3378-projected-k3-z216-gcd24-L129360-row94-one-high-torsion-closure
script: 04-computation/lrc14_reflected_residue_phase_stability_thm3381.py
output: 05-knowledge/results/lrc14_reflected_residue_phase_stability_thm3381.out
script_sha256: fedaa46451ef1582e9d6a472e4996380ad08060810c1cf2642c4b9121ce3f383
output_sha256: b41ef513589888543c3899ab0e11b684bbd0f3075f802c88adbdab09d223e4ef
hash_basis: LF-normalized bytes
---

# THM-3381 -- reflected residues need a located phase budget

**PROVED analytic lemma + VERIFIED-EXACT controls + independently audited.**

## 1. Exact affine transport

Let `L>0`, `0<=j<L`, `z>0`, and let `h` be an integer with `z+h>0`.  On the
real line define the open periodic danger comb and its clipped cell clause by

```text
B_z(j)={u in R : ||z(j+u)/L||<1/14},
A_z(j)=B_z(j) intersect [0,1].                           (1)
```

For any integer `m`, put

```text
Delta=hj-mL,
phi_h(u)=(zu-Delta)/(z+h).                              (2)
```

Then

```text
(z+h)(j+phi_h(u))/L=z(j+u)/L+m.                         (3)
```

Since `z+h>0`, `phi_h` is an affine bijection.  Equation `(3)` proves the
exact untruncated transport law

```text
B_(z+h)(j)=phi_h(B_z(j)).                               (4)
```

The intersection with `[0,1]` is not preserved by `(4)` and must still be
computed or bounded literally.  This is the first scope boundary.

## 2. Phase-aware symmetric-difference bound

Choose `m` in `(2)` so that `Delta` is a centered residue modulo `L`, and set

```text
eta=(|Delta|+|h|)/L.                                    (5)
```

At the same cell coordinate `u in [0,1]`, the new and old circle phases differ
modulo an integer by

```text
(Delta+hu)/L,
```

whose absolute value is at most `eta`.  If the two radius-`1/14` indicators
differ, the old phase lies in the `eta`-neighbourhood of one of the two danger-
arc boundary points.  This boundary neighbourhood has circle measure at most
`min(1,4eta)`.

The old phase traverses an interval of length `alpha=z/L` as `u` traverses
the cell.  Every complete turn contributes the circle measure of the boundary
set and the final partial turn contributes at most one more copy.  Therefore

```text
sigma(z,h):=mu(A_z(j) symmetric_difference A_(z+h)(j))
 <=min(1,4eta(1+1/alpha))
 = min(1,4eta(1+L/z)).                                  (6)
```

This proof uses neither a canonical residue nor a small value of `|h|` by
itself.  The located centered phase `Delta` is load-bearing.

## 3. Frozen-tree stability

Freeze a spanning tree `T` on six clauses and define its Hunter margin

```text
M_T(A)=sum_(ij in T) mu(A_i intersect A_j)
       -[sum_i mu(A_i)-6/7].                            (7)
```

Let `A_i'` be perturbed clauses and
`sigma_i=mu(A_i symmetric_difference A_i')`.  The elementary inequalities

```text
mu(A_i' intersect A_j')
 >=mu(A_i intersect A_j)-sigma_i-sigma_j,
mu(A_i')<=mu(A_i)+sigma_i                               (8)
```

sum to

```text
M_T(A')
 >=M_T(A)-sum_i (deg_T(i)+1)sigma_i.                    (9)
```

Combining `(6)` and `(9)` gives a rigorous phase-aware perturbation chamber
for any already-certified fixed packet and fixed tree.  Several coordinates
may move simultaneously.  The theorem does not say that the old tree remains
optimal, only that its displayed lower bound remains valid.

## 4. Minimal residue-blind hostile

Use the current `649`-body bank member

```text
E=(1,2,3,4,6,12),        (L,j)=(168,90),               (10)
```

and its canonical high `3:5` pair on labels `(12,4)`:

```text
z=3L-12=492,              w=5L-4=836.                  (11)
```

Exact open-interval arithmetic gives

```text
mu(A_z intersect A_w)=6/209>Dmax/5.                    (12)
```

Move only the first drift by the smallest nonzero integer displacement,
`z'=491`.  Then

```text
mu(A_z' intersect A_w)=0.                               (13)
```

For `h=-1,m=-1`, equation `(2)` gives `Delta=78`, a phase jump
`78/168=13/28`.  Thus `(13)` is not a continuity paradox.  A one-unit drift
can be a near-half-turn in the selected cell.  This witness is minimal in
changed-coordinate support and absolute nonzero integer displacement.

It refutes every proposed extension of THM-3360's canonical pair floor which
retains only a level ratio, denominator/unit passport, or `|h|` and forgets
the located phase.  It is a counterexample to a sufficient pair certificate,
not to LRC(14), and it does not alter THM-3360's canonical statement.

## 5. A strict positive noncanonical chamber

Take the distinct-drift repeated-level packet

```text
E=(1,2,3,8,9,14),       (L,j)=(7056,3780),
q=(3,3,3,5,5,5).                                       (14)
```

The exact maximum cross-`K_(3,3)` tree is

```text
T=((1,4),(1,5),(0,5),(2,3),(0,3))                      (15)
```

with canonical margin

```text
71440713312252560278 / 527394888495258135905 > 0.       (16)
```

Move only `z_0=3L-1` by `h=28`.  Here `28j=15L`, so `Delta=0` and
`eta=1/252`.  The literal and analytic bounds are

```text
sigma_0       =67424/16616095,
sigma_0 bound =28223/1333521.                           (17)
```

Vertex zero has tree degree two.  Even substituting the coarser value from
`(17)` into `(9)` leaves

```text
113864784228062404699 / 1582184665485774407715 > 0.     (18)
```

The actual perturbed margin is also positive.  Thus the theorem certifies a
genuinely noncanonical packet; the repaired chamber is not vacuous.

## 6. Typed contract, boundary, and audit

```text
source:    a reflected periodic danger comb at one selected LRC cell
target:    the perturbed comb and a lower bound for one frozen Hunter tree
map:       phi_h with centered defect Delta=hj-mL
preserved: comb order, exact untruncated endpoints, fixed tree and labels
destroyed: clipping under the affine map, tree optimality, cell selection,
           canonical singleton formulas and physical owner/entry data
sidecars:  Delta, h, z, literal singleton masses and tree degrees
tests:     the one-unit hostile (10)--(13) and positive packet (14)--(18).
                                                                    (19)
```

The recovered standard-library companion uses exact `Fraction` interval
arithmetic, checks both controls and every displayed rational identity, and
contains no optimization-dependent `assert` gate.  Normal and optimized runs
byte-match the stored output.  An independent line audit rederived `(3)`,
`(6)`, and `(9)` and found no type or inequality-direction defect.

Reproduce with

```text
python3 04-computation/lrc14_reflected_residue_phase_stability_thm3381.py
python3 -O 04-computation/lrc14_reflected_residue_phase_stability_thm3381.py
```

The all-body coherent-shift banks on the stale source branch remain
`FINITE-EXACT` sidecars, not consequences of this analytic lemma.  Large
shifts can invalidate a restricted cross-tree even when another full-tree
certificate survives.  Arbitrary residue packets, changing cells and trees,
arbitrary `k=1`, projected-to-physical transport, entry, the rung, and LRC(14)
remain **OPEN**.

**QED.**
