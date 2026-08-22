---
id: THM-3679
title: "LRC p-adic scalar seam and total-speed checksum"
status: >
  PROVED; PENDING INDEPENDENT HOSTILE AUDIT.  At every depth 13^n, four
  packet charts for each of two distinct blocker sources force a jointly
  target-zero address to be globally scalar.  On the all-coordinate-unit
  relation fibre this scalar hostile exists exactly when 13^n divides the
  total speed sum.  Hence for every positive integer row it dies by depth
  nu_13(sum w_i)+1.  THM-2334 realizes each higher difference pushforward by
  exact coordinate twists at a sufficiently delayed word clock.  No common
  grouped current across two legitimate owner-source semantics, phase
  noncancellation, covering-row exclusion, or LRC(14) conclusion is proved.
source: kps-s195 / THM-3676 Bockstein continuation, 2026-08-21
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-3676-lrc-cross-source-blind-intersection-and-scalar-seam
related:
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-3658-lrc-mod169-carry-fourier-block-intertwiner
  - THM-3676-lrc-cross-source-blind-intersection-and-scalar-seam
---

# THM-3679 -- every positive row eventually fails the scalar checksum

**PROVED; PENDING INDEPENDENT HOSTILE AUDIT.**  This theorem lifts the exact
equality graph of THM-3676 from `F_13` to every finite ring
`Z/13^n Z`.  It identifies a finite depth at which the last mod-thirteen
scalar hostile must disappear.  It does not prove that the semantic LRC
current is transported to that depth.

## 1. Higher packet differences

Retain the nine labelled coordinates

```text
U={0,1,2,3,4,5},                 B={1',2',3'},       (1)
```

and let `w in Z_(>0)^9` be the fixed speed row.  The six labels in `U` are
units modulo thirteen and the three blocker speeds are divisible by thirteen,
as in THM-3676.  For `n>=1` put

```text
A_n=Z/13^n Z,
K_n={r in A_n^9:sum_i r_i w_i=0 in A_n}.            (2)
```

Fix a selected source `j in B`, write `{a,b}=B minus {j}`, and for distinct
`k,l in U` define the full depth-`n` packet difference

```text
pi^j_(k,l,n)(r)=(r_a-r_k,r_b-r_l) in A_n^2.         (3)
```

At `n=1`, this is the target-difference map underlying THM-3671/3676.  At
larger `n`, equation (3) is first a relation-address refinement; Section 5
records its exact coordinate-twist realization and the remaining semantic
qualification.

Use the four charts

```text
C={(0,1),(2,3),(4,5),(1,0)}.                       (4)
```

Their equality graph is connected on `U union {a,b}` over every coefficient
ring.  Consequently

```text
pi^j_(k,l,n)(r)=0 for every (k,l) in C

iff

r_i=c_j for every i!=j,                            (5)
```

where the source coordinate `r_j` remains free.

## 2. Two sources force the global scalar diagonal

Let `j,j' in B` be distinct and impose (5) for both source-dependent packet
families.  A unit coordinate belongs to both outside-source sets, so their
two common values agree.  The coordinate `j` belongs to the outside-source
set for `j'`, and conversely for `j'`.  Hence every coordinate agrees:

```text
intersection_(s in {j,j'}) intersection_((k,l) in C)
  ker(pi^s_(k,l,n))

 ={c 1:c in A_n}.                                  (6)
```

This identity is exact over `A_n`; no field rank or division is used.  Adding
the third source family changes neither side.

## 3. The relation equation becomes one checksum

Intersect (6) with `K_n`.  A scalar address `c 1` is a relation exactly when

```text
c W=0 mod 13^n,                 W=sum_i w_i.         (7)
```

Now impose the all-coordinate-unit condition inherited from the all-`91`-unit
address bank:

```text
r_i not congruent 0 mod 13 for every i.             (8)
```

On the scalar diagonal, (8) says that `c` is a unit of `A_n`.  Multiplication
by `c` is then invertible, so (7) is equivalent to

```text
13^n divides W.                                    (9)
```

Therefore the exact jointly blind all-unit fibre is

```text
{c 1:c in A_n^times},       if 13^n divides W,

empty,                      otherwise.             (10)
```

At `n=1`, (10) recovers precisely THM-3676's zero-sum scalar seam.  At
`n=2`, it says that the first divided-difference layer kills that seam unless

```text
169 divides sum_i w_i.                              (11)
```

## 4. Bockstein tower and finite death depth

Suppose a packet difference is zero modulo `13^h`.  Its next divided
difference

```text
(pi^j_(k,l,h+1)(r)/13^h) mod 13                    (12)
```

is well-defined, and vanishes exactly when the packet difference is zero
modulo `13^(h+1)`.  Thus simultaneous vanishing of the ordinary target and
the first `n-1` divided-difference layers is exactly the left side of (6) at
depth `n`.  Equations (9)--(10) are therefore a complete scalar Bockstein
tower, not merely a mod-169 coincidence.

Because every LRC speed is positive,

```text
0<W<infinity.
```

Put `v=nu_13(W)`.  Then the scalar all-unit hostile exists through depth `v`
and is empty at depth

```text
n_*=v+1.                                            (13)
```

No positive speed row can remain jointly scalar-blind at every 13-adic
depth.

## 5. Exact twist realization and semantic boundary

THM-2334 permits reduction of the exact relation-address current modulo every
integer `N`.  At `N=13^n`, push that current through (3).  Fourier inversion
identifies its two-dimensional transform with the coordinate twists

```text
ell=s(e_a-e_k)+t(e_b-e_l),          s,t in A_n,      (14)
```

which translate the nine base coordinates by `ell_i/13^n`.  If the delayed
word clock is `R=13^K` with `K>=n`, its translated coordinate is shifted by

```text
R ell_i/13^n in Z,
```

so every transported word factor is unchanged.  THM-2365 permits the fixed
positive word clock to be chosen arbitrarily large.  Thus each higher packet
difference has an exact physical coordinate-twist profile on any one fixed
owner current; no new carry-convolution identification is needed.

What remains missing is the cross-source semantics.  THM-2305 does not yet
transport literally one nonzero grouped all-unit residue current, in the same
labelled coordinates and phase gauge, through two distinct legitimate
owner-source packet families.  Merely applying the second set of coordinate
twists to the first owner's current does not prove that it is the second
owner-aligned packet current.  The new exact interface is

```text
common all-unit current at depth n_*
  + sourcewise noncancellation
  + two source packet actions
  -> some higher target or divided-difference target is nonzero.   (15)
```

Equation (13) proves that there is no algebraic scalar hostile left at that
depth, while (14) realizes each action separately.  It does not supply the
common-current or noncancellation premises in (15).  No scalar row is excluded
and LRC(14) remains open.  **QED.**
