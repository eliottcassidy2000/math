---
id: THM-2970
title: "Located low-cell torsion tail closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The exact projected k=2 atlas has nineteen forced-high rows on
  1680<=z1<=1742.  Exact unit-ray envelopes force exactly one high later
  label; sixteen rows are scalar-empty and every one of the thirteen
  remaining literal classes is killed by a whole-cell order-two (also
  order-seven) torsion pair.  The theorem-local residual has thirty-nine
  all-high rows; the canonical THM-2941 exact-descent lineage independently
  closes all fifteen atlas rows, including those nine rows at height 1736.
source: codex-lrc14-located-low-cell-torsion-tail-2026-07-30
audit: >
  An independent hostile auditor replayed the settled normal and optimized
  scripts byte-for-byte against the stored transcript and independently
  checked 86,664 singleton values, 5,871 unit-ray recurrences, all thirteen
  literal classes, 52 order-two/order-seven lifts, and 4,872 unit characters.
  The universal two-high invoice, strict margins, LF hashes, q=7 seam, and
  q=8 hostile boundary all passed.
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
script: 04-computation/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.py
output: 05-knowledge/results/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.out
script_sha256: eeba4fa6968365269b776832cbcbd3753affcac10915726eb28fb7a817c9eec4
output_sha256: 8f1dd4702b64caef2b94cfdf8e474bba137900ef48200bf5d9aabb2e440e454f
hash_basis: LF-normalized bytes
---

# THM-2970 -- projected `k=2` forced-high exact-ray / located-torsion closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This note is an independently reproducible package for the
nineteen `HIGH-TAIL` rows in the projected `k=2`, `1680<=z_1<=1742` scalar
atlas.  It does not touch the other thirty-nine rows.

## Immutable inputs and evidence

The exact referee is

```text
04-computation/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.py
SHA256 eeba4fa6968365269b776832cbcbd3753affcac10915726eb28fb7a817c9eec4
```

and its stored transcript is

```text
05-knowledge/results/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.out
SHA256 8f1dd4702b64caef2b94cfdf8e474bba137900ef48200bf5d9aabb2e440e454f
```

The script pins and verifies these two inputs before doing mathematics:

```text
04-computation/lrc14_j7_k2_scalar_band_1680_1742_thm2941.py
  1224d5594571f21c91c55fe3ab165c4fc34ba7968719862d12660d24efac919d
05-knowledge/results/lrc14_j7_k2_scalar_band_1680_1742_thm2941.out
  c20607cb478ed491d7000f2b8a49213f57d1606a5152700ac3b50c836e2dc66c
```

Reproduction:

```powershell
python 04-computation/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.py --output .scratch/thm2970.normal.out
python -O 04-computation/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.py --output .scratch/thm2970.opt.out
```

Normal, optimized, and stored transcripts agree byte for byte.  The semantic
digest of all thirteen surviving literal upper classes, including their
located torsion witnesses, is

```text
1642cd1b5f7bcca7c79abf37850e7eee3c7b4aede31a17dada82dc1cac6dc427.
```

## Candidate theorem

Take any of the nineteen `HIGH-TAIL` scalar rows in the exact all-body atlas
on `1680<=z_1<=1742`.  Under the inherited projected high-wall condition,
there is no literal packet of four later distinct nonaligned labels satisfying
the projected `k=2` obligation `sum_i Delta(z_i)>=h/91`.

More precisely:

1. the exact unit-ray envelope excludes two or more later high labels on all
   nineteen rows;
2. the projected high wall supplies at least one later high label, so exactly
   one is high and the other three lie in the finite interval
   `z_1<z<floor(13L/150)+1`;
3. exact unit-ray maxima make sixteen rows scalar-empty;
4. the remaining three rows contain exactly thirteen
   `(high denominator, finite low triple)` upper classes;
5. every one of those thirteen classes has two complete fixed-safe `L`-cells
   separated by order two, and also by order seven, in the high denominator;
   the located unit-ray torsion lemma makes the projected safe residual have
   full mass one.

Hence all nineteen forced-high rows are empty, uniformly over the unbounded
height on every unit ray.  No label horizon is used as an exhaustion.

## Proof

### 1. Exact ray law must come first

Let `C_E` be the body carrier, let `L=14 lcm(E)`, and write a nonaligned label
on a denominator-`d` unit ray as

```text
z=(L/d)u+mL,   gcd(u,d)=1.                              (1)
```

Direct periodic integration gives

```text
z Delta_E(z)=A_E((L/d)u).                               (2)
```

Thus a positive ray decreases from its first point at or above the wall, a
zero ray stays zero, and a negative ray is bounded above by the unattained
supremum zero.  The referee reconstructs the carrier from the defining strict
danger combs on the master ruler

```text
14 lcm(1,...,14)=5,045,040,
```

scans every nonzero residue modulo `L`, and obtains the exact maximum over all
unbounded high labels.

This order is load-bearing.  The atlas's analytic `HIGH-TAIL` number is a
slotwise upper bound, not an attained ray.  If `H_E` is that analytic bound,
the tempting duplicate-permitting estimate

```text
Delta(z_1)+2H_E+(two best unrestricted lows)             (3)
```

closes only seventeen of the nineteen rows.  It fails on exactly

```text
z_1=1722, E=(2,8,9,10,12,14),
z_1=1736, E=(1,4,10,12,13,14).                           (4)
```

Replacing `H_E` by the exact all-residue unit-ray maximum `H`, and writing
the two best finite low values as `l_1>=l_2`, the universal upper bound for
every packet with two, three, or four highs is

```text
2H+max(H,l_1)+max(H,l_2).                               (5)
```

Indeed its three specializations dominate respectively
`2H+l_1+l_2`, `3H+l_1`, and `4H`.  In this exact nineteen-row universe the
referee also verifies

```text
l_2-H >= 3131/11953620 > 0,
```

so `(5)` agrees with the simpler `2H+l_1+l_2` invoice here.  It is strictly
smaller than `h/91` on all nineteen rows.  The minimum exact gap is

```text
121670443/941781292179                                  (5a)
```

at the first row in (4).  Thus `(5)` excludes every packet with at least two
high labels.

The inherited high wall applies because all nineteen rows have
`z_1<floor(13L/150)+1`; it forces at least one of the four later labels above
that wall.  Therefore there is exactly one high label and exactly three finite
low labels.

### 2. Exact finite reduction

For each row and each exact high denominator `d`, scan every distinct low
triple whose exact singleton excess plus the attained high-ray maximum reaches
`h/91`.  Sixteen rows have no such triple.  Their closest strict scalar-empty
gap is

```text
1323977953/569777681768295,
```

at

```text
z_1=1694, E=(2,8,9,10,12,14).
```

Only three rows survive this exact scalar step:

| `z_1` | body `E` | exact upper classes |
|---:|---|---:|
| 1708 | `(2,8,9,10,12,14)` | 2 |
| 1722 | `(2,8,9,10,12,14)` | 10 |
| 1736 | `(1,8,10,12,13,14)` | 1 |

There are nine distinct row/low packets and thirteen denominator-decorated
classes.  They are, verbatim:

| # | `z_1` | `E` | `d` | finite lows | high representative | unit `u` |
|---:|---:|---|---:|---|---:|---:|
| 1 | 1708 | `2,8,9,10,12,14` | 588 | `1722,1790,1810` | 3180 | 53 |
| 2 | 1708 | `2,8,9,10,12,14` | 588 | `1790,1810,1946` | 3180 | 53 |
| 3 | 1722 | `2,8,9,10,12,14` | 504 | `1790,1810,1946` | 3710 | 53 |
| 4 | 1722 | `2,8,9,10,12,14` | 504 | `1790,1810,1836` | 3710 | 53 |
| 5 | 1722 | `2,8,9,10,12,14` | 504 | `1790,1810,2142` | 3710 | 53 |
| 6 | 1722 | `2,8,9,10,12,14` | 588 | `1790,1810,1946` | 3180 | 53 |
| 7 | 1722 | `2,8,9,10,12,14` | 588 | `1790,1810,1836` | 3180 | 53 |
| 8 | 1722 | `2,8,9,10,12,14` | 588 | `1790,1810,2142` | 3180 | 53 |
| 9 | 1722 | `2,8,9,10,12,14` | 588 | `1810,1836,1946` | 3180 | 53 |
| 10 | 1722 | `2,8,9,10,12,14` | 588 | `1810,1946,2142` | 3180 | 53 |
| 11 | 1722 | `2,8,9,10,12,14` | 588 | `1810,1836,2142` | 3180 | 53 |
| 12 | 1722 | `2,8,9,10,12,14` | 2520 | `1790,1810,1946` | 3514 | 251 |
| 13 | 1736 | `1,8,10,12,13,14` | 196 | `1836,2004,2340` | 13260 | 17 |

The representative is the attained maximum on its denominator class.  The
located-cell conclusion is uniform in every other unit and in every later
height `m` on the ray.

### 3. Located torsion cells

Fix one class, and call the packet consisting of `z_1` and its three finite
lows the fixed packet.  An `L`-cell `j` is retained only when it is a complete
body-carrier cell and, for every fixed label `a`,

```text
14(aj mod L)>=L,
14((aj mod L)+a)<=13L.                                  (6)
```

The weak inequalities in (6) are intentional: danger arcs are strict-open.
They say that the whole cell, not merely one root or one sample point, is safe
from every fixed label.

For each class the exact census finds retained cells `j,k` such that

```text
k-j=d/q mod d                                            (7)
```

for the following numbers of classes:

```text
q=2: 13,  q=3: 12,  q=4: 13,
q=5:  1,  q=6: 12,  q=7: 13.                            (8)
```

In particular order two alone closes all thirteen, and order seven also
closes all thirteen.  From (1) and (7), the two high-label arguments differ by

```text
z(k-j)/L=u/q mod 1.                                      (9)
```

Because `gcd(u,q)=1`, their circular separation is at least `1/q>=1/7`.
Two open danger arcs of radius `1/14` cannot overlap.  At every projected
phase at least one of the two complete cell lifts is therefore safe from both
the fixed packet and the high label.  The projected safe residual has full
mass one, contradicting the inherited completion cap (indeed, any cap
strictly below one).

### 4. Sharp endpoint and hostile order

The order-seven endpoint is real.  Take

```text
L=14, d=7, q=7, u=1, m=1, z=16, j=0, k=1, theta=13/16.
```

The two high arguments have circular norm exactly `1/14`.  Their closed
danger arcs meet, but both strict-open dangers omit the seam.  Any executable
generalization must therefore use `>=1/7`, not `>1/7`, together with the
strict-open typing.

Order eight is the first failure.  Take

```text
L=56, d=8, q=8, u=1, m=1, z=63, j=0, k=1, theta=5/6.
```

The arguments are `15/16` and `1/16`; both have norm `1/16<1/14`, so both
danger indicators fire.  The two arcs overlap in length

```text
1/7-1/8=1/56.
```

This is a genuine pointwise hostile, not merely a comparison of abstract arc
lengths.

## Exact ledger and boundary

Before this reduction the scalar height histogram is

```text
1683:1, 1694:10, 1702:3, 1708:14,
1722:11, 1724:2, 1732:2, 1736:15.                       (10)
```

Removing the nineteen forced-high rows leaves thirty-nine rows:

```text
1694:7, 1702:3, 1708:9, 1722:7,
1724:2, 1732:2, 1736:9.                                 (11)
```

Within the relative residual of this reduction, the top height is `1736`, with
nine rows, and the lowest occupied height is `1694`.  The canonical THM-2941
exact-descent lineage independently closes the entire `1736` shell; later
all-body descents continue to the current global projected `k=2` cap `1656`.
Neither later step alters the exact census here.
Of the thirty-nine rows, thirty-eight
have `L=11760` and high floor `1020`; one, namely

```text
z_1=1722, E=(1,6,9,10,12,14),
```

has `L=17640` and high floor `1529`.  Every suffix certificate is four
literal `EXACT` values.  More importantly, every one of these rows satisfies
`z_1>=high_floor`.  Since later labels are ordered strictly above `z_1`, all
four later labels are high.  The first implication of the present proof—one
high plus three finite lows—therefore fails structurally, not numerically.

This theorem does **not** close its thirty-nine-row relative residual, does not
change the LRC ledger by itself, and does not assert that a pointwise or
root-only safe pair suffices.  Its exact scope is the nineteen `HIGH-TAIL` rows
and complete whole-cell fixed-packet safety.  Current navigation must compose
it with, rather than overwrite, THM-2941's canonical exact-descent and
lower-atlas closures.  The earlier `z_1=1736` hybrid is superseded evidence
under MISTAKE-331/332 because its digest included a solver-selected basis.

## Correction boundary

The atlas transcript writes `gap=upper-lower`.  During hostile audit it was
briefly read with the reverse sign, which falsely made the crude analytic
two-high estimate appear positive on all nineteen rows.  Equations (3)--(5)
are the repair: the crude estimate has the two honest hostiles (4), while the
exact all-residue unit-ray estimate closes both.  No statement in this note
uses the reversed-gap calculation.

**QED (finite exact forced-high scope).**
