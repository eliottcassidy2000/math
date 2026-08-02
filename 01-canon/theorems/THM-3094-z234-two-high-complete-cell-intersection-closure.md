---
id: THM-3094
title: "z234 two-high complete-cell intersection closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Exact whole-cell
  cardinalities close all 108,966,498 finite literal candidates in
  THM-3087's two scalar-unbounded bodies.  The other two THM-3078 rows
  remain.  No projected-cap, ledger, or LRC claim is made.
source: codex-lrc-z234-boundary-scout-2026-08-02
audit: >
  Two independent hostile audits rebuilt both literal body carriers and all
  64,740 label/body columns without importing the claimed extrema.  They
  reproduced the unique minima, second minima, strict pigeonhole slacks,
  direct extremal intersections, 223/431 mask bank, and the exact split of
  108,966,498 candidates; they accepted the strict-open endpoint,
  distinct-label, full-cell projection, and aligned-cap typing.  Root replay
  caught and repaired a stale declared script hash and an ill-typed low-label
  referent.  Fresh normal and optimized transcripts match stored output.
depends_on:
  - THM-3087-z234-two-high-common-torsion-and-finite-height-reduction
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
related:
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
script: 04-computation/lrc14_j7_k3_z234_two_high_complete_cell_intersection_closure_thm3094.py
output: 05-knowledge/results/lrc14_j7_k3_z234_two_high_complete_cell_intersection_closure_thm3094.out
script_sha256: c456e28c88cac8d0b687c1636463aa79642b2fb39ea187eafc836f46e3ab3e0e
output_sha256: 6752dbac509aba6c739eb9ce3cd59da85d39ae47b6f4f60b1657ecaf4cf36a95
semantic_sha256: 7c3fa428224d59b717770b936bb2d35aea2642dec83c77df8395c538a8015682
hash_basis: LF-normalized bytes
---

# THM-3094 -- z234 two-high complete-cell intersection closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

[THM-3087](THM-3087-z234-two-high-common-torsion-and-finite-height-reduction.md)
leaves `108,966,498` ordered literal candidates in the two scalar-unbounded
THM-3078 bodies.  Every such candidate is impossible.  Therefore the rows

```text
(1,5,9,11,12,14),
(1,9,10,11,12,14)                                      (1)
```

are empty in the pinned THM-2941 projected `k=3`, `z_1=234` necessary atlas.
The exact residual of the four-row THM-3078 boundary is reduced to

```text
(1,5,6,9,12,14),
(2,5,9,11,12,14).                                      (2)
```

This theorem does not close either row in `(2)`, lower the projected cap,
change the canonical ledger, or prove LRC(14).

## 2. Literal complete-cell sets

Fix one body `E` in `(1)`.  Its body period is `L=194040`, its pinned first
drift is `234`, and its fixed low drift is respectively

```text
243 or 260.                                             (3)
```

Let `J` be the set of complete body `1/L` cells which are safe from both
`234` and the respective label in `(3)`.  For a high label `z`, define

```text
A_z={j in J : [j/L,(j+1)/L] is disjoint from D_z}.       (4)
```

Because `D_z` is strict-open, membership in `(4)` is exactly the pair of
weak integer inequalities

```text
r=(jz mod L),       14r>=L,       14(r+z)<=13L.          (5)
```

No centered or translated-band approximation enters `(5)`.

THM-3087 proves that every remaining actual pair consists of two distinct
integer labels

```text
19111<=z,w<=51480.                                     (6)
```

The exact all-label evaluation of `(5)` gives

```text
E=(1,5,9,11,12,14), low=243:
  |J|=44062,
  unique minimum |A_48510|=22020,
  second minimum |A_48512|=25066;

E=(1,9,10,11,12,14), low=260:
  |J|=44616,
  unique minimum |A_48510|=22308,
  second minimum |A_48512|=24816.                       (7)
```

The second entry in each row of `(7)` means the least cardinality after the
unique minimizer, not merely the value at one chosen control.

## 3. Pigeonhole closure and the lossless projection

For any two distinct labels in `(6)`, `(7)` gives

```text
|A_z intersect A_w|
 >= |A_z|+|A_w|-|J|
 >= 22020+25066-44062 = 3024
```

in the first body, and

```text
|A_z intersect A_w|
 >= 22308+24816-44616 = 2508                         (8)
```

in the second.  Choose a common cell.  The complete interval over that cell
is safe from the body, `234`, the fixed low label, and both high labels.
Consequently it lies in `S_(E,Z)`.  The map

```text
phi_L(t)=Lt mod 1                                      (9)
```

maps any complete `1/L` cell surjectively onto the circle, so

```text
P_(E,Z)=phi_L(S_(E,Z))=T.                              (10)
```

But THM-2941 `(25g)` says that a literal completion would force
`P_(E,Z) subset U_A`, while THM-1166 bounds the inherited three-aligned union
by

```text
mu(U_A)<=36/91<1.                                      (11)
```

Equations `(10)--(11)` are contradictory.  This closes every pair in `(6)`
and hence all `108,966,498` candidates.

## 4. Hostile controls and first failed extension

The actual minimum/second-minimum intersections in the two bodies contain

```text
11798 and 11792 cells,                                 (12)
```

with common witness cell `14850`; the cardinality bounds `(8)` are therefore
conservative.  The proof nevertheless uses only the universal lower bounds.

Distinctness is load-bearing for that universal argument: duplicating the
unique minimum would make the first displayed cardinality sum nonpositive.
Actual LRC speeds are distinct, and THM-3087 explicitly removes equality in
the same-denominator case.  This caveat prevents the closure from being
misread as a statement about arbitrary repeated-label multisets.

MISTAKE-334's translated-band hostile is also inapplicable.  It refutes use
of a centered lattice capacity after an unknown translation; `(5)` instead
tests the literal label on the literal body cells.

The first failed extension is

```text
the THM-3087 height contract for its two scalar-unbounded bodies
  ==> the same label box for the other two THM-3078 rows.                (13)
```

Those rows have different scalar thresholds and are not inspected here.

## 5. Exact evidence

Run

```text
python 04-computation/lrc14_j7_k3_z234_two_high_complete_cell_intersection_closure_thm3094.py
python -O 04-computation/lrc14_j7_k3_z234_two_high_complete_cell_intersection_closure_thm3094.py
```

The companion pins THM-3087's source, output, and `1699`-mask bank.  It
rebuilds both literal body carriers, evaluates every one of the `64,740`
label/body columns in `(6)` by exact `int64` arithmetic, checks the two
extremal columns again by a scalar implementation, freezes both complete
count-vector hashes, and exhibits the intersections and witness in `(12)`.
Normal and optimized transcripts byte-match the stored output.

The count-vector hashes are

```text
9149b88558e7220fcf60aa1d265470cf7f1e34808b99efc1b00a48ff4a8b7c96,
190ecd77e97f4a59d83d9052a13559a124ea7991423e06fa972ae1475f3ef509. (14)
```

The independent hostile audit rebuilt the body-safe ranges directly from
the six combs and reproduced `(7)--(12)` without importing the claimed
extrema.

## 6. Scope

This theorem empties exactly the two rows in `(1)` inside the pinned
projected necessary atlas.  The two rows `(2)`, the projected cap `234`, the
canonical ledger, physical reconstruction outside this projected implication,
`k<=1`, and LRC(14) remain open.

QED.
