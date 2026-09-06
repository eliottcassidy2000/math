# Independent referee: the universal degree-zero LRC14 network certificate

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS, 2026-09-06.**
Reviewed source surface: `origin/main` at
`2e1ef24c4aa31ebde37de46f2d14bf16f4ed9620`. The native height-601 code and
its frozen output were read in full; this referee did **not** rerun the
1,317,935-row native census. Its completed finite output is an explicit
dependency, not newly generated independent evidence from this referee.

At the reviewed commit THM-4434 is still a **RESERVED / UNPROVED EMPTY
STUB**. This audit supports promotion of the assembled proof candidate; it
does not silently treat the reservation itself as proved canon.

## Exact audited consequence

For every primitive sorted distinct positive odd ternary-unit triple
`w=(a,b,c)`, the complete degree-zero raw projection networks of THM-4414
satisfy

```text
min_i E_i(w)<=6/77,
equality iff w=(1,5,11).
```

If the triple has no signed full-support relation of coefficient-magnitude
pattern `(1,1,2)`, then **each** projection is strictly below `6/77`.
The norm-four exception cannot be deleted: `(1,5,7)` has
`E=(8/245,6/49,4/35)`. Entry, synchronization, arbitrary body Haar floors,
and LRC(14) are not consequences of this local certificate alone.

## Sources reviewed and analytic checks

The proof sources are the
[full-support affine-slice theorem](lrc14_empty_core_certificate_sep06.md),
[support-two addendum](lrc14_pair_relation_empty_core_certificate_sep06.md),
[global slope and height reduction](lrc14_global_slope_empty_core_certificate_sep06.md),
[complete coefficient-box report](lrc14_coefficient_box_empty_core_three_ray_sep06.md),
and [native finite-head report](lrc14_universal_literal_empty_core_sep06.md).
The coefficient compiler and native C++ verifier were also read in full.

The affine address retains all defects, raw multiples, and owner residues.
The error cube really maps onto the strict raw carrier body: the three
phase intervals have a common interior point precisely when their pairwise
overlaps are strict, by the interval Helly property. Primitivity makes the
fixed-defect carrier lattice an integer torsor parallel to the relation.
The owner gate leaves two residues modulo three for a ternary-unit relation,
and one for a relation with one zero residue. Actual zero integer
coefficients are handled separately rather than divided through.

The planar zonotope map has generator-pair determinants `+/-a,+/-b,+/-c`.
Its area and section-width integral are exactly `9(a+b+c)/49`. Evenness and
concavity make the section width nonincreasing on the positive half-line.
Setting endpoint widths to zero correctly enforces the strict defect list,
including the possible vertical-edge boundary for actual-zero coefficients.
The symmetric rectangle comparison has error at most the central width;
both owner-sampling formulas consequently give

```text
F_v(w)/c <=6/49+4/(7M),       M=max_i|v_i|.
```

For `M>=19`, this is `142/931=15/98-1/1862`. The finite complement is
exactly the coefficient box `M<=18`, not a bounded-speed or bounded-norm
extrapolation. Its compiler enumerates all 308 primitive magnitude patterns
of even norm and support at least two, with at most one zero residue, after
excluding impossible `(0,1,1)` and inherited norm-four `(1,1,2)`.
Its rational clipping, isolated-sign representatives, coordinate
permutations, complete defect lists, and convex speed-vertex maximization
are correct. The declared cube-edge comparison supplies a second geometric
path for full support, and the rectangle formula checks support two.

One non-load-bearing wording correction belongs in the owner note:
coefficients at most eighteen imply norm **at most** 54, not that a
primitive pattern can attain 54. In this filtered primitive box the largest
norm is 52, for example `(17,17,18)`. The enumerated box and defect range
are already correct; no missing case or false bound results from the wording.

The strict intercept estimate is
`N<15c/98+2S/7+4/3`, and solving the count gate gives exactly
`c>=(308/31)S+4312/93`. The projected integer relation lattice has
determinant `c`; the explicit projected `l1` ball has area strictly larger
than `3L^2/4`. The fundamental-domain argument therefore gives a primitive
relation with `S<4sqrt(c/3)`. The norm is even because all speeds are odd.
The cases `S<=56` and `S>=58` prove the stated tail for `c>=603`:
`56056/93<603` and `g(58)=3023/372>0`, with `g` increasing thereafter.
No odd-norm gap is overlooked. The remaining eligible universe is precisely
`c<=601`; 602 is even and 603 is a nonunit.

## Native finite head: independence and implementation audit

The source
[lrc14_universal_literal_empty_core_sep06.cpp](../../04-computation/lrc14_universal_literal_empty_core_sep06.cpp)
does not import a carrier congruence, direction classifier, or raw roof
formula. It independently scans all six native sheet assignments using
exact interval endpoints on denominator `42abc`.

The residue cursor generates every interval in sorted order. Sheet zero
has the two clipped seam pieces; the other sheets lie strictly inside the
cut circle. At each positive triple contact, the contribution to projection
`i` is the full other-pair intersection length capped by the full omitted
interval length, exactly the required edge-minimum consumer. Advancing all
cursors tied at the least right endpoint visits each contact once and
excludes zero-length contacts. These operations are correct for speed one,
seam splits, and simultaneous endpoints. The integer products fit the
declared 64-bit counters on the stated head; cross-multiplied comparisons
use 128-bit integers.

Every eligible row is evaluated **before** norm-four classification. The
classification identities `c=2a+b`, `c=a+2b`, `2b=a+c` exhaust the signed
norm-four pattern in the sorted positive domain. The saved native output
reports 1,317,935 rows, unique minimum equality `(1,5,11)`, and strict
every-projection bounds on all 1,308,734 non-norm-four rows. Its eight named
controls include the norm-four hostile and wide cases outside the head.
The separately recorded H79 replay uses both literal Python and raw
implementations, while the sanitizer replay checks native low-level safety.
The native-contact count is correctly kept distinct from the raw count.

Thus the completed native artifact discharges the finite head left pending
in the earlier prose of the global note. The remaining local-proof work at
the reviewed surface is honest status promotion and dependency integration,
not an unperformed mathematical tail or head computation. The noncryptographic
row checksum is correctly labelled and is not substituted for the frozen
SHA-256 artifact hashes.

## Independent rediscovery and the changed frontier

During this audit the present referee independently recovered the full-
support `(3,2,1)` near-corner triangle, its two owner-affine lines, the
exact chord width, and the threshold `c>=25`. `observer_collision` then
independently passed that derivation. The newly read full-support source
already contains **the same mechanism in Section 5**. This is corroborating
rediscovery, not a new theorem or an additional claimed family. No duplicate
package or scarce ID is created.

The [body-specific Haar consumer](lrc14_universal_haar_consumer_empty_core_certificate_sep06.md)
has a sound compact/open equality argument and common-divisor normalization.
It needs the **actual** body `G_C`, either its measure or its componentwise
intersection with the tail-safe lifts, and an entry map that really produces
the required odd ternary-unit triple. Those native objects are the next
research obligation. Repeating the universal triple-network census is no
longer the appropriate anchor after the local theorem is promoted.
