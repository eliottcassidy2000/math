# Classifying the exterior routes of the smallest interleaved rectangles

The [original interleaving](../chanpak_mixed_defects_20261115/INTERLEAVED.md)
gave 721,723,725,727 with two disjoint, additive switches. Background
research found all such small rectangles. We independently replayed the
reduction, and classify *all* fixed-exterior calibrations of those cores.

## Six cores and a cheap continuation filter

At width 20, bucket all words by exact H and look from each base in the
H+2 and H+4 buckets. The XOR differences must be disjoint, and their
joint XOR must land at H+6. Trim only common leading zero and trailing
one pads from all four settings. All shorter words embed in this width,
so the result covers every width through 20. The exact results are
56 marked rectangles at width 20 and six reduced ones; no disjoint
additive +2,+4,+8 cube satisfies the extra four exact joint tests.
This is a finite statement, not an all-width exclusion.

| Width | Base (binary integer) | +2 mask | +4 mask | Base H |
|---:|---:|---:|---:|---:|
| 11 | 1158 | 32 | 220 | 721 |
| 11 | 1270 | 32 | 472 | 721 |
| 11 | 1244 | 478 | 32 | 959 |
| 11 | 1574 | 988 | 32 | 959 |
| 13 | 4614 | 956 | 64 | 2881 |
| 13 | 5110 | 1976 | 64 | 2881 |

For a fixed-width core switch let δ_i be its signed bit change, indexed
from the least significant bit. Let U=Σ 2^i δ_i, V=−Σ 2^(w−1−i)δ_i.
They are the changes in the ordinary and reverse-complement addresses.
If the endpoint bits are fixed, the original step changes by aV+bU
for ordinary left address a and reverse-complement right address b;
the seam parity terms are unchanged. Both numbers U,V are nonzero.
Their valuations follow from the unique lowest occupied power:
v₂(U)=min supp δ, and v₂(V)=w−1−max supp δ. Hence preservation in any
nonzero common exterior requires all switch supports to have the same
midpoint (and their full positive ratio −V/U to agree). The two
midpoints differ in each of the six rectangles. A third switch
entirely in the exterior cannot extend these cores while preserving
both original steps in every setting. A switch altering the core is
not ruled out.

## A moving modulus can make a physical four-count route

A fixed exterior does not change the **mixed** defect, so let P,Q be
the two signed steps after padding. Four distinct values of this
additive square are equally spaced exactly when either Q=±2P
or P=±2Q: translate by its minimum and take absolute step sizes.
The step along their increasing order is the smaller of |P|,|Q|.
Successive tournaments differ by the smaller control, the union of
both controls, and the smaller control again. This is a physical
three-edge route of mark changes, not a claim that each edge is a
single arc or a single control.

It suffices by reverse-complement symmetry to use the cores
1158,1244,4614; their mates interchange a,b. Here is the complete
nonnegative exterior solution (integer t≥0):

| Core | Signed steps (P,Q) | Exterior (a,b) | Modulus M |
|---:|---|---|---|
| 1158 | (2−32a+32b, 4+56a−44b) | (9t,10t) | 2+32t |
| 1158 | same | (5t+1,2t) | 30+96t |
| 1244 | (2−44a+38b, 4−32a+32b) | (11t,14t) | 2+48t |
| 1244 | same | (12+13t,9+10t) | 92+96t |
| 4614 | (2+88a−76b, 4−64a+64b) | (9t,10t) | 2+32t |
| 4614 | same | (7+11t,9+14t) | 66+96t |

To derive the table, substitute the signed steps into the four linear
equations Q=±2P, P=±2Q. Each omitted equation has its coefficient
gcd not dividing the constant. For each surviving equation, the
displayed (a,b) are one particular solution plus its primitive
kernel vector; nonnegativity is exactly t≥0. Each core begins
with 1 and ends with 0 in all settings. The same is true after
any nonzero exterior addresses and neutral trimming, so all
four witnesses are strong.

This produces infinitely many exact routes. It also gives a
sharper mismatch with the old target clock: **none of these six
cores, in any fixed exterior, realizes a four-term equally spaced
route of step 3·2^k for any k≥1**. The cases with valuation one
have M≠6; the only valuation-two family is 92+96t≡28 mod32 and
never equals 12. All these moduli and the exact placement of
their targets are tied to the same exterior line. This does not
rule out joining a pair of witnesses by a different core or
changing the reader.

Replay: compile and run `rectangles20.cpp` and run
`python verify_calibration.py`. The latter independently counts
Hamilton paths on the reduced cores and tests the exact affine
calibration by direct integer evaluations.
