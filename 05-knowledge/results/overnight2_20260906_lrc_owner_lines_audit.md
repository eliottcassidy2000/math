# Independent audit of complete owner-line count safety

**Verdict: PASS, analytically and on the complete required finite head.**
The [owner-line theorem](overnight2_20260906_lrc_owner_lines.md) correctly
closes every full carrier dictionary that spans the rank-two kernel plane
and is collinear in each owner color. No repair is needed. In fact, its
existing argument gives the strict conclusion

```text
N<2c/11,       S_i<6/77 for all three i.               (1)
```

There is no equality case in this class. This does not assert a universal
count or projection inequality. LRC(14) and physical entry remain open.
The present report is an independent audit artifact; it does not promote
any theorem namespace or edit the producer's claim.

## 1. Fullness, color, and saturation

The filters are exactly those of
[THM-4422 / projection-deficit-and-beatty-row-reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md):
`w=(a,b,c)` is primitive, sorted, distinct, positive, odd, and ternary-unit;
`L={C in Z^3:C dot w=0}`; and the full dictionary consists of all kernel
points with all coordinates nonzero modulo three inside the three strict
roofs `14|C_i|<3(w_j+w_k)`.

The admitted lattice `Gamma` has index three in `L`, its owner map
`phi:Gamma->F_3` has kernel `3L`, and the full dictionary is exactly the
convex open roof region intersected with `Gamma\3L`. These are the
same color/lattice identities used in the audited
[cap proof](overnight_20260906_lrc_cap_carriers.md). Reduction of `L` modulo
three is the full two-dimensional kernel: any modular kernel vector can
be corrected by a multiple of three using an integer Bezout row for `w`.

Write the owner sets as `A,-A`. Rank two implies `|A|>=2`, and their two
affine supporting lines are distinct and avoid the origin. Let `r` be a
primitive vector of `L` along the owner line. If `r` lay in `Gamma`, then
either `phi(r)=0`, contradicting primitivity because `r in 3L`, or a segment
between two same-color points would contain an intermediate point of the
opposite nonzero color on the same line. The latter contradicts the
distinct opposite owner line. Thus `r` is outside `Gamma`.

Consequently `3r` is primitive in `Gamma` and belongs to `3L`. Fullness and
convexity make `A` consecutive with step `3r`. Extend this primitive vector
to a basis of `Gamma`, giving integral transverse height. The owner map
vanishes on the horizontal vector and is nonzero on the other basis vector.
Thus height one is always owner-live.

Suppose the two owner lines had heights `+-k`, `k>=2`. The symmetric hull
of two consecutive points and their negatives is a parallelogram whose
horizontal section at height one has **closed length exactly one** in
units of `3r`. It contains a lattice point. All of this closed hull remains
inside the **strict** roofs, because it is the convex hull of finitely many
points satisfying each roof strictly. The new point has nonzero color and
lies on neither permitted line, contradicting the full dictionary.

Hence `k=1`; the transverse pair is a basis of `Gamma`, and

```text
u cross (3r)=+-3w,       u cross r=+-w.                (2)
```

This audit checks the important type distinction: `r` is an **invisible**
primitive direction, whereas the other generator `u` is live. Replacing
either lattice by an arbitrary selected subset would invalidate the
height-one insertion.

An immediate corollary is that **every raw carrier in this class is
primitive**, and there are exactly `N/2` primitive unoriented origin-rays.
Every point `C=u+3jr` on one line has `C cross r=+-w`. A common divisor
of its coordinates would divide all coordinates of primitive `w`, so is
one. Distinct points on that line cannot be scalar multiples, because
their cross products with `r` are identical and nonzero. This explains
precisely how the class permits arbitrarily many primitive directions.

## 2. Diameter and all seven small directions

For `n=|A|=N/2`, the endpoint displacement is `3(n-1)r`. In a coordinate
with `M=max|r_i|`, the strict diameter bound gives

```text
3(n-1)M<3(w_j+w_k)/7<6c/7,
N<4c/(7M)+2.
```

For `M>=4` this is `N<c/7+2`.

If `M<=3`, invisibility and the mod-three kernel force at least one
coordinate to vanish modulo three. A zero integer coordinate would leave
two same-parity coefficients: the primitive possibilities `(1,1)` and
`(1,3)` force equal speeds and a speed divisible by three, respectively;
`(2,2)` and `(3,3)` are not primitive. If all coordinates are nonzero,
parity and the mod-three relation leave exactly the magnitude pattern
`(1,2,3)`. Two coordinates of magnitude three cannot cancel the remaining
unit coordinate modulo three.

The seven signed/order cases in the producer's table are exhaustive. I
independently intersected each signed-permutation relation with the edges
of the closed ratio triangle `0<=alpha<=beta<=1`, retaining only cases
with an interior point satisfying the strict speed order. This gives
exactly the same seven beta intervals. To check the chord table in a
different affine chart, I used

```text
u=(0,-1/r_1,beta/r_1),       u cross r=(alpha,beta,1),
alpha=-(r_2 beta+r_3)/r_1.
```

The producer uses a chart with third coordinate zero. The two differ by a
multiple of `r`, which shifts every parameter interval together and leaves
`U_i-L_j` unchanged. Direct exact endpoint arithmetic in the new chart
reproduces all seven chosen affine upper bounds:

| r | Closed beta interval | Certificate |
|---|---|---|
| `(1,-3,2)` | `[2/3,1]` | `U_3-L_1=beta/7<=1/7` |
| `(1,3,-2)` | `[1/2,2/3]` | `U_1-L_2=2/21` |
| `(2,-3,1)` | `[1/3,1]` | `U_3-L_1=beta/7<=1/7` |
| `(2,3,-1)` | `[1/5,1/3]` | `U_1-L_2=1/21` |
| `(3,-2,1)` | `[1/2,1]` | `U_3-L_1=2 beta/21<=2/21` |
| `(3,1,-2)` | `[1/2,1]` | `U_1-L_2=2/21` |
| `(3,2,-1)` | `[1/5,1/2]` | `U_1-L_2=1/21` |

Each difference is affine on the whole ratio interval, so its endpoint
maximum is a proved uniform bound. No sampled-speed or piecewise-LP
inference is needed. If the common chord interval is nonempty, its length
`ell` is at most every chosen `U_i-L_j` and hence at most `1/7`.

Normalization is consistent: coordinates have been divided by `c`, so
the same-owner parameter spacing is **`3/c`**, not three. Thus

```text
3(n-1)/c<ell<=1/7,
N<2c/21+2<=c/7+2.
```

The combined uniform bound gives strict count safety for every `c>=53`.
The inherited odd-speed restriction makes `c<=51` the entire complementary
head; no admissible `c=52` is omitted.

## 3. Complete independent finite head

The new verifier enumerates all eligible speed triples with `c<=51` before
any support filter. For each it enumerates `C_2,C_3`, solves `C_1` exactly,
and uses literal strict integer roof comparisons. This differs from the
producer's first/third-coordinate scan solving the middle coordinate.
It then checks the entire owner class for collinearity and plane span.

```text
eligible primitive sorted odd ternary-unit triples: 678
rank-two owner-line dictionaries: 139
cardinalities: 137 of size four, 2 of size six
maximum N/c: 4/25
all maximizers in this finite head: (17,23,25)
```

Every one satisfies `N/c<2/11` strictly. Together with the strict tail,
this proves (1), including the absence of equality, at every height.
Since each raw projection summand is at most `q=3/(7c)`, the three strict
network bounds follow from `S_i<=qN<6/77`. The finite maximum is not
asserted to be the all-height maximum.

Included controls are `(17,23,25)`, `(41,47,49)`, and the complete flat
ten- and thirty-carrier dictionaries `(1,137,277)`, `(1,499,1001)`.
The noncollinear owner triangles at `(19,23,29)` and `(23,29,37)` are
explicitly excluded. Every selected row also checks (2), the invisible
primitive direction, all consecutive `3r` steps, and the strict diameter
or small-direction bound. All 3,880 explicit gates pass; normal and
optimized Python streams byte-match.

## 4. Exact connection and remaining class

The earlier
[midpoint-deficit result](overnight_20260906_lrc_midpoint_deficit.md)
exhibited affine deficits with zero second difference. This theorem closes
that whole owner-line geometry by recovering its minimum transverse lattice
height, not by claiming the missing baseline was curvature. Its overlap
with the two-ray class is exactly `N=4`, and its overlap with an exact
three-ray class is exactly `N=6`; larger owner-line dictionaries are a
distinct unbounded-ray sector.

Incoming
[THM-4431 / colored-lattice-basis-and-three-direction-lrc-network-closure](../../01-canon/theorems/THM-4431-colored-lattice-basis-and-three-direction-lrc-network-closure.md)
was inspected as **RESERVED / UNPROVED EMPTY STUB**. Its announced
colored-basis/three-direction work is a related research route, not a
dependency or an audited three-ray result in this report. The present
proof uses the invisible primitive `r` and allows arbitrarily many rays.

Combining the audited class closures leaves a possible network failure
with genuinely noncollinear points in each owner color, an affine collinear
triple somewhere in the full dictionary, at least three primitive
directions, and a positive count surplus unpaid by the proved midpoint
certificate. That is a remaining class, not an asserted counterexample.

The connection preserves every point by the map
`C=u+3jr -> (u,r,j)` and pays the raw projection count through its exact
interval. It loses no carrier. Fullness supplies the height-one test,
owner color supplies the step three, and the determinant retains the line
offset. Selected subsets, origin-ray counts alone, and curvature without
baseline do not preserve these premises.

## Reproduction

```text
python 04-computation/overnight2_20260906_lrc_owner_lines_audit.py
python -O 04-computation/overnight2_20260906_lrc_owner_lines_audit.py
```

Companions:
[independent verifier](../../04-computation/overnight2_20260906_lrc_owner_lines_audit.py),
[frozen output](overnight2_20260906_lrc_owner_lines_audit.out).
The all-height lattice proof is analytical; the finite replay supplies
exactly the complementary head, and the seven-case chord identities are
symbolic rational identities on complete ratio intervals.

```text
SHA-256, LF bytes:
source bd93b68bd6351485d1e2a75548108dbb0ffe451213e9cecc597f1bae0e203296
output 213f988777d45b2d767bb7caedfb93e0d4968eb7dcd53993fbd7b7f579c0a8e4
```
