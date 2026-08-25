---
id: THM-4030
title: "LRC(14) d=4 affine defect-lattice boundary"
status: >
  PROVED RELATIVE TO THM-4024/4004 + VERIFIED-EXACT + INDEPENDENTLY
  REFEREED. At the remaining d=4,c_4=8 equality boundary, full spoilage of
  the four labelled lifts is equivalent to a bounded support-three odd defect
  relation with exact gcd and strict-window sidecars. Certificate-negative
  rows close. Certificate-positive rows remain unresolved because the test
  forgets the divided pack's safe-phase set. This does not prove LRC(14).
source: root + prime_power_audit + d4_referee / sequence continuation, 2026-08-24
audit: >
  PASS. Direct rational wall-cell subdivision, bounded affine-centre
  enumeration, and bounded defect enumeration agree on 270 profiles. The
  audit checks 4,153 semantic gates, exact endpoint controls, the minimal
  distinct-exception hostile, and a typed ten-pack hostile. Normal and
  optimized outputs are byte-identical.
depends_on:
  - THM-4024-lrc14-complete-divisor-incidence-envelope
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
related:
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
script: 04-computation/lrc14_d4_affine_defect_lattice_boundary_thm4030.py
output: 05-knowledge/results/lrc14_d4_affine_defect_lattice_boundary_thm4030.out
script_sha256: dddd440b92b92c0826348515f9dd5bf54985ea6da456b0beda56c23cc1b36385
output_sha256: baafad92df16b304eae64e8d6e95ee1c3c44e3cbd9fd3812d2836e83dd0a89d4
hash_basis: raw LF bytes
---

# THM-4030 -- the d=4 affine defect-lattice boundary

**PROVED RELATIVE TO THM-4024/4004 + VERIFIED-EXACT + INDEPENDENTLY
REFEREED.** The theorem converts the continuous selector obstruction at the
last composite equality modulus into a finite exact Diophantine test for
each exception triple. It is necessary for a counterexample, not sufficient.
LRC(14) remains open.

## 1. Typed boundary and inheritance

Retain THM-4024's `d=4,c_4=8` boundary in the rank-eleven `11+2` branch.
The eight body speeds divisible by four and the divided pair give a cited
ten-speed pack. The three body exceptions must be

```text
e=2r, a, b,                  r,a,b positive and odd,   (1)
```

with exactly one exception of 2-adic valuation one. All thirteen LRC speeds
are pairwise distinct.

For a pack phase `y in R/Z`, its four lifts are

```text
x_j=(y+j)/4,                 j=0,1,2,3.                (2)
```

Call `y` fully spoiled by the exceptions when every label `j` satisfies at
least one of

```text
||e x_j||<1/14,  ||a x_j||<1/14,  ||b x_j||<1/14.     (3)
```

The closest proved mechanism is THM-4024's exact bad-count equality
`2+1+1=4`. Its canonical hostile is the typed `(2,9,11)` selector. The
least-used sidecar is not another count but the affine centre lattice of the
three open danger intervals.

## 2. Exact equivalence

The following are equivalent.

1. Some phase `y` is fully spoiled.
2. For one of the two orientations `(a,b)` or `(b,a)`, there are integers
   `A,B,C` such that, with

   ```text
   N_ab=2bA-a(2B-b),
   N_ar=4rA-a(2C-r),
   N_br=4rB-2bC-br,                                  (4)
   ```

   one has

   ```text
   7|N_ab|<a+b,
   7|N_ar|<a+2r,
   7|N_br|<b+2r.                                      (5)
   ```

3. For one orientation there are nonzero odd integers
   `N_ab,N_ar,N_br` satisfying `(5)`,

   ```text
   gcd(a,b)|N_ab,  gcd(a,r)|N_ar,  gcd(b,r)|N_br,     (6)

   (2r)N_ab-bN_ar+aN_br=0.                            (7)
   ```

The cited ten-speed result supplies at least one safe phase of the divided
pack. If the full row were non-lonely, the three exceptions would have to
spoil all four lifts of that phase. Consequently every hypothetical
non-lonely row at this boundary forces the genuine support-three relation

```text
e N_ab-bN_ar+aN_br=0,
|N_ab|+|N_ar|+|N_br|<2(a+b+e)/7.                      (8)
```

Absence of a defect triple satisfying all of `(5)--(7)` closes the row.
Existence says only that some abstract selector phase can be spoiled; it
does not imply that a safe phase of the divided pack is spoiled, and it does
not imply non-loneliness.

## 3. From four labels to three intervals

An odd exception has order four on the lifts and can spoil at most one label.
The speed `e=2r` has order two and, when dangerous, spoils exactly one parity
pair `{j,j+2}`. Full spoilage therefore forces the three bad-label sets to
be disjoint of sizes `2,1,1`.

Choose the label killed by `a` as the origin and swap `a,b` if necessary.
For one real lift `z`, the partition is

```text
a at z,
b at z+1/2,
2r at z+1/4 and z+3/4.                                (9)
```

The last two phases agree because their difference after multiplication by
`2r` is the integer `r`. With

```text
D_w={z:||wz||<1/14},
```

full spoilage is therefore equivalent to

```text
z in D_a intersect (D_b-1/2) intersect (D_(2r)-1/4). (10)
```

Conversely, `y={4z}` reconstructs a fully spoiled four-lift phase; the
integer part of `4z` only rotates the labels.

## 4. Interval-to-lattice proof

Lift the three selected open components in `(10)` to the real line. Their
centres and radii are

```text
c_a=A/a,        R_a=1/(14a),
c_b=B/b-1/2,    R_b=1/(14b),
c_r=C/(2r)-1/4, R_r=1/(28r).                          (11)
```

Three real intervals intersect iff they intersect pairwise. Because they are
open, all three inequalities are strict. Direct subtraction gives

```text
c_a-c_b=N_ab/(2ab),
c_a-c_r=N_ar/(4ar),
c_b-c_r=N_br/(4br),                                  (12)
```

so the pairwise radius tests are exactly `(5)`. This proves `1 iff 2` and
also proves that equality in any line of `(5)` is safe: two open components
only touch.

Expanding `(4)` proves `(7)`. Every defect is odd, since its first terms
are even and its final product is odd, and `(6)` is immediate. Hence
`2 implies 3`.

For the converse, set

```text
K_ab=(N_ab-ab)/2,
K_ar=(N_ar-ar)/2,
K_br=(N_br+br)/2.                                     (13)
```

The desired centres solve

```text
bA-aB=K_ab,
2rA-aC=K_ar,
2rB-bC=K_br.                                          (14)
```

All relevant gcd moduli are odd. Thus the divisibilities of the `N` values
in `(6)` pass through division by two in `(13)` and give the three pairwise
gcd conditions for `(14)`. Equation `(7)` becomes

```text
2rK_ab-bK_ar+aK_br=0.                                 (15)
```

It remains to choose `A` with

```text
bA=K_ab (mod a),       2rA=K_ar (mod a).              (16)
```

For a prime power `p^alpha||a`, put `beta=v_p(b)` and
`rho=v_p(2r)`. Each congruence is separately soluble. If `beta<=rho`, the
first is stronger; multiplying the error in the second by `b` and using
`(15)` shows divisibility by `p^(alpha+beta)`, so division by `b` gives
the second congruence modulo `p^alpha`. The case `rho<beta` is symmetric.
CRT supplies `A`; then

```text
B=(bA-K_ab)/a,         C=(2rA-K_ar)/a
```

are integers, and `(15)` gives the third equation. Thus `3 implies 2`.

## 5. Finite decision and cheap gates

For a fixed exception triple, it is enough to enumerate

```text
0<|N_ab|<=(a+b-1)//7,
0<|N_ar|<=(a+2r-1)//7,
0<|N_br|<=(b+2r-1)//7,                                (17)
```

retaining oddness, `(6)`, and `(7)`. Equivalently, simultaneous translation

```text
(A,B,C)->(A+a,B+b,C+2r)
```

is a complete gauge, and pairwise overlap allows the finite box

```text
0<=A<a, 0<=B<=2b, 0<=C<=3r.                           (18)
```

Since all defects are nonzero and gcd-divisible, a fully spoiled phase
requires the quick gates

```text
a+b>7gcd(a,b),
a+2r>7gcd(a,r),
b+2r>7gcd(b,r).                                       (19)
```

They are not sufficient. The first small false converse is `(e,a,b)=(2,7,11)`:
the windows force three signs `+/-1`, but
`2N_ab-11N_ar+7N_br=0` has no such solution.

This is finite for each triple, not uniformly finite over the unbounded
`d=4` branch. A height, ratio, or owner bound is still needed for that.

## 6. Sharp controls

- `(2,3,5)` passes the old odd-pair gate `3+5>7`, but its quarter-shift
  gate fails.
- `(2,5,7)` has `5+2=7`; exact endpoint contact is safe because the danger
  intervals are open.
- In the pairwise-distinct LRC exception setting, the smallest possible
  maximum speed of a fully spoiled triple is nine:

  ```text
  (e,a,b)=(2,7,9), y=6/49,
  bad labels=({0,2},{1},{3}),
  (A,B,C)=(2,7,1), (N_ab,N_ar,N_br)=(1,1,1),
  2-9+7=0.                                             (20)
  ```

  Below nine, `e=2` cannot leave two distinct odd speeds passing `(19)`;
  for `e=6`, only `{a,b}={5,7}` remains, and its unit-sign relation is
  impossible.

The canonical selector hostile retains the divided pack. Take

```text
s=1,t=4,(p,q)=(1,4),
u=(8,12,20,24,28,32,36,40,2,9,11).                   (21)
```

The divided ten-speed pack is `H={1,...,10}`. Its safe phase `y=1/11` has
exception bad-label partition

```text
({0,2},{3},{1}).                                      (22)
```

Yet the same full row has clearance `1/11` at `x=21/22`. Thus `(22)` is
a hostile to arbitrary selector choice, not an LRC counterexample.

## 7. Preservation and scope

The map from four labelled lifts to a defect triple preserves label parity,
the two orientations, strict endpoints, and exact common-phase compatibility.
It destroys the safe-phase set `G(H)` of the divided pack, the rest of the
body geometry, height filters, and which cited pack witness was selected.

The needed next sidecar is the intersection of the certificate-positive
spoiled-phase set with `G(H)`. Certificate absence is a loneliness
certificate; certificate presence has no non-loneliness consequence.

This theorem is relative to THM-4004's cited ten-speed pack and THM-4024's
`d=4,c_4=8` reduction. It does not close every triple uniformly, the
remaining `d=2,3` equality boundaries, or LRC(14).

## 8. Replay

```text
python3 -B 04-computation/lrc14_d4_affine_defect_lattice_boundary_thm4030.py
python3 -B -O 04-computation/lrc14_d4_affine_defect_lattice_boundary_thm4030.py
```

Both runs reproduce the frozen raw-LF output. **QED.**
