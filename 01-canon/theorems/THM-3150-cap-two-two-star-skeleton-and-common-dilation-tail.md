---
id: THM-3150
title: "Cap-two two-star skeleton and common-dilation tail"
status: >
  PROVED + VERIFIED-EXACT.  For the reflected body H=(1,2,3,4,6,12) at
  upper-median cell 90, every six-distinct cap-two level assignment has at
  least three high primitive channels in the two-star on labels 1 and 2.
  Hence the nine primitive overlap skeletons have total mass at least 1/35.
  Along every fixed integer common-dilation ray, the exact reflected
  two-star certificate closes from dilation 62 onward.  This is not an
  arbitrary-gcd cone theorem, a physical-survivor classification, or LRC(14).
source: root/frontier-synthesis/2026-08-02
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-3135-directed-cycle-weak-order-lane-cover-and-reflected-h-boundary
script: 04-computation/lrc_H_cap_two_two_star_skeleton_thm3150.py
output: 05-knowledge/results/lrc_H_cap_two_two_star_skeleton_thm3150.out
script_sha256: ce1af624a78791b5329ec4da8d33e849d74d5ca7d510522d19199fd88676e966
output_sha256: 3adc10e2afe0f87669e556eebe0a683262498194183fc2d66143a45df9fc4e96
hash_basis: LF-normalized bytes
---

# THM-3150 -- cap-two two-star skeleton and common-dilation tail

**PROVED + VERIFIED-EXACT.**

This theorem turns the finite two-star matching signal into one analytic
tail, while exposing the coordinate that still prevents an arbitrary-level
cone theorem.

## 1. Statement and inherited geometry

Fix

```text
H=(1,2,3,4,6,12),             L=168,             cell j=90.       (1)
```

Slots `0,1`, carrying labels `1,2`, are pivots.  Let `E_*` be the nine
edges incident to at least one pivot:

```text
01,02,03,04,05,12,13,14,15.                                (2)
```

Let `q_0,...,q_5` be six distinct positive integers satisfying

```text
max_i q_i <=2 min_i q_i.                                    (3)
```

Write an edge channel in lowest terms as `(P,Q)`.  Call it **high** when
`P+Q>=8`.  Then:

1. at least three edges in `(2)` are high;
2. the sum of the nine primitive periodized intersection skeletons is at
   least `3/105=1/35`;
3. if the pivot skeleton is blind, at ratio `1/2` or `3/2`, at least four
   edges are high and the skeleton sum is at least `4/105`.

Now fix any six distinct positive integer base levels `x_i` satisfying `(3)`
and put `q_i=n x_i`.  For every integer `n>=62`, one edge `ij` in `(2)` has
exact reflected intersection mass exceeding the full singleton debt

```text
Delta(q)=sum_i H_i/[7(q_iL-H_i)].                           (4)
```

Pairwise Bonferroni therefore closes the fixed cell on every such tail.

## 2. The finite ratio graph

A distinct cap-two channel is low exactly when its oriented ratio lies in

```text
S={1/2,2/3,3/4,4/3,3/2,2}.                                (5)
```

Normalize the pivot levels to `(1,t)`.  A leaf `x` is low to the first pivot
iff `x in S`, and low to the second iff `x/t in S`.  A leaf low to both
pivots forces `t in S/S`.  Within `[1/2,2]`,

```text
S/S={1/2,9/16,2/3,3/4,8/9,1,9/8,4/3,3/2,16/9,2}.          (6)
```

If `t` is outside `(6)`, the pivot edge is high and every leaf has low
degree at most one.  Hence at most four of the nine edges are low and at
least five are high.

Inside `(6)`, distinctness removes `t=1`, and inversion reduces the exact
finite check to five rows:

| `t` | maximum low leaf incidences | pivot low? | low-edge upper bound | high-edge lower bound |
|---:|---:|:---:|---:|---:|
| `1/2` | `4` | yes | `5` | `4` |
| `9/16` | `5` | no | `5` | `4` |
| `2/3` | `4` | yes | `5` | `4` |
| `3/4` | `4` | yes | `5` | `4` |
| `8/9` | `6` | no | `6` | `3` |

For each row, every positive-degree leaf belongs to
`(S union tS) minus {1,t}`.  Enumerating subsets of at most four such values,
with the cap retained, gives the displayed exact maxima.  Degree-zero leaves
cannot raise them.  The count three is sharp: at `t=8/9`, the leaves

```text
2/3, 3/4, 32/27, 4/3                                      (7)
```

give six low and three high edges, and the full six-level cap is exactly two.

THM-2941 proves that every high channel has primitive periodized overlap at
least `1/105`, independently of its phase.  Intersections are nonnegative,
so the first two assertions follow.

For the pivot labels, their cell residues are

```text
r_0=1*90 mod 168=90,             r_1=2*90 mod 168=12.       (8)
```

Checking `(5)` in the exact primitive phase condition shows that the pivot
skeleton is blind only at `t=1/2,3/2`.  The corresponding rows of the table
and its inverse have four high edges, proving assertion three.

## 3. Reflected transport and the common-dilation tail

For an edge `ij`, put

```text
g_ij=gcd(q_i,q_j),          q_i=g_ij P, q_j=g_ij Q.          (9)
```

The proved THM-2941 primitive-fiber transport gives, on every high edge,

```text
I_ij >=1/105-4(H_i+H_j)/(g_ij L).                          (10)
```

The map sends the exact reflected arc pair to its periodized primitive tent
fibre.  It preserves a lower overlap bound and forgets the two coefficient
shears; `g_ij,H_i,H_j,L` are the sidecar paying the loss in `(10)`.

All low-edge intersections remain nonnegative.  Subtracting the losses on
all nine edges is therefore conservative, and

```text
sum_(ij in E_*)(H_i+H_j)=65.                               (11)
```

On the common-dilation ray `q_i=n x_i`, every `g_ij>=n`.  Thus

```text
sum_(ij in E_*) I_ij >=1/35-65/(42n),
max_(ij in E_*) I_ij >=(1/9)[1/35-65/(42n)].                (12)
```

Since every `x_i>=1`,

```text
Delta(nx)<=B(n):=sum_(e in H)e/[7(nL-e)].                  (13)
```

At the exact threshold,

```text
1/35-65/(42*62)=47/13020,

(1/9)(47/13020)-B(62)
 =6445273605987946151/383853218066532083354940>0.           (14)
```

The left side of `(12)` increases with `n`, while every summand of `B(n)`
decreases.  Hence `(14)` proves every `n>=62`.  At `n=61`, the same coarse
sufficient margin is

```text
-25676690027983908791/734017257871022240108370<0.           (15)
```

Equation `(15)` is the exact failure boundary of this bound, not a
counterexample to the actual two-star certificate.

## 4. Exact referee and failure boundary

The companion independently constructs `(5)--(7)`, checks the pivot-blind
rows, verifies `(11)`, the monotonicity identities, and the exact fractions
`(14)--(15)`.  Ordinary and optimized Python transcripts byte-match.  Run

```text
python3 04-computation/lrc_H_cap_two_two_star_skeleton_thm3150.py
python3 -O 04-computation/lrc_H_cap_two_two_star_skeleton_thm3150.py
```

For arbitrary cap-two levels, the same argument gives only

```text
sum_(ij in E_*) I_ij >=1/35
 -(4/L)sum_(ij in E_*)(H_i+H_j)/gcd(q_i,q_j).                (16)
```

Large absolute levels can have small pairwise gcds, so `(16)` need not
improve with height.  The skeleton theorem is uniform, but this shear
transport is not.  A full `m>=D` cone theorem needs a sharper transport or a
new sidecar controlling the nine gcd/shear terms.  No arbitrary-gcd cone,
physical-survivor classification, or proof of `LRC(14)` is claimed.  QED.
