---
id: THM-3045
title: "K4 edge isotypic binary-ternary integral clutch"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  The six-edge permutation lattice of K4 has rational S4 decomposition
  1+[22]+[31].  The sum of the three integral intersection lattices has
  index 24 and Smith form 1,1,1,2,2,6.  Its exact equivariant quotient is
  F2[three matchings] plus one trivial F3 scalar: the binary part is the
  opposition clutch of THM-2756, while the ternary part is the failure of
  the matching permutation lattice to split integrally into its constant
  and A2 augmentation summands.  This is a finite representation-lattice
  theorem, not an identification of the binary/ternary trees and not a
  Keller, affine-owner, tournament, or LRC exclusion.
source: codex-k4-edge-integral-clutch-2026-08-01
depends_on:
  - THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch
related:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2774-tree-path-smith-index-ladder-and-binary-ternary-lattice-defects
  - THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding
  - THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary
script: 04-computation/k4_edge_isotypic_binary_ternary_clutch_thm3045.py
output: 05-knowledge/results/k4_edge_isotypic_binary_ternary_clutch_thm3045.out
script_sha256: af89e4060f293da3222db3e85f02d1279fa217a63957add2b7f6afd746722e2d
output_sha256: ae4f76b13bc060bb653c42a8176c08b53dbb91678a523588c8f272029c882c9e
hash_basis: LF-normalized bytes
---

# THM-3045 -- the six-edge lattice carries separate binary and ternary clutches

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

## 1. Inheritance and exact statement

[THM-2756, the opposite-edge integral
clutch](THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch.md)
splits the rational six-edge module of `K4` into two rank-three opposition
blocks and proves that their integral intersection lattices miss the ambient
edge lattice by an equivariant `F2^3`.  Its positive block is the permutation
module on the three perfect matchings,

```text
Q[M] = 1 + [22],          |M|=3.                           (1)
```

The integral constant/augmentation split inside `(1)` has not yet been
computed.  It contributes an independent ternary clutch.

Let `X={0,1,2,3}`, `E=binom(X,2)`, and

```text
M={01|23,02|13,03|12}.                                    (2)
```

Put `L=Z[E]`.  For `m={e,e_bar}` write

```text
p_m=e+e_bar,                 n_m=e-e_bar,                 (3)
```

where the sign of `n_m` is an arbitrary orientation gauge.  Define

```text
L_0  = Z c,                   c=sum_(e in E)e;
L_22 = {sum_m a_m p_m : sum_m a_m=0};
L_31 = direct_sum_m Z n_m.                                 (4)
```

Then these are exactly the intersections of `L` with the three rational
isotypic summands `1,[22],[31]`, and

```text
L_Q = (L_0)_Q direct_sum (L_22)_Q direct_sum (L_31)_Q.     (5)
```

The integral sum has

```text
L/(L_0+L_22+L_31)  ~=  F2[M] direct_sum F3_triv,          (6)
[L:L_0+L_22+L_31]  = 24,
Smith = diag(1,1,1,2,2,6).                                (7)
```

The isomorphism is explicit and `S4`-equivariant:

```text
Phi(x)=(kappa(x),tau(x)),
kappa(x)_m=x_e+x_e_bar mod2,
tau(x)=sum_(e in E)x_e mod3.                               (8)
```

`S4` permutes the three binary coordinates through its matching quotient
`S4/V4=S3`, while it fixes the ternary scalar.  Consequently the only primes
obstructing this canonical integral isotypic split are `2` and `3`.

## 2. The rational decomposition

Let `O(e)=e_bar=X\e`.  It is central in the six-edge permutation action.
The opposition projectors

```text
P_plus=(1+O)/2,                 P_31=(1-O)/2               (9)
```

have images spanned by the `p_m` and `n_m`, respectively.  Within the
positive block put

```text
P_0=(1/6) J_6,                 P_22=P_plus-P_0.            (10)
```

These are pairwise orthogonal idempotents of ranks `1,2,3` summing to the
identity.  The first image is the constant representation, the second is the
zero-sum matching module `[22]`, and THM-2756 identifies the third with the
faithful standard sheet module `[31]`.  This proves `(5)` and the descriptions
in `(4)`.

In the ordered edge basis

```text
01,02,03,12,13,23,                                      (11)
```

take the integral rows

```text
c,
p_01|23-p_03|12,  p_02|13-p_03|12,
n_01|23,          n_02|13,          n_03|12.             (12)
```

Their three orthogonal Gram blocks have determinants

```text
6,                         12,                         8. (13)
```

The determinant of `(12)` therefore has absolute value `24`; exact integer
row reduction gives `(7)`.

## 3. The quotient map and its kernel

Both components of `(8)` kill every row of `(12)`: an `n_m` has opposite
coordinates, an augmentation difference has total matching coefficient zero,
and `c` has total edge sum six.  Conversely, suppose `Phi(x)=0`.  The three
binary equations say that the coordinates on every opposite pair have the
same parity, so

```text
a_m=(x_e+x_e_bar)/2,       b_m=(x_e-x_e_bar)/2            (14)
```

are integers and

```text
x=sum_m a_m p_m + sum_m b_m n_m.                          (15)
```

The ternary equation says `2 sum_m a_m=0 mod3`, hence
`sum_m a_m=0 mod3`.  Put `t=(sum a_m)/3`.  Then

```text
sum_m a_m p_m = t c + sum_m (a_m-t)p_m,                   (16)
```

and the second term belongs to `L_22`.  Thus `ker Phi` is precisely the
integral sum in `(6)`.

The map is surjective.  An individual edge supplies a unit binary matching
coordinate, while `p_m` has zero binary image and ternary image `2`, a
generator of `F3`.  This proves `(6)` and `(7)` without relying on the Smith
calculation.

Equivariance is equally direct.  A sheet permutation permutes opposite edge
pairs, so it permutes `kappa`; the total coordinate sum `tau` is invariant.
The kernel of the matching action is exactly `V4`.  Hence `(6)` is an
`S4`-module statement, not merely an abstract abelian-group coincidence.

## 4. Why the primes two and three are separate

There is a canonical two-step filtration

```text
L_0+L_22+L_31  subset  L_plus+L_31  subset  L,            (17)
```

where `L_plus=direct_sum_m Z p_m`.  Its successive quotients are

```text
(L_plus+L_31)/(L_0+L_22+L_31) ~= Z/3,
L/(L_plus+L_31)                    ~= F2[M].               (18)
```

The second is exactly THM-2756's opposition clutch.  The first is the familiar
index-three failure of the permutation lattice `Z^3` to equal its constant
line plus its `A2` root lattice.  Since the orders are coprime, primary
decomposition canonically turns `(18)` into the direct sum `(6)`.

Equivalently, all three rational projectors `(9)--(10)` send `x in L` back to
`L` exactly when

```text
kappa(x)=0 and tau(x)=0.                                  (19)
```

After localization away from `6` the decomposition is integral.  At prime
`2` only the matching-permutation clutch survives; at prime `3` only the
trivial scalar clutch survives; at every prime at least five there is no
defect.  This is the precise sense in which `2` and `3` are special for this
one six-edge carrier.

## 5. Relation to the modular, quartic, and tree pictures

The result sharpens, but also limits, the motivating synthesis.

- [THM-2606](THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin.md)
  proves that `C2*C3 -> S3` acts on the three matching directions and that the
  six edges form the octahedral carrier.  Equation `(6)` adds the exact
  integral lattice left invisible by that finite-set action.
- [THM-3037](THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary.md)
  computes the genuinely modular cohomology bit in `H^1(C2*C3,V4)`.  That bit
  is not the ternary scalar in `(6)` and must not be conflated with it.
- [THM-2993](THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding.md)
  uses the same three matchings to assemble the four root-stars and four
  false-owner triangles.  The present quotient does not choose one of those
  eight transversals or an affine owner.
- [THM-2774](THM-2774-tree-path-smith-index-ladder-and-binary-ternary-lattice-defects.md)
  shows that longer tree paths carry every later cyclic index.  Thus `2` and
  `3` are the only primes in this canonical isotypic split, not the only
  primes that can occur in related frame lattices.

In particular, the binary fraction tree and ternary Pythagorean tree are not
literally the two free factors on this lattice.  The proved co-occurrence is
more precise and less sweeping: opposition contributes a binary gluing module,
while the constant/augmentation split of the three matching channels
contributes a ternary scalar gluing module.

For a quartic cover, `[22]` is the matching-resolvent channel and `[31]` is the
standard sheet channel.  Equations `(6)--(19)` do not prove that a graph order
contains a singleton idempotent, identify a source-coordinate trace, or turn
the equality of quartic and resolvent discriminants into an affine section.
Those are the separate order/owner boundaries of THM-3038 and THM-3042.

## 6. Exact verification and scope

Run

```bash
python 04-computation/k4_edge_isotypic_binary_ternary_clutch_thm3045.py
python -O 04-computation/k4_edge_isotypic_binary_ternary_clutch_thm3045.py
```

Both executions byte-match the stored transcript.  The companion uses an
explicit `require` function and no truth-bearing Python assertions.  It checks
the three rational projectors, their ranks and orthogonality; the Gram blocks,
index and Smith form; all `6^6=46,656` residue vectors for `(19)` and the exact
24-element quotient; and every sheet permutation for equivariance, the `V4`
kernel, and the order-two/order-three matching controls.

```text
PROVED HERE:       rational 1+[22]+[31] decomposition;
                   exact integral intersection lattices;
                   index 24 and Smith 1,1,1,2,2,6;
                   explicit F2[M]+F3 quotient and projector criterion;
                   S4 equivariance and prime-local split.

NOT PROVED:        literal equality of binary/ternary trees with C2,C3;
                   a new quartic discriminant or resolvent identity;
                   integral singleton or affine-source ownership;
                   a six-vertex tournament orientation;
                   a physical LRC packet/current;
                   Keller, G1, JC(2), DC(2), or LRC(14).                 (20)
```

QED.
