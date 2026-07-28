---
id: THM-2775
title: "Modular S4 to Weyl-D3 generator frame and affine-parity blindness"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  The modular (2,3,4) quotient generators act on the three quartic
  opposite-edge coordinates by explicit signed-permutation matrices
  generating W(D3)=S4.  The half-Hadamard root states recover (12),(234),
  while the (2,3,3) binary generator is the diagonal V4 element (12)(34).
  The order-four face square is diag(-1,-1,1), exactly the 110 inertia word
  realized by THM-2769.  Thus the modular/Weyl co-occurrence explains the
  binary-ternary frame but cannot force affine quasi-etaleness.  No Keller
  monodromy identification, JC(2), DC(2), graceful, or LRC result follows.
source: root/modular-weyl-d3-generator-frame-2026-07-28
depends_on:
  - THM-2758-quartic-pair-sum-sextic-resolvent-pullback-and-discriminant-square
  - THM-2766-quadratic-cubic-pullback-even-sign-kummer-plane-and-weyl-d3-s4
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
related:
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy
script: 04-computation/modular_s4_weyl_d3_generator_frame_thm2775.py
output: 05-knowledge/results/modular_s4_weyl_d3_generator_frame_thm2775.out
script_sha256: 758edc501f1f0cd6939e9e634611043cb42670e92286ba12c5b657c98f085e84
output_sha256: 4c62a2409bd5ecfec67200d820bb03840ab5f4d650214e784a075e898fae01b5
hash_basis: LF-normalized bytes
---

# THM-2775 -- the modular `S4` quotient in opposite-edge coordinates

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

The two recurring binary/ternary constructions now meet on one explicit
representation.  The meeting is exact, but it also displays the missing
affine coordinate: the same order-four face whose square creates the quartic
`V4` kernel admits nontrivial divisor inertia in that kernel.

## 1. Opposite-edge coordinates

Let `r_1,...,r_4` be four marked quartic roots and put

```text
d_1=r_1+r_2-r_3-r_4,
d_2=r_1+r_3-r_2-r_4,
d_3=r_1+r_4-r_2-r_3.                                  (1)
```

These are THM-2758's oriented opposite-edge differences.  On column vectors
`d=(d_1,d_2,d_3)^T`, define

```text
X_S = [ 1  0  0 ]       X_A = [ 1  0  0 ]
      [ 0  0 -1 ]             [ 0 -1  0 ]
      [ 0 -1  0 ],            [ 0  0 -1 ],

Y_S = [ 0  1  0 ]       Y_A = [ 0  0 -1 ] .             (2)
      [ 0  0  1 ]             [ 1  0  0 ]
      [ 1  0  0 ]             [ 0 -1  0 ]
```

Thus

```text
X_S(d_1,d_2,d_3)=(d_1,-d_3,-d_2),
X_A(d_1,d_2,d_3)=(d_1,-d_2,-d_3),
Y_S(d_1,d_2,d_3)=(d_2,d_3,d_1),
Y_A(d_1,d_2,d_3)=(-d_3,d_1,-d_2).                      (3)
```

The subscripts distinguish the octahedral and tetrahedral markings.  Direct
substitution in `(1)` shows that `(X_S,Y_S)` is the pullback action of
`((12),(234))`, while `(X_A,Y_A)` is that of `((12)(34),(123))`, exactly as
in THM-2768, with the standard left pullback convention.

Every matrix in `(2)` is a signed permutation with an even number of sign
flips.  In particular it preserves

```text
d_1d_2d_3=T,                                             (4)
```

the oriented product whose square gives the pair-sum Kummer product.

## 2. The two triangle relations become literal matrices

Matrix multiplication gives

```text
X_S^2=Y_S^3=I,           order(X_S Y_S)=4,
(X_S Y_S)^2=diag(-1,-1,1),                           (5)

X_A^2=Y_A^3=I,           order(X_A Y_A)=3.            (6)
```

Here `X_S Y_S` means first `Y_S`, then `X_S`; reversing the convention conjugates
the product and does not change its order.  Hence `(5)` and `(6)` are exactly
the `(2,3,4)` and `(2,3,3)` face relations of THM-2768.

The coordinate-permutation images of `X_S,Y_S` are a transposition and a
three-cycle, so they generate `S3`.  The element `(X_S Y_S)^2` is a nontrivial
even diagonal sign flip.  Its three `S3` conjugates generate the full
diagonal `V4`.  Therefore

```text
<X_S,Y_S>=V4 semidirect S3=W(D3),          order 24.    (7)
```

Similarly `X_A` is already a nontrivial diagonal `V4` element, its `Y_A`
conjugates generate `V4`, and

```text
<X_A,Y_A>=V4 semidirect C3,                order 12.    (8)
```

The determinant is `+1` throughout `(8)`, while `X_S` has determinant `-1`.
Under the four-state action below, `(7)` and `(8)` are `S4` and `A4`.

## 3. Half-Hadamard states recover the four roots

Define the four centered root forms

```text
h_1= d_1+d_2+d_3,       h_2= d_1-d_2-d_3,
h_3=-d_1+d_2-d_3,       h_4=-d_1-d_2+d_3.              (9)
```

They satisfy `h_i=4r_i-(r_1+r_2+r_3+r_4)` and sum to zero.  Pullback by
the matrices `(2)` gives

```text
X_S: (h_1 h_2),                Y_S: (h_2 h_3 h_4),
X_A: (h_1 h_2)(h_3 h_4),       Y_A: (h_1 h_2 h_3).     (10)
```

Thus the Weyl isomorphism `W(D3)=S4` is not merely an order coincidence:
`(9)` conjugates the signed three-coordinate action to the literal four-root
permutation action.  It also makes the quotient distinction exact:

```text
S4/V4=S3:  the binary and ternary generators both survive;
A4/V4=C3:  the binary generator lies in V4 and disappears. (11)
```

The modular kernels remain those of THM-2768:

```text
ker(C2*C3 -> A4)=F3,       ker(C2*C3 -> S4)=F5.          (12)
```

In the `S4 -> S3` matching quotient both finite-factor orders survive, and
the index-six kernel is `F2`.  In the `A4 -> C3` quotient the binary factor
collapses; its preimage is the `C2*C2*C2` carrier recorded in THM-2768.

## 4. The same frame exposes the affine blindness

Put `tau_i=d_i^2/4`.  At a prime divisor `D` of the cubic matching
normalization, the Kummer boundary row is

```text
nu_D=(v_D(tau_1),v_D(tau_2),v_D(tau_3)) mod 2.          (13)
```

The square-product condition `(4)` forces only

```text
nu_D in C_even={000,110,101,011}.                        (14)
```

Signed flips do not alter valuations, while the coordinate `S3` permutes the
three entries.  Consequently `000` is fixed and the three nonzero even words
form one allowed orbit.  Neither the triangle relations nor the invariant
product upgrades membership in `(14)` to the quasi-etale condition
`nu_D=000`.

This boundary is sharp.  THM-2769 gives a one-parameter quartic with generic
full `S4` group and, after relabelling, the divisor row

```text
nu_D=110.                                                (15)
```

Local Kummer inertia then flips the first two square roots and fixes the
third, so its signed matrix is

```text
diag(-1,-1,1)=(X_S Y_S)^2.                              (16)
```

Thus the obstructing inertia is not outside the modular/Weyl frame; it is
the square of its order-four face.  The missing predicate is geometric:
whether a `V4` element occurs as inertia at an affine divisor.  A finite
group presentation, square discriminant, or matching quotient does not
record that placement.  THM-2685's normalized valuation rows, half-divisors,
units, and class-group two-torsion are the required sidecar.

## 5. Exact verification and scope

Run

```bash
python 04-computation/modular_s4_weyl_d3_generator_frame_thm2775.py
python -O 04-computation/modular_s4_weyl_d3_generator_frame_thm2775.py
```

The no-`assert`, integer-only companion verifies every matrix relation,
enumerates the `24`- and `12`-element images, identifies all signed
permutations and their coordinate quotients, checks the Hadamard actions,
the kernel-rank ledger, and the complete even-code orbit.  Normal and
optimized runs byte-match the stored transcript.

```text
PROVED HERE (candidate):  explicit modular generators in D3 coordinates;
                          (2,3,4) and (2,3,3) relations;
                          W(D3)=S4 and orientation A4 images;
                          literal four-root Hadamard action;
                          V4 matching-quotient distinction;
                          face-square/boundary-word identity (16);
                          exact affine-parity stopping mechanism.

NOT PROVED:               a modular marking on every quartic extension;
                          a Keller/Jelonek geometric realization;
                          vanishing of affine boundary inertia;
                          exclusion of A4 or S4 Keller monodromy;
                          JC(2), DC(2), Graceful Tree, or LRC(14).          (17)
```

QED (candidate).
