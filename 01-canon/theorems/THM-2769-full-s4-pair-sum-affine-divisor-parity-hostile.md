---
id: THM-2769
title: "Full-S4 pair-sum affine divisor-parity hostile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The one-parameter quartic
  F_t(Y)=Y^4-2Y^2-8tY+1-4t has irreducible S3 squared-pair-sum cubic
  R_t(U)=U^3-4U^2+16tU-64t^2, square root product (8t)^2, Kummer rank two,
  and generic full S4 splitting group.  Its six-root pair-sum discriminant is
  an exact square, but over t=0 the three Kummer root valuations are (1,1,0),
  so the V4 cover has the nonzero even boundary word 110 and ramifies.
  Therefore full field monodromy and the square product do not supply the
  quasi-etale/unit/class-group carrier.  This is an explicit non-Keller
  affine boundary hostile, not a Keller map, JC(2), or DC(2) result.
source: root/thm2769-full-s4-divisor-parity-2026-07-28
depends_on:
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2766-quadratic-cubic-pullback-even-sign-kummer-plane-and-weyl-d3-s4
related:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2758-quartic-pair-sum-sextic-resolvent-pullback-and-discriminant-square
  - THM-2762-quartic-opposite-sum-wall-imprimitive-d4-and-keller-exclusion
script: 04-computation/jacobian_full_s4_affine_divisor_parity_hostile_thm2769.py
output: 05-knowledge/results/jacobian_full_s4_affine_divisor_parity_hostile_thm2769.out
script_sha256: 01d6b992cf58e1900eec0f88a22bff0ef32682954a45a00314459029c02c0fe2
output_sha256: 1f2a73699c26ea65386b81a18890a8e9fdc8a7e16bf20e1f5ca0c2e4280f922b
hash_basis: LF-normalized bytes
---

# THM-2769 -- a full-`S4` field packet can still ramify in the `V4` layer

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The square product in THM-2766 removes odd sign flips.  It does not remove a
double sign flip at a divisor.  The following one-parameter quartic already
realizes that distinction with genuine generic `S4` monodromy.

## 1. The quartic and its pair-sum cubic

Work over

```text
K=C(t)
```

and put

```text
F_t(Y)=Y^4-2Y^2-8tY+1-4t.                              (1)
```

For a depressed quartic

```text
Y^4+pY^2+qY+r                                           (2)
```

let `s_1,s_2,s_3` be one choice from each opposite pair of complementary
pair sums, and put `tau_i=s_i^2`.  Direct symmetric expansion gives

```text
product_i(U-tau_i)
 =U^3+2pU^2+(p^2-4r)U-q^2.                             (3)
```

Conversely, a monic cubic

```text
U^3+AU^2+BU+C,
```

with `-C=q^2` reconstructs exactly the two depressed quartics obtained from

```text
p=A/2,                 r=(p^2-B)/4,             q=+/-sqrt(-C). (3a)
```

Thus this is an exact two-branch inverse in the displayed pair-sum
coordinate, not only a field-level group construction.

For `(1)`, `p=-2`, `q=-8t`, and `r=1-4t`, so `(3)` is

```text
R_t(U)=U^3-4U^2+16tU-64t^2.                            (4)
```

Indeed `(3a)` applied to `(4)` gives `p=-2`, `r=1-4t`, and
`q=+/-8t`; the negative branch is `(1)`.

Consequently the six pair sums have polynomial

```text
H_t(V)=R_t(V^2)=V^6-4V^4+16tV^2-64t^2,                (5)
```

and

```text
tau_1 tau_2 tau_3=(8t)^2.                              (6)
```

In the invariant notation of THM-2758, `(1)` has `e_1=0`, `e_3=8t`, hence

```text
T=e_1^3-4e_1e_2+8e_3=64t.                              (7)
```

Thus the generic point is on the live `T!=0` side, not on THM-2762's
imprimitive wall.

## 2. The cubic quotient is genuinely `S3`

The cubic `(4)` is irreducible over `C(t)`.  If it had a root in `C(t)`,
monicity would make that root integral over the integrally closed ring
`C[t]`, hence an element of `C[t]`.  It would divide `64t^2`, so it would be

```text
c t^j,                 c in C*,             j in {0,1,2}. (8)
```

Substitution closes the three cases:

```text
j=0: coefficient of t^2 is -64;
j=1: coefficient of t^3 is c^3;
j=2: coefficient of t^2 is -64.                         (9)
```

Therefore there is no rational root, and a cubic with no rational root is
irreducible.

The exact discriminants are

```text
disc(R_t)=disc(F_t)
         =-2^12 t^2(27t^2-14t+3).                      (10)
```

The remaining quadratic factor has discriminant

```text
(-14)^2-4*27*3=-128,                                   (11)
```

so it is squarefree and is not a square in `C(t)`.  Hence `(10)` is
nonsquare.  Since `(4)` is irreducible and separable in characteristic zero,
its splitting field `E/K` has group

```text
Gal(E/K)=S3.                                            (12)
```

## 3. The divisor over `t=0` has word `110`

Let `A` be the normalization of `C[t]` in `E`.  The Newton points of `(4)`
at the `t`-adic valuation are

```text
(0,2),(1,1),(2,0),(3,0).                               (13)
```

Their lower hull has slopes `-1,0` and horizontal lengths `2,1`; therefore
the root valuations are

```text
(1,1,0).                                                (14)
```

Here there is no hidden ramification rescaling.  Put `U=tZ`.  Then

```text
R_t(tZ)=t^2(tZ^3-4Z^2+16Z-64).                         (15)
```

At `t=0`, the two finite `Z` roots satisfy

```text
Z^2-4Z+16=0,                 discriminant=-48,          (16)
```

and are distinct over `C`; the third root lifts from the simple root `U=4`,
where the derivative is `16`.  Thus all three labelled roots lie in the
completed local field `C((t))`, and there is a prime divisor `D` of `A`
with normalized valuation row

```text
nu_D=(v_D(tau_1),v_D(tau_2),v_D(tau_3)) mod 2=110       (17)
```

after relabelling.

Let

```text
W=<[tau_1],[tau_2],[tau_3]> subset E*/E*2.              (18)
```

Relation `(6)` gives an `S3`-equivariant surjection

```text
M=F2^3/<(1,1,1)> -> W.                                 (19)
```

The two-dimensional module `M` is irreducible over `F2`: equivalently, its
dual is the even-weight code

```text
{000,110,101,011},                                      (20)
```

and `S3` acts transitively on the three nonzero words.  A proper nonzero
submodule would be a line and would have a unique nonzero vector fixed as a
set, contradicting that transitivity.

Row `(17)` says that `(19)` is not the zero map.  Its kernel is an `S3`
submodule of the irreducible module `M`, so the kernel is zero.  Therefore

```text
dim_F2 W=2.                                             (21)
```

## 4. The generic quartic group is full `S4`

Equations `(6)`, `(12)`, and `(21)` meet exactly the full-cubic branch of
THM-2766.  The splitting field of `H_t` therefore has group

```text
V4 semidirect S3=W(D3)=S4.                              (22)
```

For completeness, this is also the splitting field of `(1)`, not merely an
abstract sextic.  Choose `s_i^2=tau_i` with

```text
s_1s_2s_3=8t.                                           (23)
```

The four even-sign half-sums of the `s_i` have elementary symmetric data

```text
e_1=0,
e_2=-(tau_1+tau_2+tau_3)/2=-2,
e_3=s_1s_2s_3=8t,
e_4=((sum tau_i)^2-4 sum_(i<j)tau_i tau_j)/16=1-4t.     (24)
```

They are precisely the roots of `(1)`.  Conversely their complementary pair
sums recover the six roots of `(5)`.  The two splitting fields coincide, so

```text
Gal(F_t/K)=S4.                                          (25)
```

## 5. What the square discriminant forgets

THM-2766's pullback formula and `(10)` give

```text
disc(H_t)
 =2^36 t^6(27t^2-14t+3)^2
 =(2^18 t^3(27t^2-14t+3))^2.                           (26)
```

Thus the six-edge discriminant is an exact square, just as field-level
even-sign monodromy predicts.  Nevertheless `(17)` is nonzero.  On the
normal affine curve `Spec(A)`, the normalization in

```text
E(sqrt(tau_1),sqrt(tau_2))                              (27)
```

ramifies at `D` by THM-2685.

The precise loss map is

```text
partial: W -> direct_sum_D F2,
[f] |-> (v_D(f) mod 2)_D.                               (28)
```

The square product `(6)` preserves only

```text
nu_D in {000,110,101,011}.                              (29)
```

It forgets the distinction between `000` and the three nonzero even words.
A nonzero word is a double sign flip, hence an even permutation in the
quartic and six-edge actions; discriminant parity cannot see it.  Quasi-etale
extension requires the strictly stronger condition

```text
nu_D=000 for every D.                                   (30)
```

The sharp control is the ramified even base change `t=s^2`: it doubles the
row `(1,1,0)` to `(2,2,0)`, hence changes `110` to `000`.  This does not repair
the original model; it confirms THM-2685's warning that even ramified base
change can erase the obstruction.

## 6. Keller/Jelonek consequence and strict scope

This family is deliberately not Keller.  Its special fibre is

```text
F_0(Y)=(Y^2-1)^2,                                      (31)
```

and the Kummer layer ramifies divisorially.  The divisor `D` above is not
being asserted to be the Jelonek divisor of a polynomial Keller map.

The result instead gives a concrete necessary test for any proposed
`A4/S4` Keller resolvent.  On the normalization of the affine target in the
cubic matching field, compute the normalized root-valuation row `(17)` at
every prime divisor.  Any odd--odd--even row excludes the model immediately.
Raw Newton slopes on the base are sufficient only after accounting for the
ramification index; in this example `(15)`--`(16)` prove that the splitting
normalization is unramified over `t=0`.  If every row is zero, one must still
form the half-divisors and locate the standard plane in global units or
`Cl[2]`, as required by THM-2655/2685.

Accordingly, `(1)` proves that all of the following field data can coexist
with failure of the affine gate:

```text
irreducible cubic quotient with group S3;
rank-two Kummer plane;
full quartic group S4;
base-field square product;
square six-edge discriminant.                           (32)
```

It does **not** construct or exclude a degree-four Keller map, classify
Jelonek divisors, or prove `JC(2)`, `DC(2)`, or any higher-dimensional
Jacobian statement.

## 7. Exact companion

Run

```bash
python 04-computation/jacobian_full_s4_affine_divisor_parity_hostile_thm2769.py
python -O 04-computation/jacobian_full_s4_affine_divisor_parity_hostile_thm2769.py
```

The no-`assert`, no-floating-point companion verifies `(3)`--`(7)`, including
the two inverse branches `(3a)`, the
three rational-root cases `(9)`, both discriminants, the exact square in
`(26)`, the Newton hull and local controls `(13)`--`(16)`, all `S3`-invariant
subspaces of `(20)`, the full 24-element even-sign action on four states, the
boundary code, the even-base-change control, and `(31)`.  Normal and optimized
runs LF-normalized-byte-match the stored 21-line transcript.

QED.
