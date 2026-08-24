---
id: THM-3909
title: "Equianharmonic E6-star marked zero-contact shell"
status: >
  PROVED + VERIFIED-EXACT.  A second independent hostile reconstruction of
  THM-3888 confirms its binary-quartic invariants, birational map and inverse,
  divisor of T, II^4+IV fibre packet, Mordell--Weil rank six and trivial
  torsion, and preserves the repaired one-way polynomial-section statement.
  The full Mordell--Weil height lattice is E6*.  Its 126 nonzero sections
  disjoint from O split into 6 constant-u, 48 linear-u and 72 quadratic-u
  sections.  Exactly nine are also disjoint from both Q+ and Q-: P0 and eight
  same-component sections.  The latter eight have nonpolynomial inverse
  coordinates, so strict-boundary avoidance is provably not a converse.
  This is a finite geometric shell, not a Keller counterexample; JC(2) is OPEN.
source: jc_elliptic_s_integral / independent post-THM-3888 audit, 2026-08-23
audit: >
  HOSTILE RECONSTRUCTION PASS.  The companion does not import either THM-3888
  script.  It symbolically rebuilds the quartic/cubic algebra and independently
  projects all 240 E8 roots away from an A2 subsystem.  It exhausts all 540
  ordered marked minimal-vector pairs, with the same 9=(8+1) answer in every
  case, and derives the complete separable degree-one parameter quartic.
  Normal and optimized streams byte-match the frozen output.  The finite
  lattice enumeration has 1,123 active gates and no Python assert.
depends_on:
  - THM-3888-f-zero-equianharmonic-jacobian-and-two-section-integrality
  - THM-3900-f-zero-generic-y-polynomial-root-color-response-classification
related:
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
  - THM-3895-f-zero-quartic-covariant-and-high-y-degree-emptiness
  - THM-3897-f-zero-residual-all-degree-global-emptiness
script: 04-computation/jc2_equianharmonic_e6_marked_zero_contact_shell_thm3909.py
output: 05-knowledge/results/jc2_equianharmonic_e6_marked_zero_contact_shell_thm3909.out
script_sha256: 219f7e5fcef35f6c63254a524b63ab85a3f030fd06b2f7393978eb1dfa4229ef
output_sha256: e0a418cf923011306760629b44b7cac6f87f4d1e5ed4f62022f3ae321706c263
semantic_sha256: af29393a69c844737e4f2ac293f6a117ca80b33b31c808ca02a02c4ef70d7947
hash_basis: raw LF bytes
---

# THM-3909 -- the marked elliptic shell is E6-star

**PROVED + VERIFIED-EXACT.**  Work over an algebraically closed field `k` of
characteristic zero and put

```text
C=overline{k(x)},                  a=x+1,
L=9x+4,                            F=15x^2+15x+4,
K=y^2-F,                           Delta=a^3L^2-K^2.       (1)
```

Retain the THM-3888 quartic and normalized cubic

```text
G^2=L^4-6aL^2T^2-8KT^3-3a^2T^4,                         (2)
v^2=K^2+L^2(u^3-a^3),
T=(v-K)/(u^2+au+a^2).                                   (3)
```

Let `S -> P^1_y` be the minimal rational elliptic surface of `(3)`, with
origin `O`, the other `T=0` section `P0`, and the two quartic boundary
sections `Q+`,`Q-`.  Then:

1. the Mordell--Weil height lattice is `E6*`;
2. the nonzero sections `R` with `R.O=0` form a 126-element shell, split by
   normalized coordinate degree as

   ```text
   deg_y(u)=0: 6,             deg_y(u)=1: 48,
   deg_y(u)=2: 72;                                      (4)
   ```

3. exactly nine of those sections also satisfy

   ```text
   R.Q+=R.Q-=0.                                         (5)
   ```

   They are `P0` and eight sections on the same `IV` component as `Q+`,`Q-`;
4. only `P0` among these nine has polynomial inverse coordinates `T,G in
   C[y]`.  Thus the other eight prove that even **total** disjointness from the
   strict transforms of `Q+ union Q-` is not sufficient for polynomiality of
   the inverse quartic coordinates.

The final assertion sharpens, rather than reverses, THM-3888's repaired
one-way statement.

## 1. Independent hostile audit of the elliptic input

For the binary quartic `(2)`, in the invariant convention used by THM-3888,
direct calculation gives

```text
I=0,
J=1728L^4Delta,
disc_T=-110592L^8Delta^2.                                (6)
```

Starting from `(2)` rather than importing the old map, define

```text
u=(G+L^2-aT^2)/(2T^2),
D=u^2+au+a^2,
v=K+TD.                                                   (7)
```

Reduction modulo `(2)` gives `(3)`, while the two inverse identities are

```text
T=(v-K)/D,                    G=(a+2u)T^2-L^2.           (8)
```

This independently confirms both directions of the dense birational map.  At
`T=0`, the two branches `G=+/-L^2` are smooth and `T` is a local parameter.
After putting `z=1/T` and `w=G/T^2`, the two infinity branches are

```text
w=+sa,                         w=-sa,       s^2=-3,       (9)
```

and are again smooth with `z` a local parameter.  Hence the four orders are
all one and

```text
div(T)=O+P0-Q+-Q-.                                      (10)
```

The polynomial `Delta(y)` has four simple roots because

```text
disc_y(Delta)=256a^6L^4(1-a)^3(9a-4)^2 !=0.             (11)
```

The short model `Y^2=X^3-64L^4Delta` therefore has four `II` fibres.  Its
minimal coefficient at `y=infinity` has order two, giving one `IV` fibre and

```text
II^4+IV,                         Euler number 12.         (12)
```

Thus `S` is a rational elliptic surface and `rho(S)=10`.  Shioda--Tate gives
rank `10-2-rank(A2)=6`.  Torsion injects into the `IV` component group
`Z/3`: the formal group and additive identity component are torsion-free in
characteristic zero.  The three-division polynomial is

```text
psi_3=3X(X^3-256L^4Delta).                              (13)
```

Its two branches would require `Delta` to be respectively a square or a
cube in `C(y)`, both excluded by the simple valuations in `(11)`.  Hence the
torsion subgroup is zero.  This checks every load-bearing geometric assertion
used below.

The polynomial implication also survives exactly in its repaired direction:
a pair `T,G in C[y]` gives a section with no finite contact with `Q+ union
Q-`.  No converse is used.  Section 5 proves that the converse is false even
after strengthening finite avoidance to total intersection zero.

## 2. The height lattice is `E6*`

The essential negative-definite lattice of a rational elliptic surface is
`E8`.  The sole reducible fibre contributes a primitive `A2` root lattice;
primitivity follows here from the already-proved trivial Mordell--Weil
torsion.  Its orthogonal root lattice is `E6`.  The unimodular gluing

```text
A2 direct-sum E6  subset E8,                 index 3       (14)
```

projects `E8` onto the full dual `E6*`.  In the positive Shioda height
normalization,

```text
MW(S/C(y)) = E6*,                    det=1/3.             (15)
```

This also follows from the discriminant formula
`disc(NS)=disc(A2)disc(MW)/|MW_tor|^2`: the left side has absolute value one,
`disc(A2)=3`, and the torsion is zero.

For an exact coordinate model, take the 240 roots of `E8` in `Q^8`, choose

```text
alpha_1=e1-e2,                    alpha_2=e2-e3,           (16)
```

and project by replacing the first three coordinates with their mean.  The
distinct projected vectors, grouped by squared norm, are

```text
norm 0:       1,
norm 2:      72,                       the E6 roots,
norm 4/3:    54,                       the two minuscule shells. (17)
```

Their multiplicities among the 240 `E8` roots are respectively `6,1,3`.
Conversely, every norm-two vector in the zero discriminant class and every
norm-`4/3` vector in a nonzero class lifts, after adding the matching `A2`
weight of norm `2/3`, to an `E8` root.  Thus `(17)` is exhaustive, not only a
list of examples.

## 3. The 126 sections disjoint from `O`

For a section `R!=O`, Shioda's self-height formula is

```text
<R,R>=2+2(R.O)-contr_IV(R).                              (18)
```

The `IV` correction is zero on the identity component and `2/3` on either
nonidentity component.  Consequently `R.O=0` gives norm two in the identity
class or norm `4/3` in a nonzero class.  Conversely each of the `72+54`
vectors in `(17)` forces `R.O=0` through `(18)`.  Hence there are exactly

```text
72+54=126                                                (19)
```

nonzero sections disjoint from `O`.

Such a section has polynomial normalized coordinates with

```text
deg_y(u)<=2,                      deg_y(v)<=3.            (20)
```

This is the usual degree-`1` fundamental-line-bundle bound, and it is also
visible directly after the infinity scaling `X_infinity=t^2X`,
`Y_infinity=t^3Y`.  In the first `IV` blow-up chart, a nonidentity section has
`X_1=tX` finite and `Y_1=t^2Y=+/-8L^2` at `t=0`.  Therefore

```text
nonidentity component  iff  deg_y(u)<=1.                (21)
```

A quadratic `u` has cubic `v` and specializes to a smooth nonzero point of
the cuspidal identity component.

### 3.1. Six constant sections

If `u` is constant, `(3)` gives a constant product

```text
(v-K)(v+K)=L^2(u^3-a^3).                                (22)
```

A nonzero constant product would make both factors units in `C[y]`, contrary
to their nonconstant difference `2K`.  Thus one factor vanishes and

```text
u^3=a^3,                         v=+K or -K.             (23)
```

There are exactly `3*2=6` such sections.

### 3.2. Forty-eight linear sections

Every linear section is obtained as follows.  Put

```text
H=a^3L^2-F^2=x^3(9x+5)^2                              (24)
```

and choose a root `R` of

```text
P(R)=-3R^4+8FR^3+6HR^2+H^2=0.                          (25)
```

This quartic is separable:

```text
disc_R(P)=-110592 x^12 a^6 L^4(9x+5)^8 !=0.            (26)
```

It also has no zero root.  Choose

```text
r^2=R,
Z=9r+H/r^3,
alpha^3=Z/L^2.                                          (27)
```

Here `Z!=0`; otherwise the coefficient equation forces `R=-F/3`, while

```text
P(-F/3)=-(F^2-3H)(F^2+H)/3 !=0.                         (28)
```

Then

```text
u=alpha(y+r),
v=+/-(y^2+(Z/2)y+3r^2)                                 (29)
```

satisfies `(3)`.  Conversely, after changing the sign of `v`, write a linear
solution as

```text
u=A y+B,                         v=y^2+c_1y+c_0.         (30)
```

Coefficient comparison first gives

```text
r=B/A,       Z=L^2A^3,       c_1=Z/2,       c_0=3r^2,  (31)
```

then `(27)` and finally `(25)`.  Thus `(29)` is complete.  The four roots
`R`, two choices of `r`, three choices of `alpha`, and two signs of `v` are
all distinct, giving

```text
4*2*3*2=48.                                              (32)
```

Together `(23)` and `(32)` account for all `54` nonidentity minimal vectors.
The remaining `72` sections are therefore exactly the quadratic-`u` identity
shell.  This proves the split `(4)` without solving a seventy-two-branch
coefficient ideal.

## 4. The marked zero-contact shell has size nine

The two boundary sections lie on the same nonidentity `IV` component.
THM-3888's local chart gives

```text
Q+.Q-=1.                                                  (33)
```

Their self-heights and mutual height are therefore

```text
<q+,q+>=<q-,q->=2-2/3=4/3,
<q+,q->=1-1-2/3=-2/3.                                   (34)
```

This is consistent with `P0=Q++Q-`, whose norm is again `4/3`.

For any `R` in the 126-element shell, Shioda's bilinear formula reads

```text
R.Qepsilon=1-<r,qepsilon>-c(r,qepsilon),                 (35)
```

where `c` is zero in the identity class, `2/3` when `r` lies in the same
nonzero discriminant class as `qepsilon`, and `1/3` in the opposite class.
Thus `(5)` asks for the following two simultaneous height pairings:

```text
class of r        required (<r,q+>,<r,q->)      count
identity                         (1,1)              0
same nonzero class               (1/3,1/3)          8
opposite nonzero class           (2/3,2/3)          1.    (36)
```

The last vector is uniquely `q++q-=p0`.  The exact companion obtains `(36)`
by projecting all 240 `E8` roots and checking every ordered pair of norm
`4/3` vectors with inner product `-2/3`.  There are `540=54*10` such ordered
pairs, and every one has exactly the same `0,8,1` profile.  The universe,
pair condition, discriminant class, local correction, and consequence
intersection numbers are all explicit rational arithmetic.  This proves

```text
#{R!=O : R.O=R.Q+=R.Q-=0}=9.                            (37)
```

## 5. Eight exact hostiles to the missing converse

THM-3888 proves the safe direction

```text
T,G in C[y]  ==>  no finite contact with Q+ union Q-.    (38)
```

It deliberately does not claim the converse from strict-transform avoidance,
because a rational inverse can acquire vertical basepoint debt.  The shell
`(37)` shows that this warning is necessary, not merely cautious.

THM-3900 independently classifies every polynomial quartic point over `C[y]`:
besides `O,P0`, there are only

```text
T_*=-2K/(3a^2),
G=+/-G_*,                G_*=4K^2/(3a^3)-L^2.            (39)
```

Only `P0` belongs to `(37)`:

- `O` is the excluded origin itself;
- the `+G_*` section is disjoint from `O` but meets each of `Q+`,`Q-` once
  at infinity—their `IV` `Y_1` addresses agree and their `X_1` slopes differ;
- the `-G_*` section meets `O` at the two simple roots of `K`, because there
  `T_*=0` and `-G_*=+L^2`.

Therefore the eight same-component sections in `(36)` have no polynomial
inverse `T,G`, despite being totally disjoint from `O,Q+,Q-`.  The implication
`(38)` is genuinely one-way.  Any converse must retain the total transform of
the rational-map base locus or an equivalent vertical exceptional sidecar;
the three strict horizontal sections do not carry enough information.

## 6. Scope and reproduction

This theorem gives a complete finite first height shell and an exact hostile
to an overstrong `S`-integrality rephrasing.  It does not enumerate sections
with `R.O>0`, impose `k[x]` descent by itself, construct a Keller pair, or
settle `JC(2)`.  THM-3897 and THM-3900 already close the polynomial `f=0`
lane by different arguments; the object-level Jacobian frontier remains the
nonzero-sidecar lane.

Reproduce with

```bash
python3 04-computation/jc2_equianharmonic_e6_marked_zero_contact_shell_thm3909.py
python3 -O 04-computation/jc2_equianharmonic_e6_marked_zero_contact_shell_thm3909.py
```

Both streams must byte-match
`05-knowledge/results/jc2_equianharmonic_e6_marked_zero_contact_shell_thm3909.out`.
