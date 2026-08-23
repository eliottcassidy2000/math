---
id: THM-3772
title: "Variable-flank three-charge rational-exact classification"
status: >
  PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.  Over an
  algebraically closed characteristic-zero field, every smooth polynomial
  Q=XA(XT)+chi(XT)+TB(XT) with a rational constant-Jacobian mate is either
  linear, a pure one-flank rational-exact near miss, or one of two dual
  affine-flank near-miss families.  In the mixed case rational exactness is
  equivalent to chi=h+gz and AB=p+g^2z/4 with p nonzero.  Every nonlinear
  case fails polynomial regularity by an explicit vertical principal-part
  mismatch.  Thus a polynomial mate exists exactly when Q is linear.  This
  does not prove JC(2), and the result has no canonical force before
  independent audit promotion.
source: root + jc_quartic_c3_construct / 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The exact companion verifies the universal
  generic quadratic equation and Jacobian sign, repeated-root eliminant,
  both mixed factorizations and rational primitives, sharp singular and
  split-fibre boundaries, pure-flank primitives and smoothness boundaries,
  all 729 linear-profile triples by direct Groebner smoothness, and bounded
  polynomial-mate hostiles.  Normal and optimized runs byte-match the frozen
  transcript; independent hostile audit remains open.
depends_on: []
related:
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3758-quadratic-radial-carrier-rational-exact-split-fibre-nonentry
  - THM-3765-normalized-three-consecutive-charge-radial-keller-classification
  - THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate
script: 04-computation/jc2_variable_flank_three_charge_classification_thm3772.py
output: 05-knowledge/results/jc2_variable_flank_three_charge_classification_thm3772.out
script_sha256: 8875a354d84a024f5e6af83aeb751a077ac72f97c2a847d482aad4453511073c
output_sha256: ebe2977416c9e8f480cf1700a6beb54367f0d3404cee80c532401ad594c01a51
semantic_sha256: 44f6673c2b293374f21a4d888a85ff37747cd546b54c42c03854a91d4088c2bb
hash_basis: raw LF bytes
---

# THM-3772 -- varying both flanks still leaves a vertical invoice

**PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.**  THM-3765
fixes the positive Euler flank to `X`.  The present theorem varies both
flanks and closes the resulting full three-consecutive-charge radial ansatz.
The generic fibre sees only the product of the two flank profiles.  Polynomial
factorization then forces one flank back to a constant on every mixed
rational-exact boundary.  The information lost by that product is recovered
by the special-fibre principal parts.

Work over an algebraically closed field `k` of characteristic zero.  Put

```text
z=XT,                       A,B,chi in k[z],
Q(X,T)=X A(z)+chi(z)+T B(z).                            (1)
```

Assume `Q` has no critical point.  Then `Q` has a rational mate
`P in k(X,T)` with `J(P,Q) in k*` exactly in the following cases.

1. **Mixed flank.**  Both `A` and `B` are nonzero, and for some
   `h,g,p in k` with `p!=0`,

   ```text
   chi=h+gz,                    A(z)B(z)=p+(g^2/4)z.   (2)
   ```

   If `g=0`, both flanks are constant and `Q` is linear.  If `g!=0`,
   polynomial factorization gives exactly two orientations:

   ```text
   A=a,       B=(p+g^2z/4)/a,       a in k*,          (3a)
   B=b,       A=(p+g^2z/4)/b,       b in k*.          (3b)
   ```

2. **Pure positive flank.**  `B=0`, `chi=h` is constant, and

   ```text
   A(0)!=0,                     gcd(A,A')=1.           (4a)
   ```

3. **Pure negative flank.**  `A=0`, `chi=h` is constant, and

   ```text
   B(0)!=0,                     gcd(B,B')=1.           (4b)
   ```

In all three cases every rational mate is written explicitly below.  A
polynomial constant-Jacobian mate exists if and only if `Q` is a nonconstant
linear polynomial.  Thus neither varying the normalized flank nor passing
through the pure-flank boundary produces a planar Jacobian counterexample.

## 1. The product quotient and the exact time form

Let `Lambda` be transcendental over `k`, set

```text
F=zAB,                  W=Lambda-chi,
Y=XA-TB.                                                 (5)
```

On the generic fibre `Q=Lambda`, one has

```text
W=XA+TB,
Y^2=Delta(z),           Delta=(Lambda-chi)^2-4F,        (6)
X=(W+Y)/(2A),           T=(W-Y)/(2B),                   (7)
```

whenever both flanks are nonzero.  The reconstruction in `(7)` takes place
in the function field and its product is `z`, because `W^2-Y^2=4zAB`.
Direct differentiation gives the decisive sign

```text
J(Q,z)=Y.                                                (8)
```

Consequently a rational mate with `J(P,Q)=c in k*` must satisfy

```text
dP=-c dz/Y                                               (9)
```

on the generic fibre.

Put `G=chi'`.  The repeated-root eliminant of `(6)` is

```text
R=F'^2-FG^2.                                             (10)
```

For completeness, smoothness forces `R!=0`.  If `R` vanished identically,
choose `alpha in k` away from the finitely many zeros of `zABG`.  Then

```text
W_0=-2F'(alpha)/G(alpha),       W_0^2=4F(alpha),        (11)
X_0=W_0/[2A(alpha)],            T_0=W_0/[2B(alpha)]     (12)
```

give `X_0T_0=alpha`.  In the etale chart `(X,z)`, equations `(11)` say
`Y=0`.  On that locus the exact identity

```text
XA (Q_z|X)=XA G+F'                                   (12a)
```

turns `(11)` into `Q_z|X=0`, so `(X_0,T_0)` is a critical point.  The only
apparent exception `G=0` would make `F'=0`; since `F(0)=0`, this makes
`F=0`, contrary to the mixed hypothesis.

It follows that `Delta` is squarefree in `k(Lambda)[z]`.  Indeed, a common
root `alpha` of `Delta,Delta'` is a root of the fixed nonzero polynomial
`R`, hence `alpha` is algebraic over `k` and belongs to `k`.  But
`W(alpha)^2=4F(alpha)` then makes `W(alpha)` algebraic over `k`, contradicting
the transcendence of `Lambda=chi(alpha)+W(alpha)`.  Thus `(6)` is a
geometrically integral quadratic generic fibre and its rational constant
field is exactly `k(Lambda)`.  All rational mates differ by `k(Q)`.

## 2. Rational exactness forces one constant flank

Let `d=deg_z Delta`.  If `d>=3`, `dz/Y` is a nonzero holomorphic
differential on the smooth projective hyperelliptic model of `(6)`, so it is
not a rational derivative.  If `d=2`, its two infinity residues are
`+/-1/sqrt(lc(Delta))`, again nonzero.  Therefore `(9)` requires

```text
deg_z Delta<=1.                                         (13)
```

The coefficient of `Lambda` in `Delta` first forces

```text
chi=h+gz.                                               (14)
```

Because `F(0)=0`, cancellation of every quadratic and higher coefficient
then gives exactly

```text
F=pz+(g^2/4)z^2,
A(z)B(z)=p+(g^2/4)z.                                   (15)
```

Conversely `(14),(15)` make `Delta` linear and integrate `(9)`.  Define

```text
D=g(Q-h)+2p.                                            (16)
```

At fixed `Q`, `Delta_z=-2D`, hence `Y_z=-D/Y`, and

```text
P_0=cY/D,                    J(P_0,Q)=c.                (17)
```

This proves rational existence and completeness.

The right side of `(15)` has degree at most one.  If it is nonconstant, the
domain `k[z]` forces one of `A,B` to be a unit, giving `(3a)` or `(3b)`.
If it is a nonzero constant, both are units.  The only remaining product-zero
case belongs to the pure analysis in Section 4.

## 3. Smoothness and the mixed vertical mismatch

In orientation `(3a)`, put

```text
s=1+gT/(2a).
```

Then

```text
Q=h+(p/a)T+aXs^2,
Q_X=as^2,                    Q_T=p/a+gXs.              (18)
```

The dual formulas hold in `(3b)` after swapping `X,T` and `a,b`.  If
`g!=0,p=0`, the line `s=0` is critical.  If `p!=0`, any zero of `Q_X`
leaves `Q_T=p/a`, so no critical point exists.  If `g=0`, the mixed product
is the nonzero constant `p`, both flanks are constant, and `Q` is linear.
Thus the exact mixed smoothness boundary is `p!=0`.

Suppose now `g,p!=0`.  Put

```text
lambda_0=h-2p/g,             w_0=lambda_0-h=-2p/g.     (19)
```

In orientation `(3a)` the exceptional fibre factors as

```text
Q-lambda_0
 =s[2p/g+aXs].                                          (20)
```

The two factors are distinct irreducibles: the first is linear in `T`; the
second is primitive linear in `X` over `k[s]` with nonzero constant term.
They are reduced because `Q` is smooth.  On them the numerator in `(17)`
takes the values

```text
s=0:                 Y=-w_0=+2p/g,
2p/g+aXs=0:          Y= w_0=-2p/g.                    (21)
```

Since `D=g(Q-lambda_0)`, the simple principal coefficients of `P_0` are
`+2cp/g^2` and `-2cp/g^2`.  A target correction `H(Q)` has the same Laurent
coefficient on both components, so it cannot cancel both.  A higher pole
only creates an uncancellable higher pole.  The dual orientation has the
same two opposite values with the component labels exchanged.  Hence no
nonlinear mixed member has a polynomial mate.

If `g=0`, `(15)` reads `AB=p in k*`, so `A=a,B=b` are constants and

```text
Q=h+aX+bT,
P=c(aX-bT)/(2ab)                                      (22)
```

is a polynomial mate.  This completes the mixed classification without
using the pending general equalizer theorem THM-3770.

## 4. Pure flanks are rational-exact but vertically unequalized

Let `B=0` and `A!=0`.  On the generic fibre,

```text
Y=XA=Lambda-chi(z),
dP=-c dz/[Lambda-chi(z)].                              (23)
```

If `chi` is nonconstant, the generic denominator is separable and every root
has a nonzero logarithmic residue.  Thus rational exactness forces `chi=h`
constant.  Conversely

```text
Q=h+XA(z),                  P_0=-cT/A(z)               (24)
```

satisfy `J(P_0,Q)=c`.  The source chart `(X,z)` shows that `Q` is smooth
exactly when

```text
A(0)!=0,                     gcd(A,A')=1:              (25)
```

the origin/`X=0` axis gives the first condition, while a torus critical point
is exactly a common root of `A,A'`.  The rational constant field is
`k(Q)`, because the generic fibre is the rational line with coordinate `z`.

If `A` is nonconstant under `(25)`, choose any root `alpha` of `A`.  It is
nonzero and simple.  The fibre

```text
Q-h=XA(z)                                               (26)
```

has the component `X=0` and the distinct component `z-alpha=0`.  The
primitive `(24)` is regular along `X=0` and has a simple pole along
`z-alpha=0`.  Any `H(Q)` regular at `h` cancels neither pole.  Any `H(Q)`
with a pole at `h` creates a pole along `X=0`, where `(24)` has none.  Thus
the complete torsor `P_0+k(Q)` contains no polynomial.  If `A` is constant,
`Q` is linear and `(24)` is already polynomial.

The case `A=0,B!=0` is the exact dual:

```text
Q=h+TB(z),                  P_0=cX/B(z),
B(0)!=0,                    gcd(B,B')=1.               (27)
```

It has a polynomial mate exactly when `B` is constant.  Finally, if both
flanks vanish, `Q=chi(XT)` is critical at the origin and never enters the
smooth classification.  Sections 1--4 prove the theorem.  **QED conditional
on audit promotion.**

## 5. Exact controls and the design consequence

Reproduce with

```bash
python3 -B 04-computation/jc2_variable_flank_three_charge_classification_thm3772.py
python3 -B -O 04-computation/jc2_variable_flank_three_charge_classification_thm3772.py
```

The assertion-free companion verifies the universal identities, both mixed
orientations and their opposite vertical values, the pure mate signs and
smoothness boundaries, and all `729` triples of linear profiles over
`{-1,0,1}`.  Direct affine Groebner ideals find `48` smooth members; all `48`
are rational-exact, and every mixed smooth member has a squarefree generic
quadratic.  Twenty finite polynomial-mate systems guard the two nonlinear
mechanisms through total mate degree nine.  These are hostile controls; the
proof is all-degree.

The construction lesson is sharper than “add another charge.”  The generic
quadratic forgets the ordered pair `(A,B)` and retains only `AB`; rational
exactness makes that product linear, so unique factorization collapses one
flank.  On the special fibre, the forgotten ordering returns as two opposite
principal values.  The next viable design must therefore leave this
one-product quadratic quotient, or add a response channel capable of
distinguishing its vertical components.
