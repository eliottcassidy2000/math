---
id: THM-3074
title: "C3 two-pole binomial cancellation and first key-form depth lattice"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let a polynomial
  planar Keller pair have a tame ramification-index-three, residue-degree-one
  valuation over the target coordinate line t=0.  If both affine source
  coordinates have poles p,q>=1, their first wedge coefficient must vanish;
  with h=gcd(p,q), the leading coefficients consequently satisfy the
  primitive binomial Y^(p/h)-cX^(q/h)=0.  A unimodular Laurent change gives a
  value-zero toric coordinate M, a complementary coordinate R of value h,
  and Z=M/c-1 of depth ell.  Necessarily 1<=ell<=D=p+q+3.  If ell=D, its
  first logarithmic wedge coefficient is exactly the Keller coefficient; if
  ell<D, that coefficient vanishes and m^h/R_0^ell is constant, forcing a
  second key-form resonance.  Together with THM-3070, the same result holds
  for every one- or two-coordinate escape branch in this coordinate-line
  scope, with D=p+3-r in the one-pole case.  If a polynomial target is
  detected without cancellation by the first (R,Z)-data, its valuation lies
  in hZ+ellZ; hence detection of t of value three implies gcd(h,ell)|3.
  This is not a claim about the full value semigroup: an exact local
  symplectic hostile with h=2, ell=4 has the polynomial x(y-x)-1 of value
  three only after its first-lattice terms cancel.  No full C3, A4/S4, or
  planar Jacobian exclusion is asserted.
source: codex-jc-resolvent-bridge-2026-08-01
depends_on:
  - THM-3070-polynomial-c3-one-face-escape-leading-cancellation-gate
related:
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift
  - HYP-9070-jc2-leading-form-circuit-and-the-euclidean-depth-search-order
script: 04-computation/jc_c3_two_pole_first_key_lattice_thm3074.py
output: 05-knowledge/results/jc_c3_two_pole_first_key_lattice_thm3074.out
script_sha256: 5bfc52f67d25e9c53dd2d951ecd4ae6e20fb03a21bc0bbfefb110defda2db757
output_sha256: 2b760bdb62b8a81bd67fd408da83f40127112788f983b3b639b20b48094225d4
hash_basis: LF-normalized bytes
---

# THM-3074 -- two-pole C3 escape starts a toric key-form tower

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Inheritance and result

[THM-3070](THM-3070-polynomial-c3-one-face-escape-leading-cancellation-gate.md)
proves that a tame `C3` valuation over a target coordinate line cannot have
exactly one source coordinate escape through a leading-wedge-nondegenerate
face.  Every surviving one-pole branch already has a primitive binomial
relation in its associated graded source.  The remaining elementary escape
geometry is the **two-pole case**, where both affine source coordinates tend
to infinity.

This theorem proves three statements in that remaining geometry.

1. The leading two-pole face always cancels and is primitive binomial.
2. A unimodular Laurent chart turns the next Jacobian obligation into an
   exact depth invoice.  Either the first key pays at its maximal permitted
   depth, or it too is resonant and a later key is compulsory.
3. The first key sees only the lattice `h Z+ell Z`.  This yields a useful
   divisibility gate for an **uncancelled first-key detection**, but not for
   the full value semigroup.  An exact hostile isolates the distinction.

The coordinate-line hypothesis is load bearing.  Nothing below straightens
an arbitrary Jelonek component by a polynomial constant-Jacobian target
change, and nothing globalizes the local Laurent chart.

## 2. Coordinate-line setup

Let

```text
P,Q in C[x,y],                    u=P, t=Q,
det partial(P,Q)/partial(x,y)=kappa in C*.               (1)
```

Let `w` be a divisorial valuation of `C(x,y)` above `(t)` satisfying

```text
w(t)=3,                  w(C(u)*)=0,
residue_field(w)=C(u).                                   (2)
```

Choose a coefficient-field uniformizer `s` in the completed valuation field,
so that

```text
t=tau(u)s^3+O(s^4),                  tau in C(u)*.        (3)
```

The Keller identity is then

```text
dx wedge dy=kappa^(-1)du wedge dt
 =3 kappa^(-1)tau(u)s^2 du wedge ds+O(s^3).              (4)
```

This theorem treats the two-pole expansion

```text
x=A(u)s^(-p)+O(s^(-p+1)),
y=B(u)s^(-q)+O(s^(-q+1)),             p,q>=1,            (5)
```

with `A,B in C(u)*`.  Section 6 states the uniform version obtained by
combining it with the one-pole result of THM-3070.

## 3. The two-pole leading face must be binomial

Differentiating `(5)` gives the first possible physical wedge term

```text
dx wedge dy
 =(p A B'-q A'B)s^(-p-q-1)du wedge ds
   +O(s^(-p-q))du wedge ds.                              (6)
```

Its exponent is strictly below two.  Comparison with `(4)` forces

```text
p A B'-q A'B=0,
(B^p/A^q)'=0.                                           (7)
```

Put

```text
h=gcd(p,q),             a=p/h,             b=q/h.       (8)
```

The constants of `C(u)` are `C`, and a rational function whose `h`-th power
is constant is itself constant.  Thus `(7)` is equivalent to

```text
c:=B^a/A^b in C*.                                      (9)
```

If `X=A s^(-p)` and `Y=B s^(-q)` denote the leading forms, `(9)` is the
primitive toric binomial

```text
Y^a-c X^b=0.                                           (10)
```

This is a necessary associated-graded relation, not a polynomial identity
between `x` and `y`.

## 4. Unimodular toric coordinates and the depth invoice

Define the value-zero Laurent monomial

```text
M=y^a/x^b.                                             (11)
```

Choose integers `c_1,d_1` with

```text
a c_1+b d_1=-1                                         (12)
```

and put

```text
R=x^(c_1)y^(d_1).                                      (13)
```

The exponent matrix

```text
[ -b    a ]
[ c_1 d_1 ]
```

has determinant one.  Consequently

```text
C[x,x^(-1),y,y^(-1)]=C[M,M^(-1),R,R^(-1)],
w(M)=0,                  w(R)=h,
dlog M wedge dlog R=dlog x wedge dlog y.                (14)
```

Write

```text
M=c(1+m(u)s^ell+O(s^(ell+1))),
R=R_0(u)s^h(1+O(s)),
m,R_0 in C(u)*.                                        (15)
```

The integer `ell` is finite and positive: if `M=c` identically, then `x,y`
would satisfy a nontrivial algebraic relation and `dx wedge dy` would vanish.
Set

```text
Z=M/c-1,                    ell=w(Z),
D=p+q+3.                                                (16)
```

From `(15)`,

```text
dlog M wedge dlog R
 =[h m'-ell m R_0'/R_0]s^(ell-1)du wedge ds
   +O(s^ell)du wedge ds.                               (17)
```

Since `dx wedge dy=xy dlog M wedge dlog R` and
`xy=ABs^(-p-q)(1+O(s))`, the first possible exponent in `(17)` is

```text
ell-p-q-1=ell+2-D.                                     (18)
```

Equation `(4)` now gives the exact trichotomy

```text
1<=ell<=D;                                             (19)

ell=D  ==>  AB(h m'-ell m R_0'/R_0)
                    =3 kappa^(-1)tau;                  (20)

ell<D  ==>  h m'-ell m R_0'/R_0=0
         ==>  m^h/R_0^ell in C*.                       (21)
```

Indeed, `ell>D` would put every term of the physical wedge above order two.
At `ell=D`, the displayed coefficient is the unique possible order-two
coefficient and cannot vanish.  At `ell<D`, it lies below order two and must
vanish.  The last implication in `(21)` is logarithmic differentiation.

Thus a strict-depth first key is itself resonant: its leading coefficient is
a rational power of the leading coefficient of `R`, and further Laurent
data must supply the Jacobian.  This is the precise, proved content of
“a second key form is required.”  No termination or Euclidean bound for the
resulting tower is asserted.

## 5. The first-key lattice and its exact scope

The unimodularity in `(14)` gives a clean test for what the pair `(R,Z)` can
detect before any further cancellation.  Let `F in C[x,y]` be nonzero and
write its unique Laurent expansion as

```text
F=sum_(j in J) R^j f_j(M),
J finite,                    f_j in C[M,M^(-1)].        (22)
```

For every nonzero `f_j`, expand at the nonzero point `M=c`:

```text
f_j(c(1+Z))=alpha_j Z^(n_j)+O(Z^(n_j+1)),
alpha_j in C*,                       n_j>=0.             (23)
```

Define the first-key predicted value

```text
mu(F)=min_(j in J)(j h+n_j ell).                       (24)
```

The coefficient predicted at that value is

```text
I_(R,Z)(F)=sum_(j h+n_j ell=mu(F))
              alpha_j R_0^j m^(n_j) in C(u).           (25)
```

Call `F` **first-key detected** when `I_(R,Z)(F)` is nonzero.  Direct
substitution in `(22)--(23)` then proves

```text
F first-key detected  ==>  w(F)=mu(F)
                             in h Z+ell Z.              (26)
```

In particular, if the target polynomial `Q=t` is first-key detected, then

```text
3=w(Q) in h Z+ell Z=gcd(h,ell)Z,
gcd(h,ell) divides 3.                                  (27)
```

The contrapositive is often more useful:

```text
gcd(h,ell) does not divide 3
 ==> the first-key initial coefficient of Q cancels.   (28)
```

Neither `(26)` nor `(28)` says that every later value belongs to this lattice.
Subleading coefficients of `R` and `Z`, or cancellation of `(25)`, can expose
arbitrary later powers of `s`.  Section 8 gives an exact counterexample to
the stronger value-semigroup reading.

## 6. Uniform one-pole/two-pole statement

Suppose instead that exactly one source coordinate escapes, in the notation
of THM-3070:

```text
x=A s^(-p)+...,
y=B s^r+...,                         p>=1, r>=0.         (29)
```

That theorem proves

```text
r<=p+2,
r A'B+p AB'=0.                                           (30)
```

Put `h=gcd(p,r)`, `a=r/h`, and `b=p/h`; this includes `r=0`, when
`(a,b)=(0,1)`.  Then

```text
M=x^a y^b,                 w(M)=0,
leading relation X^a Y^b-c=0.                          (31)
```

Choose `c_1,d_1` with `a d_1-b c_1=1` and let
`R=x^(c_1)y^(d_1)`.  Again the exponent matrix has determinant one,
`w(R)=h`, and the calculation `(15)--(28)` is unchanged after replacing

```text
D=p+q+3       by       D=p+3-r.                        (32)
```

Therefore every escaping branch covered by the coordinate-line setup--one
or both affine source coordinates having poles--starts with a primitive
binomial cancellation and obeys the same first-key depth/lattice gate.

## 7. Equality-depth two-pole control

The maximal-depth lane in `(20)` is nonempty at the exact local symplectic
level.  On the ramified cover `t=s^3`, set

```text
x=s^(-1),
y=s^(-1)+3u s^4.                                       (33)
```

Then

```text
dx wedge dy=3s^2 du wedge ds=du wedge dt.              (34)
```

Here `p=q=h=1`, `D=5`, and one can take

```text
M=y/x=1+3u s^5,                  R=1/x=s,
ell=5=D.                                                (35)
```

In fact `(33)` is the inverse of the punctured rational Keller pair

```text
P=x^4(y-x)/3,                 Q=x^(-3),
det partial(P,Q)/partial(x,y)=1.                        (36)
```

Thus the equality-depth differential invoice is sharp.  The second target
coordinate in `(36)` is Laurent, not polynomial; polynomiality of `Q` remains
a load-bearing global condition.

## 8. Strict-depth hostile to the value-semigroup overclaim

There is also an exact strict-depth packet showing why `(27)` cannot be
promoted to a statement about all polynomial values.  Put

```text
R=u s^2,
M=1+R^2-3u^3 s^7,
x=R^(-1),                    y=M R^(-1).                (37)
```

Direct differentiation gives

```text
dx wedge dy=3s^2 du wedge ds=du wedge d(s^3).          (38)
```

Both `x` and `y` have pole order two, so

```text
p=q=2,                  h=2,             D=7,
Z=M-1=u^2s^4-3u^3s^7,  ell=4<D.                        (39)
```

The first coefficient is `m=u^2`, while `R_0=u`, and indeed

```text
h m'-ell m R_0'/R_0=2(2u)-4u^2/u=0.                   (40)
```

Thus the first lattice is `2Z`, which does not contain three.  Nevertheless
the honest polynomial

```text
T(x,y)=x(y-x)-1                                        (41)
```

satisfies exactly

```text
T=Z/R^2-1=-3u s^3,
w(T)=3.                                                (42)
```

At first-key level the two weight-zero contributions in `(42)` are
`m/R_0^2=1` and `-1`; they cancel.  The order-seven correction in `Z`, after
division by `R^2` of order four, exposes the off-lattice value three.  This
packet is a local rational symplectic hostile, not a polynomial Keller pair
with both target coordinates.  Its role is exactly to refute the invalid
implication

```text
first lattice dZ   ==>   full polynomial value semigroup lies in dZ.       (43)
```

## 9. Exact verification and boundary

The companion script, replayed byte-identically under normal and optimized
Python during an independent line audit, checks symbolically:

- the two-pole wedge coefficient `(6)`;
- the two unimodular Bezout charts, including noncoprime orders and `r=0`;
- the arithmetic identity `hZ+ellZ=gcd(h,ell)Z` on a finite hostile grid;
- the unit Jacobian and depth-five equality packet `(33)--(36)`;
- every identity in the strict-depth/off-lattice packet `(37)--(42)`; and
- the one-pole equality packet from THM-3070 as an independent cross-control.

Reproduce with

```bash
python3 04-computation/jc_c3_two_pole_first_key_lattice_thm3074.py
python3 -O 04-computation/jc_c3_two_pole_first_key_lattice_thm3074.py
```

The theorem is confined to a residue-degree-one tame `C3` valuation over the
actual target coordinate line `Q=t=0`.  It supplies necessary binomial,
depth, and first-detection constraints.  It does not prove that the key-form
tower terminates, that later values remain in the first lattice, that an
arbitrary Jelonek curve can be polynomially straightened, or that `C3`,
`A4/S4`, or the planar Jacobian conjecture is excluded.
