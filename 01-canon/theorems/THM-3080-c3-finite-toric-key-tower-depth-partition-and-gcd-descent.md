---
id: THM-3080
title: "C3 finite toric key-tower depth partition and gcd descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let K=C(x,y) carry a
  normalized residue-degree-one divisorial valuation w with w(K*)=Z, and
  suppose the initial toric coordinates generate the full field,
  C(M_0,R_0)=K.  If a nonzero two-form has a
  toric presentation Omega=U_i dlog(M_i) wedge dlog(R_i), where M_i has
  constant nonzero residue, R_i has positive value g_i, and Omega is the
  inverse Keller form above a tame target parameter of ramification index E.
  If e_i is the value of M_i/res(M_i)-1 and B_i=E-w(U_i), then
  1<=e_i<=B_i.  Equality is terminal and pays the exact target-wedge
  coefficient.  Strict inequality forces a primitive binomial cancellation;
  a unimodular Laurent normalization produces the next toric pair with
  B_(i+1)=B_i-e_i and g_(i+1)=gcd(g_i,e_i).  Hence the differential
  cancellation tower always terminates, its positive depths partition B_0,
  and it has at most B_0 stages.  At the terminal stage the nonzero wedge
  makes the primitive coefficient ratio nonconstant, so no equal-weight
  Laurent initial sum can cancel.  Normalization of the divisorial value
  group then forces gcd(g_N,e_N)=1.  For the coordinate-line C3 branches of
  THM-3074 this gives sum e_i=p+q+3 in the two-pole case and sum
  e_i=p+3-r in the one-pole case, with a terminal primitive lattice.  The
  normalized keys are Laurent key forms in the completed field, not asserted
  polynomial MacLane keys; no polynomial globalization, arbitrary-Jelonek
  straightening, full C3, A4/S4, or JC(2) exclusion follows.
source: codex-jc-resolvent-bridge-2026-08-01
depends_on:
  - THM-3074-c3-two-pole-binomial-cancellation-and-first-key-form-depth-lattice
related:
  - THM-3070-polynomial-c3-one-face-escape-leading-cancellation-gate
  - HYP-9070-jc2-leading-form-circuit-and-the-euclidean-depth-search-order
script: 04-computation/jc_c3_finite_toric_key_tower_thm3080.py
output: 05-knowledge/results/jc_c3_finite_toric_key_tower_thm3080.out
script_sha256: f813fefec6e81ff1ba2db767a9699929f1d9fab0d5514e1c9da1c742cc7072a9
output_sha256: 9a07e8b2fc7292044b1e6e0d837a6f7d6fb4db003f1a16772d319f6bc975b907
hash_basis: LF-normalized bytes
---

# THM-3080 -- the C3 toric Jacobian tower has a finite depth budget

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. The missing step after the first key

[THM-3074](THM-3074-c3-two-pole-binomial-cancellation-and-first-key-form-depth-lattice.md)
proves that every one- or two-coordinate escape branch in its precise
coordinate-line `C3` scope begins with a primitive toric relation.  It also
proves a first-key dichotomy: the first deviation either supplies the
Jacobian at one exact depth or has another primitive leading cancellation.

Iterating that calculation reveals a finite invariant which was not visible
at the first stage.  Each strict cancellation multiplies the two-form
prefactor by the current deviation.  Its positive valuation is therefore
spent permanently.  The remaining Jacobian budget decreases by exactly that
amount, so an infinite cancellation tower is impossible.

The one-step calculation is independent of the number three.  We first prove
it for any tame ramification index `E`, then specialize to `C3`.  This
separates two roles:

- `E` controls the transverse differential budget;
- gcd and Bezout normalization control the Euclidean lattice descent.

The objects constructed below are **Laurent toric key forms in the completed
function field**.  Calling them keys does not assert polynomiality, the full
MacLane axioms, or algebraization on `A^2`.

## 2. Abstract ramified toric stage

Work in characteristic zero in the completed divisorial field of

```text
K=C(x,y),                    w(K*)=Z,                    (0)
```

with coefficient field `C(u)` and uniformizer `s`.  Require at stage zero

```text
C(M_0,R_0)=K.                                            (0a)
```

The unimodular transformations below preserve this full-field equality.
This condition is load bearing; algebraic independence alone would allow a
proper finite-index subfield and would not imply terminal primitivity.

Let a tame target parameter satisfy

```text
t=tau(u)s^E+O(s^(E+1)),              E>=1,
tau in C(u)*.                                             (1)
```

For a planar Keller pair of Jacobian `kappa in C*`, its inverse two-form is

```text
Omega:=dx wedge dy=kappa^(-1)du wedge dt
 =E kappa^(-1)tau s^(E-1)du wedge ds+O(s^E).              (2)
```

Suppose at some stage `i` there are full-field coordinates `M_i,R_i`, so
`C(M_i,R_i)=K`, and a nonzero prefactor `U_i` such that

```text
Omega=U_i dlog(M_i) wedge dlog(R_i),                     (3)

U_i=u_i(u)s^(sigma_i)(1+O(s)),
M_i=c_i(1+m_i(u)s^(e_i)+O(s^(e_i+1))),
R_i=r_i(u)s^(g_i)(1+O(s)),                              (4)

u_i,m_i,r_i in C(u)*,             c_i in C*,
e_i,g_i>=1.                                               (5)
```

Set

```text
Z_i=M_i/c_i-1,
B_i=E-sigma_i.                                           (6)
```

The number `B_i` is the remaining differential budget.

## 3. One-step depth dichotomy

Direct logarithmic differentiation of `(4)` gives

```text
Omega
 =u_i[g_i m_i'-e_i m_i r_i'/r_i]
    s^(sigma_i+e_i-1)du wedge ds
   +O(s^(sigma_i+e_i))du wedge ds.                       (7)
```

Comparison with `(2)` proves

```text
1<=e_i<=B_i.                                             (8)
```

There are exactly two cases.

### Terminal case

If `e_i=B_i`, the displayed exponent in `(7)` is `E-1`.  No other term can
pay the leading target coefficient, so

```text
u_i[g_i m_i'-e_i m_i r_i'/r_i]
              =E kappa^(-1)tau !=0.                     (9)
```

Call this a **terminal stage**.

### Strict case

If `e_i<B_i`, the displayed exponent is below `E-1` and its coefficient must
vanish:

```text
g_i m_i'-e_i m_i r_i'/r_i=0.                            (10)
```

Put

```text
d_i=gcd(g_i,e_i),
alpha_i=g_i/d_i,                  beta_i=e_i/d_i.        (11)
```

Equation `(10)` says

```text
c_(i+1):=m_i^(alpha_i)/r_i^(beta_i) in C*.              (12)
```

The constants are again exactly `C`; taking primitive powers in `(11)` is
why no root-choice ambiguity remains in `(12)`.

## 4. Unimodular normalization of a strict stage

Choose integers `gamma_i,delta_i` satisfying

```text
alpha_i delta_i+beta_i gamma_i=1.                        (13)
```

Define

```text
M_(i+1)=Z_i^(alpha_i)/(c_(i+1) R_i^(beta_i)),
R_(i+1)=Z_i^(gamma_i) R_i^(delta_i).                    (14)
```

The exponent matrix in the ordered coordinates `(Z_i,R_i)` is

```text
[ alpha_i  -beta_i ]
[ gamma_i   delta_i ]                                   (15)
```

and has determinant one.  Therefore

```text
C(M_i,R_i)=C(M_(i+1),R_(i+1)),
dlog(M_(i+1)) wedge dlog(R_(i+1))
                  =dlog(Z_i) wedge dlog(R_i),           (16)

w(M_(i+1))=0,             residue(M_(i+1))=1,
w(R_(i+1))=d_i.                                        (17)
```

The new value-zero coordinate is not identically one: otherwise `(14)`
would be an algebraic relation between the independent coordinates
`Z_i,R_i`.  Consequently there is a finite positive integer `e_(i+1)` with

```text
M_(i+1)=1+m_(i+1)s^(e_(i+1))+O(s^(e_(i+1)+1)).          (18)
```

Since

```text
dlog(M_i)=Z_i/(1+Z_i) dlog(Z_i),                         (19)
```

equations `(3)` and `(16)` give the next toric presentation

```text
Omega=U_(i+1)dlog(M_(i+1)) wedge dlog(R_(i+1)),
U_(i+1)=U_i Z_i/(1+Z_i).                                (20)
```

Its two decisive ledgers are

```text
sigma_(i+1)=sigma_i+e_i,
B_(i+1)=B_i-e_i,                                        (21)

g_(i+1)=d_i=gcd(g_i,e_i).                               (22)
```

Thus strict depth is consumed while the complementary value undergoes gcd
descent.

The corresponding unnormalized primitive key form is

```text
K_(i+1)=Z_i^(alpha_i)-c_(i+1)R_i^(beta_i).              (23)
```

From `(14)` and `(18)`, its exact value is

```text
w(K_(i+1))
 =lcm(g_i,e_i)+e_(i+1).                                 (24)
```

So `e_(i+1)` is precisely the excess above the primitive binomial's common
monomial value.

## 5. Finite termination and the depth partition

Begin at stage zero and repeat Section 4 whenever the strict case occurs.
At every strict step, `(8)` gives

```text
1<=e_i<B_i,
B_(i+1)=B_i-e_i>=1.                                    (25)
```

The positive integer budget strictly decreases.  Hence after finitely many
strict steps there is a terminal stage `N`.  Iterating `(21)` and using
`e_N=B_N` yields the exact partition

```text
e_0+e_1+...+e_N=B_0=E-sigma_0.                          (26)
```

In particular,

```text
number of stages=N+1<=B_0,                              (27)

g_(i+1)=gcd(g_i,e_i) divides g_i.                       (28)
```

This proves termination of the **differential toric cancellation tower**.
It does not say that a global polynomial approximate-root algorithm
terminates, because the transformations `(14)` allow negative powers and
are made only in the completed function field.

### 5.1 The terminal lattice is primitive

There is one more consequence at the terminal stage.  Write

```text
g=g_N,               e=e_N,               d=gcd(g,e),
G=g/d,                H=e/d,
theta=m_N^G/r_N^H in C(u)*.                              (29)
```

The terminal coefficient `(9)` is nonzero, so

```text
theta'/theta
 =(g m_N'/m_N-e r_N'/r_N)/d !=0.                        (30)
```

Thus `theta` is a nonconstant rational function of `u`, hence
transcendental over `C`.

Let `F` be a nonzero Laurent polynomial in `M_N,R_N`.  Group it uniquely as

```text
F=sum_j R_N^j f_j(M_N),          f_j in C[M_N,M_N^(-1)]. (31)
```

If

```text
f_j(c_N(1+Z_N))=a_j Z_N^(n_j)+higher,
a_j in C*,                                           (32)
```

then its least predicted weight is the minimum of `jg+n_j e`.  For two
pairs on the same least-weight line,

```text
(j-j_0)g+(n_j-n_0)e=0
```

implies, for a unique integer `k`,

```text
j-j_0=Hk,                  n_j-n_0=-Gk.                  (33)
```

After factoring one nonzero leading monomial, the complete least-weight
coefficient is therefore a nonzero Laurent polynomial in

```text
r_N^H/m_N^G=theta^(-1).                                 (34)
```

It cannot vanish because `theta` is transcendental over `C`.  Consequently

```text
w(F)=min_j(jg+n_j e) in dZ                              (35)
```

for every nonzero Laurent polynomial `F` in the terminal coordinates.
Since `M_N,R_N` generate `C(x,y)`, every nonzero field element is a quotient
of two such Laurent polynomials.  Hence the entire value group would lie in
`dZ`.  The divisorial valuation was normalized to have value group `Z`, so

```text
gcd(g_N,e_N)=1.                                         (36)
```

This is precisely where the early-lattice hostile of THM-3074 ceases to be
an obstruction: early resonant ratios are constant and permit off-lattice
values after cancellation; the terminal ratio is nonconstant and makes the
last lattice monomially injective.

The full-field hypothesis in `(0a)` cannot be weakened to algebraic
independence.  In `C(u)(s)`, the proper-subfield packet

```text
M=1+u s^2,              R=s^2,              U=M
```

has `U dlog(M) wedge dlog(R)=du wedge d(s^2)` and terminal data
`(E,g,e)=(2,2,2)`.  Its gcd is two because
`C(M,R)=C(u,s^2)` is the index-two subfield, rather than the full normalized
field `C(u,s)`.  This exact hostile isolates the missing hypothesis without
touching the `C3` specialization, whose unimodular charts generate `K`.

## 6. Coordinate-line C3 specialization

For the `C3` branches of THM-3074, take `E=3`.  Its initial unimodular chart
has

```text
Omega=xy dlog(M_0) wedge dlog(R_0),                     (37)
```

so `U_0=xy`.

In the two-pole case

```text
w(x)=-p,             w(y)=-q,
sigma_0=-(p+q),       B_0=p+q+3.                        (38)
```

In the one-pole case inherited from THM-3070,

```text
w(x)=-p,             w(y)=r,
sigma_0=r-p,          B_0=p+3-r.                        (39)
```

Therefore the respective key-depth partitions are

```text
two poles:       e_0+...+e_N=p+q+3,
one pole:        e_0+...+e_N=p+3-r,
both cases:      gcd(g_N,e_N)=1.                        (40)
```

The first summand `e_0` is THM-3074's `ell`; `(40)` closes the previously
unbounded **local differential** depth lane in that exact scope.  It neither
bounds global polynomial degrees nor proves that every Laurent key in the
partition is represented by a polynomial in `x,y`.

## 7. Exact one-, two-, and three-stage controls

All controls use `t=s^3` and satisfy

```text
dx wedge dy=3s^2 du wedge ds=du wedge dt.               (41)
```

### 7.1 One stage: `5`

The equality packet from THM-3074 is

```text
x=s^(-1),
y=s^(-1)+3u s^4.                                       (42)
```

Here `p=q=1`, `B_0=5`, and

```text
R=s,                M=y/x=1+3u s^5,
(e_0)=(5).                                                (43)
```

### 7.2 Two stages: `4+3`

Put

```text
R=u s^2,
M_1=1-3u s^3,
M_0=1+R^2 M_1,
x=R^(-1),                 y=M_0 R^(-1).                 (44)
```

Then `p=q=2`, `B_0=7`, and

```text
Z_0=M_0-1=R^2M_1,
M_1=Z_0/R^2,
(e_0,e_1)=(4,3).                                        (45)
```

This is THM-3074's off-first-lattice hostile written in its normalized
second-key coordinates.

### 7.3 Three stages: `4+4+3`

A new exact packet shows that more than one strict normalization can really
occur.  Put

```text
R=u s^4,
M_2=1+3u s^3,
M_1=1+R M_2,
M_0=1+R M_1,
x=R^(-1),                 y=M_0 R^(-1).                 (46)
```

Both source coordinates have pole order four, so `B_0=11`.  The normalized
keys are exact:

```text
M_1=(M_0-1)/R,
M_2=(M_1-1)/R,
(e_0,e_1,e_2)=(4,4,3),          4+4+3=11.               (47)
```

Moreover

```text
dx wedge dy
 =R^(-3)dM_0 wedge dR
 =R^(-2)dM_1 wedge dR
 =R^(-1)dM_2 wedge dR
 =3s^2du wedge ds.                                      (48)
```

Thus the finite bound is not merely a formal stopping argument: strict
stages can stack, their prefactor valuations telescope exactly, and only the
last stage need carry the Jacobian coefficient.  These are local rational
symplectic packets, not polynomial Keller maps.

## 8. Structural consequences and remaining obstruction

The proved recursion refines the open `C3` problem to a finite family of
positive compositions of `D`, decorated by the gcd chain `(28)`, the
constant primitive ratios `(12)`, and the terminal primitivity `(36)`.  For
fixed pole orders, the local differential anatomy is therefore finite-depth
even though the coefficients remain functional.

Since in both C3 geometries the initial complementary value `g_0=h` divides
`p` and the other displayed source order, one has

```text
D congruent to 3 modulo h.                              (49)
```

When `3` does not divide `h`, the partition sum already prevents any **fixed
nontrivial divisor of `h`** from dividing every depth.  Individual depths may
still meet different proper factors of `h`.  When `3|h`, an all-`3`-divisible
partition is arithmetically possible, but `(36)` forbids the tower from
remaining in that inertia-aligned lattice: some key depth must break the
surviving factor three.  This is a precise local sense in which the prime
three is exceptional in the `C3` budget.

This also explains the resemblance to the Euclidean leading-form tower in
[HYP-9070](../../05-knowledge/hypotheses/HYP-9070-jc2-leading-form-circuit-and-the-euclidean-depth-search-order.md):
both repeatedly replace a pair by primitive coprime exponents.  The present
result is not a bridge to that homogeneous-degree tower.  Here the proved
operation is gcd descent on local valuations plus subtraction from a
ramification budget; HYP-9070 studies subtraction in global degree pairs.

The highest-leverage remaining question is no longer whether a local tower
can be infinite.  It is whether simultaneous polynomiality of `P` and `Q`
can realize any of the finitely many decorated depth partitions, especially
after conjugating among the three `C2` views in an `A4` quartic.  No such
cross-chart physical intertwiner is supplied here.

## 9. Exact evidence and boundary

The companion, replayed byte-identically under normal and optimized Python
during an independent line audit, checks:

- every Bezout sign, logarithmic wedge, and gcd update for
  `1<=g_i,e_i<=24`;
- all `2^18-1` positive compositions through total budget eighteen;
- the primitive parametrization of every equal-weight terminal line on the
  same `24 x 24` grid;
- the exact one- and two-stage controls inherited from THM-3074;
- the new three-stage identities `(46)--(48)`; and
- the telescoping prefactors at every stage of that packet.

Reproduce with

```bash
python3 04-computation/jc_c3_finite_toric_key_tower_thm3080.py
python3 -O 04-computation/jc_c3_finite_toric_key_tower_thm3080.py
```

The theorem proves a completed-field differential recursion and its finite
budget.  It does not prove polynomiality of the Laurent keys, constrain an
arbitrary non-coordinate Jelonek component, couple conjugate quartic views,
or exclude `C3`, `A4/S4`, `G1`, `JC(2)`, or `DC(2)`.
