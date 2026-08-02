---
id: THM-3280
title: "FC cyclic-module simplicity and ordinary-basis repair through degree nine"
status: >
  PROVED + VERIFIED-EXACT.  For the rank-(d-1) cyclic connection attached to
  `At^d+Bt` with `A*B!=0`, every differential endomorphism over `Qbar(s)` is
  scalar.  For `3<=d<=9` the module is simple: a proper rank-m submodule would
  produce a primitive polynomial Pluecker vector of degree
  `N=sum(S)/d-m/2`, while `N<=(d-1)^2/(8d)<1`; degree zero would require a
  proper coordinate subset invariant under one nonzero cycle.  Simplicity and
  the scalar endomorphism ring classify the two marked triangle packets as
  independent copies or proportional copies modulo endpoint exponentials.
  In the proportional case the exact splitting theorem makes their difference
  a constant critical-exponential particular.  In either case the period
  scalar belongs to a linearly independent derivative-closed E-function vector
  with denominator only `s`.  Beukers Corollary 1.4 applies at `s=1`, proving
  that every non-pure affine binomial phase of degree `3<=d<=9` has simplex
  exponential period neither zero nor `1/2`.  In particular the non-pure cubic
  gap in THM-3252 is closed.  At `d=10`, the profile
  `m=4,S={6,7,8,9},N=1` only defeats this degree bound; it is not a submodule or
  counterexample.
source: codex-2026-08-02-fc3-cubic-frontier
depends_on:
  - THM-3253-fc-affine-binomial-cyclic-residue-exclusion
related:
  - THM-3250-fc3-noncollinear-pure-power-turn-current-exclusion
  - THM-3251-fc3-collinear-pure-power-spline-residue-exclusion
  - THM-3252-fc3-depressed-cubic-bessel-marked-copy-splitting-gate
external:
  - "F. Beukers, A refined version of the Siegel--Shidlovskii theorem, Annals of Mathematics 163 (2006), 369--379, Corollary 1.4."
script: 04-computation/fc_cyclic_module_ordinary_basis_repair_thm3280.py
output: 05-knowledge/results/fc_cyclic_module_ordinary_basis_repair_thm3280.out
script_sha256: 187f9befa3eb072f5499a488c554fd50439bea896b024214c4a31b1280f372d5
output_sha256: bd19bb7734ed9f901ec257e15dec64db5d3338d0bb64d4b79ee8a792bcc9d3da
---

# THM-3280 -- cyclic-module simplicity repairs specialization through degree nine

**PROVED + VERIFIED-EXACT.**

## 1. The result

Let

```text
Delta={(u,v):u>=0,v>=0,u+v<=1},
ell in Qbar[u,v] affine and nonconstant,
3<=d<=9,
A,B,C,h in Qbar,                         A*B!=0,
q=A(ell-h)^d+B(ell-h)+C,
K(s)=integral_Delta exp(sq)du dv.                         (1)
```

Then

```text
K(1)!=0,                              K(1)!=1/2.           (2)
```

The structural input from THM-3253 is valid for every `d>=3`; only the final
ordinary-point specialization was missing.  This theorem supplies that step
for `d<=9`.  At `d=3`, exact cubic depression plus THM-3250/3251 for the pure
branch therefore proves the two forbidden-value statements for **every**
algebraic affine-coordinate cubic.

Put `r=d-1`.  The homogeneous module at the heart of the proof is

```text
Y'=(Ccyc+D/s)Y,                                           (3)
D=diag(-1/d,-2/d,...,-r/d),

(Ccyc)_(j,j+1)=a=B(d-1)/d,              1<=j<r,
(Ccyc)_(r,1)=c=-B^2(d-1)/(d^2A).                         (4)
```

Both `a` and `c` are nonzero, so `Ccyc` is one weighted directed cycle.  Its
eigenvalues are distinct.  If `rho^r=-B/(dA)`, the corresponding right and
left eigenvectors may be chosen as

```text
q_rho=(1,rho,...,rho^(r-1))^T,
l_rho=(1,rho^-1,...,rho^(-(r-1)))^T,                     (5)
```

with `l_rho^Tq_sigma=r` for `rho=sigma` and zero otherwise.

The proof has two algebraic parts.  First, the differential endomorphism ring
of (3) is scalar in every degree.  Second, a Pluecker-degree estimate makes
the module simple through degree nine.  Section 6 then converts those facts
into the exact independent ordinary E-function vector missing from
MISTAKE-356.

## 2. Every rational differential endomorphism is scalar

Let `R in Mat_r(Qbar(s))` be a horizontal endomorphism of (3).  It obeys

```text
R'=[Ccyc,R]+[D,R]/s.                                     (6)
```

There is no pole at a finite nonzero point: differentiation would raise its
order while the right side would not.  If `R` had a pole of order `m>=1` at
zero with leading coefficient `R_(-m)`, the coefficient of `s^(-m-1)` would
give

```text
-m R_(-m)=[D,R_(-m)].                                    (7)
```

On matrix coordinate `(i,j)`, every eigenvalue of `ad(D)` is

```text
D_i-D_j=(j-i)/d,                     absolute value <1.   (8)
```

It cannot equal the negative integer `-m`; hence `R_(-m)=0`, a
contradiction.  Thus `R` is polynomial.

Suppose its degree is `N`, with top coefficient `R_N`.  The coefficient of
`s^N` in (6) gives

```text
[Ccyc,R_N]=0.                                             (9)
```

The coefficient one degree lower is

```text
N R_N=[Ccyc,R_(N-1)]+[D,R_N].                            (10)
```

Pair (10) by trace with every element `S` of the centralizer of `Ccyc`.
Because `S` and `R_N` commute,

```text
tr(S[Ccyc,R_(N-1)])=0,
tr(S[D,R_N])=0.                                           (11)
```

The trace pairing is nondegenerate on the centralizer of a regular semisimple
matrix.  Therefore `N>0` would force `R_N=0`.  Hence `N=0`.  Equation (6)
now says that `R` commutes separately with `Ccyc` and `D`.  The distinct
diagonal entries of `D` make `R` diagonal, and the nonzero cycle makes all its
diagonal entries equal.  Consequently

```text
End_diff(M)=Qbar                                           (12)
```

for every `d>=3`.

## 3. A proper submodule would have a polynomial Pluecker line

Assume that (3) has a differential submodule `W` of rank
`1<=m<r`.  Choose a rational basis matrix `P` for `W`.  There is an
`m`-square rational matrix `B_W` such that

```text
P'=(Ccyc+D/s)P-PB_W.                                     (13)
```

Let `w` be the wedge of the columns of `P`.  Clear denominators and divide
the polynomial gcd of its Pluecker coordinates.  The resulting nonzero
decomposable vector

```text
w in wedge^m(Qbar^r)[s]                                  (14)
```

is primitive: its coordinates have no common zero at any finite point.  For
the additive exterior actions `C^[m]` and `D^[m]`, it satisfies

```text
w'=(C^[m]+D^[m]/s)w-b(s)w                                (15)
```

for a rational scalar `b`.

At a finite nonzero point, some coordinate of `w` is nonzero and (15) shows
that `b` is regular.  At zero it has at most a simple pole.  If the polynomial
part of `b` had positive degree, its product with the top coefficient of `w`
would have degree larger than every other term of (15).  Therefore

```text
b(s)=beta+gamma/s.                                        (16)
```

Write

```text
w=w_0+w_1s+...+w_Ns^N,          w_0*w_N!=0.               (17)
```

The bottom and top equations of (15) are

```text
(D^[m]-gamma I)w_0=0,
(C^[m]-beta I)w_N=0.                                     (18)
```

Because the limits of a decomposable Pluecker vector are decomposable,
`w_0` and `w_N` represent invariant `m`-planes for `D` and `Ccyc`,
respectively.  The first is a coordinate plane

```text
W_0=span{e_j:j in S},                 |S|=m,              (19)
gamma=-sum_(j in S)j/d.
```

The second is spanned by `m` of the eigenvectors (5).

## 4. The exact degree identity

The coefficient of `s^(N-1)` in (15) is

```text
Nw_N=(C^[m]-beta I)w_(N-1)+(D^[m]-gamma I)w_N.            (20)
```

Let `T` be the `m`-element set of eigenlines spanning the top plane and pair
(20) with the corresponding left eigen-wedge.  The first term on the right
vanishes.  It remains to evaluate the normalized expectation of `D^[m]`.

For each eigenline in (5), the spectral rank-one projector has every diagonal
entry equal to `1/r`.  Hence the projector onto the `m` selected eigenlines
has diagonal `m/r`.  The additive exterior expectation is its trace against
`D`, so

```text
<D^[m]>= (m/r)tr(D)
         =(m/r)(-r(r+1)/(2d))
         =-m/2,                         because r+1=d.     (21)
```

Equations (19)--(21) give the exact integer identity

```text
N=sum_(j in S)j/d-m/2.                                    (22)
```

For an `m`-subset of `{1,...,d-1}`, the largest possible sum is

```text
m(2d-m-1)/2.
```

Therefore

```text
0<=N<=m(d-m-1)/(2d)<=(d-1)^2/(8d).                        (23)
```

For `3<=d<=9`, the last quantity is strictly less than one.  Since `N` is a
nonnegative integer, `N=0`.

Then `w` is constant and (18) says its plane is invariant under both `D` and
`Ccyc`.  A `D`-invariant plane is a coordinate subset.  No nonempty proper
coordinate subset is invariant under the single directed cycle (4).  This
contradiction proves

```text
M is simple over Qbar(s),                    3<=d<=9.       (24)
```

The numerical endpoint is honest.  At `d=10`,

```text
m=4,                  S={6,7,8,9},
sum(S)=30,            N=30/10-4/2=1.                       (25)
```

Thus (23) no longer forces a constant Pluecker vector.  Formula (25) is only
a hostile to this proof mechanism.  It neither solves the remaining
coefficient recurrences nor constructs a decomposable invariant submodule.
Simplicity for `d>=10` remains open here.

## 5. Copy dichotomy modulo endpoint exponentials

Let `mathcal E` be the `Qbar(s)`-span of the distinct endpoint exponentials,
`1`, and the extra zero-source exponential `exp(-Cs)` needed for the value
`1/2`.  It is derivative-stable.  THM-3253 proves that adjoining the extra
exponential creates no new rational split.

If a packet `U` satisfies the inhomogeneous connection of THM-3253, its class
`u` modulo `mathcal E` satisfies (3).  If `u!=0`, simplicity implies that its
coordinates are independent: the rational row relations form a differential
submodule of the dual of `M`, so a nonzero relation would kill every
coordinate.  Thus the coordinate span of `u` is a simple copy of the
coordinate/dual module.  Simplicity and the scalar endomorphism statement pass
from `M` to its dual.

For two nonzero packets `U,V`, their coordinate spans modulo `mathcal E` are
either disjoint or equal.  In the equal case, the coordinatewise map

```text
u_j |-> v_j                                                (26)
```

is a differential endomorphism of one copy.  Equation (12) makes it scalar:

```text
v=lambda u,                         lambda in Qbar.         (27)
```

Equivalently, `V-lambda U` splits over `mathcal E`.  The exact splitting
theorem in THM-3253 is stronger than mere rational splitting: every
noncritical coefficient vanishes, every critical coefficient is constant,
and every adjoined zero-source coefficient vanishes.  Hence

```text
V-lambda U = a constant Qbar-linear combination of
             critical endpoint exponentials.              (28)
```

This constant-particular conclusion is what prevents an apparent singularity
from reappearing at `s=1`.

## 6. Construction of the ordinary independent vector

THM-3253 has two geometric branches.

### 6.1 Three distinct endpoint marks

The triangle weights give packets `U,V` and scalar

```text
F=U_2-V_1,
Omega K(s)=exp(Cs)F(s),                  Omega!=0.          (29)
```

The collision audit proves `U notin mathcal E^r`, so `u!=0`.

If `v=0`, (28) with `lambda=0` writes `V` as a constant exponential
particular.  The vector consisting of the `r` coordinates of `U` and the
distinct exponential basis is linearly independent over `Qbar(s)`, and a
constant basis change includes `F`.

If `u` and `v` generate disjoint copies, all `2r` packet coordinates together
with the exponential basis are independent.  The nonzero constant row
`U_2-V_1` extends to a constant invertible basis change containing `F`.

If the copies agree, (27)--(28) give

```text
F=U_2-lambda U_1-(constant exponential particular)_1.      (30)
```

The row `(-lambda,1,0,...,0)` on the coordinates of `U` is nonzero.  Thus a
constant invertible change of the independent vector `(U,exponentials)` again
contains `F`.

In all three cases the chosen vector is derivative-closed: the packet equation
has coefficients `Ccyc+D/s` and constant source rows divided by `s`, while
each exponential satisfies a constant first-order equation.  Its common
denominator is exactly a power of `s`, so `s=1` is ordinary.

### 6.2 One doubled endpoint mark

Here THM-3253 gives the nonsplit packet `T=Y_a-Y_b` and

```text
G=-bT_1+T_2,
(b-a)^2K(s)=exp(Cs)G(s).                                  (31)
```

Simplicity makes the coordinates of `T` independent modulo `mathcal E`.
The coefficient row `(-b,1,0,...,0)` is nonzero, so the same constant basis
change gives an independent derivative-closed vector containing `G`, again
ordinary at one.

## 7. Beukers at one and the two values

Every primitive packet and algebraic exponential above is an E-function by
THM-3253.  Section 6 supplies exactly the hypothesis missing in MISTAKE-356:
the **independent functions themselves** form a homogeneous first-order
rational system ordinary at `s=1`.  Beukers Corollary 1.4 therefore makes
their values at one linearly independent over `Qbar`.

If `K(1)=0`, (29) or (31) gives `F(1)=0` or `G(1)=0`.  If `K(1)=1/2`, it gives

```text
F(1)-(Omega/2)exp(-C)=0
```

or

```text
G(1)-((b-a)^2/2)exp(-C)=0.                                (32)
```

Each is a forbidden algebraic linear relation in the constructed ordinary
vector.  This proves (2).

## 8. Scope and the new frontier

The theorem proves the non-pure affine binomial branch only through degree
nine.  It does not prove simplicity in degree ten or higher, and (25) is not
evidence of reducibility.  Three next attacks are now sharply separated:

1. solve the degree-one Pluecker recurrence at the `d=10,m=4` hostile;
2. prove irreducibility of the pulled-back generalized-Bessel/Kloosterman
   operator by an exact differential-Galois theorem; or
3. bypass simplicity by proving local saturation at `s=1` only for the two
   triangle packets.

The proof also does not cover an intermediate monomial, a genuinely
multivariate phase, the projective leading-form bridge, or `FC(3)`.

## 9. Reproduction contract

Run

```bash
python3 04-computation/fc_cyclic_module_ordinary_basis_repair_thm3280.py
python3 -O 04-computation/fc_cyclic_module_ordinary_basis_repair_thm3280.py
```

The exact verifier checks the scalar simultaneous commutant and centralizer
trace pairing for `3<=d<=10`, enumerates every bottom residue profile through
degree nine, checks the sharp Pluecker bound, verifies that no proper
coordinate subset is cycle-invariant, and records both degree-one arithmetic
hostiles at `d=10`.  These are regression controls for
the symbolic proof, not a finite replacement for it.
