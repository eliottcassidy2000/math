# Consecutive Keller fibres, curl allocation, and the plaquette factor gate

**Status: PROVED + INDEPENDENTLY HOSTILE-AUDITED elementary UFD/curl and
coefficient-span gates; VERIFIED-EXACT controls; OPEN counterexample
construction.**  No counterexample to `JC(2)` and no proof of `JC(2)` is
claimed.  The proved outcome is a substantially narrower passport for a
hypothetical counterexample in the derivative-table representation.  In
particular, the two consecutive matching fibres must contain at least four
distinct irreducible factors, equality must be balanced `2+2`, a one-toric-
ray carrier is impossible in every degree, coefficient span at most two is
tame, and the rank-three constant-trace gauge has the stronger distinct-factor
floor five.

The exact companion is
[`keller_consecutive_factor_integrability_scout.py`](../04-computation/keller_consecutive_factor_integrability_scout.py),
with stored
[`output`](../05-knowledge/results/keller_consecutive_factor_integrability_scout.out).
It runs identically with `python3` and `python3 -O`.

## 1. Inheritance pass and concept board

The closest proved mechanism is the conductance-shadow gate in
`THM-3548-planar-keller-conductance-shadow-gates.md`: the four derivative
entries form the same `K_2,2` plaquette whose squared magnitudes survive
dephasing, while determinant information lives in the closing interference
phase.  The canonical hostile object for the present lane is the polynomial
`SL_2` matrix

```text
z=xy,                  [ z+1  z+2 ]
                       [  z   z+1 ].                    (1)
```

Its determinant is one and its two matching products are consecutive, but
both row-curl defects are `x-y`; it is not a Jacobian matrix.  Thus UFD
factorization without integrability is genuinely too weak.

The corrected near miss is `THM-3551-one-ray-planar-jacobian-mate-no-go.md`.
That theorem constrains a one-ray *output polynomial*.  The present theorem
instead constrains a one-ray *matching product* of four derivatives.  Neither
statement contains the other.  The least-used decisive sidecar is
`THM-2063-one-fiber-linear-planar-keller-pairs.md`: every time a curl or a
rank-one observer produces a constant directional derivative, that theorem
closes the cell tamely.  The rank-three boundary routes to
`THM-3367-berggren-spinor-pencil-hessian-gauge-and-affine-line-keller-closure.md`
and the still-open binary constant-Hessian program
`HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs.md`.

The live board at the end of the session is:

| lane | object | representation | surviving invariant | status |
|---|---|---|---|---|
| Anchor | planar Keller derivative | `[[a,b],[c,d]]` | `ad-bc=kappa` plus both curls | OPEN globally |
| Niche | matching pencil | `T=bc`, `T+kappa=ad` | two reducible consecutive fibres | PROVED factor gates |
| Niche | coefficient image | nonconstant matrix span `S_F` | rank-one observer or trace gauge | PROVED router |
| Wildcard | toric carrier | `T=h(x^p y^q)` | exponent transport through curls | PROVED impossible |
| Wildcard | scalar carrier | `T=h(R)` | irreducibility of selected `R`-levels | PROVED conditional gate |
| Generator | Gauss chart | `E_+(V)E_-(U)E_+(W)` | balanced factors plus two curl PDEs | OPEN construction cell |
| Hostile control | three-shear map | balanced fibres, gradient gcd one, span four | explicit polynomial inverse | VERIFIED-EXACT tame |

## 2. The derivative plaquette and what its shadow forgets

Let

```text
F=(P,Q) in C[x,y]^2,       Jac(P,Q)=kappa in C*,
a=P_x,  b=P_y,             c=Q_x,  d=Q_y.               (2)
```

The determinant and the two integrability equations are

```text
ad-bc=kappa,               a_y=b_x,       c_y=d_x.      (3)
```

Put

```text
T=bc.                                                   (4)
```

Then the two perfect matchings of the plaquette are the consecutive
polynomials

```text
bc=T,                     ad=T+kappa.                   (5)
```

The map from `(a,b,c,d)` to `(T,T+kappa)` preserves the determinant and the
two matching-factor partitions.  It destroys which factors were allocated to
which edge and, more seriously, whether either row is exact.  The missing
sidecar is precisely the pair of curls in `(3)`.  Matrix `(1)` is the cheapest
decisive hostile for this information loss.

There is a literal link to the dephasing picture.  The resistor shadow keeps
the four edge magnitudes, and the first `K_2,2` Wilson correction compares the
two matching channels.  Algebraically, `(5)` keeps those channels before
taking magnitudes.  Integrability is extra global information not present in
either a conductance table or an arbitrary determinant-one plaquette.

## 3. Elementary unit and carrier lemmas

All factor counts below are over `C[x,y]`.  Write `Omega(f)` for the number of
irreducible factors counted with multiplicity and `omega(f)` for the number of
distinct nonassociate irreducible factors.

### 3.1 Cross ideals and constant edges

Equation `(3)` immediately gives

```text
(a,b)=(a,c)=(d,b)=(d,c)=C[x,y].                         (6)
```

Indeed `kappa=ad-bc` belongs to each displayed ideal.  If any of `a,b,c,d` is
constant, `F` is tame.  For example, if `a=alpha!=0`, then

```text
P=alpha x+p(y).
```

In source coordinates `(U,V)=(P,y)`, the determinant equation says

```text
partial_V Q=kappa/alpha,
```

so `Q=(kappa/alpha)V+H(U)`.  If `a=0`, the determinant makes `b` and `c`
units and gives the other triangular form.  This is the smallest special
case of THM-2063.

Consequently, in a hypothetical nonautomorphism all four edges are nonzero
nonunits.

### 3.2 Irreducible-carrier lemma

Let `h in C[x,y]` be irreducible and let `A,B in C[t]` be nonconstant.  If

```text
A(h)_y=B(h)_x,                                           (7)
```

then `h` is affine linear and a nonzero linear combination of `A(h)` and
`B(h)` is constant.

To prove this, write

```text
A'=g A_1,                B'=g B_1,       gcd(A_1,B_1)=1.
```

The curl becomes

```text
A_1(h) h_y=B_1(h) h_x.                                  (8)
```

Bezout in `C[t]` remains Bezout after substituting `h`, so
`gcd(A_1(h),B_1(h))=1`.  Hence `A_1(h)` divides `h_x` and `B_1(h)` divides
`h_y`.  Neither partial can vanish under the nonconstant hypotheses.  A
nonconstant `A_1(h)` has degree at least `deg h`, while `deg h_x<deg h`; the
same applies to `B_1`.  Therefore `A_1,B_1` are constants.  Equation `(8)`
makes `h_x,h_y` proportional, so `h` is a polynomial in one linear form.
Irreducibility over `C` makes it degree one in that form.  Finally `A'` and
`B'` are proportional, so `B-lambda A` is constant.  THM-2063 closes the
resulting Keller map.

This degree-divisibility proof replaces a more geometric generic-fibre
argument and has no hidden smoothness assumption.

## 4. The consecutive-factor theorem

### Theorem

For `(2)--(5)`, if `F` is not a polynomial automorphism, then

```text
Omega(T)>=2,       Omega(T+kappa)>=2,                   (9)
Omega(T)+Omega(T+kappa)>=4,                            (10)
omega(T)+omega(T+kappa)>=4.                            (11)
```

Moreover, equality in `(11)` is necessarily balanced:

```text
omega(T)=omega(T+kappa)=2.                             (12)
```

The multiplicity and distinct-factor assertions are different.  Equation
`(10)` is immediate because `b,c,a,d` are nonunits.  Equations `(11)--(12)`
are the new UFD-plus-curl content.

### 4.1 Excluding total distinct count at most three

Suppose one matching side has only one distinct prime; after swapping source
coordinates if necessary, write

```text
T+kappa=u h^n,                                          (13)
```

where `u` is a unit and `h` is irreducible.  Since both `a,d` are nonunits,
`n>=2`.  The roots `r_1,...,r_n` of `u z^n=kappa` are distinct and nonzero,
so

```text
T=u product_j (h-r_j).                                  (14)
```

The factors in different levels are coprime, hence `omega(T)>=n`.  If the
total distinct count were at most three, then `n=2` and each `h-r_j` would
have only one distinct prime.  Such a level cannot be a proper power: if

```text
h-r_j=v q^m,          m>1,
```

then `h=vq^m+r_j` splits over `C` into `m` nonconstant shifts of `q`,
contrary to the irreducibility of `h`.  Thus both levels are irreducible.
Every divisor of `(13)` or `(14)` now lies in `C[h]`; in particular
`a=A(h),b=B(h)` with `A,B` nonconstant.  The irreducible-carrier lemma and
THM-2063 make `F` tame.  This proves `(11)`.

### 4.2 The multiplicity-safe unbalanced equality cell

It remains to exclude distributions `(1,3)` and `(3,1)` when the total in
`(11)` is four.  Continue with `(13)`.

If `n=3`, the three level polynomials in `(14)` each have one distinct factor
and the same proper-power argument makes all three irreducible.  Again all
edges lie in `C[h]`, and the carrier lemma closes the map.

The only remaining case is `n=2`.  One of the two opposite levels splits
into two distinct primes, with arbitrary multiplicities, and the other level
is irreducible.  After a constant target gauge and rescaling, the derivative
identities take the exact form

```text
BC=R(R+2),                B_x=R_y,       C_y=R_x,       (15)
```

while the other two edges are scalar multiples of `R+1`.  Choose the sign of
`R` so that `R+2` is the irreducible level.

If `R+2` is assigned to `C`, write

```text
B=S,                  C=E(R+2),          R=SE.          (16)
```

No squarefreeness is assumed: `S,E` are complementary divisors carrying all
multiplicities.  Expanding the second curl and using the first gives

```text
(R+2)E_y=S E_x.                                           (17)
```

Since `gcd(R+2,S)=1`, equation `(17)` makes `R+2` divide `E_x`.  If `S` is
nonconstant,

```text
deg(R+2)=deg S+deg E>deg E_x,
```

so `E_x=E_y=0`: this is the whole-level allocation `B~R,C~R+2`.  If `S` is
constant, `B` is already a constant edge.  When `R+2` is assigned to `B`,
write

```text
B=S(R+2),             C=E,               R=SE.          (18)
```

The symmetric expansion is

```text
(R+2)S_x=E S_y,                                           (19)
```

and gives the same dichotomy.  A whole-level allocation makes `R_x,R_y`
proportional; equivalently, a constant linear combination of an edge
`R` and an edge `R+1` appears.  THM-2063 closes it.

This proves `(12)` with arbitrary factor multiplicities.  The squarefree
six-allocation enumeration is only a finite shadow of `(16)--(19)` and is
not used in the theorem.

### 4.3 Sharp tame control

Put

```text
v=x+y,
Q=y+v^2/2,                 P=v+Q.                        (20)
```

Then

```text
DF=[ v+1  v+2 ],          det DF=1,
   [  v   v+1 ]

T=v(v+2),                 T+1=(v+1)^2.                  (21)
```

Thus `Omega(T)+Omega(T+1)=4` and
`omega(T)+omega(T+1)=3`.  The inverse is

```text
v=P-Q,       y=Q-v^2/2,       x=v-y.                    (22)
```

So the multiplicity floor is numerically sharp among Keller maps, while the
distinct-factor strengthening necessarily uses the hypothetical
nonautomorphism assumption.

## 5. Scalar carriers and the all-degree toric closure

### 5.1 Selected-level scalar carrier

Suppose

```text
T=h(R)                                                   (23)
```

and every `R-r` is irreducible for every distinct root `r` of
`h(t)(h(t)+kappa)`.  Then `F` is tame.

If `h` is constant, then `bc` is constant.  A nonzero constant makes both
factors units, while zero makes one factor the zero constant; either case is
already closed by Section 3.1.  Thus assume `h` is nonconstant, so a selected
root exists.

UFD factorization places every one of `a,b,c,d` in `C[R]`.  Choose any
selected root `r_0`; then `R-r_0` is irreducible and generates the same
univariate ring.  The irreducible-carrier lemma applied to the first curl
gives a constant directional derivative, and THM-2063 closes the map.

Therefore a scalar-carrier counterexample generator needs at least one
reducible distinguished fibre of `R`.  This is a necessary gate, not a
classification when a selected fibre is reducible.

### 5.2 Matching toric-ray theorem

Let `p,q` be coprime positive integers and `z=x^p y^q`.  If either

```text
T=h(z)              or              T+kappa=h(z)        (24)
```

for `h in C[t]`, then `F` is tame.  This closes exactly the first reducible-
carrier exception to the scalar theorem in all degrees.

It is enough to treat `T=h(z)`.  If `h` is constant, `bc` is a unit and a
constant edge appears when `h!=0`; when `h=0`, the domain property makes
`b=0` or `c=0`, again a constant edge.  Suppose it is nonconstant.

If `h(0)` is neither `0` nor `-kappa`, every factor of `h(z)` and
`h(z)+kappa` is a binomial `z-rho` with `rho!=0`.  Primitivity of `(p,q)`
makes each binomial irreducible, so all four edges lie in `C[z]`.  The first
curl has the form

```text
q x A'(z)=p y B'(z).                                    (25)
```

The left monomials have exponents `(pr+1,qr)` and the right monomials have
exponents `(ps,qs+1)`; no pair can agree.  Hence both derivatives vanish and
a constant edge appears.

Now suppose `h(0)=0`.  The nonzero-level factors still put

```text
a=A(z),       d=D(z),
b=x^r y^s B(z),          c=x^u y^v C(z),                (26)
```

with `B(0)C(0)!=0`.  The two curls force positive integers `delta,epsilon`
such that

```text
r=p delta+1,       s=q delta-1,
u=p epsilon-1,     v=q epsilon+1,                        (27)

q A'=z^(delta-1)[(p delta+1)B+pzB'],
p D'=z^(epsilon-1)[(q epsilon+1)C+qzC'].                (28)
```

Let

```text
K=delta+deg B,                 L=epsilon+deg C.          (29)
```

The top coefficient ratios are

```text
lc(A)/lc(B)=(pK+1)/(qK),
lc(D)/lc(C)=(qL+1)/(pL).                                 (30)
```

But cancellation of the top term in `ad-bc` would require their product to
be one, whereas

```text
(pK+1)(qL+1)-pqKL=pK+qL+1>0.                            (31)
```

This contradiction also covers the exponent boundary cases in `(27)`.
Finally, if `h(0)=-kappa`, swap the two source coordinates.  The new matching
product is the old `T+kappa`, the new determinant is `-kappa`, and the
zero-level argument applies.  The second alternative in `(24)` is the same
source swap at the outset.

The obstruction is therefore not merely that toric factors look sparse.  The
curls shift the two allocations in opposite lattice directions, and their
top ratios miss determinant closure by the strictly positive gap `(31)`.

## 6. Coefficient-span observers

Fix a base point and let

```text
S_F=span_C {DF(z)-DF(z_0): z in C^2}.                   (32)
```

Equivalently, `S_F` is the span of the nonconstant coefficient matrices of
`DF`.  Use the Frobenius pairing

```text
<W,M>=tr(W^T M).                                        (33)
```

If `W=u v^T` has rank one and lies in `S_F^perp`, then

```text
u^T DF v                                                (34)
```

is constant.  Complete `u^T` to an invertible target change and `v` to an
invertible source basis.  Equation `(34)` becomes a constant derivative
entry, so THM-2063 makes `F` tame.

Every complex matrix subspace of dimension at least two contains a nonzero
singular matrix: restrict the quadratic determinant to a projective line and
take a root.  Therefore

```text
dim S_F<=2        implies        F is tame.              (35)
```

For `dim S_F=3`, the annihilator is a line.  A nonautomorphic candidate must
have an invertible generator `W`; a singular generator would trigger `(34)`.
Apply the source change

```text
z -> W^T z.                                             (36)
```

Cyclicity of trace and `(33)` make the transformed Jacobian have constant
trace `tau`.  Its determinant is the transformed constant

```text
kappa_tr=(det W)kappa.
```

For Sections 7--9, rename this nonzero transformed constant `kappa`; every
matching product and factor count there is in this fixed trace gauge.  With
`lambda=tau/2`, polynomial exactness gives a potential
`H` such that

```text
F=(lambda x+H_y,       lambda y-H_x),
det DF=lambda^2+det Hess(H).                            (37)
```

This is a router, not a solution of the binary Hessian problem.  In
particular, no equivalence between invertibility of `F` and invertibility of
`grad H` is asserted.  THM-3367 closes the affine-line coefficient image and
the top homogeneous one-ray statement; HYP-8905 records the remaining
constant-nonzero binary Hessian lane.

## 7. Rank three has a five-factor floor in its trace gauge

There is nevertheless one further unconditional squeeze.  In the trace
gauge `(37)`, put

```text
R=H_xy,                 mu=det Hess(H)=kappa-lambda^2.
```

Then

```text
a=lambda+R,       b=H_yy,       c=-H_xx,       d=lambda-R,
H_xx H_yy=R^2+mu,                                             (38)

T=bc=-(R^2+mu),             T+kappa=lambda^2-R^2.             (39)
```

Thus the open Hessian lane is exactly a scalar matching carrier whose
distinguished `R`-levels are the roots of `z^2+mu` and `z^2-lambda^2`.

### Theorem

If a rank-three coefficient-span Keller map is nonautomorphic, then in its
constant-trace gauge

```text
omega(T)+omega(T+kappa)>=5.                              (40)
```

The phrase "in its constant-trace gauge" is load-bearing: the individual
matching product `P_yQ_x` is not invariant under an arbitrary source change.

By the global theorem, an equality-four survivor would be balanced `2+2`.
There are three cases.

1. If `lambda mu!=0`, the four selected values of `R` in `(39)` are distinct;
   a coincidence would give `kappa=0`.  Equality four makes every selected
   level have one distinct prime.  A level `R-r=u f^m` with `m>1` makes every
   other level `u f^m+constant` split into `m` distinct shifts of `f`, a
   contradiction.  All four levels are therefore irreducible, and the scalar-
   carrier theorem closes the map.
2. If `lambda=0`, then `mu=kappa!=0`.  Equality four says `omega(R)=2` and
   each of `R-alpha,R+alpha`, `alpha^2=-mu`, has one distinct prime.  The same
   multiplicity point needs care.  Write

```text
R-alpha=u f^m,                 R+alpha=v g^n
```

   with irreducible `f,g`.  If `m,n>=2`, choose a constant derivation `D`
   with `Df,Dg!=0` and differentiate their nonzero constant difference.
   Coprimality gives

```text
(m-1)deg f<=deg g-1,          (n-1)deg g<=deg f-1,
```

   an immediate contradiction.  If, say, `m=1<n`, then
   `f` is a nonzero scalar multiple of `g^n-constant`, which splits into
   `n` nonconstant shifts of `g` over `C`; the other asymmetric case is the
   same.  Hence `m=n=1`, so both levels really are irreducible.  Since `b,c`
   are nonunits, each takes one whole opposite level.  The curl
   `b_x=a_y` makes `R_x,R_y` proportional and produces a constant derivative;
   THM-2063 closes the map.
3. If `mu=0`, then `lambda!=0`, `bc=-R^2`, and equality four makes `R` have
   exactly two distinct primes.  Set

```text
B=b=H_yy,                 C=-c=H_xx.                    (41)
```

   Then

```text
BC=R^2,                  B_x=R_y,       C_y=R_x.        (42)
```

   Write

```text
R=rho f^m g^n,
B=beta R L,              C=beta^(-1)R/L,
L=f^r g^s,               -m<=r<=m,  -n<=s<=n.          (43)
```

   The two curls in `(42)` express `R_y` in two ways and give the rational
   Burgers equation

```text
L_y=beta L L_x.                                           (44)
```

   If `L` is constant, `(42)` makes `R_x,R_y` proportional and THM-2063
   closes.  If a prime `h` among `f,g` is mixed (`h_x h_y!=0`) and has
   nonzero exponent `e` in `L`, then

```text
ord_h(L_y)=e-1,             ord_h(L L_x)=2e-1,          (45)
```

   because an irreducible polynomial in characteristic zero cannot divide
   either nonzero partial.  Equation `(44)` would force `e=0`, a
   contradiction.  Hence every mixed prime has exponent zero.  One mixed and
   one pure prime leaves a nonzero power of a one-variable polynomial in
   `L`, impossible in `(44)`.  Two pure primes in the same variable make
   `(44)` force `L` constant.

   The only remaining case has `f=f(x)` and `g=g(y)`.  Irreducibility over
   `C` makes both affine linear.  Exponent comparison in `(44)` forces, up to
   swapping, `(r,s)=(1,-1)`.  The original two curls then give

```text
beta(m+1)f_x=n g_y,          beta m f_x=(n+1)g_y,        (46)
```

   which would imply `mn=(m+1)(n+1)`, impossible for positive `m,n`.

All equality-four cases are tame, proving `(40)`.  This proof is specific to
the equality-four factor allocation; it does not assert the general
nonhomogeneous Hesse-zero classification that is absent from THM-3367.

## 8. Exact balanced controls and the meaning of gradient gcd one

The factor floors are necessary, not close to sufficient.  Two explicit tame
maps occupy the first surviving balanced cell.

### 8.1 Two shears: balanced and rank-three observable

Put

```text
u=x+y^2/2,             v=y+u^2/2,
Q=v,                   P=u+v.                           (47)
```

Then

```text
a=1+u,                 b=1+y+uy,
c=u,                   d=1+uy,                          (48)

T=u(1+y+uy),           T+1=(1+u)(1+uy).                (49)
```

All four displayed factors are distinct and irreducible.  Also
`gcd(T_x,T_y)=1`, so `T` cannot be a nontrivial polynomial composite
`H(R)`; otherwise `H'(R)` would divide both partials.  However, this is only
**gradient-coprime**, not no-critical-points: the exact gradient ideal has
Groebner basis

```text
[x+y/2+1, y^2-1].                                      (50)
```

The coefficient span has dimension three and Frobenius annihilator

```text
W=[-1 0; 1 0]=u_0 v_0^T,              c-a=-1.          (51)
```

Thus the rank-one observer exposes the constant derivative.  The inverse is

```text
u=P-Q,       y=Q-u^2/2,       x=u-y^2/2.                (52)
```

### 8.2 Three shears: balanced, gradient-coprime, and full span

Keep `u,v` from `(47)` and put

```text
Q=v,                   P=u+v^2/2.                       (53)
```

Then

```text
a=1+uv,                b=y+v(1+uy),
c=u,                   d=1+uy,                          (54)

T=u[y+v(1+uy)],        T+1=(1+uv)(1+uy).               (55)
```

Again the two fibres have factor pattern `(2,2)`, and independent SymPy and
python-FLINT computations give `gcd(T_x,T_y)=1`.  Yet the critical ideal is
`(x,y)`: gradient-coprime must not be confused with the unit gradient ideal.

This time the coefficient span is all of `M_2(C)`.  The coefficient columns
at monomials `y,y^2,y^3,y^6`, in entry order `(a,b,c,d)`, are

```text
(0,2,0,0),
(0,0,1/2,0),
(1/2,0,0,1/2),
(1/16,0,0,0),                                      (56)
```

and have determinant `-1/32`.  The explicit inverse is

```text
v=Q,       u=P-Q^2/2,       y=Q-u^2/2,       x=u-y^2/2. (57)
```

Therefore balanced reducible fibres, noncomposite `T`, and full coefficient
span four can all coexist tamely.  They are passport fields, not a
counterexample certificate.

## 9. A sharply parameterized balanced generator

The three-shear example is one point in a useful three-element Gauss chart;
this chart is not asserted to cover every polynomial `SL_2` table.  For
polynomials `U,V,W`, let

```text
E_+(s)=[1 s;0 1],             E_-(s)=[1 0;s 1].
```

Then

```text
M=E_+(V)E_-(U)E_+(W)
 =[1+UV,  V+W(1+UV);  U, 1+UW],                         (58)
det M=1,

T=U[V+W(1+UV)],             T+1=(1+UV)(1+UW).          (59)
```

The determinant and balanced factorization are automatic.  The whole global
problem in this chart is the pair of curl PDEs

```text
U_y=(UW)_x,
(1+UV)_y=[V+W(1+UV)]_x.                                 (60)
```

Equivalently, for potentials with `Q_x=U,Q_y=1+UW`,

```text
(partial_y-W partial_x)Q=1,
(partial_y-W partial_x)P=V.                             (61)
```

The three-shear control is the polynomial solution

```text
U=x+y^2/2,             V=y+U^2/2,        W=y.           (62)
```

This suggests a bounded counterexample generator: solve `(60)` outside the
obvious sequential-shear orbit, insist that the four factors in `(59)` are
distinct, retain coefficient span four, and then impose a genuine collision
and the known all-line/nonproperness gates.  The ansatz is not sufficient and
may itself ultimately be tame; its value is that it leaves only two explicit
polynomial PDEs after paying determinant and the minimal factor passport.

There is also an exact positive boundary bridge to the nodal construction.
For the normalized boundary frame

```text
M_k(u)=
 [2u,       1+2ku;
  3u^2-1,   3u/2+k(3u^2-1)],                            (63)
```

one has

```text
M_k=E_-(3u/2-1) E_+(1) E_-(2u-1) E_+(k).               (64)
```

The branch holonomy is the final right shear
`E_+(k_+-k_-)`.  This is only a boundary-frame factorization, not a global
integrable lift, but it shows that the balanced elementary language and the
nodal collision language are compatible rather than disjoint.

## 10. Failure frontier and next decisive tests

After these gates, a derivative-table counterexample must avoid all of the
following exits simultaneously:

```text
two matching fibres reducible with total omega>=4;
if total omega=4, the allocation is exactly balanced (2,2);
no constant edge or constant source/target directional derivative;
not a scalar carrier with every selected level irreducible;
not a positive toric matching ray in any degree;
coefficient span at least three;
if span three, the unique annihilator is invertible and the trace-gauge
  matching fibres have total omega at least five;
if span four, no linear coefficient observer remains;
plus the independent nonproperness, all-line, and collision requirements.
```

The smallest genuinely hostile construction cell is therefore a balanced
`(2,2)` table with both curls, no rank-one observer, and preferably full span
four.  The three-shear map proves that this cell is populated but tame.  The
next cheap tests for `(58)--(60)` are:

1. quotient constant affine gauges and the evident sequential-shear family;
2. solve bounded-support instances of the two curl PDEs exactly;
3. test whether every solution acquires a low-fibre-degree output-pencil
   member, invoking THM-2063;
4. impose the nodal boundary frame `(63)` and see whether a fourth elementary
   layer or a variable pivot is forced globally;
5. for any survivor, compute coefficient span, the two factor counts, the
   gradient ideal (not just its gcd), and an explicit collision before raising
   the degree cap.

The conceptual squeeze is productive even though it points toward a proof:
the resistor/dephasing analogy selects the two matching channels, UFD exposes
their consecutive fibres, and curls convert factor allocation into rigid
degree transport.  What remains is no longer an arbitrary determinant-one
matrix and not yet the whole Jacobian conjecture.  It is a small, explicit
integrability problem at the balanced plaquette boundary.

## 11. Reproduction and exact scope

Run

```bash
python3 04-computation/keller_consecutive_factor_integrability_scout.py
python3 -O 04-computation/keller_consecutive_factor_integrability_scout.py
```

The companion checks the sharp factor-floor automorphism, the nonintegrable
toric `SL_2` hostile, the symbolic leading gap and the full integer box
`1<=p,q,K,L<=7`, both multiplicity-safe allocation identities, the universal
Gauss chart, the nodal boundary factorization, both tame balanced controls and
their inverses, the rank-three/rank-four coefficient witnesses, two exact
gradient-gcd paths, and the Hamiltonian matching identities `(38)--(39)`.

The universal all-degree conclusions are proofs in this reflection, not
finite extrapolations from the script.  The script does not certify
irreducibility in arbitrary degree, global nonproperness, injectivity, or the
existence of a counterexample.  The rank-three theorem is a factor floor in
the constant-trace gauge, not a closure of HYP-8905.  The general balanced
`(2,2)`, span-four integrability cell remains open.
