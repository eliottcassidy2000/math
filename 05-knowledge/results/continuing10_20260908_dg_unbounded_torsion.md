# Graph-complement submersions with unbounded unit torsion and an actual affine-chart switch

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is an audited corollary/generalization of the named current suppliers.
It claims no literature priority. JC(2) remains OPEN.
Unboundedness below varies the surface index m; no unboundedness on a fixed
W2 is asserted. A rational mate is not a polynomial mate on the original
source, or a pair of global functions on the entire surface.

## 1. Inheritance, definitions, and complete statement

Work over C. The actual charts are inherited from
`05-knowledge/results/planar_jc48_sep08_dg_genus.md`:

```
W_m=(P1_x x P1_z) minus {z=x^m},              m>=1,
U0=A2_(x,t),       Uinf=A2_(r,b),
x=1/r,   t=-r^m-r^(2m)b,   b=-x^m-x^(2m)t,
D=W_m minus U0={r=0},
omega=dx wedge dt=r^(2m-2) dr wedge db.
```

The m2 antecedent is
`05-knowledge/results/continuing10_20260907_dg_linear_carrier.md`.
The torsion and canonical derivative supplier is
`05-knowledge/results/planar_jc48_sep06_torsion.md`, which already contains
unbounded order r-1 controls `P=x+lambda*x^r*y` on the affine plane.
The present contribution is their explicit realization on actual W_m
surfaces by globally smooth first functions, with an exact chart switch
and a matching boundary ramification degree. The incoming actual
boundary-linear order-two family is
`05-knowledge/results/planar_jc48_sep08_unit_torsion.md`.

Fix `a,h in C*`, `B0 in C`, and define

```
e=2m-1,
Q_m(x)=sum_(j=0)^m binom(2m,j)(-h)^j x^(m-j),
A=a(x-h)^(2m),             B=a Q_m+B0,
F=A t+B,                  f0=B(h),               g=F-f0,
z=x-h,                    K=z^e t+[Q_m(h+z)-Q_m(h)]/z,
G0=1/[a e (x-h)^e].
```

Q_m is the entire polynomial part of `(x-h)^(2m)/x^m`, including its
constant coefficient. Thus

```
f0=B0+a(-1)^m binom(2m-1,m)h^m,
kappa=Q_m'(h)=(-1)^(m-1)binom(2m-2,m-1)h^(m-1) !=0.
```

**Theorem.** For every stated m,a,h,B0:

1. F is a global regular function on W_m, has no critical point anywhere
   on W_m, and has boundary restriction `F|D=B0-a b`.
2. Every fibre other than f0 is A1 on W_m. The special fibre is the
   disjoint union of two reduced A1 components `E1={x=h}` and the closure
   E2 of `K=0`. Their restrictions to U0 are respectively A1 and G_m.
   Generic fibres on U0 are G_m, not A1.
3. `u=1/(x-h)` and F identify the actual open surface
   `W_m minus E1` with `A2_(F,u)`. In this chart
   `omega=a^-1 u^(e-1) dF wedge du` and `G0=u^e/(a e)`.
   Consequently `(F,G0)` is a finite flat map of degree e on this
   alternative affine chart. For m>=2 its ramification locus is exactly
   the old boundary D={u=0}, with index e. For m1 it is an unramified
   degree-one coordinate rescaling.
4. In the ORIGINAL source module
   `C_F=C[x,t]/D_F(C[x,t])`, `D_F=F_x partial_t-F_t partial_x`,
   the unit class theta has the complete annihilator
   `Ann_(C[F])(theta)=((F-f0)^e)`.
   Its j-th canonical transverse-connection derivative has exact
   annihilator `((F-f0)^(e+j))` for every j>=0.
5. In particular F has a rational mate G0 but no polynomial mate in
   C[x,t], of any degree, and hence no mate globally regular on W_m.

The response module in (4) is not a quotient of O(W_m). For m>=2 the
Hamiltonian derivation obtained by dividing by omega need not preserve
O(W_m), because omega vanishes on D. Global smoothness of F is a
geometric sidecar, not permission to replace the coefficient ring.

## 2. Globality, constants, and submersion

Substitute the actual second chart and cancel the polynomial part:

```
F_inf=B0-a b(1-hr)^(2m)
       -a sum_(j=m+1)^(2m) binom(2m,j)(-h)^j r^(j-m).       (1)
```

This is polynomial everywhere on Uinf, including its entire boundary.
It proves the boundary formula and `partial_b F_inf|D=-a !=0`.
On U0, `F_t=a(x-h)^(2m)` can vanish only at x=h, where
`F_x=B'(h)=a kappa !=0`. These two tests cover all W_m.

For clarity, the two binomial evaluations used above follow from
`sum_(j=0)^k (-1)^j binom(n,j)=(-1)^k binom(n-1,k)`.
The derivative sum is

```
sum_(j=0)^m (m-j)(-1)^j binom(2m,j)
 =m(-1)^m binom(2m-1,m)-2m(-1)^m binom(2m-2,m-1)
 =(-1)^(m-1)binom(2m-2,m-1).
```

Thus there is no generic-parameter or numerical smoothness assumption.
For m2, Q2=x^2-4hx+6h^2. The earlier source-linear note used
`a(x^2-4hx)+B0_old`; the two conventions agree after
`B0_old=B0+6ah^2`. This changes no function or theorem and prevents an
incorrect identification of the special value or boundary intercept.

## 3. An actual alternative affine chart, not only a function field

On W_m minus E1, u=1/(x-h) is globally regular. Its second-chart
expression is `u=r/(1-hr)`, regular on D and zero there.
Construct an inverse from the affine plane with coordinates (c,u).

On its principal open u!=0 use the original source chart:

```
x=h+1/u,
t=T_c(u):=u^(2m)[(c-B0)/a-Q_m(h+1/u)].                    (2)
```

T_c(u) is a polynomial in c,u, since deg Q_m=m<=2m. On the principal
open 1+hu!=0 use the boundary chart:

```
r=u/(1+hu),
b=b_c(u):=(B0-c)/a*(1+hu)^(2m)
 -sum_(j=m+1)^(2m) binom(2m,j)(-h)^j
       *u^(j-m)*(1+hu)^(3m-j).                           (3)
```

All powers in b_c are nonnegative, and b_c is polynomial in c,u.
The two principal opens cover A2, since u and 1+hu have no common zero.
Equations (2) and (3) obey the exact chart transition and satisfy F=c;
they are obtained by solving the displayed affine-linear expressions for
t and b. They avoid E1: x=h+1/u cannot equal h, and r=u/(1+hu) cannot
equal 1/h. Conversely their composites with (F,u) are identities on
both original charts. They therefore give an isomorphism of varieties.
No source translation has been called a global surface automorphism.

At u=0 the inverse is the actual D point `r=0,b=(B0-c)/a`. The value
u=-1/h, where the r formula is unavailable, is covered by (2) and has
x=0 and b=0. These are the two cheap endpoint tests for the chart switch.

For c!=f0 the fibre has no point on E1, so the new chart identifies it
with the entire affine u-line. For c=f0, the exact factorization
`g=a z K` has `K mod z=kappa !=0`. Thus its two component ideals are
comaximal and both components are reduced. The residual K is irreducible:
its t-coefficient z^e is coprime to its constant coefficient, which is
nonzero at z=0. It is linear in t and primitive. In the new chart its
closure E2 is exactly the line c=f0. E1 itself is the original affine
line x=h. This proves every fibre statement, including the boundary fill.

## 4. The finite power map and its original-source failure

Direct differentiation gives

```
D_F G0=1,
J_(x,t)(F,u)=a(x-h)^(2m-2),
omega=(u^(2m-2)/a) dF wedge du,
G0=u^e/(a e).                                             (4)
```

The last expression identifies the new-chart map with
`(c,u) -> (c,u^e/(a e))`. Its global coordinate ring is free over
`C[c,u^e/(a e)]` with basis `1,u,...,u^(e-1)`. It is finite flat of
degree e, and its Jacobian is u^(e-1)/a. For e>1 the only ramification
is u=0 with the declared index; e1 is an isomorphism.

This map is NOT defined regularly on all W_m: G0 has a pole of order e
along E1, which belongs to the original U0. The new affine chart removes
E1 and retains D; the original source removes D and retains E1. The
precise exchanged divisors matter. A finite map on the new chart is not
a finite Keller envelope for the original source. The field degree
`[C(x,t):C(F,G0)]` is e, but this numerical identity alone omits those
regularity and volume data.

## 5. Complete source-module annihilator, with no mate-degree bound

The field identity
`C(x,t)=C(F)(x)`, obtained by solving t=(F-B)/A, is exact. Holding F
fixed turns the derivation into `-a(x-h)^(2m) partial_x`. Its rational
constant field is exactly C(F). Therefore every rational solution of
`D_F P=q(F)`, for q in C[T], is

```
P=q(F)G0+H(F),          H in C(T).                         (5)
```

There is a polynomial witness at the stated order:

```
P_e=g^e G0=(a^(e-1)/e)K^e in C[x,t],
D_F P_e=g^e.                                             (6)
```

This also supplies every multiple of g^e. To prove necessity, let
`ell=ord_(T=f0) q`, with q nonzero, and suppose ell<e. At E1, g is a
uniformizer because K(0)=kappa is nonzero. Consequently q(F)G0 has pole
order e-ell there. Any cancellation by H(F) in (5) requires H to have
a pole at f0 of precisely that order. At E2, g is also a uniformizer,
whereas G0 is regular (z is invertible). Thus q(F)G0 is regular on E2,
and the same H(F) introduces an uncancelled pole. No polynomial P can
exist. This proves the full annihilator ideal and the exclusion of
polynomial mates, without restricting a proposed degree or coefficient.

The top scalar principal-part coefficient on E1, in the actual parameter
g, is `(a kappa)^e/(a e) !=0`; on E2 the principal part is zero. Lower
coefficients are unnecessary for the annihilator proof. The full
component-labelled torsion theorem applies because F is smooth on the
affine source and its rational constants are C(F). The unit class has
support only at f0. Its canonical connection differentiates scalar
principal parts, increasing the highest pole order by one with nonzero
factor in characteristic zero. This proves the entire derivative claim.

## 6. Sharp failure boundary and controls

The hypothesis h!=0 matters for m>=2. At h=0 the same formulas become

```
F=a x^m(1+x^m t)+B0,      g=F-B0,
G0=1/[a(2m-1)x^(2m-1)].
```

Now F is critical along x=0, and that fibre component has multiplicity
m. The other component `1+x^m t=0` remains reduced. The same direct
two-component argument, without invoking the smooth torsion theorem,
gives exact annihilator `(g^2)` for every m>=2: the pole order of q(F)G0
on x=0 is `2m-1-m*ord_(B0)q`, and cancelling it by a fibre-constant
correction necessarily introduces a pole on the other component.
The witness g^2 G0 is polynomial. Thus multiplicity, not simply the
denominator's exponent, controls the repair threshold after this
degeneration. For m1, h0 is harmless and the exact exponent is still1.

At m1 in the main family, E1/E2 are still distinct, the unit order is1,
and the new-chart finite map is degree1. This is a compulsory control
against calling every rational mate nonbirational. For every m, the
affine polynomial pair (t,-x) is a positive Jacobian-one control; its
second component does not extend across D.

## 7. Connection contract and decisive next question

The five live concepts are source polynomiality, global submersion,
component-labelled repair, the alternative affine chart, and boundary
ramification. The construction sends the complete polynomial part of
`(x-h)^(2m)/x^m` to a genuine global submersion, then uses u=1/(x-h) to
identify its alternative source chart. It preserves the actual two-form
and component labels. A field-only view loses which affine line was
removed. The equal integer e measures the boundary power-map degree and
the unit's repair order only after those distinct geometric roles are
retained; this note asserts no universal equality for other functions.

This is a scoped realization on W_m, not a classification of all its
submersions or polynomial endomorphisms. A remaining useful question is
whether an actual moving-source response maps into this fixed-F
component module while preserving D_F and C[x,t]. No such entry is
implied by an equal pole order or response dimension.

## 8. Finite exact reproduction

The companion source imports no prior producer. Its declared universe
includes both integer binomial identities through m=40, complete symbolic chart
identities for finitely many m, independent literal inverse-chart
parameter banks, the reduced/multiple special-fibre controls, exact
polynomial witnesses, and the actual power-map Jacobian. The unbounded
claims use the proofs above. The adjacent JSON records the complete
finite bank and every exponent; no floating arithmetic or solver status
is proof authority.

```
python -B 04-computation/continuing10_20260908_dg_unbounded_torsion.py
python -B -O 04-computation/continuing10_20260908_dg_unbounded_torsion.py
```

Both normal and optimized runs pass **546 always-active exact gates**,
with byte-identical actual LF output. The certificate regenerates unchanged.
The primary source writes its standard `_certificate.json` under results
when relocated to `04-computation`, or beside itself in this outside bundle.
Frozen SHA256:

```
source 2c6e575521e44b971877455d1984e4286611c0df2e15bec73ef17941d6a6bd49
output 35272f408c20715cd0b84df71ecc6a48bf92c841a5a2c103f0077ed0755cfd53
certificate 3517b35f819399b7aaca886bc58cbac6653506f9b387c2409d708c3772fa6ef9
```

The [independent audit](continuing10_20260908_all_m_torsion_audit.md) accepts the complete analytic theorem and its sharp parameter boundary. Filing changes only status and routing; source, output and certificate bytes are unchanged.
