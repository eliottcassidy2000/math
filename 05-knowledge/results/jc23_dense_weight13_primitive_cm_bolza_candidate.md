# Dense exact-weight-thirteen primitive-CM plus Bolza exclusion

**Status: PROOF CANDIDATE / VERIFIED-EXACT / NOT CANON / NO THEOREM ID.**

This file records a complete candidate proof for independent audit.  It does
not reserve a theorem namespace and is not a proved dependency until promoted
through the repository validity gate.

## 1. Candidate statement and firewall

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam of
THM-3992/3997.  Put

```text
A=[p^5y]H,       B=[p^2y^3]H,
U=[p^6]H,        Z=[y^4]H.                              (1)
```

> **Candidate theorem.** The complete exact-weight-thirteen locus
>
> ```text
> A*B*U*Z*(A+B) != 0                                   (2)
> ```
>
> contains no nonautomorphic planar Keller pair.

Every other lower coefficient is arbitrary subject only to the inherited
normal form.  No further collision unit is needed: `A=B` is safe, as are all
relations among lower coefficients.  The five walls in `(2)`, other cells,
entry into the inherited seam, `JC(2)`, and `DC(2)` remain open.

The inheritance pass is:

- closest mechanism: THM-4222's primitive-CM total-`Hom` obstruction;
- hostile: THM-4218's `j=0` positive-genus tail, which shows that a genus
  inventory cannot be replaced by the main face alone;
- corrected near miss: the THM-4222 internal slope probes cannot be reused;
  on the new `BD` edge the primitive transverse step is `(3,2)`, and the
  correct chain has two intermediate components, not five;
- least-used sidecar: THM-4012's exact Bolza splitting and CM mismatch.

## 2. Complete source and the three lower planes

Use `s=XT`, `p=T+s^2`, `y=sp`, `t=p-s^2=T`.  With arbitrary lower
coefficients renamed only to keep the four units in `(1)` visible, the full
polynomial is

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y+zeta_3 y^3
  +upsilon_5 p^5+xi_10 p^2y^2+alpha_11 p^4y+beta_11 py^3
  +U p^6+omega_12 p^3y^2+Z y^4+A p^5y+B p^2y^3,

K=2848/45-(7/6)Delta.                                  (3)
```

Indeed, enumerating `0<2i+3j<=13` and deleting the forbidden `y,py` rows
gives exactly eighteen monomials; the only weight-thirteen rows are `p^5y`
and `p^2y^3`.

For `Q=q^-1`, set

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2.                              (4)
```

On its torus `p=s^2` is impossible, since `(4)` would become
`-Qs^2/2!=0`; hence the standard inverse through `t=p-s^2`, `X=s/t`
recovers the actual source function field.

A monomial `p^iy^j` contributes endpoints

```text
(j+2,i+j,1),                  (j,i+j+1,1).             (5)
```

The three lower planes are

```text
M: nu=(r+2l-2)/13,
V: nu=(l-1)/6,
T: nu=(r+l-2)/8.                                     (6)
```

For the two endpoints in `(5)`, their respective gaps above these planes are

```text
M: ((13-2i-3j)/13, (13-2i-3j)/13),
V: ((7-i-j)/6,       (6-i-j)/6),
T: ((8-i-2j)/8,      (9-i-2j)/8).                     (7)
```

All are nonnegative through weight thirteen.  Equality occurs exactly as
follows:

```text
M: p^5y,p^2y^3;       V: p^6,p^5y;       T: y^4,p^2y^3. (8)
```

Thus the four required rows already own every face.  Even deleting the
coincident midpoint `(3,6,1)` leaves all three planes.  Adding or cancelling
any point with positive gap in `(7)` cannot alter the lower hull.  This is the
all-lower-coefficient proof; the exhaustive `2^14` lower-row and `2^8`
collision audits are hostile controls, not the logical substitute for `(7)`.

The eight coincident support points and their coefficients are

```text
(2,3,1): K+1376/135,          (2,4,1): Theta-Delta,
(2,5,1): xi_10-upsilon_5,     (2,6,1): omega_12-U,
(3,4,1): zeta_3-eta,          (3,5,1): beta_11-alpha_11,
(3,6,1): B-A,                 (4,5,1): Z-omega_12.     (9)
```

Only `(3,6,1)` lies on a lower face, and it is the midpoint of a retained
edge.  Hence none of the eight differences is a unit condition.  The sole
additional collision unit is `A+B`, arising from a repeated edge root and
the movement of the main intersections to the boundary below.

## 3. Polygons, faces, and every edge scheme

The face polygons and Pick ledgers are

| face | polygon | `(2Area,boundary,interior)` |
|---|---|---:|
| `M` | `(0,1),(2,0),(5,5),(1,7)` | `(39,5,18)` |
| `V` | `(0,1),(1,7),(0,7)` | `(6,8,0)` |
| `T` | `(2,0),(6,4),(5,5)` | `(8,6,2)` |

Their equations, up to nonzero monomials, are

```text
g_M=(S^2-P)(1-A S P^6-B S^3P^5)=R*C,
g_V/P=-1+P^6(U+A S),
g_T/S^2=1-Z S^4P^4-B S^3P^5.                          (10)
```

The vertical component is rational.  The tail is smooth on the torus because,
after `x=SP`, differentiation of

```text
1-Zx^4-Bx^3P^2                                          (11)
```

with respect to `P` gives `-2Bx^3P`, which is nonzero there.

The global polygon is

```text
(0,1),(2,0),(6,4),(5,5),(1,7),(0,7),
(2Area,boundary,g)=(53,15,20).                         (12)
```

The intrinsic outer packet, omitting the affine `s=0` divisor, is

```text
(12,12,8,6,2,2,2,2,1),
sum=47,                         defect=38=2*20-2.       (13)
```

Write the vertices in `(12)` as `A_0,...,F_0`.  The six outer and two
internal edge schemes, up to units, monomials, and reversal, are

```text
A_0B_0: X-1,                    B_0C_0: 1-ZX^4,
C_0D_0: -Z-BX,                  D_0E_0: (X-1)(AX+B),
E_0F_0: A+UX,                   F_0A_0: U-X^6,
A_0E_0: AX-1,                   B_0D_0: 1-BX.           (14)
```

Their nonconstant discriminants are

```text
-256 Z^3,                 (A+B)^2,                 46656 U^5. (15)
```

Consequently `(2)` is exactly the reduced, corner-avoiding edge gate.  In
particular there is no `Z+B`, `A+U`, or `A-B` wall.  When `A=B`, the fourth
scheme is `A(X-1)(X+1)`, so the cancelled midpoint in `(9)` is harmless.
The `Z+B` boundary warning in THM-4012 concerned intersections of its Bolza
factor with a rational factor on one highest face.  Here the Bolza curve is a
separate tail face and its unique internal attachment is the reduced scheme
`1-BX`; thus that warning does not create a hidden `Z+B` unit in this model.

The only main-face torus singularities are intersections of `R` and `C`:

```text
P=S^2,                         1-(A+B)S^13=0.           (16)
```

Their gradient determinant is `-13(A+B)S^12`, giving thirteen distinct
transverse points.  The `C` torus gradients have exponent determinant
`1*5-3*6=-13`, so `A*B!=0` makes `C` smooth.  Equations `(11)`, `(14)`, and
the derivative `A P^6` on `V` cover all remaining torus and compactified
boundary points.

## 4. Exact regular model

Use the common base change

```text
Q=sigma^312.                                             (17)
```

The integral lower heights and primitive normals are

```text
M:24r+48l-48,           (24,48,-1),
V:52l-52,               (0,52,-1),
T:39r+39l-78,           (39,39,-1).                    (18)
```

Every face therefore has multiplicity one.  The six outer three-dimensional
edge gcds equal the planar lengths `1,4,1,2,1,6`; both internal graph edges
are primitive.  The outer slope ledger is

```text
24,-39,-39,-24,-52,0,                                  (19)
```

and adds no chain.

On `A_0E_0`, the slopes are `-48>-49>-50>-51>-52`, giving three
intermediate rational multiplicity-one components.  On `B_0D_0`, use the
primitive transverse functional

```text
L*=-5r+3l+10,
L*(B_0)=0,                     L*(3,2)=1.               (20)
```

The two adjacent slopes are `120` and `117`, so the determinant-one sequence

```text
120>119>118>117                                         (21)
```

has two intermediate rational multiplicity-one components.  The tempting
probe `(1,-1)` has `L*=2` and would incorrectly double the transverse step.

In the main chart set

```text
s=sigma^-24 S,        p=sigma^-48 P,
H_M=sigma^312 H(sigma^-24S,sigma^-48P),
U_0=S^2-P,            V_0=(1-H_M)/S^2.                 (22)
```

Multiplying `(4)` by `sigma^48` gives the exact equation

```text
U_0 V_0=sigma^312/2.                                   (23)
```

Thus each point in `(16)` is an `A_311` smoothing.  Resolving it inserts only
rational multiplicity-one components.

The side charts are

```text
V: s=S,             p=sigma^-52P;
T: s=sigma^-39S,    p=sigma^-39P.                      (24)
```

After multiplying `(4)` by `sigma^52` and `sigma^78`, respectively, their
special fibres are exactly the last two faces in `(10)`.  For the tail, the
entire scaled equation is

```text
(S^2-sigma^39P)(1-H_T)-sigma^312S^2/2=0,
H_T=sigma^312H(sigma^-39S,sigma^-39P).                 (25)
```

The face smoothness, eight reduced edge schemes, primitive denominators,
and explicit resolutions meet Dokchitser's face/edge regularity gate through
the audited THM-4045 interface.

The two moving outer restrictions are

```text
q-1/2=K W^2+zeta_3 W^3+Z W^4,
q=H(p,0),                  deg_p H(p,0)=6.              (26)
```

They are squarefree over `C(q)`: a common derivative root is constant over
the coefficient field and cannot have the transcendental value `q`.  All
other generic edge restrictions are linear or the separated scheme in
`(14)`.  Hence the generic completion is smooth and is the actual source
completion by the inverse after `(4)`.

## 5. The two positive-genus components

### 5.1 Primitive cyclotomic main component

The nonrational main normalization is

```text
C:1-A S P^6-B S^3P^5=0.                               (27)
```

Put `x=ASP^6`.  Then

```text
P^13=(B/A^3)x^3/(1-x).                                 (28)
```

This is a cyclic degree-thirteen cover with branch residues
`(3,12,11)` at `0,1,infinity`, so Riemann--Hurwitz gives `g(C)=6`.
The Newton triangle `(0,0),(1,6),(3,5)` has interiors

```text
(1,2),(1,3),(1,4),(1,5),(2,4),(2,5).                  (29)
```

For a primitive thirteenth root `zeta_13`, the deck action

```text
tau:(S,P) |-> (zeta_13^7 S,zeta_13 P)                 (30)
```

has holomorphic character set

```text
Phi={5,6,9,10,11,12}.                                  (31)
```

Its negatives are `{1,2,3,4,7,8}`, so `H^1(C)` contains every nontrivial
character once and `Q[tau]=Q(zeta_13)` has degree `12=2g`.  Exact
multiplication by all twelve units gives

```text
{u in (Z/13)^*:uPhi=Phi}={1}.                          (32)
```

Thus the CM-pair is primitive.  Milne, *Complex Multiplication*, Sections
1.9--1.10 and Proposition 3.13, implies that `J(C)` is simple.  Since its
dimension is six,

```text
Hom(J(C),E)=0                                           (33)
```

for every elliptic curve `E`.

### 5.2 Bolza tail

For the tail in `(10)`, put

```text
x=SP,                         Y=x^2P.                  (34)
```

These are birational torus coordinates, and the normalization becomes

```text
B Y^2=x-Zx^5=x(1-Zx^4).                                (35)
```

The discriminant of the quintic is `-256 Z^3`, so it is squarefree under
`(2)` and has genus two.  Choosing `u,v in C*` with
`Z=-u^4`, `B=v^2/u` identifies `(35)` with

```text
y^2=x(x^4+1).                                          (36)
```

THM-4012 proves exactly that this Bolza Jacobian is isogenous to
`E_8000^2` and that `Hom(E_8000,E_0)=0` for the target curve
`E_0:Y^2=X^3+1`.  Hence

```text
Hom(J(Bolza),E_0)=0.                                   (37)
```

This is the decisive difference from THM-4218's dangerous `j=0` tail.

## 6. Genus completeness and degree zero

After contracting rational paths, the core has vertices

```text
R,C,Bolza,V,
```

thirteen `R--C` paths, one `C--Bolza` path, and one `C--V` path.  Therefore

```text
b_1=15-4+1=12,
g_special=6+2+12=20.                                  (38)
```

This equals the global Pick genus `(12)`, proving component and packet
completeness.  Every component outside `C` and the Bolza tail is rational.

Since `q=sigma^-312`, scale the inherited target by

```text
A_target=sigma^-104X,                C_target=sigma^-156Y. (39)
```

Its exact equation becomes

```text
Y^2=X^3+1-(3a^2/4)sigma^208X-(a^3/4)sigma^312,         (40)
```

with smooth special fibre `E_0:Y^2=X^3+1`.

A hypothetical Keller pair gives a finite nonconstant generic morphism to
this target.  Resolve its rational extension over the regular source model.
Every exceptional component is rational.  Riemann--Hurwitz makes every
rational component map constantly; `(33)` and `(37)` do the same for the two
positive-genus components.

For a relative target line bundle `L` of degree one, proper-flat degree
conservation and the multiplicity-one ledger give

```text
deg(phi_generic)=sum_i m_i deg(phi_i^*L)=0.             (41)
```

This contradicts generic finiteness and proves the candidate statement,
relative to the inherited seam and cited dependencies.

## 7. Verification and status discipline

The exact certificate checks:

1. all eighteen monomials, `2^14=16384` lower-row subsets, all eight
   collision points and `2^8=256` collision patterns;
2. the required-skeleton and analytic simultaneous-cancellation proof;
3. all face/global Pick ledgers, packet, and eight edge schemes;
4. exact unit gate `A*B*U*Z*(A+B)`, including the `A=B` hostile;
5. primitive normals, edge denominators, corrected slopes and chain counts;
6. exact main, vertical, and tail charts and thirteen `A_311` models;
7. cyclic identity, CM characters and stabilizer, Bolza identity and
   squarefreeness;
8. generic carrier squarefreeness and target scaling; and
9. complete component/graph genus inventory.

Artifact hashes:

```text
script sha256 280d2c2515fe4d1b3f5f229947319d0922b3982d83ad88fd39cd498b78e8404b
output sha256 fd9b297486e923fe6946eb7a98fb5bd09f7d22080fc354c32202f993c4d1dc96
independent audit script sha256 0eac79869bfbef796774531a1cc3607e1badbef8b165db473b18bf6b2a90280c
independent audit output sha256 774e328d0fda469e51adc956daef4df4e72799ae8700e6890f63e2601ed16d58
```

Replay with

```bash
python3 -B 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate.py
python3 -B -O 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate.py
PYTHONHASHSEED=4213 python3 -B 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate.py

python3 -B 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate_independent_audit.py
python3 -B -O 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate_independent_audit.py
PYTHONHASHSEED=4213 python3 -B 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate_independent_audit.py
```

The corresponding streams must byte-match the two frozen outputs.  The
independent audit is standard-library-only and reconstructs the monomial,
gap, collision, polygon, edge-coefficient, packet, lattice-slope,
Chevalley--Weil, Bolza, valuation, and genus ledgers without importing SymPy
or the primary implementation.

The exact certificate does **not** computationally prove Milne's theorem,
THM-4012's Hom calculation, the inherited Keller-map interface, or
proper-flat degree conservation.  Those are the explicit proof/citation
layer.  Promotion remains forbidden until an independent audit accepts this
full regular-model and collision firewall.
