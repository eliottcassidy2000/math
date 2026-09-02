---
id: THM-4334
title: "Beta-zero exact-weight-twelve endpoint-wall extinction"
status: >
  PROVED RELATIVE TO THM-4292/4297/4327 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. In the inherited reduced (2,3), exact-weight-twelve seam, the
  endpoint stratum Z=beta_11=0 and U*W*zeta_3!=0 is extinct with Lambda=U+W
  arbitrary. The beta-zero face is a connected genus-three cyclic-nine
  Kummer curve; together with the central genus-three component and graph
  b1=11 it accounts for the full genus seventeen. Their good-form orders are
  34 and 27. At Lambda=0 the sole A23 contact splits through C2 or C3 before
  the t^6 correction, so every component is Keller-constant. THM-4337 closes
  the complementary zeta_3=0 wall away from beta_11=0; seam entry, JC(2),
  and DC(2) remain open.
source: root + jc boundary agents / planar-Jacobian continuation session, 2026-09-02
depends_on:
  - THM-4292-lambda-zero-repeated-face-keller-extinction
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
related:
  - THM-4292-lambda-zero-repeated-face-keller-extinction
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4328-seam-covariant-student-stein-face-visibility
  - THM-4337-zeta-zero-exact-weight-twelve-endpoint-wall-extinction
mistake_firewall:
  - MISTAKE-487
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_z0_beta0_endpoint_extinction_thm4334.py
primary_output: 05-knowledge/results/jc2_m12_z0_beta0_endpoint_extinction_thm4334.out
primary_script_sha256: 3d6dce090ce3d7424427e09b3b76aa50385eee9e030b670be5f1b3595fbd32a2
primary_output_sha256: 912ce9c8718913468641691cf37ec67a6d8e21bc5001210749f15e7c9ae72ea8
independent_audit_script: 04-computation/jc2_m12_z0_beta0_endpoint_extinction_independent_audit_thm4334.py
independent_audit_output: 05-knowledge/results/jc2_m12_z0_beta0_endpoint_extinction_independent_audit_thm4334.out
independent_audit_script_sha256: e566081ddfc552a3b65b465741f14191070e5efbffcacce9b6e1611c62736da2
independent_audit_output_sha256: 07aabbd8e7b49ad5a761b6cbabe91fc33a6d5667c88bb1a1e6bf6679250da988
hash_basis: raw LF bytes
audit: >
  PASS. The primary exact hull path exhausts all 16,384 optional-row and
  hostile aggregate-cancellation states, checks both face fields, Pick and
  graph genus, the complete edge list, good-form orders, and the shortened
  beta-zero A23 ladder. A standard-library clean-room path independently
  combines a dual-slack/uniquely-owned-anchor certificate with its own hull,
  Kummer, edge, contact, packet, and hostile-boundary checks. A separate
  hostile referee audited the birational recoveries, toric edge orbits,
  local/global base conversion, and proper-flat bridge. Normal and optimized
  streams byte-match both frozen outputs.
---

# THM-4334 -- Beta-zero exact-weight-twelve endpoint-wall extinction

**PROVED RELATIVE TO THM-4292/4297/4327 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THE `Z=beta_11=0`, `U*W*zeta_3!=0` ENDPOINT STRATUM IS EXTINCT FOR
ARBITRARY `Lambda`; SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement

Work over an algebraically closed field of characteristic zero in the reduced
`(2,3)`, exact-weight-twelve source family inherited from THM-4230 and
THM-4327.  Impose

```text
Z=0,                 beta_11=0,                 U*W*zeta_3!=0.       (1)
```

Retain every other allowed lower coefficient, and put `Lambda=U+W`.  The
theorem is:

> No nonautomorphic planar Keller pair lies on `(1)`, for arbitrary
> `Lambda`.  Equivalently, every component of the audited special source
> fibre carries a constant map to the good elliptic target, contradicting
> proper-flat conservation of the positive generic map degree.

This is a relative endpoint theorem.  It does not prove that a hypothetical
counterexample enters this exact seam, and it does not prove `JC(2)` or
`DC(2)`.

The inheritance pass is:

- closest proved mechanism: THM-4327's good-target differential extinction
  and proper-flat degree conservation;
- canonical hostile: MISTAKE-531's `Lambda=0` double top edge, across which a
  simple-root response packet may not be transported;
- corrected near miss: MISTAKE-540, which requires genus to be computed in
  the invariant Kummer function field rather than on a ramified root cover;
- least-used relevant sidecar: the labelled toric edge orbit, which separates
  the old `A_23` contact from a new adjacent-face node.

The live concept board was

```text
{owner deletion, lower hull, invariant Kummer field, edge orbit,
 good differential, repeated splitter, arithmetic-genus ledger}.           (2)
```

## 2. The exact two-face lower complex

The complete labelled row universe consists of the sixteen allowed
`p^i y^j` with `0<2i+3j<=12`, excluding the inherited forbidden rows
`(i,j)=(0,1),(1,1)`.  A row contributes lifted points

```text
(j+2,i+j,1),                    (j,i+j+1,1),               (3)
```

together with `(2,0,0),(0,1,0),(2,0,1)`.  On `(1)`, the rows for
`Z=(0,4)` and `beta_11=(1,3)` are absent.  The rows for
`p,p^2,p^3,U=(6,0),W=(3,2),zeta_3=(0,3)` are required, leaving eight
optional rows.

To overcount every possible aggregate-coefficient cancellation, independently
delete each active point among

```text
(2,3,1),(2,4,1),(2,5,1),(2,6,1),(3,4,1),(3,5,1).          (4)
```

Notice that `(4,5,1)` is no longer deletable: after `Z=0` it is uniquely
owned by `W!=0`.  Exact enumeration of all

```text
2^8 * 2^6 = 16,384                                             (5)
```

hostile states gives the same two lower planes in every case:

```text
M =(1/12,1/6,-1/6),             N=(2/9,1/9,-4/9).          (6)
```

The independent certificate does not rely on that exhaustive hull routine.
It checks the exact slack of all twenty-six possible support points against
both planes, exhibits five uniquely owned rank-two anchors which survive
every deletion in `(4)`, and proves that their projected cells cover the
lower support.  Thus optional rows and aggregate cancellations can neither
delete `M,N` nor insert a third face.

The potential equality supports are

```text
M: (0,1,0),(0,7,1),(2,0,0),(2,6,1),(4,5,1),
N: (2,0,0),(4,5,1),(5,3,1).                               (7)
```

The point `(2,6,1)` may cancel, but the other vertices retain the same
rank-two `M` polygon.  With global base `Q=Sigma^36` (`L=36`), the scaled
slope triples are

```text
36 M=(3,6,-6),                       36 N=(8,4,-16).        (8)
```

The height coefficient in each three-dimensional normal is `-1`, so both
normals are primitive and both face multiplicities are one.  The base is also
divisible by six, making the good elliptic target integral.

## 3. Face normalizations and invariant genus

Up to torus monomials, the face equations are

```text
G_M=(S^2-P)C,             C=1-UP^6-WS^2P^5,
G_N=S^2 B,                B=1-WS^2P^5-zeta_3 S^3P^3.      (9)
```

Write `R:S^2=P`.  The Newton polygons of `M,C,N` and of the complete lower
boundary have Pick ledgers

```text
M:      (2Area,boundary,interior)=(36,10,14),
C:                                  (12, 8, 3),
N:                                  ( 9, 5, 3),
global:                             (45,13,17).             (10)
```

For `C`, put `x=P^(-1)` and `y=WS/P`.  Equation `(9)` becomes

```text
y^2=W x(x^6-U).                                             (11)
```

Here `P=x^(-1)` and `S=y/(Wx)`, so this is a birational recovery of the
actual central face field, not a cover of it.  The degree-seven right side is
squarefree because `U*W!=0`, so `C` is a connected smooth genus-three
normalization.

For `N`, remain in the invariant face field and set

```text
x=WS^2P^5,                    1-x=zeta_3 S^3P^3.
```

Eliminating `S` gives the degree-nine Kummer equation

```text
P^9=(zeta_3^2/W^3) x^3/(1-x)^2.                            (12)
```

Conversely,

```text
S=W(1-x)P^2/(zeta_3 x),                                    (12a)
```

so `k(N)=k(x,P)` with relation `(12)`.  This is a birational presentation of
the actual face field, not merely a quotient relation or a ramified ambient
coordinate.

Its valuation vector at `x=0,1,infinity` is `(3,-2,-1)`.  Since the gcd of
these valuations with nine is one, the cover is connected.  Riemann--Hurwitz
gives

```text
2g(N)-2=-18+(9-3)+(9-1)+(9-1)=4,
g(N)=3.                                                       (13)
```

This is the invariant quotient computation required by MISTAKE-540; no root
coordinate is counted as part of the source curve.

The exponent determinants for `C` and `N` are respectively

```text
det((0,6),(2,5))=-12,             det((2,5),(3,3))=-9.      (14)
```

Thus their logarithmic gradients have no common torus zero.  At their generic
points the exact Euler remainders are

```text
P C_P-5C=-(UP^6+5),
P B_P-3B=-(2WS^2P^5+3).                                    (15)
```

The first is nonzero in `k(C)`--at `S=0,C=0` it equals `-6`.  If the second
vanished in `k(N)`, then both `WS^2P^5=-3/2` and
`zeta_3S^3P^3=5/2` would be constant; the nonzero determinant in `(14)` would
make `S,P` constant, contradicting that `N` is a curve.  On the central
branch,

```text
(G_M)_P|C=(S^2-P)C_P,
```

and `S^2-P` is nonzero in `k(C)` because `C` is a distinct irreducible
factor.  On the merged branch,

```text
(G_N)_P|N=S^2B_P,
```

with `S` a torus unit.  Therefore `(15)` proves that the full face derivative
`G_P` is nonzero at both positive-genus generic points.

## 4. Complete edges and the off-contact graph

The complete nonmonomial edge schemes are

```text
X-1,                         1-zeta_3 X^3,
zeta_3+W X,                  (X-1)(UX+W),
U-X^6,                       1-WX.                         (16)
```

For example, the top terms of `G_M` are

```text
S^4P^5[-W+(W-U)X+UX^2]
 =S^4P^5(X-1)(UX+W),                 X=P/S^2,              (17)
```

and the shared `M--N` edge is `S^2(1-WS^2P^5)`.  The remaining two
`N` edges follow directly from
`S^2-WS^4P^5-zeta_3S^5P^3`.

Every edge in `(16)` is active and squarefree under

```text
U*W*zeta_3*Lambda!=0.                                      (18)
```

Indeed the top quadratic has discriminant `(U+W)^2`; the other derivatives
give only the already excluded gates `U=0`, `W=0`, or `zeta_3=0`.  The common
`M--N` segment from `(2,0)` to `(4,5)` is primitive, so `1-WX` gives one
transverse `C--N` node.  The initial form of `R` on that divisor is a torus
monomial, hence `N` does not meet `R` there.

When `Lambda!=0`, substituting `P=S^2` into `C` gives

```text
C|R=1-Lambda P^6.                                          (19)
```

Thus `R` and `C` meet transversely in twelve points.  The component graph has
three vertices, thirteen edges, and first Betti number eleven.  Consequently

```text
g(C)+g(N)+b_1=3+3+11=17,                                  (20)
```

equal to the global Pick genus in `(10)`.  No positive-genus component is
missing from the off-contact inventory; all toric subdivisions add only
rational chains.  In branch parameters each `R--C` smoothing has the exact
form

```text
u_branch v_branch=Sigma^36*(unit),                          (20a)
```

the `L=36` instance of THM-4327 equation (9d).  The primitive `C--N` edge is
the ordinary-node chart supplied by the same audited face/edge interface.

## 5. The `Lambda=0` collision survives `beta_11=0`

MISTAKE-531 forbids transporting the off-contact response packet, so no such
packet is used.  Put `W=-U`.  In the top-infinity variables of THM-4327,

```text
Q=sigma^12,                 sigma=Sigma^3,
t=sigma b,                  s=v(sigma)>0,
nu=v(b)>0,
q=r-1,                      b=1/S,
r=P/S^2.                                                    (20b)
```

The main closure is

```text
Ctilde=b^12-U r^5(r-1),        partial_r Ctilde(0,1)=-U.   (21)
```

Hence `R` and `C` have one length-twelve intersection at `r=1`, a two-branch
`A_23` contact.  In THM-4297's notation

```text
A=2U+W=U,                       U_eff=A/2!=0.               (22)
```

It remains to check the local tail argument directly, because THM-4297's
global theorem statement assumes a different coefficient gate.  THM-4292
equations (3)--(5) retain the complete local term

```text
t(alpha_11 r^5+beta_11 r^4).                               (23)
```

The prepared equation has the form

```text
F=q ell+q^2 V(q,t)-t^12/2,                 V(0,0)=A=U!=0.  (23a)
```

Before a repeated discriminant root, the calculation uses only this nonzero
quadratic coefficient and gives positive form order at least `3s+5nu`.  In
the range `nu<s`, any coefficient/`b^12` cancellation ends in a quadratic
whose discriminant is `x^2+2A`, with simple roots because `A!=0`.  On the
balanced face `nu=s`, write `c` for the `t^6` coefficient of `ell` after the
lower rows have vanished.  Its discriminant is `(c-X^6)^2+2A`: a nonzero
multiple root would again force `A=0`, and the possible root at `X=0`, where
`c^2+2A=0` and hence `c!=0`, is exactly the passage to `nu>s`.  Thus only
that last range needs a repeated-tail audit.

More explicitly, THM-4297 equations (18)--(19) write the only difference
from the `W=0` quadratic model as

```text
Delta F=-(W/2)r^4q^3,                    q=r-1.             (24)
```

After `q=t^6y` and division by `t^12`, `(24)` begins at `t^6`.  Through
order four, the exact moving-critical coefficients remain

```text
C_1=alpha_11/c^2,
C_2=upsilon_5/c^2+8/(3c),
C_3=eta/c^2,
C_4=(Delta+32/9-3c)/c^2.                                  (25)
```

There is no `t^5` coefficient.  If all four entries vanish, the terminal
`b^12q` splitter still precedes `(24)`, since for positive valuations
`nu>s`

```text
6(nu-s)<6(s+nu).                                           (26)
```

The exhaustive Keller-form lower orders are therefore exactly

```text
6s+2nu, (5s+9nu)/2, 2s+4nu, (3s+7nu)/2, s+3nu.            (27)
```

Every coefficient in `(27)` is positive.  In fact the present coefficient
gate makes the deepest part of the ladder strictly shorter.  Repetition in
the range `nu>s` requires

```text
c_1=alpha_11+beta_11=0,             c_3=eta+zeta_3=0.       (27a)
```

Because `beta_11=0`, the first equation gives `alpha_11=0`, so `C_1=0`.
Because `zeta_3!=0`, the second gives `eta=-zeta_3!=0`, and hence
`C_3=eta/c^2!=0`; here the repeated discriminant has `c!=0` because `A!=0`.
Therefore `C_2` splits first if nonzero, and otherwise `C_3` necessarily
splits.  Their gaps `2(s+nu)` and `3(s+nu)` are both strictly earlier than
the `t^6` correction, and their positive form orders are `2s+4nu` and
`(3s+7nu)/2`.  The `C_4` and terminal rows in `(27)` are valid consistency
controls but are not needed for extinction on `(1)`.  Since `L=36` is three
times the local base twelve (`sigma=Sigma^3`), every local order is merely
multiplied by the positive ramification index three.

The new `N` face creates no additional collision.  Its shared edge with `M`
is primitive with simple scheme `1-WX`; only `C`, not `R`, has a root there.
The `A_23` root lies in the open orbit of the distinct top-edge divisor.
Distinct toric divisors meet only at a fixed point, whereas both roots are
finite and nonzero, so the two events cannot coincide.  The arithmetic-genus
ledger is now

```text
0+3+3 + delta(A_23)+delta(node) - 3 + 1
 =6+12+1-3+1=17.                                           (28)
```

Thus the contact and the `C--N` node exhaust the global genus, while
THM-4292/4297 make every possible positive-genus exceptional tail over the
contact Keller-constant.

## 6. Positive good-form orders and extinction

THM-4327's face-scaling formula gives, for a supporting plane
`h=a_s r+b l+c`,

```text
ord_face(phi^*eta_0)=L(5/6-a_s-b-c).                       (29)
```

Equations `(6)` and `(8)` yield

```text
ord_M(phi^*eta_0)=27,                 ord_N(phi^*eta_0)=34. (30)
```

Equations `(14)--(16)`, the primitive height normals, the node charts in
Section 4, and Section 5's sole repeated-contact audit verify the
endpoint hypotheses of THM-4327's proper-flat face/edge interface.  By
`(15)`, the relative differential generator does not vanish generically on
either positive-genus component.  Their maps to the characteristic-zero
good elliptic curve are therefore constant.  The component `R`, every toric
chain, and every ordinary resolution component are rational and hence also
constant; Section 5 handles every exceptional tail.

After finite base change and regularization, let `M_phi` be the actual
pullback of a positive-degree line bundle on the good target.  Proper-flat
intersection, retaining every fibre multiplicity, gives

```text
deg(M_phi|generic fibre)=sum_i m_i deg(M_phi|X_i)=0.        (31)
```

A nonautomorphic Keller pair has positive generic response degree, contrary
to `(31)`.  This proves the theorem.

## 7. Sharp boundary, connections, and next tasks

The mechanism first fails at `U=0`, where `(22)` loses its Weierstrass
quadratic coefficient and `(21)` loses its transverse derivative.  At `W=0`
or `zeta_3=0`, the lower complex and active edge boundary change.  In
particular this theorem does not cover the other `Z=0` owner stratum
`zeta_3=0`, any `U=0` intersection, either repeated `U=0` cusp, seam entry,
`JC(2)`, or `DC(2)`.

The first generated task--recomputing the `zeta_3=0` lower hull before edge
specialization--is completed by THM-4337.  Two tasks survive the
concept-board comparison:

1. for the open repeated `(2,5)` and `(2,3)` cusps, recenter at the unique
   critical section and compute the invariant critical-value series, stable
   tails, and conductor/Kahler sidecar;
2. retain raw odd source-normal flags
   `beta_11=h_1-5h_0` and `zeta_3=h_1-4h_0`; THM-4328's scalar Student
   observer kills those odd faces and cannot distinguish the owner walls.

These are proposals, not consequences of the theorem.

## 8. Reproduction and honest scope

From the repository root, run the primary and clean-room exact certificates
in normal and optimized modes.  The paired streams must byte-match their
frozen outputs.

```bash
python3 -B 04-computation/jc2_m12_z0_beta0_endpoint_extinction_thm4334.py
python3 -B -O 04-computation/jc2_m12_z0_beta0_endpoint_extinction_thm4334.py
python3 -B 04-computation/jc2_m12_z0_beta0_endpoint_extinction_independent_audit_thm4334.py
python3 -B -O 04-computation/jc2_m12_z0_beta0_endpoint_extinction_independent_audit_thm4334.py
```

What is proved here is only the relative extinction of `(1)`.  The source
infinity packet is a consistency check off contact and is never transported
through `Lambda=0`.  THM-4337 closes the complementary `zeta_3=0` owner wall
away from `beta_11=0`; the residual owner intersections, exact-seam entry,
`JC(2)`, and `DC(2)` remain open. **QED.**
