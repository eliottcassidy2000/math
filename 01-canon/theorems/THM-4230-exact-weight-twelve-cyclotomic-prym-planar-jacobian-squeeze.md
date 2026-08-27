---
id: THM-4230
title: "Exact-weight-twelve cyclotomic Prym planar Jacobian squeeze"
status: >
  PROVED RELATIVE TO THM-3992/3997/4012/4045/4103/4120/4122/4147/4218
  + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the inherited b=d=0 reduced
  (2,3) seam, the complete exact-M=12 gate UZ(W^2-4UZ)(U+W+Z)!=0 has one
  genus-seven positive-genus component and response degrees 42/34. Outside
  a countable proper primitive-Prym hidden-j=0 Hom locus, the specialized
  degree is divisible by four, so no nonautomorphic planar Keller pair exists
  there. The hidden Hom locus is nonempty at W=0, where an explicit rank-two
  Eisenstein sublattice has Gram determinant 24; degree 34 is absent and all
  192 degree-42 vectors fail attachment collapse inside this sublattice.
  Its saturation in the full hidden Hom lattice, exclusion on the locus, the
  whole M=12 gate, seam entry, JC(2), and DC(2) are OPEN.
source: codex-jc-lrc-niche-crossfeed-20260826
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
related:
  - THM-4222-dense-weight-eleven-primitive-cm-planar-jacobian-exclusion
  - THM-4226-dense-weight-thirteen-primitive-cm-bolza-planar-jacobian-exclusion
  - THM-4232-weight-eleven-u-zero-primitive-cm-planar-jacobian-exclusion
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14, through THM-4045, supplies
  the face/edge model and rational toric chains. Standard characteristic-zero
  inputs below are Poincare complete reducibility, connectedness of the norm
  kernel for a ramified double cover, the rank-one CM criterion and induction
  of imprimitive CM types, finiteness of CM types for a fixed CM field,
  finiteness of automorphism groups in genus greater than one, countability
  of principally polarized members of a fixed complex isogeny class, and
  Torelli. Milne's "Complex Multiplication" notes (2020), Sections 1.9--1.10,
  supply the CM-type input. None supplies the exact in-repo arithmetic.
script: 04-computation/jc23_exact_weight12_cyclotomic_prym_squeeze_thm4230.py
output: 05-knowledge/results/jc23_exact_weight12_cyclotomic_prym_squeeze_thm4230.out
script_sha256: e1d84a4399473f5a39f35b01892f0344d83e22db3c270f640ca5dcd621d56aaa
output_sha256: 776542de91432f0217f471d8e914d0fa9fefee4160c4167b6eeac542e576a743
independent_audit_script: 04-computation/jc23_exact_weight12_cyclotomic_prym_squeeze_thm4230_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_exact_weight12_cyclotomic_prym_squeeze_thm4230_independent_audit.out
independent_audit_script_sha256: 4647d5dc58e35e53bfec5aa2e9996384786c47af0e89803fe0a9ee39347ca49b
independent_audit_output_sha256: f2f377ee8cec8424754e3c31b8e857aad2ada8e2f8c9ef96e106e6b754e2cfda
visible_saturation_script: 04-computation/jc23_exact_weight12_visible_e0_saturation_thm4230.py
visible_saturation_output: 05-knowledge/results/jc23_exact_weight12_visible_e0_saturation_thm4230.out
visible_saturation_script_sha256: bc205c520a8830f1a77e14d2fc6fee8c716f15d8ad2878690f43e7eb7320b507
visible_saturation_output_sha256: 39c88771b90bf22a094af3eaaf611e9c5ca0d0369864a42b5a78bc341c33c073
explicit_w0_map_script: 04-computation/jc23_weight12_w0_hidden_e0_explicit_map_thm4230.py
explicit_w0_map_output: 05-knowledge/results/jc23_weight12_w0_hidden_e0_explicit_map_thm4230.out
explicit_w0_map_independent_audit_script: 04-computation/jc23_weight12_w0_hidden_e0_explicit_map_thm4230_independent_audit.py
explicit_w0_map_independent_audit_output: 05-knowledge/results/jc23_weight12_w0_hidden_e0_explicit_map_thm4230_independent_audit.out
explicit_w0_map_script_sha256: dd8556b045ac3e9f2ad240c2a9e0bc3ad193105d559b0f922a0f3ec5283ed32d
explicit_w0_map_output_sha256: 05cda15ec6e7b9c384899f649d7d2444ccd9b192040b34b3212191582c849e23
explicit_w0_map_independent_audit_script_sha256: e648e25d9cc8860fcab911dded2eb12d8eca8b19de7ff4d133f36a9f1193ae3a
explicit_w0_map_independent_audit_output_sha256: 3853c7cd9971065c13e981324793f72fa70c677440e6f13c621fc577059cdc2c
explicit_w0_lattice_script: 04-computation/jc23_thm4230_w0_hidden_e0_degree42_lattice.py
explicit_w0_lattice_output: 05-knowledge/results/jc23_thm4230_w0_hidden_e0_degree42_lattice.out
explicit_w0_lattice_script_sha256: e33c8c8edba33e132abcfd8691be5df6d2fc918328d911275fe05642cfa82935
explicit_w0_lattice_output_sha256: 5d6bcc4fc2ec94e52a7a528c0593401629d4543bd875cd66569beede3c67cedc
hash_basis: raw LF bytes
audit: >
  PASS. The primary SymPy certificate executes 16,548 exact checks over all
  16,384 lower-row subsets and 128 coincident-support cancellations. A
  standard-library clean-room audit independently reconstructs the support,
  including the fixed residual point omitted by the primary support helper,
  hull, genera, character packet, response degrees, carrier wall, and both
  hostiles. A separate four-element F_4 certificate proves integral saturation
  of the visible j=0 Hom lattice. Two further exact paths verify the explicit
  W=0 degree-six hidden-E_0 map: SymPy factorization and a standard-library
  quotient-ring reconstruction. A further exact certificate computes its
  rank-two Eisenstein Gram, exhausts the 192 degree-42 vectors, checks 16
  attachment-orbit representatives in four good finite-field embeddings,
  and proves the rank-one norm-seven obstruction by a nonzero resultant.
  Normal, optimized, and fixed-hash-seed runs byte-match the frozen outputs.
---

# THM-4230 -- exact-weight-twelve cyclotomic Prym squeeze

**PROVED RELATIVE TO THM-3992/3997/4012/4045/4103/4120/4122/4147/4218
+ VERIFIED-EXACT + INDEPENDENTLY AUDITED. THIS EXCLUDES ONLY THE COMPLEMENT
OF A COUNTABLE PROPER NONEMPTY HOM LOCUS; `M=12` AND `JC(2)` REMAIN OPEN.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam. Put

```text
U=[p^6]H,       W=[p^3y^2]H,       Z=[y^4]H,
D=W^2-4UZ,      Lambda=U+W+Z,       kappa=W^2/(UZ).     (1)
```

Let `A_12(kappa)` be the primitive `Phi_12` fourfold constructed below and

```text
H_0={kappa: Hom(A_12(kappa),E_0)!=0},
E_0:Y^2=X^3+1.                                         (2)
```

> **Theorem.** On the complete exact-weight-twelve gate
>
> ```text
> U*Z*D*Lambda != 0,                                   (3)
> ```
>
> every hypothetical nonautomorphic planar Keller pair has
> `kappa in H_0`. The set `H_0` is countable and proper, but nonempty:
> `0 in H_0`, witnessed by the gate-interior locus `W=0`. Hence the
> complement of `H_0` contains no such pair.

This proves neither exclusion on `H_0` nor exclusion of the whole gate.
Walls, other cells, seam entry, `JC(2)`, and `DC(2)` remain open.

The closest proved mechanism is THM-4222/4226's complete lower-Newton and
cyclotomic analysis. THM-4218's hidden `j=0` tail is the canonical hostile.
The corrected near miss is that a rational character decomposition neither
excludes hidden factors nor saturates an integral Hom lattice. The decisive
sidecars are the `Phi_12` Prym, graph gluing on `E_0[2]`, and the twelve
labelled attachments.

## 2. Exact source, lower hull, and regular model

Use `s=XT`, `p=T+s^2`, `y=sp`, and `t=p-s^2=T`. The complete source is

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y+zeta_3 y^3
  +upsilon_5 p^5+xi_10 p^2y^2+alpha_11 p^4y+beta_11 py^3
  +U p^6+W p^3y^2+Z y^4,
K=2848/45-(7/6)Delta.                                  (4)
```

All displayed lower coefficients are arbitrary. Enumerating
`0<2i+3j<=12` and deleting only the forbidden `y,py` rows gives sixteen
monomials. For

```text
F_Q=(s^2-p)(1-QH)-Q s^2/2,                             (5)
```

the unique lower plane and equality set are

```text
nu(r,l)=(r+2l-2)/12,
(0,1,0),(2,0,0),(0,7,1),(2,6,1),(4,5,1),(6,4,1).      (6)
```

The endpoint coefficients `U,Z` own the hull vertices, so `(6)` survives all
`2^14` optional-row deletions and `2^7` coincident-support cancellations.
In particular, cancellations at `W-U=0` or `Z-W=0` require no extra gate.
The primary certificate's `expanded_support` helper omits the fixed residual
point `(2,0,1)` contributed by `-Qs^2/2`; the independent audit includes it
and verifies its gap from the unique lower plane is exactly `1`. The primary
scaled-model identity `(9)` also retains the term. Thus the omission changes
neither the hull nor any theorem conclusion, but the independent path is the
load-bearing full-support check for this point.

The global polygon, Pick ledger, and face factorization are

```text
P_global=conv((0,1),(2,0),(6,4),(0,7)),
(2Area,boundary,interior)=(48,14,18),
g_M=(S^2-P)(1-U P^6-W S^2P^5-Z S^4P^4)=R*C.          (7)
```

The four edge schemes are

```text
X-1,       1-ZX^4,       (X-1)(UX^2+WX+Z),       U-X^6. (8)
```

Their nonconstant discriminant/resultant factors are exactly `Z,D,Lambda,U`.

Take `Q=sigma^12`, `s=sigma^-1S`, and `p=sigma^-2P`. The exact integral
model is

```text
(S^2-P)(1-H_sigma)-sigma^12 S^2/2=0.                  (9)
```

The Newton triangle of `C` is

```text
conv((0,0),(0,6),(4,4)),
(2Area,boundary,interior)=(24,12,7),                   (10)
```

so `C` has genus seven. On `R`, the other factor is `1-Lambda S^12`.
Therefore `R,C` have twelve distinct contacts

```text
P=S^2,                    Lambda S^12=1,               (11)
```

with gradient determinant `-12Lambda S^11`. Locally `(9)` is
`uv=sigma^12/2`, giving an `A_11` rational path. The main-face normal has
last coordinate one, so the face components have multiplicity one. After
contracting rational subdivisions, two vertices are joined by twelve paths:

```text
b_1=12-2+1=11,             g_special=7+11=18.          (12)
```

This matches global Pick and proves that `C` is the sole positive-genus
component.

## 3. Packet and composite quartic carrier

The three non-affine edges have `(length,index)` equal to `(1,1),(4,2),
(3,11)`. The vertical `s=0` edge is affine and omitted. Hence

```text
packet=(11,11,11,2,2,2,2,1),             sum=42.      (13)
```

The length-four edge is one irreducible separable degree-four closed point
over `K_0=C(q)`:

```text
q=1/2+K w^2+zeta_3 w^3+Z w^4.                         (14)
```

The polynomial map has degree four, proving `[C(w):C(q)]=4`; characteristic
zero gives separability. This extension can have a quadratic intermediate.
Depressing the quartic gives the exact composition wall

```text
zeta_3*(zeta_3^2-4KZ)=0.                              (15)
```

For example `K=Z=1,zeta_3=2` gives
`q=1/2+(w(w+1))^2`. Thus a finite carrier is not uniformly birational of
degree four.

THM-4120 gives `E_q(K_0)={O}`. Conjugacy makes all four index-two punctures
map simultaneously either to `O` or away from it. The inherited response
interfaces therefore give exactly two **origin-fibre degrees**:

```text
full: n=42,                    finite: n=42-4*2=34.     (16)
```

A quadratic intermediate changes finite-carrier residue/meridian data but
not `(16)`.

## 4. Cyclotomic and elliptic quotients

The automorphism

```text
tau:(S,P)->(xi_12 S,xi_12^2P)                          (17)
```

has holomorphic character multiset

```text
{5,7,8,9,10,11,11}.                                   (18)
```

It comes from the seven interior points
`(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(3,4)`, where the
residue differential at `(r,l)` has character `r+2l mod 12`. With conjugates,

```text
J(C) ~ A_12 x E_(Phi6) x E_(Phi4) x E_(Phi3),
dim A_12=4,       j(E_(Phi6))=j(E_(Phi3))=0,
j(E_(Phi4))=1728.                                      (19)
```

Quotient by `tau^6:S->-S` and put

```text
u=1/P,                   x=S^2/P,                v=W+2Zx.
```

The genus-two quotient and its elliptic maps are

```text
B: v^2=D+4Zu^6,
B->E_a: (X,Y)=(u^2,v),        Y^2=D+4ZX^3,
B->E_b: (X,Y)=(u^-2,v/u^3),   Y^2=4Z+DX^3.             (20)
```

Each map `C->E_a,E_b` has degree four and `j(E_a)=j(E_b)=0`. The remaining
elliptic quotient is

```text
T=SP,       V=2UP^3+WT^2,       V^2=4U+DT^4,          (21)
```

of degree three and `j=1728`, hence Hom-orthogonal to `E_0`.

## 5. Saturated visible Hom lattice and connected Prym

Over `C`, scale `(20)` to

```text
H_2:y^2=x^6+1,                    E_0:Y^2=X^3+1,
f_1(x,y)=(x^2,y),                 f_2(x,y)=(x^-2,y/x^3). (22)
```

Let `O=Z[xi_3]=End(E_0)`. The map

```text
psi:E_0^2->J(H_2),       (P,Q)|->f_1^*P+f_2^*Q         (23)
```

is a `(2,2)`-isogeny. Indeed `a:(x,y)->(-x,y)` fixes `f_1` and negates
`f_2`, so the cross Hom terms vanish; the pulled-back principal polarization
is diagonal of degree two on each factor and `deg psi=4`.

For `T_r=(r,0)`, `r^3=-1`, the divisors `f_1^*(T_r-O)` and
`f_2^*(T_(r^-1)-O)` differ by

```text
(0,1)+(0,-1)-infinity_+-infinity_-=div(x).             (24)
```

Thus `ker psi` is the graph of `iota(T_r)=T_(r^-1)`. Identify
`E_0[2]=O/2O=F_4`. Then `iota(t)=t^2`, Frobenius rather than an `O/2O` scalar.
A row `(a,b) in O^2` descends precisely when

```text
a t+b t^2=0                    for all t in F_4.        (25)
```

At `t=1` this gives `a+b=0 mod 2`; at a primitive element it gives
`a=b=0 mod 2`. Conversely `2O^2` kills the graph. Since the two norm maps
pull back through `psi` to `([2],0)` and `(0,[2])`, respectively,

```text
Hom(J(H_2),E_0)=O f_(1*) direct-sum O f_(2*),
deg_H2(alpha f_1+beta f_2)=2(N(alpha)+N(beta)).         (26)
```

The cover `C->B` has total ramification eight by Riemann--Hurwitz. Its norm
kernel is connected and equals the Prym. A homomorphism vanishing on that
Prym therefore factors integrally through `J(C)->J(B)`, not merely after
multiplication by two. Composing `(26)` with the degree-two cover gives

```text
deg_C=4(N(alpha)+N(beta)).                              (27)
```

## 6. Countable proper nonempty hidden-Hom locus

The Prym decomposition is

```text
Prym(C/B) ~ A_12 x E_1728.                             (28)
```

Let `L=Q(xi_12)=Q[tau]`. It acts faithfully on `A_12`, with

```text
H^1(A_12,Q)=L^2,
(dim H_k^(1,0))_(k=1,5,7,11)=(0,1,1,2).               (29)
```

Suppose `Hom(A_12,E_0)!=0`. The maximal `E_0`-isotypic subvariety is
canonical up to isogeny and therefore `L`-stable. The field `L` acts
faithfully on its nonzero homology. If it is isogenous to `E_0^m`, then
`[L:Q]=4` divides `2m`, so `m=2` or `4`. For `m=4`, `A_12` lies in the fixed
isogeny class of `E_0^4`. For `m=2`, its complementary surface is rank one
over `L`, hence is a CM surface with one of finitely many CM types of `L`.
After adjoining the three fixed elliptic factors in `(19)`, every hidden-Hom
Jacobian lies in a finite union of fixed complex isogeny classes.

Write the roots of `U+Wx+Zx^2` as `alpha,beta`. The cyclic presentation

```text
S^12=x^6/(U+Wx+Zx^2)                                  (30)
```

has branch residues `(6,11,11,8)` at `(0,alpha,beta,infinity)`. The unique
ramification-index-two and index-three points distinguish zero and infinity;
the two index-twelve points are unordered. Their marked invariant is

```text
alpha/beta+beta/alpha=kappa-2.                          (31)
```

A fixed genus-seven curve has finite automorphism group and only finitely
many marked order-twelve actions, so the map from the open `kappa`-curve to
curve moduli is nonconstant. A nonconstant morphism of algebraic curves has
finite fibres. A fixed complex isogeny class has only countably many
principally polarized members, and Torelli is injective. Therefore `H_0` is
countable. It is proper because the complex parameter curve is uncountable.

It is not empty. At `W=0`, with `UZ(U+Z)!=0`, the extra involution

```text
rho:(S,P)->(S,-P)                                      (32)
```

commutes with `tau` and acts on the differential indexed by `(r,l)` as
`(-1)^l`. Its primitive eigenspaces have CM types

```text
rho=+1: {5,11},                  rho=-1: {7,11}.       (33)
```

The first is stabilized by `{1,7}` in `(Z/12Z)^*`, with fixed field
`Q(xi_3)`; the second is stabilized by `{1,5}`, with fixed field `Q(i)`.
The imprimitive-CM decomposition yields

```text
A_12(W=0) ~ E_0^2 x E_1728^2.                          (34)
```

Thus `0 in H_0`. Equations `(29)--(34)` prove that the named locus is
countable, proper, and nonempty; they do not classify it.

### 6.1. An explicit hidden degree-six quotient at `W=0`

The abstract imprimitive-CM witness in `(34)` can be made explicit. Scale the
`W=0` positive-genus component to

```text
C_0: x^6+y^4=1,
x=U^(1/6)P,                  y=Z^(1/4)SP.              (34a)
```

Put

```text
t=(1+y^2)/x^3,
x^3=2t/(t^2+1),             y^2=(t^2-1)/(t^2+1).      (34b)
```

Let `a` be any root of

```text
a^4-2a^3-2a+1=0,                                      (34c)
```

and choose `s` with

```text
s^6=4/(2a^3+3a^2-1).                                  (34d)
```

The denominator in `(34d)` is nonzero: its resultant with `(34c)` is
`-108`. Define

```text
f_a:C_0->E_0,
X=(s^2/2)x(t^2-a^2)/t,
Y=(s^3/2)y(t^2+a^3)/t.                                (34e)
```

Writing `m=t^2`, the load-bearing identity is

```text
(m-1)(m+a^3)^2-(m-a^2)^3
 =(2a^3+3a^2-1)m(m+1),                                (34f)
```

because the difference between the two sides factors as

```text
m(a+1)^2(a^4-2a^3-2a+1).                              (34g)
```

Substitution of `(34b)` and `(34d)` into `(34e)` now gives exactly
`Y^2-X^3=1`.

This map has degree six. Over `t=0` there are two points, each with
`ord(t)=3`, `ord(x)=1`, and an order-two pole of `X`; the same holds over
infinity. Over each of `t=i,-i` there is one point with `ord(t-i)=6`,
`ord(x)=-2`, again giving an order-two pole of `X`. There are no others.
The resultant of `(34c)` with `a^2+1` is `4`, so no listed pole cancels.
Thus the pole divisor of `X` has degree

```text
2*2+2*2+1*2+1*2=12=2 deg(f_a),
deg(f_a)=6.                                             (34h)
```

Let

```text
iota:(x,y)->(x,-y)=tau^6,       sigma:(x,y)->(xi_3x,y).
```

Then

```text
f_a o iota=[-1] o f_a,          f_a o sigma=[xi_3] o f_a. (34i)
```

The first identity places `f_a` in the Prym. The `E_1728` factor in `(28)`
is Hom-orthogonal to `E_0`, so `f_a` is an explicit hidden quotient of
`A_12`, rather than one of the visible maps through `B`.

The same formula also reaches the actual attachment sidecar. On the twelve
contacts `P=S^2`,

```text
c=y^2/x^3=sqrt(Z/U),
c=(t-t^-1)/2.                                          (34j)
```

Thus they split into the two six-point fibres at `t` and `-1/t`, exchanged
by `tau`; within either fibre, `sigma` and `iota` generate the three CM
rotates and their negatives. Here `t=+-1` would force `Z=0`, so the gate makes
both fibres genuine and disjoint. If `f_a` took one common value on a six-point
fibre, that value would lie in both `E_0[3]` and `E_0[2]`, hence would be the
origin. But `(34e)` is affine there: `t=0,infinity` would force a vanished
gate endpoint, while `t=+-i` is exactly `U+Z=0`, excluded by `(3)`. Therefore

```text
f_a alone does not collapse the twelve attachments.                   (34k)
```

For `alpha in O=Z[xi_3]`, postcomposition has

```text
deg(alpha f_a)=6N(alpha).                               (34l)
```

Consequently this explicit rank-one line contains no degree-34 map. It does
contain degree-42 candidates arithmetically because `N(3+xi_3)=7`. Those
candidates can also be excluded from attachment collapse, and the second
explicit hidden direction can be retained.

Let `omega=xi_3`, let `T` act on maps by precomposition with `tau`, and choose
one root `a` in either reciprocal pair of `(34c)`. After harmless `O`-unit
choices, put

```text
f=f_a,                         g=Tf,
T^2=-omega.                                             (34m)
```

The four roots of `(34c)` give, up to `O`-units,

```text
f,             g,             f+g,             g-omega f. (34n)
```

An exact elliptic addition calculation finds

```text
deg(f)=6,             deg(f+g)=6,             deg(f-g)=18.
```

The cancellation in `f+g` is load-bearing: after imposing
`a^2-(1+sqrt(3))a+1=0` and `xi_12^4-xi_12^2+1=0`, the determinant of the two
linear `t^2` coefficients in the addition slope is exactly zero; for `f-g`
it is nonzero and creates twelve additional double `X`-poles. Hence the
`O`-lattice

```text
L_exp=O f direct-sum O g                              (34o)
```

has Hermitian degree Gram matrix, with the form linear in its first slot,

```text
             [ 6,       -4-2omega   ]
H_exp=       [ -4-2omega^2,    6   ],
det(H_exp)=24.                                         (34p)
```

In particular, every degree in `L_exp` is divisible by six, so degree `34`
does not occur. The positive-eigenvalue bound for `(34p)` reduces degree
`42` to Eisenstein coefficients in `[-6,6]^4`. Exact enumeration gives

```text
192 degree-42 vectors,
16 orbits under scalar mu_6 and T,
determinant norms 13,49,61,73, each with multiplicity 48. (34q)
```

For a representative `h=af+bg`, attachment collapse is equivalent to
`h(A)=Th(A)=O`. Write `X_h=xA_h(t)` and let `d_h(t)` be the reduced
denominator of `A_h` over `K=Q(a,s,omega)`. This would force a common root of

```text
d_h(t)             and             t^deg(d_h)d_h(-1/t). (34r)
```

There is a characteristic-zero degree bound missing from a bare modular
test. Both coefficient pairs `(X/x,Y/y)` for `f,g` lie in
`t^-1 K[t^2]`, and `c=y^2/x^3=(t^2-1)/(2t)` is odd. In the coefficient-form
addition law

```text
A_(P+Q)=c((B_Q-B_P)/(A_Q-A_P))^2-A_P-A_Q,
B_(P+Q)=((B_Q-B_P)/(A_Q-A_P))(A_P-A_(P+Q))-B_P,
```

odd pairs map to odd pairs; doubling, negation, and multiplication by
`omega` preserve this. Hence `A_h(-t)=-A_h(t)`. On the smooth projective
normalization of `C_0`, both points over `t=0` map to `O` under `f` and `g`,
so every nonzero degree-42 `h` also maps them to `O`. Thus `A_h` has a pole at
zero. Coprimality of its numerator and denominator now forces `d_h` itself to
be odd, so `deg d_h` is odd.

If a finite root of `d_h` has multiplicity `e`, it contributes at least `6e`
to the pole divisor of `X_h=xA_h`, except that `t=0` contributes `6e-2`;
at `t=+-i` the contribution is `6e+2`. Indeed `t=+-1` has three points of
ramification two, `t=+-i` has one point of ramification six and
`ord(x)=-2`, and every other fibre has six unramified points. Therefore the
finite poles contribute at least

```text
6 deg(d_h)-2 <= deg((X_h)_infinity)=2 deg(h)=84.
```

So `deg d_h<=14`, and oddness sharpens this to `deg d_h<=13`.

In each of four exact good reductions

```text
(q,xi_12,a,s)=(313,29,135,21), (349,24,246,28),
               (373,69,297,33), (397,157,161,27),     (34s)
```

the reduced denominator has degree `13` and the gcd in `(34r)` has degree
zero. The embeddings satisfy the cyclotomic, quadratic-root, and sixth-root
equations with nonzero derivatives. Specialization cannot increase reduced
denominator degree, so the characteristic-zero bound forces `deg d_h=13` and
proves degree preservation. A characteristic-zero common factor would make
the resultant zero in any good specialization where these coefficients are
defined, contradicting the displayed gcd-one reductions. One prime suffices
logically; four are hostile controls. This is a finite-exact exclusion of
attachment collapse for every degree-42 vector in `L_exp`.

On the rank-one line there is also a symbolic control. The two norm-seven
CM kernels have `X`-coordinates satisfying

```text
7X^6-4X^3+16=0.                                      (34t)
```

Writing `r=X^3` and `m=t^2`, the two attachment fibres would require

```text
(m-a^2)^3-r(2a^3+3a^2-1)m(m+1)=0,
-(a^2m-1)^3-r(2a^3+3a^2-1)m(m+1)=0.                  (34u)
```

Successively eliminating `m`, the quartic root `a`, and
`7r^2-4r+16` gives the nonzero total resultant

```text
2^72 * 3^54.                                          (34v)
```

Thus degree `42` cannot collapse the attachments even on any single
`O f_a` line. Equations `(34p)--(34v)` exclude degrees `34` and `42` only in
the explicit rank-two sublattice `L_exp`. They do **not** prove that
`L_exp=Hom(A_12,E_0)`, nor do they settle the integral gluing with the visible
lattice. Saturation of the hidden lattice and the full attachment degree form
are the precise surviving `W=0` problems.

## 7. Relative exclusion and hostiles

Take `kappa notin H_0`. Equations `(28)`, the `j=1728` mismatch, and the
definition of `H_0` give `Hom(Prym(C/B),E_0)=0`. Connectedness of the norm
kernel and `(26)--(27)` make every nonconstant specialized component-map
degree divisible by four. Every other component is rational and constant.
Proper-flat degree conservation gives

```text
4 | deg(phi_generic).                                  (35)
```

But the complete responses `(16)` are `42=2 mod 4` and `34=2 mod 4`, a
contradiction. This proves exclusion only on the complement of `H_0`.

Two gate-interior hostiles fix the boundary. First, `W=0` gives `(34)`, so
the visible factors do not exhaust `Hom(J(C),E_0)` everywhere. Second, at

```text
(U,W,Z)=(2,-2,1),                 (D,Lambda)=(-4,1),   (36)
```

all attachments have `x=1,u^6=Lambda,v=W+2Z=0`. Under both maps `(20)` they
are the three nonzero two-torsion points. Composing either degree-four map
with `[2]` gives a degree-sixteen nonconstant map sending all twelve
attachments to the origin. Attachment equality therefore does not imply
constancy. Since `16` is neither `34` nor `42`, this is not a Keller map.

The exact strongest survivor is

```text
hypothetical exact-M12 Keller parameter => kappa in H_0,
H_0 countable, proper, nonempty.                        (37)
```

Classification of `H_0`, its node-annihilating degree lattice, the walls,
all-gate `M=12`, seam entry, `JC(2)`, and `DC(2)` remain **OPEN**.

## 8. Verification and firewall

Replay with

```bash
python3 -B 04-computation/jc23_exact_weight12_cyclotomic_prym_squeeze_thm4230.py
python3 -B 04-computation/jc23_exact_weight12_cyclotomic_prym_squeeze_thm4230_independent_audit.py
python3 -B 04-computation/jc23_exact_weight12_visible_e0_saturation_thm4230.py
python3 -B 04-computation/jc23_weight12_w0_hidden_e0_explicit_map_thm4230.py
python3 -B 04-computation/jc23_weight12_w0_hidden_e0_explicit_map_thm4230_independent_audit.py
python3 -B 04-computation/jc23_thm4230_w0_hidden_e0_degree42_lattice.py
```

The primary certificate covers every variable support row, lower-row subset
and collision, hull, edge schemes, Pick/packet ledgers, scaled model, nodes,
characters, elliptic quotient identities, and `W=-2Z`; its support helper has
the fixed-residual limitation recorded in Section 2. The independent
standard-library audit reconstructs the full support, geometry, carrier wall,
responses, characters, `W=0`, `W=-2Z`, and the `F_4` residue result without
importing the primary implementation. The dedicated saturation certificate
exhausts all sixteen coefficient pairs on all four elements of `F_4`. The two
explicit-map paths respectively use exact SymPy elimination and a hand-written
quotient ring. The hidden-lattice certificate checks the characteristic-zero
addition/Gram identities and pole-budget ledger, the complete bounded
degree-42 census, four good finite-field attachment resultants, and the
independent symbolic norm-seven resultant. Its finite-field statement is
**FINITE-EXACT**; the degree bound, map, and rank-one identities are exact.

The computations do not prove the stated standard CM/Prym/Torelli inputs,
classify `H_0`, handle walls, prove entry, or prove `JC(2)`. These are explicit
theorem boundaries. **QED.**
