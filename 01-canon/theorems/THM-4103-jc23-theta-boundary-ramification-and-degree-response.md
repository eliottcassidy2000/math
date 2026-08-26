---
id: THM-4103
title: "JC(2) theta-only boundary ramification and degree response"
status: >
  PROVED RELATIVE TO THM-3992/4053 + VERIFIED-EXACT + INDEPENDENTLY
  VERIFIED-EXACT on the smooth nonresonant theta-only maximum-weight-eight
  survivor; JC(2) remains OPEN. The complete generic Newton boundary forces
  seven source-infinity punctures with ramification indices
  {1,2,2,3,3,3,7}. They contribute all fourteen Riemann--Hurwitz units.
  The index-one puncture is forced to target infinity, while the two
  index-two punctures form one quadratic Galois orbit and respond together.
  The generic Keller fibre degree is therefore one of {7,12,21}, with only
  3 coarse, 4 edge-refined, and 5 fully labelled target-infinity responses.
  The index-seven branch also forces deg_t(A),deg_t(C)>=7 and total degrees
  at least 15,16. These constraints do not themselves exclude the survivor;
  THM-4130 later does so by critical monodromy.
source: codex-arithmetic-boundary-breakthrough-20260825
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
related:
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
script: 04-computation/jc23_theta_boundary_response_thm4103.py
output: 05-knowledge/results/jc23_theta_boundary_response_thm4103.out
independent_audit_script: 04-computation/jc23_theta_boundary_response_thm4103_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_theta_boundary_response_thm4103_independent_audit.out
script_sha256: 2ad7fe9eb5ab3a683c4c3ea552e4c88f989cc654133e91df014340ed8398a71c
output_sha256: dd387a9865cb4a0a397e510b8e7c6a64750b1bf9c64c8e2907dc2d6b24cc29af
independent_audit_script_sha256: 682b14ecfbd932dc0d6f15b2b0fd9c42541eba2b4a509f03e2b0f747f0a1d3c7
independent_audit_output_sha256: e3fb544a4d0e28c8b827889efe1670e5bf7fae0b2fd14b8bf4749e917bc36dfb
semantic_sha256: d1c19dc8f5a787dfec43297a86de061153b3ca543b6ab87fae4f4db1f8497c14
independent_semantic_sha256: ae508cbce4205834f0e9504ffba43430836a5954d62e6d559e31ebd29ea412d0
hash_basis: raw LF bytes for files; Python repr/JSON for semantic ledgers
audit: >
  PASS. The primary checks all eight lambda/alpha/phi presence patterns,
  including six counterfactual lambda/alpha-deletion hostiles, by convex hull, then verifies edge distances, Pick genus,
  ramification, the Eisenstein intersection, literal and grouped response
  profiles, and a support-deletion hostile. The independent standard-library
  path expands the polynomial term by term, uses primitive supporting
  inequalities instead of a hull routine, includes the permitted phi=0
  branch, and cross-checks norms, raw and Galois-refined profiles, and the
  index-seven support floors by independent paths. Normal and optimized streams byte-match both
  frozen outputs.
---

# THM-4103 -- the theta-only boundary is a complete ramification response

**PROVED RELATIVE TO THM-3992/4053 + VERIFIED-EXACT + INDEPENDENTLY
VERIFIED-EXACT on one live reduced seam; JC(2) OPEN.** Work over an
algebraically closed field `k` of characteristic zero, put
`Kbar=overline{k(q)}`, and understand every generic point and puncture below
geometrically over `Kbar`. Assume a hypothetical
Keller pair `(A,C)` lies in THM-4053's smooth nonresonant theta-only survivor

```text
delta=0,                theta!=0,                Delta_V!=0.       (1)
```

This theorem uses the entire generic Newton boundary, rather than only its
genus or its positive-genus component. It localizes the full ramification
divisor and turns the surviving degree question into a finite response
universe. It does not prove that any response is realizable.

## 1. Generic source, target, and Keller differential

Retain THM-3992's coordinates and centered target polynomial

```text
s=xt,                 p=t+s^2,                 y=sp,
E(U,V)=V^2-U^3+(3a^2/4)U+a^3/4.                         (2)
```

On the theta-only seam,

```text
E(A(x,t),C(x,t))=gamma*s^2/t+H(p,sp),
H=lambda*p+alpha*p^2+epsilon*p^3+kappa*s^2*p^2
  +phi*s*p^3+theta*s^2*p^3.                             (3)
```

The inherited rows force

```text
lambda=3a/(2gamma)!=0,       alpha=8/(3a^7)!=0,
epsilon*kappa*theta*(kappa-epsilon)!=0.                 (3a)
```

Only `phi` may be zero. Let `C_q` be the smooth projective completion
over `Kbar` of the source fibre `E(A(x,t),C(x,t))=q`, let `E_q` be the
projective completion of the target curve `{(U,V):E(U,V)=q}`, and let

```text
varphi_q:C_q -> E_q                                               (4)
```

be the extension of the polynomial map `(A,C)`. Put `Q=q^-1`. Multiplication
by `Qt`, together with `ds wedge dp=t dx wedge dt`, gives the exact source
equation

```text
F_Q(s,p)=(s^2-p)(1-QH(p,sp))+gamma*Q*s^2
        =Q*t*(E(A(x,t),C(x,t))-q).                       (5)
```

> **Lemma 1.1 (Keller residue identity).** With compatible residue
> orientation,
>
> ```text
> varphi_q^*(dU/(2V))=dA/(2C)
>   =Res(dx wedge dt/(E(A,C)-q))
>   =Q Res(ds wedge dp/F_Q)
>   =Q ds/(F_Q)_p.                                      (6)
> ```

### Proof

On the source fibre take the tangent vector

```text
v=(E_t,-E_x).                                           (7)
```

The chain rule and the Keller identity `J_(x,t)(A,C)=1` give

```text
dA(v)=A_x E_t-A_t E_x
     =E_C(A_x C_t-A_t C_x)=2C.                         (8)
```

Thus the first two differentials in `(6)` agree on the generic tangent. The
coordinate Jacobian `ds wedge dp=t dx wedge dt` and `(5)` give the middle
equality. The last equality is the usual affine residue chart. QED.

The differential `dU/(2V)` is nowhere zero on the smooth projective elliptic
curve `E_q`. Hence

```text
div(varphi_q^*(dU/(2V)))=div(dA/(2C))=R_varphi,         (9)
```

the ramification divisor of `varphi_q`.

## 2. Complete Newton polygon and edge orders

After coincident monomials are combined, the six mandatory vertices of the
Newton polygon of `(5)` are

```text
(0,1), (2,0), (4,2), (4,3), (2,4), (0,4).              (10)
```

The forced `lambda` and `alpha` terms also supply `{(2,1),(0,2)}` and
`{(2,2),(0,3)}`, while the forced shared coefficient `kappa-epsilon` supplies
`(2,3)`. The only optional pair is

```text
{(3,3),(1,4)} from phi.                                  (11)
```

They do not change the polygon. In particular, `phi=0` is included. In
counterclockwise order,

```text
Delta=conv{(0,1),(2,0),(4,2),(4,3),(2,4),(0,4)}.        (12)
```

For an edge `e`, let `nu_e` be its primitive inward normal, `ell_e` its
lattice length, and `d_e` the lattice distance from the residue monomial
`(1,1)` to its supporting line. The toric order of `(6)` at every boundary
point on `e` is `d_e-1`. Indeed, in a one-edge toric chart with boundary
parameter `z` and edge coordinate `w`, write

```text
ord_z(s^i p^j)=<nu_e,(i,j)>,
F_Q=z^h(f_e(w)+O(z)),                                  (12a)
```

where `h` is the supporting constant. At a simple geometric root of `f_e`,
use `ds wedge dp=sp dlog(s) wedge dlog(p)` in `(6)`. Taking the residue in
`w` leaves a unit times

```text
z^(<nu_e,(1,1)>-h-1) dz=z^(d_e-1) dz.                 (12b)
```

The complete ledger is

| edge | `nu_e` | `ell_e` | `d_e` | residue order | `(ord s,ord p; ord t,ord x)` |
|---|---:|---:|---:|---:|---:|
| `(0,1)--(2,0)` | `(1,2)` | 1 | 1 | 0 | `(1,2;2,-1)` |
| `(2,0)--(4,2)` | `(-1,1)` | 2 | 2 | 1 | `(-1,1;-2,1)` |
| `(4,2)--(4,3)` | `(-1,0)` | 1 | 3 | 2 | `(-1,0;-2,1)` |
| `(4,3)--(2,4)` | `(-1,-2)` | 1 | 7 | 6 | `(-1,-2;6,-7)` |
| `(2,4)--(0,4)` | `(0,-1)` | 2 | 3 | 2 | `(0,-1;-1,1)` |
| `(0,4)--(0,1)` | `(1,0)` | 3 | 1 | 0 | `(1,0;0,1)` |

The edge polynomials are squarefree under `(1)`: the only nontrivial top
quadratic has discriminant `Delta_V`, the long `s=0` edge is squarefree by
THM-4053's generic-completion gate, and the remaining edges are linear or
the visibly separable quadratic joining `(2,0)` to `(4,2)`. Thus an edge of
length `ell_e` gives exactly `ell_e` boundary points with the displayed
order.

For the first edge, its equation gives
`p~(1+gamma Q)s^2`, hence `t~gamma Qs^2`. On the index-seven edge,
`ord(H)=-8`, and

```text
t(1-QH)=gamma Qs^2                                     (13)
```

forces `ord(t)=6`. These are the only two cancellation-sensitive rows in the
last column. The other rows follow directly from `t=p-s^2` and `x=s/t`.

Finally,

```text
2 Area(Delta)=24,       |boundary lattice steps|=10,
# interior lattice points=(24-10+2)/2=8.               (14)
```

This recovers the genus-eight ledger of THM-4053 from the generic polygon.

## 3. The boundary saturates Riemann--Hurwitz

The last edge in the table is `s=0`. Its three squarefree nonzero roots have

```text
(x,t)=(0,p),                                               (15)
```

so they are affine source points and have ramification index one. Every
point on the other five edges is at source infinity. Reading index as one
plus the residue order gives exactly

```text
source-infinity indices: {1,2,2,3,3,3,7},
affine s=0 indices:      {1,1,1}.                        (16)
```

Consequently the seven source-infinity points contribute

```text
sum_(P at infinity)(e_P-1)=14.                          (17)
```

On the other hand, Riemann--Hurwitz for a map from genus eight to genus one
gives

```text
deg R_varphi=2*8-2=14.                                  (18)
```

Equations `(9)`, `(17)`, and `(18)` show that `(16)` is the entire
ramification divisor. In particular, there is no hidden affine or interior
ramification.

## 4. Exact degree sieve and two forced responses

Let `O` be the target point at infinity and `n=deg(varphi_q)`. Since `A,C`
are polynomials, no affine source point maps to `O`; hence

```text
n=sum_(P maps to O)e_P <=1+2+2+3+3+3+7=21.             (19)
```

The index-seven point exists whether it maps to `O` or to a finite target
point. A local degree never exceeds the global degree, so

```text
7<=n<=21.                                               (20)
```

THM-4053 proves independently that `n` is an Eisenstein norm
`a^2-ab+b^2`, equivalently that every prime congruent to `2 mod 3` has even
valuation. Intersecting that criterion with `(20)` leaves exactly

```text
n in {7,9,12,13,16,19,21}.                             (21)
```

This is a strict finite sharpening of the unbounded norm gate in THM-4053.

Two boundary sidecars make the sieve much smaller.

### 4.1 The index-one puncture is forced to target infinity

On the first edge take `z=s`. Its edge equation and `(13)` give

```text
p=(1+gamma/q)z^2+...,       t=(gamma/q)z^2+... .        (21a)
```

THM-3992's exact extreme Laurent rows are

```text
A=gamma^2 s^2 t^-2+b(s)t^-1+...,
C=gamma^3 s^3 t^-3+d(s)t^-2+...,
b(s) in s^2 k[s],          d(s)=(3/2)gamma s b(s).      (21b)
```

All rows after the displayed extreme have strictly larger `z`-order on
`(21a)`. Therefore there is no leading cancellation and

```text
A~q^2 z^-2,                 C~q^3 z^-3.                (21c)
```

This geometric index-one puncture maps to `O`.

### 4.2 The index-two pair is one quadratic closed point

On the second edge use `Y=sp`. Its exact edge equation is

```text
(1+gamma/q)-(kappa/q)Y^2=0,
equivalently             kappa Y^2=q+gamma.            (21d)
```

The right side divided by `kappa` has valuation one at `q=-gamma`, so it is
not a square in `k(q)`. Thus `(21d)` is irreducible and separable over
`k(q)`: the two geometric index-two punctures form one degree-two closed
point. Both `varphi_q` and `O` are defined over `k(q)`, so Galois conjugacy
forces these two punctures either both to map to `O` or both to have finite
target image.

## 5. Raw response universe and Galois sharpening

Group the seven source-infinity points by indices `1,2,3,7`. Let

```text
(a,b,c,d)=numbers of those points having finite target image,
0<=a<=1, 0<=b<=2, 0<=c<=3, 0<=d<=1.                   (22)
```

Then the target-infinity fibre equation is exactly

```text
a+2b+3c+7d=21-n.                                       (23)
```

For the seven values in `(21)`, the complete coarse profile table is

| `n` | finite-image profiles `(a,b,c,d)` | fully labelled subsets |
|---:|---|---:|
| 7 | `(0,2,1,1)`, `(1,0,2,1)`, `(1,2,3,0)` | 7 |
| 9 | `(0,1,1,1)`, `(1,1,3,0)`, `(1,2,0,1)` | 9 |
| 12 | `(0,0,3,0)`, `(0,1,0,1)`, `(1,1,2,0)` | 9 |
| 13 | `(0,1,2,0)`, `(1,0,0,1)`, `(1,2,1,0)` | 10 |
| 16 | `(0,1,1,0)`, `(1,2,0,0)` | 7 |
| 19 | `(0,1,0,0)` | 2 |
| 21 | `(0,0,0,0)` | 1 |

Thus there are exactly 16 coarse profiles and 45 fully labelled subsets.
The generating certificate is

```text
(1+z)(1+z^2)^2(1+z^3)^3(1+z^7).                       (24)
```

Splitting the three index-three points into the singleton vertical-edge
point and the two horizontal-edge points gives exactly 23 edge-refined
profiles. This is the response universe before using Sections 4.1--4.2.

Equation `(21c)` forces `a=0`, and `(21d)` forces `b in {0,2}`. Filtering the
raw table leaves exactly

| `n` | necessary profile `(a,b,c,d)` | fully labelled | edge-refined |
|---:|---|---:|---:|
| 7 | `(0,2,1,1)` | 3 | 2 |
| 12 | `(0,0,3,0)` | 1 | 1 |
| 21 | `(0,0,0,0)` | 1 | 1 |

Hence every hypothetical survivor satisfies

```text
n in {7,12,21},                                           (24a)
```

with exactly three coarse, four edge-refined, and five fully labelled
necessary responses. A Galois-aware generating certificate before the
Eisenstein filter is

```text
(1+z^4)(1+z^3)^3(1+z^7),                                (24b)
```

where the `z^4` block is the Galois-locked pair of index-two
punctures. These counts are necessary response types, not symmetry orbits and
not realizability claims. They still forget finite target values, collisions
among those values, leading coefficients, and global compatibility of `A,C`.

## 6. Hostile control and exact audit

A highest-face or genus-only observer may erase the forced `kappa` vertex
`(4,2)`. The counterfactual polygon becomes

```text
conv{(0,1),(2,0),(4,3),(2,4),(0,4)}.                   (25)
```

It still has genus eight and total canonical degree fourteen:

```text
2 Area=22, boundary=8, interior=8,
boundary indices={1,1,1,1,3,3,5,7}.                   (26)
```

But the original packet `[2,2,3]` has been replaced by `[5]`. Therefore
genus, Riemann--Hurwitz total, and the highest face do not determine the
response; complete mandatory support is load-bearing.

The primary certificate checks all eight `lambda/alpha/phi` presence
patterns (the six deleting `lambda` or `alpha` are robustness hostiles), the
complete polygon, raw and Galois-refined profile ledgers, two independent
Eisenstein criteria, the support floors below, and `(25)--(26)`. A separate
no-SymPy audit reconstructs the polynomial term by term, certifies primitive
supporting half-spaces, and checks both the 45 raw and five surviving labelled
responses literally. Normal and
optimized executions reproduce the frozen outputs byte for byte.

## 7. The index-seven branch forces high polynomial support

Take `z=1/s` at the edge of index seven. The source equation and `(6)` give

```text
T=-gamma/theta,
t~T z^6,              x~T^-1 z^-7,
varphi_q^*(dU/(2V))~-(1/theta)z^6 dz.                  (27)
```

If this puncture maps to `O`, use the target parameter `tau=U/V`. At `O`,

```text
U=tau^-2+...,          V=tau^-3+...,
dU/(2V)=-d tau+... .                                    (28)
```

Equations `(27)--(28)` force

```text
tau~z^7/(7theta),
A~49theta^2 z^-14,     C~343theta^3 z^-21.              (29)
```

If the puncture has finite target image, instead `ord_z(A),ord_z(C)>=0`.

Now use only the exact Laurent-depth boundary from THM-3992. Apart from the
unique extreme terms `gamma^2 x^2` in `A` and `gamma^3 x^3` in `C`, every
monomial `x^m t^n` satisfies

```text
m<=n+1 in A,           m<=n+2 in C.                     (30)
```

Along `(27)`, the extreme terms are

```text
gamma^2x^2~theta^2 z^-14,
gamma^3x^3~-theta^3 z^-21.                              (31)
```

If no other `A` monomial has order at most `-14`, `(31)` is incompatible
both with finite image and with the coefficient `49theta^2` in `(29)`.
Therefore some other monomial satisfies

```text
-7m+6n<=-14,           m<=n+1.                          (32)
```

Thus `14<=7m-6n<=n+7`, so `n>=7`; the sharp first pair is `(m,n)=(8,7)`
and `m+n>=15`. The identical argument for `C` gives

```text
-7m+6n<=-21,           m<=n+2,                          (33)
```

whence `21<=7m-6n<=n+14`, `n>=7`, with sharp first pair `(9,7)` and
`m+n>=16`. Consequently every hypothetical survivor in this seam obeys

```text
deg_t A>=7,       deg_t C>=7,       deg A>=15,       deg C>=16. (34)
```

This does not assert that the sharp monomials themselves occur. A deeper
DE-weight face can vanish at its initial value and re-enter at orders
`-14,-21`; resolving that iterated initial-form cancellation is part of the
remaining sidecar.

## 8. Precise next frontier and representation firewall

At this stage the next decisive object was not another genus count, but the
four unfilled geometric rows of the seven-row Puiseux evaluation table

```text
P | edge | e_P | ord_P(A) | ord_P(C) | leading target value.          (35)
```

Sections 4.1--4.2 fill the AB row and lock the two BC rows. Completing the
CD, DE, and two EF rows selects one value in `(24a)` and one of the five
labelled responses. The present index multiset cannot do
that: it preserves local degree and total ramification but destroys target
image and edge-leading-coefficient data.

Accordingly there is no intrinsic tournament here. Equal indices create
ties, and arbitrarily orienting them loses the target fibre `(23)`. A
Kakeya-style direction-only shadow has the same defect: it retains boundary
directions but discards multiplicity and target image. Those representations
may become useful only with the full edge/Puiseux sidecar.

THM-4120 later completes the labelled response, and THM-4130 then excludes
the smooth theta-only survivor by critical monodromy. THM-4134/4138 later
exclude all four degree alternatives `20,19,16,15` on `Delta_V=0`, emptying
that wall. The other two walls and `JC(2)` remain **OPEN**.

**QED.**
