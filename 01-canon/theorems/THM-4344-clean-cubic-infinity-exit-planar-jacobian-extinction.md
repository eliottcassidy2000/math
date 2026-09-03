---
id: THM-4344
title: "Clean cubic infinity-exit planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327/4339 + VERIFIED-EXACT + TWICE
  INDEPENDENTLY HOSTILE-AUDITED. In the inherited reduced (2,3), exact-
  weight-twelve seam, Z=beta_11=zeta_3=W=0 and U*K*xi_10!=0 are extinct.
  The cubic degree drop exposes a genus-two D6 replacement face. Its sole
  variable infinity collision has one critical-value series; every finite
  split order is a rational, elliptic, or genus-two tail with two proved
  attachments and positive Keller-form order when nonrational. The finite
  quadratic collision is independently a rational A11 bridge or a horizontal
  node. Exact coupled support gives 192/192 M-D6-T fans; a 512-state Boolean
  calculation is only a conservative over-atlas. All complementary charts
  and simultaneous incidences are explicit. K=0, xi_10=0, U=0, seam entry,
  JC(2), and DC(2) are not claimed here.
source: root + degree-drop/local/global hostile-referee agents / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4339-clean-interior-cubic-edge-planar-jacobian-extinction
related:
  - THM-4341-odd-self-similar-cusp-reciprocal-tail-duality
  - THM-4342-clean-cubic-zero-exit-planar-jacobian-extinction
  - THM-4343-disjoint-a23-cubic-contact-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_thm4344.py
primary_output: 05-knowledge/results/jc2_m12_clean_cubic_infinity_exit_extinction_thm4344.out
primary_script_sha256: c65a442f0eecdde92d403c88aac58790e24f16387bf0c7869bf6daf13a9e0f5a
primary_output_sha256: dd9101db7ee1eff85e09adbe64fec2e164097641ffdfd5f6d9786ad9caf3373e
local_referee_script: 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_local_referee_thm4344.py
local_referee_output: 05-knowledge/results/jc2_m12_clean_cubic_infinity_exit_extinction_local_referee_thm4344.out
local_referee_script_sha256: eb282dad789b9e1a46183a656f47aea71ef7ee5519649ffee39585c158351e7e
local_referee_output_sha256: 763a16da5fa70fefd6f8c5b90ca8f3bcc10106bdf7bcf9c336d94034dc224b41
global_referee_script: 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_global_referee_thm4344.py
global_referee_output: 05-knowledge/results/jc2_m12_clean_cubic_infinity_exit_extinction_global_referee_thm4344.out
global_referee_script_sha256: cb6f7a3e11ee3e63992a07be738142c860dcd1253a2f7a984913d3f626356610
global_referee_output_sha256: 0d4cbd9f652f0e3cafa28e3ab51ddddef75c977e3caeb63826ef3453c17a3bd7
hash_basis: raw LF bytes
audit: >
  PASS AFTER COUNT AND CHART REPAIRS. The primary pins only THM-4327's hull
  engine and checks a 512-state conservative fan over-atlas, exact faces,
  polygons, orders, the root-at-infinity chart, all five tails, the finite
  collision, and an r=3 hostile. The import-free local referee reconstructs
  the D6 chart, proves the exceptional valuation primitive, resolves the
  even A5 contact, and checks 110 exact identities. The import-free global
  referee replaces the false "512 genuine masks" description by a 192-state
  coupled atlas, supplies the exact A11 infinity complement and all three
  finite blowup charts, and verifies simultaneous 3+3 and 1+1 sign incidence.
  Pick and Riemann--Hurwitz are checks; completeness comes from fan exhaustion
  plus chart coverage. Normal and optimized streams byte-match all three
  frozen outputs.
---

# THM-4344 -- Clean cubic infinity-exit planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327/4339 + VERIFIED-EXACT + TWICE
INDEPENDENTLY HOSTILE-AUDITED. THE DISPLAYED `W=0` GATE IS EXTINCT.
`K=0`, `xi_10=0`, `U=0`, SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam. Retain the
complete sixteen-term source of THM-4230 and the exact relations

```text
e=-1376/135,                    K=2848/45-(7/6)Delta.       (1)
```

Impose

```text
Z=beta_11=zeta_3=W=0,                 U*K*xi_10!=0.         (2)
```

Then no nonautomorphic planar Keller pair lies on `(1)--(2)`. Every other
lower coefficient is arbitrary, subject only to the inherited seam. The
claim is relative to the inherited toroidal and proper-flat target
interfaces; it proves neither entry into this seam nor `JC(2)` nor `DC(2)`.

The inheritance pass was:

- closest proved mechanism: THM-4339's Laurent root tax, reciprocal finite
  collision chart, and good-differential extinction;
- canonical hostile: a leading coefficient can vanish by exposing a new
  positive-genus face rather than by harmlessly deleting one root;
- corrected near miss: the initially reported 512 fan states independently
  toggle coupled coefficients and are an over-atlas, not 512 realizable
  specializations;
- least-used sidecar: THM-4341's warning that an even hyperelliptic tail has
  two infinity points. Here those points create a graph cycle and repair the
  otherwise missing unit of genus.

The live board was

```text
unique owner | replacement face | reciprocal root | even A5 contact
critical series | attachment graph | residue buffer | proper-flat degree. (3)
```

## 2. Exact support and the replacement fan

Put `xi=xi_10`, `u=upsilon_5`, and `alpha=alpha_11`. On `(2)`, the three
load-bearing support coefficients are uniquely owned:

```text
(2,6,1)=-U,                 (4,4,1)=-xi,
(4,2,1)=-K.                                             (4)
```

They cannot cancel. The three multiply-owned coefficients are

```text
c1=K+1376/135,             c2=Theta-Delta,
c3=xi-u.                                                (5)
```

The seam makes the coefficient statuses nonindependent:

```text
K=0  iff Delta=5696/105,
c1=0 iff Delta=3968/63.                                 (6)
```

The statuses `(Delta,c1,Theta,c2)` have eight realizable classes, the pair
`(u,c3)` has three, and `Phi,eta,alpha` contribute three independent bits.
Exact rational representatives therefore give

```text
8*3*2^3=192,
192/192 realizable support patterns: M,D6,T.             (7)
```

Independently toggling the six optional source rows and the three aggregate
deletions gives a conservative `2^6*2^3=512` atlas. Each exact pattern maps
to the over-atlas state having the same six row presences and actual three
aggregate cancellations, hence identical active support. The stronger check
is

```text
512/512 conservative states: M,D6,T.                    (8)
```

The distinction between `(7)` and `(8)` is load-bearing truth discipline:
the second count is not a coefficient realization census.

The exact face equations, up to torus monomials, are

```text
G_M =(P-S^2)(UP^6-1),

G_D6=S^2(1-UP^6-alpha*S*P^5-xi*S^2*P^4),

G_T =S^2(1-S^2*P^2 A(P)),
A(P)=K+Theta*P+xi*P^2.                                  (9)
```

Thus the exited cubic root opens the intervening `D6` face; it does not
merely lower the degree of `T`.

## 3. Edges, genera, and primary differential orders

The face and global polygon ledgers are

| object | vertices | `(2Area,boundary,interior)` |
|---|---|---:|
| `M` | `(0,1),(2,0),(2,6),(0,7)` | `(24,14,6)` |
| `D6` | `(2,0),(4,4),(2,6)` | `(12,10,2)` |
| `T` | `(2,0),(4,2),(4,4)` | `(4,6,0)` |
| global | `(0,1),(2,0),(4,2),(4,4),(2,6),(0,7)` | `(40,14,14)` |

The six outer boundary schemes, followed by the two internal schemes, are

```text
X-1,
1-KX^2,
K+Theta X+xi X^2,
U+alpha X+xi X^2,
X-1,
U-X^6;

M--D6: 1-UX^6,                 D6--T: 1-xi X^2.         (10)
```

The last outer scheme is extracted from the `M` initial form after retaining
the sigma height; projecting the full source first would incorrectly retain
higher affine rows. Under `(2)`, all fixed schemes in `(10)` are reduced.
Only

```text
D_fin=Theta^2-4Kxi,             D_inf=alpha^2-4Uxi       (11)
```

can vanish.

The global source-infinity packet is

```text
(11,6,6,3,3,2,2,1),            sum(e_i-1)=26=2*14-2.    (12)
```

Generically, `M` is its rational component `P=S^2` plus the six rational
components `UP^6=1`. Their twelve intersections, six `M--D6` attachments,
and two `D6--T` attachments give

```text
V=9,                   E=20,                   b1=12.    (13)
```

Section 4 identifies the remaining face genus as two, so `2+12=14`, exactly
the global Pick genus. The supporting-plane formula gives the positive
primitive orders

```text
M:9 on Q=sigma^12,       D6:5 on Q=sigma^6,
T:8 on Q=sigma^6.                                      (14)
```

Pick and `(12)` are independent checks; fan exhaustion and the complete
charts below prove component completeness.

## 4. Exact root-at-infinity chart

At the primitive `D6` scale put

```text
Q=sigma^6,       s=sigma^-1 S,       p=sigma^-1 P,
G=sigma^2 F_Q.                                           (15)
```

Define

```text
H6=UP^6+alpha SP^5+xi S^2P^4,
H5=uP^5+eta SP^4+Theta S^2P^3,
H4=Delta P^4+Phi SP^3+K S^2P^2.                         (16)
```

Expansion of the complete source, with no ellipsis, gives

```text
G=(S^2-sigma P)
  (1-H6-sigma H5-sigma^2H4-e sigma^3P^3
     -(8/3)sigma^4P^2+3sigma^5P)
  -sigma^6S^2/2.                                        (17)
```

Use the unimodular torus coordinates

```text
x=P^-1,                 v=S/P,                 rho=sigma x.
```

With

```text
A_inf(v)=U+alpha v+xi v^2,
B_inf(v)=u+eta v+Theta v^2,
C_inf(v)=Delta+Phi v+K v^2,                             (18)
```

multiplication by `x^8` transforms `(17)` exactly into

```text
Phi_inf=(v^2-rho)
 [x^6-A_inf-rho B_inf-rho^2C_inf-e rho^3
       -(8/3)rho^4+3rho^5]
 -v^2rho^6/2.                                           (19)
```

At every root of `A_inf`, `v!=0` because `U!=0`; hence `v^2-rho` is a unit.
Division yields one complete series equation

```text
x^6=D(v,rho),

D=A_inf+rho B_inf+rho^2C_inf+e rho^3+(8/3)rho^4-3rho^5
  +v^2rho^6/[2(v^2-rho)].                               (20)
```

On the face, put `Y=2xi v+alpha`. Then

```text
Y^2=4xi x^6+D_inf.                                      (21)
```

If `D_inf!=0`, the sextic is squarefree. Its smooth projective model has
genus two; its order five in `(14)` makes its map to the good elliptic target
constant.

## 5. The even A5 collision and its complete complement

Suppose `D_inf=0`. There is a nonzero `a` with

```text
U=xi a^2,              alpha=-2xi a,
A_inf(v)=xi(v-a)^2.                                    (22)
```

At `(v,rho)=(a,0)`, `D_vv=2xi` and the prefactor in `(19)` equals `a^2`.
The parameter-preserving formal Morse lemma therefore gives a unique
critical-value series `psi(rho)`, `psi(0)=0`, and a completed equation

```text
y^2=x^6-psi(sigma x).                                  (23)
```

This is an even `(2,6)` contact. It is not covered by THM-4341's odd-cusp
statement.

For example, let `B_1=eta+2Theta*a` and `C_1=Phi+2K*a`. The first three
critical coefficients are

```text
[rho]psi   =u+eta*a+Theta*a^2,
[rho^2]psi=Delta+Phi*a+K*a^2-B_1^2/(4xi),
[rho^3]psi=e+Theta*B_1^2/(4xi^2)-C_1*B_1/(2xi).         (24)
```

The allowed specialization

```text
xi=U=1, a=-1, alpha=2,
u=eta=Theta=Delta=0,             Phi=K=2848/45          (25)
```

has critical packet `(0,0,-1376/135)`. Thus the positive-genus `r=3` tail
really occurs, and no early-Hasse nonvanishing shortcut is valid.

Let `r=ord_rho psi`, `1<=r<6`, and write `psi=rho^r C(rho)`, `C(0)=c0!=0`.
Make the honest base change and weighted chart

```text
sigma=tau^[2(6-r)],        x=tau^[2r]X,
y=tau^[6r]Y.                                           (26)
```

After removing the square `X^[2 floor(r/2)]`, with `epsilon=r mod 2`, the
actual exceptional function field is

```text
Y_0^2=X^epsilon(X^(6-r)-c0).                            (27)
```

The valuation vector `(1,2r,6r)` is primitive because its first coordinate
is one. Together with the explicit `tau`-chart, this proves that `(27)` is
the invariant exceptional component, not a ramified root cover of the kind
prohibited by MISTAKE-540.

The complete tail table is

| `r` | tail equation | tail genus | persistent delta | global genus before `epsilon_fin` | `tau` form order |
|---:|---|---:|---:|---:|---:|
| 1 | `Y_0^2=X(X^5-c0)` | 2 | 0 | 14 | 56 |
| 2 | `Y_0^2=X^4-c0` | 1 | 1 | 13 | 52 |
| 3 | `Y_0^2=X(X^3-c0)` | 1 | 1 | 13 | 48 |
| 4 | `Y_0^2=X^2-c0` | 0 | 2 | 12 | 44 |
| 5 | `Y_0^2=X(X-c0)` | 0 | 2 | 12 | 40 |

The affine chart alone does not prove the attachment count. Put

```text
rho=sigma*x=tau^12 X,       z=1/X=tau^12/rho,
h=y/x^3.                                                (28)
```

The exact complementary chart is

```text
z*rho=tau^12,
h^2=1-z^(6-r)C(rho).                                    (29)
```

At `z=rho=tau=0`, the two points `h=+1,-1` are distinct and etale because
`2h` is a unit. They lie over the toric `A11` surface `z*rho=tau^12`;
ordinary resolution inserts two disjoint rational `A11` chains, one at each
sign, and omits no carrier.

Choose `lambda^2=xi`. The repeated `D6` face factors as

```text
1-xi P^4(aP-S)^2
 =[1-lambda P^2(aP-S)][1+lambda P^2(aP-S)].             (30)
```

Each factor is rational, receives three of the six `M--D6` attachments and
one of the two `D6--T` attachments, and carries one sign of the A5 contact.
Normalizing that contact changes `(13)` to `V=10,E=20,b1=11`; before the
independent finite-edge correction, the tail adds one vertex and two paths
and restores `b1=12`. Equivalently, every finite row satisfies the exact
local conservation law

```text
floor(r/2)+g_tail+1=3=delta(A5).                         (31)
```

The final `1` is graph genus from the two attachments. It is precisely the
information lost by retaining only the normalized tail polynomial.

The differential calculation is also exact. From `Phi_inf=x^8G`,
`P=x^-1`, and `S=v/x`,

```text
G_S=x^-7(Phi_inf)_v,        -dP/G_S=x^5 dx/(Phi_inf)_v.
```

After Morse preparation,

```text
phi^*eta_0=unit*sigma^5 x^5 dx/y.                       (32)
```

Its sigma-order is `5+3r/(6-r)`; multiplying by the base index
`2(6-r)` gives exactly `56,52,48,44,40`. In particular every nonrational
tail (`r=1,2,3`) is Keller-constant.

If `r>=6` or `psi=0`, `(23)` is `x^6` times a unit. The integral element
`y/x^3` and a unit square root separate two regular horizontal sheets. The
normalization is finite, proper, and DVR-flat; its persistent delta is three,
there is no vertical tail, and the normalized global genus before the
independent finite-edge correction is eleven.

## 6. The independent finite quadratic collision

The other variable edge is

```text
A(P)=K+Theta P+xi P^2.                                  (33)
```

When `D_fin!=0`, `z^2=A(P)` is a smooth rational `T` normalization. If
`D_fin=0`, then `Kxi!=0` gives a nonzero `a` with

```text
K=xi a^2,                 Theta=-2xi a,
A(P)=xi(P-a)^2.                                      (34)
```

In THM-4339's exact reciprocal chart use the doubled base
`Q=tilde_sigma^12`, put `delta=tilde_sigma^6` and `x=P-a`, and write
`B_a=B(a)` for `B(P)=P^3(Phi+eta P+alpha P^2)`. This `tilde_sigma` is not
the primitive `D6` parameter of `(15)`. The degree-two tangent is

```text
Q=b^2-a^2xi x^2-delta B_a b.                            (35)
```

The restrictions of the three strict-transform charts to the exceptional
divisor of the ordinary blowup are

```text
delta-chart: b=delta Y, x=delta X:
             Y^2-a^2xi X^2-B_aY=0;

b-chart:     delta=bD, x=bX:
             1-a^2xi X^2-B_aD=0;

x-chart:     delta=xD, b=xY:
             Y^2-a^2xi-B_aDY=0.                         (36)
```

If `B_a!=0`, these are the three affine pieces of one smooth projective
conic. It is a rational `A11` bridge, meets the two horizontal signs once
each, has good-form order 28 in `tilde_sigma`-units (equivalently 14 on the
primitive `T` base), and restores the graph cycle lost when `T` splits.

If `B_a=0`, formal preparation instead has tangent

```text
q b^2-delta B'(a)xb-a^2xi x^2,

Disc_b=x^2[delta^2B'(a)^2+4q a^2xi].                    (37)
```

Here `q` is a DVR unit with `q(0)=1`. The bracket is therefore a DVR unit.
A unit square root gives two regular horizontal branches and no vertical
carrier; the normalized genus drops by one.

## 7. Simultaneous collisions and the complete genus ledger

The infinity and finite collision centers lie in different open outer-edge
orbits, but disjointness alone does not prove their genus corrections add.
When both discriminants vanish, label the split faces by signs. On their
shared edge the coordinate

```text
t=SP^2=v/x^3=P/zeta,                 zeta=1/(SP)         (38)
```

takes the values `+/-1/lambda`. Thus `D6+` joins `T+` and `D6-` joins `T-`.
The six `M--D6` nodes split `3+3`. Since the seven-component `M` graph is
connected, the twice-split graph is connected and has

```text
V=11,                  E=20,                   b1=10.    (39)
```

Each two-ended infinity tail or finite conic bridge adds one vertex and two
edges, hence exactly one graph cycle. A persistent infinity normalization or
finite horizontal node adds neither. The complete normalized genus formula is

```text
g_norm=14-delta_inf-epsilon_fin,                         (40)

delta_inf=0                         if D_inf!=0,
          floor(r/2)                if 1<=r<6,
          3                         if r>=6 or psi=0;

epsilon_fin=0                       if D_fin!=0, or if D_fin=0 and B_a!=0,
            1                       if D_fin=0 and B_a=0.
```

For a concrete nonempty simultaneous wall, take

```text
xi=K=1, Theta=-2, U=4, alpha=-4, Delta=5606/105.        (41)
```

It satisfies `(1)` and both discriminant equations.

## 8. Proper-flat extinction

After finite base change, toric regularization, and the displayed
normalizations, every special component is one of:

1. a rational component from `M`, `T`, a split face, a toric chain, a finite
   bridge, or an `r=4,5` infinity tail;
2. the smooth genus-two `D6` component, of order five;
3. an `r=1,2,3` genus-two or genus-one tail, with positive integral order;
4. a regular horizontal normalization branch.

Thus every component map to the good characteristic-zero elliptic target is
constant. Retaining the actual positive component multiplicities after a
common dominating base, the inherited proper-flat identity gives, for a
positive-degree target line bundle `L`,

```text
deg(phi_generic^*L)=sum_Gamma m_Gamma deg(phi_Gamma^*L)=0. (42)
```

This contradicts the positive generic response degree of a nonautomorphic
Keller pair and proves `(1)--(2)`. **QED.**

## 9. Unexpected connections and lawful compression

The root-exit reciprocity from THM-4339 survives only with an address
sidecar. For

```text
A(P)=K+Theta P+xi P^2+WP^3,
A^vee(P)=P^3A(P^-1),                                    (43)
```

`W=0` is a zero exit of `A^vee`. But inversion of the cubic polynomial does
not carry the ambient Keller fan to the `K=0` fan: the adjacent normal layer
here exposes the new `D6` polynomial `U+alpha v+xi v^2`. The lawful packet is

```text
(infinity tax, zero tax, collision length;
 adjacent normal polynomial, owner divisor, residue buffer).             (44)
```

The even-cusp calculation also sharpens THM-4341's parity boundary. The
Newton slopes `r/(6-r)` and `(6-r)/r` are reciprocal, and the differential
excesses

```text
E_r=3r/(6-r), 1<=r<=5,        satisfy E_r E_(6-r)=9.     (45)
```

This is an excess law in common sigma units, not a product law for the full
orders `5+E_r`. Unlike the odd case, tail genus plus persistent delta is one
short; the two-ended graph cycle in `(31)` supplies exactly the missing unit.

Within the fixed exponent-six `(2,6)` block, the five oriented split types
are indexed by the natural number `r`. Quotienting by reciprocity gives

```text
h=min(r,6-r) in {1,2,3},
orientation=sign(2r-6), with h=3 fixed.                 (46)
```

The odd square `(2h-1)^2` may be compressed to `h`, but only after retaining
the cusp-order block, the orientation bit, the infinity-edge owner, and the
two-attachment flag. A bare natural number loses exactly the graph term in
`(31)`.

Finally, valuation comparison is a preorder, not a tournament. Away from a
tie, the lower-valuation monomial wins. At equality, the initial polynomial
`(27)` carries genus and two attachment addresses. Arbitrarily orienting that
tie destroys the carrier rather than resolving it.

The reusable organizer is the replacement-face principle: after a highest
support vertex vanishes, pin the first nonzero unique owner in each adjacent
normal direction and recompute the hull before compressing the root packet.
Outside the present seam this remains a research method, not a theorem.

## 10. Reproduction and stopping boundaries

Run from the repository root:

```bash
python3 -B 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_thm4344.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_thm4344.py
python3 -B 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_local_referee_thm4344.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_local_referee_thm4344.py
python3 -B 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_global_referee_thm4344.py
python3 -B -O 04-computation/jc2_m12_clean_cubic_infinity_exit_extinction_global_referee_thm4344.py
```

The three normal/optimized pairs byte-match their frozen outputs. What is
proved is exactly `(1)--(2)`, relative to the inherited seam and proper-flat
interfaces. The next sharp gates are

```text
xi_10=0:    the D6 unique owner vanishes and a second infinity fan appears;
K=0:        infinity and zero exits meet;
U=0:        the root chart's unit and endpoint owner fail.                (47)
```

None is imported into this theorem.
