---
id: THM-4045
title: "Live reduced 2:3 max-seven hidden elliptic-tail no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the live b=d=0 seam
  of the oriented reduced (2,3) planar-Jacobian cell, exact total residual
  weight seven is impossible. The complete Q=q^-1 lower Newton subdivision
  has three faces: a pair of rational components meeting at seven nodes, a
  unique elliptic tail with j=1728, and a rational vertical face. After
  Q=sigma^84, every other boundary and resolution component is rational,
  while the target has good j=0 reduction. The CM-field mismatch forces the
  entire connected special map to be constant, contradicting conservation
  of the positive generic degree. Combined with THM-4012, this raises the
  unconditional floor on this seam from M>=7 to M>=8. The maximum-weight-
  seven hypothesis is load-bearing: a p*y^2 term of weight eight destroys
  the elliptic side facet, so weight eight and JC(2) remain open.
source: jc2-merged-frontiers-0824 / complete lower-Newton boundary audit, 2026-08-24
audit: >
  PASS. The primary certificate reconstructs the expanded residual support,
  enumerates exactly three lower faces, verifies all eight edge restrictions,
  the seven transverse main nodes, their exact A_83 smoothing charts, the
  tail quartic invariants, Pick and dual-graph ledgers, inner chains, target
  scaling, and the first weight-eight hostile. A SymPy-free audit independently
  expands and combines support, enumerates the lower hull, runs polynomial
  Euclidean squarefreeness checks, and reproduces the genus, slope, local,
  target, and hostile ledgers. Normal and optimized streams match both frozen
  outputs.
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
related:
  - THM-4008-pure-p-residual-totally-degenerate-generic-fibre-no-go
  - THM-4011-companion-observer-kernel-etale-log-rh-and-endpoint-gate
  - THM-4017-sharp-weight-eight-specialization-obstruction-and-newton-ledger
  - THM-4020-live-two-three-fourth-row-max-seven-exclusion
  - THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14. The general model theorem
  supplies the face/edge model and rational toric chains; the seven singular
  points of the main face are resolved by the exact local calculation here.
script: 04-computation/jc2_max7_hidden_elliptic_tail_thm4045.py
output: 05-knowledge/results/jc2_max7_hidden_elliptic_tail_thm4045.out
script_sha256: 7e5746f62400331d305e41619f3a804818a7dd2ce016ebb3e7b0808882112bcc
output_sha256: 2a6b519227d3ca6a796cfd646b89096dad4ece0fb9d3b619c01d8b711545480f
independent_audit_script: 04-computation/jc2_max7_hidden_elliptic_tail_thm4045_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_max7_hidden_elliptic_tail_thm4045_independent_audit.out
independent_audit_script_sha256: ac664386c78644f87d16bfedce224234d4b45aa84e55b4ca1d2dbb0428a420ce
independent_audit_output_sha256: 230fbfd2965646a950a039f35d4ef9c4fb26a9dc107e3178d08192e873fd41b7
independent_audit_semantic_sha256: 2b4745d3cff4ae93dd75308f5edf89f77c9c4412c0e91be3fc016d14d42b30ca
hash_basis: raw LF bytes
---

# THM-4045 -- the hidden weight-seven elliptic tail closes the lane

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. This theorem concerns
one precisely inherited lane of the planar Jacobian conjecture; it does not
claim that the reduced cell is empty.

## 1. Exact maximum-seven support

Use THM-3992's normalized reduced `(2,3)` cell:

```text
s=xt,                 p=s^2+t,                 y=sp,
G=gamma*s^2/t+H(p,y),                            gamma!=0.       (1)
```

Write `wt(p)=2`, `wt(y)=3`, and let `M` be the maximum weight of the
**entire polynomial** `H`. On the live THM-4007 seam

```text
b=[y](R/gamma)=0,          d=[py](R/gamma)=0,           (2)
```

THM-4012 proves `M>=7`. Suppose for contradiction that `M=7`. There is only
one monomial of weight seven, so the complete polynomial is

```text
H=lambda*p+alpha*p^2+epsilon*p^3+kappa*y^2+phi*p^2*y,
phi!=0.                                                       (3)
```

No ellipsis is present in `(3)`. In particular, the weight-eight monomials
`p^4` and `p*y^2` are zero. THM-4007's forced rows give, with `A5=a^5`,

```text
epsilon/gamma= 2752/(135 A5^3),
kappa/gamma  =-5696/( 45 A5^3),
kappa-epsilon=-3968 gamma/(27 A5^3).                    (4)
```

Thus `epsilon`, `kappa`, `phi`, and the expanded shared coefficient
`kappa-epsilon` are all nonzero. This prevents any silent deletion from the
Newton support below.

## 2. The complete three-face lower subdivision

Put `Q=q^-1`. The generic source fibre of THM-4012 is equivalently

```text
F_Q(s,p)=(s^2-p)(1-QH(p,sp))+gamma*Q*s^2=0.             (5)
```

Take coefficientwise `Q`-adic valuations: in particular, the coefficient
`1+gamma*Q` of `s^2` is one unit of height zero, not two support points.
After every other coincident monomial is combined, exact lower-hull
enumeration gives the six projected vertices

```text
A=(0,1), B=(2,0), C=(4,2), D=(3,3), E=(1,4), F=(0,4), (6)
```

and exactly three lower faces:

```text
face       vertices     height function nu(i,j)
M          A B D E      (i+2j-2)/7
T          B C D        (i+j-2)/4
V          A F E        (j-1)/3.                        (7)
```

The coefficients on these faces give, after removal of a monomial factor
where necessary,

```text
g_M=(S^2-P)(1-phi*S*P^3),
g_T=1-kappa*S^2*P^2-phi*S*P^3,
g_V=-1+epsilon*P^3+phi*S*P^3.                           (8)
```

This is the complete lower subdivision of `(5)`, not a truncation or a
highest-face proxy.

## 3. The three normalized faces

### 3.1 Main face: rational components and six graph cycles

The two factors of `g_M` are rational. They meet where

```text
P=S^2,                    phi*S^7=1.                    (9)
```

There are seven such points. The determinant of the two gradients on
`P=S^2` is

```text
-7 phi S^6 !=0,                                           (10)
```

so all seven intersections are transverse. Normalizing therefore gives two
rational components. Moreover,

```text
abs det((2,-1),(1,3))=7,                                (10a)
```

so these torus nodes saturate the binomial intersection number and there is
no further boundary intersection between the two completions. Their
seven-node dual graph has

```text
b_1=7-2+1=6.                                             (11)
```

The six interior lattice points of the polygon `ABDE` are accounted for by
exactly these six graph cycles; there is no abelian component on this face.

### 3.2 Tail face: the missing elliptic component

On `g_T=0`, put

```text
T_0=SP,                    W=2*kappa*T_0+phi*P^2.        (12)
```

Then

```text
W^2=phi^2 P^4+4 kappa.                                  (13)
```

This quartic is squarefree because `phi*kappa!=0`. Its binary-quartic
invariants are

```text
I=48 phi^2 kappa !=0,                J=0,               (14)
```

so its smooth projective normalization is an elliptic curve

```text
E_1728,                         j(E_1728)=1728.          (15)
```

The two torus-gradient equations reduce to
`2*kappa*S+phi*P=0` and `2*kappa*S+3*phi*P=0`; their determinant is
`4*kappa*phi!=0`. Thus the tail face itself is smooth, not merely its
normalization.

This is the positive-genus boundary tail lost by the singleton highest-face
observer in THM-4012.

### 3.3 Vertical face: rational

On the torus,

```text
partial g_V/partial S=phi*P^3 !=0,
S=(1-epsilon*P^3)/(phi*P^3).                            (16)
```

Thus `V` is smooth and rational. Its polygon has no interior lattice point.

## 4. Boundary and total-space validity gate

Make the single ramified base change

```text
Q=sigma^84.                                              (17)
```

The three height functions in `(7)` become the integral affine functions

```text
M: 12i+24j-24,       T:21i+21j-42,       V:28j-28.      (18)
```

Hence every face and edge denominator `delta_F,delta_L` is one, and all face
multiplicities are one. The eight edge restrictions,
up to nonzero monomials and constants, are

```text
AB -1+T       BD 1-phi*T          BC 1-kappa*T^2
CD kappa+phi*T                    DE -1+T
AE -1+phi*T   EF epsilon+phi*T    FA -1+epsilon*T^3.    (19)
```

They are all squarefree in characteristic zero. For the primitive inner edge
`BD`, Definition 3.12 with `L*=6-3i+j`, `P0=B`, `P1=(2,1)` gives slopes

```text
24 > 23 > 22 > 21,                                     (19a)
```

so its chain has length two. For `AE`, use `L*=3i-j+1`, `P0=A`,
`P1=(0,0)`; the determinant-one sequence is

```text
-24 > -25 > -26 > -27 > -28,                           (19b)
```

so its chain has length three. Every displayed denominator is one, hence all
chain multiplicities are one. Each outer slope is integral; Definition 3.12
sets its second slope to `floor(s_1-1)=s_1-1`, so `r=0` and no outer chain is
inserted.

There is a separate generic-completion gate. On the torus curve `(5)`,
`p-s^2` cannot vanish, since substitution would give
`F_Q=gamma*Q*s^2!=0`. Hence

```text
t=p-s^2,                    x=s/t                         (19c)
```

recovers `k(x,t)=k(s,p)`, so the torus curve is a dense open of THM-4012's
smooth generic source `C_q`. Indeed, the Jacobian of `(x,t)->(s,p)` is `t`,
so the displayed inverse is an isomorphism on these opens. Its six
**generic** outer-edge restrictions are

```text
AB -1+(1+gamma*Q)T       BC (1+gamma*Q)-Q*kappa*T^2
CD kappa+phi*T            DE -1+T
EF epsilon+phi*T          FA -1+Q*lambda*T+Q*alpha*T^2
                              +Q*epsilon*T^3.            (19d)
```

The first five are visibly squarefree over `k((Q))`. The discriminant of the
last cubic has `Q`-adic leading term

```text
-27 epsilon^2 Q^2 !=0.                                  (19e)
```

Thus the toric generic completion is smooth. By uniqueness of the smooth
projective model of a function field, it is the actual `C_q` in the
Keller-induced morphism, not an auxiliary singular completion.

The only failure of Delta-v regularity is the seven-node locus `(9)`. In the
main chart

```text
s=sigma^-12*S,             p=sigma^-24*P,               (20)
```

write

```text
H_sigma=phi*S*P^3
       +sigma^12(epsilon*P^3+kappa*S^2*P^2)
       +sigma^36 alpha*P^2+sigma^60 lambda*P,

U=S^2-P,                 V=(1-H_sigma)/S^2.             (21)
```

Transversality `(10)` makes `U,V,sigma` local coordinates. Dividing the
exact scaled equation by the unit `S^2` gives

```text
U*V=-gamma*sigma^84.                                    (22)
```

Thus each node is an `A_83` smoothing. Its regular resolution inserts 83
rational curves of multiplicity one and no positive genus.

Dokchitser's general face/edge theorem (Definitions 3.7, 3.9, 3.12 and
Theorem 3.14 in [*Models of curves over DVRs*,
arXiv:1807.00025v2](https://arxiv.org/abs/1807.00025)) applies globally and
supplies the proper flat model, multiplicities, and chains. Its regularity
clauses cover every smooth face, edge, and boundary locus. Here `T,V` are
smooth, both main branches are smooth off `(9)`, and every edge scheme in
`(19)` is smooth. Thus the seven torus nodes are its only nonregular points,
and `(22)` resolves all of them.

For an independent genus checksum, Pick's theorem gives

```text
polygon       twice area     boundary points     interior points
ABCDEF             21               9                   7
ABDE                14               4                   6
BCD                  4               4                   1
AFE                  3               5                   0.              (23)
```

The generic audit `(19c)--(19e)` makes Baker's polygon-genus equality apply,
so the global count in `(23)` is exact rather than merely an upper bound.

After normalization and regular resolution, the complete positive-genus
inventory is therefore

```text
one elliptic component E_1728;
every other irreducible component rational;
b_1(dual graph)=6;
g_generic=1+6=7.                                       (24)
```

Both internal face edges are primitive, so they attach the tail and vertical
face by one path each and create no additional graph cycle. This also proves
that the special fibre is connected.

## 5. Good target reduction and the no-Hom contradiction

The same base change is `q=sigma^-84`. In THM-4012's target equation, set

```text
A=sigma^-28*X,                 C=sigma^-42*Y.           (25)
```

The exact target model is

```text
Y^2=X^3+1-(3a^2/4)sigma^56 X-(a^3/4)sigma^84,          (26)
```

with good special fibre

```text
E_0:Y^2=X^3+1,                         j(E_0)=0.         (27)
```

In characteristic zero, an elliptic curve with `j=1728` has rational
endomorphism algebra `Q(i)`, while one with `j=0` has rational endomorphism
algebra `Q(sqrt(-3))`. Any nonconstant morphism becomes an isogeny after
translation by the negative of the image of the origin; a nonzero
homomorphism is therefore an isogeny and identifies the rational
endomorphism algebras. Consequently every morphism `E_1728 -> E_0` is
constant and

```text
Hom(E_1728,E_0)=0.                                      (28)
```

Base-change the finite nonconstant generic morphism

```text
C_q --> E_q                                             (29)
```

forced by a hypothetical Keller pair. Resolve its rational extension to the
regular source model. Every additional exceptional component is rational.
All rational components map constantly to `E_0` by Riemann--Hurwitz, and
`(28)` makes the unique elliptic component constant as well. Constants agree
across the connected special fibre.

Let `L` be relatively ample on the smooth target model. With actual positive
special-fibre multiplicities `m_i`, constancy would give

```text
deg(phi_generic^*L)
 =sum_i m_i deg((phi_i)^*(L|E_0))=0.                    (30)
```

But the generic map is finite and nonconstant, so the left side is positive.
This contradiction proves

```text
b=d=0 and M=7 is impossible.                            (31)
```

Combining `(31)` with THM-4012's unconditional `M>=7` gives the strengthened
live-seam floor

```text
b=d=0  ==>  M>=8.                                       (32)
```

## 6. Sharp boundary, old reservation, and connection contract

The first uncontrolled higher monomial is `p*y^2`, of weight eight. Its
expanded lifted endpoint `(4,3,1)` has gap

```text
1-(4+3-2)/4=-1/4                                       (33)
```

below the tail plane in `(7)`. It destroys the `j=1728` facet, exactly as
the THM-4017 hostile warns. Therefore `(32)` does not exclude weight eight,
does not extrapolate from a finite Hasse packet, and does not prove `JC(2)`.

THM-4020 had honestly reserved the same desired endpoint for an unaudited
fourth-row calculation. No step of that route is used here. It is now marked
`SUPERSEDED / UNPROVED RESERVED ROUTE`: the endpoint is proved by the
complete lower-Newton boundary mechanism above, while its proposed
fourth-row derivation remains outside the proof graph.

The typed cross-frontier connection is:

```text
source:       complete exact max-seven residual H;
target:       regular special-fibre component inventory;
map:          Q-adic lower Newton model after Q=sigma^84;
preserved:    proper fibre degree and all positive-genus components;
destroyed:    coefficients above weight seven and generic attachment labels;
sidecar:      squarefree edges plus the seven exact UV node charts;
decisive test: Hom(E_1728,E_0)=0;
hostile:      weight-eight p*y^2 lies below the tail plane.              (34)
```

This also explains the relation to THM-4044's sixty-clock firewall: a torus
observer can miss a boundary normal coordinate, while a single leading face
can miss a lower boundary tail. In both cases the repair is not more labels;
it is the missing boundary chart.

## 7. Reproduction

Run both independent exact paths in normal and optimized modes:

```bash
python3 -B 04-computation/jc2_max7_hidden_elliptic_tail_thm4045.py
python3 -B -O 04-computation/jc2_max7_hidden_elliptic_tail_thm4045.py
python3 -B 04-computation/jc2_max7_hidden_elliptic_tail_thm4045_independent_audit.py
python3 -B -O 04-computation/jc2_max7_hidden_elliptic_tail_thm4045_independent_audit.py
```

All four streams match their frozen LF outputs. **QED.**
