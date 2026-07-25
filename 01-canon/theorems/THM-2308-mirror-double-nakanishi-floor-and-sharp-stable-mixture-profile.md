---
id: THM-2308
title: "Mirror-double Nakanishi floor and sharp stable mixture profile"
status: >
  PROVED + CITED APPLICATION. For every maximal ideal of the Laurent
  coefficient ring, the residue-field dimension of a knot's Alexander
  module is additive under connected sum and is bounded above by the
  Nakanishi index and hence by unknotting number. If a knot has nontrivial
  Alexander polynomial, one maximal-ideal fiber is simultaneously nonzero
  for the knot and its mirror. Therefore every connected sum of P copies of
  the knot and Q copies of its mirror has unknotting number at least P+Q; in
  particular every n-fold mirror double costs at least 2n. For the
  Brittenham--Hermiller T(2,7)-mirror stable plane this raises the diagonal
  stable norm from the previous interval [1,5] to [2,5]. The full positive
  mixture defect is a symmetric homogeneous concave superadditive function,
  equivalently a one-variable concave profile. Mirror-averaged dual
  calibration and half-signature give optimal pointwise envelopes, and two
  explicit unconditional norms attain the two envelopes. Thus the present
  knot input determines the diagonal and chamber walls exactly but cannot
  determine any further interior ray. The exact ordinary and homogenized
  unknotting numbers remain open, and no positive Gordian catalyst is
  produced.
source: codex-2026-07-25-knot-stable-refinement
depends_on:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2292-common-catalytic-section-and-helly-calibration-nerve
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
  - THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
  - "Zhiqing Yang, Unknotting number of the connected sum of n identical knots, Journal of Knot Theory and Its Ramifications 17 (2008), 253--255, DOI 10.1142/S0218216508006130."
---

# THM-2308 -- mirror-double floors and the stable mixture profile

THM-2292 proves that the first Brittenham--Hermiller counterexample spans a
genuine rank-two stable plane and confines its unknown geometry to two
same-sign chambers. Two questions remain:

```text
how large must the mirror-paired diagonal be?
what is the exact functional form of the residual same-sign uncertainty?
```

The first question has a stronger Alexander-module answer than applying
Yang's identical-copy theorem only to the composite mirror double. The
second has a sharp convex-geometric answer: the unknown object is one
concave projective response profile, and the currently known knot data have
optimal lower and upper envelopes.

## 1. Maximal-ideal fibers pack under connected sum

Put

```text
Lambda=Z[t,t^(-1)].
```

For a knot `J`, let

```text
A_J=H_1(E_tilde(J);Z)
```

be its Alexander module over `Lambda`, and let `m_N(J)` be its Nakanishi
index, the least number of `Lambda`-module generators of `A_J`. Yang recalls
the standard lower-bound chain

```text
u(J)>=m_N(J).                                        (1)
```

Fix a maximal ideal `q` of `Lambda`, write

```text
k_q=Lambda/q,
nu_q(J)=dim_(k_q)(A_J tensor_Lambda k_q).            (2)
```

The dimension is finite because an Alexander module is finitely presented.

> **Common-fiber packing theorem.** For all knots `J,L`,
>
> ```text
> nu_q(J#L)=nu_q(J)+nu_q(L),                         (3)
> nu_q(J)<=m_N(J)<=u(J).                             (4)
> ```
>
> Consequently, for any finite packet `J_1,...,J_s`,
>
> ```text
> u(J_1#...#J_s)>=sum_i nu_q(J_i).                  (5)
> ```

### Proof

A connected-sum Seifert surface has block-diagonal Seifert matrix. Its
Alexander presentation matrix is therefore block diagonal, giving

```text
A_(J#L) isomorphic to A_J direct_sum A_L.            (6)
```

Tensoring (6) with `k_q` proves (3).

If `A_J` has `r` module generators, their images span
`A_J tensor k_q` as a vector space. Hence

```text
nu_q(J)<=r.
```

Minimizing over generating sets and then using (1) proves (4). Iterating
(3) and applying (4) proves (5). QED.

The point is not merely that each summand has nontrivial Alexander
polynomial. The same maximal ideal must see every summand. This common
fiber is the coordinate which is lost when lower bounds are applied to the
summands separately.

## 2. A two-per-copy floor for mirror doubles

Let `K` have nontrivial Alexander polynomial `Delta_K`. A square Seifert
presentation of `A_K` has determinant equal to `Delta_K` up to a unit of
`Lambda`. Since `Delta_K` is a nonunit, choose a maximal ideal

```text
q superset (Delta_K).                                (7)
```

After reduction modulo `q`, the presentation matrix is singular. Therefore

```text
nu_q(K)>=1.                                         (8)
```

The mirror polynomial is

```text
Delta_(Kbar)(t)=Delta_K(t^(-1))
```

up to a unit. Alexander symmetry makes this associate to `Delta_K(t)`.
Thus the same `q` also makes a square presentation of `A_(Kbar)` singular,
and

```text
nu_q(Kbar)>=1.                                      (9)
```

Apply (5) to `P` copies of `K` and `Q` copies of `Kbar`.

> **Mirror-mixture floor.** If `Delta_K` is nontrivial, then for all
> integers `P,Q>=0`,
>
> ```text
> u((#^P K)#(#^Q Kbar))>=P+Q.                       (10)
> ```
>
> In particular, for `X=K#Kbar` and every `n>=1`,
>
> ```text
> u(#^n X)>=2n,
> u_hash(X)>=2.                                     (11)
> ```

Equation (10) is a repo-derived common-fiber refinement. Yang's published
theorem gives at least `n` for `n` identical copies of any one knot with
nontrivial Alexander polynomial; it does not state the two-per-copy
mirror-double conclusion (11).

The hypothesis is sharp for this mechanism. If `Delta_K=1`, there is no
maximal ideal of the form (7), and the Alexander module can vanish. Nor does
(10) say that Nakanishi index, unknotting number, or their stable limits are
additive in general: it is one additive lower certificate on the particular
common fiber `q`.

## 3. The improved Brittenham--Hermiller stable plane

Now take

```text
K=T(2,7),                  Kbar=mirror(K),
a=[K],                    b=[Kbar]                  (12)
```

in THM-2191's real stable knot space, and write

```text
p=u_hash,
e=a+b,
o=a-b,
r=p(e),
delta=6-r.                                          (13)
```

THM-2292 and the Brittenham--Hermiller certificate give

```text
p(a)=p(b)=3,
p(o)=6,
r<=5.                                               (14)
```

The Alexander polynomial of `T(2,7)` is nontrivial, so (11) improves the
former Yang-only diagonal floor:

```text
2<=r<=5,
1<=delta<=4.                                        (15)
```

For a completely explicit common fiber, use

```text
Delta_K=t^6-t^5+t^4-t^3+t^2-t+1,
q=(2,t^3+t+1).                                      (15a)
```

Modulo two,

```text
Delta_K=(t^3+t+1)(t^3+t^2+1).
```

The cubic `t^3+t+1` is irreducible over `F_2`, so `q` is maximal and
`k_q` is the eight-element field. This is a concrete witness for (8)--(9),
not a computation of either fiber's full dimension.

More generally, (10) and homogenization give the physical-cone floor

```text
p(Pa+Qb)>=P+Q                  for P,Q>=0.           (16)
```

For integer coefficients this follows by applying (10) to every repeated
sum and dividing by the copy number. Positive homogeneity gives rational
coefficients, and continuity of the finite-dimensional norm gives all real
coefficients in (16).

The half-signature functional `S` and mirror involution `R` obey

```text
S(e)=0,          S(o)=6,
R(e)=e,          R(o)=-o.                           (17)
```

The seminorm is nondegenerate on this plane by THM-2292 and is invariant
under `R`. Norm symmetry together with `R` makes it unconditional in the
`(e,o)` coordinates:

```text
p(s e+t o)=p(|s|e+|t|o)                             (18)
```

where (18) denotes equality after the evident independent sign changes, not
an order relation between vectors.

For

```text
x=s+t,                  y=s-t,
s e+t o=x a+y b,                                      (19)
```

the condition `xy<=0` is exactly `|t|>=|s|`.
THM-2292's opposite-sign chamber law therefore becomes

```text
p(s e+t o)=6|t|                    if |t|>=|s|.      (20)
```

All remaining uncertainty lies in the central double cone
`|s|>|t|`.

## 4. The sharp norm envelope

Choose a dual unit functional `phi` which calibrates `e`:

```text
phi(e)=p(e)=r.                                      (21)
```

The mirror average

```text
G=(phi+phi after R)/2                               (22)
```

still has dual norm at most one and satisfies

```text
G(e)=r,                   G(o)=0.                   (23)
```

The dual unit ball also contains the half-signature `S` from (17). Hence
for all real `s,t`,

```text
p(s e+t o)>=max(r|s|,6|t|).                         (24)
```

For the reverse bound, unconditional symmetry reduces the central cone to
`s>=t>=0`. There

```text
s e+t o=(s-t)e+2t a,                                (25)
```

so (14) and the triangle inequality give

```text
p(s e+t o)<=r(s-t)+6t
             =r s+(6-r)t.                          (26)
```

Together with the exact outer law (20), this proves the global envelope

```text
max(r|s|,6|t|)
 <=p(s e+t o)
 <=r max(|s|,|t|)+(6-r)|t|.                        (27)
```

This both sharpens and symmetrizes THM-2292's chamber bound.

The envelope is optimal using only the displayed stable-plane data. For
each fixed `0<r<=6`, define two norms on `R^2` by

```text
p_low(s,t)=max(r|s|,6|t|),                          (28)

p_high(s,t)=r max(|s|,|t|)+(6-r)|t|.                (29)
```

The first is a weighted `l_infinity` norm. The second is the sum of
`r||(s,t)||_infinity` and the seminorm `(6-r)|t|`, hence is also a norm.
Both are unconditional and satisfy

```text
p_*(1,0)=r,
p_*(0,1)=6,
p_*(1/2,1/2)=p_*(1/2,-1/2)=3,                      (30)

p_*(s,t)=6|t|                    when |t|>=|s|.      (31)
```

Thus both models have the same values on `e,o,a,b`, the same mirror
symmetry, the same entire opposite-sign chamber, and, for `r>=2`, the same
physical floor (16). They attain respectively the lower and upper sides of
(27) at every point.

For the actual range `2<=r<=5`, the two envelopes differ strictly whenever

```text
|s|>|t|>0.                                          (32)
```

Indeed, if the maximum in (28) is `r|s|`, the gap is
`(6-r)|t|>0`; if it is `6|t|`, the gap is
`r(|s|-|t|)>0`. Therefore, within the still-unknown central cone, the
collected input forces exact values only on the diagonal ray `t=0` and its
walls `|s|=|t|`; it cannot force an additional central interior ray. The
entire outer cone is already exact by (20). The models are abstract norms,
not claims of knot realizability.

## 5. The mixture defect is one concave projective profile

For `P,Q>=0`, define the stable mixture defect

```text
D(P,Q)=3(P+Q)-p(Pa+Qb).                             (33)
```

For integer `P,Q`, this is the asymptotic saving rate

```text
D(P,Q)
 =lim_(n->infinity)
   [3n(P+Q)-u((#^(nP)K)#(#^(nQ)Kbar))]/n.           (34)
```

The function `D` is:

```text
nonnegative;
symmetric in P,Q;
positively homogeneous;
continuous and concave on R_(>=0)^2;
superadditive:
  D(P_1+P_2,Q_1+Q_2)>=D(P_1,Q_1)+D(P_2,Q_2).        (35)
```

Only concavity needs comment: `3(P+Q)` is linear and the restriction of the
norm `p` to the positive cone is convex. Superadditivity follows either
from (33) and the triangle inequality or from concavity plus positive
homogeneity.

Its boundary data are exact:

```text
D(P,0)=D(0,P)=0,
D(P,P)=delta P.                                     (36)
```

Write `P>=Q`. Substituting

```text
s=(P+Q)/2,                t=(P-Q)/2
```

into (27) gives the optimal defect sandwich

```text
delta Q
 <=D(P,Q)
 <=min(6Q, delta(P+Q)/2).                           (37)
```

The lower envelope is attained by `p_high`; the upper envelope is attained
by `p_low`.

Equivalently, put

```text
h(z)=D(1,z),                  0<=z<=1.              (38)
```

Then

```text
D(P,Q)=P h(Q/P)                 for P>=Q, P>0,
h is concave,
h(0)=0,                       h(1)=delta,            (39)

delta z
 <=h(z)
 <=min(6z, delta(1+z)/2).                            (40)
```

Concavity and `h(0)=0` imply that

```text
h(z)/z is nonincreasing on (0,1],
delta<=h(z)/z<=6.                                  (41)
```

Thus the stable saving **per minority summand** can only decrease as the
mixture approaches balance, and it reaches exactly `delta` on the balanced
ray. This is the functional form hidden by the binary symbiont relation.

Extending `h(z)=D(1,z)` to all `z>0`, mirror symmetry gives

```text
h(z)=z h(1/z).                                      (42)
```

Consequently, if `h` is differentiable at `1`, then

```text
h'(1)=delta/2.                                      (43)
```

Without differentiability, its one-sided derivatives satisfy

```text
h'_-(1)+h'_+(1)=delta,
h'_-(1)>=delta/2>=h'_+(1).                          (44)
```

Equations (38)--(44) reduce every stable mixture question in this
two-generator cone to one concave function on the compact interval
`[0,1]`.

## 6. Relation and information ledger

The operation-ready hierarchy on this plane is now:

```text
raw connected-sum costs
  -> homogenized norm p
  -> concave projective defect profile h
  -> signed pair gain
  -> unweighted symbiont edge.                      (45)
```

The map `p -> h` preserves every same-sign stable mixture magnitude and
uses homogeneity, mirror symmetry, and the known one-body cost to remove
scale and label redundancy. It loses the opposite-sign formal group
coordinates, but (20) reconstructs those exactly.

The map `h -> signed pair gain` remembers only which of the two relative
signs is additive. It loses `r`, `delta`, all interior values of `h`, and
the per-minority law (41). There is no intrinsic orientation: the honest
binary object is still THM-2292's switched gain edge. The magnitude sidecar
on the two-type quotient is the concave profile `h`, not a tournament.

The decisive next knot-specific tests are therefore:

```text
compute or improve r=u_hash(K#Kbar);
find one recognizable additive calibrator stronger than the
  common-Alexander-fiber floor on e;
or evaluate one interior rational slope h(Q/P).     (46)
```

Any one interior value strictly between the envelopes would rule out one or
both extremal models and constrain all other slopes by concavity. Conversely,
no argument using only mirror symmetry, half-signature, the common
Alexander-fiber floor, the Brittenham--Hermiller five-change certificate,
and norm axioms can improve (27) or (37).

This theorem does not compute `u(K#Kbar)` or `u_hash(K#Kbar)`, does not
promote the stable lower bound (10) to additivity, and does not produce a
positive Gordian catalyst. QED.
