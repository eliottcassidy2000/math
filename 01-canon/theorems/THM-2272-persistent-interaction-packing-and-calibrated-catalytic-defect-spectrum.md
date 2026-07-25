---
id: THM-2272
title: "Persistent interaction packing and the calibrated catalytic defect spectrum"
status: >
  PROVED (abstract commutative nonexpansive metric-monoid theorem) + CITED
  APPLICATION (Brittenham--Hermiller and the classical signature bound).
  The infimum over diagonal contexts of a finite packet's interaction defect
  is itself a normalized superadditive Boolean game. Consequently every
  disjoint hyperedge packing gives a simultaneous lower bound, and the
  weighted pair shadow gives a maximum-weight-matching certificate. In an
  integer metric its zero complex is exactly the some-context zero complex.
  For a packet whose atoms are rigid under every common context, the entire
  context spectrum is the root defect translated by the catalytic-saving
  spectrum of the merged composite; its minimum is the root defect, its
  supremum is the catalytic defect, and replicated spectra collapse
  uniformly to one stable rate. For individually additively calibrated knot
  atoms, Hahn--Banach identifies that rate with the exact deficit between
  independent atomic calibrations and one coherently signed common
  calibration. For K=T(2,7) and its mirror, p
  positive and q negative occurrences have context-uniform defect at least
  min(p,q), more sharply (6-u(K#mirror(K)))min(p,q). The pair-spectrum span
  of K,mirror(K) is exactly the catalytic capacity of their composite.
  This computes neither u(K#mirror(K)) nor a positive Gordian catalyst.
source: knot-refinement-probe-2026-07-25
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
  - THM-2259-boolean-continuation-hasse-field-and-signed-interaction-dividends
related:
  - THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
---

# THM-2272 -- persistent interaction packing and calibrated spectra

THM-2259 retains a separate Boolean defect game over every diagonal
continuation. There are two useful ways to compress that stalk, and they
answer different questions:

```text
infimum over contexts: interaction which cannot be removed by a context;
range over contexts:   interaction revealed by catalytic contraction of the
                       merged composite.
```

The first compression remains superadditive and therefore admits packing
certificates. Under an exact calibration hypothesis, the second has a closed
form and its width is precisely catalytic capacity.

## 1. The persistent continuation-floor game

Let `(M,+,0)` be a commutative monoid with a metric `d` satisfying joint
nonexpansivity

```text
d(a+c,b+e) <= d(a,b)+d(c,e).                         (1)
```

Put

```text
ell(x)=d(x,0),
rho_z(x)=d(x+z,z).                                   (2)
```

Fix a labelled packet `x_1,...,x_r`. For `S subset [r]`, write

```text
x_S=sum_(i in S)x_i,

Delta_z(S)
 =sum_(i in S)rho_z(x_i)-rho_z(x_S).                 (3)
```

THM-2259 proves that every `Delta_z` is nonnegative, normalized on the empty
set and singletons, and superadditive on disjoint coalitions. Define its
continuation floor by

```text
underlineDelta(S)=inf_(z in M) Delta_z(S).           (4)
```

The infimum is finite: simultaneous translation gives
`0<=rho_z(x)<=ell(x)`.

> **Persistent-game theorem.** `underlineDelta` is a normalized,
> nonnegative, superadditive set function. In particular it is monotone and
> its zero sets form an abstract simplicial complex.

Indeed, for disjoint `A,B`, every context satisfies

```text
Delta_z(A union B)>=Delta_z(A)+Delta_z(B)
                  >=underlineDelta(A)+underlineDelta(B).       (5)
```

Taking the infimum on the left gives

```text
underlineDelta(A union B)
 >=underlineDelta(A)+underlineDelta(B).             (6)
```

Normalization and nonnegativity are inherited from every `Delta_z`.
Equation (6), applied to `B=T\A`, proves monotonicity under `A subset T`.

Let

```text
K_floor={S:underlineDelta(S)=0},

K_some={S:Delta_z(S)=0 for at least one z}.          (7)
```

Always

```text
K_some subset K_floor.                              (8)
```

If `d` is integer-valued, every set of values
`{Delta_z(S):z in M}` is a nonempty subset of the nonnegative integers, so
its infimum is attained. Therefore

```text
K_floor=K_some                  in every integer metric.        (9)
```

For the Gordian metric, (4) is thus a single integer-valued potential for
THM-2259's some-context complex. Its Hasse differences are nonnegative and
curl-free, but they need not be the contextwise infima of THM-2259's Hasse
edge weights: different faces may have different minimizing contexts.

### The integer and calibration hypotheses are real

The equality in (9) can fail for a real metric even on a cancellative
monoid. Take `M=Z_(>=0)` and

```text
d(m,n)=|2^(-m)-2^(-n)|.                             (10)
```

THM-2191 verifies joint nonexpansivity for this metric. On the labelled
packet `(1,1)`,

```text
rho_z(1)=2^(-z-1),
rho_z(2)=3*2^(-z-2),
Delta_z({1,2})=2^(-z-2)>0.                          (11)
```

Thus the pair lies outside `K_some`, but its continuation floor is zero.
Its root defect is `1/4`, so the root also need not minimize the context
defect without the rigidity hypothesis in Section 3.

## 2. Hyperedge packing and the pair-matching shadow

Let `S subset [r]`, and let `A_1,...,A_m` be pairwise disjoint nonempty
subsets of `S`. Iterating (6), and then using monotonicity for their union,
gives

```text
underlineDelta(S)
 >=sum_(j=1)^m underlineDelta(A_j).                 (12)
```

More generally, suppose a finite hyperedge ledger assigns weights

```text
w(A)<=underlineDelta(A),             w(A)>=0.       (13)
```

Then every disjoint hyperedge packing `P` in `S` certifies

```text
Delta_z(S)
 >=underlineDelta(S)
 >=sum_(A in P)w(A)                    for every z. (14)
```

Maximizing the final sum gives a weighted set-packing lower bound.

The two-body specialization is exact graph language. Give the complete
graph on `[r]` the symmetric edge weights

```text
w_ij=underlineDelta({i,j}).                         (15)
```

For every matching `Q` in the induced graph on `S`,

```text
Delta_z(S)>=sum_({i,j} in Q)w_ij,                   (16)
```

and hence the maximum-weight matching is the strongest certificate obtained
from this disjoint-pair packing rule.

This is a legitimate role for a binary graph, but it is one-sided. The two
translation-invariant word metrics in THM-2259 have the same context-indexed
weighted pair graph, with edge weights `1,0,1`, while their triple defects
are respectively `1` and `3`. Their continuation floors equal those
constant defects. Thus the matching bound is exact in the screening model
and strict in the reinforcement model. A pair graph does not classify the
total persistent game, and orienting its symmetric weights as a tournament
cannot restore the missing higher dividend.

## 3. Context-rigid atoms and the exact catalytic spectrum

Before imposing rigidity, there is an exact alignment ledger. Put

```text
a_z(x)=rho_z(x)-ell_cat(x)>=0,                       (16a)

beta_cat(S)
 =sum_(i in S)ell_cat(x_i)-ell_cat(x_S).             (16b)
```

Then direct substitution gives

```text
Delta_z(S)
 =beta_cat(S)
  +sum_(i in S)a_z(x_i)-a_z(x_S),                   (16c)

underlineDelta(S)
 =beta_cat(S)
  +inf_z [sum_(i in S)a_z(x_i)-a_z(x_S)].           (16d)
```

The bracket in (16d) is the common-context alignment charge. It records
whether the atomic and merged diagonal minima occur in compatible contexts;
it has no asserted sign. This is the exact coordinate lost when one
separately minimizes every diagonal response.

Call an element `x` **context-rigid** if

```text
rho_z(x)=ell(x)                         for every z in M.        (17)
```

An additive real-valued `d`-1-Lipschitz invariant which calibrates `x`
implies (17), by THM-2176. The calibrating invariant may depend on the
packet entry.

Define the directional catalytic saving and its capacity by

```text
C_z(x)=ell(x)-rho_z(x)>=0,

kappa(x)=sup_z C_z(x)
        =ell(x)-ell_cat(x),                         (18)
```

where `ell_cat(x)=inf_z rho_z(x)` is THM-2191's catalytic root length.

Suppose every atom `x_i` in the packet is context-rigid. Then for every
subset `S` and every context `z`,

```text
Delta_z(S)
 =sum_(i in S)ell(x_i)-rho_z(x_S)

 =Delta_0(S)+C_z(x_S).                              (19)
```

This is the **calibrated spectrum law**. It has four exact consequences.

First, the root is the minimum context:

```text
underlineDelta(S)=Delta_0(S).                       (20)
```

Second, the complete value set is a translate of the catalytic-saving
spectrum of the merged composite:

```text
{Delta_z(S):z in M}
 =Delta_0(S)+{C_z(x_S):z in M}.                     (21)
```

Therefore

```text
inf_z Delta_z(S)=Delta_0(S),

sup_z Delta_z(S)=Delta_0(S)+kappa(x_S).              (22)
```

When `d` is integer-valued, the supremum is a maximum because
`ell_cat(x_S)` is attained.

Third, since context-rigidity gives `ell_cat(x_i)=ell(x_i)`, the localized
interaction defect is exactly the upper endpoint:

```text
beta_cat(S)
 :=sum_(i in S)ell_cat(x_i)-ell_cat(x_S)

 =Delta_0(S)+kappa(x_S)
 =sup_z Delta_z(S).                                 (23)
```

Thus the continuation floor and the catalytic defect are respectively the
lower and upper endpoints of one spectrum. Separately minimizing every
term is not valid in general; (23) holds because every atomic diagonal is
already rigid.

There is also a useful poison-cone lift. Let `<=_0` be THM-2176's rooted
geodesic order. If

```text
a<=_0 x,              b<=_0 y,                      (23a)
```

and `x,y` are context-rigid, then THM-2176's defect monotonicity and (20)
give

```text
underlineDelta({x,y})
 =sigma(x,y)
 >=sigma(a,b).                                      (23b)
```

Thus a rooted descendant edge becomes persistent as soon as the two
descendants have their own calibrations. The calibration is the sidecar
which upgrades a root poison cone to an all-context poison cone.

Fourth, the context dependence of THM-2259's Hasse field and signed
dividends is exactly catalytic inclusion-exclusion. With
`C_z(x_emptyset)=0`,

```text
c_z(i|A)-c_0(i|A)
 =C_z(x_(A union {i}))-C_z(x_A),                    (24)

h_z(A)-h_0(A)
 =sum_(B subset A)(-1)^(|A|-|B|)C_z(x_B).           (25)
```

No sign is asserted for the difference in (24), even though both Hasse
fields are nonnegative. A context can redistribute shortcut mass among
partial composites; its total effect on a face is the nonnegative scalar in
(19).

## 4. Replicated packets have one context-uniform stable rate

Continue under the context-rigidity hypothesis. Put

```text
A=sum_(i=1)^r ell(x_i),
X=x_1+...+x_r.                                      (26)
```

Take `n` labelled copies of the whole packet and let
`Delta_z^[n]` denote the defect of all `nr` occurrences. Equation (19)
becomes

```text
Delta_z^[n]
 =nA-rho_z(nX)
 =[nA-ell(nX)]+C_z(nX).                             (27)
```

Consequently

```text
inf_z Delta_z^[n]=nA-ell(nX),

sup_z Delta_z^[n]=nA-ell_cat(nX),                   (28)
```

and the spectral width is

```text
kappa(nX)=ell(nX)-ell_cat(nX).                      (29)
```

THM-2191 gives `ell(nX)/n -> ell_hash(X)`, while THM-2220 gives

```text
kappa(nX)/n -> 0.                                   (30)
```

Therefore the whole normalized spectrum collapses:

```text
sup_(z in M)
 |Delta_z^[n]/n-[A-ell_hash(X)]|
 ->0.                                               (31)
```

In particular (31) remains true when the context `z=z_n` changes with `n`.
Finite contexts can expose catalytic variation, but no moving sequence of
contexts changes the stable interaction rate

```text
A-ell_hash(X).                                      (32)
```

### Stable duality is a common-calibrator theorem

Specialize now to knots, so `ell=u`, and let `L` be the class of additive
real-valued Gordian-1-Lipschitz invariants. THM-2191's rational
Hahn--Banach theorem says that the supremum below is attained:

```text
u_hash(X)
 =max_(phi in L)|phi(X)|
 =sup_(phi in L)|sum_(i=1)^r phi(x_i)|.             (32a)
```

Consequently the stable replicated defect rate in (32) has the exact dual
form

```text
sum_(i=1)^r u(x_i)-u_hash(X)

 =sum_(i=1)^r u(x_i)
  -sup_(phi in L)|sum_(i=1)^r phi(x_i)|.            (32b)
```

If every atom is individually additively calibrated, then THM-2191 also
gives

```text
u(x_i)=max_(phi in L)|phi(x_i)|,
```

and (32b) becomes

```text
sum_i max_(phi in L)|phi(x_i)|
 -max_(phi in L)|sum_i phi(x_i)|.                   (32c)
```

This is exactly a **common-calibrator incompatibility**. The first term may
choose a different functional and a different absolute-value sign for every
atom. The second must choose one functional and one final sign for the
whole packet.

The equality case is precise. The rate in (32c) vanishes if and only if
there are `phi in L` and one sign `epsilon in {+1,-1}` such that

```text
epsilon phi(x_i)=u(x_i)                 for every i. (32d)
```

Indeed,

```text
|sum_i phi(x_i)|<=sum_i|phi(x_i)|<=sum_i u(x_i).    (32e)
```

Equality throughout forces every atom to be calibrated by the same `phi`
and all nonzero values `phi(x_i)` to have one sign; the converse is
immediate. Thus separate calibrators are not enough. Stable additivity asks
for one coherently signed common calibrator.

## 5. The Brittenham--Hermiller mirror pair

Let

```text
K=T(2,7),                 Kbar=mirror(K),
X=K#Kbar,
q_X=u(X),                 g=6-q_X.                  (33)
```

Brittenham--Hermiller prove

```text
u(K)=u(Kbar)=3,             q_X<=5,                 (34)
```

so `g>=1`. Half-signature is additive and Gordian-1-Lipschitz, and it
calibrates `K` and `Kbar` with opposite signs. Hence, for every knot `J`,

```text
d_G(K#J,J)=d_G(Kbar#J,J)=3.                         (35)
```

Both atoms are context-rigid. The calibrated spectrum law gives the exact
pair identity

```text
Delta_J({K,Kbar})
 =6-d_G(X#J,J)
 =g+C_J(X).                                         (36)
```

Thus

```text
minimum pair defect:       6-u(X)=g>=1,

maximum pair defect:       6-u_cat(X),

pair-spectrum span:        u(X)-u_cat(X)=kappa(X).  (37)
```

A context increases the pair defect beyond its root value if and only if it
is a positive Gordian catalyst for the composite `X`. This statement does
not assert that such a context exists. It converts that search into the
question whether the persistent Brittenham--Hermiller edge can be
strictly amplified.

Now take `p` labelled copies of `K` and `q` labelled copies of `Kbar`.
Every cross-chirality pair has continuation floor `g`, so the persistent
pair graph contains a weighted `K_(p,q)`. Its maximum matching has size
`min(p,q)`. Equations (16) and (35) yield, for every knot `J`,

```text
Delta_J(pK,qKbar)>=g min(p,q)>=min(p,q),             (38)

d_G(pK#qKbar#J,J)
 <=3(p+q)-g min(p,q)
 <=3(p+q)-min(p,q).                                 (39)
```

For `p=q=n`, this is the context-uniform amplification

```text
d_G(nX#J,J)<=n q_X<=5n                 for every J. (40)
```

The same matching argument propagates through calibrated Gordian
descendants. Suppose

```text
K<=_0 L_i                    for 1<=i<=p,
Kbar<=_0 R_j                 for 1<=j<=q,           (40a)
```

and every `L_i,R_j` is context-rigid. Equation (23b) gives persistent
cross-edge weight at least `g`, so for every `J`,

```text
d_G(
  L_1#...#L_p#R_1#...#R_q#J,
  J
)
 <=sum_i u(L_i)+sum_j u(R_j)-g min(p,q).             (40b)
```

This applies in particular to the calibrated positive-torus/mirrored
Gordian-adjacent families audited in THM-2176. It makes no new adjacency or
unknotting-number claim about those descendants.

The stable refinement is structurally sharper than merely repeating the
five-change certificate. For every sequence of context knots `J_n`,

```text
[6n-d_G(nX#J_n,J_n)]/n
 ->6-u_hash(X),                                     (41)
```

uniformly in the choice of `J_n`. The exact values of `q_X`,
`u_cat(X)`, and `u_hash(X)` remain open within the bounds used here.

The common-calibrator form is

```text
6-u_hash(X)
 =6-max_(phi in L)|phi(K)+phi(Kbar)|
 >=1.                                               (41a)
```

Half-signature supplies separate calibrations with values `3` and `-3`, so
its two atomic calibrations have opposite signs. The inequality in (41a),
which follows from `u_hash(X)<=u(X)<=5`, proves that no other additive
Gordian-1-Lipschitz invariant coherently calibrates both atoms either. The
positive stable defect is exactly the cost of that common-calibration
failure.

## 6. Relation ledger and honest frontier

The exact connection to graphs, tournaments, and prime decomposition is:

```text
source:
  the diagonal continuation stalk z |-> Delta_z;

compression:
  take its pointwise infimum on every Boolean face;

preserved predicate:
  in integer metrics, whether some context kills the face defect;

binary shadow:
  persistent pair weights, used safely through disjoint matching;

destroyed information:
  higher reinforcement/screening and the identity of minimizing contexts;

sidecar:
  the persistent Hasse field or full context spectrum, and ultimately the
  off-diagonal min-plus continuation kernel;

cheapest knot test:
  first determine q_X=u(K#Kbar), then seek a knot J with
  d_G(X#J,J)<q_X.                                    (42)
```

Schubert prime decomposition keeps `K` and `Kbar` as rigid prime atoms in
THM-2242's free commutative kernel image. The bipartite matching in (38) is
therefore an occurrence-level certificate, not an identification of the two
prime types. Its weights are symmetric and permit zeros in general; there
is no intrinsic tournament orientation. The Boolean Hasse DAG remains the
exact finite directed carrier for the chosen subset-defect table.

No unknown unknotting number, positive knot catalyst, or higher-order knot
minimal nonface is proved here. The new content is the exact separation

```text
persistent floor  = unavoidable interaction;
spectrum width    = catalytic capacity of the merged composite;
stable rate       = context-independent homogenized interaction;
pair matching     = valid lower certificate, not a classifier.             (43)
```

QED.
