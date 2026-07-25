---
id: THM-2292
title: "Common catalytic section, rank-sharp calibration nerve, and signed gain code"
status: >
  PROVED. For every finite labelled packet in a commutative nonexpansive
  integer-valued metric monoid, one context simultaneously attains the
  catalytic length of every subset sum. Hence the complete Boolean defect
  table of the catalytic length is one honest diagonal continuation slice,
  not a coordinatewise splice; the persistent continuation floor is bounded
  above by this realized catalytic defect. In a finite-dimensional
  seminormed space, the zero-defect complex is exactly the nerve of the
  atoms' exposed calibration faces. More sharply, every minimal nonface of
  arity at least three is linearly independent in the effective quotient;
  its arity is therefore at most the effective dimension, with two-body
  nonfaces the only exceptional case. The sharp uniform bound is
  max(2,d), and stable pair data are complete for zero-versus-positive
  interaction in effective rank at most two. Signed calibration faces give
  an antipodal code whose unique pair signs form a switched Z/2 gain graph.
  Triangle balance is exact in rank at most two but is insufficient from
  rank three onward without the signed nerve sidecar. The
  Brittenham--Hermiller T(2,7)-mirror packet spans an exact stable rank-two
  plane, has a two-word chirality calibration code, and obeys a quantitative
  code-distance/matching defect bound. Conversely, the l_infinity boundary
  simplex attains the rank bound and gives the first balanced-pair failure
  in rank three. No unknown ordinary or homogenized unknotting number is
  computed, and no positive Gordian catalyst is produced.
source: codex-2026-07-25-catalytic-calibration-nerve
depends_on:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2259-boolean-continuation-hasse-field-and-signed-interaction-dividends
  - THM-2272-persistent-interaction-packing-and-calibrated-catalytic-defect-spectrum
  - THM-2281-common-optimal-context-for-finite-catalytic-families
related:
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
  - "Zhiqing Yang, Unknotting number of the connected sum of n identical knots, Journal of Knot Theory and Its Ramifications 17 (2008), 253--255, DOI 10.1142/S0218216508006130."
---

# THM-2292 -- catalytic sections and rank-sharp calibration codes

There are two different compatibility questions behind a packet of
composable objects:

```text
context compatibility:
  can one continuation attain all requested catalytic minima?

calibrator compatibility:
  can one additive dual functional certify all requested stable equalities?
```

The first has finite Helly number one because optimal-context loci are
additive upper ideals. The second is a convex-intersection problem whose
effective dimension controls the first possible missing face. Confusing the
two is precisely how a separately minimized scalar table can be mistaken
for a compositional state.

## 1. One context realizes the whole catalytic Boolean table

Let `(M,+,0)` be a commutative monoid with an integer-valued metric `d`
satisfying joint nonexpansivity

```text
d(a+c,b+e)<=d(a,b)+d(c,e).                           (1)
```

Put

```text
ell(x)=d(x,0),
rho_z(x)=d(x+z,z),
ell_cat(x)=min_z rho_z(x).                           (2)
```

Fix a labelled packet `x_1,...,x_r`, and for `S subset [r]` write

```text
x_S=sum_(i in S)x_i.                                 (3)
```

Apply THM-2281 to the finite family of all `2^r` subset sums `x_S`.
There is one context `z_*` such that, simultaneously for every subset,

```text
rho_(z_*)(x_S)=ell_cat(x_S).                         (4)
```

Define the catalytic Boolean defect

```text
Beta_cat(S)
 =sum_(i in S)ell_cat(x_i)-ell_cat(x_S).             (5)
```

At the single context in (4), THM-2259's diagonal defect is

```text
Delta_(z_*)(S)
 =sum_(i in S)rho_(z_*)(x_i)-rho_(z_*)(x_S)
 =Beta_cat(S)                         for every S.   (6)
```

Thus the entire table `(Beta_cat(S))_(S subset [r])` is one actual stalk
of the continuation field. It is not assembled from incompatible
face-dependent contexts.

Several consequences are immediate.

1. `Beta_cat` is normalized, nonnegative, monotone, and superadditive on
   disjoint subsets, because the same is true of the one context table
   `Delta_(z_*)`.
2. Its zero sets form an abstract simplicial complex.
3. Its Hasse increments and signed Möbius dividends are exactly the Hasse
   field and dividends of one context:

   ```text
   c_cat(i|A)
    =ell_cat(x_i)+ell_cat(x_A)-ell_cat(x_(A union {i}))
    =c_(z_*)(i|A),                                  (7)

   h_cat(A)
    =sum_(B subset A)(-1)^(|A|-|B|)Beta_cat(B)
    =h_(z_*)(A).                                    (8)
   ```

4. If `underlineDelta` is THM-2272's persistent continuation floor, then

   ```text
   underlineDelta(S)
    =min_z Delta_z(S)
    <=Beta_cat(S).                                  (9)
   ```

Equation (9) is only an upper bound. The context which minimizes all
individual and merged diagonal responses need not minimize their
difference. No equality between persistent and catalytic defect is claimed.

For knots, (4) says that every finite packet has one knot `J_*` satisfying

```text
d_G(K_S#J_*,J_*)=u_cat(K_S)             for every S. (10)
```

Hence all catalytic subset defects, Hasse increments, and dividends of the
packet occur simultaneously at one genuine Gordian context.

## 2. Stable additivity is a calibration-face nerve

Let `(V,p)` be a finite-dimensional real seminormed vector space. Quotient
by `ker(p)` when necessary. Fix labelled vectors `x_1,...,x_r`, let

```text
Y=span_R{x_1,...,x_r},             d=dim(Y/ker(p|Y)), (11)
```

and let `B_*` be the unit ball of the dual of the resulting normed quotient.
For each atom define its exposed calibration face

```text
F_i={phi in B_*:phi(x_i)=p(x_i)}.                   (12)
```

Hahn--Banach makes every `F_i` nonempty; if `p(x_i)=0`, the condition is
automatic on the quotient.

For `S subset [r]`, put

```text
Delta_p(S)=sum_(i in S)p(x_i)-p(x_S).                (13)
```

> **Calibration-nerve theorem.**
>
> ```text
> Delta_p(S)=0
>   iff
> intersection_(i in S)F_i is nonempty.             (14)
> ```

If `phi` lies in the intersection, then

```text
p(x_S)>=phi(x_S)=sum_(i in S)p(x_i)>=p(x_S),        (15)
```

so equality holds. Conversely, suppose `Delta_p(S)=0`. Choose
`phi in B_*` calibrating `x_S`; it may be chosen of norm one when
`p(x_S)>0`, while `phi=0` is allowed when the calibrated sum vanishes. Then

```text
sum_i p(x_i)
 =p(x_S)
 =sum_i phi(x_i)
 <=sum_i p(x_i).                                    (16)
```

Every summand in the final inequality is termwise bounded, so equality of
the sums forces `phi(x_i)=p(x_i)` for every `i in S`. This proves (14).

Consequently the stable zero-defect complex is literally the nerve of the
finite family of compact convex faces `(F_i)`.

## 3. The rank-sharp calibration gate

The faces in (12) live in the `d`-dimensional real vector space dual to the
effective span in (11). Helly's theorem and (14) first give the baseline:
if `Delta_p(A)>0`, some

```text
B subset A,              |B|<=d+1                  (17)
```

still has positive defect. The special fact that the sets are calibration
faces of a symmetric dual ball sharpens this by one in every genuinely
higher-order case.

> **Rank-sharp minimal-nonface theorem.** Let `A` be a minimal nonface of
> the stable zero-defect complex. If `|A|>=3`, then the images of
> `(x_i)_(i in A)` in `Y/ker(p|Y)` are linearly independent. Consequently
>
> ```text
> |A|<=max(2,d).                                    (17a)
> ```
>
> In particular, in effective dimension at most two the zero-defect
> complex is flag: every positive defect contains a positive pair.

### Proof of linear independence

Relabel `A=[q]`, where `q>=3`. No `p(x_i)` can vanish: deleting such an
atom leaves both the sum and the defect unchanged, contrary to minimality.
For every `i`, the proper face `A\{i}` has zero defect. By (14), choose
`phi_i in B_*` satisfying

```text
phi_i(x_j)=p(x_j)                    for every j!=i. (17b)
```

The full intersection is empty, so

```text
s_i:=p(x_i)-phi_i(x_i)>0.                           (17c)
```

Suppose there is a relation in the effective quotient,

```text
sum_i a_i x_i=0.                                    (17d)
```

Put `P=sum_i a_i p(x_i)`. Applying `phi_i` to (17d) and using (17b)
gives, for every `i`,

```text
P=a_i s_i.                                          (17e)
```

If some `a_i` vanished, (17e) would give `P=0`, and then the strict
positivity of every `s_j` would force every `a_j=0`. Thus a nonzero
relation has no zero coefficient. Equation (17e) also forces all
coefficients to have one sign. Reverse the relation if necessary and
assume `a_i>0`.

Set

```text
b_i=a_i p(x_i)>0,                 B=sum_i b_i.       (17f)
```

Applying `phi_i` once more gives

```text
B-b_i=-a_i phi_i(x_i)
      <=a_i p(x_i)=b_i.                             (17g)
```

Summing (17g) over `i` yields

```text
(q-1)B<=B,                                          (17h)
```

which is impossible for `q>=3`. Hence the packet is linearly independent.
The bound (17a) follows. Two-body nonfaces really are exceptional: the
opposite vectors `1,-1` in the absolute-value norm give a dependent
minimal nonface in effective dimension one.

Thus the exact stable witness boundary is

```text
effective rank 0: no positive defect;
effective rank 1 or 2: pair defects are complete;
effective rank d>=3: a first interaction has arity at most d.  (17i)
```

No statement about recovering the **magnitude** of a larger defect from
its small faces is asserted.

For the Gordian metric, THM-2191 puts the homogenized length `u_hash` on
the rational Grothendieck group and identifies its dual unit ball with
additive Gordian-1-Lipschitz calibrators. Restricting to the finite rational
span of a knot packet and extending scalars to the corresponding real
finite-dimensional quotient gives (14)--(17i) for

```text
Delta_hash(S)
 =sum_(i in S)u_hash(K_i)-u_hash(K_S).               (18)
```

If every atom is individually stable-tight,
`u_hash(K_i)=u(K_i)`, then (18) is exactly THM-2272's common-calibrator
stable interaction rate on every subset. A genuinely `q`-body stable knot
interaction with `q>=3` and every proper subpacket additive would therefore
force effective stable rank at least `q`.

## 4. The signed calibration code and its gain-graph shadow

The unsigned nerve tests whether all atoms can be calibrated with one
common sign. Retaining every possible relative calibration sign gives an
exact antipodal refinement. For a sign word
`epsilon in {+1,-1}^S`, define

```text
Delta_p^epsilon(S)
 =sum_(i in S)p(x_i)-p(sum_(i in S)epsilon_i x_i),  (S1)

F_i^epsilon
 ={phi in B_*:epsilon phi(x_i)=p(x_i)}.             (S2)
```

Applying (14) to the signed atoms `epsilon_i x_i` gives

```text
Delta_p^epsilon(S)=0
 iff
 intersection_(i in S)F_i^(epsilon_i) is nonempty. (S3)
```

The feasible signed subsets form an antipodal simplicial complex: negating
every sign is realized by `phi |-> -phi`. Its full-face words form the
**signed calibration code**

```text
C_p(x_1,...,x_r)
 ={epsilon in {+1,-1}^r:Delta_p^epsilon([r])=0}.    (S4)
```

The pair shadow is naturally four-valued. For `i!=j`, put

```text
G_ij
 ={s in {+1,-1}:
   p(x_i+s x_j)=p(x_i)+p(x_j)}.                     (S5)
```

The set `G_ij` may be empty, a singleton, or both signs. Those are honest
obstruction, unique-gauge, and tie states; none should be forced into a
tournament orientation.

Suppose now that every `G_ij` is a singleton `{g_ij}`. The labels form a
`Z/2` gain graph. Replacing `x_i` by `eta_i x_i` switches the incident
labels by

```text
g_ij |-> eta_i eta_j g_ij.                          (S6)
```

A sign word is pairwise feasible exactly when

```text
epsilon_i epsilon_j=g_ij              for all i,j. (S7)
```

Such a word exists if and only if every triangle has trivial holonomy,

```text
g_ij g_jk g_ki=+1.                                  (S8)
```

Indeed, necessity follows by multiplying (S7). Conversely, fix
`epsilon_1=1` and set `epsilon_j=g_(1j)`; the triangles through vertex
one give (S7). The pairwise word is unique up to global sign.

Pairwise balance is not generally sufficient for membership in the full
code (S4). The rank-sharp theorem gives its exact positive boundary:

```text
d<=2:
  triangle balance iff C_p is the two-word code {epsilon,-epsilon};

d>=3:
  triangle balance is necessary, while signed minimal nonfaces of
  arity at most d are the missing sidecar.                          (S9)
```

There is also a quantitative pair certificate. For a fixed sign word put

```text
w_ij(epsilon)
 =p(x_i)+p(x_j)-p(epsilon_i x_i+epsilon_j x_j).
                                                               (S10)
```

For every matching `Q` on the packet, subadditivity after grouping the
matched pairs gives

```text
Delta_p^epsilon([r])
 >=sum_({i,j} in Q)w_ij(epsilon).                   (S11)
```

Hence the maximum-weight matching is the strongest disjoint-pair
certificate. In the singleton balanced case, let `chi` be the pairwise
word and suppose every wrong-relative-sign pair has defect at least
`lambda>0`. If

```text
D={i:epsilon_i!=chi_i},
```

then the wrong pairs form the complete bipartite cut
`D x ([r]\D)`. A maximum matching in that cut has size

```text
min(|D|,r-|D|)
 =d_H(epsilon,{chi,-chi}).                          (S12)
```

Equations (S11)--(S12) prove the code-distance bound

```text
Delta_p^epsilon([r])
 >=lambda d_H(epsilon,{chi,-chi}).                  (S13)
```

This relation has a precise tournament-analysis ledger:

```text
vertices:             labelled stable atoms;
pair observable:      allowed relative calibration signs G_ij;
gauge:                independent sign switching at vertices;
preserved target:     pairwise stable additivity;
lost data:            higher calibration-face intersections and magnitudes;
sidecar:              the signed calibration nerve through effective rank;
intrinsic direction:  none.                         (S14)
```

## 5. Boundary-simplex defects occur in honest normed spaces

The Helly bound does not make binary data sufficient. Fix `m>=2`, take

```text
V=R^m,                   p(x)=||x||_infinity,

mathbf1=(1,...,1),       x_i=mathbf1-e_i,
                                      1<=i<=m.       (19)
```

Every atom has norm one. If `S` is proper and `s=|S|`, then

```text
sum_(i in S)x_i
 =s*mathbf1-sum_(i in S)e_i.                        (20)
```

The coordinates indexed by `S` equal `s-1`, while at least one coordinate
outside `S` equals `s`. Therefore

```text
p(x_S)=s=sum_(i in S)p(x_i)             if S proper. (21)
```

For the full packet,

```text
x_[m]=(m-1)mathbf1,

p(x_[m])=m-1,
Delta_p([m])=1.                                    (22)
```

Thus its zero-defect complex is exactly the boundary of the
`(m-1)`-simplex. Every proper subpacket is additive, while the full packet
has a one-unit saving. For `m>=3`, this attains the rank-sharp bound:
the unique minimal nonface has arity exactly the effective dimension `m`.

The dual unit ball is the `l_1` ball. The calibration faces make the
mechanism visible:

```text
F_i=conv{e_j^*:j!=i},

intersection_(i in S)F_i
 =conv{e_j^*:j notin S}                 for S proper,

intersection_(i=1)^m F_i=empty,                     (23)
```

where the middle identity is for nonempty proper `S` (the intersection of
the empty family is the whole dual unit ball).

Give `V` its translation-invariant metric

```text
d(x,y)=p(x-y).                                      (24)
```

Then for every context `z`,

```text
rho_z(x)=p(x),

ell_cat(x)=ell_hash(x)=p(x).                        (25)
```

The missing full face in (22) is therefore simultaneously:

```text
a root interaction;
a persistent all-context interaction;
a catalytic interaction;
a homogenized stable interaction.                  (26)
```

There is a literal same-shadow comparator. Keep the labelled packet size
`m` but put

```text
y_i=e_1                         for every i.         (27)
```

Then every atom again has norm one and

```text
p(y_S)=|S|                      for every S,         (28)
```

so every defect of the `y` packet is zero. The `x` and `y` packets
therefore have identical defect tables on every proper subset, but their
full defects are respectively one and zero.

Choosing `m>k` proves that no fixed `k`-ary truncation determines stable
interaction in general. In particular, for `m>=3` the two packets have the
same empty pair-defect graph. No tournament or binary relation constructed
only from that pair shadow can recover which packet has the unique
`m`-body defect.

The signed pair shadow also collides. For distinct `i,j` and `m>=3`,

```text
p(x_i+x_j)=2,                  p(x_i-x_j)=1,
p(y_i+y_j)=2,                  p(y_i-y_j)=0.         (28a)
```

Thus both packets have the same unique gain label

```text
G_ij={+1}                     for every pair.        (28b)
```

Every triangle holonomy is trivial. Nevertheless the `x` packet has empty
full signed code: its unique pairwise candidates are the two constant sign
words, and (22) excludes both. The `y` packet has exactly those two
constant words in its code. At `m=3`, this is the first possible
rank boundary:

```text
rank at most two: balanced pair gain data are complete;
rank three:       balanced pair gain data can miss a triple obstruction.
                                                               (28c)
```

This strengthens the abstract stopping boundary in THM-2248 from an
integer word metric to a genuine normed and already homogenized object. It
does not realize the boundary simplex by knots.

## 6. The Brittenham--Hermiller stable plane and its code

Let

```text
K=T(2,7),                   Kbar=mirror(K),
a=[K],                     b=[Kbar]                 (BH1)
```

in the rational stable knot group of THM-2191, and write `p=u_hash` for
its real seminorm extension. Orient half-signature so that

```text
(sigma/2)(a)=3,             (sigma/2)(b)=-3.         (BH2)
```

Half-signature is additive and Gordian-1-Lipschitz. It therefore gives

```text
p(a)=p(b)=3,
p(a-b)=6.                                          (BH3)
```

The first two equalities also use the obvious upper bounds from
`u(K)=u(Kbar)=3`; the last uses the triangle upper bound
`p(a-b)<=p(a)+p(b)=6`.

Brittenham--Hermiller exhibit at most five crossing changes for
`K#Kbar`, so

```text
p(a+b)<=u(K#Kbar)<=5.                              (BH4)
```

No equality in (BH4) is known or used.

The stable span is genuinely two-dimensional. Yang's published theorem
states that if a knot `J` has nontrivial Alexander polynomial, then

```text
u(#^n J)>=n                         for every n>=1. (BH5)
```

The Alexander polynomial of `X=K#Kbar` is the product of two nontrivial
mirror Alexander polynomials and is nontrivial. Applying (BH5) to `X` and
passing to the stable limit gives

```text
p(a+b)=u_hash(X)>=1.                               (BH6)
```

Put `e=a+b` and `o=a-b`. If `alpha e+beta o` lies in `ker(p)`, then
half-signature gives `6 beta=0`; after `beta=0`, (BH6) gives
`alpha=0`. Thus

```text
dim_R(span{a,b}/ker p)=2.                           (BH7)
```

The rank-two flag theorem now says that pair data are complete for
zero-versus-positive stable interaction throughout this plane. For a
labelled packet made of occurrences of `a` and `b`, its unique pair gains
are explicit:

```text
two a-occurrences:     G_ij={+1},
two b-occurrences:     G_ij={+1},
one of each type:      G_ij={-1}.                   (BH8)
```

Indeed, positive same-type sums have norm six by homogeneity while their
negative differences vanish. The cross-type negative sum has norm six by
(BH3), while the positive sum has norm at most five by (BH4).

Let `chi` be `+1` on every `a` occurrence and `-1` on every `b`
occurrence. The gain graph in (BH8) is balanced and has pairwise words
`{chi,-chi}`. Both are actual full calibration words because
half-signature calibrates every signed atom in the `chi` word and its
negative calibrates `-chi`. Conversely every full word restricts to every
pair, so (BH8) forces one of these two patterns. Hence

```text
C_p(packet)={chi,-chi}.                             (BH9)
```

This is an exact stable calibration code, not a claim that formal negative
coefficients are represented by physical connected sums.

There is a quantitative consequence. Put

```text
delta=6-p(a+b)=6-u_hash(K#Kbar).                    (BH10)
```

Equations (BH4)--(BH6) give

```text
1<=delta<=5.                                        (BH11)
```

A wrong relative sign on a same-type pair has defect six; a wrong relative
sign on a cross-type pair has defect `delta`. Therefore every wrong pair
has defect at least `delta`, and (S13) gives, for a packet of `P`
occurrences of `a` and `Q` occurrences of `b` and every coefficient-sign
word `epsilon`,

```text
3(P+Q)-p(sum_i epsilon_i x_i)
 >=delta d_H(epsilon,{chi,-chi}).                   (BH12)
```

For the physical all-positive connected sum, the distance from the
all-positive word to `{chi,-chi}` is `min(P,Q)`. Thus

```text
u_hash((#^P K)#(#^Q Kbar))
 <=3(P+Q)
   -(6-u_hash(K#Kbar))min(P,Q)
 <=3(P+Q)-min(P,Q).                                 (BH13)
```

The first inequality can also be seen by grouping
`min(P,Q)` copies of `K#Kbar`; the gain-code proof identifies the
coefficient as the exact wrong-pair stable defect and the minority count as
the Hamming distance to coherent calibration.

Equations (BH7)--(BH13) refine the operation-ready description of this
counterexample but compute neither `u(K#Kbar)` nor `u_hash(K#Kbar)`.
They produce no positive Gordian catalyst and do not promote a stable
inequality to an ordinary unknotting lower bound.

## 7. Information ledger and scope

The three finite carriers now separate cleanly:

```text
catalytic contexts:
  additive upper ideals;
  every finite family has one common section;

stable calibrators:
  exposed convex faces;
  their nerve is exact and obeys the rank-sharp max(2,d) gate;

signed pair shadow:
  a switched gain graph, exact in effective rank at most two;
  incomplete from rank three without the signed nerve sidecar;

tournament shadow:
  no intrinsic orientation; forcing one loses ties and switching gauge.
                                                               (29)
```

The source-to-target maps are:

```text
all subset sums x_S
  -> one common optimal context z_*
  -> one realized catalytic Boolean table;

stable atoms x_i
  -> calibration faces F_i
  -> their intersection nerve
  -> zero stable-defect complex.                    (30)
```

The first map loses the size and identity of a smallest context. The second
loses the full off-diagonal continuation kernel and retains only stable
additivity. The gain graph further loses higher face intersections and
defect magnitudes. The rank-sharp gate and the Brittenham--Hermiller plane
recover exactly the zero-versus-positive target in effective rank at most
two; they do not classify raw unknotting number. No positive knot catalyst,
unknown unknotting number, exact value of `u_hash(K#Kbar)`, or higher-order
Gordian minimal nonface is produced. QED.
