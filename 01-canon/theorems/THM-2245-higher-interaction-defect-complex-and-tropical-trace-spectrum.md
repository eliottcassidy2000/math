---
id: THM-2245
title: "Higher interaction defect complex and tropical trace spectrum"
status: >
  PROVED. For every finite labelled family in a commutative nonexpansive
  metric monoid, the subsets on which root length is additive form an
  abstract simplicial complex. Its minimal nonfaces are irreducible
  multiway interaction certificates, and they can have arbitrary arity even
  in translation-invariant integer word metrics; hence the binary
  interaction graph, and a fortiori any pairwise orientation of its
  one-skeleton, is not a complete compositional carrier. For the faithful min-plus
  continuation kernel, raw root length, catalytic root length, and
  homogenized root length are respectively a pinned diagonal, the infimal
  diagonal, and the asymptotic infimal-diagonal slope. Every fixed diagonal
  has the same asymptotic slope. Applied to knots, this gives a canonical
  higher-interaction complex on any labelled prime-summand packet and a
  tropical trace hierarchy for unknotting number; it does not exhibit a
  higher-order knot symbiont or compute a new unknotting number.
source: klein-2026-07-25-higher-interaction-defect-complex
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2235-response-antipode-barriers-for-lrc-sheets-tournament-cycles-and-knot-kernels
---

# THM-2245 -- the interaction object is a complex, not a graph

THM-2176 replaces unknotting number by its continuation kernel and records
the pair interaction

```text
sigma(x,y)=ell(x)+ell(y)-ell(x+y).
```

The support of `sigma` is useful, but it is still only binary data. This
theorem identifies the exact finite higher-order shadow and shows why no
binary orientation can replace it.

## 1. The higher defect

Let `(M,+,0)` be a commutative monoid with a metric `d` satisfying

```text
d(a+c,b+d)<=d(a,b)+d(c,d),                          (1)
```

and put `ell(x)=d(x,0)`. Fix a finite labelled family

```text
x_1,...,x_r.
```

Labels are retained even when two entries of the family are equal. For
`S subset [r]`, write

```text
x_S=sum_(i in S)x_i,

Delta(S)=sum_(i in S)ell(x_i)-ell(x_S).             (2)
```

The empty sum is zero, so `Delta(emptyset)=0`.

Subadditivity of `ell` gives

```text
Delta(S)>=0.                                        (3)
```

More strongly, if `S subset T` and `W=T\setminus S`, then

```text
Delta(T)-Delta(S)
 =ell(x_S)+sum_(i in W)ell(x_i)-ell(x_S+x_W)

 =Delta(W)+sigma(x_S,x_W)>=0.                       (4)
```

Here `sigma` is THM-2176's nonnegative pair defect. Thus `Delta` is
monotone under inclusion.

## 2. The zero-defect simplicial complex

Define

```text
K_ell(x_1,...,x_r)
 ={S subset [r]:Delta(S)=0}.                        (5)
```

Equation (4) proves:

> **Defect-complex theorem.** `K_ell` is an abstract simplicial complex.
> Equivalently, the positive-defect subsets form an inclusion upper set.

Indeed, if `T` has zero defect and `S subset T`, then

```text
0<=Delta(S)<=Delta(T)=0.
```

Every vertex belongs to `K_ell`, since `Delta({i})=0`. A minimal nonface
`A` therefore has at least two vertices and satisfies

```text
Delta(A)>0,
Delta(B)=0                    for every proper B subset A.      (6)
```

It is an **irreducible multiway interaction**: no smaller subpacket records
any saving. For every `i in A`, equation (4) specializes to

```text
Delta(A)
 =sigma(x_(A\setminus{i}),x_i)>0.                  (7)
```

Thus an arity-`m` minimal nonface is visible as a pair interaction only
after `m-1` objects have already been merged. It need not be visible on any
pair of original vertices.

For an ordering `i_1,...,i_m` of `S`, put
`S_k={i_1,...,i_k}`. Telescoping gives the exact path formula

```text
Delta(S)
 =sum_(k=2)^m sigma(x_(S_(k-1)),x_(i_k)).           (8)
```

The summands depend on the order, but their total does not. THM-2176's
cocycle identity is precisely the elementary adjacent reassociation that
preserves (8).

## 3. Minimal nonfaces have unbounded arity

The higher complex is not determined by its graph. Fix any `r>=3`. Let
`M=Z^r`, let `e_1,...,e_r` be its standard basis, and put

```text
g=e_1+...+e_r.
```

Give the Cayley generators `+/-e_i` cost one and `+/-g` cost `r-1`.
Let `ell` be the resulting weighted word length and

```text
d(a,b)=ell(a-b).                                    (9)
```

This is a translation-invariant integer metric and satisfies (1).

For `S subset [r]`, write `e_S=sum_(i in S)e_i` and `s=|S|`.
If a word representing `e_S` has net `g`-coefficient `k in Z`, the cheapest
coordinate correction has cost

```text
F_s(k)
 =|k|(r-1)+s|1-k|+(r-s)|k|.                        (10)
```

Conversely, the corresponding `g`-letters and coordinate corrections
realize this cost, so

```text
ell(e_S)=min_(k in Z)F_s(k).                        (11)
```

When `0<=s<r`, the minimum is `F_s(0)=s`. For `k=1`,

```text
F_s(1)=2r-s-1>s,
```

and `k<=-1` or `k>=2` is still larger. When `s=r`,

```text
F_r(1)=r-1<F_r(0)=r,
```

and the other integers are larger. Hence

```text
ell(e_S)=|S|                  for every proper S,
ell(e_[r])=r-1.                                    (12)
```

For the labelled family `(e_1,...,e_r)`,

```text
Delta(S)=0                    for every proper S,
Delta([r])=1.                                      (13)
```

Therefore `K_ell` is the boundary of an `(r-1)`-simplex and its unique
minimal nonface has arity `r`.

This proves two sharp limitations:

```text
pairwise zero defect does not imply total zero defect;
no fixed arity truncation determines all interaction defects
in abstract metric monoids.                         (14)
```

The example is a hostile abstract control. It is not asserted to occur in
the Gordian metric.

## 4. Kernel formula for the defect

THM-2176 associates to `x in M` the min-plus continuation kernel

```text
P_x(a,b)=d(x+a,b)
```

and proves

```text
P_x tensor P_y=P_(x+y).                             (15)
```

Consequently the higher defect is exactly the loss at one pinned entry
after composing the full kernels:

```text
Delta(S)
 =sum_(i in S)P_(x_i)(0,0)
  -(tensor_(i in S)P_(x_i))(0,0).                  (16)
```

Equation (16) explains both the success and the limitation of the defect
complex. It detects every shortcut visible at the chosen root, but it
forgets all other input-output entries of the faithful kernel.

## 5. The three min-plus traces

For a kernel `Q` on `M`, define its infimal diagonal

```text
tr_min(Q)=inf_(a in M)Q(a,a).                       (17)
```

For the continuation kernel, THM-2191 gives

```text
P_x(0,0)=ell(x),
tr_min(P_x)=ell_cat(x),                             (18)
```

where `ell_cat` is the common-context localization of the metric. From
(15),

```text
P_x^(tensor n)=P_(nx),

tr_min(P_x^(tensor n))=ell_cat(nx).                 (19)
```

Define the asymptotic min-plus trace slope by

```text
lambda_min(P_x)
 =lim_(n->infinity)
    tr_min(P_x^(tensor n))/n.                       (20)
```

The limit exists, and THM-2191's exact commutation of localization and
homogenization proves

```text
lambda_min(P_x)=ell_hash(x),                        (21)
```

where `ell_hash(x)=lim_n ell(nx)/n`.

THM-2220 supplies the stronger fixed-diagonal statement. For every fixed
context `a`,

```text
lim_(n->infinity) P_x^(tensor n)(a,a)/n
 =ell_hash(x).                                      (22)
```

Thus the three exact levels are

```text
ell(x)       =one pinned diagonal entry;
ell_cat(x)   =the one-step infimal diagonal;
ell_hash(x)  =the asymptotic infimal-diagonal slope
              =the slope of every fixed diagonal.  (23)
```

Moving the minimizing context with `n` can change finite values, but it
cannot change the slope. The full kernel retains the finite interaction
that all three scalar traces forget.

## 6. Knot and tournament interpretation

For oriented knots, take `M` under connected sum and `d=d_G`. Then

```text
ell(K)=u(K),
ell_cat(K)=u_cat(K),
ell_hash(K)=u_hash(K).                              (24)
```

Schubert prime decomposition gives a canonical multiset of prime summands up
to permutation. After choosing labels for its occurrences, (5) defines an
**unknotting interaction complex** canonical up to relabelling and
simplicial isomorphism. A
Brittenham--Hermiller symbiont pair is a two-vertex minimal nonface.
Equations (6)--(14) show what the binary symbiont relation omits: a
higher-order minimal nonface could have every pair additive and become
nonadditive only after a composite context is formed.

No such higher-order knot example is claimed here. The theorem says that
pairwise data cannot exclude it for formal metric-monoid reasons. The
correct search record for a finite knot dictionary is therefore

```text
labelled prime packet;
weighted defect Delta(S) on tested subsets;
minimal-nonface antichain;
continuation-kernel witnesses for every positive defect;
pinned, catalytic, and asymptotic traces.            (25)
```

A tournament records an orientation of every pair. It can encode neither
an unoriented minimal nonface of arity greater than two nor the composite
vertex `x_S` appearing in (7) unless those data are added as sidecars.
The natural finite shadow is therefore a weighted simplicial complex or
hypergraph. A tournament becomes legitimate only after an intrinsic
asymmetric comparator on its faces has been specified.

## 7. Scope and failure boundary

1. The interaction complex depends on the labelled decomposition packet;
   it is not a complete invariant of the total knot.
2. The full continuation kernel remains the faithful operation carrier.
   Equations (16) and (23) are compressions of it.
3. Arbitrary-arity minimal nonfaces are proved for abstract weighted word
   metrics, not for knots.
4. Neither the complex nor the trace hierarchy computes an unknown
   unknotting number or produces a positive Gordian catalyst.
5. The theorem gives no reason to orient symmetric defect data as a
   tournament.

The exact reframe is:

```text
scalar value        -> pinned kernel entry;
binary symbiosis    -> missing edges of a zero-defect complex;
multiway symbiosis  -> minimal nonfaces;
catalysis           -> one-step infimal diagonal;
stable behavior     -> asymptotic min-plus trace slope;
exact composition   -> the full continuation kernel.             (26)
```

QED.
