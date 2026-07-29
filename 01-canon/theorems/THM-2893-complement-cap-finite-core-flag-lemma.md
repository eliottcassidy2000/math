---
id: THM-2893
title: Complement-cap finite-core flag lemma
status: PROVED
source: root-2026-07-29
depends_on: []
related:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-753-safe-peel-reduction-to-irreducible-cores
  - THM-2888-eight-body-first-apex-global-pair-cap-atlas
  - THM-2894-unmarked-residual-semilattice-order-and-group-clutch-no-go
---

# THM-2893 -- complement-cap finite-core flag lemma

## 1. Abstract carrier and complement cap

Let `(Omega,mu)` be a finite measure space, let `C` be a measurable carrier
with

```text
h=mu(C)>0,
```

and let `(D_v)_(v in V)` be a family of measurable danger sets.  For every
finite `S subset V`, put

```text
U_C(S)=mu(C intersect union_(v in S) D_v),     c_C(v)=U_C({v}).       (1)
```

Fix integers `1<=s<k`.  Suppose that a number `B` is a uniform
`(k-s)`-set cap:

```text
U_C(Q)<=B       for every Q subset V with |Q|=k-s.                    (2)
```

Put

```text
theta=h-B,
E_s={S subset V: |S|=s and U_C(S)>=theta},
H_s={v in V: c_C(v)>=theta/s}.                                      (3)
```

Call `E_s` the heavy `s`-graph and `H_s` its high-coverage core.

## 2. Cover-to-flag theorem

If `K subset V`, `|K|=k`, covers `C`, then:

1. every `s`-subset of `K` lies in `E_s`;
2. `H_s` meets every edge of `E_s`; and consequently
3. `|K intersect H_s|>=k-s+1`.

Equivalently, every hypothetical `k`-cover induces the complete
`s`-uniform hypergraph on its vertices, inside the heavy hypergraph, and all
but at most `s-1` of its vertices lie in the high-coverage core.

### Proof

Take `S subset K` with `|S|=s` and put `Q=K minus S`.  Subadditivity and
`(2)` give

```text
h=U_C(K)<=U_C(S)+U_C(Q)<=U_C(S)+B.
```

Thus `U_C(S)>=h-B=theta`, proving the first assertion.

For every `s`-set `S`,

```text
U_C(S)<=sum_(v in S)c_C(v).                                         (4)
```

If `S` missed `H_s`, every summand in `(4)` would be strictly below
`theta/s`, so `U_C(S)<theta`; hence no heavy edge misses `H_s`.  Finally,
if `K minus H_s` had at least `s` vertices, those vertices would contain an
`s`-edge of the complete heavy hypergraph on `K`, contradicting that `H_s`
meets every heavy edge.  Therefore `|K minus H_s|<=s-1`, which is the third
assertion.  This also audits the equality convention: heavy and high both
use `>=`, while the contradiction uses the strict complement `<`.  QED.

## 3. The strict discrepancy cutoff

Now let the vertices be positive integers, up to any fixed finite set of
excluded labels, and suppose that for all `w>=w_0`

```text
c_C(w)<h/7+gamma/w,       gamma>=0.                                 (5)
```

If

```text
B<(7-s)h/7,                                                     (6)
```

then `H_s` is finite.  More exactly, with

```text
beta=theta/s-h/7>0,       N=ceil(gamma/beta)-1,                    (7)
```

every `w in H_s` with `w>=w_0` satisfies `w<=N`.

Indeed, membership in `H_s` and `(5)` imply

```text
theta/s<=c_C(w)<h/7+gamma/w,
```

so `beta<gamma/w` and therefore `w<gamma/beta`.  For integral `w`, this is
exactly `w<=ceil(gamma/beta)-1`.  Condition `(6)` is equivalent to
`beta>0`.

Since a genuine cap has `B>=0`, `(6)` can hold only for `s<=6`.  Its
strictness is structural, not cosmetic.  At the feasible equality boundary
the core threshold is only `h/7`; estimate `(5)` permits infinitely many
labels with coverage exactly `h/7`, so it supplies no finite core.  For
example, repeat one set of measure `h/7` at infinitely many labels and
retain the non-sharp cap `B=(7-s)h/7`.

### Finite first-apex gate before any complement cap

Estimate `(5)` also gives a finite hitting statement without `(2)`.  For
any number of remaining labels `1<=p<=6`, put

```text
W_p=max(w_0, floor(7p gamma/((7-p)h))+1).                     (7a)
```

If all `p` labels satisfy `w>=W_p`, then

```text
U_C(K)<=sum_(w in K)c_C(w)
       <p(h/7+gamma/W_p)<h.                                  (7b)
```

Thus every `p`-cover meets the finite set of allowed labels below `W_p`.
If every nonempty literal residual again obeys an estimate of the form `(5)`—as
the LRC interval carriers do—then the same argument applies after each
possible first apex.  Since the depth is at most six, this yields a finite
rooted decision tree for every fixed initial carrier.  It is only a
finiteness theorem: successive cutoffs may be enormous.  Complement caps
and heavy flags compress that tree into the much smaller objects used below.

### Ordered first-apex suffix refinement

There is a lossless refinement whenever a finite hitting set is retained as
an ordered sidecar.  Choose any ordering

```text
A=(a_1,...,a_K)
```

of a set that meets every `p`-cover, and place the labels outside `A` after
this ordered prefix.  A putative cover has a unique least-indexed member
`a_r` of `A`.  After subtracting `D_(a_r)`, its other `p-1` labels all lie
in the strict suffix

```text
V minus {a_1,...,a_r}.                                      (7c)
```

Consequently it is enough, for each `r<=K`, to exclude `(p-1)`-covers of
the literal residual by that suffix; allowing the earlier prefix labels
again is unnecessary.  The order need not be decreasing coverage.  That is
one canonical choice in the LRC application, but an elimination order may
instead be chosen to reduce the later residual workload.

The index `r` and suffix in `(7c)` are retained branch data.  They are not
reconstructed from the unmarked residual, which by THM-2894 has forgotten
the order of deletions.  Thus this refinement uses precisely the marked
sidecar that the residual semilattice itself lacks.

### Monotone scalar bootstrap on the apex gate

The order search has a finite monotone formulation.  For `a in A`, put

```text
C_a=C minus D_a,       h_a=mu(C_a),
c_a(v)=mu(C_a intersect D_v).
```

For `P subset A minus {a}`, define

```text
S_a(P)=sup {
  sum_(v in Q)c_a(v) :
  Q subset V minus (P union {a}), |Q|=p-1
}.                                                          (7d)
```

Take the supremum to be zero if there are fewer than `p-1` allowed labels.
If

```text
S_a(P)<h_a,                                                 (7e)
```

then subadditivity excludes every `(p-1)`-cover of `C_a` from that suffix.
Moreover `S_a(P)` is nonincreasing as `P` grows, so this scalar certificate
is monotone under further prefix deletion.

Now reserve an ordered seed `H=(b_1,...,b_m)` whose branches are
independently certified on their actual suffixes
`V minus {b_1,...,b_i}`.  Starting with `P={b_1,...,b_m}`, repeatedly append
any `a in A minus P` satisfying `(7e)` and add it to `P`.  If this process
reaches `P=A`, the seed order followed by any activation order satisfies
the ordered suffix refinement and excludes every `p`-cover.

Equivalently, the inflationary map

```text
Phi(P)=P union {a in A minus P : S_a(P)<h_a}
```

has a least fixed point above the seed.  Finiteness and monotonicity show
that fair sequential activation and simultaneous rounds reach the same
fixed point.  Finding a smallest seed is therefore a target-set-selection
problem on a directed threshold hypergraph.  A seed is only a list of
branches still needing nonscalar certificates; its membership alone proves
nothing.  This is the proof-obligation analogue of THM-753's safe-peel core,
not a tournament obtained by forcing pairwise orientations.

### Longest-component singleton seal

For the LRC danger comb

```text
D_w={t in R/Z : ||wt||<=1/14},
```

each tooth is a closed interval of length `1/(7w)`.  Let `R` be a nonempty
finite union of intervals and let `lambda` be the length of its longest
connected component.  If one comb covers `R`, then it covers that entire
component.  A connected interval contained in the disjoint union `D_w`
lies in one tooth, so necessarily

```text
lambda<=1/(7w),       w<=floor(1/(7 lambda)).             (7f)
```

Thus the one-label leaves of the residual ladder have an exact geometric
cutoff independent of discrepancy.  It is the literal toothpick scale:
a tooth cannot cover a longer residual needle.

For an interval component `(a,b)` in a fixed lift of the circle, containment
in `D_w` is equivalent to the existence of an integer `q` with

```text
wb-1/14<=q<=wa+1/14.
```

Consequently the exact noncontainment test is

```text
ceil(wb-1/14)>floor(wa+1/14).                            (7g)
```

Equations `(7f)` and `(7g)` turn a singleton residual obligation into a
finite integer ledger.  They supplement rather than replace the
discrepancy cutoff: discrepancy controls high cores and multi-label tails,
while the longest-component law is sharper at the final one-label leaf.

## 4. Literal residual recursion

Assume `(6)`, and choose an integer

```text
1<=ell<=k-s+1.                                                     (8)
```

For `L subset H_s`, define the literal residual

```text
C_L=C minus union_(v in L)D_v.                                    (9)
```

Suppose every `ell`-set `L subset H_s` whose every `s`-subset is heavy
(with this flag condition vacuous when `ell<s`) has the following property:

```text
C_L is not covered by any k-ell allowed labels outside L.          (10)
```

Then `C` has no `k`-cover.

For if `K` were a `k`-cover, Section 2 would give at least `k-s+1`
vertices of `K` in `H_s`.  Choose `ell` of them as `L`.  Every
`s`-subset of `L` is heavy, because every `s`-subset of `K` is heavy, and
the remaining `k-ell` vertices of `K` cover `C_L`, contradicting `(10)`.
Thus a global covering problem reduces to finitely many literal residual
problems indexed by heavy flags in `H_s`.

The word *literal* is essential: the cap and the flag remember only union
sizes.  They do not preserve the pointwise overlap pattern needed to decide
whether the remaining labels cover `C_L`.

## 5. The two LRC rungs

For the lonely-runner danger combs, `(5)` is the standard strict
`h/7+gamma/w` carrier discrepancy form.

- **Four slots after one apex (`j=5`).**  Take `(k,s)=(4,2)`.  A global
  pair cap `B_2<5h/7` forces a heavy `K_4`, at least three of whose vertices
  lie in the finite core.  Taking `ell=3` reduces the proof to literal
  residuals behind finite heavy triangles, with one label left.  This is
  the mechanism isolated in THM-2888.
- **Five slots after one apex (`j=6`).**  Take `(k,s)=(5,2)`.  A global
  triple cap `B_3<5h/7` forces a heavy `K_5`, at least four of whose vertices
  lie in the finite core.  One may take `ell=3` and close residual pairs, or
  `ell=4` and close residual singletons.  Alternatively, `(k,s)=(5,3)`
  shows that a stronger pair cap `B_2<4h/7` forces a finite heavy
  three-graph and at least three core vertices.

This is the precise self-similarity of the covering ladder: cap a
complementary block, enumerate a finite heavy flag, subtract it literally,
and apply the same operation to the smaller residual.  The binary `s=2`
case is an intrinsic undirected relation—vertices are allowed speeds and an
edge records the symmetric observable `U_C({x,y})>=theta`.  No orientation
is present, so replacing this graph by a tournament would discard ties and
add non-mathematical data.

## 6. Scope

The cover-to-flag theorem is all-`k` under its stated cap; the cap-free
first-apex corollary uses the LRC denominator seven and `p<=6`.  The
theorem proves neither the required complement caps nor any residual inequality.
THM-2888 supplies the pair cap on its stated eight-body carriers; a
separate exact census is required to close their heavy-triangle residuals.
For the seven-body `j=6` rung, the global triple cap and subsequent residual
caps remain separate obligations.  In particular, this theorem alone does
not prove either rung and does not prove LRC(14).
