---
id: THM-4463
title: "Homogeneous tournament contacts and the exact pinned scale floor"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED.
  Zero contacts are precisely matching swaps in the path graph of homogeneous
  pairs. The bounded quotient-forced envelope branch has exact floor one;
  the pair-free branch has floor at least N. A fixed integral pinned context
  has min(N,c)<=phi_D(N)<=c, hence exact eventual floor c, where
  c=min_homogeneous_pair(1+2Duv). Strong-component and rooted pair-expansion
  descriptions are retained. This does not determine the full internal
  reversal kernel, unpinned exchanges, actual Hamilton response maximizers,
  or H>=disc.
source: tournament-continuation-20260908
depends_on:
  - THM-2249-directed-triangle-forced-quotient-frustration
  - THM-2256-automorphism-contact-dichotomy-for-quotient-frustration
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-4162-rooted-pair-mixed-two-ear-tensor-and-enumeration-free-johnson-cosets
  - THM-4163-order-eleven-homogeneous-pair-johnson-centrality
  - THM-4460-adverse-budget-continuation-state-and-shear-composition-collapse
script: 04-computation/tournament_continuation_metric_20260908.py
output: 05-knowledge/results/tournament_continuation_metric_20260908.out
script_sha256: c48db128df96ba6390ba382d734a18e6e9ca3c3e693020ab425d55d811b44e7b
output_sha256: 20a0b3e17dea47a292e74c610122e352b6b371f065609ad2c3e80799ac20a7fe
hash_basis: raw LF bytes
audit: >
  Independent root and two peer reviews accept the external-row identity,
  endpoint metric repair, homogeneous-pair path classification, matching
  contact count, strong-component criterion, sharp automorphism floor, and
  pinned integral-pseudometric invoice. One peer independently ran -O.
  All labelled tournaments of orders two through five give 124468 permutation
  probes; the complete order-three pinned bank checks 329308 nonfree Hall
  multisets. Normal and optimized output bytes agree with every explicit
  exception check active. No Lean or literature-priority claim.
---

# THM-4463 -- homogeneous-pair contacts determine tournament transport scale

**Status: PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED;
promoted on 2026-09-08 after root and two peer reviews.** This concerns the
quotient-forced reversal envelope. The full internal response, unpinned
exchange, and the live strong-tournament `H>=disc` problem remain OPEN.

## Inheritance and the decisive reframing

The closest proved mechanism is
[THM-2256, automorphism-contact dichotomy](../../01-canon/theorems/THM-2256-automorphism-contact-dichotomy-for-quotient-frustration.md).
It classifies the forced quotient envelope as bounded or linear according
to a finite but permutation-sized zero-contact test. Its construction comes
from [THM-2249, forced quotient frustration](../../01-canon/theorems/THM-2249-directed-triangle-forced-quotient-frustration.md).
The least-used sidecar is the actual external vertex comparison in each
ordered interaction term. Looking at that coordinate reveals an elementary
graph on quotient vertices and fully resolves the bounded branch.

The Anchor is this transport-scale classification inside actual tournament
substitutions. The Niche is the strong-component decomposition of its
homogeneous-pair obstruction. The Wildcard asks whether a pin that separates
all such pairs creates a linear tax. A singleton pin fails that cheap hostile:
it pays a constant invoice, however large the block size becomes.

The canonical controls are
[THM-2183, order-join metric product](../../01-canon/theorems/THM-2183-order-join-is-an-exact-tournament-metric-product.md),
[THM-2195, substitution product boundary](../../01-canon/theorems/THM-2195-transitive-quotients-exactly-control-universal-substitution-products.md),
and [THM-2221, pinned cut-metric response](../../01-canon/theorems/THM-2221-tournament-context-cut-metric-and-pinned-transport-response.md).
Their distinctions remain load-bearing: fixed block orders, quotient
automorphisms, pins, and the internal kernel cannot be discarded.
MISTAKE-249's one-sided versus two-sided continuation correction is respected:
permutation multiplication here is noncommutative, and zero contacts are not
called a congruence. The related
[THM-4460, adverse-budget continuation state](../../01-canon/theorems/THM-4460-adverse-budget-continuation-state-and-shear-composition-collapse.md)
supplies the warning that present zero cost need not be a removable state.
Here the exact source of failure is pair-dependent endpoint deletion.

| Live concept | New comparison |
|---|---|
| Quotient interaction | Zero contact reduces to matchings of homogeneous pairs |
| Strong components | Internal homogeneous pairs plus adjacent singleton components |
| Full Hamming metric | Restoring the omitted endpoints adds one per moved vertex |
| Pinned cut metric | Distinguishing pairs and paying a scale-sized bill are separate |
| Full reversal kernel | Quotient envelope is only a literal subset of its reversed pairs |

The source object is the Hall-layer transport of a labelled quotient
tournament `R`. The target is its graph of homogeneous vertex pairs; the map
expands each two-way layer interaction into external row disagreements. It
preserves all zero contacts, their diagonal reversal costs, and the bounded/
linear forced-envelope distinction. It does not preserve arbitrary positive
interaction values after taking only the zero graph, nor internal block
orientations. The full interaction matrix and the core kernel `G(X)` are the
respective restoration sidecars.

The intrinsic directed arcs remain the arcs of `R`. The homogeneous-pair
graph, cut metric, and layer interaction are symmetric objects; no tournament
orientation or artificial tie-breaking is imposed on them.

## 1. Contact is external row separation

Let `R` be a tournament on `C`, `|C|=q>=2`, with adjacency `A` and diagonal
zero. Every transport matrix `X` below is a **nonnegative integer** `q x q`
matrix with every row and column sum equal to the same integer `N>=1`.
Use THM-2249/2256's definitions

```
f(sigma,tau) = #{(i,k): i->k and tau(k)->sigma(i)},
w(sigma,tau)=f(sigma,tau)+f(tau,sigma),
d(tau)=f(tau,tau).
```

For distinct vertices define

```
h(u,v)=#{x outside {u,v}: A(u,x)!=A(v,x)},  h(u,u)=0.
```

Thus `h(u,v)=0` means `{u,v}` is a homogeneous pair: every exterior vertex
sees its two members identically. This definition does not assert that their
full adjacency rows are equal; their mutual arc prevents that.

**PROVED.**

```
w(id,tau) = sum_v h(v,tau(v)).                             (1)
```

Indeed, for each `v`, `f(id,tau)` counts paths `tau(v)->x->v`, and
`f(tau,id)` counts paths `v->x->tau(v)`. These are exactly the two ways
an exterior vertex distinguishes the pair. If `sigma` is an automorphism,
relabel the target by `sigma^{-1}` to obtain

```
w(sigma,tau)=sum_v h(v,(sigma^{-1} tau)(v)),
d(tau)=d(sigma^{-1} tau).                                  (2)
```

There is a useful endpoint repair:

```
Hamming(A(u,-),A(v,-)) = h(u,v)+1_[u!=v].                  (3)
```

Exactly one of the two endpoint coordinates differs. The full-row Hamming
distance is a metric. The endpoint-deleted `h` need not be a pseudometric:
in the ordered transitive triple, `h(0,1)=h(1,2)=0` but `h(0,2)=1`.
Consequently its zero relation cannot be collapsed as though it were an
equivalence. This is a precise lost-coordinate mechanism, not shared syntax
with an unrelated metric problem.

## 2. Complete zero-contact classification

Let `J(R)` be the undirected graph whose edges are homogeneous pairs, oriented
temporarily by their existing tournament arcs.

**PROVED.** Its connected components are consistently oriented paths. Their
vertex sets are exactly the maximal transitive modules of `R`.

For the first assertion, a vertex has at most one incoming and at most one
outgoing homogeneous-pair edge. If `v` beat two such neighbors `a,b`, say
`a->b`, then `a` would distinguish `v,b`; the incoming case is dual.
Thus each nontrivial component is a directed path or directed cycle.
A cycle `v0->v1->...->v(k-1)->v0` is impossible: successive homogeneity of
`{v1,v2},...,{v(k-2),v(k-1)}` transports `v0->v1` to `v0->v(k-1)`, a
contradiction. The same argument along a path proves every earlier vertex
beats every later one. An outside vertex has constant incidence along all
edges of the path, so the component is a module. Conversely, the consecutive
vertices of any transitive module form homogeneous pairs in the whole
tournament, so that module lies in one component. This proves maximality.

**PROVED (all zero contacts).** For any automorphism `sigma`,

```
w(sigma,tau)=0
 iff sigma^{-1} tau is a product of disjoint transpositions
     along edges of J(R).                                (4)
```

By (2), every nonfixed point must be sent to a graph neighbor. A permutation
cycle of length at least three would give a cycle in the path forest, which
is impossible. Its nontrivial cycles are therefore disjoint edge swaps.
The converse follows directly from (1).

If that matching has `k` edges, then

```
d(tau)=k.                                                (5)
```

Each swapped homogeneous pair reverses its own arc and preserves every
other incidence. Disjoint swaps therefore reverse exactly their `k` arcs.
In particular, a nonempty matching is never an automorphism.

If the component orders of `J(R)` are `l_1,...,l_b`, the exact number of
nontrivial zero-contact partners of each automorphism is

```
product_j Fibonacci(l_j+1) - 1.                          (6)
```

This is the elementary matching recurrence for a path; isolated vertices
contribute `Fibonacci(2)=1`. The total number of ordered
automorphism/nonautomorphism contacts multiplies (6) by `|Aut(R)|`.

Zero-contact permutations still do not form a monoid. In the transitive
triple, swaps `(01)` and `(12)` each have zero contact with the identity,
but their composition has contact one. The needed continuation sidecar is
the path position and matching, not only membership in one transitive module.

## 3. Strong components and the exact unpinned envelope dichotomy

The homogeneous pairs of `R` are precisely:

1. homogeneous pairs inside a single strong component; and
2. pairs of adjacent singleton components in the linear strong-component order.

Outside a strong component, incidence is uniform, so its internal pair test
agrees with the whole tournament. If a homogeneous pair belongs to different
components, an intervening component would distinguish it. Each of its two
components must be a singleton: otherwise homogeneity would make the chosen
vertex an internal source or sink, impossible in a nontrivial strong component.
The converse for adjacent singleton components is immediate.

Strong connectivity alone therefore does not decide this obstruction. The
strong four-vertex substitution `C3[T2,1,1]` contains a homogeneous pair,
whereas the nonstrong order-join `C3 join 1` contains none.

There is a direct connection to the existing Hamilton-response mechanism:
by [THM-4145, rooted homogeneous-pair expansion](../../01-canon/theorems/THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula.md),
contracting a homogeneous pair of a strong tournament gives a rooted strong
quotient `(Q,r)`, and re-expansion `P_r(Q)` recovers the tournament. Thus the
strong bounded-contact class is exactly the class of such rooted pair
expansions. The same map satisfies `H(P_r(Q))=F_Q(N_Q^+(r))` in that
theorem's ear convention, while its child continuation requires the mixed
two-ear field. The pair is cheap for the quotient-forced reversal observer,
but is not removable from the Hamilton response.
[THM-4162](../../01-canon/theorems/THM-4162-rooted-pair-mixed-two-ear-tensor-and-enumeration-free-johnson-cosets.md)
supplies the exact child tensor, and
[THM-4163](../../01-canon/theorems/THM-4163-order-eleven-homogeneous-pair-johnson-centrality.md)
already closes order-eleven Johnson support-floor centrality on this pair
branch. The older open wording in THM-4145 is not inherited here. None of
these connections extends that finite result to actual response maximizers
or proves `H>=disc`.

Let `phi_R(N)` be THM-2256's minimum forced energy outside the scaled
automorphism axes. **PROVED, all `N>=1`:**

```
J(R) has an edge  => phi_R(N)=1;
J(R) has no edge => phi_R(N)>=N and phi_R(N)=Theta_R(N).    (7)
```

The first claim uses `X=(N-1)I+P_(uv)` for any homogeneous pair: its cross
interaction is zero and diagonal cost is one. Positivity gives equality.
For the second, (4) rules out all zero contacts, and THM-2256 already supplies
the upper linear bound. Its lower bound strengthens from `N-1` to `N` using
the following sharp refinement of its automorphism interaction lemma:

```
sigma,tau distinct automorphisms
 => w(sigma,tau)>=|supp(sigma^{-1}tau)|>=3.                (8)
```

To prove (8), set `alpha=sigma^{-1}tau`. The number `h(v,alpha(v))` is
constant along each `alpha`-orbit. On a nontrivial orbit it cannot be zero:
otherwise that orbit would trace a cycle in `J(R)`. Automorphisms cannot
have two-cycles because they would reverse the mutual arc. Thus every moved
vertex contributes at least one in (1), and a nonempty support has size at
least three. Rotations of `C3` attain equality three.

For a Hall decomposition with `k` nonautomorphism layers, their diagonal
energy is at least `k`. If both types of layers occur and no zero contact
exists, the cross contribution is at least `k(N-k)`, giving at least `N`.
If all layers are nonautomorphisms the diagonal alone gives `N`. If all
layers are automorphisms but the transport is nonfree, (8) gives at least
`3(N-1)>=N` for `N>=2`; this case is absent at `N=1`.

This classifies the formerly unspecified bounded branch by a vertex-level
test, and gives its exact value. It does not say the full kernel `G_R(X)`
has a bounded floor: `F_R` counts only pairs with distinct source and target
quotient labels. THM-2183's strong order-join metric theorem is compatible
with `phi_(Tq)=1` because it retains the omitted internal-pair response.

## 4. The exact eventual floor with a fixed pinned context

Let `D` be any nonnegative integral cut semimetric supplied by THM-2221.
More generally the proof only needs an integral pseudometric. Define

```
H_D(X)=F_R(X)+<X,D>,
A_D={sigma in Aut(R): <P_sigma,D>=0},
phi_(R,D)(N)=min {H_D(X): X has row/column sums N,
                          X notin {N P_sigma:sigma in A_D}}.
```

For each fixed transport `X`, `H_D(X)` is a lower bound for its actual pinned
cost `G_R(X)+<X,D>`. The minimum `phi_(R,D)` constrains only the nonfree
transport branch; free-axis internal costs require a separate check. The
free axes are exactly the zero set of `H_D`, and the identity is always in
`A_D`.

If `J(R)` has an edge, put

```
c_D=min_({u,v} in E(J(R))) [1+2D(u,v)].
```

**PROVED.** For every integer `N>=1`,

```
min(N,c_D) <= phi_(R,D)(N) <= c_D.                        (9)
```

Consequently `phi_(R,D)(N)=c_D` for every `N>=c_D`. If `J(R)` is edgeless,
then `phi_(R,D)(N)>=N` and, for a fixed `D`, it remains `Theta_(R,D)(N)`.

The upper bound in (9) is the same one-swap transport as before, whose
exterior cost is exactly `2D(u,v)`. For the lower bound, take any Hall
decomposition and let `k` be the number of nonautomorphism layers and
`a=N-k` the number of automorphism layers.

- If `a=0`, diagonal energy is at least `N`.
- If at least two automorphism types occur, their mutual energy is at
  least `3(a-1)` by (8), so `H_D>=k+3(a-1)>=N`.
- Suppose exactly one automorphism type `sigma` occurs. If some bad layer
  has positive contact with it, then `F_R>=k+a=N`. If
  `<P_sigma,D>>0`, integrality gives `H_D>=k+a=N` instead.
- In the remaining case `sigma in A_D` and every bad layer has zero contact
  with it. By (4), each relative permutation is a nonempty matching.

In this last case `D(v,sigma(v))=0` for every `v`. The triangle inequality
therefore gives `D(v,sigma(x))=D(v,x)`. For a relative matching `rho`,

```
d(sigma rho)+<P_(sigma rho),D>
  = sum_({u,v} in rho) [1+2D(u,v)] >= c_D.                (10)
```

The Hall formula has nonnegative cross coefficients and a diagonal
`d(tau)n_tau^2>=d(tau)n_tau`. Hence its cost is at least the sum of the
individual matching invoices (10), and at least `c_D` if a bad layer occurs.
If none occurs the transport is an excluded free axis. This proves (9).
When `J(R)` is edgeless the last nonfree case cannot occur, giving the
stronger lower bound `N`. A fixed nonautomorphism layer amid identity layers
gives an upper bound linear in `N`, completing the claimed fixed-context scale.

Integrality is used to replace a positive contact/context invoice by at
least one. Arbitrary positive real pin weights require their minimum positive
invoice in the lower bound; they are not silently included in (9).

## 5. Why pair-separating pins need a scale invoice

A single pinned exterior vertex can distinguish every homogeneous pair:
choose its incidence word alternately along each component path of `J(R)`.
This is tournament-realizable because exterior incidences may be prescribed
freely. Then `D(u,v)=1` for every homogeneous pair, so `c_D=3` and

```
phi_(R,D)(N)=3 for N>=3                                (11)
```

whenever a homogeneous pair exists. Thus removing all *zero* exterior
separations still leaves a constant-size one-vertex exchange at arbitrarily
large block order. Treating the pin as a new quotient block of order `N`
would change its mass and hence the actual transport problem.

The sharp three-vertex control is the transitive triple and word `(0,1,0)`.
Its floor is three also for `N=1,2`, checked by complete Hall enumeration;
together with (11), this gives the all-size constant law in that control.
With the same pin repeated `N` times, the exact finite minima for `N=1..7`
are instead `3,5,7,9,11,13,15`. No all-size equality is inferred from these
seven values.

There is an all-size scale guarantee: apply (9) pointwise to a varying
integral context `D_N`. If every homogeneous pair satisfies
`D_N(u,v)>=kappa N` for a fixed `kappa>0`, then

```
phi_(R,D_N)(N)>=min(N,1+2kappa N).
```

Thus a linear amount of distinguishing context restores a linear tax in the
bounded branch. The sidecar is its actual exterior mass, not merely whether
its word separates labels. Pins remain fixed by admissible bijections;
unpinned exterior/core exchange is outside every claim here.

## Exact audit and remaining work

Run:

```
python 04-computation/tournament_continuation_metric_20260908.py
python -O 04-computation/tournament_continuation_metric_20260908.py
```

[Source](../../04-computation/tournament_continuation_metric_20260908.py) and
[frozen output](../../05-knowledge/results/tournament_continuation_metric_20260908.out).

Normal and optimized output bytes agree. Raw LF SHA-256 pins are:

```
script c48db128df96ba6390ba382d734a18e6e9ca3c3e693020ab425d55d811b44e7b
output 20a0b3e17dea47a292e74c610122e352b6b371f065609ad2c3e80799ac20a7fe
```

The complete labelled universe is every tournament of orders `2..5`, every
permutation, and every automorphism-relative contact pair. Direct directed-
path counting is compared with row separation, graph matchings, Fibonacci
counts, and strong-component decomposition. This gives 124468 permutation
probes and 2054 ordered nontrivial zero contacts. The recovered bounded/linear
census is `(2,0),(6,2),(48,16),(680,344)`, matching THM-2256's existing audit.
At order five, strong quotients split `(240,304)` and nonstrong quotients
split `(440,40)`; no inference to higher-order class counts is made.

The pinned bank exhausts every order-three quotient, all eight binary context
words, multiplicities `0,1,2`, and `N=1..7`. It enumerates 329308 nonfree Hall
multisets including the extra pin hostile controls. Every integer transport
has a Hall decomposition, so the minimum over all these multisets is the
minimum over all transports, even though different multisets may encode the
same matrix. Checks use explicit exceptions and remain active under `-O`.

Searches covered zero contact, homogeneous/twin pairs, path components,
bounded frustration, transposition kernels, and related correction entries.
The generic old Hall expansion and dichotomy are credited; the new scoped
content is the full matching classification, sharp automorphism floor,
strong-component criterion, and pinned-context invoice. No literature novelty
claim is attached to the elementary modular-decomposition fact.

The actual full reversal kernel `G_R(X)` can exceed `F_R(X)` and can depend
on internal orientations at precisely the scale omitted here. The next useful
question is a lower bound for that residual on a contact swap in a strong
quotient, with internal block markers retained. Neither the contact test nor
the pin invoice proves `H>=disc` or an LRC statement.
