---
id: THM-2221
title: "Tournament context cut metric and pinned transport response"
status: >
  PROVED. For equal-order tournament substitution blocks over one labelled
  quotient, with block-uniform pinned exterior context and unit Hamming
  reversal cost, the whole incidence-word histogram collapses exactly for
  every prescribed exterior transport-cost probe to an integer cut
  semimetric D_mu. For arbitrary split block transport X the exterior cost
  is <X,D_mu>, and the full pinned response is min_X[G(X)+<X,D_mu>].
  Whole-block transpositions recover every entry of D_mu. Bypass-
  suppression loci are additive upper ideals with finite antichain bases.
  Vertex marginals do not determine D_mu, while D_mu itself discards
  higher correlations and phase which can affect general whole-fiber complex
  amplitudes; no Krenn matching realization is asserted. The core kernel G
  and unpinned core/context exchange remain open.
source: klein-2026-07-24-tournament-context-cut-transport
depends_on:
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2216-residual-capacity-hinge-gram-law
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
external:
  - "Mario Krenn, Xuemei Gu, and Daniel Soltesz, Questions on the Structure of Perfect Matchings Inspired by Quantum Physics, arXiv:1902.06023."
---

# THM-2221 -- tournament context cut metric and pinned transport response

Assume throughout: one common labelled quotient tournament, equal positive
order `N` in every corresponding core block, unit arc-reversal/Hamming cost,
and block-uniform pinned exterior adjacency identical label by label in
source and target.  Admissible bijections may transport vertices arbitrarily
among the core blocks but may not exchange them with pinned exterior
vertices.  This is not an unrestricted formula for `d_iso`; it supplies the
exact exterior term of the transport kernel, while the core-pair kernel
remains a necessary sidecar.

## 1. Inheritance and the missing coordinate

The closest proved result is
`THM-2195-transitive-quotients-exactly-control-universal-substitution-products`.
For a prescribed whole-block permutation it writes the exterior cost as a
weighted Hamming derivative of the exterior incidence words.  Its sharp
hostile example is the cyclic transport around an induced directed triangle.
The theorem explicitly leaves an unmarked block-transport matrix as the next
object.

The whole-fiber warning comes from the Krenn--Gu--Soltesz inherited-color
model: an induced coloring must be aggregated over its entire perfect-matching
fiber before cancellation is assessed.  A pairwise orientation is not a
faithful replacement.  Here the coefficient algebra is more special:
nonnegative reversal cost with a Hamming kernel.  That special form permits an
exact, rigorously delimited collapse of the whole incidence-word fiber.

The repaired object is not a tournament.  It is an integer cut semimetric on
the core blocks, paired with a min-plus core transport kernel.

## 2. Pinned contexts and their color fibers

Let `C` be a finite set of core block labels, and let `R` be a tournament on
`C`.  As assumed above, every core block has the same order `N>=1`.  Let

```text
P = R[T_c : c in C],
Q = R[S_c : c in C],                 |T_c|=|S_c|=N.       (1)
```

A pinned exterior context consists of labelled vertices `t` which occur
identically in the source and target and which every admissible bijection
fixes.  Each `t` has uniform adjacency to all vertices inside a displayed
core block, with that block-label incidence identical in source and target.
Arcs among exterior vertices are also identical and hence cost zero.
The incidence color of `t` on the core is the binary word

```text
w_t in {0,1}^C,       w_t(c)=1_[c -> t].                    (2)
```

More generally an exterior block of order `m_t`, fixed as a whole and matched
internally at zero cost, contributes `m_t` copies of this word.  Thus the
whole-fiber context ledger is the histogram

```text
mu(w)=sum_(t:w_t=w) m_t,       mu in N^({0,1}^C).           (3)
```

Every such histogram is tournament-realizable: prescribe the core--context
arcs word by word and orient the still-free context--context pairs
arbitrarily.

Define

```text
D_mu(c,d)
  =sum_w mu(w)|w(c)-w(d)|
  =total exterior weight which separates c and d.          (4)
```

Each word `w` defines the elementary cut semimetric

```text
delta_w(c,d)=|w(c)-w(d)|.                                  (5)
```

Hence

```text
D_mu=sum_w mu(w) delta_w.                                  (6)
```

In particular `D_mu` is a symmetric integer-valued pseudometric: symmetry
and the zero diagonal are immediate, and the triangle inequality holds for
each `delta_w`, hence for their nonnegative sum.  Its zero relation has the
literal meaning

```text
D_mu(c,d)=0
 iff every positive-weight context word has w(c)=w(d).      (7)
```

Thus exact ties are observational equivalences.  There is no intrinsic
orientation to impose on them.

The image of `mu -> D_mu` is exactly the nonnegative integral cut semigroup on
`C`.  Constant words map to zero, and `w` and `1-w` generate the same cut.
Equivalently, if

```text
m(c)=sum_w mu(w)w(c),
m(c,d)=sum_w mu(w)w(c)w(d),
```

then

```text
D_mu(c,d)=m(c)+m(d)-2m(c,d).                         (7a)
```

Thus one-coordinate exterior scores lose precisely the pairwise
co-incidence needed even before any higher-order color question is asked.

### Signed Gram gauge and negative type

The incoming signed-coordinate viewpoint of THM-2195 identifies an
equivalent semidefinite chart.  Put

```text
M_mu=sum_w mu(w),
s_w(c)=2w(c)-1,

Gamma_mu=sum_w mu(w) s_w s_w^T.                    (7b)
```

Then `Gamma_mu` is positive semidefinite, has constant diagonal `M_mu`,
and

```text
D_mu(c,d)=[M_mu-Gamma_mu(c,d)]/2.                  (7c)
```

Indeed,

```text
|w(c)-w(d)|=[1-s_w(c)s_w(d)]/2
```

word by word.  Consequently `D_mu` is of conditional negative type:
whenever `sum_c a_c=0`,

```text
a^T D_mu a=-1/2 a^T Gamma_mu a<=0.                 (7d)
```

For every admissible transport matrix `X`, whose total mass is
`N|C|`,

```text
<X,D_mu>
 =[M_mu N|C|-<X,Gamma_mu>]/2.                      (7e)
```

Thus the cut-cost minimization is equivalently a signed-Gram trace reward.
Adding a constant context word changes `(M_mu,Gamma_mu)` by
`(1,mathbf1*mathbf1^T)` and leaves (7c)--(7e) invariant.  The cut metric is the
gauge-free observable; the raw Gram matrix is not.

## 3. Exact split-transport formula

Let `pi:V(P)->V(Q)` be an arbitrary bijection between the core vertices.  It
need not preserve the displayed blocks.  Its block transport matrix is

```text
X_pi(c,d)=#{x in T_c : pi(x) in S_d}.                       (8)
```

Every row and column sum of `X_pi` is `N`.  Conversely every nonnegative
integer matrix with these margins is realized by at least one core bijection.

Let `g(pi)` be the number of reversed pairs whose two vertices both lie in
the core, and define the exact core kernel

```text
G(X)=min_{pi:X_pi=X} g(pi).                                 (9)
```

> **Theorem 1 (cut-metric transport law).** After extending `pi` by the
> identity on the pinned context, its total reversal cost is
>
> ```text
> g(pi)+sum_(c,d) X_pi(c,d)D_mu(c,d).                       (10)
> ```
>
> Consequently the exact pinned-context response is
>
> ```text
> R_P,Q(mu)
>   =min_X [G(X)+<X,D_mu>],                                 (11)
> ```
>
> where the minimum is over the integer transportation matrices with all
> row and column sums `N`.

**Proof.**  Vertex pairs split into core--core, context--context, and
core--context pairs.  The first class costs `g(pi)`, and the second costs
zero.  Fix a source vertex `x in T_c` which maps into `S_d`, and an exterior
vertex of word `w`.  The source arc is encoded by `w(c)` and the target arc
by `w(d)`, so this pair reverses exactly when `w(c)!=w(d)`.  Summing over the
whole word fiber gives `D_mu(c,d)` for this moved core vertex.  Summing over
the `X_pi(c,d)` vertices moved from `c` to `d` proves (10).  Minimizing first
inside each transport fiber and then over `X` proves (11).  QED.

This is the precise first part of the block-transport object requested by
THM-2195.  The exterior term is linear in `X` and depends on the context
only through an `l1` cut metric.  The remaining core term `G(X)` is not
generally determined by `X` alone without performing the minimization in
(9): it retains internal tournament arcs and pairwise coupling of moved
vertices.

## 4. Recovery of the THM-2195 derivative

For a permutation `sigma` of `C`, let `P_sigma` denote its permutation
matrix with the convention that `c` maps to `sigma(c)`.  On the whole-block
face

```text
X=N P_sigma.                                                (12)
```

The core kernel is then

```text
G(NP_sigma)
 =sum_c d_iso(T_c,S_(sigma(c)))
  +N^2 #{ {c,d} :
       (c ->_R d) xor (sigma(c) ->_R sigma(d)) }.           (13)
```

Indeed, internal block minimizations are independent, while a quotient-arc
disagreement reverses all `N^2` cross pairs.  The exterior term becomes

```text
<NP_sigma,D_mu>
 =N sum_c D_mu(c,sigma(c))
 =N sum_w mu(w) delta_sigma(w),                             (14)
```

where

```text
delta_sigma(w)=sum_c |w(c)-w(sigma(c))|.                    (15)
```

Equations (13)--(15) recover the prescribed-block law of THM-2195.  If
`sigma` is an automorphism of `R[C]`, the quotient-arc term in (13) vanishes.
If `sigma` is a cycle, (15) is exactly the number of binary transitions
around that cycle.

## 5. Exact classification of prescribed exterior-cost probes

Write

```text
E_mu(sigma)=N sum_c D_mu(c,sigma(c)).                       (16)
```

This is the exterior reversal cost of the whole-block permutation `sigma`.

> **Theorem 2 (prescribed exterior-cost classification).** For two context
> histograms `mu,nu`, the following are equivalent:
>
> 1. `D_mu=D_nu`;
> 2. `E_mu(sigma)=E_nu(sigma)` for every permutation `sigma` of `C`;
> 3. `<X,D_mu>=<X,D_nu>` for every admissible split transport matrix `X`.

**Proof.**  The first statement immediately implies the other two.  Conversely,
take the transposition `tau=(c d)`.  Only `c,d` move, so

```text
E_mu(tau)=2N D_mu(c,d).                                    (17)
```

Equality of every permutation response therefore recovers every entry of
`D_mu`.  This proves `2=>1`; `3=>2` follows by restricting to the whole-block
matrices.  QED.

Thus `mu -> D_mu` is not merely a convenient bound.  It is the exact additive
quotient of the whole incidence-word ledger for all exterior Hamming
transport responses:

```text
D_(mu+nu)=D_mu+D_nu.                                       (18)
```

This is the operation-ready classification.  It should not be enlarged into a
tournament: its observable is symmetric, weighted, and legitimately tied.
It is not an injectivity statement for the single minimized scalar
`R_P,Q(mu)`: after pinning has forced identity transport, different larger
singleton pinning levels have different `D_mu` but the same saturated value
`R_P,Q=A_0`.

## 6. Min-plus context profile and finite pinning bases

Extend `mu` to nonnegative real weights.  Since the transportation set is
finite, (11) is a minimum of finitely many affine functions of `mu`.
Therefore:

```text
R_P,Q is piecewise-linear and concave;
mu<=nu coordinatewise implies R_P,Q(mu)<=R_P,Q(nu).         (19)
```

The monotonicity follows also directly from (10): adding pinned context adds
nonnegative tests to every transport.

Let

```text
A_0=G(NI)                                                  (20)
```

be the decomposition-respecting cost.  The identity transport has zero
exterior cost, hence `R_P,Q(mu)<=A_0`.  Define the block-bypass saving

```text
S_P,Q(mu)=A_0-R_P,Q(mu)
 =max_X [A_0-G(X)-<X,D_mu>].                               (21)
```

Then `S_P,Q` is nonnegative, piecewise-linear, convex, and coordinatewise
nonincreasing.  For every integer `s>=0`, the suppression locus

```text
J_s={mu in N^({0,1}^C): S_P,Q(mu)<=s}                      (22)
```

is an additive upper ideal:

```text
mu in J_s  =>  mu+nu in J_s for every nu.                  (23)
```

It is nonempty.  Let `mu_K` consist of `K` copies of every singleton cut
word `1_{c}`.  Their sum gives

```text
D(c,d)=2K for c!=d.                                        (24)
```

Put

```text
q(X)=sum_(c!=d)X(c,d).
```

Then the exact added exterior cost is `2Kq(X)`.  Every nonidentity integer
transport with equal row and column margins has `q(X)>=2`: its nonzero
off-diagonal flow has zero divergence and therefore contains a directed
cycle of length at least two.  Identity
transport is optimal exactly when

```text
2Kq(X)>=A_0-G(X)             for every X!=NI,        (24a)
```

and it is uniquely optimal if all these inequalities are strict.  Since
`G(X)>=0`, the crude but explicit sufficient condition

```text
4K>A_0                                                     (24b)
```

forces every `X!=NI` to cost more than identity.  Hence `R_P,Q=A_0` and
`S=0`.  No sharpness is claimed for the sufficient bound (24b).

Because the word alphabet is finite, Dickson's lemma applies to
`N^({0,1}^C)`: every upper ideal has finitely many minimal elements.
Consequently each `J_s` has a finite antichain basis of minimal pinning
contexts.

This is the operation-sensitive dual of the knot context response in
THM-2191:

```text
knot common context:       distance response decreases,
                           saving loci are upper ideals;

tournament pinned context: distance response increases,
                           bypass-suppression loci are upper ideals.          (25)
```

The reversal of monotonicity is structural.  A knot summand creates possible
common intermediates; a pinned tournament context creates additional
incidence tests.

## 7. Two hostile examples and the exact loss ledger

### 7.1 Vertex marginals and cosmetic orientations fail

Let `C={0,1,2}` and compare the two context ledgers

```text
mu_const = [000]+[111],
mu_split = [100]+[011].                                   (26)
```

Both have the same one-coordinate marginals:

```text
sum_w mu(w)w(c)=1            for c=0,1,2.                  (27)
```

Thus every ranking or binary orientation made only from exterior win counts
sees three exact ties in both cases.  Nevertheless,

```text
D_const=0,

D_split(0,1)=2,   D_split(1,2)=0,   D_split(2,0)=2.         (28)
```

For the cyclic permutation `rho=(0 1 2)`,

```text
E_const(rho)=0,
E_split(rho)=4N.                                          (29)
```

If the core is the directed triangle, `rho` is precisely the partial
automorphism used in THM-2195.  In the first context it extends across the
exterior; in the second it does not.  Two exterior vertices are the first
possible size for this marginal-collision mechanism: with one exterior
vertex its marginal vector is the word itself.

The first failed implication is therefore

```text
same exterior score/marginal data
   does NOT imply same block-transport response.            (30)
```

The repaired invariant is the cut metric `D`, equivalently the pairwise
separation counts, with zero distances retained as ties.

### 7.2 The cut metric is exact for Hamming transport but not for arbitrary
whole-fiber amplitudes

Now compare

```text
mu_even = [000]+[011]+[101]+[110],
mu_odd  = [111]+[100]+[010]+[001].                         (31)
```

Both have coordinate marginals two and

```text
D(c,d)=2             for every distinct c,d.                (32)
```

By Theorem 2, every exterior block-transport response is identical.  Yet the
three-way parity functional distinguishes the whole fibers:

```text
sum_w mu_even(w)(-1)^(w(0)+w(1)+w(2)) =  4,
sum_w mu_odd (w)(-1)^(w(0)+w(1)+w(2)) = -4.                (33)
```

Thus the map `mu -> D_mu` destroys higher color correlation.  That loss is
harmless for reversal cost because every kernel (15) is a sum of two-point
cut functions.  It can matter for a general sum-product or complex-amplitude
whole-fiber observable, which requires the full weighted ledger (including
phase), not merely `D`.  The parity witness is an information-theoretic
boundary; no realization as a Krenn perfect-matching instance is asserted.

## 8. Connection contract

```text
source:
  exterior incidence-word fibers mu on {0,1}^C;

target:
  the integer cut semimetric D_mu and the min-plus response
  min_X [G(X)+<X,D_mu>];

map:
  mu |-> sum_w mu(w) delta_w;

preserved:
  every pinned exterior reversal cost, for arbitrary split vertex transport,
  and hence the exact pinned context response once G is retained;

destroyed:
  constant context words, complement orientation, and higher-order
  color correlations/complex phase;

needed sidecars:
  G(X) for core-pair reversal cost;
  the full complex fiber ledger for general whole-fiber cancellation;
  a larger membership/transport kernel if exterior vertices themselves may
  exchange with core vertices;

cheapest decisive tests:
  a transposition recovers D(c,d);
  (26) refutes marginal/tournament compression;
  (31)--(33) exhibit the exact higher-order information discarded.           (34)
```

## 9. Honest frontier

This theorem supplies the exact exterior functional for the unmarked
block-transport program proposed by THM-2195, including transports which
split the core blocks.  It does not compute `G(X)`, prove that an unrestricted
optimal isomorphism preserves the core/context partition, or collapse
complex whole-fiber amplitudes.

The next high-leverage tournament problem is now sharply typed:

```text
classify or bound the core kernel G(X)
subject to the integer transport margins,
then minimize G(X)+<X,D>.                              (35)
```

For a transitive quotient, THM-2183's uncrossing says the unique relevant
transport is the decomposition-respecting one.  Directed triangles create
low-`G` cyclic transports.  The cut metric `D` measures exactly how much
pinned exterior context suppresses those transports.  This is an
assignment/transport problem on a symmetric semimetric, not a tournament on
contexts or on knots.
