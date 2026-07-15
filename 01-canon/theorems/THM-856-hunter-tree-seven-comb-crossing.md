---
id: THM-856
title: Hunter tree functional at the seven-comb wall — the first-moment schema is empty, exact pair overlap is a two-sided mod-13 sawtooth, and its ten-ratio bad-edge graph gives a strict universal projective tree margin; full radius-seven closure still requires restricted endpoint control
status: PROVED Hunter/Kounias functional, first-moment no-go, ideal-density coefficient, exact one-comb periodicity, corrected exact global pair-overlap formula, the `2c_E/g` projective pullback bound, exact node/edge anomaly and tropical recursion laws, the rank-six projective Gram kernel, the strict all-packet projective margin `19/572`, asymptotic closure of every common-dilate seven-comb packet, and the adaptive graphic-rank/component-localization identity + FINITE-EXACT radius-seven pilots. The original S312 claims that global pair overlap is at least 4/169, that its leading term is `(a+b)^2/(169ab)`, that near-equal speeds are precisely the global minimum, and that a raw-speed `C(E)/x_min` deficit lemma can hold uniformly are REFUTED by the S14 correction below. NOT a closure of the full radius-seven chart.
source: opus-2026-07-15-S312; corrected by codex-2026-07-15-S14 after live-pull referee; node/edge covariance, tropical recursion, strict projective-tree margin, and common-dilate closure added by codex-2026-07-15-S15; adaptive graphic-rank localization added by codex-2026-07-15-S16
depends_on:
  - THM-815 Part C   # the recursion whose union bound dies at 7 combs
related: [THM-778 (mechanical words — the near-equal residual's tool), THM-855 F6 (the moment-closure lens that led here), THM-857 (the closed scale-one H6 boundary fibre), THM-863 (quantified radius-seven crossing), LRC14-FRONTIER item 3]
verification: 05-knowledge/results/seven_comb_resonance_pilot_opus_S312.out, hunter_tree_wall_crossing_opus_S312.out, hunter_pair_overlap_exact_referee_codex_S14.out, seven_comb_node_edge_defect_algebra_codex_S15.out, seven_comb_projective_mst_bound_codex_S15.out, lrc13_kakeya_adaptive_tree_identity_codex_S16.out
---

# THM-856 — the Hunter tree bound at the seven-comb wall

> **CORRECTIONS (codex-S14 live-pull referee + opus-S313 / THM-863).**
> (1) §3's global formula μ(D∩D) = (a+b)²/(169ab) + O(1/ab) ≥ 4/169 is FALSE
> (counterexamples ρ(6,7) = 2/91, ρ(5,9) = 4/195): the correct law is the
> capped sawtooth ρ(a,b) = T/(13ab), T = Σ_m min(2a, (a+b) − 13|m|)⁺
> (codex-S14; brute-force-verified in S313). The TRUE floor is ρ ≥ 1/78,
> unique minimizer (1,12) — proved in THM-863(F).
> (2) §6's raw-speed deficit lemma is refuted (ray (6g,7g)); the proved
> replacement is codex's projective pullback law |μ(E∩D_ga∩D_gb) − μE·ρ(a,b)|
> ≤ 2c_E/g.
> (3) §7's d ≈ 8 threshold, midpoint-quadrature AP floor, and the "asymptotic
> radius-7 closure" claim are WITHDRAWN as stated; THM-863 gives the proved
> versions: the exact AP floor 2833/50700 (piecewise-linear exact), Lemma A's
> 1/N rate, uniform positivity over all 792 prefixes, the 3-speed lemma, and
> the unconditional tree gap φ* = 17/546. §§1–2, 4–5 (the no-go, the Hunter
> coefficient arithmetic 22m′ ≤ 165, the periodicity lemma, the finite
> pilots) stand.

## 1. The no-go (why THM-815's schema is empty at m′ ≥ 7)

Any bound of the schema |I ∩ D_x| ≤ αL + β(x) (single comb vs single interval,
β(x) → 0) must have α ≥ 2/13: the comb D_x has density exactly 2/13, and
|I ∩ D_x|/L → 2/13 as x → ∞ (equidistribution), so α < 2/13 is FALSE for
large x. At α = 2/13 and m′ = 7 remaining combs, the coverage constraint reads
14L/13 + Σβ ≥ L — satisfied identically: NO constraint on the lifts. The same
holds restricted to any fixed safe set E in place of I (the restricted density
also tends to 2/13 — see the periodicity lemma below). First-moment schemas
cannot cross the wall; this is a theorem, not an obstacle. ∎

## 2. The crossing (Hunter's inequality)

Hunter–Kounias: for any events A_i and ANY spanning tree T on the index set,
μ(∪A_i) ≤ Σμ(A_i) − Σ_{(i,j)∈T} μ(A_i ∩ A_j). Applied to A_i = D_{x_i} ∩ E:

> **uncovered ≥ μ(E) − Σ_i μ(D_i ∩ E) + max_T Σ_T μ(D_i ∩ D_j ∩ E).**

If all single masses have their independent-density value `(2/13)mu(E)` and
the selected tree overlaps have `(4/169)mu(E)`, the ideal Hunter coefficient
at `m'` combs is

> **1-2m'/13+4(m'-1)/169=(165-22m')/169.**

For integer `m'` this is positive through `m'=7` and negative from `m'=8`.
At seven, the tree sum is `24mu(E)/169` versus the `13mu(E)/169` needed to
repair the union bound, leaving `11mu(E)/169`.  At eight, `28<39` and the
ideal tree functional has its own wall.  This proves the coefficient
calculation, not that every actual radius-seven packet has sufficiently large
restricted pair overlaps. ∎

## 3. S14 correction: the exact global pair overlap is a trapezoid

Put `D_u={t:||ut||<=1/13}`.  Write `x_1=ga`, `x_2=gb` with `(a,b)=1` and,
without loss, `a<=b`.  Multiplication by `g` preserves Haar measure, so the
answer depends on the reduced projective pair `(a:b)`, not on the raw scale.
Define

```text
T(z)=sum_(m in Z) (z-|m|)_+,
psi(theta)=theta(1-theta),       0<=theta<1.
```

Then the exact formula is

```text
mu(D_x1 intersect D_x2)
 = [T((a+b)/13)-T((b-a)/13)]/(ab)                      (3.1)
 = 4/169
   +[psi({(a+b)/13})-psi({(b-a)/13})]/(ab).            (3.2)
```

Indeed a tooth centred at `j/a` and one centred at `k/b` have centre offset
`m/(ab)`, where `m=bj-ak`.  Coprimality makes the tooth pairs bijective with
`m mod ab`.  Since `2(a+b)<13ab`, every contributing residue has a unique
least integer representative and no contributing overlap wraps around the
residue circle.  The overlap is a **trapezoid**, not the triangle used in the
original S312 draft: after multiplication by `ab` its length is

```text
((a+b)/13-|m|)_+ - ((b-a)/13-|m|)_+.
```

Summing gives (3.1).  If `r=floor(z)` then
`T(z)=(2r+1)z-r(r+1)=z^2+psi({z})`; subtracting the two squares gives
`4ab/169` and proves (3.2).  Equivalently, with
`Q(c)=r(13-r)` for `r=c mod 13` in `{0,...,12}`,

```text
mu(D_x1 intersect D_x2)
 = 4/169 + [Q(a+b)-Q(b-a)]/(169ab),                    (3.3)
|mu-4/169| <= 42/(169ab).
```

The error has both signs.  Exact counterexamples to the original lower bound
are

```text
mu(D_6 intersect D_7)=2/91 < 4/169,
mu(D_5 intersect D_9)=4/195 < 4/169.
```

Thus `(a+b)^2/(169ab)` is not the leading term; it results from omitting the
containment plateau when the narrower tooth lies inside the wider one.
Near-equality is not a global minimizer theorem.  The sign is the mod-13
sawtooth `Q(a+b)-Q(b-a)`, and the restricted overlap inside a prefix safe set
has an additional endpoint-correlation term.

There is an exact scale-normal estimate for that term.  If `E` is a union of
`c_E` intervals and `B_(a,b)=D_a intersect D_b`, then

```text
|mu(E intersect D_(ga) intersect D_(gb))
   -mu(E)mu(B_(a,b))| <= 2c_E/g.                        (3.4)
```

On every full cell `[j/g,(j+1)/g)`, multiplication by `g` maps the pullback of
`B_(a,b)` linearly onto `B_(a,b)`, so its mass is exactly
`mu(B_(a,b))/g`.  Each component of `E` has at most two partial boundary
cells; bounding their actual and expected contributions by the cell length
proves (3.4).  Thus common scale `g`, reduced ratio `(a:b)`, and the mod-13
sawtooth are the natural edge coordinates.  Raw speeds alone are not.

## 3.5 The exact node-coloured Hunter defect algebra

The correction above admits a useful exact reorganization.  Let `e=mu(E)`
and let `x_1,...,x_m` be the remaining comb frequencies.  Give vertex `i` the
single-comb anomaly

```text
s_i = mu(E intersect D_(x_i)) - 2e/13.                  (3.5)
```

For every edge `ij`, write `x_i=g_ij a_ij`, `x_j=g_ij b_ij` with the reduced
pair coprime, and define

```text
h_ij = [Q(a_ij+b_ij)-Q(b_ij-a_ij)]/(169 a_ij b_ij),
eta_ij = mu(E intersect D_(x_i) intersect D_(x_j))
         - e mu(D_(a_ij) intersect D_(b_ij)),
c_ij = e h_ij + eta_ij.                                 (3.6)
```

Here `h_ij` is the projective mod-13 defect from independent pair density,
while `eta_ij` is the endpoint/pullback defect and satisfies
`|eta_ij|<=2c_E/g_ij`.  Substitution into Hunter--Kounias gives the exact
decomposition of its lower bound on the uncovered part of `E`:

```text
L_H(E;x_1,...,x_m)
 = e - sum_i mu(E intersect D_(x_i))
     + max_T sum_(ij in T) mu(E intersect D_(x_i) intersect D_(x_j))
 = (165-22m)e/169 - sum_i s_i + MST(c),                 (3.7)
```

where `MST(c)=max_T sum_(ij in T)c_ij`.  In particular, at the seven-comb
wall,

```text
L_H = 11e/169 - sum_i s_i + MST(c).                     (3.8)
```

Thus the residual is not a statistic of seven raw speeds.  It is a coloured
complete graph with vertex colours `s_i` and edge colours
`(a_ij:b_ij,g_ij,h_ij,eta_ij)`, evaluated by a tropical spanning-tree
character.  This is the precise sense in which pair data must remain joined
to its endpoint incidences.

There is also an exact recursive classification of the edge information that
the Hunter evaluator uses.  Let `lambda_1>...>lambda_q` be the distinct values
of `c_ij`, let `F_l={ij:c_ij>=lambda_l}`, and put
`r_l=m-kappa(F_l)`, where `kappa` is the number of connected components of
the threshold graph and `r_0=0`.  Kruskal's theorem gives

```text
MST(c)=sum_(l=1)^q lambda_l (r_l-r_(l-1)).               (3.9)
```

Only connectivity-increasing edge levels contribute, but which levels do so
depends on their full incidence pattern.  For example, a connected graph of
nonnegative `c_ij` edges implies `MST(c)>=0`; it does not follow from the
number or average of nonnegative edges.  Equations (3.7)--(3.9) turn the open
uniform step into a finite-type projective edge-classification problem once
the rational-prefix periodic tables for `s_i` and `eta_ij` are fixed.

## 3.6 The defects form a restricted covariance algebra

There is a canonical all-packet completion of the coloured graph.  Put
`p=2/13` and define centered functions on the circle by

```text
F_E=1_E-e,                  F_i=1_(D_(x_i))-p.
```

Then the vertex and global projective defects are ordinary correlations:

```text
s_i=<F_E,F_i>,              h_ij=<F_i,F_j>.                (3.9a)
```

In particular, for every packet, not only a coherent common-scale packet, the
matrix with off-diagonal entries `h_ij` and diagonal `22/169` is positive
semidefinite.  Its rank need not be six.

The additional prefix information is the centered third moment

```text
theta_(Eij)=integral F_E F_i F_j.
```

Expanding the indicators gives the exact ANOVA identities

```text
eta_ij=p(s_i+s_j)+theta_(Eij),
c_ij=e h_ij+p(s_i+s_j)+theta_(Eij).                         (3.9b)
```

Equivalently, the restricted Gram matrix

```text
G^E_ij=integral_E F_i F_j
      =c_ij-p(s_i+s_j),                         i!=j,
G^E_ii=p(1-p)e+(1-2p)s_i
      =22e/169+9s_i/13                                      (3.9c)
```

is positive semidefinite.  Thus valid node/edge colours lie in a constrained
Gram cone.  For every pair,

```text
|c_ij-2(s_i+s_j)/13|^2
 <=(22e/169+9s_i/13)(22e/169+9s_j/13),                     (3.9d)
```

and all higher principal-minor inequalities also hold.  This is a stronger
structural statement than the separate error bound (3.4): `s_i` is the
prefix--comb covariance, `h_ij` the comb--comb covariance, and the genuinely
new edge stalk is the prefix--comb--comb correlation `theta_(Eij)`.

For a tree `T`, (3.9b) also rewrites the defect evaluator as

```text
L_H=(165-22m)e/169
    +max_T [sum_(ij in T)(e h_ij+theta_(Eij))
             +sum_i(2 deg_T(i)/13-1)s_i].                  (3.9e)
```

Node colours therefore couple to tree degree; erasing node--edge incidence
loses even this covariance form.

## 3.7 The projective sawtooth is a six-channel odd kernel

The numerator in (3.6) has more structure than its sign table suggests.  On
`Z/13`, put

```text
H_(r,s)=Q(r+s)-Q(r-s).                                      (3.10)
```

Let `C_Q` be circular convolution by the even function `Q` and let
`(Rf)(s)=f(-s)`.  Then

```text
H=C_Q(R-I).                                                 (3.11)
```

Consequently `H` is symmetric, every row sums to zero, and it kills the
seven-dimensional reflection-even subspace.  On the odd basis
`u_i=e_i-e_(-i)`, `1<=i<=6`, its Gram matrix is

```text
u_i^T H u_j = 8 min(i,j)(13-2 max(i,j)).                    (3.12)
```

The six normalized leading principal minors are

```text
11, 117, 1183, 10985, 85683, 371293,
```

so the odd restriction is positive definite.  Hence `H` is positive
semidefinite of rank six and its kernel is exactly the even sector.
Equivalently,

```text
H_(r,s)=2 sum_(k=1)^6
  sin(2 pi k r/13) sin(2 pi k s/13) / sin(pi k/13)^2.       (3.13)
```

This follows directly by Fourier transforming `Q`; (3.12) and its integer
minors give a rational certificate independent of trigonometric numerics.

If a packet has one common scale `x_i=g a_i` and the `a_i` are pairwise
coprime, its global projective defects are therefore off-diagonal entries of a
six-dimensional Gram matrix:

```text
h_ij=<z(a_i),z(a_j)>,
z_k(a)=sqrt(2) sin(2 pi k a/13)/(13a sin(pi k/13)).          (3.14)
```

This is a sufficient coherent-reduction condition for the explicit rank-six
completion, not a characterization of every packet.  When pairwise gcd
reductions disagree, the general `L^2` Gram completion from (3.9a) still
exists, but these six residue coordinates do not assemble directly into one
packet.  Even under (3.14), the restricted credits `c_ij` need not be positive
semidefinite: the third moments in (3.9b) are additional edge data.  Thus the
residue sawtooth is a six-channel reflection-odd current, not an arbitrary
13-by-13 sign table, but it is only one layer of the Hunter state.

There is also a sharp all-packet scalar floor.  For every reduced coprime pair,

```text
h_ij >= -11/1014,                                           (3.14a)
```

with equality at projective ratio `1:12`.  Indeed the numerator in (3.3) is
at least `-42`; this proves the claim when `a_ij b_ij>=23`, and the 31 coprime
pairs with product at most 22 have unique minimum
`[Q(a+b)-Q(b-a)]/(ab)=-11/6` at `(1,12)`.  A seven-vertex tree has six edges,
so for every speed packet

```text
11/169+MST(h)>=0.                                           (3.14b)
```

Thus projective defects alone can consume, but cannot cross, the ideal
seven-comb margin if every edge is bounded independently.  This scalar-floor
argument is deliberately non-strict: `1:12` edges can occur along
multiplicative chains.  The complete graph nevertheless forces alternative
ratios, and the next lemma converts that consistency into strict uniform
progress.

## 3.7a The ten bad ratios force a strict projective tree margin

Put

```text
lambda_0=-3/676.
```

For a reduced pair `a<=b`, the strict inequality `h_(a,b)<lambda_0` is
equivalent to

```text
[Q(a+b)-Q(b-a)]/(ab)<-3/4.                                (3.14c)
```

The numerator is at least `-42`, so only `ab<=55` can occur.  The complete
coprime table is

```text
reduced ratios                 169 h_(a,b)
1:9, 2:9, 4:9                 -10/9
1:10, 3:10                    -7/5
1:11, 2:11                    -18/11
3:11                           -28/33
1:12                           -11/6
1:25                           -22/25.                    (3.14d)
```

Include both orientations of these ten ratios.  Their exact multiplication
table has no `r,s,r/s` all in the set.  Therefore, on any set of distinct
positive integer speeds, the graph whose edges satisfy `h<lambda_0` is
triangle-free: a bad triangle, ordered as `x<y<z`, would put the two ratios
`y/x,z/x` and their quotient `z/y` in the table.

Let `G_+` be the complementary graph of edges with `h>=lambda_0`.  It has at
most two connected components, since representatives of three components
would form a bad triangle.  If it is connected, one of its spanning trees
gives

```text
MST(h)>=6lambda_0=-9/338.                                  (3.14e)
```

If it has two components, good spanning trees inside them use five edges.
The cross cut has at least six pairs.  The unique minimum edge value is
`m_1=-11/1014`, attained exactly at ratio `1:12`; such edges form disjoint
subgraphs of the chains `x--12x` and have maximum degree two.  A complete
cross cut on seven vertices has a vertex of degree at least four, so not all
its edges can have value `m_1`.  The second distinct global edge value is

```text
m_2=-18/1859,
```

attained at `1:11` and `2:11`.  The same `-42` cutoff reduces its verification
to `ab<=25`.  Hence some cross edge has value at least `m_2`, and

```text
MST(h)>=m_2+5lambda_0=-237/7436.                           (3.14f)
```

The disconnected bound is the weaker of (3.14e)--(3.14f), so every seven
distinct positive integer speeds satisfy

```text
11/169+MST(h)>=19/572>0.                                  (3.14g)
```

This is a strict all-packet projective theorem, not a common-scale statement.
It does not say `MST(h)>=0`: the exact packet
`{4,9,21,32,70,170,189}` has all 21 edges negative and
`MST(h)=-505/1447992`.  The finite ratio table, both edge minima, and this
packet are independently checked in the projective-MST referee.

## 3.8 Exact recursion and the edge-order tournament

Write `tau(G)=MST(c)` for a finite weighted graph, with `tau=-infinity` when
no spanning tree exists and `tau(K_1)=0`.  Splitting spanning trees according
to whether they use a nonloop edge `e` gives tropical deletion--contraction:

```text
tau(G)=max(tau(G without e), c_e+tau(G/e)).                 (3.15)
```

Parallel labelled edges created by contraction must remain separate when
tree multiplicity or tie structure matters; for the value alone they may be
merged by taking their maximum.  A bridge makes the deletion branch
`-infinity`, and a loop is simply deleted.

There is also an exact node-insertion recursion.  Add a new comb vertex `v`
to the complete weighted graph on `V`, and write `d_i=c_(iv)`.  Then

```text
tau(G+v)=max_(pi partition of V)
  [sum_(B in pi) tau(G[B]) + sum_(B in pi) max_(i in B)d_i]. (3.16)
```

Indeed, deleting `v` from a spanning tree leaves a forest whose components
are the blocks `B`; each block has exactly one edge back to `v`.  Conversely,
optimal trees inside the blocks plus one maximizing `v`-edge per block form a
spanning tree.  For a noncomplete graph, restrict to blocks for which `G[B]`
is connected and which meet the star of `v`.

Put

```text
Delta_G(pi)=tau(G)-sum_(B in pi)tau(G[B]).
```

Then the insertion surplus is

```text
Gamma_v(G)=tau(G+v)-tau(G)
 =max_pi [sum_(B in pi) max_(i in B)d_i-Delta_G(pi)],        (3.17)
```

and the Hunter defect updates exactly by

```text
L_(m+1)=L_m-22e/169-s_v+Gamma_v(G).                         (3.18)
```

The full induced-tree profile `S -> tau(G[S])`, equivalently its partition-
defect profile, is therefore a recursively sufficient node colour for adding
one comb at fixed `E`.  This is not a finite-state or minimality theorem, and
if the prefix set changes then `e`, `s`, and `c` must also be transported.

Tournament Analysis becomes exact only after changing its vertex set and
retaining a sidecar.  Put the graph edges at tournament vertices, orient
`e -> f` when `c_e>c_f`, and declare a fixed endpoint gauge inside ties.  A
tie Hamiltonian path is a compatible descending order `e_1,...,e_N`.  With
`rho` the graphic-matroid rank, define

```text
k_j=rho({e_1,...,e_j})-rho({e_1,...,e_(j-1)}) in {0,1}.     (3.19)
```

Greedy accepts exactly the edges with `k_j=1`, so

```text
sum_j k_j=m-1,            tau(G)=sum_j c_(e_j)k_j.          (3.20)
```

Individual `k_j` values inside a tie level depend on its gauge.  The invariant
block count at level `lambda` is

```text
k_lambda=rho(E_(>=lambda))-rho(E_(>lambda)),                (3.21)
```

which is the rank jump already used in (3.9).  After contracting the
components formed by higher-credit edges, a level-`lambda` edge that becomes a
loop lies in no maximum tree, a bridge of the level quotient lies in every
maximum tree, and every other nonloop level edge is optional.  This gives an
exact recursive three-colour classification of edge roles.

The edge-order tournament plus the graphic circuit/rank sidecar preserves the
Hunter evaluator.  The tournament alone does not: away from ties it is always
transitive, and attaching the same ordered values to a different incidence
pattern changes (3.19).  Moreover, the present `(lambda,k_lambda)` vector
evaluates `tau(G)` but does not determine future node insertions; (3.16) needs
the induced-subgraph profile.  A useful coarse signed grade is

```text
delta_+(c)=m-1-rho({e:c_e>0})=kappa({e:c_e>0})-1,           (3.22)
```

the exact number of nonpositive edges Kruskal must still add.  Positive
credits contain a spanning tree exactly when every nontrivial vertex cut has
a positive crossing credit.

## 3.9 Every common-dilate family closes asymptotically

The strict projective margin gives a larger uniform subcase of the open
radius-seven step.  Let `E` be a fixed union of `c_E` circle intervals, of
measure `e`, and let `x_1,...,x_7` be distinct positive integer speeds.  Among
the trees maximizing `sum_T h_ij`, choose one minimizing the reciprocal-gcd
cost, and put

```text
Gamma_h(x)=min_(T: sum_T h=MST(h))
             sum_(ij in T) 1/gcd(x_i,x_j).                 (3.23)
```

The full-cell/partial-cell proof of (3.4), applied to one comb and to every
pair, gives

```text
|s_i|<=2c_E/x_i,
|eta_ij|<=2c_E/gcd(x_i,x_j).                               (3.24)
```

Evaluate `MST(c)` on the tree chosen in (3.23), then use (3.14g).  The result
is the all-packet sufficient criterion

```text
L_H >= 19e/572
       -2c_E[sum_i 1/x_i+Gamma_h(x)].                      (3.25)
```

This formula still sees incoherence: a large raw speed does not make an edge
endpoint error small when its pair gcd is small.  It nevertheless closes every
common-dilate ray.  If

```text
x_i=G y_i,             i=1,...,7,                          (3.26)
```

with arbitrary distinct positive integers `y_i`, then every pair gcd is at
least `G`.  After sorting, `y_(i)>=i`, so no pairwise-coprime hypothesis is
needed and

```text
sum_i 1/x_i+Gamma_h(x)
 <=(H_7+6)/G=1203/(140G),
L_H >= 19e/572-(1203/70)c_E/G.                            (3.27)
```

Hunter therefore proves noncoverage whenever

```text
G > (1203/70)(572/19)c_E/e.                               (3.28)
```

This leaves only finitely many common dilates for every fixed projective
packet.  The earlier rank-six sine completion remains a useful exact
description when the reduced coordinates are pairwise coprime, but it is not
needed for (3.27).  The full frontier consists of the finite initial dilates
and packets whose reciprocal-gcd tree cost in (3.25) does not decay.

## 4. The periodicity lemma (the finite table for E-restricted masses)

For a prefix safe set E with rational endpoints and a scale-one lift
x = r + 13h: **x·(|E ∩ D_x| − (2/13)μ(E)) is EXACTLY periodic in h** (period
dividing an explicit Λ(E, r); verified: prefix {1,...,5}, r = 6: period 60,
exact). The per-comb data of the infinite radius-7 chart is a finite table of
rationals. *Proof:* the anomaly depends only on the positions of E's endpoints
in the comb's tooth-coordinate, i.e. on x·(endpoint) mod 1, which is periodic
in h with period = denominator/gcd data. ∎

## 5. Pilot verification (prefix {1,2,3,4,5}, exact ℚ)

μ(E) = 7/15, 10 components. Four random radius-7 packets (lifts h ∈ [2,40]):
Hunter bound = +0.033, +0.034, +0.036, +0.038 — ALL COERCIVE (non-coverage
PROVED per packet by one inequality; actual uncovered ≈ 0.13–0.14, i.e. the
independence prediction (11/13)⁷·μE = 0.145 is accurate to 9%). Consecutive
packets {499..505}, {32..38}: Hunter −0.001, −0.017 while actual uncovered =
0.055, 0.040 — still far from tight.  These are exact finite witnesses that
the functional can succeed and fail.  In view of (3.3), they do not prove
that consecutive packets are the unique failure locus.

The S15 node/edge replay checks all six packets directly with exact interval
arithmetic.  The four random packets have connected positive-credit graphs,
and their Kruskal trees use six positive edges.  In both consecutive packets
all 21 endpoint discrepancies and all 21 total credits are negative, even
though the projective terms `h_ij` are positive on 14 and 18 edges respectively.
Their failure is therefore prefix-coupled, not visible in the global sawtooth
sign.  Conditional-overlap tournaments on the seven proof events are
transitive in all six packets (score histogram `0,...,6`, singleton SCCs, no
directed triangles, one Hamiltonian path); this fingerprint occurs on both
sides of Hunter positivity and is a deliberately lossy quotient.  The
historical high consecutive stress packet `{499,...,505}` has residues
`5,...,11`, not the scale-one chart `6,...,12`, and is retained with that scope.

## 6. What remains for a radius-7 closure (named residuals)

1. **Incoherent endpoint/gcd cost.**  A bound on each edge by raw speed alone
   is false.  On the legal scale-one ray

   ```text
   g=1+13k,       x=6g=6+13(6k),       y=7g=7+13(7k),
   ```

   the reduced pair is always `(6,7)` and
   `mu(D_x intersect D_y)=2/91<4/169`.  For every fixed finite union `E` of
   intervals,

   ```text
   mu(E intersect D_(6g) intersect D_(7g))
      -> mu(E) 2/91,
   ```

   because `D_(6g) intersect D_(7g)` is the pullback of
   `D_6 intersect D_7` by multiplication by `g`.  Hence the originally
   proposed `C(E)/min(x,y)` deficit from the ideal `4mu(E)/169` has a
   nonzero limiting left side and a zero right side.  The ten-ratio lemma now
   resolves the joined projective tree uniformly, and (3.27) closes this
   common-dilate example for large `g`.  What remains is the endpoint field
   when the reciprocal-gcd cost `Gamma_h(x)` in (3.25) does not decay.
2. **Boundary third moments at fixed reduced ratio.**  For `x=ga,y=gb`, use the
   proved split of the restricted edge weight into

   ```text
   mu(E) mu(D_a intersect D_b)
   + [endpoint discrepancy of the g-fold pullback].
   ```

   Equation (3.4) bounds the second term by `2c_E/g`; for rational endpoints
   its scaled value is periodic.  The projective margin is already `19/572`.
   What remains is a finite operation-stable language for the joined node
   covariances and third moments on small-gcd edges, including how its
   maximizing tree changes across endpoint cells.
3. **Correlated/AP-window residual.**  The two consecutive pilot packets evade
   the tree certificate but remain far from covering.  Once the projective
   densities and boundary errors are separated, use an AP-window/mechanical-
   word argument (THM-778) on any genuinely consecutive residual rather than
   attributing it to a nonexistent global Bezout minimum.

**The surviving replacement potential is Hunter's weighted-tree functional.
Its ideal-density coefficient crosses the first-moment wall at seven combs,
and the exact pilots show that it can certify real packets.  Its uniform input
is not a scalar raw-speed deficit: it is the reduced-ratio, mod-13 sawtooth,
and endpoint-discrepancy packet above.**

## 7. Exploratory cluster pilots from S312 (not an asymptotic theorem)

The companion cluster script adds three useful finite observations on
`E_[5]`, but its original “complete asymptotic schema” language is not proved.
In particular, the script labels a calibration as such and estimates one
integral by midpoint sampling; neither step is an exact all-packet argument.

1. At `x=500`, differences `d=1,2,3,5` give restricted overlaps below the
   independent baseline, while the sampled `d=8,13,40,200` rows are closer to
   it.  This supports a beat-window diagnostic.  It does **not** prove a
   threshold at `d=8`: equations (3.3)--(3.4) show that raw difference alone
   cannot determine the edge weight.
2. Three displayed multi-cluster packets have positive Hunter certificates
   when a chosen tree uses inter-cluster edges.  This is a three-packet exact
   pilot, not a proof that every large raw separation gives baseline overlap;
   individual projective terms can remain negative.  Section 3.9 now absorbs
   every sufficiently large common-dilate ray, but not arbitrary small-gcd
   endpoint packets.
3. The seven-consecutive-speed avoid function has positive sampled integral
   over `E_[5]` (about `0.0558`) and the recorded spot values expose narrow
   denominator-seven windows.  Turning this into a uniform floor still needs
   exact integration over its Farey cells and a quantitative skew-
   equidistribution error.

Thus the cluster experiment names plausible mechanical-word sublemmas but
does not close the asymptotic radius-seven chart or reduce the remaining work
to finite bookkeeping.  Its stored numbers remain reproducible finite data.

## 8. S16 adaptive graphic-rank localization

For a measurable prefix-safe set `E`, let

```text
A_i=E intersect D_(u_i),             S(t)={i:t in A_i}.
```

The graphic-matroid rank of the clique induced by the active set is
`r(t)=max(|S(t)|-1,0)`.  Pointwise,

```text
1_(S(t)=empty)=1-|S(t)|+r(t),
```

and hence

```text
mu(E minus union_i A_i)
 =mu(E)-sum_i mu(A_i)+integral_E r(t)dt.                 (8.1)
```

This is an exact all-order identity.  A fixed spanning tree contains at most
`r(t)` edges induced by `S(t)`, so Hunter is precisely the lower bound obtained
by replacing the spatially adaptive graphic basis in (8.1) by one global tree.
The discarded nonnegative quantity is the tree-adaptivity gap

```text
Gamma_E=integral_E r(t)dt-max_T sum_(ij in T)mu(A_i intersect A_j). (8.2)
```

There is a monotone exact localization hierarchy.  For a finite measurable
partition `P` of `E`, let `MST_C` be the maximum spanning-tree weight from the
pair overlaps restricted to a cell `C`, put

```text
H_C=mu(C)-sum_i mu(C intersect A_i)+MST_C,
L(P)=sum_(C in P) max(0,H_C).                            (8.3)
```

Then

```text
0<=L(P)<=mu(E minus union_i A_i).
```

If `P'` refines `P`, then `L(P')>=L(P)`; on the Boolean activation-atom
partition equality holds.  Indeed `MST(w_1+w_2)<=MST(w_1)+MST(w_2)` under a
cell split, and on an atom with nonempty active set `S` a tree takes exactly
`|S|-1` positive internal edges, whereas the empty atom contributes its full
mass.  This proves the claims without a generic-position or independence
hypothesis.

The exact-Fraction audit on the twelve first proper Hamming-one packets

```text
([12] minus {r}) union {r+13},          1<=r<=12,
```

uses the five smallest speeds as `E` and the other seven as needles.  One
global Hunter tree is positive in only `3/12` rows; literal-component-local
Hunter and a rooted six-comb localization are positive in `12/12`.  For the
row missing `6`, the three values are respectively

```text
-634/9009,       1/8151,       1/286,
```

and the rooted value equals the exact uncovered mass.  The bank was already
closed by stronger Hamming recursion: this is a new general certificate and
method test, not global sporadic-branch closure.

Tournament Analysis uses the 21 pair-overlap obligations rather than runners.
Every component order is transitive, yet the missing-`6` row realizes five
local edge orders and five local Kruskal trees, with 77--135 flips from the
global order.  Thus a tournament retains local Kruskal priority only together
with the original `K_7` edge-incidence sidecar; activation atoms remain the
exact carrier.  In the Kakeya language the periodic combs are fixed-intercept
needles, and (8.1)--(8.3), rather than a Euclidean dimension estimate, are the
part of incidence geometry that transfers. ∎
