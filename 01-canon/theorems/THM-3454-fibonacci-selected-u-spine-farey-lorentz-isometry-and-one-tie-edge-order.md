---
id: THM-3454
title: "Fibonacci-selected U-spine Farey--Lorentz isometry and one-tie edge order"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The fixed-cusp
  fraction leaf chain, parabolic Berggren U-spine, and Lorentz chord metric
  have the same integral separation.  Four consecutive Fibonacci spine indices
  (at rooted depths one lower) give
  a transitive rooted T4, while their six edge separations form a five-level
  total preorder whose only stable tie is the adjacent pair 02/23.  Two marked
  edge equalities characterize additive four-windows, and the inherited
  Cassini/Pell unit then characterizes the Fibonacci windows.  The quadratic
  branch labels sampled at Fibonacci spine indices obey a minimal order-six linear
  recurrence.  No Farey-graph, full-tree, LRC, or Jacobian equivalence follows.
source: codex-2026-08-15-fibonacci-u-spine-distance-transplant
audit: >
  independent algebra, chart-sign, small-index, recurrence, Pell, tournament,
  incidence-orbit, Farey-distance, source-dependence, and repo-novelty audits;
  exact normal/optimized/stored replay and dependency/hash gates
depends_on:
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
related:
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost
  - THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase
  - THM-3379-fibonacci-ray-local-t4-bit-is-mod3-boolean-projection
  - THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses
script: 04-computation/fibonacci_u_spine_farey_lorentz_one_tie_thm3454.py
output: 05-knowledge/results/fibonacci_u_spine_farey_lorentz_one_tie_thm3454.out
script_sha256: fadfe97113ea1d38472face781a39f3c5f615f9aea093025269b3c72e03372fa
output_sha256: 35c964e0635e07f16146b3c439d0eee32a0d465009ceea17e1d71729d3f68b07
semantic_sha256: ff0262595d3dea442d874294e5166c95c47c5b7514375d8b13ab604cdbd6b4f9
hash_basis: LF-normalized bytes
---

# THM-3454 -- Fibonacci-selected U-spine Farey--Lorentz isometry and one-tie edge order

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is a repository synthesis and corollary theorem, not a
literature-priority claim.  Its new content is the joined metric statement,
the six-distance preorder and sharp boundary, the marked recurrence converse,
and the Fibonacci-sampled branch-label recurrence.  The ambient Farey,
Lorentz, Berggren, and Pell mechanisms are inherited.

## 1. Inheritance and validity gate

The closest proved mechanisms are:

1. [THM-3333](THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md),
   which proves the Lorentz cross-determinant identity for the ordered
   Gaussian-square lift;
2. [THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
   which identifies the entire fixed-cusp fraction fan and parabolic Berggren
   `U`-spine; and
3. [THM-3339](THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md),
   which proves the Cassini/Pell classification of consecutive Fibonacci
   parameters and a different, product-weighted strict `T6`.

The canonical hostile is
[THM-3345](THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost.md):
away from one branch, Berggren ancestry costs are source-dependent.  The two
edges of one XOR colour can have different path lengths.  Therefore the
metric equality below is a same-`U`-ray theorem, not a global tree law.

[MISTAKE-394](../MISTAKES.md#mistake-394-2026-08-15-fibonacciberggren-17-adic-torsor-scope----a-parameter-norm-square-was-misnamed-and-the-tied-root-entered-a-t6-support-claim)
is the corrected near miss: a tied root may not be silently admitted to a
tournament.  Here ties are part of the theorem and the six-cost object is not
called a tournament.

The least-used decisive sidecars are:

- the marked spine index (equivalently rooted depth plus one), because
  pairwise distances forget a common translation; and
- edge incidence in `K4`, because `S4` preserves whether two edges meet.

The exact objects are:

| role | declaration |
|---|---|
| four vertices | four marked nodes on one rooted `U`-ray |
| six vertices | the six labelled unordered pairs `01,02,03,12,13,23` |
| pairwise observable | rooted tree separation, equivalently fan cross-determinant magnitude |
| strict orientation | smaller separation to larger separation |
| equality | retained as a missing strict comparison, or as one bidirected weak comparison |
| preserved target | the five distance levels and the marked additive seams |
| lost data | common index/depth translation, Cassini sign, product order, owner, phase, and off-ray ancestry word |
| required sidecar | one absolute spine index (equivalently one rooted depth) plus the labelled edge-incidence structure |

## 2. The fixed-cusp leaf chain is the `U`-spine

For every integer `t>=1`, define the reduced fraction vector, its swapped
Gaussian spinor, and its Pythagorean point by

```text
r_t=(t,t+1),                 rho_t=t/(t+1),
s_t=(t+1,t),                 P_t=Phi(s_t),                (1)
```

where the repo chart is

```text
Phi(m,n)=(m^2-n^2,2mn,m^2+n^2).                          (2)
```

The swap in `(1)` is load-bearing: `r_t` is the numerator/denominator vector,
whereas `s_t` is the positive Gaussian chart.  Directly,

```text
P_t=(2t+1,2t(t+1),2t(t+1)+1).                            (3)
```

The Berggren matrix

```text
U=[1 -2 2; 2 -1 2; 2 -2 3]                              (4)
```

satisfies `U P_t=P_(t+1)` and `P_1=(3,4,5)`.  Hence

```text
P_t=U^(t-1)(3,4,5).                                      (5)
```

This is THM-3334's parabolic spine with its index shifted by one.

Every `rho_t` lies in `(0,1)` and is Farey-adjacent to the fixed cusp `1/1`,
because

```text
det(r_t,(1,1))=-1.                                       (6)
```

For `1<=s<t`,

```text
det(r_s,r_t)=s-t,
det(s_s,s_t)=t-s.                                        (7)
```

Thus the graph on the leaves `{rho_t:t>=1}` with Farey adjacency inherited
from the full Farey graph is an infinite path: two leaves are adjacent iff
their indices differ by one.  Let `d_leaf` be distance in this leaf path.

On the rooted Berggren tree, `P_s` is the ancestor of `P_t` on the same
`U`-ray, at depths `s-1` and `t-1`.  Therefore

```text
d_B(P_s,P_t)=t-s.                                        (8)
```

Finally THM-3333 gives

```text
<Phi(u),Phi(v)>_L=2 det(u,v)^2,                          (9)
```

with `<X,Y>_L=X_C Y_C-X_A Y_A-X_B Y_B`.  Combining
`(7)--(9)` proves the exact joined metric law

```text
d_leaf(rho_s,rho_t)
=d_B(P_s,P_t)
=|det(r_s,r_t)|
=sqrt(<P_s,P_t>_L/2)
=t-s.                                                    (10)
```

Equation `(10)` is **not** a Farey-graph isometry.  In the full Farey graph,
nonadjacent leaves have distance two through `1/1`.  For example, `rho_1=1/2`
and `rho_4=4/5` have cross-determinant and `U`-distance three but Farey graph
distance two.

## 3. The hidden quadratic label and its Fibonacci sampling

Write `P_t=(a_t,b_t,c_t)`.  To distinguish the user's parameter index from
THM-3334's depth index, denote the scalar here by lowercase `q_t`:

```text
q_t=a_t^2+2=2c_t+1=4t(t+1)+3=(2t+1)^2+2
   =Q_(t-1) in THM-3334's depth notation.                (11)
```

The sequence beginning at `t=1` is

```text
11,27,51,83,123,171,227,...                              (12)
```

and has constant second difference eight.  The underlying hypotenuse norm

```text
c_t=t^2+(t+1)^2=2t^2+2t+1                               (13)
```

is the discriminant-`-4` quadratic carrier.  The affine label `q_t=2c_t+1`
itself has polynomial discriminant `-32`; these two discriminants must not be
conflated.  The exceptional scalar `2` is not a member of `(11)`.

Now sample `(11)` at Fibonacci spine indices.  With `F_0=0,F_1=1`, put

```text
R_n=q_(F_n)=Q_(F_n-1)=4F_n(F_n+1)+3,        n>=1.        (14)
```

Then

```text
R_0,R_1,...=3,11,11,27,51,123,291,731,1851,... .         (15)
```

Here `R_0=3` is the polynomial extension `q_0`; geometrically it would have
depth `-1` and is not a rooted positive Berggren triple.  The geometric
sequence starts at `n=1`, with its sole initial collision `R_1=R_2=11` coming
from `F_1=F_2=1`.

This is a sparse sample of `(12)`, not the same constant-second-difference
sequence.  With the polynomial boundary retained, it obeys the minimal
order-six recurrence

```text
R_n=4R_(n-1)-2R_(n-2)-6R_(n-3)
    +4R_(n-4)+2R_(n-5)-R_(n-6),             n>=6.        (16)
```

On the canonical geometric sequence `n>=1`, the first self-contained use of
the same recurrence is at `n=7`.

Indeed, if `alpha,beta` are the roots of `x^2-x-1`, then

```text
F_n^2=(alpha^(2n)+beta^(2n)-2(-1)^n)/5.                 (17)
```

Consequently `(14)` is a nonzero linear combination of the six distinct
modes

```text
alpha^2, beta^2, -1, alpha, beta, 1.                     (18)
```

Their characteristic polynomial is

```text
(x^3-2x^2-2x+1)(x^2-x-1)(x-1)
=x^6-4x^5+2x^4+6x^3-4x^2-2x+1.                         (19)
```

All six mode coefficients in `(14)` are nonzero, so no proper factor of
`(19)` annihilates the sequence.  This proves both `(16)` and minimality over
the rationals.

## 4. Four Fibonacci spine indices and six common separations

For `k>=4`, put

```text
(e_0,e_1,e_2,e_3)=(F_(k-1),F_k,F_(k+1),F_(k+2)),
X_i=P_(e_i),                                             (20)
```

Here `e_i` is the spine index.  The actual rooted depth of `X_i` is
`h_i=e_i-1`.  Thus all metric differences below are also rooted-depth
differences, while an additive recurrence written in the `h_i` coordinates
acquires the affine `+1` term made explicit in Section 6.

For `i<j`, define the common separation

```text
d_ij=d_B(X_i,X_j)
    =|det(r_(e_i),r_(e_j))|
    =sqrt(<X_i,X_j>_L/2)
    =e_j-e_i.                                            (21)
```

The four nodes are strictly ordered by rooted ancestry, so orienting an
ancestor toward a descendant gives the transitive `T4`

```text
X_0 -> X_1 -> X_2 -> X_3.                                (22)
```

The six edge separations are

| edge | exact cost |
|---|---:|
| `01` | `F_(k-2)` |
| `12` | `F_(k-1)` |
| `02` | `F_k` |
| `23` | `F_k` |
| `13` | `F_(k+1)` |
| `03` | `2F_k` |

For `k>=4`, positivity and the Fibonacci recurrence give

```text
d_01 < d_12 < d_02=d_23 < d_13 < d_03.                  (23)
```

Therefore the strict comparison relation on the six edge labels is the
transitive orientation of

```text
K_6 minus {{02,23}}.                                     (24)
```

It has `14` strict arcs.  Equivalently, `d(E)<=d(E')` is a five-level total
preorder; after deleting reflexive loops, its weak arc encoding is a
semicomplete digraph with the one digon `02<->23`.  Collapsing the tied pair
gives a transitive `T5`.  There are exactly two strict `T6` refinements,
obtained by orienting the tied pair either way, but neither refinement is
intrinsic.

### Sharp lower boundaries

| `k` | window | six-cost levels | consequence |
|---:|---|---|---|
| `1` | `(0,1,1,2)` | `0 < {1,1,1,1} < 2` | `P_0=(1,0,1)` is not a positive descendant; duplicates |
| `2` | `(1,1,2,3)` | `01 < {02,12,23} < {03,13}` | vertex duplicate; `11` strict arcs |
| `3` | `(1,2,3,5)` | `{01,12} < {02,23} < 13 < 03` | first rooted `T4`, but two cost ties; `13` strict arcs |
| `k>=4` | strictly increasing | `01 < 12 < {02,23} < 13 < 03` | exactly one tie; `14` strict arcs |

Thus `k=3` is the first four-node tournament and `k=4` is the first stable
one-tie six-edge preorder.  This explicitly respects MISTAKE-394.

## 5. Why this is not THM-3339's product `T6`

THM-3339 compares the six products `w_ij=e_i e_j`.  For `k>=3` they are all
distinct and form a genuine transitive `T6`; its only parity-changing
comparison is the opposite-edge pair `03/12`.

The present observable compares distances and forces equality on the
adjacent-edge pair `02/23`.  The `S4` action on the six edges of `K4`
preserves incidence.  Among the `15` unordered edge-pairs, the adjacent orbit
has size `12` and the opposite orbit has size `3`.  Hence

```text
0 of 24 K4 relabelings send {02,23} to {03,12}.           (25)
```

The two six-edge objects are not relabelings of one another.  Nor are they
THM-3364's local `T4` on a parent and its three children, or THM-3379's six
orders of three channels.  Cardinality alone supplies no identification.

## 6. The sharp marked recurrence converse

Let `x_i` denote positive spine-index coordinates, and put `h_i=x_i-1` for
their actual rooted depths.  Thus let

```text
0<x_0<x_1<x_2<x_3,          c_ij=x_j-x_i.                (26)
```

Then

```text
c_12=x_0  and  c_02=c_23
iff
x_2=x_0+x_1  and  x_3=x_1+x_2.                          (27)
```

Equivalently, in actual rooted-depth coordinates,

```text
c_12=h_0+1  and  c_02=c_23
iff
h_2=h_0+h_1+1  and  h_3=h_1+h_2+1.                     (27a)
```

The first equality in `(27)` is exactly `x_2-x_1=x_0`.  Substituting this
into `x_2-x_0=x_3-x_2` gives `x_3=x_1+x_2`.  The reverse implications are
immediate.

Under `(27)`, the six labelled costs, in the order
`01,12,02,23,13,03`, are

```text
x_1-x_0, x_0, x_1, x_1, x_0+x_1, 2x_1.                 (28)
```

Consequently the stable order `(23)` holds exactly when

```text
x_1<2x_0.                                                (29)
```

The boundary `x_1=2x_0` creates the additional tie `01=12` seen at `k=3`.

The appearance of `x_0` in `c_12=x_0` is essential.  Pairwise distances
recover only relative positions and are invariant under

```text
(x_0,x_1,x_2,x_3) -> (x_0+h,x_1+h,x_2+h,x_3+h).         (30)
```

Thus `(27)` uses one absolute spine-index sidecar, equivalently one absolute
rooted depth together with the known `+1` chart shift, not edge order alone.

Finally define the Cassini/Pell norm

```text
D(x_0,x_1)=x_1^2-x_0 x_1-x_0^2.                         (31)
```

THM-3339 proves by strict descent that, for positive `x_0<x_1`,

```text
|D(x_0,x_1)|=1
iff
(x_0,x_1)=(F_r,F_(r+1)) for a unique r>=2.               (32)
```

Equivalently,

```text
(2x_1-x_0)^2-5x_0^2=+-4.                                (33)
```

Combining `(27)` and `(32)` gives the exact classification

```text
(27) and |D|=1
iff
(x_0,x_1,x_2,x_3)=(F_r,F_(r+1),F_(r+2),F_(r+3)).        (34)
```

Adding `x_0>=2`, equivalently the strict lower comparison `c_01<c_12`, gives
the stable range `k>=4`.  In the notation `(20)`,

```text
D(F_(k-1),F_k)=(-1)^(k-1).                               (35)
```

### 6.1 Long-window gap transplant

The four-point converse is the first instance of an all-length statement.
Let `m>=2`, let

```text
0<x_0<x_1<...<x_m,          g_i=x_(i+1)-x_i.             (35a)
```

Then

```text
x_(i+1)=x_i+x_(i-1) for 1<=i<=m-1
iff
g_1=x_0 and g_i=g_(i-1)+g_(i-2) for 2<=i<=m-1.          (35b)
```

In words: a spine-index window is Fibonacci-recurrent exactly when its
consecutive gap sequence is Fibonacci-recurrent and one marked seam aligns
the two sequences.  In the forward direction `g_i=x_(i-1)`.  Conversely, `g_1=x_0`
first gives `x_2=x_1+x_0`; then simultaneous induction using the gap
recurrence gives

```text
g_i=x_(i-1),              x_(i+1)=x_i+x_(i-1).           (35c)
```

For rooted depths `h_i=x_i-1`, the same statement is affine:
`h_(i+1)=h_i+h_(i-1)+1`, and the marked seam is `g_1=h_0+1`.

For four vertices, `g_1=x_0` is `c_12=x_0`, while
`g_2=g_1+g_0` is exactly `c_23=c_02`.  Thus the adjacent-edge tie is the
first gap-recurrence equation, not a numerical accident.  Adding the Pell
unit `(32)` to `(35b)` classifies every finite consecutive Fibonacci window,
not only the four-point case.

## 7. Hostile anatomy and strongest survivors

| overclaim | minimal or sharp witness | first failed implication | survivor |
|---|---|---|---|
| cross-determinant is Farey graph distance | `rho_1,rho_4` | determinant `3`, graph distance `2` via `1/1` | leaf-chain and `U` distances equal determinant |
| one recurrence seam suffices | `(1,2,3,4)` | `c_12=x_0`, but `c_02!=c_23` | both seams are iff |
| the tie seam alone suffices | `(2,3,4,6)` | `c_02=c_23`, but `c_12!=x_0` | both seams are iff |
| recurrence forces Fibonacci | `(1,3,4,7)` | recurrence holds, `|D|=5` | add the Pell unit |
| stable order plus recurrence forces Fibonacci | `(3,4,7,11)` | stable one-tie order holds, `|D|=5` | add the Pell unit |
| six labelled costs plus `|D|` identify the index window | `(2,3,5,8)` and `(1,2,4,7)` | both have costs `(1,2,3,3,5,6)` and `|D|=1`; only the first is recurrent | retain one absolute index/rooted-depth coordinate and both seams |
| Pell plus stable weak order suffices even with `x_0>=2` | `(2,3,6,10)` | `|D|=1` and the weak order is stable, but recurrence fails | order is not equation data |
| same metric law holds across a Berggren fibre | THM-3345's equal-colour costs | source-dependent words have unequal lengths | restrict to one marked ray |

The sixth row is the clean quotient-loss witness: the edge metric forgets a
common translation, and the absolute Pell norm can accidentally survive that
translation.  The missing coordinate is an origin, not another orientation.

## 8. Boolean subsets and the harmonic-series loss wall

This elementary corollary records the user's subset viewpoint; no novelty is
claimed over the repo's existing Boolean/harmonic carriers.  For every
`A subset N_(>=1)`, put

```text
R(A)={rho_t:t in A},
B(A)={P_t:t in A},
H(A)={1/t:t in A}.                                       (36)
```

Each is a set-valued injection, and each preserves union, intersection,
complement relative to its ambient image, and symmetric difference.  Thus
every subset of the natural numbers is literally a labelled subset of the
fixed-cusp fraction fan, of the `U`-spine, and of the harmonic series.

The scalar harmonic valuation is not faithful.  Already

```text
sum_(t in {2}) 1/t = 1/2 = sum_(t in {3,6}) 1/t.         (37)
```

So the term subset is the Boolean carrier; its total mass is only a lossy
invariant.

A faithful scalar-valued repair is not one harmonic evaluation but the full
Dirichlet germ

```text
D_A(s)=sum_(t in A) t^(-s),                 real s>1.     (38)
```

If `A!=B`, let `m` be the least element of their symmetric difference.  After
choosing the sign of the coefficient at `m`,

```text
m^s(D_A(s)-D_B(s)) -> +1 or -1              as s->infinity,  (39)
```

because the remaining tail is bounded by
`sum_(n>m)(m/n)^s`, which tends to zero.  Hence the full Dirichlet germ
recovers the labelled subset exactly.  This is the same sidecar lesson as the
metric theorem: one scalar value loses the kernel; a sufficiently rich germ
restores it.

The Fibonacci lane is even sharper.  For `J subset {2,3,...}`, the labelled
subseries

```text
sum_(n in J) 1/F_n                                       (40)
```

always converges.

Starting at `n=2` removes the duplicated value `F_1=F_2=1`, so `(40)` is a
support subset rather than the indexed term multiset forbidden by
MISTAKE-209.  Retaining both indices would be a different, explicitly
labelled object.  Indeed `F_(n+2)>=2F_n`, so grouping even and odd indices
gives the uniform majorant

```text
sum_(n>=2) 1/F_n <= 2 sum_(m>=1) 2^(-(m-1))=4.           (41)
```

Hence the convergence/divergence bit collapses every Fibonacci-indexed
Boolean subset, even though the labelled term set remains faithful.

## 9. Exact companion

The exact companion checks, using only integral and rational arithmetic:

- all `40` fan points and all `780` pairs among them;
- the matrix word, fixed-cusp, determinant, Lorentz, null, `q`, and
  second-difference identities;
- the Fibonacci-sampled order-six recurrence through index `40` and its
  independent `6 x 6` Hankel determinant `393216`;
- every Fibonacci window for `2<=k<=20`, including strict-arc and linear-
  extension counts `11/12`, `13/4`, and `14/2` at the three boundaries;
- the complete `S4` adjacent/opposite pair orbits and all `24` relabelings;
- all `3060` increasing quadruples in `[1,18]` for the marked recurrence iff;
- `1666` recurrence-order boundary rows;
- every positive Pell pair through `500`, with exactly the twelve Fibonacci
  pairs listed in the transcript; and
- the translation, partial-recurrence, non-Pell, Farey-distance, harmonic,
  and convergence hostiles above.

Normal and optimized runs reproduce the stored transcript byte for byte after
LF normalization:

```text
python3 -B 04-computation/fibonacci_u_spine_farey_lorentz_one_tie_thm3454.py
python3 -B -O 04-computation/fibonacci_u_spine_farey_lorentz_one_tie_thm3454.py
```

The script freezes the LF-normalized hashes of THM-3333, THM-3334, and
THM-3339 and uses explicit exceptions rather than optimization-sensitive
assertions.

## 10. Scope and non-consequences

| field | exact boundary |
|---|---|
| source | one marked parabolic `U`-ray and its fixed-cusp leaf chain |
| target | integral line metric, Lorentz chord metric, and labelled edge preorder |
| preserved | spine-index/rooted-depth differences, edge incidence, additive seams, and `Q` labels |
| destroyed by edge metric | absolute spine index/rooted depth, Cassini sign, products, phase, owner, and off-ray word |
| needed sidecar | one absolute spine index (or rooted depth) and the labelled `K4` incidence structure |
| cheapest positive | the window `(2,3,5,8)` |
| cheapest translation hostile | `(1,2,4,7)` |

Nothing here identifies the full Berggren tree with the Farey or
Stern--Brocot tree.  It does not turn the six costs into an intrinsic `T6`,
does not identify the distance preorder with THM-3339's product `T6`, and
does not construct an LRC current, phase, owner, spectral closure, Keller
map, Jacobian counterexample, or degree-spectrum exclusion.  LRC(14) and
`JC(2)` remain open.
