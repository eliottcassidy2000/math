# Session140: proof-preserving representations, not numerical analogies

**Status.** Elementary statements below are **PROVED** as written, with
same-agent exact controls. Older mechanisms retain their original status and
credit. The explicit Collatz macro-tree gives a useful geometric obstruction;
exponential inverse-tree growth is classical, so no historical novelty claim
is made. Neither Collatz convergence nor a new Hadwiger–Nelson bound follows.
Source baseline: `00f55b6c70a4` (the complete pinned hash is in SOURCE_PINS.json).

**Concurrent credit.** The Collatz lane in `../collatz_hex/PROOF.md` proves
the stronger equal-five-step binary supply for every root prime to three.
The eight-step increasing odd hierarchy here is an independently derived
redundant control, not the session's sharp clock. The HN lane in
`../hadwiger_nelson/PROOF.md` supplies the induced continuum realization
and exact cycle-displacement lift used as complementary constructions.

## Immediate result: a fixed hexagonal mesh cannot carry exact Collatz locally

Write the **standard**, unaccelerated map as

```
T(n)=n/2 if n is even, and T(n)=3n+1 if n is odd.
```

There is no injective map of its positive-integer graph into either the
honeycomb vertices or the triangular lattice of hexagon centers with a
uniform bound on the lattice distance between the images of `n` and `T(n)`.
Thus allowing any fixed number of routing edges per arithmetic step does
not repair the proposed exact fixed-mesh picture. The proof is elementary:
the concurrent Collatz lane constructs `2^k` different ancestors at exactly
`5k` standard steps, whereas either lattice has only quadratically many sites
within that radius.
This concerns a discrete fixed-spacing lattice, not arbitrary curves in the
continuous plane or scale-changing drawings. The independently concurrent
five-step construction strengthens the radius constant. Details are Principle 3.

## 1. A quotient is exact only for its declared consumers

### Theorem Q: factorization and operation descent

Let `q:X->Y` be surjective and `p:X->Z` a target observable. A map
`pbar:Y->Z` with `p=pbar o q` exists **iff** `p` is constant on every fibre
of `q`. Let `F:D->X` be a partial operation. It descends to an exact partial
operation on Y, including its legality, **iff**

```
qx=qy => [x in D iff y in D],
qx=qy and x,y in D => qF(x)=qF(y).
```

**Proof.** Necessity follows by applying the descended function. For
sufficiency choose a representative of each fibre; the two conditions make
both the domain and output independent of that choice. The same argument
works for multiple input slots, with a common tuple of witnesses.

For an operation-closed language C of admissible one-hole continuations,
put `x~y` when every continuation has the same terminal response on x,y,
with an explicit illegal symbol when legality is observable. Composition
shows that `~` is an operation congruence: a future context separating the
outputs would compose with the operation to separate the inputs. This is
classical continuation/Nerode reasoning, already explicitly recovered in
`05-knowledge/results/clock_even_guarded_composition_20261026/THEORY.md`.
It is not a newly claimed abstract theorem.

### A useful linear strengthening: exact boundary rank

Let a piece have a finite boundary B and response vector `w in K^B`, over
`K=R` or `Q`. Suppose gluing to a legal context with response `c` gives
`<w,c>`. A linear compression `L` preserves every such response **iff**

```
ker L subset span{c : c is a legal context}^perp.
```

Indeed every difference erased by L must annihilate every consumer, and
that condition defines the descended linear forms. Hence the minimum
linear dimension for all responses is the dimension of the legal context
span. In particular, if legal pinned contexts supply all delta vectors,
no proper linear compression of the full table can work. **Realizability
matters:** formal delta vectors do not prove that unprescribed geometric
colouring gadgets realize those pins.

For two colour-constrained graph pieces with common named boundary B,
disjoint private interiors, and no cross-edges between those interiors,
the number of extensions after gluing is

```
sum_b w_G(b) w_H(b).
```

Products are taken on the **same boundary assignment** before summing.
Additional cross-edges require more boundary data, not multiplication of
marginal counts. Boolean feasibility uses existential rather than counting
semantics; the two quotients are not interchangeable.

### Typed ledger for Principle 1

| Source -> target; map | Preserved predicate | Lost information / repair | Decisive hostile or test |
|---|---|---|---|
| Collatz word -> diagonal/trace of its affine matrix | counts and multiplicative scale | carry and arithmetic source guard; keep the full affine map | Words `(1,7)` and `(7,1)` have matrices with diagonal `(9,256)` but carries 5 and 131. Sources 483 and 469 merge at 17. |
| Colouring piece -> proper boundary marginals | each recorded marginal | full-support interaction; keep the complete joint relation or the proven legal context quotient | Even/odd subsets of `{0,1}^3` have identical proper counting marginals; pinning `000` separates them. |
| Tournament -> ordinary scalar H or zeta | only the named scalar/spectral observation | continuation response, order, path-cover incidence | All transitive L_n have H=1 and zeta=1; actual context `C3[L_n,1,1]` has H=`1+2^n`. |
| Independent torus -> physical one-parameter curve | restrictions along that curve | relation-lattice character ideal; retain independent endpoint response | THM-2615 computes `ker R_w=<chi_a-1:a.w=0>`; further samples on the same curve cannot recover the kernel. |
| Nodal differentials -> normalization | branchwise differential data | conductor-supported class `x dy` | THM-3955 proves this nonzero class vanishes after normalization; branchwise vanishing alone is insufficient. |
| Ordered prime-wheel gaps -> length and total | those two statistics | collision partition of all partial sums modulo the next prime | Actual Q=30, p=7 words `(2,4,2,4,6)` and `(4,2,4,6,2)` have the same length/total but different next length distributions. |

The exact all-modulus tournament carrier is an instructive **positive**
example: the factorial transform `Psi_T(x)=sum_c (c-1)! p_T(c)(x)_c`
represents the complete declared homogeneous-context quotient modulo M.
This is a proved source/target theorem, unlike an informal assertion that
all spectra classify all future operations. It does not preserve integer H,
order, arbitrary marked ports, real inequalities, or other context languages.

### Session139 distinguishes three different meanings of “all contexts”

The fixed weighted five-vertex pair A,B is stronger than the transitive
order example. It has `B=U A U^T` with rational orthogonal U fixing the
all-ones vector. For every binary **homogeneous one-port** exterior context
Q, its cross-blocks are zero or all-one vectors, so block conjugacy gives

```
Q[B]=diag(U,I) Q[A] diag(U^T,I).
```

Thus all those substitutions retain exactly the adjacency/skew spectra,
determinant zeta, and **all-one** marked resolvent. This does not claim
equality of every individually vertex-marked response. Nevertheless

```
H(Q[B])-H(Q[A])=H(Q)/3456.
```

The lost first path-cover coordinate survives every such exterior context.
At order five M23 (disjoint digon/complement-triangle incidence) restores
H. The same pair also has identical zero mean four-face cocycle but
different second moment, whose sidecar is M22. These are different losses.

This is compatible with the positive all-modulus **H-response** quotient:
the latter is defined by the Hamiltonian observable, not by spectra. An
all-context theorem for one observer is not an all-observer theorem. The
actual binary simulators at primes 3,5,7,11 preserve modular H responses,
not exact real H, order, zeta, or the rational pair's contextual cospectrality.
Source: `session139.../all_moduli/WEIGHTED_BINARY_SIMULATION.md`, section 3,
plus the independently audited weighted zeta and cocycle proofs. These
reconstruction hostiles do not exclude spectral sufficient inequalities.

## 2. Path pictures need an exact clock and a target-reflecting map

### Theorem S: exact first-return lifting, including recurrence

Let `F:X->X`, `R:V->V`, and let `J:X->V` be injective with image B. Suppose
for every x the **first positive return time**

```
tau(x)=min{t>=1:R^t J(x) in B}
```

is finite and `R^tau(x) J(x)=J(F(x))`. Set
`a_k=sum_(i<k) tau(F^i x)`. Then:

1. `R^a_k J(x)=J(F^k x)`, and these are exactly the visits of the R-orbit
   to B. Thus `F^k x=y` for some k iff `R^t J(x)=J(y)` for some t.
2. For every target `A subset X`, first hitting of J(A) from J(x) occurs
   at `a_h`, where h is the first F-hit of A; if one hit does not exist,
   neither does the other.
3. x is periodic under F iff J(x) is periodic under R. If its least
   F-period is p, its least R-period is
   `sum_(i=0)^(p-1) tau(F^i x)`.
4. x is eventually periodic under F iff J(x) is eventually periodic under
   R. If every v in V eventually enters B, no R-cycle or forward orbit can remain forever
   outside the designated states.

**Proof.** Induct at the successive first returns. Their times strictly
increase and tend to infinity because each tau is a positive integer.
There are no visits to B between successive a_k by definition. This proves
1 and 2. An R-return to J(x) must therefore be a designated return, and
injectivity of J makes it an F-return to x; the first such return proves 3.
An eventual R-cycle encountered by J(x) contains a later designated state,
since designated returns continue forever. That state is periodic for F
by 3, proving the nontrivial direction of 4. The reverse direction follows
by expanding the finite F-cycle. The final assertion follows from the
extra entry hypothesis. QED.

**Canonical construction.** For any F and any positive integer tau, use

```
V={(x,j):0<=j<tau(x)},
R(x,j)=(x,j+1) if j+1<tau(x),
R(x,tau(x)-1)=(F(x),0),       J(x)=(x,0).
```

The hypotheses hold, giving an exact variable-time suspension. A geometric
subdivision with marked original vertices is an instance when route
interiors contain no other marked vertex and its labelled outgoing ports
agree. Clock values may be unbounded across X; every particular finite
cycle still has a finite summed period. Reachability and recurrence are
preserved, while native periods, step lengths and local metric balls are
not. This is the precise positive repair of an overly local lattice model.

**Endpoint-only hostile.** Let F fix two states 0 and 1, let R swap them,
and take J to be identity and proposed route lengths 2. Then
`R^2 J(x)=J(F(x))` for both x, but R goes from J(0) to J(1) whereas F does
not go from 0 to 1. The endpoint equations alone are insufficient: a route
passes through another designated state before its claimed endpoint.
First-return/port information is necessary. Separately, if routes share a
site but require different next moves, that site cannot be one functional
state without an extra lane label. Such labels are legitimate but change
what “one state per lattice node” means.

### Positive exact-clock model: spanning sites by nonlocal jumps

Here C denotes the same standard Collatz map T defined above.

Choose an explicit bijection `ell:N+->Lambda` to the triangular lattice of
hexagon centers: label the origin first, then the six axial sides of each
radius-r ring in a fixed orientation. The ring has 6r sites and the ball
has `1+3r(r+1)` sites. Define

```
F_Lambda=ell o C o ell^(-1).
```

This is an **exact conjugacy**, not a quotient or time change. It occupies
every lattice site and preserves one-step arithmetic dynamics, indegrees,
components, hitting and periods exactly. Its arrows are nonlocal jumps,
not necessarily lattice edges and not a single spanning orbit. The chosen
ring order is a convention, not an intrinsic unique geometry.

The radius of the site ell(n) is

```
r(n)=ceil((sqrt(12n-3)-3)/6)=sqrt(n/3)+O(1).
```

For even n, the lattice distance of the jump is at least
`r(n)-r(n/2)=(1-1/sqrt(2))*sqrt(n/3)+O(1)`, so this explicit model pays
unbounded jump lengths. It provides a precise positive alternative:

- **Site conjugacy:** all sites, original clock, nonlocal jumps.
- **Honeycomb subdivision:** local lattice edges, marked original vertices,
  auxiliary sites, variable clock, no claim of a spanning site bijection.
- **Free induced Euclidean realization:** unit edges and exact simple
  adjacency, but no uniform separation or canonical coordinates.

These alternatives preserve different predicates; none is rejected merely
because it is not the others. The explicit site model is credited to the
root coordinator, and the honeycomb/Euclidean existence constructions to
the concurrent Collatz/HN lanes respectively.

### Theorem H: when a quotient genuinely reflects hitting

Suppose `qF=Gq`. If `q^{-1}(B)=A`, then for each x and each k,

```
F^k x in A iff G^k qx in B.
```

Thus exact first hitting time transfers. If `q^{-1}(B)` strictly contains
A, a sufficient repair is to prove separately that **every** point in
`q^{-1}(B)` eventually reaches A. Semiconjugacy plus quotient convergence
then implies source convergence. Without that repair the implication is
false: collapse the translation `F(n)=n+1` on Z to one fixed point. The
quotient is at its target forever; the source orbit from 1 never hits 0.

For forward-invariant A, universal finite hitting is equivalent to the
existence of a function `r:X->Z_{>=0}` with `r=0` on A and
`r(Fx)<r(x)` for x outside A. One direction is well-ordering; the other
uses the actual first hitting time. This equivalence is not an effective
construction of r. THM-4507 already excludes a broad class of global
one-step logarithmic potentials with finitely many valuation corrections;
finite features and a planar picture do not supply the missing rank.

### Typed ledger for Principle 2

- **Source:** the actual arithmetic functional graph, including legal branch.
- **Target:** a routed graph, macro-transition system, or residue machine.
- **Map:** `J` with a specified time cocycle, or semiconjugacy q.
- **Preserved:** hitting of the corresponding marked target under the stated
  hypotheses, and exact step count only when the clock is retained.
- **Lost:** intermediate states, parity/valuation guards, time, and possibly
  root-fibre membership.
- **Sidecar:** arithmetic label, branch guard, elapsed time, landmark marker,
  and a proof that target fibres contain no hidden nonterminating dynamics.
- **Hostile:** constant quotient of integer translation; shortcut Collatz
  and standard Collatz have different native clocks and cycle lengths.
- **Cheapest test:** commute each actual edge with the representation,
  check legality on every proposed quotient fibre, and compare first hits.

This recovers the older exact precision arrow
`Z/(qm)->Z/m, [x]->[floor((ax+b)/q)]`. Integer division generally does not
factor through `Z/m`. The shortcut Collatz first k parity observations need
`2^k` residue states; a fixed finite clock is a bounded-horizon observer,
not an all-time arithmetic state.

## 3. Geometry has independent growth and colouring capacity constraints

### Theorem G: growth obstructs bounded edge stretch

Let G,H be graphs, and `f:V(G)->V(H)` injective. If every edge of G has
image endpoints at H-distance at most L, then

```
|B_G(v,r)| <= |B_H(f(v),Lr)|.
```

**Proof.** Map each path of length at most r; the triangle inequality puts
its endpoint within Lr, and injection prevents cardinality collapse.
Consequently a source with exponential ball growth cannot map this way
into a target whose balls have a uniform polynomial upper bound.
No planarity assumption or edge-disjoint routing is used.

### Apply the stronger concurrent Collatz certificate, without reproving it

The concurrent Collatz lane (`../collatz_hex/PROOF.md`, Theorems 1--3)
constructs, for every root n prime to three, `2^k` distinct predecessors
at **exactly 5k standard steps**. Hence any injective routing into a lattice
with ball upper bound `(2R+1)^2` must have maximum edge stretch through
that finite witness satisfying

```
D_k >= (2^(k/2)-1)/(10k).
```

No uniform bound is possible, in any component. The lane's proof also
handles bounded fibre multiplicity. We do not duplicate that sharper proof
here. `verify.py` retains a separate, slower increasing odd-ancestor
construction from root 5 as an independent arithmetic control; its constant
is eight rather than five and is not the session's main certificate.

**Continuum boundary.** Suppose an injective Euclidean realization has
edge lengths at most L and all distinct image points at least epsilon
apart. Disjoint radius-`epsilon/2` disks around the `2^k` predecessors fit
inside the radius-`5Lk+epsilon/2` disk around the root image. Therefore

```
2^k <= (10Lk/epsilon+1)^2,
```

impossible for fixed L and epsilon. In particular any exact unit-edge
plane realization of Collatz must have arbitrarily close **nonadjacent**
pairs. For the depth-k finite family its minimum separation is at most
`10Lk/(2^(k/2)-1)`. This is compatible with the HN lane's constructive
induced-unit embedding of every countable locally finite simple pseudoforest:
generic unit-circle choices avoid exact unwanted distances but give no
uniform separation. Geometric existence and discrete fixed-capacity
realizability are different predicates. Neither tests convergence.

Exponential inverse-tree growth is classical; see Applegate--Lagarias,
*Density Bounds for the 3x+1 Problem, I. Tree-Search Method*,
[author manuscript](https://dept.math.lsa.umich.edu/~lagarias/doc/applegateI.pdf).
That source uses the shortcut map. The session's standard-time certificate
and geometric corollaries make no priority claim for inverse-tree growth.

| Representation | What is required | What it does not automatically preserve |
|---|---|---|
| Planar topological drawing | vertices/edges with no unwanted intersections | unit distance, fixed mesh, native dynamical time |
| Unit-distance drawing | graph edges have Euclidean length one | crossing-freeness, absence of extra unit edges, separation |
| Induced unit-distance drawing | adjacency iff Euclidean distance one | crossing-freeness or a positive minimum vertex separation |
| Subdivision in a fixed lattice | original edges become finite lattice paths with marked endpoints | one-step time, original unit metric, chromatic lower bounds |
| Bounded-stretch injection into a fixed lattice | uniformly bounded route lengths and distinct lattice sites | impossible for the actual whole Collatz component by growth |

### Independent chromatic and cycle obstructions

1. Standard Collatz has the directed triangle `1->4->2->1`. Honeycomb is
   bipartite, so even a one-edge injective realization of this tiny part
   already fails. The growth argument is stronger because it rules out
   every fixed routing length, not merely length one.
2. Every graph homomorphism `G->H` pulls a k-colouring of H back to G.
   Hence direct unit-neighbour realization in the triangular lattice forces
   `chi(G)<=3`; in honeycomb it forces `chi(G)<=2`. This cannot carry a
   known 4-chromatic unit-distance obstruction, let alone de Grey's
   5-chromatic examples. The latter are primary-source established in
   [de Grey's paper](https://arxiv.org/abs/1804.02385).
3. Subdivision is not a chromatic-preserving repair. Once subdividing
   every edge makes **every** loopless graph bipartite, with original
   vertices on one side and subdivision vertices on the other. K4 versus
   its subdivision is a tiny combinatorial hostile, not a planar unit
   realization of K4. The same loss applies to subdividing a genuine HN
   witness: the routed graph is no longer the original unit-distance graph.
4. HN's plane colouring is not just a colouring of lattice sites or tile
   contacts. It constrains every pair of plane points at Euclidean distance
   one. A hexagonal tiling used as a colouring partition has a different
   source and target from the honeycomb adjacency graph.

### An explicit three-colouring of the actual standard Collatz graph

The root coordinator observed the simple exact colouring

```
c(n)=floor(log_2 n) mod3.
```

This is independently checked here. An even step lowers the integer
bitlength exponent by one. If `2^k<=n<2^(k+1)`, an odd output obeys
`2^(k+1)<=3n+1<2^(k+3)`, so the exponent rises by one or two. Every
standard edge therefore changes colour. The triangle 1,4,2 uses all three
colours, hence the underlying simple Collatz graph has chromatic number
**exactly three**, unconditionally and constructively. This does not produce
a lattice embedding: colouring capacity and spatial growth capacity are
independent necessary tests. Nor does it colour all unit-distance pairs
in R^2; only the actual Collatz edges are the source relation.

### The degree condition alone says neither convergence nor uniqueness

`outdegree=1` and `indegree<=2` is satisfied by integer translation
`n->n+1` on Z, which has no periodic orbit, and by disjoint finite directed
cycles, which have arbitrarily many basins. These have different dynamics
but the same local degrees (indeed all indegrees equal one). On the
honeycomb itself, choosing a consistent direction along its disjoint
infinite zigzag rows gives a one-in/one-out map with no finite cycle.
Thus even the mesh does not force convergence. At a branching node the
whole functional graph is not a single simple directed path.

The older THM-4523 backward-tree separation statement applies to vertices
**prime to three**. Every positive multiple of three has only a doubling
backward ray, so the pure backward trees of 3,6,9,... are isomorphic.
Whole-graph rigidity uses forward enrichment: fix the 3-units, then for an
odd multiple m of three the fixed vertex `3m+1` has fixed even predecessor
`2(3m+1)` and the remaining predecessor m, which is therefore fixed; the
doubling chain above m is then fixed. This credits the inherited separation
lemma, not a reproof of it. Rigidity explicitly does not prove connectivity.
Unique arithmetic labelling, a unique geometric drawing, a single basin,
and convergence are four different predicates.

### Constructive lifting is stronger than merely matching local edges

The HN lane proves a complementary exact lift criterion. Given an abelian
group quotient `q:A->B`, distinct vertex labels y_v in B, and selected
allowed step lifts s_uv satisfying `q(s_uv)=y_v-y_u`, there are positions
x_v in A with `x_v-x_u=s_uv` iff every cycle sum of the lifted s is zero
**in A**, not merely modulo ker q. A spanning tree constructs the positions;
each remaining edge is certified by its fundamental cycle sum. Distinct
quotient labels make the resulting positions injective. Thus this is a
constructive finite certificate, not only an information-loss warning.

The explicit hostile uses `A=Z[omega]`, `omega^2+omega+1=0`, unit steps
`+-1,+-omega,+-(1+omega)`, and quotient `a+b omega -> a+2b mod7`. The
source triangular lattice is 3-colourable; the quotient adjacency is K7.
The quotient triangle `0->1->2->0` uniquely lifts to steps `1,1,-omega`,
whose sum `2-omega` is nonzero. Every local edge is correct, but the
closed walk does not lift. Consequently the quotient's chromatic lower
bound cannot be transferred back. The exact cycle displacement is the
missing sidecar. This example and proof are credited to the concurrent
`../hadwiger_nelson/PROOF.md`, not independently rediscovered here.

### Typed ledger for Principle 3

- **Source:** arithmetic graph or unit-distance colouring graph.
- **Target/map:** injective bounded-edge-stretch drawing; for colouring,
  a graph homomorphism into the proposed unit-neighbour substrate.
- **Preserved:** enough adjacency to transport paths/ball size; homomorphism
  pulls target colourings backward.
- **Lost if relaxed:** spatial scale, collision identity, native time,
  original unit distances, or original chromatic obstruction.
- **Sidecar:** explicit routing lengths and lane/landmark labels for dynamics;
  exact coordinates, squared-distance certificates and full port-colouring
  relation for HN. These are different sidecars, not one interchangeable one.
- **Cheapest tests:** trivial Collatz cycle; binary ancestor count versus
  lattice balls; inherited 3-colouring; subdivision hostile.

## 4. Quantized certificates must survive the actual allowed corrections

### Exact repair quotient

Let A:U->V be a homomorphism of abelian groups and let K be a desired set
of corrected defects in V. There exists u with `v-Au in K` **iff**

```
[v] in q(K),          q:V->V/im A.
```

This is the definition of a coset, but its coefficient ring is load-bearing.
For a finite graph Laplacian over R, the quotient is represented by one
mean per connected component. Over Z the same matrix may leave torsion.
For the triangle Laplacian

```
L=[[2,-1,-1],[-1,2,-1],[-1,-1,2]],       v=(1,-1,0),
L(v/3)=v,                              L(v)=3v.
```

An integral solution of Lz=v would imply `3(z_0-z_1)=2`, impossible.
Thus zero mean proves real repairability but not integral repairability.
This is the concrete older THM-3990 real/integral firewall. At a node,
normalization similarly discards a conductor class rather than repairing
it; THM-3955 supplies the explicit local kernel.

### A genuinely useful positive certificate: the full arithmetic jet

The new remote all-modulus tensor statement was read directly, including
its scope warning. Let m>=2, a nonzero modulo m, and
`delta=min_(t in Z)|a+mt|`. If an integer multiaffine P has top coefficient
in `a+mZ` and every proper coefficient in mZ, then for `|x_i|<=1`,

```
sum_(S subset [r]) (partial_S P(x))^2
    >= delta^2 product_i (1+x_i^2).
```

This is `05-knowledge/results/tensor_primorial_ports_20261102/TENSOR_COSETS.md`,
an elementary all-modulus extension of the inherited parity tensor lemma,
with same-agent proof audit in its source. It is not a generic consequence
of a nonzero scalar residue. In the coefficient-reversed induction, write
the final split as `(u,v+xu)`. The lower-dimensional coefficient vectors
for u and u±v lie in the same permitted coset. If their norm floor is L,
write `A=||u||^2,C=||v||^2,D=<u,v>`; then

```
A>=L,       A+C+2D>=L,       A+C-2D>=L,
A+||v+xu||^2
 >=(1-|x|+x^2)A+(1-|x|)C+|x|L
 >=(1+x^2)L.
```

All three nonnegative terms are needed. The base is the integer-coset
norm floor delta. This independently checks the short inductive mechanism.

A scalar value may vanish: `3x-2=0` at x=2/3, despite valid parity
coefficients. Drop divisibility of proper coefficients and `P=x-1` at
x=1 has full jet energy 1<2. Drop `|x|<=1` and `P=x+2` at x=-2 likewise
has energy 1<5. Thus the retained jet, residue class and bounded parameter
cube each perform a proof step.

For tournament Pfaffians the available hypothesis is parity. Proving the
norm theorem for moduli 3,5,11 or composites does **not** prove that actual
Pfaffian coefficients occupy those stronger cosets. The source correctly
keeps this application boundary. Likewise the incoming order-42 capacity
proof is marked same-agent audited; we read its full proof and audit, but
do not upgrade it to independent acceptance in this synthesis.

### Typed ledger for Principle 4

- **Source:** a defect plus its allowed correction image, or a coefficient
  vector plus its proved integral residue class.
- **Target/map:** the actual repair quotient; alternatively the full
  derivative-jet map at x and its squared norm.
- **Preserved:** exactly the coset obstruction; under the tensor hypotheses,
  a quantitative nonzero norm floor after real evaluation.
- **Lost:** component identity, torsion, transverse derivative energy, or
  divisibility when passing too early to means, real scalars, or values.
- **Sidecar:** coefficient ring and integral image; full jet and actual
  coefficient congruences. A spectral trace is not a substitute.
- **Hostiles:** triangle integral repair; `3x-2`; `x-1`; `x+2` outside the cube.
- **Cheapest test:** compute the actual image/kernel, test its smallest
  nontrivial coset, and verify every application hypothesis before borrowing
  the norm inequality.

## Recovered history: what is genuinely connected, and what is not

The following **12 representative mechanisms across seven subject families**
were inspected, not counted from titles or asserted as an exhaustive repo
census. Their statements and precise source paths are in SOURCE_PINS.json.
“Progress” here means the specific proved mechanism or obstruction, not a
claim of external originality.

| # | Mechanism recovered | Role in this session |
|---:|---|---|
|1| Collatz backward-tree rigidity (THM-4523) | Distinguishes global arithmetic shape from local degree and from connectivity. |
|2| Finite-valuation potential obstruction (THM-4507) | Blocks a tempting finite-feature monotonicity proof; preserves quantifier order. |
|3| Collatz affine carry/trace dictionary (THM-4604) and guarded precision clocks | Identifies what arithmetic and exact-time representations must retain. |
|4| Seven-image unit cube classification / hexagonal flower | Exact metric quotient with a specified collision fibre, not a chromatic analogy. |
|5| Full colour-boundary table and its lost Walsh modes | Gives exact gluing and explicit pair-marginal hostiles. |
|6| All-modulus homogeneous tournament contexts | Positive example of a proved operation-relative sufficient quotient. |
|7| Weighted tournament zeta/mean blindness | Negative examples: spectral and first-moment information miss different incidence coordinates. |
|8| Tensor matching reserve and incoming capacity potential | A successful quantitative certificate that retains compensation, not just an obstruction. |
|9| Prime-wheel collision partitions and marked-port recurrences | Counts compose only with the right internal boundary data; counts still lose location. |
|10| LRC physical-diagonal toric kernel (THM-2615) | More measurements inside one restricted carrier never recover its kernel. |
|11| Node differential conductor class (THM-3955) | A local normalization quotient can erase a real obstruction. |
|12| Componentwise Poisson / integral repair quotient (THM-3990) | Defines the actual correction-invariant object and its ring-dependent boundary. |

The seven families are Collatz dynamics, metric colouring, tournament
invariants/stability, sieve arithmetic, LRC harmonic/toric observations,
singularity algebra, and Poisson/integral gluing. This is a selected verified
connection inventory, **not** a count of all novel results in repository
history. The common theorem-level discipline is to specify the consumer,
operation and loss kernel before choosing a representation. It does not
identify all seven mathematical objects with one another.

## Proof targets that remain after the synthesis

1. **Collatz:** a valid arithmetic/ranking argument excluding other cycles
   and nonperiodic positive components. A fixed local hex mesh cannot do this
   by exact bounded-stretch injection; a useful replacement must state the
   unbounded scale/extra-state cost and preserve source guards and time.
2. **HN:** realizable unit-distance gadgets with complete legal boundary
   relations that compose to an unprescribed colouring obstruction. Pinned
   palettes, lattice site colouring, graph minors and subdivisions do not
   by themselves supply it.
3. **Shared method:** seek the smallest *realizable* family of contexts
   separating the target-bearing kernel, or prove no proposed restricted
   family can separate it. This is sharper than enriching an invariant
   without naming which lost distinction the enrichment restores.

Reproduction:

```
python3 verify.py > normal.json
python3 -O verify.py > optimized.json
cmp normal.json optimized.json
```

The verifier independently checks modulo-9 inverse guards, 16 complete
macro-levels (131071 distinct odd labels including the root), native-time
routes, lattice balls, proper colouring, boundary parity loss, subdivision
loss, explicit carry values, integer repair torsion, and 42525 exact
rational two-variable jet evaluations. It imports no repository source.
It additionally exhausts 4330 finite map/clock systems, checks their
first-return itineraries and exact cycle-time sums, checks the endpoint-only
route hostile, enumerates 30301 ring-ordered sites, verifies 10000 exact
site-conjugate steps and 100000 bitlength-colouring edges. These are controls
for the proofs, not a replacement for their all-depth or all-parameter
reasoning.
