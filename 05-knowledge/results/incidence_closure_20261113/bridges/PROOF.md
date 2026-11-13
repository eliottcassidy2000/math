# Connected fibres, loop transport, and why sheafification takes at most two passes

**Status: PROVED IN NOTES after independent proof/code audit; see [AUDIT.md](../AUDIT.md).** Claims B1–B6
below have audited elementary proofs in their stated scopes.
Inherited results retain their own status and scope. The small controls are
**FINITE-EXACT**. No external priority claim is made. No theorem ID is reserved.

## 1. Bounded inheritance and the intended Neményi theorem

The closest proved mechanism is the tight-flat/block contraction in
[polynomial staircase diameter, Theorem 3 and Corollary 4](../../graceful_polynomial_diameter_20261019/POLYNOMIAL_STAIRCASE_DIAMETER.md).
For labels `0,...,m`, retain an acyclic choice of one edge of every difference
color `r+1,...,m`. The completions to a tree by one edge of each free color
`1,...,r` form a connected graph under changes of at most two free coordinates.
The retained edges never move. For `r>=3` the inherited diameter bound is
`r²−3r+3`; the source gives the small cases separately. This does not preserve
the abstract isomorphism type of the completed tree.

The canonical hostile is a frozen one-color reconfiguration fibre. The corrected
near miss is the requirement to retain every initially target-correct coordinate:
the known five-state example requires moving a correct buffer. These are separate
from the loop obstruction developed here. The least-used sidecar for connectivity
is actual edge transport, including its composition around loops.

For sheafification the closest proved source is
[THM-4892, tournament sites and the plus construction](../../../../01-canon/theorems/THM-4892-tournament-sites-plus-construction.md),
with [its independent audit](../../sheafification_plus_construction_20261107/audit/AUDIT.md).
The repaired near misses are the reversed finer/coarser topology relation and an
unsupported Hamiltonian-path/global-section identification (MISTAKE-731 in
[the correction ledger](../../../../01-canon/MISTAKES.md)). Hamiltonian paths do
not restrict to Hamiltonian paths on arbitrary vertex subsets, so the putative
restriction rule does not define that presheaf. The real OCF packing presheaf and
the pair-orientation presheaf do have restriction maps. The empty-set value is
the decisive sidecar for the strong-component site's second plus pass.

The accented and unaccented repository searches identify **Neményi's isothermal
principal-stress-net result**, not a theorem about graceful reconfiguration.
The exact working version is the local, no-body-force, isotropic plane-elasticity
classification in [STRESS.md](../../five_stress_connection_atlas_20261112/STRESS.md).
The scope and source limitations remain those of
[that packet's audit](../../five_stress_connection_atlas_20261112/AUDIT.md): its
local calculations are explicit; the original 1933 full text and the Sáenz erratum
were not available there. This session likewise does not import an unread erratum.

The primary [Sáenz publisher abstract](https://link.springer.com/article/10.1007/s00707-008-0118-3)
identifies a broader two-to-five-parameter result under simply connected plane
elasticity hypotheses and explicitly links an erratum. Neményi's 1933 article is
[Stromlinien und Hauptspannungstrajektorien](https://doi.org/10.1002/zamm.19330130506).
The 1933 DOI and [2010 erratum DOI](https://doi.org/10.1007/s00707-010-0370-1)
did not yield readable full text in this session. Attribution is not a substitute
for the narrower repository proof.

The live concept board is: fixed boundary; scalar observer; invertible transport;
loop holonomy; matching-family descent; marked stabilization. The anchor is the
incidence-closure problem; the niche is completion-groupoid transport; the wildcard
is the actual second derivative in the five-dimensional stress module.

## 2. B1: connectivity gives scalar descent; transport needs holonomy

Let `q:X→Y` be a surjection and let each fibre carry a connected undirected graph
of allowed moves. For any set `Z`, a map `f:X→Z` factors uniquely through `q` iff
its values agree across every allowed edge. Necessity is immediate. For
sufficiency, follow an edge path within a fibre; the common value defines the
factor map. Surjectivity gives uniqueness. This forgets every within-fibre state
coordinate, so it applies only to move-invariant observers.

Now let `G` be one connected graph, give each vertex `v` a set `E_v`, and each
oriented edge `v→w` a bijection `t_vw:E_v→E_w`, with `t_wv=t_vw^{-1}`.
A section is a family `s_v` satisfying `s_w=t_vw(s_v)` on every edge.
Equivalently, the transports extend to a functor from the graph's path groupoid,
where immediate edge/reverse-edge backtracks cancel. No further loop relations
are imposed just by knowing the move graph.

**B1 (proved).** Evaluation at a root `o` identifies the section set with
the elements of `E_o` fixed by every closed-walk transport based at `o`.

**Proof.** Every section is determined along paths from its root value. Returning
along a closed walk must preserve that value. Conversely, choose a root element
fixed by all closed walks and transport it to each vertex. Two paths to the same
vertex differ by a closed walk after reversing one path, so give the same value.
This defines a section and is inverse to evaluation. Empty fibres are allowed;
the resulting sets can be empty. QED.

Consequently all root states extend iff every loop transport is the identity.
For vector spaces and linear isomorphisms, the same proof gives the invariant
subspace of the loop representation. Mere connectedness gives neither trivial
loop transport nor a nonempty torsor section. If moves have noninvertible maps,
the path-groupoid theorem is inapplicable; use the full diagram limit instead.

## 3. B2: an actual graceful fibre with eight possible obstruction classes

Retain the edge `03` of color 3 on labels `0,1,2,3`, and release colors 1 and 2.
Writing only the free edges, the four completions are

| State | Color 1 | Color 2 |
|---|---|---|
| A | 01 | 02 |
| B | 12 | 02 |
| C | 12 | 13 |
| D | 23 | 13 |

All four become trees after adjoining `03`. The two other choices form a triangle
and leave one vertex isolated. Any two completions differ in at most two free
colors, hence the width-two exchange graph is exactly `K4`.

Give every exchange edge a freely specified bit `s_e`, and let that edge transport
a two-element torsor by `b→b+s_e` over `F2`. No canonical choice of these bits is
inferred from the graceful trees. This is an explicit extra-data experiment.

**B2 (proved).** For a finite graph with `v` vertices, `e` edges and `c`
connected components, changing a bit at each vertex acts on edge bits by

    s_vw → s_vw + x_v + x_w.

There are `2^(e−v+c)` gauge classes. A compatible torsor section exists iff every
cycle has bit sum zero; then there are exactly `2^c` sections. Every other class
has none.

**Proof.** On each component choose a root and a spanning tree. The root bit is
free; the tree equations determine every other vertex bit. Each remaining edge
imposes exactly the zero-sum equation around its fundamental cycle. Equivalently,
normalize all tree-edge bits to zero by vertex gauges. The remaining `e−v+c`
bits are free and classify gauge orbits. Gauges fixing all edge bits are exactly
the componentwise constant ones, which also count solutions when one exists.
QED. This is the elementary graph coboundary calculation, not a new classification
of graceful trees.

Here `K4` has cycle rank 3. Among its 64 edge-bit assignments, eight are balanced
and have two sections; the other 56 have none. The eight gauge classes each have
eight assignments. For a concrete obstruction, put bit 1 on `AB` and zero on all
other edges: the loop `A→B→C→A` reverses the torsor label.

This is also the exact mechanism in the independently owned
[signed lifting candidate](../pairings/SIGNED_LIFT.md): after a tested base
permutation preserves zero support, magnitudes and diagonal entries, ratios of
the nonzero signed entries define edge signs. A lift exists iff these ratio signs
are a coboundary. Raw entries of `K_ij=C_(i,j)−C_(i,Rj)` are not bijective edge
transports and may vanish; only the guarded ratios enter this comparison. The
completion graph and the signed support graph remain different objects, with
different retained boundaries and predicates.

## 4. B3: a sheaf can retain the obstruction completely

Build a two-sheet graph cover from the preceding bits: vertices are `(v,b)`;
each base edge `vw` has lifts from `(v,b)` to `(w,b+s_vw)`. Regard graphs as
topological realizations of their one-dimensional cell complexes. Projection
`p:G_s→G` is a local homeomorphism: an open star at each lifted vertex maps
homeomorphically onto the corresponding base star.

For every open `U` of `G`, let `F_s(U)` be the continuous sections of `p` over `U`.
Restriction is ordinary restriction of functions.

**B3 (proved).** `F_s` is a sheaf. Its global sections correspond exactly to
the torsor sections in B2. In particular a nontrivial holonomy class on `K4`
defines an already-sheafified object with no global sections.

**Proof.** Sections agreeing on overlaps glue to one function on the union.
Continuity is local on an open cover, and the equation `p∘s=id` holds locally,
so the glued function is a section, uniquely. On each edge, a section must choose
one of the two lifted edges, and its endpoint choices must satisfy B2. Conversely
compatible endpoint bits define these lifted edges and give a continuous section.
QED.

Thus sheafification guarantees effective gluing **of matching families**. It
does not assert the existence of such a family, remove every loop obstruction,
or make every root state extend. In this example `F_s^+=F_s` already. The positive
unification is the same compatibility limit; a hypothesis supplying a compatible
family cannot be omitted.

## 5. The precise two-pass theorem and its smallest overlap repair

**CITED.** For a set-valued presheaf on a genuine Grothendieck site, plus is the
colimit of matching families over covering sieves, with families identified after
common refinement. The plus presheaf is separated; applying plus to a separated
presheaf gives a sheaf; hence `F++` is a sheaf and subsequent plus maps are
isomorphisms. The site axioms and functorial restriction maps are essential:
pullback stability supplies restriction, and transitivity combines local covers.
See [Stacks Project, Theorem 7.10.10, tag 00WB](https://stacks.math.columbia.edu/tag/00WB).
This is a general upper bound, not a claim that every presheaf needs two passes.

**B4 (proved elementary model; inherited mechanism).** On the discrete
two-point space `{a,b}`, take the constant presheaf with two-element value `A`
on **every** open, including the empty set, and identity restriction maps.

| Open | `F` | `F+` | `F++` |
|---|---:|---:|---:|
| empty | 2 | 1 | 1 |
| `{a}` | 2 | 2 | 2 |
| `{b}` | 2 | 2 | 2 |
| `{a,b}` | 2 | 2 | 4 |

The first pass makes the empty open have one section. At the full open, matching
two singleton labels originally requires equality of their restrictions to
`F(empty)=A`, so it still permits only the diagonal. In the second pass the
overlap has become a singleton, so all four pairs match. The result is the sheaf
of arbitrary functions to `A`. If `F(empty)` is normalized to a singleton before
starting, the full open obtains four sections in one pass. The same calculation
on a finite discrete space gives `A^U` in two passes, or one after normalization.
This is the explicit finite form of THM-4892(B), not a new sheafification theorem.

The checker computes the colimit over **all** covering sieves on discrete spaces
of sizes 1,2,3, including the empty covering sieve of the empty open. It constructs
all matching families, glues them along every refinement, and checks restriction
maps and the sheaf condition through three iterations. It does not hard-code the
smallest-sieve formula as the algorithm producing plus.

For comparison, finite observer refinement has no universal two-pass repair
bound. Even the deterministic transition `i→max(i−1,0)` on `{0,...,m}`, observed
initially only as `{0}` versus the rest, reveals one new distance layer each
round: block counts `2,3,...,m+1`. There are `m−1` strict refinements. The root
lane's reversible construction is stronger; this small independent control
already prevents identifying generic iterative refinement with sheaf plus.

## 6. Neményi periods are two different representation-level obstructions

**PROVED INHERITED, with the source's scope.** In the no-body-force local model,
write `h=(Z')²`, `q=(sigma_xx−sigma_yy)/2−i tau_xy`, and

    q=h(A|z|²+Bz+conj(B)conj(z)+C),
    F=h(Az+conj(B)),       G=h(Bz+C),
    p=p0−2 Re ∫F dz,

where `A,C,p0` are real and `B` is complex. This is a five-real-parameter family
on a connected simply connected chart with `Z'≠0`.
At an isotropic point `q=0`, the stress does not determine unique principal axes;
retain the prescribed net by continuation, as in the inherited theorem, or work
on a subdomain with `q≠0`.

On a multiply connected domain, the stress is single-valued exactly when
`Re ∮F dz=0`. A globally single-valued elastic displacement needs the stronger
pair of conditions

    ∮F dz=0,
    kappa ∮zF dz + conjugate(∮G dz)=0

for every loop (physical `kappa+1≠0`). The first stronger condition removes the
rotational continuation; the second removes translation. These are not the same
period copied twice, and neither condition comes from the number of plus passes.

For `Z=z+1/z`, `h=1−2/z²+1/z⁴`, the relevant residues are

    Res F=−2A,   Res zF=−2 conj(B),   Res G=−2B.

Thus every stress in the local five-parameter family is single-valued around the
hole: `∮F` is purely imaginary. But at `kappa=3`, displacement requires `A=0`
and `B=0`, leaving only `C,p0`. Setting `A=1,B=0` isolates the rotational failure;
setting `A=0,B=1` passes that test and fails the translational one. Domains also
exclude zeros of `Z'`; those additional holes contribute no residues.

The exact map to B1 is **analytic continuation acting on local potentials or
displacements**, not a claimed bijection from graceful trees to stress fields.
Changing the observer from stress to displacement changes which continuation
increments are invisible. Connectivity of the domain supplies paths, while the
chosen observable determines the required loop-invariance predicate.

## 7. B5: stabilization preserves a loop obstruction; periodicity is separate

**B5 (proved, flat transport retained).** For a finite graph local system of finite-dimensional vector
spaces over a field, and a fixed finite-dimensional vector space `W` with trivial
transport, the natural map

    sections(E) tensor W → sections(E tensor W)

is an isomorphism. Choose a basis of `W`; both sides are the same direct sum of
copies of the edge compatibility equations. Equivalently, loop operators become
`rho(g) tensor I_W`. A loop acting as `−1` over a field of characteristic not 2
still has no nonzero fixed vector after this stabilization. QED.

Here sections mean **parallel** sections for the specified flat transport. The
underlying vector bundle can already be trivial after that transport is forgotten.
For a precise hostile, identify the endpoints of `[0,1]` to make a circle and set

    L=([0,1] × C)/((1,z) ~ (0,−z)).

Flat holonomy is `−1`, so a parallel section is constant on the interval and must
satisfy `z=−z`, hence is zero. But `t↦[t,exp(pi*i*t)]` is a well-defined nowhere-zero
ordinary section: its endpoint values agree under the displayed identification.
It trivializes the complex line bundle. For sixteen real copies, the matrix with
eight diagonal rotation blocks `R(pi*t)` goes from `I_16` to `−I_16` and descends
to a global frame for the underlying real bundle. Thus even that rank-sixteen
bundle is trivial while its specified flat transport has no nonzero parallel
section. These are explicit trivializations, not an appeal to a classification.

This gives a precise representation-theoretic test for a proposed Bott bridge.
The actual classical periodicities are **CITED** in their own categories:
stable unitary homotopy has period 2 and stable orthogonal homotopy period 8
([Bott 1959, pp. 314–315](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/bott4.pdf)).
Real negative-square Clifford algebras satisfy
`Cl_(0,k+8) ≅ Cl_(0,k) tensor M16(R)`; the corresponding graded-module quotient
has the specified degree-eight periodicity
([Atiyah–Bott–Shapiro, §4 and §6](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/abs.pdf)).
No Clifford action has been supplied on these completion fibres. An actual
intertwiner and retained observer would be required to transfer a predicate.
In particular B5 does **not** claim that a nontrivial KO/KU class survives Bott
periodicity or Morita equivalence. Passing from flat transport to an underlying
bundle or its K-class has already forgotten the character in the displayed
examples. Bott's graded K-theory periodicity does not classify global parallel
sections for a prescribed connection.

The inherited [marked-periodicity note, §2](../../collatz_bott_marked_periodicity_20261007d.md)
already separates a lossless displayed tensor frame from its lossy homotopy
class. Keeping a matrix corner recovers the original operator; forgetting to a
stable class can forget arithmetic guards. Neither periodic equivalence nor
matrix amplification proves that the original dynamics returns after 2 or 8
steps. B5 similarly shows that amplification does not repair holonomy.

The stress module also supplies a literal second derivative: on its Cartesian
five-dimensional family, `D_x D_y=0`, `D_x³=D_y³=0`, while `D_x²` need not vanish.
For `sigma_xx=2y²`, `sigma_yy=−2x²`, `tau=0`, the second x derivative is
`(0,−4,0)` and the third is zero. This is a nilpotent differentiation algebra,
not an idempotent reflector. The parent lane's Gram-defect calculation is a
different, explicitly typed second-order diagnostic; detection depth does not
bound closure repair depth.

**CITED.** Quillen's topological plus starts with a connected CW complex and a
perfect normal subgroup of its fundamental group. It gives an acyclic map that
kills precisely that subgroup; acyclicity includes the stated isomorphisms in
homology with local coefficients. See [Quillen, Higher K-theory for categories
with exact sequences, §1, pp. 95–96](https://aareyanmanzoor.github.io/assets/articles/quillen-higher-K-exact-seq.pdf).

**B6 (proved).** A sign holonomy `rho:pi→C2` kills every perfect normal subgroup
`N` of `pi`, and therefore factors through `pi/N`. If `rho` is nontrivial, the
factored character is nontrivial. Indeed `rho([N,N])=0` since `C2` is abelian,
while `N=[N,N]`; the quotient map is surjective. QED. Thus the sign obstruction
in B2 cannot be removed by a Quillen plus map. This uses the actual group map,
not a coincidence between two constructions named plus.

Grothendieck presheaf plus, Quillen's topological plus, and Bott periodicity have
different source categories and maps. A single sign obstruction supplies three
precise tests: an already-sheafified section object may retain it; a perfect-normal
quotient preserves its character; and trivial tensor amplification preserves its
linear fixed-vector obstruction while flat transport is retained. Forgetting flat
transport can destroy that obstruction, as the explicit trivializations show.
None identifies those three constructions with
one another or with a double derivative.

## 8. Connection contract and reproduction

| Source | Target/map | Preserved predicate | Lost without sidecar | Decisive test |
|---|---|---|---|---|
| Fixed-tail completions | Connected exchange graph | Reachability with tail fixed | Tree type; state coordinates | The four completions above |
| Move-invariant scalar | Fibre quotient | Exact scalar value | Noninvariant data | Compare one allowed edge |
| Invertible edge transport | Root loop representation | Compatible sections | Path independence if loops discarded | Fundamental-cycle products |
| Signed adjacency plus tested base map | Ratio-sign transport | Lift existence | Zero support, magnitude, diagonal guards | Coboundary solve on support |
| Local stress potentials | Analytic continuation | Selected stress/displacement observer | Periods invisible to weaker observer | `Z=z+1/z`, two parameter choices |
| Presheaf on actual site | `F++` | Universal sheafification | Original nonseparated distinctions | Empty-overlap two-point example |
| Local system plus tensor frame | Amplified local system | Linear compatibility | Frame after forgetting it | `−1→−I_16` still obstructed |

Reproduce from the worktree root:

```sh
python3 05-knowledge/results/incidence_closure_20261113/bridges/check_bridges.py > 05-knowledge/results/incidence_closure_20261113/bridges/check.normal.json
python3 -O 05-knowledge/results/incidence_closure_20261113/bridges/check_bridges.py > 05-knowledge/results/incidence_closure_20261113/bridges/check.optimized.json
```

The two receipts are byte-identical, with 5,432 explicit checks. The program has
no removable assertions and imports no producer code. It exhausts all 6 raw
color choices in the declared graceful fibre, all 64 transports and 16 vertex
gauges, all covering sieves on the three finite sites, and the declared finite
refinement paths. These finite universes do not establish the general claims
without their proofs. A normal/optimized replay is not independent mathematical
review. The statements above were promoted after that independent review; see AUDIT.md.

Frozen initial check source SHA256:
`32eea5653748102011bded91edcd1c846ebd8aeac7b190253e2cd6e33436759d`.
Both check receipt SHA256 values:
`46bc599c3848ae69e04b71248c2bc96c4d474daa5fca54c391dfa6eb34b06d15`.

**OPEN:** identify a canonical transport or observer on a larger completion fibre
whose holonomy adds a useful invariant to incidence closure. Connectivity alone
does not choose it. The graph-coboundary bridge to guarded signed lifts is exact;
a functor to elasticity or Clifford modules has not been constructed.
