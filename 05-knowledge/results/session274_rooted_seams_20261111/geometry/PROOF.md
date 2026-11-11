# Rooted seams as a birth fan over the inherited holonomy quotient

Status: **PROVED** for the elementary statements below; finite controls are supplementary. This does not prove Chan–Pak, and no seam move is asserted to preserve an individual Hamiltonian-path count.

## 1. Inheritance and conventions

The exact staircase gauge and the complement-line metagraph are inherited from `01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md`, §0. The switching normalization is THM-474, `THM-474-tilings-are-switching-classes.md`. The cut/cycle quotient itself is **not new**: see `THM-1382-map-graph-star-flips-nested-cut-cycle.md` and both THM-1405 star-quotient files. They also explicitly show that this quotient does not descend to tournament isomorphism classes.

For legibility relabel the inherited downward path as

`P = v_1 -> v_2 -> ... -> v_n`.

Set `G_n=K_n\E(P)` and `m_n=(n-1)(n-2)/2`. A tiling `t in F_2^{E(G_n)}` has bit 0 on a forward chord `v_i -> v_j` (`i<j`) and bit 1 on its reverse. Every path edge remains forward. Distinguish:

* a single free-tile flip;
* a seam rectangle, defined below;
* the metagraph line `t <-> t+1`, which flips **all** free tiles;
* replacing the chosen Hamiltonian path in the same tournament, which is a reframing and need not preserve any of the following sector labels.

This last distinction is essential: holonomy below belongs to a framed tournament, not an unframed tournament.

## 2. Seam rectangles recover the inherited star group

For the seam between `v_i` and `v_{i+1}`, switch the cut between the two path intervals and then reverse just the path edge `v_i v_{i+1}`. The result again contains `P`. Its free bits differ from `t` precisely on the rectangle

`R_i = {v_a v_b in E(G_n): a<=i<b}`.

The prefix-cut vectors generate the vertex-cut space `C_n=cut(G_n)`: a singleton cut is the sum of its two neighboring prefix cuts, with the usual endpoint convention. Thus these path-native seam moves realize exactly the already known star group. For `n>=4`, `G_n` is connected. Consequently

`dim C_n=n-1`, `Q_n=F_2^{E(G_n)}/C_n`, `dim Q_n=q_n=(n-1)(n-4)/2`.

An orbit has `2^{n-1}` tilings. Its class is specified by the sums of bits on cycles in `G_n`; these are the inherited F_2 holonomies. Initial sector counts are `1,4,32,512` at orders `4,5,6,7`.

The all-tile vector `1` is a cut iff `G_n` is bipartite. Hence the metagraph complement stays within each seam orbit at `n=4` (where `G_4` is a tree), but crosses sectors at every `n>=5`, since `v_1,v_3,v_5` form a triangle in `G_n`. More precisely, it translates every sector by the same nonzero order-two element `[1]`. Every odd-cycle parity reverses and every even-cycle parity remains unchanged. This is a statement upstairs, before quotienting by isomorphism/converse.

## 3. Every sector has the same mean complete path-cover profile

Let `c_k(T)` count covers of the vertex set by `k` unordered, disjoint, directed paths, with singleton paths allowed. In particular `c_1=H` and `c_n=1`.

**Forest-mean theorem.** For `n>=4`, every coset `t+C_n` has the same mean `c_k`, for every `k=1,...,n`, equal to the uniform mean over the entire tiling cube. In fact the equality holds for every specified oriented linear-forest event, not merely the aggregated coefficients.

**Proof.** Draw a uniform vertex sign vector `x in F_2^n` modulo global complementation. The corresponding orbit point is `t+delta x`, where `(delta x)_{uv}=x_u+x_v` on chord edges. Fix a directed spanning path cover `F`. If `F` uses a path edge backward, it occurs nowhere in the framed cube. Otherwise its only random orientation requirements occur on `r=|E(F) intersect E(G_n)|` chords. Those chords form a forest. The equations `x_u+x_v=b_uv` on a forest have independent rank `r`, for arbitrary right-hand sides. Their probability is exactly `2^{-r}`, independently of `t`. Uniform independent free-tile bits give the identical probability. Sum over covers. QED.

Thus all-order cycle-sector information is invisible to every *first moment* of the complete cover profile. This does **not** identify the distributions, individual profiles, root responses, or higher moments. Products of cover observables can form cycles and therefore retain sector dependence. The construction is a graphic version of a limited-independence design: it fools every forest-supported test, not every test of a given cardinality.

This turns the historical cut/cycle obstruction into a useful resource statement. We may move through genuinely different sectors while keeping all these aggregate baselines identical, then use a rooted exposure to access the information the averaging lost. It does not guarantee any selected H inside any sector.

## 4. Physical rooted birth fills every compatible child sector

Let `a=v_j` be the root of the framed tournament. Add a new vertex `b` satisfying

`b -> a`, and `u -> b` for every old `u != a`.

Insert `b` immediately before `a` in `P`. This is always a Hamiltonian path of the child. Denote the resulting affine bit embedding by `B_a`. The original free-chord bits embed unchanged. If `j>1`, the old path edge `p->a`, `p=v_{j-1}`, becomes a free chord fixed forward. Every new free chord incident with `b` has its physically prescribed orientation. Their orientations are fixed independently of the source tiling.

Restriction of child chord cochains to old `G_n` induces a surjective map

`rho: Q_{n+1} -> Q_n`.

Indeed a cut restricts to a cut, and any old cochain can be extended. The dimension difference is

`q_{n+1}-q_n=n-2`.

**Birth-fan theorem.** Fix ANY parent seam orbit `t+C_n`, with `n>=4`, and ANY root position `j`. Its `2^{n-1}` physically born children meet **every** sector in `rho^{-1}([t])`, exactly twice each. There are `2^{n-2}` such child sectors.

**Proof, internal root.** Compare `B_a(t+delta x)` and `B_a(t)`. Their difference is `delta x` on old chords and zero on all new chords. It is a child cut `delta y` iff, after a global gauge adjustment, `y_u=x_u` on all old vertices; this follows because `G_n` is connected. Vanishing on the newly freed chord `pa` forces `x_p=x_a`. Vanishing on the newborn chords forces `y_b=x_u` for every old `u` outside `{p,a}`. Thus `x` must be constant on that outside set and constant on `{p,a}`. Modulo a global constant, the kernel is exactly one-dimensional, generated by the old cut of `{p,a}`. The image therefore has dimension `(n-1)-1=n-2`. It lies in a fiber of `rho` of that same dimension, so fills it with multiplicity two.

**Proof, first root.** There is no newly freed old path edge. Newborn chords join `b` to every old vertex except `a`; vanishing forces all `x_u`, `u!=a`, equal. The kernel is generated by the old cut of `{a}`, again dimension one. The same conclusion follows. QED.

No small-order extrapolation enters this proof. The statement is also independent of the actual constants chosen on the new free chords: any fixed framed one-vertex insertion has the same linear fan. The rooted-sink birth above is singled out because it is an explicit physical count operation, not just an affine embedding.

**A useful recursive interpretation.** Every old holonomy sector admits every compatible next-level extension after allowing its seam controls. The kernel records precisely the one bit of old gauge choice that becomes invisible in the child quotient. Recursive sector supply is therefore complete, although scalar H supply need not be.

## 5. A lossless birth-and-seam compiler for every framed tournament

The sector fan has a stronger coordinate-level form, valid already at `n>=1`, without a connectivity assumption.

**Paired-root compiler.** Keep the newborn `b` and old root `a` at gauge value zero. Choose arbitrary binary gauge values `y_u` on the other `n-1` old vertices. Then

`(t,y) -> B_a(t)+delta_{G_{n+1}} y`

is a bijection from the parent tiling cube times `F_2^{n-1}` onto the entire child tiling cube. The gauge operation is realized by a sequence of the seam rectangles from §2. It need not preserve H.

**Explicit inverse.** Given any child tiling `z`, compare its new free-chord bits with the fixed bits of `B_a`. Every chord `bu` for `u` outside `{p,a}` determines `y_u`, because `y_b=0`. At an internal root, the additionally freed old path edge `pa` determines `y_p`, because `y_a=0`. At a first root, all `u!=a` are determined by newborn chords. All `n-1` unknown bits are therefore recovered uniquely. Subtract their chord-cut vector from `z`; the new chord coordinates now equal the birth constants. Restriction to old chords recovers a unique `t`. This proves existence and uniqueness for every order. QED.

Using first-root births at each stage gives a unique recursive program from the singleton to **every** framed tournament. The program consumes exactly

`0+1+...+(N-2)=binom(N-1,2)`

binary controls. Choosing those controls independently and fairly samples the tiling cube exactly uniformly. The usual number `2^{binom(N-1,2)}` is inherited; this is an explicit seam-compatible compiler for it, not a new enumeration formula. Each physical rooted birth preserves a chosen Hamiltonian frame, and its subsequent paired-root seam controls insert all missing chord information.

This generator is not uniform on tournament isomorphism classes: the inherited frame-fiber weight of a class is `H(T)/|Aut(T)|`. Nor is it a bounded-state proof of Chan–Pak. Its control word grows with size, and the inverse algorithm targets a fully specified child tournament rather than a requested scalar H. The progress is a complete physical grammar whose **only remaining constructive selection problem** is how to choose seam controls for H while retaining the interfaces needed later.

Provenance scope: THM-474 already supplies unique switching normalization along a spanning path, and THM-1382/1405 already supply the star/cut geometry. This compiler makes the rooted-birth factorization and inverse explicit; it does not claim novelty for those background gauge principles or for the tiling cardinality.

## 6. Birth does not generally commute with the complement-line operation

**Scope in this section: parent order n≥4**, so the old chord graph is connected.
At smaller orders the internal-root nonzero-defect conclusion need not hold.

Write `J_n(t)=t+1`. Then

`D_a = B_a(J_n t)+J_{n+1}(B_a t)`

is independent of `t`: it equals 1 on precisely the newly introduced free coordinates and 0 on old chords. This is the exact defect between the two routes.

If `a` is first in `P`, those new coordinates form the star of `b`, so `D_a in C_{n+1}`. Hence the two routes agree at child-sector level.

If `a` has a predecessor `p`, the new coordinates consist of newborn chords **and the newly freed chord `pa`**. Were `D_a` a cut, its zero restriction to connected old `G_n` would force all old gauge values equal. Such a cut cannot take value 1 on `pa`. Contradiction. Therefore `[D_a] !=0`.

An explicit detector is any path from `p` to `a` in old `G_n`, closed by the new chord `pa`. This cycle has zero defect on its old portion and defect one on its closing chord. It witnesses precisely the old path-edge information released by an internal seam birth.

The seam birth is therefore compatible with metagraph complementation at the initial boundary but has a nontrivial holonomy defect at every internal location. It is incorrect to infer recursive metagraph line transport from an unmarked one-vertex construction alone.

## 7. Exact initial controls and remaining obligation

`python3 geometry/check.py` (from this packet's parent directory) exhausts the tiling cubes at `n=4,5,6`, partitions every seam orbit, computes every path-cover coefficient by independent subset recurrences, and checks every root position in every orbit. It confirms common profile sums, the exact child-sector fan, complement crossing, and the boundary/internal complement defect. It additionally checks the paired-root compiler bijectively for every root at n=1,...,6, reaching all 32,768 child tilings at order 7. The small enumeration is a structural control, not new evidence for ultimate Chan–Pak coverage.

The missing bridge is an **H-selective sector-and-gauge controller**: the birth fan supplies every new cycle-charge configuration, but choosing a sector still leaves many physical tournaments and no guaranteed prescribed H. A useful candidate theorem must retain the rooted response (or equivalent full interface) together with the sector and prove that recursively selecting gauge representatives hits a requested odd H. Neither forest means, complete sector supply, nor metagraph connectivity by itself supplies that selector.
