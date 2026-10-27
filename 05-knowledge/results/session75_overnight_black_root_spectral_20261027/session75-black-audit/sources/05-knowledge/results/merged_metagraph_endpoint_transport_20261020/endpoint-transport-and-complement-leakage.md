# Endpoint transport and complement leakage in the fixed-path metagraph

Status: PROVED in the precise scopes below; finite controls are separately labelled.

## Scope, inheritance, and useful terminology

The object is THM-796, `THM-796-three-sorted-recursive-tiling-line-node-incidence.md`: fixed-path tilings, their all-tile-complement lines, and converse-merged tournament nodes. It is **not** the full single-arc-flip quotient of THM-588, `THM-588-metagraph-laplacian-algebraic-connectivity-four.md`. No integer spectrum or gap 4 is asserted here.

The closest proved mechanisms are the labelled score-complement identity in THM-646 (`THM-646-line-metagraph-score-complement-skeleton.md`), the leg/triangle-current identity in THM-649 (`THM-649-flip-parity-trilogy-affine-D-binomial-law.md`) and THM-790 (`THM-790-blue-parity-law-proved.md`), and the Hamiltonian-path-orbit fibre in THM-796. The hostile is a count-preserving reassignment of two blue lines at order five. The near miss is treating a generic symmetric weighted graph's PSD tests as tournament-specific realizability tests. The least-used coordinate here is the joint score distribution of the first and last vertices of Hamiltonian paths.

Terminology firewall: current rigidity below concerns only the cyclic-triangle value transported by complement. It is not the single-observer-path condition `H(T)=|Aut(T)|`, and it is not pure-blue status. Neither converse is inferred.

Useful names:

- **Endpoint transport inventory:** Hamiltonian paths counted by the scores of their first and last vertices, optionally restricted to blue paths.
- **Triangle-current-rigid node:** its complement-neighbour cyclic-triangle count is constant across all incident half-edges.
- **Complement leakage:** the exact covariance discarded by forgetting the half-edge after one complement step.
- **Regular-shell bottleneck:** all half-edges of a regular tournament must land in the almost-regular score shell.

The main new arrangement constraints are the endpoint inventory identities and the structural classification of triangle-current rigidity. The elementary compression and rank identities are infrastructure, not a claim of a new general matrix theorem.

## 1. Conventions

Fix `n>=3`, base Hamiltonian path `n -> ... -> 1`, tile set `S={(x,y):x-y>=2}`, and `X=F_2^S`, of size `N=2^m`, `m=binom(n-1,2)`. Let `C(t)=t+1` and let `R` be the staircase reflection. The merged-node map `pi:X->M` satisfies `pi R=pi`. A blue tiling is fixed by `R`; blue and black subsets are `C`-invariant. Write `F_u=pi^{-1}(u)` and `d_u=|F_u|`.

The half-edge count matrix is

`A_uv = #{t in F_u: C(t) in F_v}`.

Thus `A` is symmetric, `A 1=d`, and a loop line contributes **two** to `A_uu`. Write `A=B+K` for its blue and black parts, with row-mass diagonals `D=diag(d)`, `D_B`, `D_K`. Let `h_u=c3(T_u)`, which is well-defined after converse merging.

If a node is represented by an ordinary tournament `T`, set `a(T)=|Aut(T)|` and `kappa(T)=1` when `T` is self-converse, `2` otherwise. For a Hamiltonian path `P=(p_1,...,p_n)`, let

`e(P)=(s_T(p_1),s_T(p_n))=(a,b)`.

A path is called blue precisely when the positional reversal `p_i -> p_(n+1-i)` is an anti-automorphism of `T`. This criterion does not require choosing an arbitrary converse isomorphism.

## 2. Endpoint transport theorem

Let `E_T(a,b)` count Hamiltonian paths with endpoint scores `(a,b)`, and `E_T^B(a,b)` count only blue paths. Define `E_T^K=E_T-E_T^B`. Then for either colour `c`, and also after summing colours,

`sum_v A^c_uv z^(h_v-h_u) = (kappa(T)/a(T)) sum_(a,b) E_T^c(a,b) z^(a-b-1)`. (2.1)

This is an identity of Laurent polynomials. In particular every coefficient, not only the total degree or a parity, is prescribed by the endpoint inventory.

### Proof

Place the observer path at `n -> ... -> 1`. THM-646 gives

`s'_n=n-s_n`, `s'_1=n-2-s_1`, and `s'_i=n-1-s_i` for `1<i<n`.

Using `c3(T)=binom(n,3)-sum_v binom(s_v,2)` and `sum_v s_v=binom(n,2)` yields

`c3(Ct)-c3(t)=s_n-s_1-1=a-b-1`. (2.2)

Equivalently this is THM-790's triangle form of the leg law. The action of `Aut(T)` on ordered Hamiltonian paths is free: an automorphism fixing the entire ordered path fixes every vertex. Each path orbit supplies one fixed-path presentation, and its endpoint-score and blue properties are invariant on the orbit. This proves (2.1) for an ordinary isomorphism fibre after division by `a(T)`.

If `T` is not self-converse, its converse supplies the second ordinary fibre. Reverse the path and converse the tournament. Endpoint scores become `(n-1-b,n-1-a)`, leaving `a-b-1` unchanged; the blue set is empty. Thus the merged fibre contributes the factor `kappa`. If `T` is self-converse there is only one ordinary fibre. This proves all cases. QED.

### Full score-profile version

For a sorted score multiset `S`, let `[S]` denote its equivalence class under `s -> n-1-s`. For endpoint scores `(a,b)`, form `Phi(S;a,b)` by reflecting all scores, then replacing one reflected start score `n-1-a` by `n-a`, and one reflected end score `n-1-b` by `n-2-b`. The two vertices are distinct, including when `a=b`.

For every folded score profile `sigma`,

`sum_(v:[scores(T_v)]=sigma) A^c_uv`
` = (kappa/a(T)) sum_(a,b:[Phi(S;a,b)]=sigma) E_T^c(a,b)`. (2.3)

The pointwise score equations above prove the formula, with the same orbit argument. It is stronger than merely asking whether a pair of score profiles can occur. For an arbitrary destination-node set `J`, it also gives the necessary capacity bound obtained by summing the right-hand side over the score profiles present in `J`. Equality holds when `J` is a union of full score-profile classes. These are necessary constraints, not sufficient conditions for actual tournament-class matching.

The inventory can be computed without the metagraph: count Hamiltonian paths of one representative by endpoints, for example by a subset dynamic program. Thus (2.1) and (2.3) are not definitions that inspect the desired neighbour row.

## 3. Exact structural classification of triangle-current-rigid nodes

Order the strongly connected components of `T` as

`T_1 -> T_2 -> ... -> T_k`,

so all arcs between components point forward. Singletons count as regular tournaments.

**Theorem.** The following are equivalent:

1. The merged node represented by `T` is triangle-current-rigid.
2. `s_T(p_1)-s_T(p_n)` is constant over all Hamiltonian paths of `T`.
3. Both `T_1` and `T_k` are regular.

If `k=1`, this says simply: a strong tournament is triangle-current-rigid if and only if it is regular.

### Proof

Equivalence of 1 and 2 is (2.1). Every Hamiltonian path concatenates one Hamiltonian path in each component, in the displayed order.

A nontrivial strong tournament has a Hamiltonian cycle. Cutting a cycle just before any prescribed vertex realizes that vertex as a Hamiltonian-path start; cutting just after it realizes it as an end. A proof of the cycle fact is recalled below.

If `k>=2`, start choices in `T_1` and end choices in `T_k` are independent: choose the respective component paths and concatenate them with arbitrary Hamiltonian paths in the middle components. Therefore their score difference is constant exactly when the source-component scores are all equal and the sink-component scores are all equal. External contributions to scores are constant within a component, so this is precisely regularity of `T_1,T_k`.

If `k=1`, take a Hamiltonian cycle `v_1 -> ... -> v_n -> v_1`. If every path endpoint difference is a constant `q`, cutting each cycle edge gives `s(v_(i+1))-s(v_i)=q` for every `i`. Summation gives `nq=0`, so `q=0`; equality around the cycle makes all scores equal. Conversely regularity makes every endpoint difference zero. QED.

**Constant current.** If the source and sink component sizes are `a,b`, every incident half-edge has

`h_v-h_u=n-1-(a+b)/2`. (3.1)

This includes the strong regular case `a=b=n`, where the value is `-1`, and the transitive case `a=b=1`, where it is `n-2`. Middle components can be arbitrarily nonregular without changing rigidity. Regularity of the entire tournament is therefore not the correct condition in the nonstrong case.

**Elementary cycle fact (Camion).** Choose a longest directed cycle in a strong tournament. An outside vertex with both incoming and outgoing arcs to the cycle can be inserted at a transition between consecutive cycle vertices. Otherwise each outside vertex either dominates the entire cycle or is dominated by it. Strong connectivity supplies an outside edge `b -> a` from a vertex dominated by the cycle to a vertex dominating it. Replacing any cycle edge `x -> y` by `x -> b -> a -> y` extends the cycle, contradiction. Thus the cycle spans. Strong tournaments of size greater than one contain a cycle, since an acyclic tournament is a transitive total order.

### Current rigidity forces an independent set

For every `n>=4`, the full set of triangle-current-rigid nodes is independent, and none of them has a complement loop. This is stronger than the regular-node statement below.

Indeed, along an edge joining two current-rigid nodes their constant currents must be negatives of each other. A strong current-rigid node is regular and has current `-1`; it can exist only at odd order. Every nonstrong current-rigid node has odd end-component sizes `a,b` and positive current `n-1-(a+b)/2`. At odd `n>=5`, `a+b<=n-1`, so that positive current is at least `(n-1)/2>=2`, never `+1`. At even order there are no strong regular nodes, so every current-rigid current is positive. A zero current is impossible in either case. Thus neither an edge between current-rigid nodes nor a loop at one can occur. At `n=3` both nodes are current-rigid and adjacent, so the lower boundary is sharp.

Consequently, for either colour `c`,

`sum_(u triangle-current-rigid) d_u^c <= (1/2) sum_u d_u^c`. (3.2)

Every half-edge at a current-rigid node must be matched to a half-edge at a current-nonrigid node. This is a capacity restriction on the placement of the formulaic colour masses, not merely their total.

The historical theorem is Camion (1959); this classical input is independently proved here. A modern primary research source identifying the same input is Kühn–Lapinskas–Osthus–Patel, https://arxiv.org/abs/1303.4213. No new claim of priority for Camion's theorem is intended.

## 4. Regular-shell bottleneck

Let `n=2r+1` and let `u` be any regular node. Then every half-edge has

`h_v=h_u-1`, `scores(T_v)={r-1,r,...,r,r+1}`. (4.1)

Indeed all path endpoint scores equal `r`, so (2.2) gives `-1`; the score-profile transformation changes only two entries, by `+1` and `-1`.

Consequences:

- Regular nodes have no complement loops.
- No two regular nodes are adjacent; the regular set is independent.
- Their full degree mass, colour by colour, is confined to the almost-regular shell.
- If `T` has scores `{r-1,r^(n-2),r+1}`, with unique high-score vertex `v_+` and low-score vertex `v_-`, then

`sum_(v regular) A_uv = (kappa(T)/a(T)) HP_T(v_+,v_-)`. (4.2)

The colour-restricted version uses only paths of that colour. Formula (4.2) holds because the transformed score profile is regular exactly when the observer path starts at `v_+` and ends at `v_-`.

For `n>=5`, every tournament in the almost-regular shell is strong and nonregular, hence triangle-current-nonrigid. To see strong connectivity directly, a sink union of `k` strong components would have total outdegree `k(k-1)/2`. But any `k<n` vertices in the shell have total score at least `kr-1`, and `kr-1>k(k-1)/2` because `k(n-k)>2` for `n>=5`. This excludes such a sink union.

This is a shell mass constraint, not a claim that every pair of nodes in the two shells is adjacent. It does not characterize pure-blue tournaments.

## 5. Complement compression, exact rank, and the missing half-edge

Let `U` be the `N`-by-`|M|` fibre-indicator matrix. Then

`D=U^T U`, `A=U^T C U`, `Pi=U D^(-1) U^T`.

`Pi` is orthogonal projection onto fibre-constant functions. Consequently

`D +/- A = U^T(I +/- C)U >=0`,

`L:=D-A D^(-1)A = U^T C(I-Pi)C U >=0`. (5.1)

The normalized matrix `D^(-1/2) A D^(-1/2)` is a self-adjoint contraction. These are generic reversible compression facts. In fact any nonnegative symmetric matrix with row sums `d` already satisfies `D +/- A>=0`; these conditions alone cannot reject a correctly weighted abstract multigraph.

For a node observable `f`,

`f^T L f = sum_u sum_(t in F_u) [f(pi(Ct)) - (Af)_u/d_u]^2`. (5.2)

Thus `L` measures the information discarded when an actual half-edge is replaced by a uniform choice within its node. True double complement is the identity, whereas the node Markov surrogate has two-step operator `(D^(-1)A)^2`.

**Full equality.** `L=0` if and only if complement descends to a deterministic node involution, equivalently each row of `A` has only one nonzero entry. Equivalently the fibre-constant function space is invariant under `C`. For a particular `f`, equality in (5.2) only says the complement value of that observable is constant in every fibre; it does not imply full lumpability.

**Exact rank.** Let `q=|M|`, `c` be the number of connected components of the support graph of `A`, and `b` the number of those components that are bipartite (loops destroy bipartiteness). Then

`rank L=q-c-b`. (5.3)

Proof: `L=D^(1/2)(I-S^2)D^(1/2)` with `S=D^(-1/2)AD^(-1/2)`. The `+1` eigenspace has one constant-on-component vector for each component. The `-1` eigenspace has one bipartition-sign vector for each bipartite component. There are no other unit-modulus eigenvectors, since equality in the convex averaging step forces equal values, respectively opposite signs, across every support edge. This proves (5.3).

### Colour-resolved leakage

Use a pseudoinverse at zero colour masses and put `L_c=D_c-A_c D_c^dagger A_c`. Each `L_c` is PSD. If `b_u=sum_v B_uv`, `k_u=sum_v K_uv`, and `v_u=B_(u,*)/b_u-K_(u,*)/k_u` whenever both denominators are nonzero, then

`L=L_B+L_K+sum_(u:b_u k_u>0) (b_u k_u/d_u) v_u^T v_u`. (5.4)

This is the matrix law of total variance. It separates hidden variation inside each colour from variation caused by forgetting the colour itself. Equation (5.4) is generic; the endpoint inventories identify its tournament-specific values.

### The triangle equality boundary

Combining (2.1), (3), and (5.2),

`h^T L h = sum_u (kappa(T_u)/a(T_u)) sum_(P in HP(T_u)) [s(P_start)-s(P_end)-mean_u]^2`. (5.5)

The constant `-1` cancels in the variance. The contribution of a node is zero **exactly** when both end SCCs are regular. This is a concrete structural meaning for one direction of the leakage tensor.

For every `n>=4`, the partition is not complement-invariant. At `n=4` the unique strong nonregular tournament gives triangle leakage `14/5`. For `n>=5`, put an ordered chain of `n-4` vertices before the strong four-vertex tournament whose row-bit adjacency is `(12,9,2,4)` in zero-based coordinates. It has five Hamiltonian paths, with final score 1 three times and final score 2 twice. Its automorphism group is trivial; after adjoining the chain the merged node has mass 10 and contributes `12/5` to (5.5). Hence

`h^T L h >= 12/5 >0` for every `n>=5`. (5.6)

The strong-four witness has arcs `0->{2,3}`, `1->{0,3}`, `2->{1}`, `3->{2}`. Its Hamiltonian paths are `(0,2,1,3)`, `(0,3,2,1)`, `(1,0,3,2)`, `(2,1,0,3)`, `(3,2,1,0)`. These explicit data also verify the stated counts without a classification assumption.

## 6. A hostile surviving counts, colours, global currents, and score support

FINITE-EXACT at `n=5`, with the reproducible node labels in `n5.json`. Code uses zero-based vertices, tiles ordered by increasing `x` then increasing `y`, and bit 1 for `x->y`. Relevant nodes are:

| node | representative mask | fibre mass | blue mass | c3 | score multiset |
|---|---:|---:|---:|---:|---|
| 0 | 0x0 | 9 | 3 | 3 | 1,1,2,3,3 |
| 2 | 0x2 | 13 | 3 | 4 | 1,2,2,2,3 |
| 3 | 0x5 | 11 | 3 | 4 | 1,2,2,2,3 |
| 6 | 0xc | 3 | 3 | 5 | 2,2,2,2,2 |

Remove one blue line of each pair `(0,2),(3,6)` and add one of each pair `(0,3),(2,6)`. Every listed pair already occurs in the true graph. Each removed pair originally has multiplicity one and disappears, while each added pair originally has multiplicity one and becomes double. **The support itself is not preserved.** No previously forbidden endpoint pair is introduced.

This hostile preserves:

- every node's total, blue and black half-edge masses;
- all black lines and their paired-reflection structure;
- every diagonal entry and loop count;
- every parity law following from those masses;
- the full global joint triangle-count distribution on the two endpoints;
- the full global joint endpoint score-profile distribution;
- consequently every global triangle-current moment, not merely its parity or variance;
- nonnegative symmetric matrix compression/PSD conditions.

It nevertheless violates (2.1). Node 2's blue row changes from

`z^(-3)+z^(-1)+z` to `z^(-3)+2z`,

and node 3's from

`2z^(-1)+z` to `3z^(-1)`.

The blue endpoint inventory of each original tournament is unchanged, so neither altered row can be its actual complement row.

There is even an abstract lift preserving the same reflection and node partition: replace complement pairs `(0x16,0x29)` and `(0x1a,0x25)` by `(0x29,0x1a)` and `(0x16,0x25)`, leaving every other tiling pair unchanged. All four tilings are blue. This gives another fixed-point-free involution commuting with reflection, so all generic involution-compression constraints still hold. It is not the all-tile complement, and the endpoint transport identity detects exactly that missing tournament coordinate.

## 7. Exact controls and reproducibility

`probe.py` independently classifies fixed-path tilings by score-refined canonical labelling, folding converse only after computing both ordinary canonical codes. `verify.py` independently enumerates all Hamiltonian paths and automorphisms of one representative per node, verifies both colour-current and folded-score row identities, checks the SCC rigidity equivalence, and compares the covariance formula with the matrix leakage. The rank lower bound is exact modular elimination after reducing rational entries modulo a prime that divides none of the fibre masses; (5.3) supplies the matching graph-theoretic upper bound, so the recorded rank is exact, not floating point.

Commands from this packet directory:

`python -u probe.py`

`python -u verify.py`

`python -u verify_hostile_and_colour.py`

Results through order seven:

| n | nodes | triangle-current-rigid rows | current-nonrigid rows | leakage rank | triangle leakage |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 0 | 0 | 0 |
| 4 | 3 | 2 | 1 | 2 | 14/5 |
| 5 | 10 | 4 | 6 | 9 | 378284/6435 |
| 6 | 34 | 6 | 28 | 33 | 4411601212005792352/2803448856786435 |
| 7 | 272 | 17 | 255 | 271 | exact rational in `audit_n7.json` (approximately 71488.4358) |

All assertions pass. `verify_hostile_and_colour.py` verifies the full rational colour-covariance matrix identity through order six, the hostile, and an exact mask-by-mask match with an independently implemented Hamiltonian-orbit census classifier through order seven. It uses an explicit tile-order permutation and bit-complement gauge rather than assuming the two coordinate systems match.

The independent audit in `independent_transport.py` separately verifies the endpoint current and score-profile rows at `n=3..6` and the hostile, using a different classifier and orientation implementation.

The scalar identity `h^T(D-A)h=N(n-1)/4` is also checked; it is a second-moment corollary of the already-proved THM-790 global current generating function, not a new result.

## 8. Preserved/lost ledger and open-question relevance

Source: Hamiltonian-path orbits above ordinary tournament classes. Target: coloured weighted complement incidence on merged nodes. Map: forget the ordered observer path, retain node and complement-neighbour node. Preserved: endpoint-current inventory, score-shell capacities, masses and colour. Lost: the particular observer path, exact endpoint identities, and neighbour identity inside a common score shell. Sidecar: endpoint-score inventories restore the displayed row identities but do not restore the entire matching.

The hostile shows exactly why global formula counts need arrangement rules: even preserving every global current moment and all permitted score pairs leaves wrong node rows. Triangle rigidity further gives a structural answer to a visual question: when an entire node's incident lines lie on one triangle-level difference, its first and last strong components must be regular.

The following remain open here:

1. Which additional low-cost path inventories determine matching inside one destination score shell?
2. Can endpoint transport and full Hamming/Walsh lift constraints jointly characterize finite candidate realizability?
3. What is the minimum continuation-stable refinement that removes complement leakage while respecting deletion/extension operations?
4. How do pure-blue node classifications sharpen the regular-shell bottleneck at higher orders?

No assertion is made about LRC(14), another prize problem, universal sufficiency of these gates, or a spectrum/gap inherited from THM-588. The contribution is a rigorous realizability obstruction and an exact equality classification for the actual fixed-path metagraph.
