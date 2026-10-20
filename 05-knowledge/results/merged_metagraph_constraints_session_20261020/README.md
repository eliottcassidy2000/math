# Merged metagraph: structural rules beyond counts

Status: PROVED for the structural deductions stated below; FINITE-EXACT where an explicit census range is given. The complete pure-blue atom alphabet is settled, but the complete metagraph arrangement is not determined by these constraints.

## 1. The object and four predicates that must stay distinct

A node is a tournament isomorphism class merged with its converse. A tiling retains a fixed directed Hamiltonian path, and a literal complement line pairs a tiling with the tiling obtained by flipping every nonpath arc. Reflection permutes tile coordinates without flipping bits. A line is blue when its tilings are reflection-fixed, black otherwise. Loops and parallel lines remain literal instances; a loop contributes two half-edges.

- **Self-converse** concerns isomorphism to the reversed tournament; “SC” must not be reused here for strong connectivity.
- **Pure-blue**, or a “point of symmetry”, means every tiling in the node is blue.
- **Path-rigid** means `H(T)=|Aut(T)|`, equivalently one ordinary tiling in the class fibre. This does not mean current-rigid.
- **Triangle-current-rigid** means every incident complement half-edge has the same cyclic-triangle change. This is controlled by the endpoint components, not by a singleton fibre.

Regular cyclic `T5` is pure-blue with fibre three and has three distinct blue neighbors. It is not path-rigid or a leaf. These distinctions are useful constraints, not terminology for an already solved classification of adjacency.

## 2. Strong core, independent periphery, and exact cuts

For every `n>=3`, the partner of a nonstrong tournament is strong: the complement operation reverses the nonpath endpoint arc and closes the retained path into a Hamiltonian cycle. Thus nonstrong nodes form an independent periphery and every loop lies in the strong core.

Let `T_n=2^binom(n-1,2)` and `B_n=2^floor((n-1)^2/4)` count all and blue tilings; write `s_n,b_n` for their strong subsets, with `T_0=B_0=1`. Unique strong-component order gives

    T_n = sum_(k=1)^n s_k T_(n-k),
    B_n = b_n + sum_(k=1)^floor(n/2) s_k B_(n-2k).

The literal periphery/core crossing counts are exactly `T_n-s_n` in total and `B_n-b_n` in blue. At order eight these are `64,648` and `70`. These are forced arrangement cuts, not just global mass totals.

The unique transitive blue pipe lands mixed at every `n>=4`, with triangle jump `n-2`; at `n=3` its partner is the pure-blue 3-cycle. Every even-order pure-blue node is a singleton-fibre blue leaf attached to a mixed node for even `n>=4`.

[Proofs and all-order recurrences](../merged_metagraph_purity_and_strong_core_20261020/README.md).

## 3. The paired-black quotient retains a charge

Reflection pairs literal black lines. Halve these pairs to obtain a multigraph `Q`, with loops retained. If a merged node has total mass `M_u` and blue mass `B_u`, then

    deg_Q(u)=(M_u-B_u)/2,
    chi_u=deg_Q(u) mod 2.

Every non-self-converse node has charge one. For any shore `S`, the literal black cut satisfies `|delta_K(S)| = 2 sum_(u in S) chi_u (mod 4)`. Both shores of a `Q` bridge have odd total charge. The cycle rank splits as `beta(K)=|E(Q)|+beta(Q)`: the extra term is the forced reflected-pair digon space, not unexplained tournament topology.

Let `a_0=1`, `a_j=0` for `j<0`, `a_n=a_(n-1)+a_(n-3)`, and let `p_n` count palindromic `{1,3}` words of weight `n`. Two explicit leaf families prove

    bridges(Q_n) >= (a_n-p_n)/2 + [n even] a_((n-4)/2),   n>=5.

At `n=4`, the two constructed leaves share one bridge; they must not be counted twice. Equality and the claim that every bridge is a leaf edge are finite-exact through eight, not all-order theorems. Active blue and black connectivity is also only finite-exact through eight.

[Quotient, bridge proof, census, and independent witnesses](../merged_metagraph_black_quotient_bridges_20261020/README.md).

## 4. Endpoint inventories constrain destination rows

The Hamiltonian-path endpoint-score inventory of one tournament prescribes, separately in each colour, the entire row of cyclic-triangle changes and folded score profiles of its complement neighbors. This is stronger than allowing a destination with a compatible score sequence: it fixes the available multiplicity in each profile.

Triangle-current rigidity holds exactly when the first and last strong components are regular. If their sizes are `a,b`, the constant current is `n-1-(a+b)/2`; a strong regular tournament has current `-1`. For every `n>=4`, all current-rigid nodes form an independent set, so at most half of each colour's total half-edge mass can lie there. The `n=3` pair is the sharp exception.

A regular node can meet only the almost-regular score shell. Exact covariance and rank identities quantify the lost half-edge information when complement moves are compressed to nodes. A five-vertex hostile preserves colours, degrees, loops, and global joint triangle/score distributions but violates an individual node's endpoint-transport row.

[Endpoint transport, current rigidity, and complement leakage](../merged_metagraph_endpoint_transport_20261020/README.md).

## 5. Distance layers must satisfy matrix moments together

For ordered fibre-to-fibre Hamming-distance counts `W_d`, the Krawtchouk transforms are integer Gram matrices:

    Q_j = sum_d K_j(d) W_d = Z_j Z_j^T >= 0.

Reflection supplies stronger rank bounds than the ambient cube alone. At the top level the entire alternating-distance matrix has rank one, `Q_m=z z^T`, where `z` is the nonzero tiling parity charge. It is determined by the transitive node's alternating-distance row and is linked explicitly to a sign-corrected signed Rédei count.

A separate five-vertex hostile preserves every layer degree, diagonal distance distribution, entry parity, pair total, and even the entire complement graph, yet a two-by-two Gram minor has determinant `-64`. Thus adding all these scalar constraints still does not make a candidate arrangement realizable. The matrix constraint detects an obstruction that those summaries miss.

[Distance-moment proof, exact audits, and hostile](../merged_metagraph_distance_moment_constraints_20261020/README.md).

## 6. The pure-blue alphabet is completely classified

The all-order rigid-strong classification gives exactly the singleton and `C3`; hence all rigid classes are ordered `{1,C3}` stacks and the rigid self-converse classes are their palindromes. This removes the historical above-sixty-two gap without the old auxiliary group lemma ML.

A different mechanism closes the pure-blue alphabet: every nontrivial pure-blue strong tournament is an odd circulant; explicit prefix reroutings and unit-step reindexings then exclude every odd order at least seven. The strong pure-blue atoms are exactly

    {single vertex, C3, regular cyclic C5}.

Regular cyclic `C5` is the old `T5`, not the unrelated exceptional tournament sometimes denoted `T(5)` in bypass literature. Its fibre is three, so the distinction between pure-blue and path-rigid remains essential even after classification.

The center law now gives the exact ordinary generating function for all pure-blue isomorphism classes, including the empty object:

    P(z)=(1+z+z^3+z^5)/(1-z^2-z^6).

With `a_0=1`, `a_j=0` for `j<0`, and `a_r=a_(r-1)+a_(r-3)`, the all-order counts are

    pure-blue(2r)=a_r,
    pure-blue(2r+1)=a_(r+2).

The even formula includes the auxiliary empty object at `r=0`. These are unlabelled class counts, not tiling masses. The old three-letter formula is now a theorem, not an extrapolation of the census. The classification and center law also force every pure-blue node's literal line-degree to be one or three in the metagraph range `n>=3`.

### The pure-blue neighborhood is now exact

Every pure-blue node at `n>=4` has either one or three **distinct mixed neighbors**, with one blue line to each. There are no pure-blue loops or pure-blue-to-pure-blue lines. At `n=3` the two pure-blue nodes instead share their single line, a necessary exception.

If `L_n` and `J_n` count the leaf and trivalent pure-blue nodes, then

    L_(2r)=a_r,                 J_(2r)=0,
    L_(2r+1)=a_r+a_(r-1),      J_(2r+1)=a_(r-2).

The three-neighbor statement is not inferred merely from degree three. In every `A -> C5 -> A^op` branch, the central five vertices of the complement partner form its unique middle-score shell. The three induced center types have Hamiltonian-path counts `11,15,13`, so their full partners cannot become isomorphic or converse-isomorphic.

Consequently, writing `beta_n=2^(floor((n-1)^2/4)-1)` for all blue lines, exactly `L_n+3J_n` are pure-blue-to-mixed and exactly `beta_n-L_n-3J_n` are mixed-to-mixed at every `n>=4`. The branch count and placement are determined even though the full mixed wiring is not. The identity `p_(2r+1)=p_(2r+4)` also has a direct normal-form suffix bijection; it is not a cross-order graph isomorphism.

[Exact counts, distinct-neighbor proof, and line budget](../pure_blue_atom_classification_20261020/PURE_BLUE_COUNTS_AND_BRANCH_GEOMETRY.md).

The [uniform prefix-rerouting proof](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md) uses explicit violated arcs, not a large finite sweep or asymptotic estimate. The [independent rigid-strong proof and finite-atom method](../merged_metagraph_finite_pure_atoms_20261020/README.md) remain useful: the latter proved an intermediate order-131 bound and exhaustive circulant exclusions through 31. Those are valid separate results, not the current classification frontier or prerequisites for the uniform exclusion.

## 7. Odd-order blue spectra come from auxiliary cubes

At every odd `n>=3`, the self-converse classes admit a fixed-point-free involution preserving the number of blue tilings and reversing cyclic-triangle parity. Therefore the blue-count spectra on the two triangle-parity sides agree at every odd order, closing the former finite-observation question.

There is a stronger structure. Delete the paired edges of an involutory anti-automorphism to obtain a **hollow frame**. Its automorphism group has odd order, so each orbit of missing pairs has odd size. Toggling these orbits yields free auxiliary cubes `(C2)^r` of self-converse classes. Blue count is constant on each cube, and triangle parity is a checkerboard. For `n=2k+1`, `r>=1` and `r=k (mod 2)`.

Thus every blue-count multiplicity is even at odd order; when `n=1 (mod 4)`, `n>=5`, every full multiplicity is divisible by four and each shore multiplicity is even. The order-eight multiplicity-three hostile is unaffected: the repaired theorem is explicitly odd-order. The proof also retains an exact affine formula for the integer triangle count in the missing-edge bits, stronger than parity alone.

The map is genuinely on unmarked classes, via conjugacy of anti-involutions, but these cubes are **not** claimed subgraphs or automorphisms of the actual blue metagraph. The construction preserves a centralizer and a fixed-reversal path set, not necessarily the full automorphism group, Hamiltonian-path count, total fibre mass, or purity. At `n=7` it maps the pure-blue `C3 -> 1 -> C3` node of mass one to a mixed node of mass 45, both with blue count one; the pair is not an actual blue complement edge.

[Odd-spectrum theorem, hollow-frame cubes, and explicit loss boundary](../odd_blue_spectrum_hollow_cubes_20261020/ODD_BLUE_SPECTRUM_AND_HOLLOW_CUBES.md).

## 8. What remains beyond these constraints

The metagraph's actual cross-class matching is not determined by counts, colours, currents, or the constraints proved here. No sufficiency theorem is claimed. Useful residual targets include:

1. A realizability criterion coupling the endpoint inventory, distance moments, and defect-labelled factor incidence, without erasing the half-edge sidecar.
2. All-order bridge and connectivity laws beyond the exact order-eight census.
3. How the auxiliary hollow-frame cube structure constrains actual complement incidence. A class-level pairing preserving blue count does not automatically preserve the blue graph.
4. The smallest continuation-equivalent state for recursive constructions. A complete static fingerprint need not be a recursive state.

The historical 3-power cap on blue fibre counts is false at order eight (`39>27`), and a direct order-nine witness with `(H,|Aut|,B)=(2659,1,39)` refutes even its odd-order-only restriction. Universal even multiplicities of blue-count values are also false there (blue count one occurs three times). These corrections preserve the original parity laws; the repaired odd-order statement is now proved by the separate auxiliary-cube mechanism.
