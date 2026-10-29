# Tournament classification: inherited frontier and precise next targets

Source snapshot: `3b6e3ef82382da222fc59aa02914eea72d5a105b` (2026-10-29).
This is a bounded source-recovery map, not a new proof or literature-priority claim. Repository theorem status is inherited, not independently re-proved here. Sources below were read by exact `git show` from the pinned snapshot; the companion JSON validates their blob bytes. No publication is performed by this recovery pass.

## 1. What “classification” currently means

There are four distinct deliverables:

1. **CITED + PROVED elementary reconstruction:** canonical strong-module decomposition gives an exact ordinary-tournament classifier, retaining ordered linear nodes and child-coloured prime quotients. The decorated tree reconstructs the tournament up to isomorphism. This is already complete as a framework; it does not supply a concise bank of every prime species.
2. **FINITE-EXACT enumeration:** literal representatives exist through order ten (9,733,056 classes), according to `hn_metagraph_information_20261026` and the routed `consecutive_flags_order10_20261028` packet. This recovery did not replay the order-ten census. A000568 coefficients through 350 are inherited aggregate counts, not representatives at those orders.
3. **PROVED predicate-specific classifications:** small Hamilton presentation count and converse-distinguishing exceptions have genuine all-order atom lists, described below.
4. **OPEN operation-compatible classification:** full general-substitution metagraph grammar, concise prime-species atlas, and smallest carrier for selected operation menus remain separate questions. A counting invariant, rank, or ordinary object classifier does not automatically classify operations on presentations.

The useful finite-simple-group analogy is “irreducible quotients plus extension/action data.” Prime tournaments are not simple groups. Their automorphism groups have odd order; their action, not only their group order, controls child attachments. The analogy is a research architecture, not an achieved unrestricted classification.

## 2. Classical inputs with checked source scopes

- **Gallai modular decomposition; cyclic-triple realization.** Boussaïri–Chergui–Ille–Zaidi, [3-uniform hypergraphs: modular decomposition and realization by tournaments](https://arxiv.org/html/1805.04607v1), recalls tournament strong-module decomposition and proves the realization machinery. Prime tournament cyclic-triple data determines its orientation up to converse (Theorem 24, crediting the earlier Boussaïri et al. result). Realizations share strong modules. Repository `tournament_composition_classification_20261026/CLASSIFICATION.md` uses Theorem 10, Theorem 24, Proposition 51 and Corollary 57, giving labelled lift count `2^(prime nodes) product(linear arity!)`; quotient by the hypergraph automorphism action, not equal orbit sizes.
- **Critical prime families.** Schmerl–Trotter, *Discrete Mathematics* 113 (1993), 191–205, [DOI](https://doi.org/10.1016/0012-365X(93)90516-V): the critical tournaments of order at least five are exactly the odd-order families `T_n,U_n,W_n`. This statement was checked also as Theorem 23 in the preceding primary-author paper. “Critical” means every one-vertex deletion is decomposable; it does not mean every prime belongs to these three families.
- **A less-used next stratum.** Belkhechine–Boudabbous–Dammak, [Les tournois (-1)-critiques](https://arxiv.org/abs/1007.2726), published 2007: classifies indecomposable tournaments with exactly one noncritical vertex; orders are odd and the abstract gives `3m-15` types at each odd `m>=7`. This is a promising existing atlas layer to import rather than rediscover. The full family definitions were not audited here.
- **Switching-native decomposition is already classical.** Bui-Xuan–Habib–Limouzy–de Montgolfier, [Unifying two Graph Decompositions with Modular Decomposition](https://arxiv.org/html/cs/0604065v3), Theorem 5, Proposition 10, Theorem 7: bijoin/umodule decomposition for tournaments; root-normalized Seidel switching converts appropriate complementary umodules to modules after deleting the root; diamond-free, locally transitive, totally bijoin-decomposable, and twin/antitwin-constructible tournaments coincide. This is the correct classical boundary for a signed rank-one cut project. Its root-deleting normalization is not automatically the same operation as keeping all vertices and replacing a skew matrix by `DSD`. The author's PDF timed out; full arXiv HTML was successfully inspected. No novelty claim is made for these equivalences.


- **Actual finite-group classification for switching classes.** Babai–Cameron, [Automorphisms and enumeration of switching classes of tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v7i1r38/pdf), *EJC* 7 (2000), R38, Theorem 5.2: a finite abstract group is the full automorphism group of a tournament switching class iff its Sylow 2-subgroups are cyclic or dihedral. Theorem 5.1: a given permutation group is contained in some switching-class automorphism group iff its Sylow 2-subgroups are cyclic or dihedral and act semiregularly. The distinction is **full abstract realization versus containment in a specified permutation action**. These are classical theorems directly realizing part of the requested group-classification ambition, not a prime-tournament atlas or a classification of every switching class. Both theorem statements were checked in the primary PDF, page 11. Existing repository use is explicit in `THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md`, audit field: Theorem 5.1 supplies semiregular 2-subgroups for its `n=2 mod4` argument. `THM-1460-babai-cameron-zero-at-all-odd-n.md` separately studies Remark 7.4 via canonical score-parity anchors; this source is not new to the repository.

## 3. Repository structural map

All paths in this section are relative to the repository root.

### Exact modular reconstruction and the hostile that prevents scalar gluing

`05-knowledge/results/tournament_composition_classification_20261026/CLASSIFICATION.md` proves, relative to the cited decomposition, the decorated-tree iff and

`Aut(Q[T_i]) = (product Aut(T_i)) semidirect K`,

where `K` is the child-colour stabilizer in `Aut(Q)` and the decomposition partition is canonical. For arbitrary substitutions this only identifies the partition-preserving subgroup until intrinsicality is proved.

The same file retains the sharp 12-vertex **81/27 hostile**: `A=C3→1`, `B=A*`, and `C3[A,A,A]` versus `C3[A,A,B]`. Both have 67 identical labelled cyclic triples and `H=462915`, but automorphism orders 81 and 27; hence path-class counts 5715 and 17145. Current selectors retain attachment-orbit counts 176/512 and transporter counts 27/0. An observer's fibre is a permutation action, not merely its cardinality.

### Genuine complete predicate strata

- `01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md`: **PROVED current audited addendum**, with Gluck and Feit–Thompson explicitly cited. No subset has trivial stabilizer under `Aut±` iff `T=X[C3^[k]]`, `X∈{C3,T5,P7}`, `k>=0`, for `|T|>=2`. The original unsupported completeness sentence is not validated retroactively; the later structural necessity proof and MISTAKE-650 lineage matter.
- `05-knowledge/results/metagraph_classification_20261021a/README.md`: **PROVED** `h=H/|Aut|=3` iff the unique strong nontrivial atom is `R5`; arbitrary tournaments are ordered SCC stacks on `{K1,C3,R5}` with exactly one `R5`. Inherited `h=1` atoms are `{K1,C3}`. This is not the classification of all Hamilton counts or all prime tournaments.
- `01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md`: a prime parent of order `q>=3` has `2^q-2-2q` prime labelled attachments, with the other patterns sources/sinks or homogeneous-pair extensions. Critical exceptions and root/card incidence are explicit. Its 7,921,973,322 rows are a presentation cover, not a new ordinary-class census.

### Signed Hamilton data and distinct sign operations

`01-canon/theorems/THM-1936-signed-redei-join-multiplicative.md`: with block-compatible labels, signed count `R=Σ sgn(path)` is multiplicative under ordered join; relabelling multiplies R by the label-permutation sign, so only `|R|` is an ordinary invariant. The strong four-tournament has `H=5,|R|=1`: do not use the corrected old H=3 sentence.

`01-canon/theorems/THM-1966-signed-redei-count-independent-invariant-n7.md`: exact order-seven witness shares adjacency characteristic polynomial and `H=81` but has `|R|=1` versus 17. Thus spectrum + Hamilton count is not a full signed classifier.

Keep these operations separate: global converse; relabelling parity; cut switching `DSD`; free-arc complement relative to a selected Hamilton path; permutation of names of disjoint runs. The path complement keeps consecutive path arcs and reverses all other arcs. It is not global converse. Signed Rédei gauge cannot be silently identified with symmetric-character colouring or a switching-invariant scalar.

### Operation-closed carriers: known sufficiency versus unresolved minima

`minimal_blue_selectors_20261028/README.md` proves the minimum for its singleton/native-complement/fresh-re-marking language: two endpoint-rooted relations plus B, whose boundary part has `4n-8` states for `n>=5`. It does not solve arbitrary strong-atom substitution.

`moving_selector_atoms_20261028/README.md` proves an all-size read-only-cursor carrier and a stronger failure of marginal scans. In child `C3[C3,1,1]`, ordered singleton runs `(3,0,1,2,4)` and `(3,0,2,1,4)` agree under all individual cursor probes but assemble to Hamilton counts 2,899,675 and 2,906,091. One joint parity bit is missing in this fibre. Persistent C3 marks can already reconstruct every old arc; full presentation reconstruction in that richer language is therefore not an open target. The narrower read-only carrier's universal reconstruction/minimal uncoloured composition remains OPEN. Over an aligned source-target pair, the faithful residual is component ordering, with `k!/|G|` orbits.

`session94_recursive_ports_crt_clocks_20261029/README.md`: full cycle indices distinguish attachment behaviour missed by every fixed-point moment. Regular `C9` and `C3²` have identical named-port moments but different three-subset orbit counts 10/12. Paired palettes must share the same acting element; nine-vertex binary paired count 192 is not independent quotient 128.

`session95_niche_recurrences_20261029/README.md`: the circulation-boundary Schur matrix is a closed substitution grammar, unlike its scalar series. Children with identical scalar flow series yield 9470 versus 9462 at degree 12 in a common C3 context. This is a concrete operation-compatible classification model, but not a Hamilton grammar.

`session98_tournament_growth_20261029/composition/QUOTIENT_EXTENSION.md`: finite-surplus cycle-index jets have exact observer quotients. At surplus 12, ambient dimension 70, formal rank 57, ordinary-tournament-palette rank 33. These are minimal *linear ambient* quotients, not minimal realizable action atlases or metagraph carriers. At surplus 3 actual palettes already impose a kernel invisible in a formal-colour calculation.

### Inherited Cayley aggregate: connectivity is already proved at the correct quotient

`05-knowledge/results/tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md`, section 8, Theorem 3, already proves that aggregating legal Hamilton-complement incidences over labelled switching classes gives the weighted Cayley matrix on `Q=X/B` with generators `[E\P]`. Each undirected Hamilton path contributes weight two; row sums are `n!`. For `n>=3` the support is connected unless `n=1 mod4`, when it has exactly two components. It is bipartite exactly in the connected `n=3 mod4` case; even-order graphs and the individual `n=1 mod4` components are nonbipartite. Distinct paths give distinct generators for every `n>=8`, proved by a cut-degree argument including the order-eight endpoint exception. At seven the inherited exact census has 2520 undirected paths, 1260 generators, and weight four.

This is an **aggregate incidence quotient, not an equitable quotient**: representatives inside one switching class can have different Hamilton counts. Consecutive quotient steps need not lift from the same intermediate tournament. It therefore proves neither native black-metagraph connectivity nor fixed-score-fibre connectivity. The sharp two-step obstruction used in the current session is an inherited corollary of this generator-injectivity theorem, not a fresh discovery. Potential current additions are the refined seven-vertex legal skeleton/colours and general word-length bounds; their proof and status belong to the session's separate computation/theory artifacts, not this provenance map.

### Higher precision, not merely rank: newest native gluing frontier

`overnight_research_20261029/THREE_RUN_GLUE.md`: the regular-S3 chart `(u,v,W)` has index 54 and Smith factors `(1,1,1,3,3,6)`; exact congruences couple sign/trivial/standard sectors. Independent rational Fourier blocks lose integral data. Run-name S3 is not a tournament automorphism group.

`four_vertex_glue_fermat_20261029/README.md`: for the literal 60-by-60 three-run matrix M of T5, native Smith form is `1^52,3^2,6^2,0^4`. The tetrahedral kernel has four directions, Gram determinant 125, and its native lattice is six times the unglued kernel. For all r>=1 the packet proves

`0 → coker(M^r) → coker(A^r) → (Z/6)^4 → 0`,

where A is the separate-sector matrix. At r=2 the dyadic extension is nonsplit. The defect carries an irreducible C5 action at both bad primes. Every-power torsion order is `2^(6r-4)3^(8r-4)5^(3r-3)` and a degree-37 minimal common matrix recurrence is supplied. Most relevant for research: the native/unglued mod-m image ratio is

`2^(4[v2(m)>=r]) 3^(4[v3(m)>=r])`.

Thus a fixed prime-field rank can be completely blind to glue first visible at digit r. Powers of a fixed linear response are not disjoint-cover composition on growing tournaments. Actual native complement exposes all four inherited T5 kernel directions; the stacked observer has rank 60. A universal fixed-kernel compression is already ruled out.

## 4. Five precise research targets

1. **Signed-cut / bijoin classification interface (anchor, root pursuing).** Identify exactly which signed cross-block rank-one cuts are classical bijoins; normalize at a named root and track both the deleted root and switching signs. Deliver a reconstruction theorem plus hostile example showing ordinary modular decomposition alone is insufficient. Compare every purported new statement with Bui-Xuan et al. The next genuinely new question is compatibility with native Hamilton-path complement or endpoint response, not rediscovering twin/antitwin equivalence.
2. **Mixed-word Smith glue (niche).** For a specified legal T5 complement partner U on the same labels, compute native `M_U M_T` and sector-chart product, including primitive left/right kernels, exact Smith factors and induced quotient action. Test whether the fixed-source cokernel exact sequence survives; if not, give the first failed map and a replacement. Controls: fixed T powers, order-four unimodularity, all four zero directions exposed under complement. Stop with a proved scoped rule or a minimal explicit obstruction, not ranks alone.
3. **Read-only versus persistent-mark boundary (anchor alternative).** Fix K_scan's exact allowed cursor/reset/assembly syntax. Determine whether the run-order parity hostile survives all permitted uncoloured futures; if not give a shortest separating context. Do not claim novelty for persistent-mark reconstruction, already proved. This asks for a concrete behavioural quotient rather than unrestricted “minimum state.”
4. **Deletion-defect prime atlas (structural niche).** Import critical and (-1)-critical family generators, verify low orders against the inherited prime census, and compute rooted augmentation stabilizers and endpoint profiles on these strata. A precise new target is whether the root-coloured response recurrence closes within a finite family of attachment types. Raw `2^q-2-2q` counts do not remove root/card multiplicities.
5. **Palette/attachment refinements beyond scalar growth (wildcard).** On a fixed prime arity, restrict Session98's observation kernels to the actually realized cycle-index span, then add one named-port or transporter test. Measure and prove which kernel directions survive. Use C9/C3² and 81/27 as hostile controls. This can yield a genuinely smaller correct query compiler without pretending to reconstruct individual metagraph neighbors.

## 5. Least-used recovered sidecars and exclusions

The old unsigned K4 clutch `THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch.md` is a genuine integer-lattice sidecar, but unsigned K4 edges are not oriented edge cochains or the six orders of three runs. Their respective gluing indices 24,16,54 must not be conflated. The June `tournament_contradiction_grammar_codex_s260.out` already organized certificate families, but is a historical computation receipt, not current proof authority; its candidate H-gap and LRC links require fresh canon checks before reuse. Low-deletion-defect classification is more actionable than re-running a generic invariant zoo.

No claim here closes unrestricted metagraph connectivity, a concise prime atlas, H>=disc, LRC(14), or a global compression theorem. All finite quantities retain their original observation universes.
