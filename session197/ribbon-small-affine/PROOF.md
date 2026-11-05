# No fixed-block ternary or base-four affine clock at span three

Session 197. This extends the binary append-block obstruction in `../ribbon-periodic/PROOF.md` to multipliers three and four. The conclusion is all-order and covers append blocks of **every length**, using a structural reduction followed by finite exact matrix certificates. It is not a proof or disproof of span-three spectrum completeness.

## 1. Statement and conventions

Use the same span-three tournament and mask-word convention as the ribbon packets. For each `m` in `{3,4}`, any fixed prefix `P`, nonempty word `W`, and fixed **nonzero** constant `c`, the sequence

`H_k=H(P W^k)`

cannot be eventually nonconstant and satisfy

`H_(k+1)=m H_k+c`

for every sufficiently large `k`.

The proof does not claim a corresponding obstruction for arbitrary state-dependent append rules, reordering, or modifications of old edges. The restriction `c!=0` is deliberate. It is exactly the relevant case for unit-offset ternary rules and a two-step cycle of unit-offset binary rules.

## 2. Small path counts force a very simple cell

**Lemma.** A tournament with at most four Hamiltonian paths is an ordered join of singleton components and at most one directed-triangle component. Its count is therefore one or three.

First, every strong tournament of order at least four contains a strong four-vertex subtournament. Start with a directed triangle. If an outside vertex is mixed toward it, adding that vertex gives such a strong four-vertex graph. Otherwise all outside vertices are uniform toward the triangle. Strongness forces a vertex `x` dominated by the triangle and a vertex `y` beating the triangle, with `x -> y`. Choose two triangle vertices `u -> v`; then `y -> u -> v -> x -> y` is a four-cycle, so those four vertices induce a strong tournament.

Every strong four-tournament has five Hamiltonian paths. This can be checked on the four choices of diagonals of a spanning four-cycle, which are isomorphic. For one choice the five paths, in zero-based cycle labels, are

`0123, 1230, 1302, 2301, 3012`.

The verifier also enumerates all 64 labelled four-tournaments and checks all 24 strong ones directly.

Induced-subtournament path monotonicity now shows that any strong component with at least four vertices contributes at least five paths. Every nontrivial smaller strong component is a directed triangle and contributes three. Counts multiply over ordered strong components, so two triangle components contribute nine. The lemma follows, without assuming any global spectral gap or general parity theorem.

## 3. Separated boundary components give exact multiplication

Let `G=T(W)` be a single cell, of length `L>=7`. Write `A` for its first three vertices and `Z` for its last three. Suppose no directed path in `G` goes from any vertex of `Z` to any vertex of `A`.

The strong components of a tournament are totally ordered. The supposition means that every component meeting `A` is strictly before every component meeting `Z`. Split their order into

`I -> M -> J`,

where `I` ends at the last component meeting `A`, `J` starts at the first component meeting `Z`, and `M` contains the intervening components. Thus `A subset I`, `Z subset J`, and `M` may be empty.

All backward arcs between consecutive cells go from the next cell's first three vertices to the preceding cell's last three vertices. Therefore all such arcs are contained inside the boundary graph

`D_j = J_j union I_(j+1)`.

Every other arc between the following factors points forward:

`(P union I_1) -> M_1 -> D_1 -> M_2 -> D_2 -> ... -> M_k -> J_k`.

This is an actual ordered-join decomposition, not a formal factorization of response coordinates. Any backward edges from the first cell into the prefix also lie inside `P union I_1`, because `I_1` contains the first three new vertices. Nonconsecutive cells have only forward cross arcs. Every `D_j` is the same labelled-interface graph up to shifting indices.

Using the empty-product convention `H(empty)=1` only for a possibly empty `M`, we obtain for `k>=1`

`H(P W^k)=H(P union I_1) H(M)^k H(D)^(k-1) H(J)`.

In particular

`H_(k+1)=q H_k`, where `q=H(M)H(D)` is a fixed positive integer.

An eventually nonconstant sequence of this form cannot also satisfy `H_(k+1)=mH_k+c` with fixed `c!=0`: the equation `(q-m)H_k=c` would force constancy unless `q=m`, in which case it would force `c=0`.

## 4. Why no long block remains to be searched

The deterministic merging of Hamiltonian paths from disjoint vertex parts, proved in `../ribbon-compression/PROOF.md`, implies

`H(P W^k) >= H(T(W))^k`.

Thus a nonconstant affine clock of either multiplier must have `H(T(W))<=4`. By the preceding lemma, its cell is a join of singletons and at most one directed triangle.

If `L>=12`, a first-three vertex and a last-three vertex differ in natural index by at least `L-5>=7`. Their direct arc is forward. If there were a directed path back from the latter to the former, they would lie in the same nontrivial strong component. That component would have to be the single directed triangle, but every triangle in a span-three ordering has index diameter at most six. This is impossible.

Hence every block of length at least twelve is excluded either by growth at least `5^k` or by the exact multiplicative decomposition of Section 3. There is no extrapolation from a finite word bound.

The same growth argument applies to every cyclic rotation of `W`: after discarding bounded end fragments, `P W^k` contains `k-O(1)` disjoint cells of any rotated phase. If a rotated cell has at least five paths, it too rules out an affine clock of either multiplier.

## 5. The remaining exact spectral certificates

The finite matrix representation is the same independently verified one as in `../ribbon-periodic/representation.json`. Its integer matrices satisfy

`M_a B = B R_a`, and `B f = F`,

where the original nonnegative frontier has 119 states and the represented output has 62 coordinates. The checker re-verifies these identities before using the representation; no minimality assertion is needed.

An eventual sequence `alpha m^k-c/(m-1)`, with both `alpha` and `c` nonzero, requires **both** one and `m` as eigenvalues of `R_W`. Indeed applying its characteristic-polynomial recurrence gives

`alpha p(m)m^k-c p(1)/(m-1)=0`

for all sufficiently large `k`, forcing `p(m)=p(1)=0`.

The finite certificate checks the two multipliers separately and proves that no remaining word has both required eigenvalues:

### Lengths one through six

Every cyclic rotation class is checked, without any path-count filter. The numbers of classes are

`8, 36, 176, 1044, 6560, 43800`.

For each of these 51,624 classes, at least one of `det(I-R_W)` and `det(mI-R_W)` is nonzero modulo one of `3,5,7,11,13`.

### Lengths seven through eleven

Only cells with `H<=3` can remain. The exact cap-four frontier automaton has 288 reachable states and recognizes all such words. Its nonnegative saturation and induced-path pruning are justified by the companion ribbon proof.

Exhaustively enumerating accepted words to these lengths, retaining cells with a path from the last three vertices to the first three, then rejecting any cyclic rotation whose cell has `H>3`, gives these counts:

- Length 7: 62,976 low-count words; 29,952 boundary-active words; 417 final rotation classes.
- Length 8: 199,616; 39,872; 1,091.
- Length 9: 604,160; 29,440; 921.
- Length 10: 1,768,448; 24,576; 675.
- Length 11: 5,147,392; 9,216; 329.

For each multiplier, all 3,433 remaining classes have at least one of the two required determinants nonzero modulo three or five. No class remains unresolved.

The enumeration includes all mask bits in the first three letters, even those ignored by a standalone cell's left padding. Those bits can affect seams between repeated cells and are **not** quotiented away.

Cyclic rotation does not change the characteristic polynomial of a product of square matrices, so one representative suffices for the determinant test. The low-count and boundary-active filters are physical graph tests, not inferences from floating-point spectral estimates.

## 6. Conclusion and controller interpretation

The structural alternatives exclude every block of length at least twelve. The exact finite tests exclude all remaining blocks. This proves both stated all-length nonzero-offset affine obstructions. In particular a fixed-block `3h-2` or `3h+2` induction is impossible at span three.

Together with the all-length binary result, it has one concrete controller consequence. Suppose a construction is meant to implement repeated `h -> 2h+1` by appending words selected by a deterministic finite control state. If that control state eventually repeats with period one or two, then one fixed append block would eventually implement respectively `2h+1` or `4h+3`, contradicting the two theorems. Thus this particular strategy needs an eventual control cycle of length at least three, or a dependence not captured by that fixed finite control cycle. This does not rule out such a richer construction.

## 7. Reproduction

Run from the repository root:

`python session197/ribbon-small-affine/verify.py`

Requirements: Python 3 with assertions enabled, and a C++17 compiler available as `g++`. No third-party Python package, network access, repository build, or external data is needed. The Python driver rebuilds the original frontier and exact low-count automaton, verifies the integer representation, generates the C++ verifier's input, compiles into its own temporary directory, and checks the complete enumeration totals. The C++ checker uses only modular integer arithmetic and literal graph reachability.

`VALIDATION.json` records the replay. The dependencies are the frozen `../ribbon/transfer.py` and `../ribbon-periodic/representation.json`. This packet does not change them or enlarge the target-value census.

Run `python session197/ribbon-small-affine/verify_physical.py` for the independent physical controls. These compare incremental reachability with Floyd closure on all 4,096 canonical order-six inputs and 600 randomized examples, and check 437 literal SCC factorizations against future path counts. `PHYSICAL_VALIDATION.json` records those checks.

## 8. General finite period reduction

The SCC argument has a useful width-independent form.

**Proposition.** Fix an integer width `w>=1` and an integer multiplier `b>=2`. If a fixed-prefix, fixed-block family `H(PW^k)` is eventually nonconstant and obeys `H_(k+1)=bH_k+c` with `c!=0`, then its append block necessarily has length

`|W| < w(b+1)`.

**Proof.** Every strong tournament on `s` vertices has a Hamiltonian cycle and therefore at least `s` Hamiltonian paths, obtained by cutting that cycle at each vertex. For completeness, the cycle exists by the elementary extension argument: start with a directed triangle; insert an outside vertex mixed toward the current cycle; if none is mixed, strongness forces an edge from an outside vertex dominated by the cycle to an outside vertex beating the cycle, and that pair can be inserted together. Repetition spans the graph. Singleton components satisfy the same bound directly.

Consequently, if a cell has `H(G)<=b`, each of its strong components has at most `b` vertices. Any strong component with `s` vertices has index diameter at most `w(s-1)`: take a simple directed path from its largest to its smallest index. It has at most `s-1` arcs; each backward step has length at most `w`, and forward steps cannot increase the total net backward drop. This argument does not require the component to be an interval in the supplied ordering.

If `L=|W|>=w(b+1)`, the distance between every first-`w` vertex and every last-`w` vertex is at least

`L-2w+1 >= w(b-1)+1`.

Their direct edges point forward, and no path can return from the latter to the former: a shared strong component would exceed the diameter bound. The same `I -> M -> J` factorization as in Section 3, with boundary size `w`, therefore makes the repeated sequence exactly multiplicative. This contradicts a nonconstant affine recurrence with nonzero offset.

If instead `H(G)>b`, disjoint-cell path merging gives growth at least `(b+1)^k`, which also contradicts such an affine recurrence. QED.

This bounds the period search for every fixed multiplier and width. It does **not** by itself decide which short blocks admit a suitable physical prefix, and the existence of an eigenvalue pair is only a necessary condition. The sharper low-cell classification in Sections 2–4 reduces the span-three multiplier-three/four checks to length eleven, as used in the completed certificates.
