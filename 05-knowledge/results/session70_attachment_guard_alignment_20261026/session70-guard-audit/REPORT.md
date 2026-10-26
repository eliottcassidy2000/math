# Session 70 independent guard audit: the full Hamilton-path language reconstructs every finite labelled tournament

**Status: PROVED**, by the argument below, using a published arc-traceability theorem whose hypotheses and proof were inspected. This is not a literature-priority claim or a formal Lean verification. **FINITE-EXACT** tests are separate from the all-order argument.

The independent attack began before reading the producer's proposed proof. The decisive mechanism is the single-vertex reachability information retained by the union of Hamilton paths inside each strong component. Arcs absent from all Hamilton paths are not counterexamples to reconstruction.

## 1. Definitions and result

For a finite labelled tournament `T`, let `L(T)` be the **set** of vertex permutations that form directed spanning Hamilton paths, and let `U(T)` be the directed union of their consecutive arcs. Then

`L(U(T)) = L(T)`.

The forward containment follows from `U(T) ⊆ T`; the reverse follows from the definition of `U(T)`. Thus the full language and its directed edge union have exactly the same information, even for arbitrary finite digraphs. The naive assertion `U(T)=T` is false, but this does not imply an information loss.

**Theorem.** If two tournaments on the same finite labelled vertex set have the same set of directed spanning Hamilton paths, they are equal. Equivalently, the map `T ↦ U(T)` is injective on finite labelled tournaments.

## 2. The imported structural fact

Busch, Jacobson and Reid, *On arc-traceable tournaments*, Journal of Graph Theory **53** (2006), 157–166, DOI [10.1002/jgt.20173](https://doi.org/10.1002/jgt.20173), prove the following necessary structure. It is also Theorem 1.1 on printed pp. 4–5 of Busch and Jacobson, *Extremal results on arc-traceable tournaments*, JCMCC **63** (2007), 3–15, [publisher PDF](https://combinatorialpress.com/article/jcmcc/Volume%20063/vol-063-paper%201.pdf).

If `T` is strong and `x→y` occurs in no Hamilton path, there is a vertex `z` such that the strong components of `T−z`, ordered as `C1 → … → Ck`, satisfy:

- `k ≥ 4`;
- `x ∈ C1` and `y ∈ Ck`;
- `C2 → z → C(k−1)`.

The proof also appears as Theorem 4.12, printed pp. 45–47, in Arthur H. Busch's 2005 doctoral thesis, *Arc-traceable tournaments*, [institutional PDF](https://digital.auraria.edu/downloads/ma7yx-pcr25/AA00005831_00001.pdf). The preceding local traceability lemma is Theorem 4.7, printed pp. 42–43. Both proofs were inspected. Only the first two displayed properties are needed in the reconstruction argument below; the third explains the bound of four components.

## 3. Two elementary consequences

### Lemma A: nonseparating vertices expose all their incident arcs

If `T−a` is strong, every arc incident with `a` belongs to a Hamilton path of `T`.

Indeed, take a Hamilton cycle of `T−a` (Camion). For an arc `a→v`, break that cycle immediately before `v` and prepend `a`. For an arc `v→a`, break it immediately after `v` and append `a`.

### Lemma B: interior vertices of a cut decomposition are nonseparating

Let `T` be strong, and write the components of `T−z` as `C1 → … → Ck`. Every vertex `a ∈ Ci` with `1<i<k` has `T−a` strong.

Strongness forces at least one arc `z→c1` into `C1` and at least one arc `ck→z` from `Ck`. Deleting an interior vertex leaves both end components intact. The first component reaches every remaining middle vertex, every such vertex reaches the last component, and the route `Ck → z → C1` closes these routes. Vertices within either end component retain their strong connectivity. Thus every remaining vertex can reach every other.

### Corollary: two traceable two-step replacements for every missing arc

For a nontraceable arc `x→y` in strong `T`, apply the structural fact and choose `a∈C2`, `b∈C(k−1)`. These vertices are distinct. Lemmas A and B show that

`x→a→y` and `x→b→y`

consist entirely of arcs of `U(T)`. These are internally vertex-disjoint, length-two replacements.

## 4. Single-vertex reachability preservation

**Lemma C.** If `T` is strong, then for every vertex `v` and every pair of remaining vertices `s,t`,

`s reaches t in T−v  ⇔  s reaches t in U(T)−v`.

The reverse direction is immediate from inclusion. In a directed `s,t` path of `T−v`, replace each arc not in `U(T)` by one of the two paths from the corollary whose interior vertex is not `v`. The resulting directed walk stays in `U(T)−v`. Erasing repeated-vertex segments gives an `s,t` path.

**Scope is essential:** this is false for a nonstrong tournament. In the transitive tournament `0→1→2` with `0→2`, the core is just `0→1→2`. Deleting vertex 1 destroys core reachability `0→2` but not tournament reachability.

## 5. Reconstruction on a strong block

Suppose strong tournaments `T,T'` have the same core `U`. Present arcs of `U` have the same orientation in both. For a missing unordered pair `{x,y}`, suppose its orientation in `T` is `x→y`. The structural fact gives a vertex `z` such that in `T−z`, `x` reaches `y` but `y` does not reach `x`. Lemma C transfers exactly this relation to `U−z`, and then to `T'−z`. Consequently `T'` cannot have `y→x`, since that edge alone would be a reverse path in `T'−z`. Hence it also has `x→y`. Every pair is recovered.

This gives a decoder: compute reachability in `U−z` for every `z`. For each missing pair, choose any deletion for which exactly one direction is reachable and orient it that way. Such a deletion exists by the structural theorem, and no deletion can prescribe the wrong direction by Lemma C.

## 6. Arbitrary tournaments: ordered SCCs

Let `D1 → … → Dr` be the ordered strong components of `T`. Every Hamilton path is a concatenation of a Hamilton path from each `Di`. Conversely, every such concatenation is legal. Therefore

`U(T)[Di] = U(T[Di])`.

Each nontrivial strong block has a Hamilton cycle; every arc of that cycle belongs to a Hamilton path of the block. Thus its core remains strong. All cross-block arcs point forward, and arbitrary ends/starts of successive strong blocks can be joined, so `U(T)` has the same ordered SCC decomposition as `T`. Singleton blocks cause no difficulty.

Recover that ordered decomposition from `U`, orient every pair in different blocks according to the order, and apply the strong-block decoder within each block. This proves the theorem for every finite order, including the singleton and empty cases under the usual conventions.

Given the core as input, direct repeated reachability gives a polynomial decoder (for example `O(n^4)` with Floyd–Warshall for every deletion). This does not claim that enumerating the full path language is efficient.

## 7. Consequence for the requested attachment carrier

In the source convention, `X_n` consists of tournaments containing the identity path and has `M=(n−1)(n−2)/2` free arc bits. The operation name `R_p` records a particular labelled permutation, and its legal/failure guard is exactly membership of `p` in `L(T)`.

The theorem therefore proves that the complete named legal re-marking guard signature separates **all** elements of `X_n`, for every `n≥1`. Any equivalence respecting every such guard is equality. Accordingly any continuation-stable quotient with these operations and observable failure has exactly `2^M` states at fixed order, and its worst-case fixed-length lossless state cost is `M` bits. One-step observations of **all named permutation guards** already force this; no multi-step argument is needed.

This upgrades the precise all-order open statement in the pinned `strong_atom_continuation_vitali_20261026/README.md`, §6/C5. It does not collapse the distinction between:

- all named permutation operations;
- a small menu of local moves;
- unnamed/fresh uniform re-marking distributions;
- scalar Hamilton-path counts;
- separate endpoint isomorphism types without alignment.

The aligned source/target/forest carrier of source C1 is closed under the specified operations and equivalent to the literal marked state; for this named protocol, no nontrivial quotient of that literal state can retain the guards. Adding attachment operations cannot weaken this lower bound. This is not a minimality theorem for every weaker observation protocol, and in particular does not settle incoming HYP-9416's fresh unnamed multiset protocol.

## 8. Hostile controls and computational evidence

`verify_reconstruction.py` is independent code, importing no producer code. It enumerates Hamilton paths by DFS, builds their directed core, decodes by SCCs and vertex-deleted reachability, and compares the decoded tournament bit-for-bit. It also checks the strong-block reachability lemma directly. `receipt-normal.json` and `receipt-optimized.json` record deterministic completed test universes, outcomes and negative controls; the older timing-bearing JSONL files remain historical scratch receipts. All correctness checks are explicit runtime checks and remain active under Python `-O`.

Reproduce:

```bash
python3 verify_reconstruction.py --max-fixed 7 --samples 500 --receipt receipt-normal.json
python3 -O verify_reconstruction.py --max-fixed 7 --samples 500 --receipt receipt-optimized.json
cmp receipt-normal.json receipt-optimized.json
```

Dense masks use lexicographic unordered pairs and bit 1 for the increasing orientation. `X_n` masks use lexicographic nonconsecutive pairs, fixing `i→i+1`.

A concrete strong hostile to naive arc-union equality is the five-vertex dense mask 10. Its nine Hamilton paths omit arc `4→1`. Deleting vertex 0 gives component order `4,3,2,1`; vertices 3 and 2 give the two traceable replacements `4→3→1` and `4→2→1`. The decoder restores the omitted arc.

Exploratory `closure.py` additionally tested a different unit-propagation decoder: a forbidden permutation with only one unresolved adjacency forces that adjacency's opposite orientation. It closed all dense tournaments through order five and seeded samples at orders six through eight. **Its all-order success is not proved here.** The published-structure decoder above does not depend on that experiment.

An early test implementation sorted SCCs using a key that read the same list during Python's in-place sort; random nonstrong controls caught the bug. The corrected implementation precomputes its keys and was rerun. This was an implementation error, not a failed mathematical case.

## 9. Source verification and publication boundary

Inspected source bytes, kept **scratch-only**, not intended for repository publication:

- 2007 publisher PDF SHA-256: `30f3fd1533d3c28f250060982981bd997ed3d70b98dc1773fafacde6dbe33a3a`.
- 2005 institutional thesis PDF SHA-256: `76c57cfb33b8ad32e7e830e811bb3a3a38d4a2eeea72416d03b45412f51cdf19`.

The web tool could not open the thesis's application/octet-stream response; an ordinary direct fetch from the same official institutional URL succeeded. The theorem in the 2007 scanned PDF was checked visually, including `k≥4`; OCR alone initially misread that symbol. Full source PDFs, images and OCR are not publication artifacts. Intended audit artifacts are this report, independent verification code and test receipts.

## 10. Cross-audit: elementary relocation guards versus their traces

The producer's one-step counterexample was checked independently. At order four, fixed-path masks `0` and `5`, with outgoing rows `(2,4,9,3)` and `(6,12,8,1)`, have identical guards for every named single-vertex relocation in either direction. Exhaustive signatures at orders one through three are injective, so this first failure order is minimal.

The two-step operation word

`p=(1,2,3,0)`, followed in the new position frame by `q=(0,2,3,1)`,

is legal at its first step for both states. Its second step fails for mask 0 and succeeds for mask 5. The accumulated original-label word is `p∘q=(1,3,0,2)`. This is a shortest separating elementary trace because no one-step guard separates the pair.

The pinned all-order dependency

`05-knowledge/results/hamilton_path_relocation_connectivity_20261020/THEOREM.txt`

was independently inspected in full. It proves connectivity of the graph of **labelled Hamilton words** under legal single-vertex relocation. The delicate common-deletion, shared-gap edge-lift, and recurrent-front arguments are valid. Its normalized-isomorphism-class corollary is not substituted for the stronger labelled-word statement needed here.

Consequently, the complete **trace** language of named leftward and rightward relocations also reconstructs the marked tournament. Every Hamilton word is reachable from the identity word by a legal relocation trace. Given the operation names, the accumulated permutation is determined independently of the tournament: update the current original-label word `w` by `w_new=w_old∘p`. Equality of successful trace languages therefore implies equality of full Hamilton-word sets, to which the new reconstruction theorem applies.

This corollary requires all finite traces, not merely the one-step guard vector. The connectivity proof uses both move directions and arbitrary spans; it does not by itself prove a rightward-only or bounded-span trace statement. No polynomial trace-length or efficient exploration bound is asserted.

## 11. Final validation controls

The deterministic normal and optimized receipts cover all **33,868** states in `X_1` through `X_7`, plus 503 explicitly seeded/structured cases at each of orders 8, 9 and 10. Negative controls require rejection of (i) applying the strong reachability lemma globally to transitive `TT3`, (ii) treating the five-vertex mask-10 core as the whole tournament without reconstruction, and (iii) a corrupted mask-10 core with the genuinely traceable arc `4→3` deleted. These controls use explicit exceptions, not removable `assert` statements.

The publication manifest lists only this own-written report, the independent verifier, deterministic normal/optimized receipts and source metadata. It excludes downloaded source bytes, source OCR, page images, exploratory programs and timing-bearing historical logs.
