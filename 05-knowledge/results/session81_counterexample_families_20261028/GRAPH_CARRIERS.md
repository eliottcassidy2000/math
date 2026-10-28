# Sparse graph carriers: the recovered trijection, a capacity obstruction, and a faithful repair

**Session 81; status: elementary PROVED statements with proofs below; FINITE-EXACT census n≤30 and codec tests n≤5. No novelty claim.** The historical results below are inherited, not newly proved literature results.

## Findings first

1. A tournament type cannot in general be encoded by a pair or triple of **independent graph types on the same n vertices with exactly n−1 edges each**. Already n=3 fails: two tournament types, only one such graph type. Further exact failures begin at n=9 for pairs and n=14 for triples in the n=4..30 census.
2. This is a capacity issue, not a mysterious failure of a proposed algorithm. For any fixed number k of sparse layers, even retaining all labels and cross-layer alignment eventually supplies too few outputs: O(kn log n) bits versus Θ(n²) bits.
3. Allow vertex expansion and the result reverses completely: **ONE sign-even simple graph with N−1 edges suffices**. The recovered Session66 port codec, padded by isolates, has

   `N = binom(n,2)+8n+1`, `E=N−1`.

   It is explicit, polynomial-time, intrinsically decodable, sign-even, and faithful on tournament isomorphism classes. Dense cyclic structure is paid for by disconnected components; `E=V−1` does not mean tree.
4. The earlier trijection concerns **marked** tournaments, not an unmarked sparse-graph class encoding. The exact equal free-coordinate count `binom(n,2)−(n−1)=binom(n−1,2)` is real; the mark's symmetry and orientations are indispensable sidecars.
5. “Even graphs” in tournament equinumerosity means **automorphism-sign-even**, not even degrees, not an even number of edges, and not E=V−1. These are three different restrictions.

## 1. Inheritance and type ledger

Closest mechanisms: `trijection_four_nodes_three_edges_20261020.md`, Theorem 1; Session66 `REPORT.md`, §§2–6, graph port codecs; THM-4531/4721 and Session72 `EVEN_GRAPH_OPERATION_ENVELOPE.md` for the signed harmonic explanation of equinumerosity. Canonical hostile: triangle plus an isolate has V=4,E=3 but is not a tree. Corrected near miss: equal raw free-bit counts do not identify orbit sets; separate graph isomorphisms forget a shared frame. Least-used sidecar: the actual marker stabilizer, or alternatively intrinsically marked center/tail/head ports.

Live board: (i) marked free-coordinate cubes, (ii) automorphism-sign-even graph modes, (iii) sparse graph capacity, (iv) connected core plus isolated vertices, (v) alignment versus independent quotients. The core-padding construction connects (iii) and (iv), and its intrinsically oriented port version also lands in (ii); the orbit-stabilizer law connects (i), (ii), and (v), but does not turn their maps into the same bijection.

### The recovered trijection, precisely

Let t(n) be the number of tournament isomorphism classes. For a tournament together with one specified coherent mark on all n vertices:

| mark | orientations | marked-type count |
|---|---|---:|
| Hamiltonian path P_n | directed along path | `2^binom(n−1,2)` |
| C_(n−1) plus one isolate | directed cycle | `(1/(n−1)) sum_(d|n−1, d odd) phi(d) 2^(binom(n−1,2)/d)` |
| star | all outward OR all inward | `2t(n−1)` |

The respective stabilizers of a fixed oriented mark are 1, C_(n−1), S_(n−1). The star has two distinct orientation orbits (source and sink). The line forgetful fibre over [T] has size `H(T)/|Aut(T)|`, not one: automorphisms act freely on Hamiltonian path words. At n=4 the three counts are 8,4,4; at n=5 they are 64,16,8. Therefore the numerical sum 8=4+4 is not an all-order law.

At n=4 the three **unoriented** sparse graph types are P4, K_(1,3), K3+K1 (12,4,4 labelled copies). Complement swaps the last two and preserves the first type. This half-edge coincidence is special among n≥2: `binom(n,2)=2(n−1)` iff n=4 (n=1 is a degenerate zero-edge exception). The historical note restricts n≥4; do not export its uniqueness claim to n=1 without that exception.

### Even-graph equinumerosity is a separate target

For any orientation of a simple graph G, an automorphism has sign (−1) to the number of edges it reverses. G is sign-even if every automorphism has positive sign; this condition is orientation-independent. Royle–Praeger–Glasby–Freedman–Devillers prove t(n) equals the number of sign-even graph types on n vertices. **CITED**, primary abstract and definition checked 2026-10-28: [Tournaments and Even Graphs are Equinumerous](https://arxiv.org/abs/2204.01947), J. Algebraic Combin. 57 (2023), 515–524.

There is no E=n−1 restriction in that theorem. For example P3 is sign-even but not degree-even; C3 is degree-even but sign-odd; the empty graph is sign-even. The inherited naturality obstructions concern equivariant maps, whereas an arbitrary type bijection only requires equal finite cardinalities. No new symmetry-preserving claim is made for the padding map below.

## 2. Fixed-carrier capacity theorem (PROVED)

Let S_N be simple graph types on exactly N vertices with N−1 edges, and s(N)=|S_N|. Let k≥1 be fixed. Ordered independent k-tuples supply exactly s(N)^k possible outputs. Joint tuples up to one common relabelling supply at most the number of labelled tuples,

`B(N)^k`, where `B(N)=binom(binom(N,2),N−1)`.

There are at least `2^binom(n,2)/n!` tournament types on n vertices, since every labelled orbit has size at most n!.

**Theorem.** If

`2^binom(n,2) > n! B(N)^k`,

no injection of tournament types on n vertices into k sparse graph layers on N vertices exists, even when the target is fully labelled. In particular, with N=n any fixed k fails for all sufficiently large n. With N allowed to expand and k fixed, any successful sequence of encodings must have

`N log(N+1) = Ω(n²)`, hence `N = Ω(n²/log(n+1))`.

**Proof.** Pigeonhole principle gives the displayed inequality. For N≥2, `B(N) ≤ (eN/2)^(N−1)` by `binom(a,b)≤(ea/b)^b`; thus log2 B(N)=O(N log N). On the source side `log2 t(n)≥binom(n,2)−log2(n!)=n²/2−O(n log n)`. Comparing proves the first necessary asymptotic. If N≥n² the second is automatic; otherwise log(N+1)=O(log(n+1)), proving it. QED.

No unproved asymptotic for unlabelled sparse graphs is needed. Independent tuples can be much smaller than joint tuples, so an independent-tuple obstruction must not be silently promoted to a joint-tuple obstruction. The elementary labelled upper bound already supplies explicit joint/labelled obstructions at n=28 for k=2 and n=43 for k=3; these are example orders, not a claim that every larger order has been separately checked by this finite computation.

### Exact independent-tuple census (FINITE-EXACT)

| n | s(n) | t(n) | t(n)≤s(n)²? | t(n)≤s(n)³? |
|---:|---:|---:|:---:|:---:|
| 1 | 1 | 1 | yes | yes |
| 2 | 1 | 1 | yes | yes |
| 3 | 1 | 2 | **no** | **no** |
| 4 | 3 | 4 | yes | yes |
| 5 | 6 | 12 | yes | yes |
| 6 | 15 | 56 | yes | yes |
| 7 | 41 | 456 | yes | yes |
| 8 | 115 | 6880 | yes | yes |
| 9 | 345 | 191536 | **no** | yes |
| 13 | 46415 | 48542114686912 | **no** | yes |
| 14 | 175652 | 28401423719122304 | **no** | **no** |

Running `python3 check_sparse_capacity.py > counts.json` generates all rows n≤30. “Yes” means an arbitrary injection exists by finite ranking, not that a natural, constructive, operation-preserving, or efficient map has been supplied. “No” forbids every injection in the specified universe. For n=3 all sparse graphs are P3, so transitive and cyclic tournaments give the smallest explicit pair/triple hostile. Taking an unordered pair/triple only reduces capacity.

**Reproduction formula.** For a partition λ of n, let z_λ be its symmetric-group centralizer order, and let ℓ range over the lengths of its orbits on unordered pairs. Then

`s(n)=sum_(λ⊢n) [x^(n−1)] product_ℓ(1+x^ℓ)/z_λ`,

`t(n)=sum_(λ⊢n; all parts odd) 2^(number of pair orbits)/z_λ`.

Inside a vertex cycle of length a, there are floor((a−1)/2) pair orbits of length a and, if a is even, one of length a/2. Between cycles a,b there are gcd(a,b) orbits of length lcm(a,b). An even vertex cycle has a power swapping endpoints of a pair, so it fixes no tournament; for odd cycles each pair orbit admits exactly two consistent orientations. These observations prove the formulas independently of stored count sequences.

## 3. Disconnected deficit completion (PROVED)

For a finite simple graph G with v vertices, e edges and c components, its cycle-space dimension is

`β(G)=e−v+c`.

Consequently `e=v−1` iff `β=c−1`. Connected graphs in this family are trees; disconnected graphs can carry arbitrary cycles, provided extra components pay the deficit.

**Padding lemma.** Suppose G has no isolated vertices and e≥v−1. Put q=e−v+1 and

`P(G)=G disjoint-union q K1`.

Then P(G) has N=e+1 vertices and N−1 edges. Deleting every isolated vertex recovers G. Hence P is faithful on isomorphism classes of such G. Also

`Aut(P(G)) ≅ Aut(G) × S_q`.

**Proof.** The size assertion is arithmetic. Isolated vertices are intrinsic and G has none, so isomorphisms restrict to and extend from G; the automorphism factors act on disjoint intrinsic sets. QED.

Moreover P(G) is sign-even iff G is sign-even: every automorphism splits as above, isolated permutations reverse no edges, and the remaining edge sign is exactly that on G.

The no-isolates condition is necessary for this stated inverse. A pair of graphs differing only in the number of pre-existing isolated vertices, with the same e and with nonnegative q, can pad to the same graph. Removing isolates is not a general inverse without that hypothesis or a sidecar.

**Procedural hostile family to “E=V−1 implies tree.”** Take any connected non-tree G and append β(G) isolates. For instance `C_r disjoint-union K1`, r≥3, gives V=r+1,E=r and one cycle; `K_r disjoint-union (binom(r,2)−r+1)K1`, r≥3, permits unbounded clique number, chromatic number, treewidth and cycle rank. Thus E=V−1 plus disconnectedness imposes no bound on those core invariants. This is a counterexample generator, not just an isolated witness.

## 4. One sign-even sparse graph faithfully encodes every tournament (PROVED)

Use the inherited Session66 graph Γ(T). For each tournament vertex v create center c_v, tail t_v, head h_v; join c_v to both ports; attach one leaf to c_v, two leaves to t_v, three leaves to h_v. For each arc u→v add undirected edge {t_u,h_v}. This produces a simple graph with

`v(Γ)=9n`, `e(Γ)=8n+binom(n,2)`.

It has no isolated vertices; every nonleaf has degree at least three. The number of leaf neighbors identifies center, tail, head intrinsically as 1,2,3. Each center identifies its two ports, and port edges recover every tournament arc. Therefore Γ(T)≅Γ(T′) iff T≅T′, with construction/decoding polynomial-time. This is the recovered result, not a new codec claim.

Apply the padding lemma with

`q=binom(n,2)−n+1=(n−1)(n−2)/2 ≥0` for all n≥1.

The resulting Ψ(T) is one simple graph with

`N=binom(n,2)+8n+1`, `E=N−1`,

and Ψ(T)≅Ψ(T′) iff T≅T′. Delete isolates, recognize the leaf-neighbor roles, and read the tournament. Graphs in the image have one connected nontrivial component plus q isolates: Γ is connected because every pair of tournament vertices is joined by one port edge between their connected local gadgets. Its cycle rank is exactly q, which pays for the q extra components.

Its automorphism group is `(S2×S3)^n ⋊ Aut(T)`, times S_q from the padding. This is **not** an automorphism-group-preserving code; it is faithful on classes. The target can in fact be restricted to sign-even sparse graphs (proof immediately below), but not to degree-even graphs: the pendant leaves have degree one. Two or three coordinates can trivially be supplied by repeating Ψ or adding a fixed sparse graph of size N, but the second and third coordinate are unnecessary for fidelity.

**Sign-even strengthening (PROVED).** Orient every pendant edge from its core endpoint to its leaf, both center–port edges from center to port, and every tail–head edge from tail to head. This is an orientation of Γ(T) fixed by **every** graph automorphism: leaf/core status and all three core roles are intrinsic. Thus every automorphism reverses zero edges, stronger than merely an even number. Γ(T), and hence Ψ(T), is sign-even. This does not orient Γ as a tournament: only existing edges receive directions. It also does not make Aut Γ odd-order; the leaf permutations preserve this sparse orientation while contributing an even-order kernel.

This is a genuine common target with the even-graph thread, but remains an **expanded-order injection**, not the same-order counting bijection of Royle et al. The sign-even subset at the expanded size N is large enough to contain the constructed image; no claim of surjectivity onto it is made.

The general expansion lower bound Ω(n²/log n) and this explicit O(n²) construction leave a logarithmic gap. No optimality claim or efficient near-optimal construction is made.

## 5. General mark law and what actually unifies the themes

Fix any **oriented** simple carrier F on n vertices with r edges, and let A=Aut(F) be its orientation-preserving automorphism group. Tournaments extending F form an affine cube of dimension binom(n,2)−r. Marked tournament isomorphism classes with a mark isomorphic to F are precisely the A-orbits of that cube. This is simply transport to a fixed mark followed by the stabilizer action; it does not assume F is connected.

For a∈A, let o(a) be its number of orbits on the unmarked pairs. The number of fixed completions is zero if some orbit returns with its endpoints swapped, and otherwise it is 2^o(a). Burnside therefore gives the exact orbit count. This general rule recovers line, cycle-plus-isolate and oriented star. The number r=n−1 identifies the **dimension** with binom(n−1,2), but an actual correspondence to tournaments on n−1 vertices requires a compatible pair-coordinate action. Equal dimensions alone do not supply it.

The unifying law is an information ledger:

- marking removes free bits but retains a frame and its stabilizer;
- quotienting independent layers discards alignment;
- sign-even graph equinumerosity counts stabilizer-compatible Fourier modes;
- sparse carrier bounds restrict output information only if vertex count is controlled;
- isolate padding changes the carrier size while retaining all core information.

These are rigorous compatible mechanisms, not a single unsupported grand bijection. The new repair says exactly where to pay the missing capacity: in carrier expansion and disconnected components.

## 6. Verification and remaining scope

Run from this directory:

```sh
python3 -B check_sparse_capacity.py > counts.json
python3 -B -O check_sparse_capacity.py > counts-optimized.json
cmp counts.json counts-optimized.json
python3 -B check_codec.py > codec.json
python3 -B -O check_codec.py > codec-optimized.json
cmp codec.json codec-optimized.json
```

Both scripts use only Python's standard library, exact integers/rationals, and explicit raises rather than removable asserts. Capacity formula is checked against direct graph canonicalization for n≤5 in the separate codec script. Codec tests enumerate every labelled tournament n≤5 and check exact decoding under three explicit relabellings per graph, all sizes, and isolate counts. These tests support, but do not replace, the all-order structural proofs.

Open useful next question: obtain a practical, intrinsically decodable O(n²/log n)-vertex E=V−1 codec without canonicalizing the tournament first, or prove a stronger lower bound under an explicit naturality/local-gadget requirement. The present result does not solve that problem. It does settle the ambiguity between same-n and expanded-n pair/triple proposals.

## Independent audit

A separate implementation formed the actual unordered-pair permutations for each partition, reproducing the displayed order-9 and order-14 capacity failures. Both public checkers passed normal and optimized runs with identical output. For exact preservation of the source automorphism group, see [Sparse symmetry](SPARSE_SYMMETRY.md); the present smaller codec preserves classes and sign-evenness, not the full group.
