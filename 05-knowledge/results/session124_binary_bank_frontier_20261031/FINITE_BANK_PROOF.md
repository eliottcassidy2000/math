# Independent finite binary-bank audit

Session 124, 2026-10-31.

**Verdict: the balanced depth-two obstruction through parent order nine passes an independent exhaustive implementation audit. The explicit three-level phase bank also passes.** Neither statement is an all-order obstruction, an infinite compiler, or a proof of the unrestricted Hamiltonian-count spectrum.

## Exact finite statement and quantifiers

For a tournament `T`, let `H(T)` count directed Hamiltonian paths as ordered vertex sequences. For every subset `S` of its vertices, `T+x(S)` means the new vertex has exactly `S` as its out-neighbors. Put

```
B(U) iff
  (there exists A subset V(U) with H(U+x(A)) = 2H(U)-1)
  and
  (there exists B subset V(U) with H(U+x(B)) = 2H(U)+1).
```

The independently checked finite theorem is:

> For every tournament `T` with `1 <= |V(T)| <= 9` and `H(T) >= 9`, it is not the case that both
> (i) some cut `S_-` has `H(T+x(S_-))=2H(T)-1` and `B(T+x(S_-))`, and
> (ii) some cut `S_+` has `H(T+x(S_+))=2H(T)+1` and `B(T+x(S_+))`.

Every existential ranges over **all** cuts of the actual tournament at that node. The minus and plus choices are independent; the two cuts at each child are independent; every eligible first child is considered, not merely the first candidate. The independent checker deliberately examines all eligible children even when it has already found a successful one-sided continuation. A parent lacking either first digit cannot be a depth-two bank and is discarded only for that necessary reason.

The result remains true under arbitrary relabeling because both Hamiltonian path count and all physical one-vertex extensions are transported by isomorphisms. No hypothesis of strong connectivity, common physical gadget, common cut, or fixed active vertices is imposed.

## Independent finite universe: generation, uniqueness and coverage

The submitted program minimizes an edge code over permutations within sorted outdegree classes. `independent_census.cpp` instead uses equitable partition refinement and individualization:

1. Start with the single cell containing all vertices.
2. Split every cell by each vertex's vector of outgoing-neighbor counts into the current ordered cells; order the new cells lexicographically by those vectors. Repeat until stable.
3. If a cell is nonsingleton, take the first such cell and branch on each of its vertices, placing the chosen vertex in a singleton immediately before the rest of that cell. Refine again.
4. At a discrete partition, encode the labeled tournament. Minimize over the leaves, and count the leaves attaining the minimum.

Every operation commutes with relabeling. Thus isomorphic tournaments have equal resulting codes. Conversely, equal codes explicitly identify isomorphic labeled adjacency matrices. This is an exact canonical form, not a hash or heuristic signature.

Different branches give different final vertex orders: the first selected vertex stays before the remainder of its cell, since refinement never merges or moves vertices across cell boundaries. Automorphisms act bijectively on the refinement search tree. For any fixed leaf code, its leaf labelings form a free transitive orbit under the automorphism group: two labelings giving the same code differ by a unique automorphism, and every automorphism transports a leaf to a leaf. Therefore the multiplicity of the minimum is exactly `|Aut(T)|`.

Starting with the unique one-vertex tournament, the checker extends every representative by every cut and canonicalizes independently. Induction covers all tournaments: deleting the last vertex of any labeled tournament gives an already represented parent, and relabeling that parent transports the deleted vertex's out-neighbor cut. Sets of exact canonical codes guarantee uniqueness.

There is also a separate, exact exhaustion certificate at each order:

```
sum over generated representatives T of n! / |Aut(T)| = 2^(n(n-1)/2).
```

The summands count disjoint labeled isomorphism orbits. Their sum exhausts all labeled orientations of the complete graph. Thus the class counts do not rely on merely agreeing with a published sequence. The original audit additionally mapped every submitted representative to this independently generated universe, rejected duplicate isomorphism classes, and required exact equality of sets (`PRODUCER_BIJECTION_RECEIPT.json`). This optional comparison is not an input dependency: the publication and actual relocated runs generate the full universe afresh without producer representative files.

A further implementation check literally enumerates all vertex permutations for every submitted representative of orders one through seven and independently counts automorphisms; its complete histograms match the refinement multiplicities. This check is in `check_automorphisms.py`.

## Independent all-cut response and direct-DP cross-checks

The independent C++ code computes start/end counts by **removing** a vertex from a subset; it does not use the submitted forward-propagation implementation. Write `E[A,u]` for paths ending at `u` on `A`, and `R[A,v]` for paths starting at `v`. For distinct vertices,

```
W_uv = sum over nonempty proper A subset V of E[A,u] R[V-A,v].
q_uv = W_uv + W_vu.
```

Deleting the new vertex from a path gives the quadratic polynomial in the submitted note. For an original arc `u->v`, concatenation identifies `W_uv` with the Hamiltonian paths using that arc. Summing endpoint divergence therefore turns the polynomial into

```
H(T+x(S)) = H(T) + sum_{u in S, v outside S, u->v} q_uv.
```

This identity is also proved and independently audited in the companion [canonical cut-response theorem](../session124_tournament_cut_response_20261031/README.md). The present implementation evaluates each cut by a direct **nonnegative outgoing-arc sum**. It does not use the submitted linear-coefficient/subtraction recurrence. The polynomial proof and this orientation convention agree; there is no missing factor of two (`q_uv=2b_uv`).

The independent implementation compares these responses with Hamiltonian subset DP on the **actual explicitly built child tournaments** for:

* every cut of every representative parent through order seven;
* every cut of the first parent and every 257th independently ordered representative at orders eight and nine;
* **every cut of every qualifying first-level child of every parent with both first digits**, through parent order nine.

The last item means that all potentially decisive grandchild counts are counted by direct DP, not accepted merely from the cut formula. Exact numbers are in `census.normal.json` and `census.optimized.json`.

`check_polynomial.py` separately imports the actual submitted `cut_polynomial.py` (a byte-preserving snapshot is included). It compares its full output to the independent program's responses after that program has directly counted every actual child. This covers all labeled parents through order four and 22 deterministic test parents at each order five through ten: **207 parents and 45,450 all-cut comparisons**, reaching child order eleven. This is an implementation cross-check, not the proof of the general identity or the exhaustion argument.

## Results

| Parent order | Isomorphism classes | H >= 9 | Both first digits | Minus continuation only | Plus continuation only | Neither continuation | Full depth two |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3 | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| 4 | 4 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5 | 12 | 6 | 2 | 0 | 2 | 0 | 0 |
| 6 | 56 | 48 | 14 | 2 | 4 | 8 | 0 |
| 7 | 456 | 446 | 45 | 6 | 6 | 33 | 0 |
| 8 | 6,880 | 6,868 | 128 | 12 | 8 | 108 | 0 |
| 9 | 191,536 | 191,522 | 301 | 20 | 10 | 271 | 0 |

At order nine, the orbit mass is **68,719,476,736 = 2^36**. The 840 qualifying children give 860,160 grandchild cuts, all directly checked. The order-nine run makes 1,242,112 direct child-DP comparisons in total, including the parent controls. Across all nine orders, the run checks 1,472 qualifying children and all their 1,113,856 grandchild cuts, making 1,565,142 direct child-DP comparisons overall. The independent counts exactly reproduce the submitted standard balanced census.

## No overflow through order eleven

All Hamiltonian counts on `m` vertices are bounded by `m!`, simply because each path is a permutation. In a fixed endpoint DP entry the stronger bound is `(m-1)!`. For a response calculation the parent has `m <= 10`; its returned child has order at most eleven.

For each ordered pair, `W_uv <= (m-1)!`: concatenate the two paths into a permutation with the consecutive block `uv`. This is injective, since the position of `uv` recovers the cut. Also

```
sum_{u != v} W_uv <= (m-1)m!,
```

by injecting an arbitrary ordered two-path cover into a permutation plus one of its `m-1` internal break positions. Thus at `m=10`, every `W_uv <= 362,880` and their total is at most **32,659,200**. Each multiplication in the convolution is bounded by `(a-1)!(m-a-1)! <= (m-2)!`, and all accumulation terms are nonnegative. In the independent implementation every response is a sum of nonnegative terms and is at most **11! = 39,916,800**.

The submitted implementation's potentially signed linear coefficient satisfies `|L_v| <= m! + (m-1)m! = m*m!`. Its recurrence first adds `L_v` to a valid previous response and then subtracts nonnegative pair terms. The intermediate magnitude is bounded, conservatively, by `11! + 10*10! = 76,204,800`; the final response is nonnegative, so partial subtraction cannot overshoot it. These bounds are far below signed 64-bit capacity. Target doubling is likewise safe. The independent program checks the factorial, pair, total-cover and response bounds with explicit exceptions.

Canonical codes for parent representatives use at most **36 bits**. The independent orbit mass uses 64-bit integers (32-bit would be insufficient for `2^36`); factorials and automorphism multiplicities also fit. Graph masks use at most eleven bits. The audited search never shifts a 64-bit word by 64. No floating-point arithmetic or probabilistic hash is used in the mathematics. A supplementary undefined-behavior/bounds-sanitized build also passes the complete census through order seven and three full order-ten response tests (all child cuts, reaching order eleven), with empty diagnostics and output identical to the normal build; see `SANITIZER_RECEIPT.json`. The sanitizer scope is not the full order-nine census.

## Positive phase certificate: stronger independent verification

`phase_bank.json` is the submitted physical tree, whose seed adjacency is `(8,1,3,6,15)`, edgecode 4 in the stated pair order. `check_phase.py` constructs every child adjacency entry directly and counts paths by **literal vertex permutations**, importing neither DP nor cut-polynomial code. It checks all 15 nodes and 14 edges, yielding

```
5
13,15
23,25,27,29
73,75,77,79,81,83,85,87
```

with common offsets `(3,-3,27)` and the `+2` companion at each depth. Five hostile mutations—wrong cut, wrong offset, missing branch, invalid tournament and nonempty terminal proof—are rejected. Python normal and `-O` receipts agree.

This proves the concrete three-level tree on this parent only. The arithmetic form of its leaves is `8H+33+2j`, `0 <= j <= 7`, at this seed. It does not prove uniform realizability on arbitrary parents or indefinite iteration. In particular it does not conflict with the standard offsets `(-1,-1)` obstruction.

## Promotion and remaining boundaries

The **standard balanced depth-two, parent n<=9, H>=9** result may be promoted from COMPUTED to independently VERIFIED finite exhaustive census. The positive three-level phase tree is independently verified even without trusting subset DP.

This audit does **not** promote the separate constant-offset depth-three census, period-two offset census, fourth-level stopping searches, strong-component counts, random probes, or moving-active-vertex implementation. Those are separate claims requiring their own scoped audit. No conclusion here quantifies over arbitrary parent order, unbounded offset sequences, or infinite depth. Nothing here establishes that `7,21` are the only forbidden tournament Hamiltonian-path counts.

## Publication replay boundary

`AUDIT_RECEIPT.json` is the new, standalone publication replay without external representatives. `PRODUCER_BIJECTION_RECEIPT.json` preserves the prior optional scratch-input comparison as historical provenance. `PORTABILITY.json` compares the full publication run with an actually relocated package run under Python `-O`; both include normal and optimized C++ censuses and all normal/optimized Python checks. Neither replay reads the original scratch representative directory.
