# Independent audit: exact intercyclic spectrum and dyadic obstruction

Status: ACCEPT. All mathematical claims in the final `intercyclic_rays.md` pass proof-level review, and both producer and independent normal/optimized controls agree. No numerical fit or finite prefix is used to close an all-order implication.

Audited final proof SHA-256: `54dba01e199b76c139ba6b58a2a2451a921f323b60e23470a38f98a5b24f4e4b`.

This supplements `AUDIT.md`, which covers the initial one-apex and bounded-run packets. Together the reports cover all three frozen packets. The global odd-spectrum conjecture remains OPEN, and no novelty claim is certified by this review.

## 1. Exact classification and its quantifiers

The previously verified Bang-Jensen--Bessy Theorem 1.8 reduces every intercyclic tournament to a two-feedback presentation or a one-apex presentation. With feedback arc x->y, the one-feedback triangles correspond exactly to inversions in the two binary coordinates. A triangle containing both feedback vertices meets every triangle, since the ordinary chain is transitive. Thus absence of two disjoint triangles is equivalent to cross-intersection of the two inversion-edge families. Empty inversion families imply an actual one-apex presentation, not merely a scalar equality.

Stripping leading 0 and trailing 3 removes actual sources and sinks. Neither can participate in an inversion edge, so the nonempty-family hypothesis and cross-intersection property survive.

Every endpoint case in the reduced-word classification is exhaustive and valid:

- For endpoints 1,2 (or 2,1), an outer inversion in one coordinate forces the other family's edge to meet an endpoint that cannot belong to such an edge. This is impossible.
- For endpoints 1,1, every x-zero and every y-one must occur at the same unique interior position. It is therefore one letter 2, with all other letters 1. The bit-swapped case is identical.
- For endpoints 1,0, every y-inversion ends at the last vertex. If there are at least two earlier y-ones, every earlier x-zero would conflict with one of them, so the prefix is 1^a3^b. If there is one y-one, it is penultimate; its value gives either that family or 1^a20. The bit swap gives the other two families.
- Coordinate complement plus word reversal handles endpoints 3,1 and 3,2. This is a symmetry of the inversion condition; the later tournament-duality map is correctly treated separately and additionally swaps the feedback labels.
- For endpoints 3,0, neither coordinate has an internal inversion. An interior 0 forces every interior letter to be 0, and an interior 3 forces every interior letter to be 3. Among letters 1 and 2, repetition creates disjoint endpoint inversions, while mixing creates an internal inversion disjoint from the other outer inversion. Only one such letter or none remains.

Conversely every listed template cross-intersects. The duplication of 30 is harmless. The classification is an all-length proof, not an extrapolated language model.

The stated deletions for the one-apex subfamilies genuinely remove every triangle. For the remaining families, duality is exactly reverse(word) with 0<->3 and 1,2 fixed after relabelling the feedback pair. It preserves ordinary H and gives the stated template pairings.

## 2. Direct cycle counts and the two-ray equality

The odd-cycle-collection formula counts collections of vertex-disjoint directed odd cycles with weight 2 per nontrivial cycle. Since the tournament is intercyclic, only the empty collection and singleton cycles occur. Directed cycles must still be counted modulo rotation, not merely by their vertex sets. The proof and verifier respect this convention.

I independently checked the source: Irving--Omar, arXiv:2412.10572v3, Corollary 20 on PDF page 14, attributes the formula to Grinberg--Stanley Theorem 1.39 and Lemma 6.5. The source's complement convention is harmless here: tournament reversal bijects directed Hamilton paths and directed cycles. Primary source: https://arxiv.org/pdf/2412.10572 . Direct PDF retrieval succeeded. The producer added this citation before the final freeze.

The single-feedback cycle columns are exactly the weighted-inversion sums. For both-feedback cycles, each cycle either uses x->y or splits uniquely as x,A,y,B,x with nonempty increasing ordinary blocks A,B. Rechecking all parity classes gives:

- For 1^a20, the x->y class has one triangle. The other classes contribute 2^(a-1), 2^(a-1)-1, and 2^(a-1), for a total of 3*2^(a-1).
- For 2^a10, both the arc-using and arc-avoiding classes contribute 3*2^(a-1)-1, for a total of 3*2^a-2.
- For 2^a1 2^b, the three contributions are 2^(a+b-1), 2^(a+b-1)-2^(a-1)-2^(b-1)+1, and 2^(a+b-1)-1.
- For 2^a3^b0, they are 2^(a-1), 2^(a+b-1)-1, and 2^(a+b-1)-2^(a-1)-2^(b-1)+1.

All subsets and parity cases are accounted for, including a=b=1. This yields exactly the table in the proof. The short words 310 and 320 have six and seven directed odd cycles in total, respectively, hence H=13 and H=15.

The reverse inclusion is not missing: 1^(t-1)20 realizes 3*2^t+1 for every t>=2, and 2^t10 realizes 9*2^t-3 for every t>=1. Their displayed spanning cycles are valid. Therefore the exact all-order identity is accepted:

S_intercyclic = S_oneapex union {3*2^t+1: t>=2} union {9*2^t-3: t>=1}.

The unions need not be disjoint. Only the values outside S_oneapex are asserted to have O(log X) counting function, which follows immediately and correctly.

## 3. H=29, exact triangle packing, and the divisibility boundary

The one-apex argument for excluding H=29 uses all orders correctly. Orders at most six have maximum 23; strong cores of order at least eight have minimum at least 45. For chain length six, fixed-content minima are 33,27,25,27,33, and every nonminimal word lies at least six above its own minimum. Thus 29 is absent. Neither dyadic ray equals 29.

The six-vertex witness has H=29 and two disjoint directed triangles. Every tournament of total prime count 29 has just one nontrivial SCC, by oddness, H=1 iff transitive, and SCC multiplicativity. The Busch bound puts that component at order at most seven, hence its triangle packing is at most two. The intercyclic exclusion supplies at least two, proving exact packing number two, including with arbitrary singleton SCC padding.

The strong one-apex H=87 witness is also correct. Its count is divisible by 29 but all cycles meet its apex. The corollary depends on the whole count being 29, not merely divisibility. The verifier and independent replay exhibit this boundary explicitly.

## 4. Dyadic counting constants

For X=2^L, a count at most X has reduced one-apex chain length m<=L+1. The full candidate bank contains X-1 words. Its complement-reversal involution has exactly 2^floor((L+1)/2)-1 fixed words.

At m=L+1, the four endpoint inversions in every word beginning 11 and ending 00 contribute 9*2^(m-4), so H>=1+9*2^(m-4)>X. This excluded family has X/8 words and is involution-stable. Counting the remaining orbits using the original fixed-point count, then adding the transitive value, gives exactly

|S_oneapex intersect [1,X]| <= 7X/16 + 2^(floor((L+1)/2)-1).

For L>=5 the first ray has L-3 terms and the second has L-4 terms below X. Adding those possible values gives the stated intercyclic bound and the resulting deficit of at least

X/16 - 2^(floor((L+1)/2)-1) - 2L + 7.

Relative to the X/2 positive odd integers, this is 1/8-o(1) along the specified dyadic cutoffs. The claim is correctly about the intercyclic spectrum, not the global tournament spectrum. It does not show the missing values are globally unattainable, and this audit does not silently promote the dyadic statement to a claimed asymptotic density at every cutoff.

## 5. Independent computational evidence

The producer verifier reran successfully under normal Python and `-O`, with outputs byte-identical to both frozen receipts.

`independent_intercyclic.py` imports only the auditor's clean-room utilities. It uses a regular-expression recognizer rather than the producer's finite template bank, actual triangles rather than just inversion sets for its small cross-condition controls, feedback-first labels, DFS directed-cycle enumeration, and a fresh memoized Hamilton-path recursion. All arithmetic is exact Python integer arithmetic.

Normal and optimized independent outputs byte-match. The controls include:

- Exact reduced grammar through chain length eight, with counts 1,14,18,24,30,36,42 at lengths two through eight.
- 765 actual-triangle versus inversion-family equivalence checks.
- All accepted reduced presentations through chain length seven checked by actual triangle intersection and Hamilton count; literal permutations also checked through chain length five.
- 28 direct DFS odd-cycle table and independent Hamilton-count checks for the four ray families.
- Literal permutation verification of H=29 and H=87, including disjoint-triangle pairs for the former and none for the latter.
- Direct one-apex order-seven spectrum {25,27,31,33,37,39,45,47,51}.
- Dyadic source-word, involution, orbit-filter, exact ray-count, and numerical-spectrum controls for L=5..14, beyond the producer's L<=12 controls.

The finite bounds are explicitly recorded; the accepted all-order claims rest on Sections 1--4, not these finite tests. No mathematical correction was required during this audit. Only the missing explicit OCF primary citation was added before final acceptance.
