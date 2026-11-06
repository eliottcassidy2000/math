# Independent audit: cover equality, sharp-defect normal form, and fixed-start sparsity

**PASS in the stated scope.** This audit covers the frozen `reviewed-inputs/PROOF.md` and `DENSITY.md`. It does not audit the separate bounded-defect extension and does not prove F(8,137) or an eventual Hamilton-spectrum threshold.

## Exact input scope

- Main proof SHA256: `680ee2d8e97f700d33b008e2a42f9e61a9287da03d25668e4b7d430e8db72ea7`.
- Density addendum SHA256: `39fc861102006e76265f703640c9ecfc0b54d2d9ad29e8d67532b7603bf22f45`.
- Producer verifier SHA256: `745fb1e78bd221d109a0a722f64af2c1a303c2a421e9a1bc144bf89e39064d93`.
- The sharp-defect dependency is the corrected nontrivial-source version, SHA256 `582687a3d966a928efb7355cca6d9f1181d76d09d4065da17a0815d1343f457d`. Its original singleton-boundary error and exact correction are preserved. That separately audited theorem is now published at https://github.com/eliottcassidy2000/math/blob/20d185c937163242289918282869826d0e1ddb75/session203/sharp-z-defect/PROOF.md (publication changes its status only).

Here H(empty)=1; path bins are labeled and may be empty unless designated as nonempty. The principal L=h theorem permits the singleton. The D=s normal form concerns nontrivial positive-start sources, so it does not identify the singleton's formula-defined D=0 with its physical Z difference 1.

## All-order argument review

1. **Two-terminal cover lemma.** If both terminal-deletion alternatives failed in a strong tournament S, the last SCCs of S-x and S-y would both survive as sink SCCs of S-{x,y}. They coincide. Strongness of S forces an arc from that common sink to x, while its sink property in S-y forbids the same arc. The common component is nonempty. The reduction from arbitrary tournaments to their ordered SCCs correctly gives two nonempty paths; omitted intermediate components are assigned to the other path and all concatenations go forward.
2. **Exact excess decomposition.** Cutting immediately after v bijects Hamilton paths with the L-covers whose second path is empty or begins at an outneighbor. Every remaining cover, after deleting terminal v, is either a Hamilton path starting at an inneighbor, or a two-nonempty-path cover with end/start in two distinct inneighbors. Every unordered pair supplies a distinct endpoint type. Consequently L-H is at least the indicated sum plus binomial(indegree,2); for indegree zero or one the two-path contribution is exactly zero. This proves the claimed iff, not just a necessary condition.
3. **Complete physical normal form.** Apply the independently audited unique-outneighbor equality test to a=the unique successor of v. The new L-equality classification leaves either no inneighbor of a or exactly one x. In the latter case end_a(U-v)=end_x(T), and the required start_x(T)=end_x(T)=0 means precisely an interior SCC. Deleting initial v,a gives start_v(U)=H(T). Conversely each described graph satisfies the equality test. T is nonempty, and the explicit reachability argument proves strongness. No additional unlisted graph form survives.
4. **Optional edge reversal readout.** In the unflipped graph, a Hamilton path lost by reversing a->x must have v,a consecutive: the other possibility would require a Hamilton path of T starting at x, which is excluded. The lost paths correspond to R_x(T). Newly gained paths have the indicated L_x(T) decomposition; its trailing T path is nonempty because end_x(T)=0. For two labeled covers, deleting a from a newly admitted x->a segment creates one designated x-ending path and two labeled ordinary bins, giving exactly 2E_3,x. The lost-cover factor 2 similarly records which bin contained a. Thus all displayed h,e,d,b and B+t formulas have literal inverse maps.
5. **Marked SCC factorization.** Restricting each labeled path to an ordered component and concatenating component pieces is an exact bijection, including empty pieces. The x-start path cannot meet P; the x-end path cannot meet Q. This gives exactly the three formulas for c,r,z, with C3(P), rather than C2(P), in z.
6. **Finite two-parameter family.** At fixed positive t, each nontrivial strong component has 3<=H<=t and order<=H. There are finitely many types and at most floor(log_3 t) nontrivial components. Arbitrarily many singleton components affect C2 and C3 only through their total numbers before and after the marked component. A marked singleton is retained, not counted in either number. Empty retained cores have cover count one. The constants alpha,beta,gamma are positive; disjoint events that either labeled bin starts at x give C2(M)>=2R_x(M), and C2(P_core)>=H(P_core), so alpha>=2gamma. Therefore no cancellation hides large k,l: B+t>=beta*2^l*3^k. This bounds the number of distinct B values by O_t(log^2(X+2)), regardless of repeated parameterizations.

The finite union over 1<=t<=K and any fixed finite shifts remains zero-density. This conclusion excludes only equality-face sources as a complete bounded-start strategy. General D>s sources, arbitrary profile occupancy, and the ordinary spectrum remain outside it. For t not in the ordinary spectrum the set is empty, consistent with the bound.

## Independent exact controls

All checks use explicit failures, so debug or optimized Python mode cannot remove validations.

- `check.cpp` enumerates actual directed path prefixes on all subsets, independently of the producer's endpoint subset recurrence. It checks all 33,867 labeled nonempty tournaments through order 6, all 202,013 roots, all 502,170 terminal pairs, and 340,390 strong terminal-deletion alternatives. There are 15,959 cover equalities and 3,670 sharp-defect equalities. It reconstructs the normal form of every sharp equality source and checks strongness.
- The same independent code constructs all 33,867 baseline sources and all 10,694 admissible flipped sources. It verifies all numerical formulas, SCC factorization, and the exact two-parameter expression. This includes 6,254 marked-singleton cases. It does not merely compare predicted counts with another evaluation of the same formula.
- `check_bijections.py` enumerates literal permutation/cover sets. It checks the exact bad-cover deletion map on all 5,405 marked inputs through order 5 (31,322 bad covers), plus 54 physical flips through body order 4. It checks equality of the lost/gained Hamilton-path sets and both labeled two-cover sets, including multiplicities and empty bins.
- `check_templates.cpp` separately checks 72 literal sources with both unmarked retained cores nonempty, marked component either a singleton or directed triangle, singleton padding placed on either side of each core, and total order up to 15. It verifies the same formula with actual graph counts, not just symbolic arithmetic.
- Independent C++ optimized/debug outputs agree; Python normal/-O outputs agree. The frozen producer optimized/debug replays also reproduce its saved receipts byte-for-byte.

## Reproduction

Run from this audit directory:

```sh
g++ -std=c++20 -O3 check.cpp -o /tmp/cover-equality-audit && /tmp/cover-equality-audit
g++ -std=c++20 -O0 -D_GLIBCXX_ASSERTIONS check.cpp -o /tmp/cover-equality-audit-debug && /tmp/cover-equality-audit-debug
python3 check_bijections.py
python3 -O check_bijections.py
g++ -std=c++20 -O3 check_templates.cpp -o /tmp/cover-equality-templates && /tmp/cover-equality-templates
g++ -std=c++20 -O0 -D_GLIBCXX_ASSERTIONS check_templates.cpp -o /tmp/cover-equality-templates-debug && /tmp/cover-equality-templates-debug
g++ -std=c++20 -O3 reviewed-inputs/verify.cpp -o /tmp/cover-equality-producer && /tmp/cover-equality-producer
sha256sum -c SHA256SUMS
```

No mathematical correction to these two frozen proof files was required. Status-only promotion with preserved originals may proceed through the separate release review. The producer input manifest and all audited inputs are retained under reviewed-inputs/.
