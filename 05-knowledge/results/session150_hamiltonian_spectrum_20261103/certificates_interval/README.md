# A verified bridge in the tournament Hamiltonian-path spectrum

**Result (finite, exact).** Every odd integer from 80,407 through 80,873 occurs as the number H of directed Hamiltonian paths of a strongly connected tournament of order 12. The packet contains one explicit labelled witness for each of the 234 values. It is an existence certificate, not a census of all order-12 tournaments.

Combining this bridge with the inherited results in `01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md` gives every odd H ≤ 84,259 except 7 and 21: THM-4104 supplies the allowed odd prefix through 80,405 and the entire odd interval 80,875 through 84,259. The omitted values 7 and 21 are known impossible (THM-1370). This finite extension does not settle the all-order spectrum conjecture.

**Encoding.** `witnesses.txt` contains one row `(H, parent, cut)` per value in ascending order. The parent has vertices 0–10. The 55 bits of `parent` correspond to pairs `(i,j)` with `0 ≤ i < j ≤ 10` in lexicographic order, least significant bit first. A set bit means `i → j`, and a clear bit means `j → i`. Add vertex 11: cut bit i set means `i → 11`, otherwise `11 → i`.

**Independent replay.** `verify.py` constructs the directed adjacency matrix directly from the bits, counts Hamiltonian paths using start-vertex subset dynamic programming, and independently checks strong connectivity using transitive closure. `check_interval.py` calls it for each record and checks exact interval coverage, uniqueness, ordering, and input ranges. Neither checker uses the search recurrence. Run `python3 check_interval.py` and `python3 -O check_interval.py` from this directory. Both produce the same PASS line; no correctness gate uses `assert`.

**Search scope.** The search sampled 15,000 labelled order-11 parents with independent edge-bit probability 0.25, then tested nonconstant one-vertex cuts through an endpoint-path and two-path-cover counting identity. The witness verifier does not depend on that search or on its formulas. This was not an exhaustive search; its random-parent horizon makes no nonexistence claim.
