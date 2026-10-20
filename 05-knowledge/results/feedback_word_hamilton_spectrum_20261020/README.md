# Session 16: a proof-based spectrum grammar and its sharp one-apex boundary

**Status: PROVED (elementary constructions and obstructions) + FINITE-EXACT independent controls.**
The assertion that every positive odd integer other than 7 and 21 is a tournament Hamilton-path count remains **OPEN**. This packet does not turn a connected metagraph, a recurrence, or a finite prefix into an existence theorem.

## Results worth carrying forward

1. A tournament whose deletion of one distinguished vertex is transitive has an exact geometrically weighted binary-inversion formula for its ordinary Hamilton-path count. The formula has a four-coordinate nonnegative integer transfer representation.
2. Adjacent exchanges in its binary-word metagraph give a strict rearrangement theorem. Among **strong** such tournaments on `n>=5` vertices, the minimum is
   
   `(2^floor((n-3)/2)+1)(2^ceil((n-3)/2)+1)`.
   
   The maximum is
   
   `(2^floor((n-1)/2)-1)(2^ceil((n-1)/2)-1)+2`.
   
   Both bounds are attained by explicit words. The minimum is 81 at order nine, whereas the unrestricted strong minimum is 75. This explains a genuine construction-class boundary, not an arithmetic obstruction for all tournaments.
3. For every `a,b>=1`, a **strong** tournament on `a+b+3` vertices has ordinary count `(2^a+1)(2^b+1)`. The scalar factorization is not an SCC factorization.
4. The natural adjacent-exchange metagraph has **no ordinary-H unit step of two** on the strong one-apex interior once `n>=6`. Its connected fixed-content word fibres therefore do not give odd-interval coverage by a discrete intermediate-value argument. The full exact jump is proved below.
5. The entire one-apex class, over all orders, omits the globally attainable value 13. Consequently one-apex gadgets, SCC products of them, and source/sink padding cannot prove the full spectrum conjecture: no product of their values can first create the prime 13.
6. For any fixed feedback-vertex bound `k`, all tournament path counts in that class are exactly the outputs of a finite, explicitly specified, nonnegative integer word automaton on `2^k` letters, with at most `k! 3^(k+1)` states. This includes all interactions and all orphan paths. At `k=2`, there are at most 54 states. It is an exact constructive grammar, not a completeness claim for its numerical image.

The full proofs are in [proofs.md](proofs.md). The exact verifier imports no repository mathematics implementation. Normal and optimized runs agree byte-for-byte. No large tournament census was run.

## Inheritance / scope

- Anchor: `THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`, the two infinite strong prime lanes, still unclosed.
- Closest mechanism: `THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md`, which correctly requires state-labelled unit moves, not merely an additive gcd.
- Hostile: the prime 13 in the one-apex grammar; the zero/unit-step failure proved here; the inherited regular seven-carousel example has a different target class.
- Corrected near miss: THM-029's historical Helly step and THM-079's historical poisoning-DAG proof. The incoming `spectrum_basins_20261020e` packet already repairs or warns on these truth surfaces. This packet makes no duplicate shared-repo edits.
- Least-used sidecar: complete gap assignment around the ordered feedback vertices. It preserves ordinary path multiplicities, unlike a selected matching or `H/Aut(T)`.

All results concern **ordinary vertex-word Hamilton paths**. A relabelling does not change their count. No regularity restriction or automorphism quotient is silently imposed.

## Literature boundary

The literature check found a directly relevant 2026 paper by Maria Abi Aad, *About the number of directed Hamiltonian paths in special tournaments*, DOI [10.1016/j.dam.2026.07.017](https://doi.org/10.1016/j.dam.2026.07.017), indexed with issue date 15 December 2026. Its visible primary-source abstract and propositions describe formulas for special reversed-arc families, extending the 2020 paper [10.1016/j.dam.2019.09.005](https://doi.org/10.1016/j.dam.2019.09.005). It does **not** assert all-odd spectrum completeness. Its inherited Proposition 1 overlaps the one-apex counting territory here; therefore no novelty claim is made for the basic counting formula. The extremal/exchange consequences and word-automaton formulation are independently proved here; no literature-priority claim is made for them either.

The author is verified by the author preprint record [SSRN 5804403](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5804403); a social-media search extraction had misleadingly associated another person with the title and is not used as authorship evidence. The publisher pages were searchable, but direct full-page opening returned 403/redirect errors. Mathematical formulas displayed in search extraction contain obvious indexing/OCR problems (one line gives `2^(n-1)+1` for a nearly transitive tournament, contradicted already by the directed triangle). They are not used as proof inputs. The present formulas are derived from path words and independently checked by subset DP and literal permutation counts.

For the current published/open status and the correct Moon–Busch exclusions, use the inherited `rotation_owners_spectrum_20261020e/spectrum.md` and the incoming correction packet. This session has not found or claimed a definitive resolution of the global conjecture.

## Reproduction

Run from the parent workspace:

```text
python session16-hamilton-spectrum/verify.py
python -O session16-hamilton-spectrum/verify.py
```

Receipts: `verify.out`, `verify_optimized.out`. Explicit exceptions, not Python assertions, enforce checks under `-O`.

Universe: 511 one-apex words of chain length at most eight (ordinary order at most nine); exact extrema and 4,097 adjacent exchanges through chain length twelve; 256 symbolic-family parameter checks; 1,272 general feedback automaton controls (`k=0,2,3`, total order at most six). Literal permutation checks are independent of subset DP. The displayed all-order results are proofs, not extrapolations from these checks.
