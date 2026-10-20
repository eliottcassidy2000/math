# Independent review: blue maxima continuation

Status: **APPROVED**. No unresolved mathematical or computational blocker was found.

Reviewed on 2026-10-20. Producer files were not edited by this review; all replays were performed in an isolated copy.

## Approved claims

- `M_6=1185` and `M_7=8295`, with a unique maximizing hollow-frame type at each order.
- Exactly 12 maximum-blue tournament classes at order 13 and 24 at order 15; triangle-parity shore splits are respectively 6+6 and 12+12.
- `liminf M_k/k! >= sqrt(e)`.
- The explicit one-sided rates `M_k/k! >= sqrt(e)(1-101/(8k))` for `k>=4`, and `M_k/k! >= sqrt(e)(1-13/(8k))` for `k=4 (mod 12)`.
- Certified lower-bound witnesses `M_8>=66123`, `M_9>=597249`, and `M_10>=5975877`. No global optimality is certified for these three orders.

## Mathematical audit

The signed-permutation randomization preserves the hollow anti-involution. Independent local enumeration confirms probabilities 1/2 for one edge, 1/4 for disjoint edges, 5/16 for a directed wedge, and 3/16 for a three-edge path. The wedge proof's excluded-path count is exactly two per unordered pair of parts, hence twelve in total.

Every intersection of an index Hamiltonian path with a four-index block is one of the five listed linear forests. Consequently its probability correction is exactly `1+X_B/4`, including the case of four consecutive indices. Independent block samples therefore produce the stated exact expectation identity. Vertex-disjoint triples have distinct parent blocks, occur once in the nonnegative expansion, and have simultaneous consecutiveness probability `6^r(k-2r)!/k!`. The greedy degree estimate and fixed-truncation limit then prove the liminf claim without exchanging an unbounded expectation and a limit.

The quantitative estimates in Section 6 are valid across their stated ranges. For the uniform rate, `C_r<k` ensures that the retained term exists and every product factor lies in [0,1]; the lower bound on `4q` is valid for `k>=13`, while the remaining finite range has a negative asserted lower bound. For admissible orders, the factorization involving `(1-1/k)(1-9j/k)` is exact. The constants 101/8 and 13/8 are the correct weighted sums of the respective quadratic errors.

The primary publisher abstract explicitly states the required existence theorem for resolvable `(v,4,1)` designs for every `v=4 (mod 12)`: [Hanani, Ray-Chaudhuri and Wilson, On resolvable designs](https://www.sciencedirect.com/science/article/pii/0012365X7290091X). This review checked the theorem statement, not the original existence proof. Only the underlying pairwise block design is used here. The transfer from `max h(D)` to maximum blue mass continues to use the inherited regular-subset theorem, as disclosed in the packet.

## Fresh complete computations

All C++ executables were recompiled from source with `g++ -O3 -std=c++17`.

- At `k=6`, the complete forward-subset-DP census, direct subcube replay, and Walsh computation agree on all 1,835,008 reduced frames. These cover all 32,768 tournament masks and all 1,073,741,824 labelled external frames.
- At `k=7`, the fresh Walsh producer and ordinary-monomial subset-zeta replay agree pointwise on all 956,301,312 reduced frames. Both exited successfully. The replay covers all 2,097,152 tournament masks and all 4,398,046,511,104 labelled external frames.
- Every non-timing field in the freshly regenerated exact JSON summaries matches the retained summaries. The representative text files agree exactly, including order, with the JSON representative lists.
- The full Python verifier, including orbit coverage, direct word counts, moments, completion orbits, and the local rational checks, was rerun successfully. Orbit coverage was not skipped.

The arithmetic and transform identities were also inspected. Integer capacities are sufficient at the supported orders; the order-15 weighted second moment uses 128-bit accumulation.

## Separately written checks

`independent_audit.py` does not import producer code. It reconstructs adjacency from the mathematical bit definition and counts words using a backward memoized recurrence. It confirms the exact witnesses, all three search witnesses, the local distribution, and the frame automorphism groups.

It directly computes triangle parities for all matching-completion orbits. It also performs exact full-tournament isomorphism searches between every pair of completion-orbit representatives: all 120 pairs at order 13 and all 496 pairs at order 15 are nonisomorphic. Thus the completion counts are not merely assumed from the group-orbit summary. Results are retained in `independent_audit.json`. Run the script with ordinary `python3`, not `python3 -O`, because this separate check uses assertions.

## Scope limits

This approval establishes the stated finite maxima and one-sided lower bounds. It does not establish a general closed form, a universal `sqrt(e)` upper bound, asymptotic equality, or optimality of the order-17, order-19, or order-21 witnesses. No literature-priority claim was reviewed or made.

## Reviewed mathematical text hashes

SHA-256 of `SQRT_E_LOWER_BOUND.md`:

`5d3118f045e1b636a245bc69604d73c14e27a7b652e1bd1d24a67597834ecb87`

SHA-256 of `BLUE_MAXIMA_K6_K7.md`:

`f3eb53c7298941fcfbf67b9de49ee8cfe9bcbf8f4b402ca8ca6f1099c2846e51`

## Final-text revalidation

The final `SQRT_E_LOWER_BOUND.md` differs from its previously reviewed version only in its approval status: restoring the earlier status string reproduces the earlier SHA-256 exactly.

A status-only comparison could not be established for `BLUE_MAXIMA_K6_K7.md`. Its exact final text was therefore reread and freshly checked against the approved computations and proofs. It is approved at the final hash above. This revalidation includes the already independently counted `k=10` witness, its lower-bound-only qualification, the finite-pattern comparison value 5,982,795, and the correct distinction between trivial matching stabilizer and a regular tournament completion. Direct degree checks also confirmed regular-completability of all displayed `k=6,7,8,9,10` witnesses. No claim that this document changed only in formatting or status is made.
