# Square sum coverage and the two feedback tournament spectrum

**Status: PROVED with finite certificates for the square-sum proof chain and the scoped algebraic results below; FINITE-EXACT for tournament coverage through 99,999; global tournament completeness OPEN. Same-agent review, with independently implemented consumers.** Ordinary interval square-sum Hamiltonicity is Gerbicz's known theorem, not a new resolution claimed here. The permanent tournament omissions 7 and 21 were already proved. This packet makes the square-sum proof locally replayable and sharpens the construction class in which to pursue the remaining tournament problem.

## Inheritance and concept board

| Lane | Closest mechanism | Hostile or corrected boundary | Sidecar and next use |
|---|---|---|---|
| Anchor | THM-4094, [two prime lanes](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md) | A finite multiplicative bank never creates a new prime | Retain the feedback word and all continuation responses |
| Anchor | [Exact feedback-gap automaton](../feedback_word_hamilton_spectrum_20261020/proofs.md) | One feedback vertex misses 13 at every order; bounded run templates miss infinitely many primes | Unbounded words on four letters, two marked vertices |
| Niche | [Seeded square-sum theorem at threshold 6242](../seeded_optimal_threshold_6242_20261019/OPTIMAL_THRESHOLD_THEOREM.md) | Individual edge optionality does not imply compatibility of several edges | Named endpoints, positional parity, and complete glue words |
| Wildcard | Reachable/observable quotient of a weighted word automaton | Scalar H and even (order,H) discard future responses | Exact rational Hankel rank and actual continuation counts |
| Shared mechanism | A complete transition table plus well-founded descent | A recurrence without a universally covered input domain is insufficient | Square-sum integer division closes; tournament inverse entry remains missing |

The [current ordinary-spectrum audit](../rotation_owners_spectrum_20261020e/spectrum.md), not a historical session log, supplies the status and correction lineage. THM-1370 / [h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) excludes the two raw counts. THM-4812 / [presentation-fibres-7-and-21](../../../01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md) concerns the distinct quotient H/Aut. Nothing here identifies these spectra.

## A complete local proof chain for ordinary square sum cycles

Let Q_N have vertices 1,...,N, with distinct vertices joined exactly when their sum is a square. There is a Hamiltonian cycle in Q_N for every N >= 32.

The proof has three explicit pieces:

1. **32 <= N <= 6242.** The new `square_bridge.json.gz` contains a 32-cycle followed by 6,210 insertion/reversal records. The independent consumer reconstructs and checks **every one of the 6,211 cycles in full**, checking the vertex set and every cyclic sum. It performs 19,531,376 checks. There are 28,834 reversals, 925 direct insertions, and no fallback reset cycles. Search used no theorem about success of local improvement; only the saved witnesses are used in the proof.
2. **6242 <= N <= 128452.** The inherited [optimal-threshold packet](../seeded_optimal_threshold_6242_20261019/OPTIMAL_THRESHOLD_THEOREM.md) supplies every intermediate cycle by a checked extension trace preserving 482–3362–6242. Its 122,211-length replay passed again, with 7,147,854 checks.
3. **N >= 128453.** The inherited [eventual-cycle proof](../seeded_eventual_cycle_theorem_20261019/SEEDED_INTERVAL_THEOREM.md) supplies a well-founded base-49 construction. Its 2,568 parity-compatible base pairs at indices 53 through 2620, all 294 pair-transition states, and all 392 final-cycle states passed replay. The induction reduces n to floor((n-24)/49), and every residue 24,...,72 and both input parities are covered. This is an infinite proof using finite certificates, rather than extrapolation from the bridge.

The finite traces need not meet at the same **cycle** at N=6242; they prove existence for overlapping integer ranges. Requiring equality of those cycle representatives would add an unnecessary condition.

Cutting one cycle edge gives a path for every N >= 32. The inherited [complete small-order inventory](../square_sum_ports_20261023/README.md) supplies paths at 25 through 31, and the other small cases. Thus, including the one-vertex convention, a square-sum Hamiltonian path exists exactly for

    N in {1,15,16,17,23} or N >= 25.

The all-size existence result is credited to Robert Gerbicz, whose own [A090461 comment](https://oeis.org/A090461) records the thresholds. The original forum link returned 403 during this session; its unread contents are not a proof dependency of the local chain above. The stronger repo theorem already preserves the exact seed at every N >= 6242, and the [all-size connected-path packet](../square_sum_all_size_connected_paths_20261019/ALL_SIZE_CONNECTED_PATHS.md) preserves 2,303 specified paths at a larger threshold. These stronger inherited results must not be replaced by the weaker ordinary theorem in future summaries.

## Two feedback vertices already cover a substantial finite interval

Fix marked vertices x,y with x -> y. The remaining vertices v_1,...,v_m form a transitive chain, v_i -> v_j for i<j. Letter a_i in {0,1,2,3} has bit r set exactly when marked vertex r points to v_i. Let H(w) count ordinary Hamiltonian vertex words of the resulting tournament.

The scanner visits all normalized prefixes through length 20 whose current H does not already exceed 100,000. It can omit a leading 0, since that chain vertex is a source. Pruning uses the insertion injection H(T) <= H(T+v). These optimizations explain the search, but **positive coverage does not rely on exhaustiveness**: every target has its own saved word.

`two_feedback_100000.json` contains **49,998 explicit words**, one for every odd number from 1 to 99,999 except 7 and 21. The independent Python consumer checks every word's count using the proved gap-assignment bijection. An independent vertex-subset DP checks all 341 words of length at most four; literal permutation counting additionally checks the 85 words of length at most three. In particular, the target

    80557 = H(10000020300312)

has an independent subset-DP check on its 16-vertex tournament. The fixed two-feedback condition is essential to the scope of this positive signal. No claim is made that these witnesses have minimum order, or that this is a record for unrestricted finite coverage.

**OPEN strengthened target:** Every positive odd number other than 7 and 21 is H(w) for this single four-letter, two-feedback language. A proof would imply the unrestricted raw-H spectrum conjecture. Failure of this stronger target would not refute the unrestricted conjecture. The current positive data justify pursuing it; they do not prove it.

## Exact minimal linear memory is seventeen

**Proposition, proved by exact finite linear-algebra certificates.** The scalar series w -> H(w) has rational Hankel rank 17. Consequently the minimum dimension of a rational weighted linear representation is exactly 17. In characteristic two its rank is one; in every odd characteristic it is 17.

For a self-contained proof, use the inherited gap semantics. For each of the two orders of x,y, a path has three chain gaps. Record empty/good-last/bad-last in the first two gaps and empty/nonempty in the last. This gives 36 states. The letter matrices M_a, initial row u and terminal column v count every gap assignment with its multiplicity, so H(w)=u M_w v.

The exact reachable row space has dimension 28; the exact observable column space has dimension 20. `minimal_feedback.py` grows each space until it is closed under all four letter matrices. It selects 17 observable columns O and 17 reachable rows B, with K=B O nonsingular, and computes

    A_a = K^(-1) B M_a O.

The checker verifies R M_a O = R O A_a on a basis R of the full reachable space. Thus, for every word, the 17 coordinates

    X(w) = (H(w s_1),...,H(w s_17))

update exactly by X(wa)=X(w) A_a. The suffix s_1 is empty, so the first coordinate is the desired H(w). The explicit suffixes, prefixes, initial vector and four matrices are frozen in `minimal_feedback.json`.

The corresponding 17-by-17 Hankel minor has determinant

    150994944 = 2^24 * 3^2.

Every linear realization factors that minor through its state space, so no dimension below 17 suffices. `audit_minimal_feedback.py` reconstructs all 289 entries with **vertex-subset DP**, without using the gap or matrix constructor, and checks the determinant independently. It also checks all 17 continuation coordinates on all words through length four.

All denominators of the displayed realization have prime support contained in {2,3}. Therefore the exact characteristic-zero identities reduce modulo any prime p>=5 and give dimension at most17; the same nonzero minor gives dimension at least17. At p=3 the direct finite-field closure computation gives rank17; at p=2 it gives rank1, in agreement with H(w) always being odd. The checked finite-field results are in `modular_feedback.json`. Extending scalars does not change these ranks.

This is a **linear-state dimension**, not a lower bound on bits, arbitrary nonlinear encodings, or graph size. The displayed observation-coordinate matrices have negative rational entries. The following integral change of basis removes their denominators, while preserving the distinction between a lattice point and an actual word. The unreduced 36-state model remains available when nonnegative path-assignment counts matter.

## One extra prefix completes the integer observation lattice

Let L be the integer span of all observation rows X(w), and L_0 the span of the 17 prefix rows used in the nonzero minor K above. These are lattices in the specified 17 continuation-count coordinates, not sets of attainable single-word states.

**Proposition, with exact certificate.** L/L_0 has order two and is generated by X(203). The words 100 and 203 both have length three and H=5, but X(100) lies in L_0 and X(203) does not. Thus the extra class is invisible to both order and scalar path count.

The certificate `integral_feedback.py` constructs an integer basis B of span_Z(L_0,X(203)) and verifies:

    |det K| = 150994944,
    |det B| = 75497472,
    X(203) not in L_0,       2 X(203) in L_0,
    B A_a B^(-1) is integral for every letter a,
    X(empty) B^(-1) is integral.

The index and membership checks independently certify the one-step enlargement, rather than trusting that a normal-form routine produced the right generated lattice. Invariance and the initial-state check imply that every reachable row lies in this lattice. Conversely, its generators are actual prefix rows, so it lies in their full integer span. Hence it is exactly L.

The Smith invariant factors of this **specified observation-coordinate lattice** are

    1, 2,2,2,2,2,2,2,2,2,2,2,2, 4,4,8,144.

The factor144 reflects this choice of continuation probes; it does not contradict full series rank17 modulo3, where other continuations distinguish the missing mode. In normalized coordinates Y=X B^(-1), the four transition matrices, initial row, and output column are all integral and are saved in `integral_feedback.json`. This gives an **integer** linear realization of the same minimal dimension17. Some coefficients are negative.

This is a precise one-bit completion of a finite-index lattice, not a one-bit completeness proof for the spectrum. The hostile is explicit: 7 X(empty) is a positive integral observation row in L, with first coordinate7, yet no tournament realizes it. The same holds for21 X(empty). Therefore integrality, all the lattice congruences, and positivity of the observed continuation counts are jointly insufficient. A successful inverse construction must remain in the actual matrix orbit of the initial row, rather than its linear span. This is the tournament analogue of retaining one connected cycle rather than only a feasible degree cover in the square-sum gluing problem.

## Every nontrivial fixed scalar context fails

The following obstruction applies at every feedback size k>=1 and every fixed tournament F on the marked vertices. Put z=0 and o=2^k-1 for the all-zero and all-one chain letters. Fix words P,Q.

**Proposition.** There is a function Phi satisfying

    H(P w Q) = Phi(H(w)) for every word w

if and only if P consists entirely of z and Q entirely of o. In that case H(PwQ)=H(w), so Phi is the identity on the attained spectrum.

**Proof.** The easy direction inserts source vertices before and sink vertices after the tournament. For the converse, first suppose P contains a letter with some bit one. Let w=z^m. All these inputs have H(w)=H(F), independently of m. Choose that prefix vertex b and its marked vertex x. The induced tournament on b,x and the m chain vertices has b -> chain -> x -> b; the chain is transitive. Its path count is 2^m+1. Hence H(Pz^mQ)>=2^m+1 by monotonicity under adding vertices, contradicting a fixed value Phi(H(F)).

If Q contains a letter with some bit zero, instead set w=o^m. Choose that suffix vertex b and a marked vertex x that it beats. The induced tournament has x -> chain -> b -> x, again with 2^m+1 paths. The same contradiction follows. Therefore only the stated padding contexts survive. The count 2^m+1 is the near-transitive formula, or follows directly by assigning the transitive block across the return arc. QED.

For pure appending this is stronger still: H(z^m)=H(o^m)=H(F) and the input orders agree, while H(z^m Q)=H(Q) is constant and H(o^m Q) grows without bound for every non-padding Q. **Even (order,H) cannot determine the appended result.** This generalizes the obstruction behind the earlier bounded short-context probes, but does not rule out controllers on proper state-restricted languages, contexts chosen from the actual source state, or higher-dimensional recurrences. Those are the intended reroutes.

## The precise bridge between the two problems

The square-sum source is a pair of actual paths with named endpoints and a shared positional-parity map. The target is a larger interval cycle; affine block lifting preserves internal square sums and glue words certify the remaining edges and single-cycle connectivity. Forgetting the parity map destroys the block-choice coverage proof.

The tournament source is a word with a marked two-vertex feedback frame. The target is its continuation-count function; the gap map preserves every Hamiltonian-path multiplicity. Forgetting that response and retaining H alone destroys composition, by the proposition above. Neither construction is a bijection between square-sum graphs and tournaments.

What square-sum coverage has and tournament coverage lacks is a **complete inverse address rule**. Division by49 reaches a certified base while retaining every port condition. The tournament matrices are exact, but there is not yet a proved terminating procedure taking each allowed target h to a realizable state with first coordinate h. Establishing that procedure, or finding its first unavoidable obstruction, is the overnight anchor.

## Reproduction

Run from this directory:

```text
python3 verify_square_bridge.py
python3 -O verify_square_bridge.py
c++ -std=c++17 -O3 two_feedback_scan.cpp -o /tmp/two_feedback_scan
/tmp/two_feedback_scan 20 100000 > /tmp/two_feedback_100000.json
python3 audit_feedback.py
python3 minimal_feedback.py
python3 audit_minimal_feedback.py
python3 modular_feedback.py
python3 integral_feedback.py
```

The bridge discovery script additionally requires OR-Tools; no square bridge consumer requires it. Exact matrix work uses SymPy. The inherited induction dependencies and replay receipts are listed in `INHERITANCE.json`; the saved finite bridge alone is not the infinite theorem. The handoff records the current live targets and the overnight cutoff.
