# The missing shared-dictionary case of the three-color triangle

Status: **PROVED** under the explicit retained-zero and deleted-zero grammars below; finite checks are supplementary. This closes the shared-dictionary question explicitly left separate in Session41, `fibonacci/BRANCH_ADDENDUM.md`, and Session42, `fibonacci-addition/RESULTS.md`. It is a theorem about marked Fibonacci representations, not a theorem that constructs arbitrary tournament Hamilton-path totals.

## 1. The recovered construction and the new quantifier

Use Fibonacci values `F_2=1,F_3=2,F_n=F_(n-1)+F_(n-2)` and three distinct unit slots `A=F_(-1),B=F_1,C=F_2`. All selected slots are distinct and nonadjacent. In grammar R the zero slot 0 stays present and unselected, so A and B are compatible. Grammar D deletes it and reconnects A to B, making that pair incompatible. The other adjacencies are unchanged.

Session41's canonical proper support of F_n expands into an ascending unit word W_n. Its reverse V_n satisfies `V_3=CA,V_4=CAB,V_n=V_(n-1)V_(n-2)`, with infinite limit V. Integer row N is the reversed prefix of length N, so its cell j has color `V[N-j+1]`: constant diagonals have constant colors. The substitution is `C→CA,A→B,B→CA`. This is the requested integer-row triangle; it is different from Pascal's triangular array or an F_2^2 additive charge.

Here a **shared dictionary** chooses exactly one legal proper support for each F_n, n≥3, and reuses its recursively expanded word at every occurrence. At each integer N the *outer* legal support may vary. All its Fibonacci summands then use the same dictionary. We ask for a coherent infinite stripe X whose prefix of length N is an expanded legal support of N for every N. This retains the intermediate synchronization quantifier that the occurrence-dependent language theorem deliberately did not classify.

The inherited complete proper-support list is particularly simple. In reverse order, the two R choices are:

* odd n: the blocks `n−1,n−3,...,4`, then `CA` or `BA`;
* even n: the blocks `n−1,n−3,...,3`, then `B` or `A`.

For D the odd choice is forced to be `CA`; the even choice remains free. Every proper support for n≥4 contains its largest slot n−1; n=3 is the explicitly displayed CA/BA boundary. Let V_n denote its expanded reverse word under the chosen dictionary, and let V_A=A,V_B=B,V_C=C.

## 2. Every dictionary works, and its stripe is unique

**Theorem 1.** Every independent sequence of the indicated legal choices produces exactly one coherent infinite stripe. Its prefix of length F_n is V_n for every n≥3.

**Existence.** For n≥4 every proper support of F_n begins with slot n−1 in reverse order, so V_n begins with V_(n−1). The blocks V_3,V_4,... therefore have a unique infinite limit X. Take a prefix of V_n of length N. Traverse its descending proper-support expansion until the prefix ends. Keep every completed block, discard all later blocks, and recursively expand only the final partial block. This process ends at unit slots and returns a legal outer support whose dictionary expansion is that prefix.

The legality assertion needs the actual slot order: when a final slot i is replaced by smaller proper slots, any preceding retained slot has index at least i+2, while every replacement index is at most i−1. No new adjacent pair can appear across that seam. Inside the replacement, legality is inherited; discarding later slots preserves it. The same argument works in R and D, including their distinct unit-slot boundaries. The completed unit word has length N, so the resulting support has value N.

**Uniqueness.** In any coherent stripe, inspect the row of length F_(n+1). Every legal outer support of that value is either its singleton or a proper support. The singleton expands to V_(n+1), and every proper support begins with F_n. Thus every such expanded row begins with the already chosen V_n. Arbitrarily large n force the entire stripe to be X. This argument allows arbitrary legal outer supports; it does not silently require ordinary Zeckendorf support.

The prefix-truncation proof is a useful constructive algorithm: it produces an actual legal support, not merely a color-count equality.

## 3. Exact information growth and owner recovery

Define

`a_n = F_n−1` for odd n, and `a_n = F_n` for even n.

The a_n are strictly increasing for n≥3. The choice for n first appears at position a_n: its two alternatives are C/B for odd n and B/A for even n. The blocks before that terminal choice depend only on smaller indices. Later blocks begin with V_n and cannot alter this position.

Consequently the number of length-N prefixes extendible under shared dictionaries is exactly

`2^r(N)`, where `r(N)=#{n≥3:a_n≤N}`, in R;

`2^d(N)`, where `d(N)=#{even n≥4:F_n≤N}`, in D.                 (1)

Thus their growth is respectively `Θ(N^(log_phi 2))` and `Θ(N^((log_phi 2)/2))`. Both spaces have continuum cardinality but zero prefix entropy. This contrasts with the inherited occurrence-dependent counts `2^(N−a(N))` in R and `2^b(N)` in D, which have positive entropy. In those latter formulas a(N),b(N) count canonical A,B letters, and are unrelated to the first-appearance address a_n.

There is an explicit **owner sidecar** that makes synchronization testable. Build the canonical V_n by the proper-support list above, tag its terminal C by owner n when n is odd, and tag its terminal B by owner n when n is even. Forced A letters have owner 0. Copy tags along with every reused lower block. The nested tagged words define a fixed owner map ρ on positive positions. Then every shared-dictionary stripe is obtained by one binary choice per owner:

* odd owner n: change *all* its C positions to B, or change none;
* even owner n: change *all* its B positions to A, or change none;
* owner 0: keep A.

For D the odd-owner choices are forced to remain C. This is an iff, because the tagged recursion is exactly the shared-dictionary recursion. It gives a finite exact selector: a proposed colored prefix is realizable iff forced positions agree and each visible owner's occurrences request one consistent choice. First-appearance positions a_n recover every visible choice.

For example positions 1,6,9,14,19,22,30,35 all have owner 3. The historical prefix with canonical colors through 34 and B at 35 is therefore impossible under *every* shared dictionary: its first position fixes owner 3 to C. The occurrence-dependent R grammar does allow this one-site change at every horizon. The missing operation is occurrence-specific renewal, not a different numerical value or another global dictionary.

## 4. An exact target selector, and the first synchronization gap

Fix a target N and count unordered pairs of distinct positive A-colored integers summing to N. Let ε_i∈{0,1} be the choice to promote B to A for each visible even owner i. Canonical C positions never become A, even in R. The exact pair count is the Boolean polynomial

`P_N(ε)=r_AA + Σ_i l_i ε_i + Σ_(i≤j) b_ij ε_i ε_j`,       (2)

with `ε_i²=ε_i`. Here r_AA counts canonical AA pairs; l_i counts AB pairs whose B endpoint has owner i; b_ij counts BB pairs with endpoint owners i,j. Pairs involving C contribute zero. The coefficients are nonnegative integers obtained directly from the tagged prefix. The attainable counts are exactly the image of this polynomial on the owner cube. Multiple target constraints retain the *same* ε variables.

The occurrence-dependent Session42 partition product counts each B position independently. Replacing position variables by their shared owner variable is the exact map to (2). This substitution preserves actual pair counts but destroys independent local choices; therefore interval filling need not survive.

**Small exact obstruction.** At N=16, the only potentially A-colored distinct pairs not blocked by C are `3+13` and `5+11`. Positions 3 and 11 are both B with owner 4; their other endpoints are forced A. Hence

`P_16=2 ε_4`,

and the shared-dictionary attainable set is `{0,2}`. Owner 6 occurs only at the excluded midpoint 8 and has no effect. Occurrence-dependent selectors attain `{0,1,2}`. Releasing just one of the two owner-4 occurrences recovers 1. The checker establishes that 16 is the first target with a noninterval image in either grammar; the displayed proof of the gap itself is all-exact and independent of that minimality check.

Even coefficient content1 does not restore filling. At target18 the exact polynomial is `P_18=1+3ε4+ε6`, with image `{1,2,4,5}`: pair5+13 is forced AA; pairs2+16,3+15,7+11 use owner4; pair8+10 uses owner6. Thus the missing3 survives despite gcd(3,1)=1. Keeping the whole selector polynomial is stronger than keeping its coefficient gcd.

There is also a positive retained result: set every even-owner bit to 1. This produces exactly the inherited A_max class (all canonical A and B positions), so the shared-dictionary family still contains the fixed additive basis representing every integer from 7 onward by two distinct A elements. The new gap concerns selectable representation *multiplicity*, not existence of a pair.

**Exact target-selection algorithm.** Given N and a desired pair count k, build the tagged prefix and (2), then enumerate its `2^d(N−1)` even-owner assignments. This decides feasibility and returns an actual shared dictionary when feasible. Tagged blocks through the first Fibonacci length ≥N have total size O(N), so the cost is `O(N+d(N−1)^2 2^d(N−1))`, with O(N) simple storage or a streamed owner parser. Here `2^d(N−1)=Θ(N^γ)`, `γ=(log_phi 2)/2≈0.7202`. This is a numerical-size algorithm, not a polynomial-time claim in log N. Several requested targets can be tested simultaneously by evaluating their polynomials on the same owner assignment; each successful finite assignment extends freely to a complete infinite dictionary. Thus the proof supplies actual target selection and reconstruction, not only an entropy obstruction.

## 5. Typed transfer to target-dependent tournament work

Source: the marked-unit dictionary and its owner map. Target: a repeated-module or marked-port compiler in which several occurrences must reuse one physical chart. Map: substitute a common selector variable for every occurrence owned by that chart. Preserved predicate: the actual target response polynomial evaluated at a legal common choice. Lost information under a positionwise relaxation: equality of choices at repeated occurrences. Necessary sidecar: an owner-incidence map together with the response polynomial, not a list of independent local response sets. Cheap hostile test: two contributions owned by one bit give `2ε`, whereas independently relaxed bits give `ε_1+ε_2`.

This is a rigorous selector operation and a concrete failure test for a proposed tournament compilation. It does not assert that colored row pairs are Hamilton paths. A physical tournament use must supply its own response polynomial and prove that releasing an owner is an allowed port or module operation. The distinction explains how a much larger target-adaptive family can exist even when every globally reused dictionary has only polynomial prefix information.

## 6. Provenance and replay

Inherited grammar and canonical grid: `05-knowledge/results/session41_owned_attachments_fibonacci_20261023/fibonacci/RESULTS.md` and `COLORED_GRID_AND_PHASE_QUOTIENT.md`. Inherited occurrence-dependent classification and explicit shared-dictionary boundary: its `fibonacci/BRANCH_ADDENDUM.md`. Inherited independent-position pair product and A_max basis: `05-knowledge/results/session42_shared_lattices_square_colors_20261023/fibonacci-addition/RESULTS.md` §§5–7. The new statements here are Theorem 1, the exact shared-dictionary prefix counts/owner classification, and their shared-owner target polynomial and gap.

Run `python3 check_shared_dictionary.py` and `python3 -O check_shared_dictionary.py` from this directory. The checker independently enumerates legal slot supports, checks every R dictionary through F_11 and every D dictionary through F_12, verifies every row prefix, exact prefix counts, owner reconstruction, the historical-35 obstruction, and the target-16 selector gap. It uses explicit RuntimeError guards, never removable assertions. Exact small computations supplement the proofs; no tournament census is performed.
