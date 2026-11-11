# Positive spine automata, exact target languages, and a physical binary certificate

**2026-11-11. PROVED elementary reductions and explicit constructions; FINITE-EXACT controls; OPEN coverage.** No claim of novelty over all historical automata work, and no Chan–Pak resolution.

## 1. Inheritance and scope

Inherited frontier: Session260 `SYNTHESIS.md`, especially physical substitution with positive Stirling coefficients, the gap between congruence and equality, and the need for realization sidecars. Parent identifies THM-4940, THM-4952, THM-5070 as existing odd-cycle/ribbon/port-chain mechanisms: this note does not reclaim them. The closest newly used historical source is Session252's five-dimensional interface, whose ambient algebra need not be physically realizable. A hostile is count29: one-apex words cannot realize it, but two-apex words can. A corrected near miss is a fixed-dimensional matrix algebra mistaken for arbitrary realizability. The least-used sidecar here is a literal word producing the complete nonnegative transfer vector.

The principal output is **not** a fixed-ray membership algorithm. It is a positive finite automaton for every arbitrary spine word around a fixed apex core, and a finite linear test that certifies a whole physically realized binary continuation cone if a suitable seed is found.

## 2. Concrete tournament grammar

Fix a labelled tournament A on k vertices a_1,...,a_k. A word w over the alphabet {0,1}^k constructs T_A(w): the other vertices are the letters in order, every spine edge points forward, and a_i points to the vertex carrying letter b precisely when b_i=1. The core edges are those of A.

This grammar is exactly the tournaments with a specified k-vertex deletion leaving a transitive tournament. Union over the finitely many cores captures deletion number at most k for orders at least k; smaller tournaments can be handled separately.

For each permutation pi of the core, a Hamilton path traverses the core in that order. Between these vertices are k+1 possibly empty spine segments. Since all spine edges point forward, each segment is an increasing subsequence. Therefore a path is uniquely encoded by pi and an assignment of each spine vertex to one segment j in {0,...,k}.

For an internal nonempty segment j, its first vertex must receive an arc from pi_j and its last vertex must point to pi_(j+1). For segment0 there is only a last condition; segment k has only a first condition. An empty internal segment is allowed precisely when pi_j points to pi_(j+1). These conditions are necessary and sufficient, including consecutive empty segments.

## 3. A nonnegative integer weighted automaton

Retain pi and one status per segment. A status is empty, nonempty with last endpoint legal, or nonempty with last endpoint illegal. The last segment never needs the third status. Thus there are

    D = k! * 2 * 3^k

states (also valid for k=0). Initial weight1 is placed on the all-empty state for each pi. Reading b chooses exactly one segment j:

* If that segment is empty and j>0, require b_(pi_j)=1.
* Set its status to legal if j=k or b_(pi_(j+1))=0, and illegal otherwise.
* Keep every other status unchanged.

Different assignments reaching the same state add their integer weights. A state accepts if no last endpoint is illegal and each empty internal segment passes its core-edge condition. Empty first and last segments impose no condition.

The encoding in section2 is a weight-preserving bijection, so the accepted weight is exactly H(T_A(w)). Every transition matrix has nonnegative integer entries. The explicit compiler in `check_automaton.py` implements precisely these rules, without any odd-cycle identity or unproved tournament input. The initial all-empty states also correctly count paths in the bare core.

## 4. Exact all-length decision, not a finite word cutoff

Fix a target N>=0. Replace all scalar arithmetic by the finite semiring

    S_N={0,1,...,N,N+1},
    x plus y = min(N+1,x+y),
    x times y = min(N+1,xy).

The truncation map from nonnegative integers is a semiring homomorphism. Consequently the capped transfer vector after any word is exactly the coordinatewise truncation of its actual vector. An output equals N in the capped computation iff its actual value equals N. This remains true even when some large intermediate coordinate later dies: its contribution is either zero or still at least N+1. There is no cancellation.

**Theorem.** For fixed A,N, the words with H(T_A(w))=N form an effectively constructible regular language. Breadth-first search on capped vectors decides existence over **all word lengths**, returning a shortest witness or a finite complete rejecting state graph. There are at most (N+2)^D capped vectors; an accepted word, if any, has length at most (N+2)^D-1.

This upper bound is intentionally crude. It is a decision theorem, not an efficient algorithm claim or a proof that every allowable N is accepted. Modulo any positive integer, the same construction gives a finite residue automaton, but residue acceptance does not imply exact acceptance.

**Regular side conditions.** Intersect the capped automaton with any finite word automaton, retaining its state alongside the count vector. Thus any separately justified forbidden-subsequence, bounded-packing, endpoint, or phase constraint can be included without changing exactness. This is the interface to the parent's packing language; none of its structural theorems is needed here.

## 5. Reproducible exact controls

`python3 check_automaton.py` compares the automaton to independent subset-DP Hamilton counting on 5,336 specified cases with k=0,1,2,3 and total order<=6. It enumerates all core orientations, exhausts each word bank of size<=500, and otherwise samples exactly500 words using the recorded seeds. These are focused tests of the new compiler, not a broad audit of inherited work.

All-length capped searches in the committed program give:

* k=1,N=7: absent,39 reachable vectors.
* k=1,N=13: absent,62 vectors.
* k=1,N=21: absent,98 vectors.
* k=1,N=29: absent,127 vectors.
* k=1,N=43: absent,187 vectors.
* k=2, with a_1->a_2, N=7: absent,13,579 vectors.
* k=2,N=21: absent,108,618 vectors (the script uses a150,000-state budget and closes before it).
* k=2,N=29: word [1,2,0,2].
* k=2,N=43: word [3,1,2,0].

Letters are integer bitmasks with bit0 attached to a_1. The two-apex absence results recover consequences of the inherited global exceptions; they are controls demonstrating that the algorithm can actually terminate negatively. The 29 and43 witnesses show the extra core vertex is productive. No minimal tournament-order claim is made from these word searches.

## 6. Productive all-size constructor rays

Use a_1->a_2 and the same bitmask convention. The following are actual tournaments, with m repeated letters, not formal arithmetic:

    H(T_A(1,2,0,2^m)) = 15*2^m-1,                    m>=1;
    H(T_A(3,1,2,0^m)) = 12*3^m+8*2^m-9,             m>=1.

Thus the first family carries the physical operation H->2H+1 starting at29. The second couples genuine binary and ternary channels. Their m=0 values are13 and15, respectively; **the displayed formulas must not be extended to m=0**.

There are also periodic-block constructions with common characteristic polynomial

    (x-9)(x-4)^2(x-1)^2.

The three literal word families and their first five counts are

    (1,2)(0,2)^m:       5,29,187,1361,10823;
    (1,2,0,2)(1,0)^m:  29,231,1825,15179,131061;
    (1)(2,0)^m:         1,13,91,673,5383.

In each case m>=0, and subsequent terms satisfy

    h_(m+5)=19h_(m+4)-123h_(m+3)+337h_(m+2)-376h_(m+1)+144h_m.

**Proof certificate, not extrapolation.** The transfer representation has dimension36. For any proposed recurrence polynomial p, its residual sequence is another scalar sequence with the same36-dimensional transfer. By Cayley–Hamilton, vanishing on36 consecutive terms beginning at the claimed start proves all later residuals vanish. `certify_rays.py` checks these36 integer residuals exactly for each family and checks the initial conditions. The first two characteristic polynomials are (x-1)(x-2) and (x-1)(x-2)(x-3), with start1. This gives finite exact certificates of all-size identities; it does not rely on fitting a handful of terms.

The second-order and third-order rays are sparse. The repeated eigenvalues of the periodic families explicitly retain polynomial-in-m factors that a scalar multiplier heuristic would erase.

## 7. A finite certificate that WOULD prove large-count completeness

Write the positive compiler in column convention v(wb)=M_b v(w), with H=Lv. Adjoin a constant coordinate, so a word block u acts by B_u=diag(M_u,1). For candidate physical digit blocks u_+,u_-, begin with the rows

    r_+ = (L M_(u_+) - 2L, -1),
    r_- = (L M_(u_-) - 2L, +1).

Close their row span under right multiplication by B_(u_+) and B_(u_-). This process stabilizes in dimension at most D+1 and uses exact rational linear algebra. Call the resulting space R.

**Binary-cylinder criterion (PROVED).** A literal seed word w satisfies

    H(w followed by arbitrary digit blocks)

with the exact rules h->2h+1 and h->2h-1 at every step iff every row of R annihilates (v(w),1).

Indeed every reachable defect is one of the initial defect rows multiplied by a word in B_+,B_-. Their span is exactly R; evaluating at the physical seed gives the iff. This includes all continuations, not only immediate digit tests or an ambient affine line that happens to have the right scalar output.

**Conditional completeness theorem.** If for every odd h in [B,2B-1], with B odd and B>=3, one supplies a literal seed of count h and a certified binary cylinder (the blocks may differ between seeds), then every odd integer >=B is realized. The seed interval contains B,B+2,...,2B-1. From any odd N>=2B+1 choose the unique sign for which (N-sign)/2 is odd; this smaller parent is at least B. Iterate to a seed. The required digits are realized by that seed's two blocks.

For example B=23 requires the twelve seeds23,25,...,45. Adding literal seeds for allowable smaller values would establish the full Chan–Pak spectrum. **No such seed bank with certified cylinders is supplied.** Existing count witnesses do not automatically satisfy the rowspace condition.

`probe_binary_cylinders.py` tests all twelve ordered distinct one-letter pairs for k=2. The exact rowspace ranks are9,9,11,9,17,5,9,17,5,11,5,5 in lexicographic ordered-pair order. A separate finite physical seed scan is reported in `binary_cylinders.jsonl`. Ambient invariant-line equations can have solutions even when this physical scan finds none; their existence is not counted as progress to coverage.

## 8. What remains open

1. Find a physically reachable two-apex seed in a productive binary cylinder, possibly with blocks longer than one letter and with finite control states between different cylinders.
2. Replace a single return space by a finite family of exact invariant affine cones or linear slices, with transitions certified in the full positive representation.
3. Decide whether two-apex spectra are already all odd values other than7,21. The all-target conjecture does not follow from the per-target decision theorem.
4. Use the parent's packing-constraint automaton to eliminate impossible physical rowspace solutions before searching for seed words.

A finite graph of control states can generalize section7: each state carries its own row constraints and digit blocks; closure includes the destination state's constraints pulled back by each physical block. A finite certificate still requires literal seed reachability and a complete inverse-digit covering argument. Neither may be replaced by a numerical matrix orbit.
