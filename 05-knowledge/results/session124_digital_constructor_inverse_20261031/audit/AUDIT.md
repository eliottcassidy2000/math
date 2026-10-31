# Independent audit: digital singleton/triangle normal form

**Verdict: the exact identities, frozen fibers and complete finite non-3 inverse are mathematically sound within the stated actual-seed constructor language.** Completeness means a decision/reconstruction algorithm over every finite word; it does not mean every target has a positive answer.

## Independent reconstruction and checks

`independent_check.py` reconstructs literal paths and ordered two-path covers for every labeled seed tournament of orders1–5 and every cut. The33,866 seed/cut profiles give exactly the273 retained states. Every stored witness is separately reconstructed, and every h is verified positive and odd. Stored t is converted to g=h−t explicitly; it is not mistaken for g.

The checker then compares actual adjacency Hamiltonian-path DP with seed profiles, singleton extensions and all four triangle cuts; tests every singleton word through length6 at every seed; checks exact Fraction-valued radical-inverse errors; and compares inverse decisions with an independently generated **finite forward closure without a word-depth limit**. This closure uses F≤M and q<4M with M=63, not the unsound arbitrary-block q≤M pruning. The proof below justifies this comparison for non-3 targets. Returned small witnesses are again checked by actual adjacency DP.

`structure_check.py` independently computes forward singleton fibers to the full source boundary, compares the exact stored boundary lists, checks descendants of every boundary cylinder, verifies all four macro defects on an integer basis, and verifies the requisite modulo3 operator products on a basis. These are exact calculations, not numerical relaxations.

Receipts record exact counts, hashes and the four large frontier-repair controls. The large repairs are checked at their known normalized final Q and legal macro placements, not by an O(target) scan; their profiles and actual graph orders are reconstructed. The independent Python audit does not repeat the large-graph DP; the curated reproduction driver additionally recompiles the separate adjacency verifier and reruns all four full-matrix counts at O0/O2. Normal/optimized independent runs and copied producer replays are compared byte-for-byte. The curated `fibers.py` requires an explicit output path; regenerated data are compared against the retained JSON. Inputs are not overwritten.

## Global proof and quantifier audit

### Digital formula and radical inverse

The chronological address has reversed convention: its least significant bit is the earliest singleton, and address1 means singleton0. Summing F increments at singleton0 steps gives the claimed formula. Trailing singleton1s can be discarded for **F**, not for the whole state.

Unrolling the even recurrence gives the nonnegative inversion sum. Expanding it yields the pair-interaction formula. The exact run error is positive; the missing bit immediately above a maximal one-run implies R<1/2. Each run term is at most3/4, giving the stated O(log n) additive error. There is no uniform positive linear lower bound: all-one addresses have f=1. The curated proof denotes this explicitly defined radical inverse by rho_2(n), not a 2-adic valuation.

### Scalar inverse and the depth-free13 obstruction

For q>0, the nonnegative formula bounds B by floor((N−F)/q). For q=0, an actual seed must have an empty external cut: any outneighbor v supplies a cover with Q=(v) and P a Hamiltonian path of A−v. Hence F=h and g=0. Initial singleton0s are actual neutral loops.

Removing trailing binary ones from B leaves an even argument2n. If f(B)=m>1, the recurrence forces `(m−1)/2≤n≤m−1`. For m=13, the exact G(n) values for n=6,...,12 are `11,9,17,15,19,15,23`, none13. This is an all-word-length negative result, not a failed search to some depth. The two possible one-vertex marked seeds are independently checked by the scalar inverse; all transitive-spine/one-ear presentations begin from one of those seeds. This is not an unrestricted-tournament exclusion of13.

### Canonical q ancestry and actual source fibers

Because h is odd, parity identifies the last singleton uniquely: even q comes from q/2 through0; odd q≥h comes from(q−h)/2 through1. Every nonzero legal predecessor is strictly smaller. The only self-loop is q=0 through0, which is removed. Odd positive q<h are genuine ancestry roots, not missing predecessors.

After a singleton, g=h(q mod2). The potential increment at0 is exactly q+h(q mod2), and at1 it is0, proving the offset invariant. For a noncanonical actual seed, its two first-child offsets are different precisely as stated. Their **source nodes** must remain attached: an offset cannot be used at a node whose ancestry misses that source.

For a nonempty finite bank at fixed h, the largest source L is at least h. Thus every node above L has a unique first ancestor crossing the boundary; it lies in(L,2L+h]. No source lies strictly above L, so no new offset can enter below that boundary in its descendant tree. This proves the infinite frozen-fiber statement. It does not prevent overlap between different q fibers or offsets created by triangles.

### Triangle normalization and defects

After r triangles, h=3^r h0 and q is divisible by3^r. Therefore Q=q/3^r obeys the same singleton rule with h0, while C_j acts as Q→4Q+jh0. Its chronological binary shadow is the ordinary two-bit expansion of j. This statement concerns the normalized q address, not all state coordinates.

The four full-state operator differences from three times their two-singleton shadows were independently checked on all four basis vectors. They are linear identities, not realizability assertions about a scaled vector. Only execution of actual constructors from a retained actual seed certifies a tournament.

### Modulo3 exclusion and the safe final bound

C0/C3 kill all coordinates modulo3. After a nonconstant triangle, h and q vanish modulo3; after any intervening singleton, g also vanishes and stays zero. A subsequent triangle then kills F as well. Three adjacent nonconstant triangles also kill the entire vector. Hence non-3 output permits at most two triangles, adjacent if two, and neither can be constant. These are operator statements valid for every starting integer vector; no word-depth assumption enters.

Remove all trailing singleton1s. An empty remaining word gives the seed F. Otherwise its last block is S0,C1 orC2. The displayed positive differences give respectively q<2F, q<2F orq<4F at that **reduced final state**. They do not justify q≤N pruning in an arbitrary intermediate block word.

For the independent forward closure with cap M, q is nondecreasing under every physical block and F never decreases. Every non-3 target≤M therefore has a reduced witness all of whose intermediate q values are<4M and F values≤M. A seed itself is separately recorded even if its q exceeds the cap. The finite state space is bounded since h≤F≤M and0≤g≤h, with bounded q. This proves that the closure comparison has no hidden depth limit or invalid intermediate q≤M pruning.

### Completeness of the macro inverse

For every reduced non-3 witness, its r∈{0,1,2} and final Q lie in the enumerated interval. Expanding triangles into their shadows gives a genuine singleton ancestry from the seed q. The only information lost by reduced ancestry consists of initial0s at q0=0. All initial singleton0s are neutral; only the0 belonging to an initial C1=01 must be retained. One optional leading0 therefore suffices, even for two adjacent macros. The checker includes an actual case whose C1 representation requires this exception.

All remaining macro positions are exactly the listed one opposite-bit pair, or two adjacent opposite-bit pairs. Executing each candidate from the actual seed and checking F=N proves soundness. The preceding reduction proves necessity. It reconstructs an actual word with the target count, not necessarily the original unreduced word or its final q.

The four repaired counts103061,594793,774073,949777 have normalized Q values40836,265882,340300,394098. Each legal candidate decoding recovers its entire reduced actual state and satisfies q_final<4F_final.

### Complexity

The final-state loop is linear in the **numeric target** N, not in its binary encoding length. An ancestry has O(log N) bits; with at most two adjacent macros there are O(log N) placements, each requiring O(log N) execution work. A conservative fixed-seed bound is O(N log²(N+1)) arithmetic operations for this straightforward implementation, with additional integer bit costs. The seed-bank factor is finite. This is a pseudo-polynomial decision reduction, not a polynomial-time algorithm in log N or a universal positive-answer theorem.

## Publication hygiene

The originally audited unused `digital.invert_singleton` helper had an inaccurate reduced-word docstring and a singleton-only pruning rule. It was not used by any theorem or checker. The curated packet removes that helper, records the original and curated source hashes in `../CURATION.json`, and reruns both independent audits. The audited `inverse_addresses` and `inverse_non3` algorithms are unchanged.

No global coverage claim, prime-lane closure, unrestricted factor rule, or new forbidden tournament count is established.
