# Independent review of information, Hodge normalization and reflection defect

2026-10-26. **ACCEPTED within the stated scopes.** Reviewer: recurrence/enumeration agent, separate from this lane's producer. No new open-problem implication is inferred.

I read `RESULT.md`, `check.py`, the repaired **THM-895**, `01-canon/theorems/THM-895-helmholtz-identity-perron-census.md`, its legacy checker, and the actual reflection coordinates in `unit_distance/RESULT.md` §4. I replayed the final information checker normally and under `-O`, including its full involution-defect polynomial gate. Both review receipts have SHA256 `d18ef3ced44b9ec8b0db40326ac694fd107e6aeef7f287c963e6a0962b18867a`. The legacy script also ran successfully; its output is `review-legacy.txt`.

## Proof and quantifier checks

1. **Chart fibers.** Equal positional words are exactly automorphism-related Hamiltonian paths. Their stabilizers are trivial because a complete ordered path lists all vertices. Thus the fiber is H/Aut, not H; converse preserves that number, so a non-self-converse merged pair has twice the fiber. The reciprocal-fiber count is exact but already requires the fibers. This makes no compressed-recursion or connectivity inference.
2. **Entropy.** The displayed Shannon identity is for uniform chart words. The conditional value at a particular class C is `H(W | C=c)=log2 h(c)`; averaging gives the stated entropy subtraction. Uniform labeled tournaments and uniform classes use the other displayed distributions. The score-information lower bound explicitly uses **uniform classes**, where `H(C)=log2 A_n`. The multiset bound `binom(2n-1,n)` and `A_n>=2^binom(n,2)/n!` give the claimed lower bound. Adding the positive integers H and Aut costs at most `2 log2(n!)` more bits. The order-32 integer test is sufficient; no claim of a sharp first order is made. These arguments obstruct reconstruction from the listed statistics, not short aggregate counting formulas.
3. **Hodge.** For n>=1, the incidence normal equation is `(nI-J)phi=d`. Since sum(d)=0, the unique zero-mean solution is d/n. Its squared norm is `sum(d_i^2)/n`; the residual has zero divergence and is orthogonal to it. K2 already rejects the old factor 1/2. The cyclic-triangle identity gives the residual's nonzero transitive baseline exactly. Filling all triangles is essential to the zero-harmonic-space statement. A one-dimensional complete graph has cycles. The claim concerns a unit flow, not merely the signs of arbitrary potential differences.
4. **Legacy numerical scope.** The legacy checker now tests integral residual divergence, orthogonality and energy, beyond its floating energy check. It reproduces Perron correlations0.9313/0.9218 and the stated level spreads. Those correlations remain numerical census evidence. The exact nilpotence/acyclicity proof of spectral-radius-zero iff transitive is separate. I raised the ambiguous historical phrase that lambda alone “refines” x; the producer repaired it to the joint invariant `(x,lambda)`. The explicit counterexample `C3⊕C3` versus `C3⊕TT3` has lambda1 but x54 versus62.
5. **Reflection.** Conjugation of the displayed hexagon fixes z0,z3 and exchanges z1,z5 and z2,z4, so it is exactly two transpositions after relabeling. For a general involution with r transpositions and f fixed vertices, r arcs are forced to disagree. The remaining pair orbits consist of `binom(f,2)` pointwise fixed pairs and `b=(binom(n,2)-r-binom(f,2))/2` two-element pair orbits. Their independent choices give `2^(r+binom(f,2)+b) t^r (1+t^2)^b`. Thus the n6 spectrum is exactly `512 t^2(1+t^2)^6`, and minimum2. This is an obstruction to a tournament relation preserved by that same permutation, not a unit-distance chromatic obstruction or an arithmetic rank.

## An independent structural separation of the order-five twins

For each vertex record `(outdegree, number of cyclic triangles containing it)`. This multiset is invariant under isomorphism. Direct extraction from upper-pair codes8 and10 gives

```
code8:  (1,1), (1,3), (2,1), (3,1), (3,3)
code10: (1,1), (1,2), (2,3), (3,1), (3,2).
```

This proves their nonisomorphism without relying on a canonical-label oracle. Their common sorted scores, H=9, Aut=1 and three cyclic triangles are independently consistent with the enumeration lane's subset-DP and automorphism results. The smallest counterexample is not “cycles versus no cycles”; it is different incidence of the same number of cycles.

## Reproduction and remaining limits

Run `python check.py` and `python -O check.py` in this folder; the review saved their exact output separately as `review-replay-normal.json` and `review-replay-optimized.json`. Run the legacy script from the worktree root. The finite universe is all labeled tournaments through five, plus all six-vertex orientations for the reflection histogram, not all-order numerical evidence. The algebraic identities are proved independently of those finite universes. No claim about the single arithmetic source7 follows from the information or reflection statements.
