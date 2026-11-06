# A unary-ray obstruction to selector-free rebuilding menus

Status: PROVED, independently audited. This is an obstruction to one proposed universal witness class, not to an eventual odd Hamiltonian-path spectrum threshold.

## Statement

Let T(w) consist of apices 0->1 and a transitive chain in word order. A letter a in {0,1,2,3} records which apices that chain vertex beats. Write H(w) for its number of directed Hamiltonian paths.

Fix any finite list of rewrites

    f_i(w) = p_i phi_i(w) q_i,

where p_i and q_i are arbitrary fixed finite words and each phi_i sends each letter to an arbitrary finite word, possibly empty. There is an integer N such that for every n>=N, no rewrite in the list satisfies

    H(f_i(1^n)) = 2 H(1^n) + 1.

Consequently no such finite list supplies the plus branch on every word with H(w)>=23. The prefix/suffix lengths need not be bounded uniformly across all possible menus; each particular menu only needs to be finite.

This does not exclude an appropriate restricted closed witness language, stateful transducers, input-dependent contexts, or freely located interior insertion. It does not exclude rebuilding a count using a different representative.

## Proof

The endpoint-forest representation used in the Session204 finite-menu bridge has 28 states. A state is an unordered multiset of endpoint pairs of directed path components, using endpoint labels 0,1,* and allowing at most one component starting at *. Ports internal to paths are absent from the endpoint list. Initially the nonzero states are {(0,1)} and {(0,0),(1,1)}, each with weight one. The terminal row selects single-component states.

When a new chain vertex z is added, choose an optional predecessor at an old component end and an optional successor at an old component start. Every bulk end is a legal predecessor; port i is a legal predecessor iff the new letter has bit i zero. Only a port whose bit is one can be a successor. Joining both ends of the same component is prohibited. Replace z by * in the resulting endpoints and discard states with more than one bulk start. Discarding is safe: subsequent chain vertices cannot enter an existing bulk-start component, so its number cannot decrease. Every retained path forest has a unique inverse on deleting the newest vertex. These rules give exact matrices A_a and H(w)=ell A_last ... A_first x.

The accompanying standalone exact certificate constructs the matrices from those rules and verifies, for every a,

    characteristic polynomial(A_a) = t^10 (t-1)^12 (t-2)^5 (t-3),
    rank(A_a) = rank(A_a^2) = 18.

Equality of ranks shows that the zero generalized eigenspace has no nilpotent block of size greater than one. Hence, for arbitrary row u and column v, and n>=1,

    u A_a^n v = P_1(n) + P_2(n) 2^n + P_3(n) 3^n,

for polynomials P_j over the rationals (some may vanish). This assertion applies to every fixed pair of contexts, without enumerating or bounding them.

The certificate also verifies the vector recurrence

    (A_1-I)(A_1-2I) A_1 x = 0

and scalar initial values H(1)=3 and H(11)=5. Therefore H(1^n)=2^n+1 for n>=1.

Consider one rewrite. If phi_i(1) is empty, its output is constant along the unary ray, and so can meet the strictly growing target 2^(n+1)+3 only finitely often. First suppose phi_i(1)=a is a single letter. Its defect is

    d_i(n) = H(p_i a^n q_i) - 2^(n+1) - 3.

For n>=1 this is an exponential polynomial with bases 1,2,3. A nonzero such expression has only finitely many positive integer zeros: choose the largest base with a nonzero polynomial coefficient, then the leading term of that polynomial eventually dominates all lower bases and its own lower-degree terms. Thus infinitely many zeros would force all coefficients to vanish, and d_i(n)=0 for every n>=1. Taking n=1 would give a literal tournament T(p_i a q_i) with H=7. This is impossible by the established all-order forbidden-seven theorem; use the repaired formulation in THM-029-alpha1-gap.md or THM-1370-h-spectrum-omits-7-21-all-n.md, not the historical false assertion that three pairwise-intersecting cycles always share a common vertex.

It remains to handle phi_i(1)=v of length k>=2. The exact certificate additionally supplies a common permutation of the 28 coordinates making every A_a lower triangular, with the same diagonal vector for each letter. Its entries are ten zeros, twelve ones, five twos, and one three. Therefore A_v, the chronological product for v, is triangular in that same order with characteristic polynomial

    t^10 (t-1)^12 (t-2^k)^5 (t-3^k).

Even without any assertion about its zero Jordan blocks, for all sufficiently large n the scalar H(p_i v^n q_i) is an exponential polynomial with bases 1,2^k,3^k. Its defect from 2^(n+1)+3 has a nonzero coefficient -2 on the distinct base 2. The same dominant-term argument shows it has only finitely many zeros. No repeated-word product rank assumption is needed.

Each fixed rewrite therefore succeeds for only finitely many unary exponents. The union of these finite sets over a finite menu is finite, proving the statement. Since H(1^n)>=23 for n>=5, this ray lies in the proposed selector-free class beyond its first few elements.

## Exact scope for the Session204 criterion

The Session204 finite-menu bridge remains logically valid as a sufficient criterion. This theorem proves that the criterion cannot be satisfied on L={w:H(w)>=23} by any finite menu of arbitrary word-valued morphisms with fixed end contexts, even when those contexts are arbitrarily long. Such a menu must either restrict the witness language, enrich the rewriting operations, or first apply normalizations outside this same class.

A finite collection of normalizations and subsequent rewrites all of the displayed form does not evade the theorem: composition is again of this form. This includes finite branches required to preserve H before applying the plus step. An unbounded normalization procedure is outside the theorem.

## Reproduction

Run `python3 verify_unary_barrier.py` (Python 3 and SymPy). The verifier checks the exact characteristic polynomials, ranks, common triangular order and diagonal, unary recurrence, and initial values. It independently counts all 341 words of chain length at most four by literal permutation enumeration. `RESULT.json` is the author run. The forbidden-seven theorem is a mathematical dependency, not a claim inferred by finite enumeration.
