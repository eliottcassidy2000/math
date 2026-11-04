# Independent audit of named defect transport and prototype census

2026-11-04. **VERIFIED proof review and independent exact census.** Audited the sibling grammar/PROOF.md, search.cpp, search_allmarks.cpp and verify.py. No mathematical errors found.

The bad-join identity is correct for every tournament and physical mark: splitting each Hamiltonian path immediately after the mark gives exactly the joinable two-slot covers, leaving precisely J(w,v) for v->w. The proof includes the empty second slot correctly.

The birth recurrence is correct coefficientwise. The old-mark endpoint sum gains z_x(1+L)R_w. Each old in-neighbor start sum gains its path-through-x correction, and the new in-neighbor x contributes the constant c_w in K_w. Multiplication and z_x^2=0 leave exactly the two displayed terms. All RHS polynomials are old-state quantities. No bounded-state, native complement legality, or affine realizability follows merely from this transfer; the packet states these limitations correctly.

Replayed grammar/verify.py with normal Python and Python -O. The stdout files grammar_audit_normal.out and grammar_audit_optimized.out are identical. This replays the author's exhaustive small-case algebra tests and permutation checks rather than claiming those implementations are newly independent.

For a genuinely independent completeness check, audit_grammar.cpp uses **forward path extension** instead of the original backward last-vertex recurrence, computes every mark in one table, and uses **Floyd transitive closure** instead of the author's iterative reachability routine. It exhausts all8192 labelled two-birth extensions of each source and counts both strong and rejected weak matches. Results in grammar_independent_census.txt exactly match every cut/mark in the author's saved all-mark census:

- Both sources have7938 strong outputs.
- C2 has six fixed-mark +2 outputs, twelve all-mark +2 outputs, and no +4 output.
- C4 has no strong affine output for either digit or any mark.
- The two fixed-mark C4 +2 scalar matches are rejected weak outputs, at cuts(57,127) and(63,57). Their final or preceding newborn is a global sink respectively.

All twelve surviving cases also have strong first-birth intermediates: each first cut is neither empty nor full, and the old source is strong. This is a supplementary check on the literal two-birth interpretation, not an expansion to edge-rewrite programs.

The source scope is explicitly the Session180 six-vertex uniform prototype outputs C2 and C4. The newer actual170807 lineage certificates are neither searched nor excluded. Exactly two births with all old induced edges retained is the complete search universe; arbitrary markings expand it but do not change its birth count or permit old-edge rewrites.

Reproduction:

    g++ -O2 -std=c++17 audit_grammar.cpp -o audit_grammar
    ./audit_grammar
    python ../grammar/verify.py
    python -O ../grammar/verify.py

The count bounds at order8 are safe in uint64_t. The author's F<=8*8! is a valid loose bound; for a fixed marked vertex each permutation determines its unique split, giving the stronger F<=8! as well. No overflow concern affects this finite audit.
