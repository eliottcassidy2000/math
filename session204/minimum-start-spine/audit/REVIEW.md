# Independent audit: minimum-entry ancestry spine

Status: PASS. Independently reviewed the complete author manuscript `session204-minimum-start-spine/PROOF.md`, including its product-charge refinement, four-exception bound, five-stretch bound, and total tail-block bound. No mathematical correction required.

## Dependency and derivation

Use the independently audited bounded-entry SCC decomposition, available locally at `session203-finite-profile-interface/dependencies/BOUNDED_ENTRY_PROOF.md`, Section 2.

Let (T,v) be strong and non-singleton, with positive integer start count k. Write the ordered SCCs of T-v as C,D1,...,Dt, put A=N_T^+(v) intersect C, m=|A|, and q=product H(Di). Then A is nonempty, C is strong (possibly a singleton), and

    k = q sum_{a in A} start_a(C).

Every entry has positive integer start count. Choose a minimizing entry a and set k'=start_a(C). Therefore

    q m k' <= k,
    1 <= k' <= floor(k/(q m)).

Repeat inside (C,a), stopping at a singleton. Order strictly decreases at every step, so the process terminates even through arbitrarily many equal-budget steps.

For successive selected steps i, multiplication gives

    product_i(q_i m_i) <= k_initial / k_terminal = k_initial,

because the terminal singleton has start count one. Define a structural exceptional junction to mean m_i>1 or q_i>1. Every such factor q_i m_i is at least two. Thus their number is at most floor(log2 k_initial), hence at most four for k_initial<=26.

In fact q_i is odd. If q_i>1, at least one nontrivial strong tail component occurs, so q_i>=3. This strengthens individual drops but is not needed for the four-junction bound.

## Classification of unexceptional steps

If q=m=1, every tail SCC is a singleton, because a nontrivial strong tournament has at least three Hamilton paths. The entry set is {a}. Construct T literally from (C,a) by introducing v beating only a and marking v (J), then introducing the singleton tail vertices in order. Each tail vertex loses to all earlier tail/core vertices except possibly v: it is an O move if v beats it and a Z move if it beats v. Thus this selected SCC step is exactly a J followed by O/Z graph word.

When C is a singleton, J first creates a rooted edge. This boundary is valid as a graph construction, but nontrivial-input transfer identities must not be applied to its singleton input. The final strongness conditions are those already imposed by the original strong T.

## Scope restrictions

1. The bound counts exceptions on this minimum-entry selected ancestry spine. It is not a bound on every possible ancestry choice or on the total branching/marked reconstruction complexity.
2. An exceptional junction is defined by its chosen SCC data. The theorem does not say that its output graph/profile has no alternative O/Z/J representation.
3. The other entry vertices of C remain part of the actual graph. Knowing only (C,a)'s scalar start count, or even its unmarked rooted profile, does not automatically retain the identities and simultaneous responses needed to reconstruct the original multi-entry junction.
4. The initial assertion is for strong rooted sources. For a general positive-start source, first isolate its initial SCC; the weak-tail factor may be handled separately, and should not silently be included in the strong-step count.
5. No numerical B-occupancy, finite binary automaton, F(8,137), or eventual Hamilton-path-count threshold follows without further construction arguments.

No exhaustive computation is necessary for this corollary: it is an exact inequality and telescoping argument from the existing all-order SCC identity.

## Review of additional author refinements

The disjoint stage classes E_t (q>1) and E_m (q=1,m>1) obey 3^E_t 2^E_m<=k_initial. Thus E_t<=2 for k_initial<=26. Every nontrivial tail block contributes a factor at least three to product q_i, so the total number of these blocks on the selected spine is at most floor(log3 k_initial), also at most two at budget 26. The tail blocks are pairwise vertex-disjoint because later cores lie inside earlier C components. At most E exceptional steps divide the ordered spine into at most E+1 (possibly empty) same-budget stretches. All these refinements are valid.
