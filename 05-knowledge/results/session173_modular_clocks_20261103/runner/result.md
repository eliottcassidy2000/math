# Composite shift masks, cycle defects, and a constructive runner certificate

## External input and conventions (source checked; global computation not rerun)

Jaan Allikvere, *Fourteen and fifteen lonely runners*, arXiv:2609.02604v2 (24 September 2026), https://arxiv.org/html/2609.02604v2, Theorem 1.1, establishes LRC(13) and LRC(14) in its convention: k nonzero relative speeds, k+1 total runners, target 1/(k+1). This includes both 14 and 15 total runners; the repository sometimes calls the 13-speed problem “LRC14.” The paper is a computer-assisted preprint, not a journal-publication claim. Its stated theorem, conventions, and Lemmas 4.5–4.6 were checked against the linked source; its global computation was not independently rerun.

The useful imported mechanisms are divisor-level lifting (Lemma 4.5) and few-exception shifting (Lemma 4.6). The latter fixes a common-divisor core while moving exceptional speeds around a finite grid, then bounds the union of forbidden shifts. Its stated prime-modulus, fourteen-speed criterion is e ceil(2d/15)<d. Its finite proof combines modular searches and a product bound; we do not repeat that computation or claim to extend its global range.

## Prior project results consulted

- THM-760: a scaled core plus one coprime exception, with an explicit one-sheet dodge. It explicitly leaves two or more exceptional classes untreated.
- THM-373: the runner-phase comparator is a finite tournament walk between rational walls.
- THM-639: Hamiltonian path frames and the distinction between a structural invariant and a lossy tournament quotient.
- Filenames LEM-044, LEM-045, THM-1217 show Hunter/tree overlap machinery already exists in the project. The contribution here is its concrete composite shift-mask specialization and exact cycle defect, not a claim to invent the general union bound.

## Proposition: composite-modulus overlap-aware sheet dodge

Let d>=2, 0<delta<=1/2, C a finite set of integer core speeds divisible by d, and E=(w_1,...,w_e), e>=1, exceptional integer speeds. Suppose t0 satisfies ||v t0||>=delta for every v in C. Define

B_i={q in Z/dZ: ||w_i(t0+q/d)||<delta}.

For any tree T on the e exceptional indices, set

U_T=sum_i |B_i| - sum_{ij in T}|B_i intersect B_j|.

Then at least max(0,d-U_T) of the d shifted times are delta-witnesses for C union E. In particular U_T<d supplies an explicit witness, found by scanning q. The best bound of this form uses a maximum-weight spanning tree with weights |B_i intersect B_j|.

Writing g_i=gcd(d,w_i), one always has

|B_i| <= g_i ceil(2 delta d/g_i).

Thus a phase-independent sufficient condition is sum_i g_i ceil(2 delta d/g_i)<d. This is valid for composite d and noncoprime exceptions; it specializes to the external paper's counting mechanism when d is prime and delta=1/15. It is not necessary: actual masks and their overlaps can be far smaller.

### Proof

Adding q/d fixes every core phase modulo one. For a fixed q, let S(q)={i:q in B_i}. The subgraph T[S(q)] is a forest. If S(q) is nonempty, |S(q)|-|E(T[S(q)])| is its number of connected components, at least one. If empty, both sides are zero. Summing over q gives |union B_i|<=U_T. Exactly d-|union B_i| shifts avoid every exception; the core already remains safe.

The exceptional phases form d/g_i equally spaced circular points, each repeated g_i times. An open circular arc of length 2delta contains at most ceil(2delta d/g_i) such points. The bad arc is open because equality at delta is allowed. This proves the count bound, including its endpoint convention. For delta=1/2 the bound is still valid, though it may be loose.

### Exact defect

The gap U_T-|union B_i| equals

sum_{q:S(q) nonempty}(components(T[S(q)])-1).

Consequently the certificate is exact precisely when each nonempty bad-index set induces a connected subtree. This is a testable structural condition, not a heuristic independence assumption.

## Cyclic upgrade, directly using an ordered cycle

For e>=3 and cyclic order 1,...,e, put B_{e+1}=B_1 and

U_cycle=sum_i|B_i|-sum_i|B_i intersect B_{i+1}|+|intersection_i B_i|.

Then |union B_i|<=U_cycle. Its gap is the sum, over shifts for which S(q) is nonempty and proper, of (number of cyclic runs of S(q))-1. The all-covered shifts are handled exactly by the full-intersection correction. This follows pointwise: on a proper subset of a cycle, selected vertices minus selected edges equals the number of runs; on the full cycle it is zero, so the correction restores one.

Deleting one edge turns the cycle into a tree. Since the full intersection has size at most every edge intersection, U_cycle is at least as strong as the best bound obtained solely by deleting one edge from this fixed cycle. A different spanning tree may be better still.

## Exact example in the specified order (4,2,5,3,7,13)

Take d=9, C=9*{1,2,...,8}, E=(4,2,5,3,7,13), delta=1/15, t0=1/81. There are fourteen moving speeds. The core margin is 1/9. The six masks, in exactly that cyclic order, are

{0,2}, {0}, {0,7}, {0,3,6}, {5}, {2,4}.

Their sizes sum to 11. Consecutive intersection sizes are 1,1,1,0,0,1 and the total intersection is empty, giving U_cycle=7<9. Their union is exactly {0,2,3,4,5,6,7}. Thus q=1 and q=8 give witnesses 10/81 and 73/81. This is constructive at the standard 15-runner threshold; it does not depend on accepting the global LRC theorem.

The same cycle at d=8 with C=8*{1,...,8}, t0=1/72 has masks {0,2,4,6},{0,4},{3},{0},{1},{3}. The prescribed cycle bound is 8 and is inconclusive, but a maximum spanning tree gives 6, exactly the union size, and witnesses exist. This is an actionable distinction between keeping a specified cycle frame and optimizing the overlap frame.

## Clocks {1,3,5,7}, 44, 527, and the 5+8=13 motif

With d=8, C=8*{1,...,7}, E={1,3,5,7,13,44,527}, t0=1/64, the masks are {0},{0},{3},{1},{6},empty,{2}. The tree bound is 5, yielding q=4,5,7 and witnesses 33/64,41/64,57/64. The phase-independent gcd bound is 16, so actual clock geometry provides substantial information.

For d=9, add exception 8 and use C=9*{1,...,6}, t0=1/63. The eight masks are {0},{0,3,6},{7},{5},{4},{6},{3},{1}; their size sum is 10, the tree bound and union size are 7, and witnesses are 5/21 and 19/21. These examples are constructed families, not evidence that 44 or 527 are universally special.

A useful exact interpretation of 5+8=13 is available on the 8-sheet clock. At t0=1/64, their q-dependent steps agree because 13=5+8, but their initial phases differ by 8t0=1/8. Hence B_13 is a translate of B_5: B_5={3} and B_13={6}. They have the same shape but are not the same forbidden event. More generally, if w'=w+ad and gcd(w,d)=1 and ad t0 is an integral multiple h/d modulo one, then B_{w'}=B_w-h*w^{-1} in Z/dZ. This is a conditional exact translation law; the congruence alone is insufficient without the phase condition.

The same experiment includes d=44 and d=527=17*31. No primality is needed by the proposition. For the cycle family, the tree/union counts are 26/26 at d=44 and 329/321 at d=527; a nonzero defect is permitted and measured.

## Spectral bridge and boundary

For f_i=1_{B_i} and unnormalized DFT hat f_i(r)=sum_q f_i(q) exp(-2 pi i rq/d), Parseval gives

|B_i intersect B_j|=(1/d) sum_r hat f_i(r) conjugate(hat f_j(r)).

Thus the mask Gram matrix supplies exactly the weighted overlaps used by the tree certificate; a Hamiltonian path in any chosen tournament on the exceptions is one admissible tree. Retaining the masks or their Fourier vectors gives a lawful bridge from modular clocks to a weighted path certificate. The eigenvalue multiset of an unrelated tournament adjacency matrix does not by itself recover these intersections. No spectral-only implication is asserted.

## Verification and scope

`test.py` uses Python Fraction throughout, checks every claimed witness against every speed, verifies the gcd capacities, and checks 345 systematic parameter families. Of these, 271 have a successful tree certificate even though the phase-independent sum bound fails. These are explicitly constructed finite families, not a classification or global LRC proof. `cycle_tests.py` additionally checks all set families on a four-point universe with four event masks, verifies the pointwise cycle defect formula, and records the prescribed-cycle examples.

## Independent audit

The explicit-check `independent_audit.py` passes normally and under `python -O`: 8,190 membership patterns for tree/cycle bounds, 22,560 signed-speed and endpoint capacity cases, four exact fourteen-distinct-speed examples, and 476 conditional translation checks. Author tests also pass in both modes with identical receipts. The proofs were checked independently; finite controls support, but do not replace, the proofs. The source comparison does not certify the external global computer-assisted theorem.
