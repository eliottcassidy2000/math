# Incoming cut-arc law and a canonical digit-one readiness bridge

**Date:** 2026-11-04. **Scope:** targeted comparison, not a replay or recertification of the incoming large censuses. This note proves the transport argument below using the explicitly pinned incoming structural theorem. No novelty against the literature or earlier repository work is claimed.

## 1. Pinned inputs and their precise status

1. [S74 audit, section A](https://github.com/eliottcassidy2000/math/blob/a0518f53cd622f71d2e268effa09519cce3bdc15/04-computation/experiments/uh_one_arc_cut_20261104_audit_S74/AUDIT.md). Local copy `S74_AUDIT.md`.
   - It proves, for every word of the recorded Shelor grammar, `k = 3 + #U - #D`, with parent offsets `j=0,1,2` giving changes `+1,0,-1` in the number of cut arcs.
   - Its all-order interpretation as the distribution of **isomorphism classes of uniquely Hamiltonian tournaments** uses the classical completeness and injectivity of that grammar. The audit explicitly labels this dependence; the primary Shelor text was inaccessible to it.
   - With `m=n-4`, its class formula is `T(n,k)=sum_(d>=0) binom(m,k-3+2d) binom(k-2+d,d)` and its series is `sum_(n>=4,k) T(n,k) X^(n-4) Z^k = Z^2 (X+(1-X)Z)/(1-2X-X(1-X)Z)`.
   - Keep the small-order ranges: `T(n,2)=2^(n-5)` for `n>=5`; `T(n,3)=n 2^(n-7)` for `n>=6`; `T(n,4)=(n^2+n-8)2^(n-10)` for `n>=7`. At order four only `T(4,3)=1` is nonzero. The binomial top-diagonal pattern stops after its first three diagonals.
   - The audit reports comparison with Shelor-free class sets through order 16 and 9,227,464 generated-word checks over orders 4–20. Those large computations were **not rerun here**.

2. [THM-4868](https://github.com/eliottcassidy2000/math/blob/c2fbfacee956e64bc1320413614c2622195509fb/01-canon/theorems/THM-4868-the-noncritical-run-of-a-uniquely-hamiltonian-tournament-local-roundness-one-run-cut-arcs-plus-noncritical-equals-n-plus-1.md) and [full S10 proof, sections 1–2](https://github.com/eliottcassidy2000/math/blob/c2fbfacee956e64bc1320413614c2622195509fb/05-knowledge/results/uniquely_hamiltonian_noncritical_run_20261104.md). Local copies `THM4868.md`, `S10_NOTE.md`.
   - These contain actual all-order proofs, not just a table fit: the noncritical vertices of a uniquely Hamiltonian tournament of order at least four form one run `y -> ... -> x`; `y` is its unique noncritical sole entry and `x` its unique noncritical sole exit; `k+r=n+1`; the admissible insertion parents and transitions are as in section 3 below.
   - This is a second, Shelor-free proof of the grammar and cut-arc recurrence, using Moon's strong-deletion theorem and the elementary insertion lemma rather than assuming Shelor's characterization. The pinned `c2fbface` change corrects priority, not the underlying proof mechanism.
   - The source reports the full census through order 18 and parent-transition checks through order 14. Its section 9 is still the placeholder `(Filled after the independent audit.)`. Thus do **not** call S10 independently audited on this evidence. This note does not independently certify its entire theorem or its finite census.
   - Primary literature bodies were not accessible to S10. Classical priority is source-reported; this note makes no independent literature-priority claim.

3. The already reviewed Session189 [correspondence repair](https://github.com/eliottcassidy2000/math/blob/80147be9d0693a3b1b0436d3fe4ac3781d492bd6/05-knowledge/results/session189_response_frontier_20261104/readiness_repair/REPAIR.md) and its [independent review](https://github.com/eliottcassidy2000/math/blob/80147be9d0693a3b1b0436d3fe4ac3781d492bd6/05-knowledge/results/session189_response_frontier_20261104/readiness_audit/REVIEW.md) are pinned at checkpoint `80147be9`. Their local copies `../session189-returns/incoming_compare/REPAIR.md` and `../session189-root/READY_REVIEW.md` match the immutable Git blobs. Its exact marked image requires a **strong deleted core**; choosing arbitrary sole-exit marks is invalid. The essential marked identity is re-proved below, so the bridge does not rely on an uncorrected marked-count assertion.

The complete immutable file identities are in `SOURCES.json`.

## 2. Definitions and the canonical bijection used here

`H(R)` counts all directed Hamiltonian paths. `hc(Y)` counts Hamiltonian cycles modulo cyclic rotation. For a tournament `U`, `e_p(U)` and `s_p(U)` count Hamiltonian paths ending and starting at `p`.

A **strong digit-one ready tournament** `R` has a sole exit `x -> q` and `e_q(R-x)=1`. Necessarily `R-x` is strong: strongness puts `q` in the first strong component of `R-x`, while a spanning path ending at `q` puts it in the last.

For a marked sole-exit tournament `(R,x)`, let `Phi(R,x)` reverse exactly the arcs with both ends in `R-x`, leaving every arc incident to `x` unchanged. It is a named involution. A Hamiltonian cycle through `x` corresponds to a Hamiltonian path of the core starting at `q`, so

    hc(Phi(R,x)) = s_q((R-x)^op) = e_q(R-x).

Conversely, if `Y` is uniquely Hamiltonian, `x` is a sole exit and `Y-x` is strong, then `R=Phi(Y,x)` is strong and digit one: the reversed core is strong, `x` reaches it via `q`, and it reaches `x` through any other core vertex. Here `n>=4`, so such other vertices exist. Also `e_q(R-x)=s_q(Y-x)=hc(Y)=1`.

There can be at most one sole exit whose deletion is strong. Indeed if two sole exits satisfy `a->b`, deleting `b` leaves `a` a sink. Every positive-ready mark has strong deleted core by the preceding argument. By THM-4868(C), every uniquely Hamiltonian `Y` of order `n>=4` has exactly one noncritical sole exit `x`. Therefore choosing this canonical `x` gives a relabeling-compatible bijection between unmarked UH classes and strong digit-one ready classes, in both directions. The order-three directed triangle is excluded: it has no strong two-vertex deletion.

## 3. Canonical ready heredity and the exact 2/3-child rule

**Theorem (transport of THM-4868(C,E,G,H)).** The strong digit-one ready classes of orders `n>=4` have the following canonical rooted tree. Its unique order-four seed is

    R4 outgoing rows = (12,5,8,2), canonical exit x=2, q=3.

For a ready parent `R`, form `U=Phi(R,x)`. Let the unique Hamiltonian cycle of `U` have noncritical run from `y` to `x`. On this cycle write `y-` and `y--` for the preceding vertices. The allowed parents are exactly

    p = y, y-, and additionally y-- if y- has indegree one in U.

For each such `p`, append a vertex `w` to `R` with

    N+(w) among the old vertices = {x,p}.

Every other old vertex beats `w`. The result is again strong digit one with the same canonical sole exit `x`. Every ready class of order `n+1>=5` occurs exactly once as one of these children, up to isomorphism. In particular every such class has a canonical induced ready parent and retains its sole-exit mark under that deletion.

**Proof.** THM-4868(G) says that every UH child is uniquely reconstructed as `U^p`, obtained by adding `w` with `p->w` and `w->u` for every other old `u`, for one of the displayed admissible `p`. Its new run starts at `w`.

None of the allowed `p` equals `x`: `y!=x` because the run has length at least two; `y-` is critical; and if `y-` is a sole entry, `y--` is critical by the source's Lemma 2. Since `U-x` is strong, `(U^p)-x` is strong: `w` has in-neighbor `p` there and an out-neighbor other than `p`, because `|U|>=4`. Thus `x` remains the unique noncritical sole exit of `U^p`. This also follows directly from (H).

Now apply `Phi` at that same physical `x`. On the old vertices it is exactly `Phi(U,x)=R`. The edge `w->x` is retained. For each old core vertex other than `p`, the edge `w->u` reverses to `u->w`; the edge `p->w` reverses to `w->p`. Thus the ready child has exactly the stated two old out-neighbors `{x,p}`. Section 2 proves strongness and digit one.

For the converse, start from a ready child, use its unique positive-ready mark `x`, and apply `Phi`. THM-4868(C,G) gives its canonical UH deletion of `w`, with the preceding argument identifying the same `x` on its parent. Applying `Phi` before or after deleting `w` therefore yields the same induced ready parent. Uniqueness of the canonical UH decomposition and THM-4868(E)'s trivial automorphism group exclude duplicate child classes from different marked admissible pairs. The source seed is `S4=(6,12,8,1)` with cycle `(0,1,2,3)` and canonical exit `2`; applying `Phi` gives the displayed ready seed. QED.

This proves all-order **digit-one** induced-deletion heredity from the incoming structural theorem. It neither assumes nor proves digit-two heredity, and does not imply fixed-order H-neutral rearming or global spectrum coverage.

## 4. A concrete path-count update, and the necessary richer state

For a tournament `R` on `V` define

    F_z(R) = sum_(A disjoint union B = V, z in B) H(R[A]) s_z(R[B]),
    G_(q,z)(R) = sum_(A disjoint union B = V, q in A, z in B) e_q(R[A]) s_z(R[B]),

with ordered parts and `H(empty)=1`. For the two-outneighbor child above, partition a Hamiltonian path by the position of `w`. If `w` is last, its predecessor is neither `x` nor `p`. Otherwise its successor is `x` or `p`; a nonempty prefix cannot end at either of these vertices. This gives the exact identity

    H(R+w) = H(R)-e_x(R)-e_p(R)+F_x(R)+F_p(R)-G_(x,p)(R)-G_(p,x)(R).

At the ready mark, `F_x=H+1`. For a direct derivation, put `Q=R-x`, with `x->q` its only exit. Partitioning a path by the position of `x` gives `H(R)=H(Q)-e_q(Q)+F_q(Q)`. In an ordered two-path cover whose second path starts at `x`, that path is either the singleton `x` or starts `x,q`; hence `F_x(R)=H(Q)+F_q(Q)`. Their difference is `e_q(Q)=1`. Thus the same law is

    H(R+w) = 2H(R)+1-e_x(R)-e_p(R)+F_p(R)-G_(x,p)(R)-G_(p,x)(R).

This is a one-step observable update, **not** a proof that these observables form a closed state under repeated growth. It suggests transporting endpoint and two-path-cover responses along the canonical grammar, rather than extrapolating spectrum density from its unweighted generating function.

**Exact obstruction to an unweighted-state shortcut.** The following two UH tournaments both have order six, run length three, four cut arcs, two sole entries, two sole exits, and even the same `H(U)=15`:

    U1 = (6,60,8,1,13,29), cycle (0,1,5,4,2,3), H(U1)=15.
    U2 = (22,12,8,1,46,15), cycle (0,4,5,1,2,3), H(U2)=15.

Their common named canonical sole exit is `x=2`, with `q=3`. But their ready images are

    R1 = (60,5,8,50,38,6), H(R1)=33, sigma=1;
    R2 = (44,53,8,50,5,20), H(R2)=29, sigma=1.

Rows encode the outgoing-neighbor bit mask on names `0,...,5`. Thus `(n,r,k,#entries,#exits,H(U))` does not determine `H(Phi(U,x))`. In particular the cut-arc/run generating function is not a ready-spectrum formula.

## 5. Small reproducible bridge check and boundaries

`check_bridge.py` is independent standard-library code. Run:

    python3 check_bridge.py > check_bridge.json
    python3 -O check_bridge.py > check_bridge.optimized.json
    cmp check_bridge.json check_bridge.optimized.json

The receipts agree byte-for-byte. It checks 88 generated UH classes at orders 4–8 and 87 canonical insertions. It uses exact subset path DP, an independent cycle count, reachability for strong deletions and cut arcs, and rotation-canonical UH codes. It verifies the canonical Phi image, persistent exit, induced ready parent, admissible parent set, run transitions, the two-outneighbor path update, and the two explicit witnesses above. This is not an independent all-tournament census or a rerun of the reported order-18/20 calculations; completeness beyond the checked generated family comes from the stated source theorem, not from these tests.

**Useful next step:** use this canonical, isomorph-free digit-one generator to organize exact ready-spectrum computations with the named endpoint/cover data needed for section 4. Separately address digit two and weighted-state closure. Nothing here improves the prime-nine named-response collision frontier, proves arbitrary ambient response injectivity, closes reusable two-branch returns, or establishes general Hamilton-spectrum coverage.
