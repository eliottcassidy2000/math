# Exact edge interfaces, orphan paths, and realized response components

Session 175. **Proved identities and finite-exact constructions. No unbounded spectrum coverage or historical novelty claim.**

## 1. Correct target and inheritance

Ordinary Hamilton-path counts 7 and 21 are impossible at every order. The open question is completeness of all remaining positive odd integers. This packet does not search for the two impossible counts. The two infinite strong-witness lanes are primes p other than 7, and 7p for odd primes p other than 3; SCC products and carry witnesses then suffice. These reductions, the orphan warning, and the unit-component criterion are inherited from:

- [Repository spectrum frontier, pinned a55a46e2](https://github.com/eliottcassidy2000/math/blob/a55a46e2a9fd3880980b26149c0c6680f3ec8909/05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md), §§3–6, and its THM-4094/THM-4118 references.
- [S174 wrapper obstruction](https://github.com/eliottcassidy2000/math/blob/a55a46e2a9fd3880980b26149c0c6680f3ec8909/05-knowledge/results/session174_palindrome_symmetry_20261103/synthesis/NOTE.md), §§3,5. Its suggestion of constructing ordinary H=7 or 21 is superseded by the earlier global exclusion, not an open construction target.
- [Busch 2006, Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) gives the strong minimum and finite cutoff underlying the exclusion.
- [Chan–Pak, §6.5, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf) states the completeness conjecture. [Grinberg's April 27, 2026 notes, Remark 1.1.11](https://darijgrinberg.gitlab.io/t/22s/lec12.pdf) still describes the higher-value problem as unknown.

The session contribution here is an explicit endpoint-factor implementation of local responses, literal nonuniform-interface witnesses, and completely enumerated response components. These finite counts were already within known global coverage; they are not new attained values globally.

## 2. Endpoint-refined edge-flip identity

For every nonempty vertex set S define E(S,u) as the number of Hamilton paths of T[S] ending at u, and B(S,u) as the number beginning at u. Both are zero when u is not in S. Singleton counts are one. The usual endpoint recurrence is

    E(S,u) = sum_{w in S\{u}, w→u} E(S\{u},w).

B is computed by applying E to the reversed tournament. Fix distinct u,v and put R=V\{u,v}. Define

    C(u,v) = sum_{S subset R} E(S union {u},u) B((R\S) union {v},v).

This definition uses only induced subtournaments missing at least one of u,v, and therefore is unchanged by flipping their mutual edge.

**Proposition.** If T has u→v and T' reverses only this edge, then

    H(T') − H(T) = C(v,u) − C(u,v).

**Proof.** Paths avoiding the mutual edge are the same in T and T'. Every path of T using u→v splits uniquely at that edge into a prefix ending at u and a suffix beginning at v. Their other vertices partition R, giving exactly C(u,v) choices. The corresponding changed paths of T' number C(v,u). All factors are genuine endpoint counts; there is no formal-profile realizability assumption. Subtract. ∎

The formula can be evaluated for every possible edge from one pair of endpoint tables. Given a neighboring interface state, all tables are recomputed; one must not reuse them after other edges change.

## 3. Exact orphan contribution under vertex insertion

Let x be new, with I the old vertices beating x and O the old vertices beaten by x. For each old Hamilton path P, let a_x(P) be the number of legal insertion positions. Every new Hamilton path either deletes to an old Hamilton path, or is an orphan whose deletion breaks the shortcut. The latter necessarily has x internally with consecutive neighbors u∈I and v∈O and old reverse arc v→u.

Writing U for the old vertex set, the exact orphan number is

    O_x = sum_{empty ≠ S proper subset U}
          sum_{u in S∩I, v in (U\S)∩O, v→u}
          E(S,u) B(U\S,v).

**Proof.** Split an orphan at u,x,v. The two remaining path segments are disjoint and cover U. Conversely every displayed prefix/suffix choice joined through x is a Hamilton path whose deletion fails at the unique new shortcut. The prefix vertex set, endpoints, and two paths are uniquely recoverable. ∎

Thus

    H(T+x) = sum_P a_x(P) + O_x,
    H(T+x)−H(T) = sum_P(a_x(P)−1) + O_x.

Concrete smallest hostile: old rows (0,1) mean 1→0. Set 0→x→1. The old path has two legal positions, but the new cycle has three paths. The missing third path is an orphan. A three-vertex transitive base with rows (0,1,3), with only vertex0 beating x, has two insertion-slot paths and three orphans, hence H=5.

## 4. A deliberately nonuniform interface

Take A={0,1,2}, B={3,4,5}, and z=6. Fix z→A and B→z. Give each of the nine A–B edges an independent bit. Bit 3i+j is one when i→3+j and zero for the reverse. The interface state is k from 0 through 511. Internal A and B are either transitive T=(6,4,0) or cyclic C=(2,4,1), in outgoing-neighbor bitmask notation. We enumerate TT, CC, and CT.

The all-one state is a uniform cyclic wrapper. Moving off it changes actual intermodule adjacency and leaves S174's uniform-wrapper hypothesis. Most states are not uniform modules. The edge-flip formula applies without modification.

For each cube build the response graph whose edges are single interface-bit flips with |ΔH|≤2. Also build its induced subgraph on strong states. Counts of strong states are 400,511,448 in TT,CC,CT respectively. There are 324,275,321 full response components, and 279,274,292 strong response components.

**Component interval lemma (inherited criterion, explicitly realized here).** The H-values in any connected response component are a solid odd interval. Indeed all H-values are odd; take a graph path between vertices of minimum and maximum H. Each step changes H by at most2, so every intermediate odd value must be visited. This does not claim the original hypercube is connected using these restricted moves.

Largest-span certified strong components:

- TT: eleven states, values 153,155,157,159; states 233,235,241,249,346,458,459,474,489,490,491.
- CC: eighteen states, values 141,143.
- CT: eight states, values 151,153,155; states 91,115,123,227,243,251,339,371.

The receipt includes every component, every state, and every H count. Fragmentation is real: these cubes do not furnish a long unit-step interval or an induction extending intervals to infinity.

## 5. Reproducible strong prime and 7p witnesses

All tuples below are outgoing-neighbor masks on vertices 0,…,6. Each is strong and lies in the TT cube:

    H=49, k=279: (62,20,32,118,100,66,7)
    H=77, k=103: (62,36,8,114,102,68,7)
    H=91, k=74:  (22,12,8,113,102,71,7)
    H=133,k=202: (22,12,24,113,98,71,7)
    H=157,k=249: (14,60,24,112,97,69,7)

The complete receipt contains 30 distinct strong lane values across the cubes. In particular 77,91,133 have 7-adic valuation one, explicitly showing that the mirrored uniform-wrapper norm obstruction does not extend to these nonuniform interfaces. This does not invalidate the norm theorem, whose hypotheses have been changed.

The finite packet does not produce witnesses beyond known coverage, prove prime infinitude in the family, or close either infinite lane. A future scalable constructor needs a recurrence retaining actual endpoint states and proving overlapping realized components across increasing orders.

## 6. Checks

Run from any directory:

    python3 check.py > certificate.json
    python3 -O check.py > certificate.optimized.json
    cmp certificate.json certificate.optimized.json
    python3 audit.py > audit.json

The main checker enumerates all labelled tournaments of orders1,…,5, independently checks every H by literal permutations, checks every edge flip, and checks every insertion into orders1,…,4 against literal orphan classification. It then enumerates all 1536 interface states and checks every one of their nine directed neighbor transitions. There are 24474 edge checks total. Explicit exceptions, not assert, enforce all checks.

The separate audit imports no counting implementation: it rebuilds all 1536 seven-vertex tournaments as sets of arcs, enumerates all 5040 permutations for each, and checks component connectivity and interval coverage. These are independent implementations, not an independent human proof review.
