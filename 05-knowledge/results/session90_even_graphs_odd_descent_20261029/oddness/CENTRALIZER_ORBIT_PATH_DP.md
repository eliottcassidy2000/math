# Centralizer-orbit compression of the folded path recurrence

Session90, 2026-10-29. **PROVED equitable-partition recurrence, FINITE-EXACT implementation.** The underlying orbit-aggregation principle is standard. This is a concrete reuse of the supplied centralizer from the paired-tournament lift, not a new A000568 formula or a claimed runtime record.

## Input contract

Input is an odd tournament W, a supplied anti-involution tau, and the **actual permutation subgroup** C=C_Aut(W)(tau). The recurrence also works for any supplied subgroup H of C, giving possibly weaker compression. Its output is still h_tau; obtaining the blue count requires division by the full |C|, not |H|. Dividing by |H| instead counts H-orbits of the fixed paths. An abstract group order is insufficient. Computing C is a separate cost; Session89's aligned pair generator already returns it as the pair automorphism group.

Tau fixes one root r and exchanges n mirror pairs. The inherited labelled recurrence has states `(S,v)`, where S is the set of used pairs and v is a member of one of those pairs. A transition appends w from an unused pair when v→w. Its terminal predicate is v→r. This counts tau-fixed Hamiltonian paths, not all Hamiltonian paths.

The companion `FIXED_ORBIT_DESCENT.md` explains this observer and credits the earlier recurrence at `05-knowledge/results/hollow_frame_blue_extrema_20261020/replay_subset_dp.cpp`.

## The weighted quotient

C permutes pairs, signed vertices, states, and transitions, preserving the level |S| and terminal predicate. For state-orbits A and B, let m_AB be the number of legal outgoing transitions from **one representative state of A** into B. This does not depend on the chosen representative: a group element transports its outgoing transitions bijectively.

Store `P_A`, the **total** partial-path count over every state in A. Initialize a level-one orbit by

`P_A = |A|`.

Process increasing levels and update

`P_B += P_A * m_AB`.

The output is the sum of P_B over terminal orbits. Divide this complete fixed-path total by |C| to obtain the number of blue ordinary presentations. Do not divide intermediate state counts by |C|.

**Proof.** For every source state s in A, let d_s be its partial-path count. Its contributions to B total d_s m_AB. Summing over s gives `(sum_s d_s)m_AB=P_A m_AB`. Levels make the state graph acyclic, so induction proves every orbit total. At level one each signed vertex contributes one path, giving |A|. The terminal predicate is constant on each orbit because C fixes r and preserves arcs. C acts freely on complete Hamiltonian paths, but need not act freely on intermediate states; this is why the division is postponed.

The computation needs one integer per state-orbit and one weighted update per nonzero m_AB. Constructing the orbit dictionary and transition weights incurs additional cost. Our audit implementation deliberately constructs all labelled states and verifies equitability at every state, so it is a correctness reference, not an optimized performance implementation.

## Hostile control: nonfree intermediate action

Use the nine-vertex tournament `C3[C3,C3,C3]`, with vertices `(i,j) in F3²`, outer arcs i→i+1 and inner arcs j→j+1. Let tau(i,j)=(-i,-j).

Its ordinary automorphisms are exactly

`(i,j) -> (i+a, j+b_i)`

with a and all b_i in F3. This has order81. The centralizer of tau has a=0, b_0=0, b_2=-b_1, hence order3. It fixes all three vertices in the i=0 block, including two nonroot vertices. Exhaustive9! permutations independently confirm the full group.

The labelled DP has64 states and96 legal arcs. The quotient has24 state-orbits and32 weighted arc entries. The state-orbit sizes are **four singleton orbits and twenty size-three orbits**. The fixed path count is33 and the blue presentation count is11. Literal Hamilton-word enumeration confirms H=3159 and exactly33 tau-fixed words.

A tempting wrong initialization gives every state-orbit mass1 and runs the same outgoing-weight recurrence. It returns19, neither33 nor11. The eight initial signed vertices form four orbits, not8/3 orbits. This explicit nonfree example prevents the accidental validity seen when every state-orbit happens to have size |C|.

## Complete small lift audit

The checker covers all531 normalized pair-matrix lifts at n=0,1,2,3. For each input it enumerates the full centralizer using all n! pair permutations. This is exact because commuting automorphisms fix the root and must preserve each normalized internal direction p_i→q_i, so they are precisely pair permutations passing the adjacency test.

Every quotient result equals the labelled pair-mask result, and every intermediate orbit total equals the sum of labelled state counts. The companion checker independently compares the latter with literal Hamilton paths on all531 lifts.

Aggregate state counts are12,420 labelled versus12,292 quotient; legal labelled arcs12,352 versus12,224 quotient weighted entries. Only eight of these531 labelled matrices have nontrivial centralizer. Each has C3, reducing24 states and24 arcs to8 state-orbits and8 weighted arcs. This is a real structural reduction but a small aggregate saving in this universe; no wall-clock improvement follows automatically.

For n=0, the single root has one empty first-half path and needs no DP states. All arithmetic is integral; no averaging or rounding is used.

## Reproduction

Run the companion files from the same directory:

```
python3 check_centralizer_path_dp.py > orbit_dp.json
python3 -O check_centralizer_path_dp.py > orbit_dp.optimized.json
cmp orbit_dp.json orbit_dp.optimized.json
```

Standard library only. Explicit active checks remain under optimization. The checker imports only the local companion's codec and basic adjacency routines; the compressed recurrence and its exhaustive centralizer tests are implemented here. It uses an independently accumulated labelled-state recurrence for all intermediate comparisons, and a literal9! path/group control for the nonfree hostile.

This algorithm does not compress arbitrary future attachment operations. It compresses the stated folded-path count and terminal predicate under a known group action. Shared ports or later guards that distinguish these state-orbits must be added to the state and stabilizer before reusing the quotient.
