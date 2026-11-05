# A surviving symmetry blocks universal two-sided two-port digits

Session 197, 2026-11-05. **All-order obstruction for an exactly specified physical future language.** This does not exclude target-dependent constructions from small asymmetric seeds and does not claim any new global Hamilton-count gap.

## 1. Exact exposed-vertex rules

The initial input is a fixed tournament T together with exactly two distinct ordered exposed old vertices. Every step:

1. adds one new vertex;
2. either makes its old out-neighbor set a subset of the two current ports (rightward birth), or makes its old in-neighbor set such a subset (leftward birth);
3. forgets one of the two old ports or the new vertex, retaining the other two, in their inherited order.

A forgotten vertex remains in the graph. No old arc changes. There are always exactly two distinct exposed ports. **An old vertex that was not initially exposed can never become exposed later**, and a forgotten vertex can never be re-exposed: only a surviving current port or the newly born vertex may be a next port.

Both birth directions, every subset including empty and both ports, and every forgetting choice are allowed at every step. Choices may depend arbitrarily on the target and the entire graph. The obstruction below does not assume a periodic word, fixed template, bounded order, or finite-dimensional numerical state.

If a different architecture allows arbitrary old vertices to be re-exposed, nonuniform incidences to non-ports, or old-arc changes, the theorem is not asserted for that architecture.

## 2. Pointwise stabilizers survive every permitted birth

**Theorem 1.** Let G be a subgroup of Aut(T) fixing both initial ports individually. Extend every g in G by fixing all subsequently born vertices individually. For every permitted word, this extension is an automorphism of the output tournament and fixes both current ports. Consequently

    |G| divides H(output).                            (1)

**Proof.** Initially the assertion holds by hypothesis. Assume it holds just before a birth. Every current port is either an initial port still retained or a later-born vertex, so every g fixes it. The new vertex distinguishes only those two fixed ports individually; all remaining old vertices receive the same orientation relative to it. Thus g preserves every new edge when extended to fix the new vertex. It preserves old edges by induction. After any forgetting choice the two surviving ports are still fixed. This proves the assertion for every finite word.

An automorphism acts on Hamilton paths by applying it to each position of the vertex word. If a path word is fixed, every one of its named vertices is fixed, so the automorphism is the identity. The action of G is therefore free and partitions the finite path set into orbits of size |G|. This gives (1). QED.

The free-action principle is also recorded in the [S196 symmetry packet](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/symmetry/PROOF.md). The present assertion applies it to the exact persistent-port language above. No equivariance assumption on the strategy is needed: every individual allowed next action preserves the fixed subgroup, however that action was chosen.

## 3. An explicit strong input defeating every initial pair

Take a directed triangle of directed-triangle modules:

    T=C3[C3,C3,C3].

The modules are `M0={0,1,2}`, `M1={3,4,5}`, `M2={6,7,8}`. Inside each module orient the cyclic order. Between modules orient `M0 -> M1 -> M2 -> M0`. Its rows are

    [58,60,57,464,480,456,135,263,71],
    H(T)=3159.                                       (2)

The graph is strong because both the quotient cycle and each module are strong. The path count is independently verified by all 9! vertex permutations and by vertex-subset endpoint DP.

For each module Mj, its cyclic rotation by one step, fixing every other vertex, is an automorphism of order 3. These three disjoint rotations commute.

Any two initial ports meet at most two modules. At least one module is untouched, and its rotation fixes both ports individually. Theorem 1 gives

    3 divides every descendant Hamilton count.        (3)

If both ports lie in the same module, the two untouched rotations generate a subgroup of order 9, so 9 divides every descendant count. If they lie in different modules, the guaranteed subgroup has order 3. No assertion that these are the whole pointwise stabilizers is needed.

In particular, starting from this strong physical input, no possible initial port pair and no later adaptive two-sided word reaches either

    3*3159+2=9479 or 3*3159+4=9481.                   (4)

More strongly, no descendant of this input can have either literal nonzero ternary child of its own current H, since that current H is divisible by 3 and the two children are not.

Thus permitting both birth directions does not yield a per-input universal digit theorem. This obstruction is independent of the infinite linear-rank result for two-sided continuation: it is a persistent physical automorphism, not a failure of a particular coordinate system.

## 4. What this does and does not rule out

The theorem excludes an induction hypothesis claiming that **every** strong input, with freely chosen two initial ports, can always reach both nonzero ternary children by this language. It does not exclude a carefully selected family of asymmetric inputs, a target-dependent construction beginning from a small seed, or a procedure that changes the permitted exposed-vertex rules.

The concrete H49 input used by the Session 197 leftward H149/H151 gates is a useful contrast. Enumerating all 7! permutations of its explicit seven vertices finds only the identity automorphism. Thus this particular surviving-symmetry obstruction does not apply there. The positive paired-fan construction and this strong hostile address different quantifiers and are consistent.

No claim is made that 9479 or 9481 are absent from the unrestricted tournament spectrum, or even from the same two-sided language starting at a different input. Nor is a literal return to a marked state inferred from any scalar attainment.

## 5. Reproduction

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The standard-library verifier uses explicit exceptions. It checks the literal nine-vertex graph, strongness, both independent H3159 counts, and the three actual order 3 rotations. For all 72 ordered initial port pairs it checks every one of the 24 first births and every one of the 24 second births, including both directions and all forgettings: 1,728 one-step and 41,472 two-step adjacency-level preservation checks. It independently checks subgroup divisibility on 72 selected one-step graphs and the rigid H49 contrast.

These finite controls are not substituted for the inductive proof of automorphism survival. The divisibility and all-order nonattainment assertions follow from Theorem 1 and the explicit module rotations. No dirty checkout or worktree was read or changed.
