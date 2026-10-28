# Session 81: counterexample families, exact observer quotients, and sparse carriers

**Status:** elementary proofs and deductions from explicitly cited structural theorems; independent proof audits; reproducible finite-exact checks. No literature-novelty claim. This session builds on [Session80's symmetry composition](../session80_tournament_symmetry_20261027/README.md), rather than treating its existing results as new.

## The coherent picture

The strongest unifying statement is not that 7,21,11,13 or A1,A2 denote interchangeable objects. It is that a representation must be judged against its **future operations and requested output**. Several apparently strong observers provably forget the very coordinate that a later symmetry calculation needs. Conversely, some quotients are exactly closed under an entire operation grammar, so collecting more of the same observations can never restore that coordinate.

This packet establishes nine complementary parts.

### 1. A complete path observer can still hide exponentially differing symmetries

[LIFTING.md](LIFTING.md) proves that two labelled tournaments have identical cyclic triples if and only if they have identical Hamiltonian-path counts on **every named induced subset**, equivalently all induced path-cover counts. A legal module reversal supplies an explicit path bijection; the converse uses cited canonical decomposition theorems.

A constructive order-4·3^h family lies in one such observer fiber and realizes **every** automorphism order 3^(3^h+j), for 0≤j≤(3^h−1)/2. The old 12-vertex orders 81 and 27 become the first case of a gapless all-height symmetry distribution, not isolated coincidences.

### 2. The entire missing fiber is an explicit symmetry action

[FIBER_ACTION.md](FIBER_ACTION.md) gives the complete realization coordinates: two orientations per prime node and a total child order per linear node. Hypergraph automorphisms act by signed transports; ordinary tournament automorphisms are exactly coordinate stabilizers. An exact subgroup-fixed-point product and subgroup Möbius inversion recover the distribution, given the full acting group and its node action.

For the 12-vertex example, **128 labelled realizations** form **four tournament isotypes** under a hypergraph group of order **1296**. Their ordinary groups have orders 81,81,27,27. The all-height polynomial actually covers every isotype of this hierarchical fiber, although its original representatives do not cover every labelled realization.

This does not extract individual groups from aggregate prime cycle indices alone. Signed node transports and subgroup data are genuine additional inputs.

### 3. Procedural counterexamples persist under the allowed grammar

[OBSERVER_GRAMMAR.md](OBSERVER_GRAMMAR.md) proves C3 is the exact H-observable quotient for the full grammar of named restrictions and tournament substitution. Thus no amount of H readouts after these operations separates the collision family. Endpoint-sensitive guards are a genuinely richer observable and can separate it.

A finite refinement algorithm automatically extracts distinguishing operation words when a proposed summary is not closed. At order five, scalar H gives 7 labelled-observation cells, one deletion-refinement round gives 107, and two give the exact 208 cyclic-triple structures.

### 4. Sparse graph pairs fail at fixed size, but expanded carriers succeed

[GRAPH_CARRIERS.md](GRAPH_CARRIERS.md) distinguishes independent graph types from jointly aligned layers, and both from marked-tournament trijections. Fixed pairs/triples of graphs with n−1 edges on the original n vertices cannot encode all tournament types. Exact independent-pair capacity already fails at n=9, and independent triples at n=14, besides the elementary n=3 obstruction. These small thresholds are not claims about jointly aligned layers.

Yet one explicitly decodable **sign-even** graph on N=binom(n,2)+8n+1 vertices with N−1 edges suffices. Dense information is stored in a core, balanced by isolates. “Sign-even” here concerns automorphism edge-reversal signs, not even vertex degrees.

### 5. Preserve the group too, and the padding changes

[SPARSE_SYMMETRY.md](SPARSE_SYMMETRY.md) replaces interchangeable leaves by rigid port arms, then appends pairwise distinct rigid trees. The resulting single graph has E=V−1, is sign-even, faithfully reconstructs the tournament, and has **exactly its automorphism group**. It uses O(n² log n) vertices.

The extra Θ(k log k) vertices for k balancing tree components are optimal **for symmetry-preserving forest padding of the fixed core**, not for all conceivable encodings. Identical isolates are sufficient for class reconstruction but add unwanted permutations.

### 6. Natural conflict graphs have a much sharper realizability law

[CONFLICT_IMAGES.md](CONFLICT_IMAGES.md) examines the natural odd-cycle intersection graph Ω, not an artificial codec. Each component is K1, K2, or contains K4. Its forests are exactly matchings plus isolates. The first nonforest Ω with E=V−1 is **K4 plus three isolates**, realized by a **14-vertex tournament with H=243**; both the graph-size and tournament-order bounds are sharp.

The infinite family tK4 ⊔ (2t+1)K1 is realized on 11t+3 tournament vertices with H=3^(4t+1). A general SCC deficit equation generates further examples from any chosen large strong blocks. This is a useful distinction between “a faithful graph encoding exists” and “this graph occurs as a specific natural invariant.”

### 7. Root debt has a real tournament connection, but not the proposed 7/21 equivalence

[ROOT_COST.md](ROOT_COST.md) proves that directed-edge charges freely generate the root e_i−e_j **iff j is reachable from i**. Hence every nontrivial strong tournament has exactly the debt function of A_(n−1): C3 has positive debt somewhere, while all strong orders at least four are free. The same underlying lattice does not suffice without strong connectivity; transitive adjacent double roots still cost 2.

The audited inherited lattice classification is precisely for equal-unit ADE root sums with roots as charges: A1 and A2 are the only nonfree cases. Reducible sums can remove debt by perpendicular catalysts. The triangle's skew operator gives sqrt(−3), while Paley7 gives sqrt(−7), with the different lattices and operators stated explicitly. Dynkin A4 already evaluates to 21 under the independence polynomial despite being free, obstructing the simplest diagram analogy.

The H and H/|Aut| exclusions of 7 and 21 remain separately attributed inherited theorems, with different proof scopes. No numerical correspondence is promoted to an isomorphism.

### 8. Aliquot fibers reveal the same need to name missing coordinates

[ALIQUOT_OBSERVERS.md](ALIQUOT_OBSERVERS.md) starts from the inherited exact fibers s⁻¹(11)={21}, s⁻¹(12)={121}, s⁻¹(13)={27,35}. It proves that **every prime-power multiplier** q makes the pair (s(n),s(qn)) injective, with a guarded inverse.

A mixed multiplier genuinely fails: q=40 sends both 22,655 and 24,275 to observations (5857,1303480), and the complete inverse algorithm certifies exactly those two sources. The valuation profile on multiplier primes restores injectivity. This is an exact shared *descent criterion* with tournament attachment, not an identification of their objects or a claim about Collatz termination or prime untouchables.

### 9. Even group order, exponent and nilpotency class can all collide

[GROUP_SHAPES.md](GROUP_SHAPES.md) strengthens the counterexample factory again. At order 4·3^h, it constructs at least h pairwise nonisomorphic automorphism groups within one full path-observer fiber, all of order 3^(3^h+1), exponent 9 and nilpotency class 3. Moving the single active symmetry node through the composition tree changes the center and lower central series without changing these coarse invariants.

At 36 vertices the two groups both have order **59,049**, but their centers have orders **2,187 and 27**. Their commutator subgroups have orders 9 and 729. Both generated permutation groups were enumerated completely; identifying them as the full automorphism groups uses the canonical decomposition theorem, not an unperformed search through 36! permutations. Thus even knowing the prime modules and total symmetry order does not replace the location of the stabilizer action in the composition tree.

## What to regard as the proposed unifying framework

For each representation, specify:

1. **Objects and alignment:** labelled structures, isomorphism classes, or marked presentations?
2. **Observer and fiber:** which distinctions does the representation identify?
3. **Allowed operations:** does each operation descend to that quotient?
4. **Restoration coordinate:** a frame, a realization choice, subgroup action, valuation profile, or constrained charge set?
5. **Fidelity target:** scalar count, class reconstruction, full permutation group, or natural-image realizability?
6. **Hostile generator:** a seed collision plus a proved operation that preserves it, or a transport that deliberately splits its fiber.

The exact theorems above instantiate this framework. A useful next theorem must supply a real map and preserved predicate, not merely align small exceptional numbers.

## Reproduce the packet

All ten Python files use the standard library and explicit failure checks. From this directory:

```sh
python3 check_lifting.py
python3 check_fiber.py
python3 check_observer_grammar.py
python3 check_sparse_capacity.py
python3 check_codec.py
python3 check_sparse_symmetry.py
python3 check_conflict.py
python3 check_root_cost.py
python3 check_aliquot.py
python3 check_group_shapes.py
```

Each was also checked under `python3 -O`, with matching output. `check_observer_grammar.py` imports the adjacent `check_lifting.py`; the others are standalone. Results are printed, not required as preexisting receipt files. The reports distinguish public-checker ranges, additional independent audits, symbolic proofs, and inherited citations. No large-class enumeration beyond the stated bounds is implied.

## Remaining open boundaries

- An efficient, intrinsically decodable sparse carrier near the capacity lower bound Ω(n²/log n), or a stronger lower bound for a specified encoding model.
- Full connected realizability of Ω beyond the K1/K2/K4 gap, and the complete Hamiltonian-path spectrum.
- Which non-prime-power multipliers have injective two-observation aliquot maps; q=40 is a counterexample, not a classification or proved smallest failure.
- Efficient recovery of individual prime actions/subgroup transport from compact aggregate family data. The new fiber formula does not remove this input requirement.
- A root-debt model with suitably retained constraints that detects tournament orbit phenomena. The unrestricted strong-charge model provably loses too much information.
