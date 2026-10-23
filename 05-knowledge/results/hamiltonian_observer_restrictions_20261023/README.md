# Hamiltonian observers, conditional kernels, and guarded return clocks

23 October 2026. Status: **PROVED** statements on the explicitly declared carriers below, with **FINITE-EXACT** controls and independent audits. The all-size native tournament-metagraph wiring rule remains **OPEN**. The native incidence coupling, the automorphism-frame graph, and the auxiliary return chain are three different graphs/operators.

## 1. A Hamiltonian path is an observer, with a retained choice

Every finite tournament has a directed Hamiltonian path. The elementary insertion proof constructs one: given a path through all but `v`, insert `v` before the first path vertex it beats; the preceding vertex, if any, beats `v`. If it beats no vertex, append it. The new consecutive arcs are directed correctly. This is the standard existence mechanism, not a new theorem of this session.

Once a path `P=(v_0,...,v_(n-1))` is supplied, it supplies an actual linear observer order. Every arc is then a forward/backward bit, and the consecutive bits are forced forward. The remaining `binom(n-1,2)` bits are lossless *with the actual list P retained*. Forgetting `P` forgets the owner-to-position alignment.

Existence does not supply a label-free unique observer. On a tournament with a nontrivial automorphism, an automorphism-equivariant choice of one ordered spanning path is impossible: the chosen path would have to be fixed by that automorphism, but fixing every position in a spanning list fixes every vertex. An external order or a distinguished reference path can make a choice, at the cost of retaining that additional input. The [conditional compiler](compiler/README.md) makes this distinction explicit for an unpointed candidate fibre versus a chart-dependent completion.

Our inheritance pass used the incoming full-cube `8/9` mechanism, the regular order-63 nine-flag hostile, the corrected native-versus-auxiliary quotient distinction, and the underused harmonic-extension sidecar. The live concept board is: **observer choice; shared input-prefix coordinates; guarded moves; conditional expectation; integer clocks; cycle transport**. These concepts are compared below through concrete maps, rather than through repeated small numbers alone.

## 2. Structural progress on the native incidence coupling

[The all-size connectivity proof](coupling/NATIVE_CONNECTEDNESS.md) shows that the native fixed-path incidence coupling has **exactly two components** at every order `n>=3`: the transitive/empty singleton and all nonempty compressed graph types.

The constructive mechanism is a backward-arc descent. A compressed vertex with degree between 2 and `n-3` provides two neighbors and a nonneighbor. They create a strong four-vertex prefix with a certified outgoing boundary arc; a different Hamiltonian presentation reduces its backward count. Matchings, spanning stars and complete graphs have explicit separate repairs. This retains the boundary condition that a naive local reorder would lose.

This result concerns the incidence relation between compressed ordinary graph classes and native tournament classes under a supplied Hamiltonian path. It does not assert connectivity of the blue/black line metagraph or of a fixed-source guarded path fibre.

The [restriction analysis](coupling/RESTRICTIONS.md) also determines why the old spectral bound fails. The isolated transitive state prevents any positive global gap. Even after removing it, order four has exact gap `7/15`, not `8/9`. Source complement pairing and target converse pairing produce a connected order-four chain with gap `8/15`, still insufficient.

For any specified weighted incidence table, the surviving gap is governed by an exact PSD/maximal-correlation criterion. Coarsening observations on the same probability space is safe. Changing positive point weights by a ratio at most `b/a` retains at least `(a/b)` times the original gap. Hard restriction requires a new argument.

## 3. Native shared-frame completion: count, observe, move

For the authenticated mixed-depth ternary cells from incoming commit `8b71b746f7a9`, each actual automorphism is parametrized by coordinates

    a_(depth, shared-cell-component, input-prefix) in F3.

An indexed owner answer fixes the coordinates along its source address. Shared ports therefore impose shared equations, not independently multiplied counts. If `D` is the total number of coordinates and `F` the consistently forced set, then

    number of remaining actual path lifts = 3^(D-|F|).

A requested tuple exposing `q` new coordinates has exactly `3^q` possible answers, with uniform remaining lift multiplicities. The minimum number of additional actual arc flags needed for required untouched bottom slots `R` is

    |R| - maximum_matching_size(J[R]),

where `J` uses consecutive arcs at their actual path clocks. This is a targeted conditional rule, extending the incoming whole-path flag theorem.

Changing one free coordinate across all cells that share it realizes a genuine subtree rotation. With `r=D-|F|`, the resulting actual frame graph is `H(r,3)`: connected, diameter `r`, and gap `1/r` for the specified random-scan heat-bath search. The group law is nonabelian even though its parameter set is a coordinate cube.

These moves preserve the ordered tournament bits: `gP` and `P` have the same backward word for `g in Aut(T)`. They move actual owner lifts and shared-port alignments within a word fibre; they do not by themselves traverse different metagraph owners.

There is now an exact conditional way to reduce a guard. If an authenticated guard depends on only `q` exposed free coordinates, with allowed set `A subset F3^q`, then its full surviving frame graph is

    H(q,3)[A] Cartesian-product H(r-q,3).

Its components and visible bottlenecks are therefore determined entirely by the first factor; all invisible kernel lifts survive in the second. This requires a proof that the real guard factors through that observation. It does not assert that the native electrical mod-8 guard is local. A tested three-frame guard already disconnects single-coordinate moves.

Complete proofs and finite universes: [compiler/README.md](compiler/README.md).

## 4. A positive gap survives arbitrary guards at return times

The compressed free bits also encode an auxiliary tournament on `n-1` vertices. Its full-cube sandwich has gap `8/9` for `n>=4`. Adding a uniformly refreshed actual owner order retains that gap.

For any nonempty accepted subset `S` of this lifted space, let `C` be its complement. The watched/return kernel is

    R = Q_SS + Q_SC (I-Q_CC)^(-1) Q_CS.

The [self-contained proof](trace/TRACE_WITH_CLOCKS.md) gives gap at least `8/9` for **every** such `S`, including arbitrary exact source, port and guard predicates. Its proof uses harmonic extension and the Poincare inequality. This is an application of the standard trace construction, not a claim of a new general Markov-chain theorem.

The retained clock is essential. From the stationary conditional distribution, the mean excursion length is `1/pi(S)`. For a fixed actual source `T` without additional guards this is

    n! 2^binom(n-1,2) / H(T).

All accepted return states preserve the requested data. Intermediate excursion states need not. To obtain a native legal walk, the remaining task is to lift return transitions to allowed guarded paths and charge their lengths and congestion. The theorem does not bypass that task.

The earlier native electrical example explains the same information issue on a different carrier: boundary energy may be preserved while integer-current parity is lost. We retain return targets, duration generating functions, harmonic extensions and original guard receipts as distinct data.

## 5. What the 13/4, mod-8/mod-9 and geometry routes actually connect

The [recovered valuation chain](history/VALUATION_CYCLE_GAIN_LIFT.md) is

    13 -> 2 -> 5 -> 3 -> 7 -> 4.

Identifying the endpoint 4 with the input 13 gives five quotient states. The reset is an identification, not an additional valuation step. The historical normalization matters: replacing the Baby Monster denominator by its actual double-cover centralizer changes this to a two-cycle after the same identification.

A new gain decoration `g_i=x_(i+1)/x_i` preserves what the quotient loses. The five gains multiply to `4/13`; a reset costs the reciprocal `13/4`. Gauge changes can move this defect to one seam but cannot erase it. With an anchored starting value and indexed gains, the original path labels are recovered exactly. A complex fifth-root gauge produces a contracting spiral, whose logarithmic branch cut changes the displayed argument but leaves the real displacement `log(4/13)` intact.

For rational gains, valuation-plus-unit registers give

    rho_(2,3)(4/13)=(2,5),    rho_(3,2)(4/13)=(0,1).

The mod-9 register hides this particular defect; the mod-8 register with its integer valuation detects it. Their pair still misses 73, so finite congruence agreement cannot certify exact closure. These registers are not substitutes for the native electrical mod-8 invariant.

The recovered spherical/Euclidean/hyperbolic mechanism uses a different indispensable sidecar: dart reversal and cyclic rotation determine faces. For cubic maps with face size `q`, the curvature factor changes sign at `q=5,6,7`. The same Heawood graph supports genus one or genus three depending on the relative rotation choices. A Hamiltonian order, a scalar cycle gain, and a rotation system retain different information; an explicit transition map is needed to connect them.

## 6. Next precise obligations

1. **Guard footprint:** determine which actual owner coordinates the original mod-8 and port guards inspect. Apply the product decomposition only after proving that factorization; retain complete original receipts otherwise.
2. **Coupled guarded blocks:** when the visible factor is disconnected under single-coordinate moves, find the smallest simultaneously allowed block that joins components, or prove a conserved obstruction. The invisible kernel factor need not be re-enumerated.
3. **Native spectrum:** use the new all-size backward-arc descent to bound the nontransitive incidence component quantitatively. Connectivity alone does not restore `8/9`.
4. **Return-path lifting:** realize auxiliary return transitions by native guarded moves, retaining integer lengths, original port ownership and harmonic-extension data. A gap in return-event time is insufficient without this cost.

## 7. Reproduction, provenance and audits

The [late incoming integration](INCOMING.md) distinguishes Session43's full backward graph from our compressed graph: their isomorphism partitions are incomparable from order five onward. It also extends the incoming native phase assembly to an exact `C_(2p)`-to-`P_(2p-1)` representative-cut construction, and connects bounded-range modular guards to the incoming bounded-degree Jacobian result without claiming a general transfer.

Run `python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/reproduce.py` from the worktree root. This runs each exact checker and verifies that frozen receipts reproduce. [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json) pins the inherited inputs. [AUDIT.md](AUDIT.md) records proof reviews, hostile controls and scope boundaries.

No theorem ID was reserved, no old frontier was silently promoted, and no proof of Collatz, LRC, Hadwiger–Nelson or the planar Jacobian conjecture is asserted. The progress is a constructive all-size incidence result and a more precise native shared-frame compiler, with the guard and clock obligations explicitly exposed.
