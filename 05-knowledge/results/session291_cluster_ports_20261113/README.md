# Session291: exchange without forgetting the interface

**Status:** independently reviewed research synthesis. Chan–Pak remains
open. Exact statements, conditional programs and analogies are separated below.

The starting suggestion was to compare twelve bases with two ports each to the
Leech lattice, and to use cluster-algebra ideas to change tournament cores. The
useful answer is not a shared number 24. It is an exact distinction between
**changing a presentation**, **changing an observable**, and **changing a physical
tournament while controlling every continuation**.

## Results and limits

1. [Paired-port arithmetic](ports/PROOF.md). The inherited 24 charts belong to
   twelve marked `(adjacency, selected-set)` sources, with only eight distinct
   literal adjacency arrays. Their eight-base exponential coefficient matrix
   has rank eight. The twelve paired differences have rank seven and index
   `2^9·3^6` in the saturated evaluation kernel. For every integer `m>=2`,

       gcd_i (U_(2i+1)(m)-U_(2i)(m)) = gcd(2^m-2,3^m-3).

   This yields an exact prime-power schedule for when port differences become
   invisible to a modular observation. All 4096 binary choices have distinct
   **formal aggregate** values for every `m>=3`; this is an all-order theorem,
   not extrapolation. The coefficient observation modulo 13 retains all choices;
   modulo 2 retains only four images. None supplies a tournament addition gadget.

2. [Cluster/exchange lane](cluster/PROOF.md). Ordinary Fomin–Zelevinsky quiver
   mutation stays inside tournaments only at sources or sinks. It cannot be
   silently used as the desired core-changing move. An alternative exact
   marked-path elimination atlas retains subset and endpoint data. Its exchange
   identities change presentations of the same inventory, not automatically
   physical tournaments. The original squarefree path matrix is inherited work.

3. [Witt-inspired carrier](witt/PROOF.md). The attached paper is Sulyma's
   *RO(G)-graded norms for de Rham–Witt forms*. Its useful methodological theme
   is to keep gradings and mixed carries before collapsing them. Classical
   primitive-necklace constructions supply positive integral corrections, but
   their scalar values for odd-prime path orbits are even, hence cannot directly
   be tournament Hamiltonian counts. A square-zero vertex-resource ring gives
   an exact determinant/odd-cycle carrier and an endpoint resolvent. A block
   Schur complement keeps reentries that scalar contraction would lose.

4. [Historical connection contracts](history/CONNECTIONS.md). The repo already
   contains a genuine Paley/Golay bridge, a score-lattice/Leech-cusp bridge, and
   endpoint/cluster work. Each uses a different map and forgets different data.
   Those contracts matter more than a resemblance between dimensions.

## What the Leech comparison actually says

The canonical paired difference code is the self-dual `[24,12,2]` doubled
repetition code. Its Construction A is a rotated integer lattice, not Leech.
A proper parity-restricted subcode cannot remain dimension twelve and self-dual.
A formal order-three operator on each pair gives `A_2^12` with roots, not Leech;
it also leaves the positive physical cone. These are exact obstructions to the
naive identification, not obstructions to every future bridge.

A meaningful new bridge must derive **cross-pair coupling from endpoint
continuation data**. Imposing Golay by an arbitrary relabeling would not explain
the tournament geometry. A single independent port swap cannot preserve any
self-dual distance-eight binary code; coupled swaps with support at least eight
are a possible test, not a proved port symmetry. No intrinsic Golay code or
Leech lattice has been extracted here.

## Next constructive question

See [the synthesis](SYNTHESIS.md): seek a quotient of the marked carrier that is
stable under chosen physical attachments, and a well-founded inverse rule that
lowers the requested target while keeping that quotient feasible. Exact modular
blindness can screen candidate moves; it cannot replace absolute anchoring or
physical realization. Neither the presentation atlas nor the pair lattice proves
an eventual threshold, eliminates all missing values except 7 and 21, or closes
Chan–Pak.

## Reproduction and receipts

Run from the repository root:

```sh
python3 05-knowledge/results/session291_cluster_ports_20261113/ports/verify.py
python3 05-knowledge/results/session291_cluster_ports_20261113/cluster/verify.py
python3 05-knowledge/results/session291_cluster_ports_20261113/witt/verify.py
```

Each lane gives its exact dependency and scope. Repeat with `python3 -O` to
check assertion-independent controls. See [independent review](audit/REVIEW.md) and its final-source hash manifest.
The review covers the new claims and stated dependency boundaries, not a global
re-audit of the repository or the attached paper.
