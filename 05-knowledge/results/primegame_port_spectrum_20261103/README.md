# Growth faithfulness, the Boolean inverse boundary, and PRIMEGAME guards

**PROVED all-order structural theorems / FINITE-EXACT constructions /
CITED inherited frontier. Same-agent audit with independent algorithms.**

The assertion that7 and21 are the only missing positive odd tournament
Hamiltonian-path counts remains **OPEN**. This session does not claim a
completed proof or an extension of the existing99,999,999 coverage.
It proves a new scoped carrier theorem and locates a stronger averaging
obstruction while testing a concrete construction grammar.

## What changed

1. **All-order growth faithfulness.** Four actual source-preserving growth
   queries recover a signed old-edge derivative. Iterating recovers the full
   Boolean derivative response. A proved nonsingular fair-point Hessian then
   shows that named growth responses distinguish EVERY labelled tournament.
   At fixed order, a carrier closed under this language needs the entire
   orientation:2^binom(n,2) states, or2^binom(n-1,2) in a fixed-path chart.
   [Proof and precise operations](GROWTH_FAITHFULNESS.md).
2. **Averaging remains blind at every attachment size.** Two different
   positive parity laws on old orientations have equal means after EVERY
   fixed source-oblivious attachment. Their shared-source moments can only
   separate them from order ceil(n/2) onward. At n3 a two-birth example
   attains that threshold: mean10 for both, second moments113 and111.
   Allowing shared copies and ordinal sums recovers all rational source
   weights by finite interpolation. An H-only inverse cannot be unimodular:
   all its query rows are identical modulo2.
3. **The81/27 hostile is retained literally.** Both inputs have67 identical
   cyclic triples and H462915. Appending a sink and then a vertex beating
   precisely old port8 gives H2390013 versus1754919.
4. **A small constructor language has a substantial exact prefix.** All
   4,037,914 recursive-tree presentations through11 vertices were checked.
   The saved11888 witnesses cover every allowed odd value through20405;
   all witnesses have a second exact counting check. A runnable constructor
   returns the actual adjacency. [Finite scope and branching obstruction](TREE_CONSTRUCTIONS.md).
5. **The numerical prompts acquire specified operators and guards.** The
   five mixed blue nodes induce P4 plus an isolate, giving golden-ratio
   eigenvalues. The eight self-converse blue nodes form a leafed C4, giving
   silver-ratio eigenvalues. PRIMEGAME's quartets have respectively1,2,0
   executable complete orderings over all positive inputs. The last
   quartet's product is integral at34, but no ordering executes ANYWHERE.
   [Spectra, recurrence and full source cells](BLUE_SPECTRA_AND_PRIMEGAME.md).

At each tested order from 3 through 7, all labelled tournaments are also
distinguished by H plus its first edge derivatives alone: 2,097,152 exact
states at n=7. This first-derivative assertion fails at n=2. Extending
this shallower bound to all n is OPEN; the all-order proof uses higher
derivatives. This is distinct from legality-only relocation observations.

## Reproduce

From the repository root:

    python -B 05-knowledge/results/primegame_port_spectrum_20261103/port_scout.py
    python -B 05-knowledge/results/primegame_port_spectrum_20261103/growth_identity.py
    python -B 05-knowledge/results/primegame_port_spectrum_20261103/spectra_primegame.py
    python -B 05-knowledge/results/primegame_port_spectrum_20261103/verify.py
    python -B -O 05-knowledge/results/primegame_port_spectrum_20261103/verify.py --output verification_optimized.json
    python -B 05-knowledge/results/primegame_port_spectrum_20261103/construct.py 20405
    python -B 05-knowledge/results/primegame_port_spectrum_20261103/replay_cpp.py
    python -B 05-knowledge/results/primegame_port_spectrum_20261103/seal.py --check

Only spectra_primegame.py needs SymPy; the other Python verifiers use the
standard library. replay_cpp.py needs a C++17 compiler and regenerates the
two full tree censuses, the complete first-derivative census and the
independent all-witness receipt. It writes temporary executables outside
the repository. [Audit](AUDIT.md) and [next decisive tests](HANDOFF.md).
