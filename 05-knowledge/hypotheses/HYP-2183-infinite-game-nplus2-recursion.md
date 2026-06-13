# HYP-2183: infinite-game n+2 recursion turns V* into a period-3 automaton and the LRC row/column program into omega^2 proof rank

**Status:** SUPPORTED synthesis and finite scout. This is not a proof of LRC n=14; it is a proof-program improvement that merges S575, HYP-1881/HYP-1891, HYP-2177, HYP-2181, and incoming S599t's partition-function/infinite-game HYP-2182.

**External prompt note:** The linked X status did not expose its text to the web fetcher. Public Hamkins pages verify the relevant background of infinite games, infinite Go, and finite positions with ordinal values such as `omega^2+k`; the mathematical import here uses only the user-supplied cue "n+2 recursion / infinite go" plus that public infinite-games context.

## Claim

The clean even-LRC `n -> n+2` recursion should be treated as a finite automaton carried by `C=2n-1`, while the older denominator recursion should be treated as a well-founded infinite-game proof rank. In particular:

1. On the clean even ladder, `n -> n+2` sends `C=2n-1` to `C+4`, so `C mod 3` advances by `+1`.
2. By HYP-2177, the doubled AP sporadic `V*=AP[(n-2)->2(n-2)]` is tight on the `(2n-1)` pinch lattice exactly when `3|C`.
3. Therefore `V*` is alive on exactly one state of a period-3 `n+2` automaton:

   `n = 8,14,20,26,...`, equivalently even `n == 2 mod 6`.

4. The HYP-1881/HYP-1891 row/column recursion has ordinal shape `omega^2`: odd-root column motion `b -> b+2` is the outer `omega` ladder and dyadic row motion `N -> 2N` is the inner row coordinate.

## Proof Program

Define the denominator rank

`rho(N) = omega*((odd(N)-1)/2) + v_2(N)`.

The intended induction has two moves:

- **row discharge:** for fixed odd root, prove that all dyadic lifts descend in finite row depth;
- **column transfer:** after row discharge, move the odd root by `b -> b+2` using a finite automaton over the residue phases of `C=2n-1`.

This imports the correct infinite-game lesson: do not try to win the proof game by one scalar invariant on ever-larger boards. Use a well-founded rank for the proof state and a finite automaton for the families that persist through infinitely many successor steps.

Incoming S599t supplies the global partition-function reading: infinite-game ordinal recursion, tournament `H` spectra, and LRC covering depth are all spectra of decomposable partition functions. Incoming S613's Theorem D scout also verifies that loose `V*` gaps center on the `(2n-1)` lattice. HYP-2183 is the local LRC specialization of both signals: the partition-function/free-energy viewpoint tells us to keep rank and decomposable side channels, while the `C mod 3` automaton tells us exactly which clean `n+2` columns keep the `V*` zero alive.

## Consequences For LRC n=14

The n=14 obstruction is no longer an isolated sporadic. It is the `C=27=3^3` member of the alive phase of the `n+2` automaton. The next transfer target should be:

- prove the period-3 automaton theorem for the clean `n+2` ladder;
- plug HYP-2177's lattice proof into the alive phase;
- use HYP-2167/HYP-2166 carry-fiber conservativity to lift from the `Res_27` quotient to integer rows;
- keep HYP-2181's carrier-side channels instead of collapsing them to raw scalar H or raw exponent numerology.

This reframes the n=14 task as "rank-decreasing row discharge plus finite column automaton" rather than as a one-off fight against `V*`.

## Assumption Challenge

The scout does not assume tournament vertices are runners or arcs. The useful vertices for this session are proof routes and recursion states:

- clean dihedral polygon ladder;
- period-3 `V*` automaton;
- `<2,-1>` shell-orbit ledger;
- `omega^2` two-mode proof rank;
- unit-distance/infinite-game carrier analogy;
- raw next-`n` brute force;
- unseen X-post text as theorem.

The quotient preserves the predicate "does this proof route retain the data needed to decide the `V*`/n+2 family?" It destroys individual lift/carry choices, non-single-swap sporadic detail, and row-specific owner certificates. That loss is intentional: those side channels must be reattached through HYP-2166/HYP-2167 rather than hidden inside the automaton.

## Evidence

S616 computes the clean even ladder through `n=30`:

- `V*` alive columns: `[8,14,20,26]`;
- automaton law: `C -> C+4`, so `C mod 3` cycles through the three states;
- route tournament: transitive, score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, no directed 3-cycles, one Hamiltonian path;
- top routes: period-3 `V*` automaton, then `omega^2` two-mode proof rank, then `<2,-1>` shell-orbit ledger.

## Related

- `04-computation/infinite_game_nplus2_recursion_s616.py`
- `05-knowledge/results/infinite_game_nplus2_recursion_s616.out`
- `07-reflections/infinite-game-nplus2-recursion-s616.md`
- HYP-1881, HYP-1891, HYP-2091, HYP-2177, HYP-2178, HYP-2181, HYP-2182/S599t, incoming S613 Theorem D
- S575 clean simplex/polygon recursion; HYP-2166/HYP-2167 `Res_27` quotient and carry-fiber conservativity.
