# Infinite-Game n+2 Recursion: Period-3 V* and omega^2 Rank

S616 was prompted by the user's X link and "n+2 recursion / infinite go" cue. The exact X status text was inaccessible, so the session used only verified public Hamkins context about infinite games, infinite Go, and ordinal game values, plus the user's mathematical phrase. That restraint was useful: the repo already had the real n+2 material in S575 and the two-mode recursion in HYP-1881/HYP-1891.

Rebasing brought in S599t's partition-function synthesis, which had already claimed HYP-2182 for the global infinite-game frame. It also brought in S613's Theorem D scout, verifying that loose `V*` gaps sit on the `(2n-1)` lattice. S616 therefore becomes HYP-2183: the local LRC automaton that sits under the partition-function view and specializes the lattice signal to the clean `n+2` ladder.

The new observation is small but sharp:

`n -> n+2` on the clean even LRC ladder sends `C=2n-1` to `C+4`.

Modulo 3, that is a three-state successor game. Since HYP-2177 says the doubled AP sporadic `V*` is tight on the lattice exactly when `3|C`, the infinite-looking family is not miscellaneous:

`V*` appears at `n=8,14,20,26,...`, exactly the alive phase `n == 2 mod 6`.

So n=14 is the `C=27=3^3` member of a period-3 automaton, not a lone exception.

The other half of the insight is the infinite-game rank. HYP-1881/HYP-1891 split denominator recursion into odd-root columns and dyadic rows. Write

`rho(N)=omega*((odd(N)-1)/2)+v_2(N)`.

Then row work is finite descent inside a column, and column work is the outer successor ladder. That is the `omega^2` import: prove that the proof game cannot cycle by assigning it a rank, then separately solve the finite automaton of persistent states.

## Tournament Analysis

The vertices were deliberately not runners. They were proof routes:

- `period-3 V* n+2 automaton`
- `omega^2 two-mode proof rank`
- `<2,-1> shell-orbit ledger`
- `clean dihedral polygon ladder`
- `unit-distance/infinite-game carrier analogy`
- `raw next-n brute force`
- `unseen X-post text as theorem`

The pairwise observable was usefulness for the LRC n+2 proof route: preserves the right predicate, keeps rank data, offers a finite automaton, retains side channels, gives an actionable next lemma, and avoids external-text risk. The tournament was transitive with one Hamiltonian path. The ranking put the period-3 automaton first, the omega^2 proof rank second, and the shell-orbit ledger third.

## What This Improves

This pushes the S615/S2181 carrier story one notch closer to a proof tool. HYP-2181 identified `gcd(3,2n-1)` as the live carrier obstruction across LRC and tournament-style gaps. S616 turns that into a successor automaton: the obstruction phase is predictable under `n+2`.

For LRC n=14, the next lemma should not be "find more sporadics." It should be:

1. prove the period-3 automaton for the clean ladder;
2. splice in HYP-2177's lattice proof for the alive phase;
3. discharge lift/carry side channels through HYP-2166/HYP-2167;
4. treat the row/column recursion as an `omega^2` proof game, so every reduction visibly lowers rank or returns to a finite automaton state.

That is a real mathematical improvement because it shrinks an infinite family to a finite successor machine plus a well-founded rank.

## Artifacts

- `04-computation/infinite_game_nplus2_recursion_s616.py`
- `05-knowledge/results/infinite_game_nplus2_recursion_s616.out`
- `05-knowledge/hypotheses/HYP-2183-infinite-game-nplus2-recursion.md`
