# What the compensated h=77 lock buys, and its next exact obstruction

The [two compensated fold arms](FOLD_FIRST_SEAM.md) have base counts 173 and 175 at ordinary h=77 and common normalized start 181=oddpart(S+h). They synchronize with count-neutral ones; a common zero then adds 181·2^(j−1). This note tests that *exact* scale rather than switching silently to a nearby slope or comparing residues.

## No third neighboring arm at the same normalized start, at any order

By the [fixed-state fibre inverse](INVERSE_STATE_FIBRE.md), enumerating all non-singleton strong SCC factors dividing 77 is finite. [Busch](https://doi.org/10.37236/1141) bounds their order by nine. The earlier [order-nine canonical generator](../chanpak_tail_return_20261115/NINE_ORDER.md) supplies a complete isomorphism list, and every mask of every qualifying strong core is visited. There are respectively 0, 32, and 2770 distinct nonsingleton local signatures of counts 7, 11, and 77, with a literal endpoint DP control for each signature. Add the two singleton signatures.

The putative neighboring base counts 171 and 177 are coprime to 77. Run the all-order normalized-start scan from the h=15 proof with total target h=77 and cap X=177. The same proof applies: N cannot grow less than by the ordinary multiplier; a future 0 or a nontrivial partially unmarked SCC is impossible once S>X+B, where B bounds the absolute bias. An all-marked nontrivial SCC would introduce a nontrivial divisor of 77 and is incompatible with a unit target; all-one singleton padding is neutral in both N and oddpart(S+h). At S=0 all histories with the same h already have the same full state. The finite replay has B=77.

The two complete normalized lists for 171 and 177 have respectively two and seven entries, **neither containing 181**. In contrast the known 173 and 175 arms both meet 181. Therefore no third adjacent base value can join this exact h=77, λ=181 return via any spine or count-neutral transitive preparation. The program also prints a bounded diagnostic count for 175; because 175 is not a unit modulo 77, its whole normalized list is not certified by this stopping rule. Its inclusion at 181 is literal independently. This excludes widening at this signature, not at other h or λ, and not by an operation that changes the absolute base readout.

## Congruence complete still means thin here

After the shared guard, S_j=181·2^j−154 and h=77, so gcd(S_j,h)=gcd(181,77)=1. The scalar terminal-zero transfer of the [fixed-cycle modular theorem](../chanpak_cross_problem_theorems_20261114/FIXED_CYCLE_CONGRUENCES.md) makes either arm reachable in every odd-compatible residue class modulo every modulus. This uses its scalar cut recurrence, not any extra left-end cycle-count claim.

For the exact integer targets through X, each guard j≥3 and each arm has output N_j+S_j t+77p(t) with t≥0 and p≥0. Thus it contributes at most 1+X/S_j presentations. Only O(log X) heights can appear, and
```
S_j ≥ (647/4)·2^j,      sum_(j≥3) 1/S_j ≤ 1/647.
```
Together both arms across every height and every common transitive tail therefore hit at most (2/647)X+O(log X) distinct absolute counts through X. This quantifies the remaining work even for a physically closed, all-height, congruence-complete seam.

Replay: compile `replay_h77_edge.cpp` with C++17. The canonical code is inherited from the order-nine census; the local signature and forward stop-rule replay are standalone here.
