# Five strong rank-preserving steps: all 512 odd classes modulo 1024

Let `p=(0,0,3,5,2,4)` be a literal span-3 ribbon word, counting its first zero as the first vertex. The seed is strong, has H=25, and its prefix product in the 61-dimensional literal model after the first vertex has **rank 4**.

**Finite exact theorem.** For every odd r modulo 1024 there is a five-letter continuation using only nonzero masks such that the final count is r modulo 1024. Every intermediate tournament remains strong and each prefix product retains rank 4. The same physical prefix works for all 512 target classes.

This is a rank-preserving statement, not just a set of endpoint counts. The certificate lists a five-letter word for each odd residue in increasing order. The verifier rebuilds the rational observable and reachable spaces O,B from the [literal model audit](https://github.com/eliottcassidy2000/math/blob/2e1da5f6e4fba682fc62e326085890e5f3cbccf4/05-knowledge/results/chanpak_target_ribbon_audit_20261116/RANK_BUDGET.md), proves the seed's exact rational rank 4, and extracts four selected columns of the compressed observable image. For each listed continuation it propagates those four raw columns and certifies a nonzero projected rank-four minor modulo the prime 10007. This is a one-sided exact proof: the nonzero modular minor implies rational rank at least 4; the rational upper bound 4 follows from multiplication by linear maps on the same quotient. Rank cannot increase, so equality at the end certifies equality at *every* intermediate time. Modular row elimination only uses small integers and modular inverses of nonzero pivots, not floating point.

Strongness at each step follows from the physical geometry: each nonzero mask makes the appended vertex beat at least one old vertex, and some older vertex more than three places back beats it; the old tournament is one SCC. Count residues use the exact frontier transfer. Seventeen residue representatives also have their literal Hamiltonian counts checked by a separate subset-path DP.

Run `python check_rank4_steering.py` (NumPy and SymPy; keep the sibling target/ribbon audit directory when checking out). This replays `RANK4_MOD1024.json` rather than trusting a sampling estimate.

## Scope

Rank retention need not return the image space, just as keeping three independent leash directions does not return to the original tangle. For example the one-letter seed fan H=25 at word (0,1,3,3,2,5) does spend rank for half of its choices; a better seed and multiple letters are needed. The new five-letter certificate makes no claim that it can be iterated from its new states, that it hits an exact prescribed integer, or that it supplies a uniform lower growth certificate. An actionable next test is an input/output family of *image states* on which one can repeat this modular control without losing the physical growth cap or the absolute target band.
