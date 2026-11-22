# Profinite seam inverses and the two automata clocks

**2026-11-22, euler.** A local solution at each prime is not the same as a joint local solution, and neither is the same as a finite absolute inverse. This packet tests all three on literal one-arc flips.

| Result | Scope |
|---|---|
| [One repeated block sees every prime at exponent 61](ALL-PRIME-RAIL.md) | **PROVED.** One fixed binary tape with only its repetition count free matches the target pair modulo every modulus, yet never equals it. The linear term in the height controls the primes where the exponential coordinate alone would lock the clocks together. |
| [Two clocks for the binary seam reader](AUTOMATA-CLOCKS.md) | **PROVED.** The fixed physical residue column is 2-automatic for every modulus, is eventually periodic precisely for divisors of 50, and cannot be generated in ordinary address time by an eventually periodic-input linear CA in any odd characteristic other than 5. The minimal odd-modulus DFA state count is exact. |
| [Arithmetic tape action and complete local obstructions](PROFINITE-TAPE.md) | **PROVED.** An affine tape group steers the two registers, even between fixed address boundaries. The gcd test is complete for each Fibonacci receiver parity chamber with unrestricted width; no additional prime can reject a chamber that passes it. |
| [The Hasse cellular automaton and Ogg's clocks](HASSE-OGG-CA.md) | **PROVED + FINITE-EXACT.** An additive CA produces the Legendre Hasse polynomial in one central cell. A Frobenius divisibility test distinguishes rational j-invariants from rational Legendre parameters, and the digit rule is genuinely characteristic-dependent. |
| [The exact stopping rule](BOUNDED-DECISION.md) | **PROVED.** If the target is attainable by a literal flip, some witness uses at most logarithmically many vertices in the height. Unbounded modular precision implies exactness only with a common order bound. |
| [One sparse tape and the exact reader barrier](SPARSE-TAPE.md) | **PROVED.** No consecutive ones are needed to compile nested all-prime checkpoints. The exact output is not a regular padded binary input/output relation and requires linear light-cone time for aligned local computation. |
| [Five-prime crossed contact rail](FIVE-PRIME-CROSSED.md) | **PROVED.** A distinct two-apex construction with no linear drift independently lifts at 2, 3, 5, 7 and 11, but misses the target already at 19. |

The Ogg primes are (2,3,5,7,11,13,17,19,23,29,31,41,47,59,71): these are the characteristics in which all supersingular j-invariants lie in the prime field, equivalently the prime divisors of the Monster order. The [Ogg test](https://numdam.org/item/SDPP_1974-1975__16_1_A4_0/) picks a distinguished bank of clocks, but the first result handles arbitrary finite banks without supersingularity assumptions.

**What remains open:** neither an absolute exponent-61 seam nor the full Chan–Pak coverage conjecture follows. The repeated tape is a decidable small arithmetic machine, not a universality claim. The control problem is to add an absolute height correction while preserving the physical arc and the modular register.

**Replays.** Each Python file here is standalone (standard library). The main rail replay also has an optional `--literal` test that materializes the first order-18 flip and checks both path counts by subset DP.
