# Session67 — enumerate tournament classes from the complete metagraph carrier

**First verified checkpoint, 2026-10-26.** The local-relocation lane is audited.
Weighted-core recurrence, information-history integration and their corrections
remain ongoing and are not results of this checkpoint.

## What is established

- **PROVED corollary of inherited all-order connectivity:** the full fixed-path
  cube, partitioned by legal single-vertex relocation, enumerates ordinary
  tournament classes without importing their number. Connected blocks plus
  closure under every legal edge give a count-free recognition certificate.
- **FINITE-EXACT through order 8:** the implementation produces every ordinary
  class representative and its normalized tiling multiplicity, then recovers
  converse-merged nodes and all coloured complement lines, including loops
  and repeated lines. Order 8 has **6,880 ordinary classes**, **3,528 merged
  nodes**, and **1,048,576 literal lines** stored as 438,374 weighted endpoint rows.
- **PROVED linear algebra:** the relocation incidence quotient has one basis
  dimension per ordinary class over every field. The invariant-function
  algebra is a separate object. The mod-two graph Laplacian is not a substitute:
  a smallest order-4 `K_(2,3)` component has incidence rank 4, Laplacian rank 2.

| Order n | Ordinary classes A000568(n) | Converse-merged nodes | Self-converse classes |
|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 |
| 1 | 1 | 1 | 1 |
| 2 | 1 | 1 | 1 |
| 3 | 2 | 2 | 2 |
| 4 | 4 | 3 | 2 |
| 5 | 12 | 10 | 8 |
| 6 | 56 | 34 | 12 |
| 7 | 456 | 272 | 88 |
| 8 | 6,880 | 3,528 | 176 |

Order zero is the empty-tournament convention, not a program run. Orders 1–2
have no complement lines under this implementation's convention.

## Evidence and reproduction

- [Statement, proof, guards and information ledger](session67-root/COUNT_FREE_ENUMERATION.md).
- [Reproduction commands](session67-root/REPRODUCE.md), [C++ source](session67-root/enumerate.cpp),
  [exact census](session67-root/output/census.json), and [artifact verifier](session67-root/verify_exports.py).
- [Independent audit](session67-local-audit/AUDIT.md).
- [Frozen-source provenance](PROVENANCE.json) and [publication gates](PUBLICATION_CHECKS.md).
- `python3 verify_packet.py` or `python3 -O verify_packet.py` verifies the exact
  sealed file set and SHA-256 hashes. `--repo /path/to/math` also verifies inherited
  source bytes against their recorded commit.

## Boundaries

The algorithm is all-order mathematically, but the supplied implementation is
explicitly bounded to n<=8 and uses exponential space. It is not a polynomial-time
class-count formula, a black-core connectivity theorem, or a claim that coarse
metagraph statistics determine classes. Multiplicity is `tc(T)=H(T)/|Aut(T)|`, not
raw Hamiltonian-path count H. No novelty claim is made for inherited connectivity.
