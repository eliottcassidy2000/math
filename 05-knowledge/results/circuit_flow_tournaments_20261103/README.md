# Tournament circuits and information flow

**PROVED scoped; FINITE-EXACT controls.** Same-agent review with separate
counting and native-move implementations. Builds on S176's actual strong
edge gates and the older path-cover substitution calculus.

The two gates now have complete three-shift operators. Their composition
order has an exact meaning:

    H(L[U[A]]) - H(U[L[A]]) = -72 pc_3(A).

Here pc_3 counts three unordered disjoint paths covering the actual input.
Order is invisible exactly for inputs of at most two vertices. Starting
with one vertex gives two nine-vertex sources with H611, connected by
four legal native complements while retaining every old name and the
source root star. A common attachment reads their hidden difference as
88693 versus86317. At the physical root the neutral reset changes F by44.

For arbitrary input tournaments, at most d future gates read exactly a
`d+1` dimensional linear space of path-slot information. For the singleton
word grammar there is also a sharply different finite quotient: **216 is
the largest modulus hiding every gate reordering**, with a minimal
19-state output clock. The mod27 clock is
`H = 2^(i+12i^2), i=#L+14#U mod18`; a parity bit gives the mod8 factor.

The 81/27 pair has equal complete profiles and remains invisible to every
uniform gate word. A common nonuniform one-vertex attachment distinguishes
it. Full source data is retained; no arbitrary joint-port quotient or full
tournament spectrum theorem is claimed.

- [Proofs and exact formulas](PROOFS.md)
- [Recovered older work and concept board](RECOVERY.md)
- [Audit universes and reproduction](AUDIT.md)
- [Next decisive test](HANDOFF.md)

From the repository root, with Python and SymPy:

```sh
python -X utf8 -B 05-knowledge/results/circuit_flow_tournaments_20261103/reproduce.py
```

Four scripts run normally and under `-O`, with identical exact receipts.
