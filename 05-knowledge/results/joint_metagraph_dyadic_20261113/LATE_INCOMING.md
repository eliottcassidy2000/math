# Final concurrent merge: slot equality and exact elimination

The last bounded incoming pass is pinned at
`5479c530c7907b5442dd885cbd6b47d76b41b8ba`. This note records incoming
interfaces and their compatibility with the present results. It does not
adopt an unreviewed larger census as a new proof in this packet.

[THM-5153, exact slot pairings and self-converse structure](../../../01-canon/theorems/THM-5153-exact-slot-pairings-of-merged-tournament-classes-and-self-converse-structure-to-order-nine.md)
and its [packet](../kps_s28_exact_pairing_sc_census_20261113/README.md)
use the word pairing for **equality of a full slot polynomial** between two
merged classes. Their source status is FINITE-EXACT through order nine;
the incoming commit is labelled pre-audit, and this packet does not rerun
that full census. This is a different relation from an involution of all
ordinary classes preserving the physical dart graph.

The slot polynomial is
`Q_T(t)=sum_j F_j(T) binom(t,j)`, where F_j counts ordered covers into j
nonempty directed paths. The inherited [uniform-context theorem](../session310_positive_clocks_20261113/NEUTRAL_TAILS.md)
says equality of Q permits replacement in every uniformly attached exterior
context by retaining the number of module runs. It explicitly does not
preserve selective internal ports. Our [independent small slot census](SLOT_WIRING.md)
tests the resulting interface with actual complement-line wiring rather
than treating the larger incoming enumeration as an audited input.

[Session316's squarefree port elimination](../session316_hinge_and_signed_repairs_20261113/SQUAREFREE_PORTS.md)
provides a useful exact counterpart to our walk-memory Schur complement.
It works over the ring `Z[x_i]/(x_i^2)` on interior vertices, with effective
port matrix

```text
E_A = A_PP + A_PI X_I (I-A_II X_I)^(-1) A_IP.
```

The nilpotent monomials retain interior support; products kill reuse of a
vertex. For distinct endpoints, each monomial counts a path with that exact
interior. The inverse is a finite formal series. Substituting every x_i=1
is not a ring map, since x_i^2=0 would become1=0. This explains the missing
coordinate in the source's positive-reachability but opposite-pairing
witnesses. In our recurrence, retaining hidden state through K(z) similarly
repairs loss under conditional averaging; the two constructions use
**different algebras and different observers**. Their common exact operation
is elimination with a retained response, not equality of their outputs.
The incoming two-sign determinant channel additionally retains cycle parity.
Its all-seams A5 condition concerns permutations of five closing seams,
not an A5 action on the tournament's vertices; it therefore does not violate
our signed-action obstruction.

[Session316's cyclic wrap](../session316_hinge_and_signed_repairs_20261113/CYCLIC_WRAP.md)
recovers the full slot operator
`Q_W(z)=z^2 Q_T(z)+z(Q_T(z+1)-Q_T(z-1))` for `W=C3[T,1,1]`.
Its incoming ray result has self-converse witnesses for every `7*3^j`
with j>=3, seeded by strong counts189 and567; ordinal propagation is not
strong. The exclusions7,21,63 have the earlier finite-target cutoff proof.
This is compatible with our H<=21 defect formula and with dyadic residue
coverage, and supplies no universal exact-target selector.

[Sequential reentry](../session316_hinge_and_signed_repairs_20261113/SEQUENTIAL_REENTRY.md)
retains a fixed pattern of revealed local probes in a neutral mirror bank.
Unrevealed completions stay neutral, with binomial fibre size, while a
specified reading order separates successive bit values. The source's
proof preserves the named apex and all already revealed positions. It
does not make its branch magnitudes arbitrary. This is a positive example
of carrying precisely the port data which a scalar count would discard.

Root read the elimination, wrap and sequential arguments in this bounded
pass. Their standalone replays and large literal mirror-bank computations
were not rerun here; their statements remain credited to the source packet.
The fresh slot/wiring census below is independently reproduced with its
own finite universe and hashes. No new all-order spectrum claim is imported.
