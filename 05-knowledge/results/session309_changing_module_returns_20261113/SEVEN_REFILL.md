# A seven-prime neutral return that increases actual trailing capacity

This is the same Π frame and outgoing-neighbor mask convention as [PHYSICAL_REFILL](PHYSICAL_REFILL.md). Unlike that witness, both replacement bases remain nontransitive, and the exact target is 7 times a prime.

| Slot base | Row masks | H | Relevant cover count |
|---|---|---:|---:|
| Old left M | 6,12,80,37,11,23,59 | 101 | P(3)=8823 |
| Old right N | 2,12,145,165,203,151,175,3 | 131 | P(2)=2864 |
| New left L | 2,20,9,35,45,7 | 29 | P(3)=1845 |
| New right R | 2,4,9,51,135,215,159,15 | 59 | P(2)=1712 |

The old right slot is N→TT₂ and the new right slot is R→TT₅, with ordinal joins. All four base modules are strong. The equality is

    60·8823·2864·4 + 101·131
    = 60·1845·1712·32 + 29·59
    = 6064590511
    = 7·866370073,

where 866370073 is prime. Both right bases are strong, so the actual available trailing sites rise from two to five. The intrinsic two-slot valuations rise from six to nine. Keep the seven exterior labels and all their mutual arcs in the rebase, replacing the two modules as units; no induced inclusion of old module interiors is claimed. The new five physical tail transfers retain both new module interiors, exactly as in the prime construction, by putting each moved singleton at the front of the left slot and reversing its arc to c. Their scalar law is H↦(3H−1711)/2. This is a one-time capacity-increasing exact-target return, not a proof that the new terminal count permits a second such rebase.

A search used the [inherited response bank](../target_collatz_transitions_20261113/response_bank.json.gz), which keeps one row per pair (H,P(2)); recomputing P(3) gives a supply of actual witnesses, not a complete (H,P(3)) bank. Only the explicit verified witness is claimed, not a negative enumeration result. The initial exploratory search incorrectly partitioned an already occupied set instead of its complement when counting P(3); an independently written verifier caught the error. No earlier exploratory hits were promoted. The included clean verifier recomputes all literal module rows and five named graph edits without depending on search code.
