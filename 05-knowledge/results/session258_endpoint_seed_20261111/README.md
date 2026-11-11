# Session258: exact endpoint constraints and odd ray seeds

**Research completed 11 November 2026. Chan–Pak and universal odd-slope seed existence remain OPEN.**

Start with [the synthesis](SYNTHESIS.md). The central results move from a gcd condition to exact constructive endpoint conditions and a complete finite decision for every fixed ray multiplier.

- [Strong-core seed classification](seed_decision/CORE_CLASSIFICATION.md): all unit-offset seeds have one strong selected core plus unselected transitive padding; fixed-multiplier search becomes finite.
- [Endpoint cut theorem](endpoints/PROOF.md): exact positive cut image, heavy-arc constraints, and the sharp one-vertex+2 obstruction.
- [Nine-state positive compiler](seed_decision/PROOF.md): complete fixed-q decision, all-length constructor, and an independent finite certificate.
- [Zeckendorf predecessor and carry bridge](grammar/HISTORY_BRIDGE.md): 27-coordinate suffix requirements, legal inverse guards, and exact native carries.
- [Positive seed probes and machine interpretation](grammar/PROOF.md): two actual tail counts replace a signed equation; least-resistance and TM/CA statements retain their precise quantifiers.
- [Carmichael, Catalan, and99·2^k+23](arithmetic/RESEARCH.md): source fibers, a whole-chamber exclusion, and the exact first25-prime phase sieve.
- [Incoming strong zero digit](INCOMING_SYNTHESIS.md): the universal+2 theorem closes the entire final one-birth repair, beyond the incoming finite cut test.
- [Open proof obligations](OPEN_QUESTIONS.md), [source pins](SOURCE_PINS.json), [reproduction](REPRODUCE.md).

The odd slope c and actual Hamilton total N are different. Slopes7 and21 are valid; totals7 and21 are excluded. A fixed multiplier q=c2^j is a third quantity. The inherited seed bank mixes selected-set tail roots and doubling vertices; absence within one mechanism does not exclude the other.

New results are labelled proofs in their notes or finite exact computations. Classical bounds and the repository's earlier constructions are inherited at their stated scopes, not re-audited wholesale. No canon promotion or literature-priority claim is made.

**Executable constructor:** `python3 construct_seed.py --q 20` emits a slope5 seed and its first actual tail tournament. `--q 7` returns a complete exact-multiplier exclusion. See [the command guide](REPRODUCE.md).
