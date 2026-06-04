# Unit-Distance Impairment Atlas S623

The useful shift in this session was to stop asking only how to improve the
unit-distance beam and instead ask how to break it. Small failures are cheaper
than large successes, and they identify the side channels that a proof or
construction has to preserve before scaling.

The inspirations came from other repo threads. LRC quotient work taught that
compression only becomes proof-relevant when owner, carry, pinch, and lift data
are reattached. Tournament work taught that scalar counts hide packet
structure. The adversarial cauldron sessions taught that changing the schedule
can reveal a resource the base game forgot. Equidecomposability sharpened the
distinction between equal count and predicate-preserving class. S617 then gave
the concrete unit-distance tool: edge counts can be made state-local through
frontier gains, but that local observable is not the whole proof.

Incoming HYP-2194 already gives a sharp Moser spectroscopy lane: direction-pair
dropout, observed support masks, and gain caps through `n=14`. This S623 atlas
keeps that as signal and asks a wider impairment question: what do width,
policy, triangular controls, and canonicalization tell us that direction
dropout alone does not?

S623 therefore builds a small impairment atlas for triangular and rank-4 Moser
carriers. It deliberately damages width, ranking policy, direction support,
gain ceilings, and canonicalization, then measures what changes.

The strongest small signal is width on the Moser carrier. At target `14`, widths
`1`, `3`, and `10` find only `28` or `29` edges, while width `30` reaches `33`.
The high-width state has a tighter span and a future frontier containing gain
`3` and gain `4` moves. Width is not merely sampling volume here; it preserves
future high-gain shape.

Policy impairment separates scalar edge count from usable geometry. On the
Moser target `12`, edge-only ranking still finds `27` edges but lets the span
expand to `10`, while sprawl-bias drops to `24`. Compactness is a retained
construction invariant. The triangular carrier does not show the same policy
sensitivity at this size, which is also useful: the atlas distinguishes robust
toy carriers from fragile Moser-like carriers.

Direction-drop jackknife is the most proof-like new tool. Dropping any
triangular antipodal direction still reaches the same best small count, but six
of nine Moser antipodal directions at target `10` lose one edge when removed.
That gives a certificate shape: dense Moser candidates should carry a
direction-support loss vector, and a proposed extension should explain which
directions are load-bearing.

Gain ceilings make high-gain packets explicit. On Moser target `12`, ceilings
`1`, `2`, `3`, and `4` give best counts `12`, `21`, `25`, and `27`. For the
S614/S617 `n=22` frontier, this is exactly the missing object: a `61`-edge graph
has to be a dense `21`-core plus a high-enough gain extension, with geometry and
obstruction side channels still attached.

Canonicalization also belongs in the ledger. The triangular target `10` test
shows D6 canonicalization reduces last-layer children from `1797` to `1428`
with the same best count. The ratio is modest at small size, but it is a
measured duplicate-work channel rather than a vague implementation preference.

The technique program coming out of the session is:

1. damage-response ledger,
2. direction jackknife certificates,
3. gain-threshold extension solver,
4. deletion-core resilience score, and
5. orbit-budget accounting.

For the `n=22` problem, the immediate next object should be a side-channel
complete extension ledger. Each retained dense `21`-core should carry its
candidate gain packet, direction-support profile, deletion-resilience vector,
canonical orbit class, and totally-unfaithful obstruction status. Then the
search question becomes narrow: does any side-channel-complete gain-`4` or
gain-`5` extension lane survive?

Tournament Analysis uses impairment lenses as vertices rather than points or
unit edges. The S623 route tournament is transitive and ranks extension gain
ledger, direction-drop jackknife, and canonical-orbit repair above raw wider
beam. That is the right bias for future work: widen only after the small
impairment atlas tells us which invariants the width is preserving.
