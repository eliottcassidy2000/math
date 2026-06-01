# Sectors as nodes: LRC as token chip-firing on a fixed circulant arena (S521)

*claudebox-2026-06-01-S521. A new family of LRC-to-tournament correspondences where
the NODES are the n evenly-spaced sectors of the circle (a fixed arena, independent
of the speeds) and the runners are a moving charge whose boundary-crossings drive
the dynamics. This is the program's addition/multiplication split made into a
single object: the sectors carry the multiplicative circulant structure, the
runners are the additive walk. Extends the encoding taxonomy.*

## The construction

Divide the circle into `n` sectors `S_k = [k/n, (k+1)/n)`. The observer sits at
the boundary between `S_{n-1}` and `S_0`. The runners occupy sectors with
occupancy `o(t) = (o_0,...,o_{n-1})`, `o_k = #{ i : v_i t mod 1 in S_k }`. The
forbidden zone (within `1/n` of the observer) is `S_0 U S_{n-1}` (up to the
boundary), so:

> **strict-LRC(n) for v  <=>  some t has `o_0 = o_{n-1} = 0`** (observer's two
> adjacent sectors cleared).

The nodes are FIXED (`n` of them, whatever the speeds); the speeds enter only
through the moving occupancy. As `t` increases a runner crossing a sector boundary
`v_i t = k/n` (`t = k/(n v_i)`) moves one token from `S_{j}` to `S_{j+1}`. So the
dynamics is **`m` tokens chip-firing forward on the `n`-cycle of sectors**
(verified: 40/49 consecutive occupancy-changes are single adjacent-token hops; the
other 9 are simultaneous crossings forced by speed resonances).

## The addition/multiplication split, as one object

- **Multiplicative arena.** The `n` sectors carry the half-turn circulant
  tournament `R_n` (the regular n-gon tournament). This is the fixed,
  multiplicatively-structured backdrop (the `A000016` / Paley / THM-369 side).
- **Additive walk.** The runners are tokens; their motion is the additive flow
  `t -> {v_i t}`. The "edges change when runners cross boundaries" = tokens hopping
  along the `n`-cycle edges.
- **Residue identity.** At `t = a/n` the occupancy is exactly the residue
  distribution `{ v_i a mod n }`, and clearing sectors `{0, n-1}` means
  `v_i a not in {0, n-1} (mod n)` for all `i` — the multiplicative-walk / THM-369
  picture, now read as token positions on the sector cycle.

So LRC becomes: **an additive token walk on the fixed multiplicative sector arena
must, at some time, vacate the observer's two forbidden nodes.**

## Faithfulness subtlety (strict vs closed)

Half-open sectors `[k/n,(k+1)/n)` capture STRICT loneliness (`||v_i t|| > 1/n`):
`o_0 = o_{n-1} = 0` iff no runner is *strictly* within `1/n`. The extremal/tight
sets (AP and friends) are lonely only at the boundary `t = 1/n` (distance exactly
`1/n`), so they are strictly NON-lonely — `o_0=o_{n-1}=0` never holds for them
(e.g. AP `(1,2,3,4)` at `t=1/5` has occupancy `(0,1,1,1,1)`: `S_4` holds the runner
at `4/5`, distance exactly `1/5`). So this encoding cleanly separates **strict
loneliness (open) from the boundary-tight extremizers** — the sector boundary is
exactly where the tight cases live. (A threshold-centred sector convention, or
tracking the boundary occupant, recovers closed loneliness.)

## Restriction

Realizable occupancy vectors per speed set: `~50/70`, `~123/252`, `~243/924`
(n=5,6,7) — moderately restricted per set, but the UNION over speed sets is
essentially all compositions of `m` into `n` parts. So unlike the gap-necklace
(globally restricted by `#tight <= n-1`), the sector-occupancy is restricted only
per-instance; its value is the FIXED-ARENA, dynamical (chip-firing) viewpoint, not
a global compression.

## A derived dynamic tournament

One can also build a genuine changing tournament on the `n` sector-nodes:
orient `S_i -> S_j` toward the less runner-crowded arc (clockwise vs ccw runner
count). This reorients as tokens move; realizable iso-classes are few (9 for n=5).
But it is a function of the occupancy and loses information, so it is not faithful
to loneliness on its own — a cautionary instance: a dynamic sector-tournament is
restricted but typically unfaithful, exactly like the residue-orbit and round
tournaments. The faithful content is the occupancy itself.

## Creative hypotheses (the multitudes)

- **HYP (chip-firing / abelian dynamics).** The token motion on the `n`-cycle is a
  continuous-time deterministic process with rates `v_i`. LRC(strict) = the
  configuration ever vacates nodes `{0, n-1}`. Is there an abelian-sandpile /
  rotor-routing invariant governing reachable configurations? A counterexample
  would be a token process that keeps `{0,n-1}` perpetually occupied — a "covering"
  of the observer nodes in token-time, the chip-firing form of the covering-system
  obstruction.
- **HYP (fixed-arena reduction).** Since the arena `R_n` is fixed and the speeds
  enter only via residues at `t=a/n` (occupancy = residues), strict-LRC at
  denominator `n` is purely the multiplicative-walk question on the `n` nodes; the
  full strict-LRC adds the inter-`a/n` token motion (the additive interpolation).
  Conjecture: strict loneliness, if it occurs, occurs arbitrarily close to a
  rational `a/q` with `q` a pairwise sum (S521 Thm A), realized as a token
  configuration vacating `{0,n-1}`.
- **HYP (boundary = tight locus).** The extremizers are exactly the speed sets
  whose token process touches but never vacates `{0,n-1}` (boundary-only), which is
  the chip-firing signature of `M = 1/n`. So the sector encoding draws the tight
  locus as a measure-zero touching set in token-space.
- **Other sector-node tournaments to mine.** (i) the boundary-crossing PARITY
  tournament (each boundary an edge whose state flips at crossings — a Z_2 spin
  system on the `n`-cycle); (ii) the flow/transfer tournament (net token flux
  across each boundary); (iii) the threshold-centred sector convention giving
  closed loneliness directly.

## Assessment

The "sectors as nodes" idea yields a genuinely new correspondence: LRC as a
**fixed-arena token chip-firing** problem (`m` runners hopping on an `n`-cycle,
clear the observer's two forbidden nodes). Its strengths: the arena is fixed and
small (`n` nodes), the addition/multiplication split is explicit (circulant arena +
additive token walk), and the residue/THM-369 picture is recovered at `t=a/n`. Its
limit: faithful only to STRICT loneliness (the boundary-tight extremizers escape),
and globally restricted only per-instance. It does not crack LRC, but it adds the
**chip-firing / token-dynamics** viewpoint to the toolkit — a setting with its own
deep invariants (abelian sandpile, rotor-routing) that the runner-based encodings
do not expose, and where the obstruction reappears as "perpetual occupation of the
observer nodes," the token form of the covering obstruction.
