# A taxonomy of LRC-to-tournament encodings: the faithfulness/restriction lattice (S521)

*claudebox-2026-06-01-S521, long creative session. A systematic survey of ways to
map the Lonely Runner Conjecture into "which iso-classes can a speed set exhibit?",
ordered by two axes — FAITHFUL (loneliness is a function of the iso-class) and
RESTRICTED (small realizable set). Computed with exact arithmetic across many
speed sets (some via parallel sub-agents). The honest meta-conclusion: the
difficulty is conserved — it migrates between the static iso-class set and the
dynamics — but the survey isolates where a proof should live.*

## The setup

Every encoding is a map `(speed set v, time t) -> iso-class of some structure`.
The *exhibitable set* `R(v) = { class(v,t) : t in [0,1) }`.  LRC(n=m+1) asks that
a designated *lonely class* lies in `R(v)` for every `v`.  Two failure modes make
an encoding **unfaithful** (unable to even state LRC):

- **Observer-blindness** — forgetting which vertex is the observer.  The runner
  half-turn tournament does this: its class is `A000016(m)` (round/locally-
  transitive tournaments), beautiful and tiny but it cannot see loneliness.
- **Single-modulus** — seeing only `t = k/n`.  The residue-orbit at `n` (the
  speed set's multiplicative orbit in `Z/n`) is `<= phi(n)` classes but ~**50%**
  of lonely sets are lonely *only* at other denominators (measured n=6,8,9), so
  it misses them.

A **faithful** encoding must determine every `||v_i t||` with the observer marked
— equivalently the full *gap structure*.  All faithful encodings are therefore
quotients of one finest invariant and refinements of one coarsest.

## The surveyed encodings (exact data, m=n-1 movers)

| encoding | faithful? | ambient | realizable / restriction | LRC criterion |
|---|---|---|---|---|
| **round tournament** (runner half-turn) | NO (obs-blind) | `A000568(m)` | `A000016(m)` (2,4,6,10,16,30,52) | — cannot state |
| **residue-orbit at n** | NO (single-modulus) | — | `<= phi(n)` | only `q=n` loneliness |
| **distance-rank tournament** (order by `||v_i t||`) | YES | `m!` | union saturates `~m!`; every runner blocks | a min-distance `>= 1/n` cell exists |
| **danger graph** (edge iff circ. dist `< 1/n`) | YES (obs-aware) | `2^{C(n,2)}` | tiny per-set (`~10^-3..10^-6`); partition-union unrestricted | observer's component is a singleton |
| **finest faithful** (round + tight-pairs + observer) | YES | huge | 40–231 (m=4..5) | observer 2-gaps loose |
| **gap-type necklace** (tight/loose of the `n` gaps, obs-marked) | YES | `2^n` | union `= 2^n-2`; per-set `~12..54` | both observer-gaps loose |
| **two-gap** `{LL,LS,SL,SS}` | YES | 4 | 4 | `LL` exhibitable |

(Loneliness for the extremal AP families is a **boundary/wall** event `t = 1/n`
with min-distance `= 1/n`; open-cell scans miss it.  All faithful counts above use
closed cells.)

## Two structural theorems found (rigorous)

1. **Danger graph = circular proper-interval graph.**  Its connected components
   are exactly the maximal arcs of circularly-consecutive points cut wherever a
   consecutive gap `>= 1/n`; within a component adjacency is "circular distance
   `< 1/n`" (a chain/indifference graph, *not* a clique — verified, 0 contiguity
   failures over 22256 cells).  **Observer lonely <=> its arc is a singleton.**
2. **Gap-necklace law: `#tight <= n-1` always**, and the only two forbidden
   length-`n` words are all-TIGHT (`n` gaps `< 1/n` cannot sum to 1) and all-LOOSE
   (forces every gap `= 1/n` exactly, measure zero).  So the realizable necklace
   union is exactly `2^n - 2`.  **LRC <=> a necklace with both observer-gaps loose
   is exhibitable.**

## The faithful-encoding lattice and the restriction<->provability tension

The faithful encodings form a lattice of quotients:
```
finest faithful  (round + tight pairs + observer)
        |                |                 |
   gap-necklace     danger graph      (other refinements)
        \________________|_________________/
                two-gap  {LL,LS,SL,SS}
```
all determined by "the cyclic order of the points and which gaps are `< 1/n`,
observer marked."  As you descend (more restriction):

- **Maximal restriction (two-gap, 4 classes)** makes LRC a one-bit statement
  ("`LL` exhibitable"), but the *static* exhibitable set is trivially
  `{LL,LS,SL,SS}` for every true-LRC set — **all the content moves into the
  DYNAMICS** (the transition structure of how the 2-gap state evolves with `t`,
  i.e. THM-386).  The encoding is faithful but the proof obligation is just
  relocated.
- **Maximal information (finest faithful)** carries everything but is too rich to
  reason about globally.
- **The sweet spot** is an encoding that is *both* restricted *and* retains a
  usable dynamics.  Two candidates:
  - **danger-graph coalescent**: the contiguous arcs merge and fragment as `t`
    varies; LRC = the observer's arc can be reduced to a singleton.  A
    fragmentation–coalescence process on a small state (the arc partition).
  - **gap-necklace flow**: the `n` gaps flip tight<->loose at walls under the
    `#tight <= n-1` constraint; LRC = both observer-gaps simultaneously loose.

So the survey answers "make the exhibitable set more meaningfully restricted":
the gap-necklace (`2^n`, `#tight<=n-1`) and danger-graph (observer-singleton) are
the most restricted *faithful* carriers, and the two-gap is the extreme — but
restriction buys crispness at the cost of pushing the difficulty into the
generating dynamics.  **Difficulty is conserved.**

## Multitudes of hypotheses (computed and conjectural)

- **HYP (gap-necklace realizability).** For a fixed `v`, the exhibitable gap-
  necklaces are exactly those consistent with the speeds' 3-distance gap spectrum;
  characterize them and the local rules that flip tight<->loose at walls.
- **HYP (coalescent monotonicity).** In the danger graph, track the observer's
  arc length as `t` varies; conjecture it always reaches a singleton (LRC) because
  the arc-endpoints (the two observer-adjacent runners) recede under the two-gap
  dynamics (THM-386).  A counterexample = a covering of `[0,1)` by "observer-arc-
  nonsingleton" cells.
- **HYP (faithful boundary).** An encoding is faithful iff its iso-class
  determines `{ i : ||v_i t|| < 1/n }` together with the cyclic order — i.e. iff
  it refines the finest faithful invariant's loneliness fibers.  Observer-
  blindness and single-modulus are the two ways to fall below this.
- **Abstract encodings to mine (uncomputed conjectures).**
  - *second-order pair-cell* tournament (vertices = unordered runner pairs;
    project's HYP-1976) — does it separate the hard locus better?
  - *even-graph dual* `E_n` (the project's dual object) of the round tournament —
    a parity encoding of loneliness.
  - *OCF conflict graph* `Omega(T)` of the round tournament — encode loneliness in
    the independence structure that computes `H`.
  - *voltage / doubling lift* (THM-378) — a signed encoding where the observer
    mark becomes a voltage.
  Each is a candidate faithful (or near-faithful) carrier; the question for each
  is its restriction and whether it exposes the blocker/fully-covered core.

## Honest meta-conclusion

No encoding cracks LRC; each reframes it, and the difficulty is *conserved* — it
sits in the static exhibitable set (rich encodings) or in the generating dynamics
(restricted encodings).  The value of the survey is to locate the sweet spot:
**a faithful, restricted encoding that still carries a dynamics** — the danger-
graph coalescent and the gap-necklace flow — where the residual obstruction is the
same one every S521 thread reached: free the observer's arc / both observer-gaps,
i.e. the blocker / fully-covered core, now dressed as a fragmentation process on a
partition (danger graph) or a constrained binary necklace flow (gap necklace).
