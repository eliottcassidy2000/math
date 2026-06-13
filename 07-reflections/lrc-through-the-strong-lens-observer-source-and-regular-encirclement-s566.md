---
source: opus-2026-06-02-S566 (remote-control)
status: SYNTHESIS — the strong (strongly-connected) lens on LRC; observer-source reformulation, generic strong encirclement, tight = regular strong block; verified
tags: [LRC, strong-tournament, strongly-connected, SCC, observer-source, Moon, Hamiltonian, regular-tournament, n14, THM-381, S525, S565]
---

# The LRC through the strong lens: the observer escapes to a source, encircled by a strong block

**Prompt (user):** focus on strongly connected tournaments and their relationship
to the LRC; see through the strong lens.

## 1. The reformulation the strong lens gives

**THM-381:** LRC(n) ⟺ for every primitive speed set, at some time the observer
(vertex 0) is a **source** of the (observer + runners) half-turn tournament. A
source is its own top SCC. So:

> **LRC ⟺ the observer can always be lifted to the UNIQUE TOP SCC (sole source)
> of the condensation at some time.**

The condensation of any tournament is a transitive chain of strong components. As
`t` varies, the observer's SCC-membership oscillates: when a runner is within `1/n`
of `0`, the observer is swallowed into a strong block with it; LRC says the
observer **escapes all strong blocks to source-hood** at least once. A
counterexample = an observer **permanently entangled** — never the sole source.

## 2. What the runner block looks like when the observer IS a source (verified)

When the observer is a source, the runner sub-tournament (S525/HYP-2000) has
**#SCC ∈ {1, m}** — nothing in between. Verified on 152 n=14 configs
(`lrc_strong_lens_s566.py`):

- **#SCC=1 (a single STRONG block): 139 (91%)** — the generic case;
- **#SCC=13 (transitive): 13 (9%)**, exactly the configs whose runners fit in a
  **semicircle** at that time (13/13);
- **intermediate: 0** — confirmed.

So **loneliness is generically realized by a single strong block**: the runners
**encircle** the observer (wrap fully around, not a semicircle), and the observer
is lonely because it sits in a **gap** of that encirclement. The minority
transitive case = runners clustered in a half-circle, observer lonely on the
empty side.

## 3. Moon's theorem: strong = the runners Hamiltonian-encircle the observer

A tournament is strong **iff** it has a Hamiltonian cycle (Moon: strong ⟹ vertex-
pancyclic). So the runner block being strong ⟺ the runners admit a **Hamiltonian
cycle in half-turn order** ⟺ they **cyclically surround** the observer. The
picture of loneliness-via-strong is exact:

> the observer is **surrounded by a Hamiltonian cycle of runners** yet sits in one
> of its gaps (a `≥ 2/n` opening).

## 4. The tight cases = the REGULAR strong block (and this closes the circle)

The AP and V* are lonely (boundary `t=1/14`) via a **strong** block that is the
**regular rotational tournament**: every runner beats the next `(m-1)/2=6` in the
cyclic order — the unique **doubly-regular** (all out-degrees equal) strong
tournament, the roots-of-unity encirclement. Generic strong configs use
**irregular** strong blocks (uneven out-degrees). So within "strong," the worry-set
is the **maximally regular** encirclement. This is exactly the dual-Burnside
**fix-side** object (S565): the regular rotational tournament is **self-converse**
and rotation-symmetric. The threads converge:

| lens | easy / IGNORE | hard / WORRY |
|---|---|---|
| measure (S564) | positive | zero |
| orbit (S563) | low resonance | resonance-maximal |
| Burnside (S565) | orbit side (round/A000016) | fix side (self-converse) |
| **strong (this)** | **transitive *or* irregular strong** | **REGULAR strong (rotational encirclement)** |

## 5. The strong-lens proof shape

> Prove: for every speed set, the observer reaches sole-source at some time —
> equivalently, the runners can always be brought to **encircle (or semicircle-
> cluster) the observer leaving a `≥2/n` gap**. The only sets where this is in
> doubt are those forced into a **regular rotational encirclement** (the tight
> family); for everything else an irregular strong block or a semicircle does it
> with room.

The hard core is thus: **can a regular rotational encirclement always leave the
observer a gap?** For the AP it does (the `1/14` wall). The conjecture is that no
speed set forces the regular encirclement to close the gap entirely — i.e. the
**regular strong block never becomes a perfect cover of the observer's
neighborhood.** This is the strong-tournament face of "the danger arcs never
completely cover" and of the dual-Burnside "self-converse boundary classes are
always source-realizable."

## 6. Honest status

The observer-source reformulation (THM-381) and #SCC∈{1,m} (S525) are known; the
verified **generic-strong / transitive-semicircle split** and the identification
**tight = regular rotational strong block = the dual-Burnside fix-side / Hamiltonian
encirclement** are the strong-lens synthesis. It unifies the measure, orbit,
Burnside, and strong views onto one object — the regular rotational encirclement —
and states the proof target in strong-tournament language. Not a proof.

**Artifacts:** `04-computation/lrc_strong_lens_s566.py` (+`.out`). Builds on
THM-381 (observer source), S525/HYP-2000 (#SCC∈{1,m}), Moon's theorem, S565 (dual
Burnside), S563/S564 (orbit/worry). New: **HYP-2089**.
