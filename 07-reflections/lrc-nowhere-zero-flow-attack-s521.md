# Attacking LRC from the nowhere-zero-flow mindset (S521)

*claudebox-2026-06-01-S521. Reformulating the Lonely Runner Conjecture in the
language of flows, tensions, divergence, and nowhere-zero flows on directed graphs.
Extends the sectors-as-nodes / chip-firing encoding. Connects to the project's
even-graph (mod-2 cycle space) framework.*

## 1. Flow/tension duality: positions are potentials, the threshold is a tension

Read the runner positions `v_i t mod 1` as **vertex potentials** on the complete
graph `K_n` (observer + runners), valued in `R/Z`. Each edge `(i,j)` carries the
**tension** `||v_i t - v_j t||`. The LRC constraint `||v_i t|| >= 1/n` is a
*lower bound on the tensions incident to the observer*. So:

> **LRC = a one-parameter family of tensions (the line `t |-> (v_i t)`) must, at
> some `t`, keep every observer-incident tension `>= 1/n`.**

Tensions are the flow-dual (Tutte): tension = coboundary of a potential, flow =
cycle in the kernel of the boundary. The positions live in the **cut space**
(potentials), the speed-relations in the **cycle space**.

## 2. The crossing-flow on the sector cycle; occupancy = divergence

On the directed `n`-cycle of sector-boundaries `b_k = k/n`, let `F_k(t)` = number
of runner-crossings of `b_k` during `[0,t]`. Conservation of runners gives the
identity
> **`o_k(t) = o_k(0) - (F_{k+1}(t) - F_k(t))`**, i.e. **occupancy = (initial) -
> divergence of the cumulative crossing-flow**; total divergence `= m` (the `m`
> runner-"sinks").

The flow `F` is **nowhere-zero** (every boundary is eventually crossed, all `F_k > 0`)
— that part is automatic. The LRC content is:

> **LRC(strict)  <=>  at some `t`, `div F = 0` at the two observer-vertices**
> (`o_0 = o_{n-1} = 0`): the runner-sinks vacate the observer's sectors, so the
> crossing-flow is *locally conserved* there.

So LRC is a flow with `m` units of prescribed divergence (the runners), asked to
become divergence-free at two designated vertices.

## 3. The observer as a nowhere-zero-flow SOURCE (the cleanest reframe)

In the THM-381 observer-marked tournament, `indeg(observer) = #{ i : ||v_i t|| < 1/n }
= N(t)` (the danger count). A nowhere-zero flow needs every vertex balanced
(in = out), so a **source** (in-degree 0) carries no circulation. Therefore:

> **The observer is a nowhere-zero-flow source  <=>  `indeg(observer) = N(t) = 0`
> <=>  the observer is lonely.**

The runner sub-tournament is bridgeless (a round tournament on `>= 3` vertices), so
it *always* carries nowhere-zero flows; the question is purely whether the observer
can be expelled from every circulation. **LRC = the observer becomes a pure
flow-source at some time.** As `t` varies the observer's flow-role cycles
source <-> internal <-> sink; LRC asks it to hit "pure source."

This is the flow form of "N(t) reaches 0", but it reframes the obstruction
crisply: a counterexample is a speed set whose observer is *perpetually on a
circulation* — never a source — i.e. the runner-flow always loops through the
observer. A "covering of the observer by circulations."

## 4. The relation lattice = the flow lattice of the speed system

The integer relations `{ a in Z^m : sum a_i v_i = 0 }` form a rank-`(m-1)` lattice
— the **flow space** of the 1-vertex multigraph with `m` loops of "lengths" `v_i`
(or the cycle space of the speed system). A **nowhere-zero** relation (all
`a_i != 0`) always exists for `m >= 2` (computed: `(1,2,3) -> (-1,-1,1)`,
`(2,3,5,7) -> (-2,2,1,-1)`). The minimal nowhere-zero relation is the "shortest
nowhere-zero flow" of the speed system; its `L1`-norm is a complexity measure of
the resonance structure. (No direct bound on lonely times found yet — flagged as a
lead: does the minimal NZ relation control the covering modulus?)

## 5. Mod-2 flows and the even-graph framework

A nowhere-zero `GF(2)`-flow = an **even subgraph** (every vertex even degree) = an
element of the cycle space. This is exactly the project's **even-graph metagraph
`E_n`**, the dual of the tournament metagraph `G_n`. So the NZF mindset is native
to this repo: the runner round-tournament's mod-2 flows are its even subgraphs, and
the OCF (which computes `H` from the odd-cycle collection) is the *odd*-flow side.
The observer-marked threshold edges are a `GF(2)`-coboundary (a cut) separating the
observer; loneliness = that cut is "all-out" (observer source).

## Honest assessment

The nowhere-zero-flow mindset gives a clean, unifying language and one crisp new
reframe — **LRC = the observer becomes a pure flow-source / the runner-sinks
clear the observer's two sector-vertices** — plus the divergence identity
(occupancy = `div` of the crossing-flow) and the connection to the project's
even-graph / cycle-space framework. It does **not** crack LRC: the surviving
obstruction is identical to all prior S521 threads — the observer's divergence
(`= N(t)`) must vanish, and a counterexample is a runner-flow that perpetually
loops through the observer (the flow form of the covering obstruction). What the
flow language *adds* is the toolkit: divergence/conservation, the cycle/cut
decomposition, the sandpile group of the sector cycle, and Tutte flow-coloring
duality.

## Creative hypotheses / leads

- **HYP (flow-coloring duality).** If a *planar* graph can be built whose
  nowhere-zero `k`-flows correspond to lonely configurations of an `n`-runner
  system, then `k`-colorability of its planar dual would give LRC. The sector
  `n`-cycle is planar but its cycle space is 1-dimensional (too thin); seek a
  richer planar carrier (e.g. the bad-arc incidence / interval structure, which is
  a circular-arc graph — when is it planar?).
- **HYP (source-expulsion).** LRC <=> for every primitive speed set the marked-
  tournament observer is expelled from all circulations at some `t`. Equivalent to
  `N(t) = 0`; the flow form suggests a Menger/connectivity attack: the observer is
  a source iff there is no directed path returning to it, i.e. no runner is within
  `1/n` — bound the "return paths."
- **HYP (minimal NZ relation controls the modulus).** The minimal nowhere-zero
  relation `sum a_i v_i = 0` has `L1`-norm `|a|_1`; conjecture the optimal lonely
  denominator (S521 Thm A: a pairwise sum) is bounded in terms of `|a|_1` / the
  flow lattice's covering radius — tying the flow lattice to the lonely time.
- **HYP (sandpile / abelian).** The runner-sinks moving on the sector cycle form a
  chip configuration; "clear the observer vertices" is a reachability question in
  the sandpile group `Z_n` of `C_n`. Is there a conserved sandpile invariant an
  adversary could exploit to keep the observer vertices charged forever?
- **HYP (mod-2 flow obstruction).** Express loneliness as a `GF(2)`-coboundary
  condition on the marked tournament (observer-cut all-out) and study it via the
  even-graph metagraph `E_n`; the odd/even-flow (OCF / even-subgraph) split may
  separate the lonely fiber from the observer-blind base.
