# Ever more complex vertices: LRC as a cutting sequence, and what the flow/Menger views add (S521)

*claudebox-2026-06-01-S521, long creative session. Pushing the "what are the
tournament's vertices?" question to richer, more precise objects, pursuing the
flow-coloring and source-expulsion/Menger threads, and connecting LRC to symbolic
dynamics. Extends the encoding taxonomy and the nowhere-zero-flow attack.*

## Source-expulsion / Menger: a clean NEGATIVE

In the THM-381 marked tournament `indeg(observer) = N(t)` (#runners within `1/n`),
and the observer is a nowhere-zero-flow SOURCE iff `N(t)=0` iff lonely. I tested
whether Menger / connectivity adds anything: it does **not**. "Observer is a
source" is the *local* degree condition `indeg = 0`; there is no global
connectivity obstruction (whenever `N >= 1` the observer trivially lies on a
directed cycle). So the source/flow language re-expresses `N(t)=0` crisply but the
Menger/connectivity machinery has no purchase — a useful dead end to log.

## The chamber sign-vector = the cutting sequence of the torus line (the rich vertex)

Make the vertex the FULL **sector-vector** `s(t) = (floor(v_1 t n), ..., floor(v_m t n))
in (Z/n)^m` — which sector each runner occupies (refining occupancy, which only
records how many per sector). The curve `t -> (v_i t)` is a line on the `m`-torus;
the `1/n`-grid cuts it into cells, and `s(t)` is the **cutting sequence** of that
line. Properties (computed, exact):

- **Polynomially thin.** The line visits `~ n*sum(v_i)` distinct sector-vectors
  inside the `n^m` grid: `30/625, 50/625, 60/7776, 144/7776, 84/117649` — ratio
  `-> 0` (polynomial cells in an exponential grid). A genuinely RESTRICTED yet
  FAITHFUL (to strict loneliness) encoding.
- **Cutting-sequence dynamics.** Consecutive cells differ by a single coordinate
  `+-1 (mod n)` (a runner crossing one sector wall) — verified 40/49 single-
  coordinate moves; the rest are simultaneous resonant crossings. So the encoding
  is a walk on `(Z/n)^m` along grid edges.
- **LRC = central box in the cutting sequence.** strict-LRC `<=>` the sector-vector
  with all coordinates in `{1,...,n-2}` (the safe band `[1/n,1-1/n]^m`) appears.
  Tight/extremal sets miss it (boundary-only), exactly as before.

This is the **symbolic-dynamics view**: LRC asks whether the cutting sequence of a
rational torus line contains the "central" letter. It ties LRC to the factor
structure of multidimensional Sturmian / billiard sequences and to discrepancy
theory (next).

## Even more complex vertices (a hierarchy)

- **Rauzy-graph / factor vertices.** Take the vertices to be the FACTORS (subwords)
  of the cutting sequence of length `k` — the Rauzy graph `R_k`. Its size is the
  factor-complexity `p(k)`, polynomial (`~ k^{m-1}` for an `m`-torus line). LRC = a
  central-box factor appears in `R_k` for some `k`. This makes the vertex set the
  symbolic *complexity* of the orbit — precise and deep; the appearance of a given
  factor is governed by the line's continued-fraction / Diophantine data.
- **Bad-arc-endpoint vertices.** The `2*sum(v_i)` endpoints of the danger arcs on
  the time-circle; the depth function `N(t)` is the covering multiplicity between
  consecutive endpoints; LRC = a depth-`0` gap. A 1-dimensional "Davenport-Schinzel"
  / arrangement structure.
- **Discrepancy / bounded-remainder vertices.** The central box is (for the right
  shape) a *bounded remainder set* of the torus rotation; the number of orbit
  points it contains is `measure * #orbit + bounded error`. LRC = the box is hit;
  the obstruction is the Diophantine discrepancy of the finite orbit
  `{ (v_i a/q) }`. Vertices = the orbit points / the discrepancy excursions.
- **Poincare-section / first-return vertices.** Cross-section the torus flow by a
  hyperplane; vertices = return-map states. The return map is an interval-exchange
  / rotation; LRC = a return into the central box.
- **Relation-lattice vertices.** Points of the rank-`(m-1)` flow lattice
  `{ a : sum a_i v_i = 0 }` (the speed system's cycle space); the minimal
  nowhere-zero relation is the shortest nowhere-zero flow (S521 NZF reflection).

The pattern: as the vertices get more precise (occupancy -> sector-vector ->
factors), the realizable set stays polynomially thin and the LRC criterion stays a
clean "central object appears," while the proof difficulty migrates into the
Diophantine/symbolic-complexity data of the line — the same conserved difficulty.

## Assessment

The richest new vertex is the **sector-vector / cutting sequence**: polynomially
restricted in an exponential ambient, faithful to strict loneliness, with a clean
single-coordinate-walk dynamics, and it places LRC inside **symbolic dynamics**
(cutting sequences, factor complexity, bounded remainder sets, discrepancy) — a
toolkit not previously brought to bear here. The Menger/source-expulsion thread is
a logged negative (source = local degree). The flow-coloring and pair/witness
tournament threads are under separate computation; their results extend this
taxonomy. Difficulty remains conserved — it now lives in the Diophantine data
controlling whether the central factor appears — but the cutting-sequence /
discrepancy framing is the most principled new home for it after the covering-
systems route.
