# Session mac-mini-2026-10-20-kepler: the lonely runner on a regular polygon, Kepler's second law, the polygonal order, and the colouring bridge

**Historical session record; corrected 2026-10-20.** The current THM-4737
uses per-runner one-sided rates at closed endpoints, allows infinite
polygonal order, and separates a real-orbit supremum from attainment.
HYP-9350's recalled fractional-colouring equality is **REFUTED** by
D={1,3,4,7}: chi_f=4 but1/kappa=5. The surviving bound retains the rotation
constraint. [Correction and exact proofs](runner_area_clock_20261020j/boundary_and_colouring.md).
The original narrative below is provenance, not the current truth source.

Original status: pre-audit. Canon: THM-4737 (PROVED elementary + FINITE-EXACT), HYP-9350 (OPEN). THM-4738 was
reserved and is ceded unused. Experiments: `04-computation/experiments/kepler_lrc_20261020/`.

## 1. The owner's prompt and the repo state

"Attempt the lonely runner conjecture in the general sense, each runner count on a loop shaped as a
convex regular polygon, using Kepler's second law creatively; keep tournaments and Hadwiger-Nelson in
view." The repo's LRC lane is the `k = 13` programme (LRC(14): THM-527/661/663, the relation-lattice
dictionary THM-1820, the runner/tournament dictionary THM-639) and polygon lenses already existed
(oracle-S517: the tight witness `{1..k}` at `t = 1/n` is the regular `n`-gon minus the observer with the
rotational tournament; codex-S529 outside/inside arcs; codex-S431 vertex-time half-turn tournaments).
Kepler's law had never been used. Literature state (CORE-PAPERS): `k <= 6` Barajas-Serra, `k <= 12`
by a computer-assisted sieve (preprint), `k = 13` open.

## 2. THM-4737, the Kepler-polygon reformulation

Loop = regular `n`-gon of perimeter 1, observer at a vertex, side = threshold `1/n`.

- **(K1)** Constant-speed motion along the sides is Kepler motion about the centre (area rate `a v/2`
  on every side), and about the observer's vertex the area rate on side `s` is `d_s v/2`,
  `d_s = R(cos(pi/n) - cos((2s+1)pi/n))`, zero exactly on the two incident sides; per lap the swept
  area about either point is the polygon area.
- **(K2) Kepler criterion.** The observer is lonely at `t` iff every runner has nonzero angular
  momentum about the observer, i.e. every observer area clock is strictly increasing; each clock is
  flat exactly `2/n` of the time; the lonely set has measure `m(D)`, positive iff `kappa(D) > 1/n`,
  finite (measure zero) exactly on tight sets. 20,000 random checks, no mismatch.
- **(K3) Vertex times.** At `t = j/N` the observer is lonely iff all `v_i j mod N` lie in
  `[N/n, N - N/n]`; for `N = n` this is "`n` divides no speed", so the `n`-gon alone proves the
  conjecture for every speed set without a multiple of `n` (exhaustively confirmed, `k = 2..6`).
  The **polygonal order** `N_min(D)` (least `N` with a lonely vertex time) is a new invariant: `n` on
  all six tight sets in range (`{1,2}`, `{1,2,3}`, `{1,2,3,4}`, `{1,3,4,7}`, `{1,2,3,4,5}`,
  `{1,3,4,5,9}`), `2` when all speeds are odd, and unbounded but slow: records on `(1, m)` are
  `3, 5, 8, 11, 15, 16, 17, 19` at `m = 2, 6, 12, 18, 72, 3960, 8640, 16560`; maxima `15, 18, 17`
  over all primitive sets with speeds `<= 200, 45, 24` for `k = 2, 3, 4`.
- **(K4) Real speeds.** The Kepler flow is Haar-equidistributed on the annihilator torus of the
  relation lattice; the conjecture for `D` is equivalent to the zonotope `M[1/n, 1-1/n]^k`
  containing an integer point; relation rank `<= 1` is settled by the range of a linear form
  (`sum|m_i| >= 3`); integer speeds are the maximal-rank case, the Kepler curve itself.

Exact `kappa(D)` over all primitive `k`-sets of speeds `<= M` (`k = 2..6`): no violation of
`kappa >= 1/(k+1)`; the tight sets in range are exactly the six listed.

## 3. HYP-9350

(a) the growth law of the polygonal order (a provisional bound `12` for two runners was refuted by the
session's own census at speeds `<= 200`; conjecture `N_min = O(log max D)`, `(1, m)` extremal, tight
sets have `N_min = n`); (b) the Kepler measure `m(D)` and its minimum over non-tight sets; (c) the
colouring bridge: `chi_f(G(Z,D)) = 1/kappa(D)` (Chang-Liu-Zhu 1999, CITED from memory, to verify) makes
the conjecture a fractional-colouring statement about one-dimensional distance graphs, and the vertex
time criterion is a multiplier statement about the chord circulant `C_N(D)`; whether the HN lane's
residue/idoneal colourings (THM-4558/4566) can certify lonely vertex times with controlled `N` is
untested; (d) what a polygon-based general proof must supply: the arithmetic of multipliers and moduli,
since the tight families have finite lonely sets and no measure argument can close them.

## 4. Honest scope

Nothing here proves the conjecture for any new `k`. The contributions are the Kepler dictionary (K1-K2),
the one-line `n`-gon settlement of all speed sets without a multiple of `n`, the polygonal order with
its growth data, the zonotope form for real speeds with the rank `<= 1` lemma, and the colouring
bridge as a research direction. The arXiv API was unreachable from this machine, so the Chang-Liu-Zhu
identity is recalled, not verified.

## 5. Files

`kepler_polygon_lrc.py` / `.out` (rates, criterion, exact kappa, lonely measure, N_min, tight family),
`polygon_orders.py` / `.out` (n-gon law, tight sets, extremal N_min, polygon-size hierarchy),
`polygonal_order_bound.py` / `.out` (maxima and histograms over larger ranges),
`pair_one_m_growth.py` / `.out` (records on the pairs (1, m), m <= 30000).
