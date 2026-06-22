---
id: HYP-2895
title: The LRC(14) covering crux decomposes -- bounded (compactness) + unbounded (equidistribution); large committed speeds equidistribute, keeping loneliness at unbounded-denominator witnesses (explains THM-566)
status: STRATEGY + verified on the adversarial family. Single-large-speed equidistribution is clean (explains THM-566); multi-large-speed reduces to mac-mini's resonance/atom analysis (HYP-+2878)
source: kind-pasteur-2026-06-22-S31o
related:
  - THM-523    # q-witness lemma: LRC(14) reduces to COVERING sets (contain a mult of every q<=14)
  - THM-566    # adversarial: no BOUNDED-denominator witness for {1..11,13,84*lcm(1..B)}
  - THM-527    # compactness reduction (bounded speeds)
  - HYP-+2878  # mac-mini covering-system / prime-basis / single-atom over-determination
  - THM-560    # structured/sporadic tiler characterization (the boundary layer)
---

# HYP-2895: the covering crux decomposes via equidistribution

## The reframe (THM-523)
LRC(14) reduces to COVERING sets: any speed set omitting all multiples of some `q<=14` is lonely at
`t=1/q` (q-witness lemma). So a counterexample must contain a multiple of EVERY `q in {2,...,14}`. The
tight locus (AP, GW; THM-560) is the trivial boundary -- both omit multiples of 14, lonely at `t=1/14`.

## The decomposition (kps-S31o, `lrc_equidistribution_covering_kps.py`, `lrc_small_covering_kps.py`)
A covering set `S` falls in one of two regimes:

**(A) BOUNDED covering sets** (all speeds `<= B`): lonely with MARGIN. E.g. `{2,...,14}` (covers all
`q<=14`) has `M = 1/8 = 0.125 >> 1/14`, witnessed at `t = 1/16`. Several small covering sets all give
`M = 1/8`. Finitely many shapes per bound -> mac-mini's compactness reduction (THM-527) closes these.

**(B) UNBOUNDED covering sets** (one speed `v >> ` the rest): the large speed EQUIDISTRIBUTES (Weyl).
`U_v = {||vt||<1/14}` has measure `1/7`, and for large `v`, `meas(U_v ∩ L) -> meas(L)/7` for the
1/14-lonely set `L` of `S\{v}`. So `U_v` covers only `~1/7` of `L`; loneliness SURVIVES on `~6/7` of
`L`, at an unbounded-denominator witness. **VERIFIED on THM-566's adversarial family**:
`{1..11,13,840}` -> `M = 14/169 = 0.0828`, witness `t = 493/845`; `{1..11,13,5040}` -> `M = 84/1009 =
0.0833`, witness `t = 2102/5045`. Both COVERING (contain a multiple of every `q<=14`, `14|v` so
`t=1/14` fails), yet `M` stays near the 12-subset's `1/12 >> 1/14`.

## This EXPLAINS THM-566
THM-566 shows no FIXED denominator `D` witnesses all covering sets (the tail `84*lcm` sits at the
observer for every numerator over `Z/D`). HYP-2895 resolves the apparent paradox: the witness EXISTS,
it just has UNBOUNDED denominator (`845`, `5045`, growing with `v`). Equidistribution guarantees the
witness; it is the bounded-denominator search that fails, not loneliness. So THM-566 is not an
obstruction to LRC -- it is an obstruction to a NAIVE finite-denominator certificate, exactly resolved
by the equidistribution (unbounded-denominator) witness.

## What is clean vs open
- CLEAN: the SINGLE-large-speed case. One speed `v -> infinity` over a fixed lonely core `S\{v}`:
  `M(S) -> M(S\{v})` and the witness denominator grows. Explains the THM-566 family directly.
- OPEN (= mac-mini's territory): MULTIPLE large speeds. For integer speeds the joint orbit
  `(v_1 t,...,v_k t)` is a closed geodesic on `T^k` (not equidistributed on the full torus -- the `v_i`
  are rationally dependent), so the joint danger measure is governed by RESONANCES/gcds, not `(1/7)^k`.
  This is exactly mac-mini's strong-component-atom / CRT-independence analysis (HYP-+2878). So HYP-2895
  gives the CONCEPTUAL reason for over-determination (large speeds can't conspire to cover the lonely
  set) and reduces the multi-speed case to the resonance count.
- OPEN: the MODERATE regime (speeds comparable, neither small-margin nor equidistribution) + the
  recursive-tightness gap (THM-560): if the lonely core `S\{v}` is itself near-tight (`meas(L) -> 0`),
  the equidistribution slack shrinks. This is the genuine residual, shared with OPEN-Q-108.

## Net
The covering crux is not monolithic: BOUNDED (compactness, lonely with margin) + UNBOUNDED
(equidistribution, lonely at large denominator) cover the extremes cleanly; the hard core is the
moderate/resonant middle = mac-mini's atom over-determination. Equidistribution gives the proof its
"why," explains THM-566, and shows the finite-denominator obstruction is not a loneliness obstruction.

-> THM-523, THM-566, THM-527, HYP-+2878, THM-560, OPEN-Q-108.
