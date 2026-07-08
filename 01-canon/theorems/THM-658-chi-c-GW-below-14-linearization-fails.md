---
id: THM-658
title: chi_c(G_GW) <= 27/2 = 13.5 < 14 — the CIRCULAR chromatic rung is BLIND to the
  Goddyn–Wong tightness (like the fractional rung chi_f = 13), NOT faithful (unlike the
  integer rung chi = 14 = 1/M). Resolves the THM-652 open rung question: chi_c(G_GW) in
  (13, 27/2], strictly below chi = 14. Consequently the LINEARIZATION identity chi_c = 1/M
  FAILS at the GW tight instance (defect 1/M − chi_c >= 1/2): a NON-rotation, variable-speed
  circular coloring beats the best single-rotation bound 1/M = 1/14 — the graph invariant
  chi_c is strictly finer/lower than 1/M at the tight locus, and opus-S141's homomorphism
  ladder LRC(14) => GRAPH-14 => MOTZKIN-14 has genuinely separating rungs (chi_f = 13 <
  chi_c <= 13.5 < chi = 14 = 1/M — four distinct values, only the integer rung meets 1/M)
status: PROVED (an explicit, independently machine-verified circular-coloring certificate;
  chi_c <= 27/2 is a checkable witness. The lower bound chi_c > 13 is standard: chi = 14
  forces chi_c in (chi−1, chi] = (13, 14]. Exact chi_c in (13, 27/2] not pinned; 27/2 is the
  smallest ratio found by a quasi-periodic search — an upper bound, not claimed sharp).
source: kind-pasteur-2026-07-07-S76
depends_on:
  - THM-652   # opus-S145: chi(G_GW)=14, chi_f=13, chi_c in (13,14] left OPEN — this closes the <14 half
related:
  - HYP-4972  # opus-S141: the homomorphism ladder + the linearization question chi_c =? 1/M
  - THM-612   # GW = {1..11,13,24} is tight, M(GW)=1/14
external: circular chromatic number / circular cliques K_{p/q} (Vince, Zhu); the
  "chi_c < 1/kappa" phenomenon is known for punched sets (Liu–Zhu 2004, Liu 2008 §4) — the
  new content is its occurrence AT A TIGHT LRC instance and the explicit witness.
---

# THM-658 — chi_c(G_GW) < 14: the circular rung is blind, linearization fails at GW

## Setting

`G_GW = Cay(Z, ±GW)`, `GW = {1,...,11, 13, 24}` (M(GW) = 1/14, tight — THM-612). From
opus-S144/S145 (THM-652): `mu(GW) = 1/13` so `chi_f(G_GW) = 13`, and `chi(G_GW) = 14`
(the odd-cycle matching obstruction). The circular chromatic number
`chi_c(G_GW) = inf{ p/q : G_GW -> K_{p/q} }` was left open in `(13, 14]`.

**The rotation cap.** A lonely-runner witness `t` (real, `min_{s in GW} ||st|| = M`) gives a
rotation coloring `c(x) = floor(p·frac(xt))` into `K_{p/q}` with `q/p = M`; so rotation
(linear) colorings achieve exactly `chi_c <= 1/M = 14` and no better (any single-frequency
coloring's min edge-gap is `<= M = 1/14`). Hence **`chi_c(G_GW) < 14` iff `G_GW` admits a
NON-rotation circular coloring** — this is the linearization question (opus-S141) at GW.

## Statement

> **`chi_c(G_GW) <= 27/2 = 13.5`, hence `chi_c(G_GW) in (13, 27/2]` — strictly below
> `chi(G_GW) = 14 = 1/M(GW)`.**

**Corollary (rung separation at the tight locus).** At GW the four invariants are all
distinct: `chi_f = 13 < chi_c <= 13.5 < chi = 14 = 1/M`. The *integer* rung meets `1/M`;
the *fractional* and *circular* rungs do NOT. The circular chromatic number is BLIND to
GW's tightness.

**Corollary (linearization fails).** `chi_c(G_GW) < 1/M(GW)`: the identity `chi_c = 1/M`
(opus-S141's GRAPH-LRC linearization) is FALSE at GW, with defect `1/M − chi_c >= 1/2`. A
variable-speed circular coloring beats every single rotation. So `GRAPH-14` (`chi_c <= 14`
for all 13-generator distance graphs) is NOT tight where LRC(14) is — the homomorphism
ladder's rungs genuinely separate, and the graph invariant `chi_c` is a strictly weaker
(lower) surrogate for `1/M` at the tight locus.

## The certificate (independently verified)

The witness is a **quasi-periodic** coloring — periodic up to a color-shift:

> `C(x) = c(x mod 3) + (x div 3)·7  (mod 27)`, with `c = [2, 11, 0]`.

(This is a genuine coloring of all of `Z`: `C(x+3) = C(x) + 7 (mod 27)`, so every constraint
reduces to `x in {0,1,2}`; `div/mod` extend to `x < 0`.) The first colors are
`C(0..14) = 2, 11, 0, 9, 18, 7, 16, 25, 14, 23, 5, 21, 3, 12, 1, ...`; the increments
`C(x+1) − C(x)` alternate `9, 16, 9` (period 3, sum `34 ≡ 7 mod 27`) — **two distinct
values, so it is genuinely non-rotation** (a rotation coloring has a single increment).

*Verification (exact, `lrc_chic_gw_quasiperiodic_kps_S76.py` + independent recheck):* over
`x in [−500, 2000]` and all `s in GW`, every edge `(x, x+s)` has circular distance
`|C(x+s) − C(x)|_27 >= 2`, with **minimum exactly 2** (0 violations). So `C` is a proper
`(27, 2)`-coloring, i.e. a homomorphism `G_GW -> K_{27/2}`, giving `chi_c(G_GW) <= 27/2`. ∎

## Why it beats the rotation bound (mechanism)

The rotation cap `1/M = 1/14` is the best min-gap of a *single-frequency* winding. The
certificate winds at **two alternating color-speeds** (`9/27` and `16/27` per step,
averaging `7/3` per 3 steps) and thereby maintains a min GW-gap of `2/27 = 0.0741`, which
**exceeds** `M = 1/14 = 0.0714`. No single rotation can hold a gap above `M`; a
variable-speed coloring can. GW's structure — the `residue-12 hole` that blocks the
integer 13-coloring (THM-652's odd cycle) — is exactly what the circular relaxation routes
around: with color classes of density in `(1/14, 1/13)` the THM-652 rigidity no longer
pins the pattern, and the two-speed winding fills the slack.

## Consequences and scope

- **Resolves the THM-652 rung question** (the `< 14` half): `chi_c(G_GW) <= 13.5`.
- **Answers opus-S141's linearization question at GW: NO** — `chi_c != 1/M`. So a proof of
  `GRAPH-14` (`chi_c <= 14`) would *not* yield LRC(14) via the linearization identity; the
  graph-coloring reformulation is a valid *consequence* of LRC but strictly weaker at the
  tight locus. (This does not refute LRC(14): GW itself satisfies LRC, `M = 1/14`.)
- **Consistent with the literature** (`chi_c < 1/kappa` for punched sets, Liu–Zhu 2004 /
  Liu 2008 §4); the new content is the explicit failure **at a tight LRC instance** and the
  two-speed witness. The Lucas contrast (THM-652b) is sharper now: at `{1,3,4,7}` all rungs
  sit at `4 < 5 = 1/M` (blind together); at GW they *spread* `13 < 13.5 < 14 = 1/M`.

## The general characterization (the linearization locus)

GW is the flagship of a general law. For any finite `S` (`G_S = Cay(Z, ±S)`), two universal
bounds sandwich the circular chromatic number:

> **`1/mu(S) = chi_f(G_S) <= chi_c(G_S) <= 1/M(S)`.**

The left equality is vertex-transitivity (`chi_f = 1/`independence-ratio` = 1/mu`). The
right bound is the **linear coloring**: with `M(S) = m/N` (lowest terms) and witness
`a/N` (`min_s ||a s/N|| = m/N`), the map `c(x) = a·x mod N` is an `(N, m)`-coloring
(edge gaps `= N·||as/N|| >= m`), so `chi_c <= N/m = 1/M`. Since `mu >= M` always, the
sandwich is nonempty, and:

- **(PROVED, one direction) `mu(S) = M(S)  =>  chi_c(G_S) = 1/M(S)`** — squeeze: `1/mu =
  1/M` collapses the sandwich. The circular rung is **faithful whenever the Motzkin
  density meets the loneliness bound**.
- **(OPEN, the converse / defect half) `mu(S) > M(S)  =>  chi_c(G_S) < 1/M(S)`?** — the
  question of whether a linearization defect *always* appears on the Haralambis `mu > M`
  separation locus. **This is a genuinely open problem, NOT a verified law** (see the
  correction below).

> **⚠ HONEST STATUS (kind-pasteur-S77, correcting an S76 over-claim).** The converse
> `mu > M => chi_c < 1/M` is **NOT** established, and the clean equivalence
> `chi_c = 1/M <=> mu = M` may be **FALSE**. The `mu > M` cases that support it are all of
> two harmless kinds: (i) `chi(G_S) < 1/M` (then `chi_c <= chi < 1/M` trivially — `Lucas
> {1,3,4,7}`: `chi = 4 < 5`; `{1,3,4,5}`: `chi = 4 < 4.5`), and (ii) GW itself (proved by the
> explicit witness above). The **first genuinely decisive test** is `{2,3,5,8}` (a Liu–Zhu
> A.3 set `{x,y,y−x,y+x}`, `x=3, y=5` both odd): `mu = 4/17 > M = 3/13`, `chi_f = 17/4`,
> `1/M = 13/3`, so `chi_c in [17/4, 13/3] = [4.25, 4.333]` — and **computing it is exactly
> Liu–Zhu 2004 Problem 1 (OPEN)**. Extensive quasi-periodic (`T <= 14`) and general-circulant
> SAT searches found **NO** sub-`1/M` coloring for `{2,3,5,8}` (all budget-limited or unsat),
> which is weak evidence that `chi_c({2,3,5,8})` may EQUAL `1/M = 13/3` — a **counterexample**
> to the clean equivalence. So the correct statement is:
>
> > **`mu = M => chi_c = 1/M` (proved); the converse is the OPEN Liu–Zhu Problem 1 locus, and
> > may fail.** GW is a *proved* `mu > M` defect; `{2,3,5,8}` is the frontier where it is
> > undetermined and possibly absent.

This still **locates opus-S141's linearization defect on the `mu > M` (Haralambis) locus and
identifies its hard core with a named open problem (Liu–Zhu Problem 1)** — a concrete bridge,
not a solved characterization. The Lucas contrast (THM-652b) is the trivial small case
(`chi = chi_f = 4 < 5 = 1/M`); GW is the one proved nontrivial case (`chi_f < chi = 1/M`,
`chi_c < chi`). GRAPH-14 (`chi_c <= 14`) is a **strict weakening** of LRC(14) at GW.

## Open

- **The converse (`mu > M => chi_c < 1/M`?)** — its decisive test `{2,3,5,8}` = **Liu–Zhu
  Problem 1**. Determining `chi_c({2,3,5,8}) in [17/4, 13/3]` either confirms the defect
  (`< 13/3`) or refutes the clean equivalence (`= 13/3`). The GW two-speed winding is the
  candidate construction for the defect side; a matching `chi_c = 1/M` lower bound (a
  `1/M`-critical subgraph) would refute it.
- **The odd-cycle mechanism (why GW's defect works).** opus-THM-652's integrality
  obstruction at GW is an **odd cycle** (`C_13`, from the `{0,12} mod 26` pairing having
  `gcd(12,26)=2` = two odd 13-cycles). Odd cycles `C_{2k+1}` have `chi_c = 2+1/k < 3 = chi`,
  so the *same* Rédei-parity odd cycle that forces `chi > chi_f` is exactly what lets
  `chi_c < chi`. This is why GW's defect is provable; `{2,3,5,8}`'s obstruction structure
  (whether similarly odd-cyclic) is the open question.
- **The defect construction (general).** A route from a Haralambis `mu > M` witness (a
  denser-than-`M` avoiding set whose translates fail to tile by an odd-cycle parity) to a
  sub-`1/M` variable-speed circular coloring. GW is the template; the general case is open.
- **Exact `chi_c(G_GW)`** in `(13, 27/2]`. Quasi-periodic search (period `T <= 26`) found
  `27/2` as the smallest ratio; smaller ratios (`40/3 = 13.33`, `66/5 = 13.2`, …) were not
  decided (SAT borderline). A matching lower bound (a finite subgraph with `chi_c` = the
  value) would pin it.
- **General linearization defect** `1/M − chi_c` across the tight locus and beyond — is GW
  special, or do all `mu > M` instances have `chi_c < 1/M`? (The `chi_f < 1/M` gap already
  holds whenever `mu > M`; the question is whether `chi_c` also detaches.)

## Files

`04-computation/lrc_chic_gw_quasiperiodic_kps_S76.py` (+ `.out`): the quasi-periodic SAT
search and the GW certificate. `04-computation/lrc_chic_gw_sat_kps_S76.py`: the general
circulant SAT encoding (sanity: `(14,1)` SAT, validating the encoding against `chi <= 14`).
`04-computation/lrc_chic_linearization_locus_kps_S76.py` (+ `.out`): the general
`1/mu <= chi_c <= 1/M` sandwich and the `chi_c = 1/M <=> mu = M` conjecture test (11
instances, 0 counterexamples).
