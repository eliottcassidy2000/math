---
id: THM-686
title: Measure floors for the residual families — (A) THE WINDOW-22 MEASURE-FLOOR THEOREM (exhaustive, exact): every covering 13-set in [1..22] has μ(S) ≥ 7/858, attained uniquely in census order at {1..14}∖{6}, with ZERO tight families (μ > 0 everywhere on all 31471); with THM-685 the whole branch is certified at every modulus by [banks q ≤ 12135] + [transfer] — the quantified-margin complement of LEM-024; (B) THE INSTRUMENT REFUTATION: truncated Bonferroni (continuum B5, chain-B3/B5) NEVER fires on the raw small-speed covering class — S₅ explodes to 156–172× iid under relation stacking (the standing law's 12th confirmation, now in the continuum; THM-671's discrete B5 domain does not port); (C) the t=2 deviation lemma |A∞(a,b) − 36/49| ≤ (24/7)/max(a,b) proved, sharp constant measured 6/49 at (1,13)
status: PROVED ((A) is a finite exact computation, exhaustively executed — 31471 families, integer-scaled rational sweeps, no floats anywhere; (B) exact ledgers; (C) Koksma proof below, constant verified over all coprime pairs to b = 40). Cross-validations: the sweep reproduces boxeph's chain ladder μ₂..μ₅ exactly (11/14, 5/7, 9/14, 33/56); the class floor 7/858 EQUALS mac-mini's k=1 engine floor in LRCWitnessFloorRepair (flagged for mac-mini to confirm the objects are one); the ≤18 sub-class has exactly 966 members = the kps bank domain.
source: klein-2026-07-10-S235 (HYP-5900; executing the S234 successor lead "measure floors for the residual families")
depends_on:
  - THM-685   # the transfer that converts these floors into certificates at every modulus
related:
  - LEM-024 / opus LRCWindow22Census (the loneliness half; this adds the quantified margin)
  - boxeph HYP-5853 (μ_L ladder = A∞ on chain directions — cross-checked exactly)
  - mac-mini LRCWitnessFloorRepair (engine floor 7/858 — the same number, flagged)
  - THM-671 (discrete B5 — its per-ruler cluster domain does NOT port to raw families; scoping)
---

# THM-686 — measure floors for the residual families

## (A) The window-22 measure-floor theorem (exhaustive, exact)

Enumerate every 13-subset of [1..22] that is covering (∀q ∈ 2..14 ∃v: q | v;
note 13 and 14 are forced elements, so primitivity is automatic): there are
**31471** (14002 with min = 1 — the census branch — and 17469 with min ≥ 2,
matching the opus/kps counts). For each, μ(S) computed EXACTLY by the
integer-scaled breakpoint sweep (THM-685 machinery; scale 28·lcm(S), no
floating point anywhere). Results:

> **μ(S) ≥ 7/858 ≈ 0.008159 for every window-22 covering family**, with the
> minimum at {1,2,3,4,5,7,8,9,10,11,12,13,14} = {1..14}∖{6} (the near-AP);
> **zero families have μ = 0** — no covering family in the window is tight;
> median μ = 0.0676, max = 0.1564.

With THM-685 (|LM(q) − qμ| ≤ Σv): every family has live multipliers at every
q > Σv/μ(S); the worst threshold over the branch is **max q\* = 12135**.
Hence the window-22 branch is certified at EVERY modulus by [finite bank
q ≤ 12135] + [transfer beyond] — the quantified-margin complement to
LEM-024's 6-witness loneliness (which gives one witness per family; this
gives witnesses at all large moduli with clearance ≥ 1/14, uniformly).

The ≤18 sub-class (exactly the 966 bank domain) has the same floor and the
same argmin; on it max q\* ≈ Σv/μ evaluated per family (worst ≈ 12135 comes
from the [1..22] tail).

## (B) The instrument refutation — truncated Bonferroni does not port

The continuum quintic Bonferroni μ(S) ≥ B5(S) = 1 − S₁ + S₂ − S₃ + S₄ − S₅
(danger sets D_l, S_k = Σ_{|U|=k} meas(∩_U D), all S_k exact via
inclusion–exclusion of line sweeps) is a valid unconditional lower bound, and
its iid reference value is +0.1221 (THM-671's number). MEASURED, exact:

- B5 = **−6.87** (deep well), **−7.89** (worst-covering {1,2,3,4,7,…,17}),
  **−0.98** (GEN), **−6.16** (DIL); on a 49-instance spread of the ≤18 class:
  **0/49 fire**, values −6.7 … −8.9.
- Cause, quantified: relation stacking explodes the deep intersection terms —
  S₅ = 11.9 (well) and 13.2 (worst-covering) against iid 0.0766: **156–172×**.
  S₂ elevates only 1.3–1.4×; the blow-up is in the tail exactly where the
  truncation cuts.
- Chain-coarsening (boxeph super-runners, exact chain-danger intersections)
  comes remarkably close on the well — chain-B5 = **−0.00086** — but still
  does not fire; on 7 chains the full depth-7 chain I-E is exact (= μ) anyway.

Scoping vs THM-671: the discrete B5's 100% covering fire rate lives on the
per-ruler CLUSTER objects (P∪L co-offset instances at moduli q ∈ (V, 2V]),
not on raw small-speed covering families; the continuum limit of the raw
families is where the coherence is fully expressed. This is the standing
law's 12th documented confirmation — truncated/absolute instruments fail on
resonance-stacked objects while the exact signed object (μ itself, one
O(Σv log Σv) sweep) is computable and fine. **For finite residual classes the
instrument IS the exact census; truncated Bonferroni has no role there.**

## (C) The t=2 deviation lemma (the a-priori uniformity input)

For coprime 1 ≤ a < b: unfolding the larger coordinate, A∞(a,b) =
(1/b)Σ_{j<b} g(j/b) with g(x) = meas{β ∈ B̃ : frac((a/b)β + x) ∈ B̃},
∫g = (6/7)² (Fubini), and g of bounded variation V(g) ≤ 24/7 (the window
sweeps ≤ 2 band translates, each a trapezoid of height ≤ 6/7). Koksma with
the exact grid (D* = 1/b):

> **|A∞(a,b) − 36/49| ≤ (24/7)/max(a,b).**

Measured sharp constant over all coprime pairs b ≤ 40: max |dev|·b = **6/49**
at (a,b) = (1,13) — i.e. the true constant is (6/7)(1/7), the band's
variance; the proved 24/7 is 28× conservative. [t ≥ 3 uniformity: the S234
exact tables show non-relation triples at the same 1/height scale; the
class-uniform t ≥ 3 lemma remains open but is NOT needed for (A) — the
census is exact.]

## Consumption and the honest state of the measure-floor program

- Bounded residual classes (window-22 today; any future bounded window the
  assembly needs): the exact census is the floor theorem — finite, decidable,
  Lean-shaped (rational interval sweeps = boxeph's μ_L evaluator stack;
  compression: most families clear μ ≥ 0.03 by wide margin, only the near-AP
  tail needs individual data — a 6-witness-style shrinkage is available).
- Unbounded residual content: unchanged by this session — it lives in the
  other assembly legs (hpartA / the analytic node), not in the window branch.
- The deep well {1..12,182} is outside every bounded window; its exact
  μ = 4637/194040 (THM-685) with q\* = 10880 certifies IT individually the
  same way.

## Verification & files

`04-computation/lrc14_measure_floors_klein_S235.py` (+ `.out`): the sweep vs
μ₂..μ₅ cross-check, exact B5/chain-B5 ledgers, the ≤18 census, the deviation
lemma scan. `05-knowledge/results/lrc14_window22_mu_census_klein_S235.out`:
the full 31471-family census (floor, argmin, zero-tight, max q\*).
