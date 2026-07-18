---
id: THM-980
title: THE LIVE-FLOOR BRIDGE AND THE REDUCTION CAPSTONE — the program's nucleus served, and the residual stated in its final form. (I) THE WEIGHT-1 SAMPLING BRIDGE (three lines, the modulus-supply's simpler sibling): the good indicator has ≤ E(V) = 2·Σv breakpoints, so liveCount(q) ≥ q·(μ₀ − E(V)/q) — hence ANY family with positive good measure fires the census live floor at the EXPLICIT modulus q₀ = ⌈2E/μ₀⌉, with liveCount ≥ q₀μ₀/2 > 0. REFEREED on the two structured controls where the B₅-route fails (geometric: B₅ < 0 yet 163,790 live points at q₀ = 589,934; deep well: 1,040 at 43,520) and on the tight control (μ₀ = 0: zero live at every q — no false floor); (II) THE CERTIFICATE PAIR (the pipeline's two explicit moduli): DISSOCIATED branch → THM-979 (B₅-route: q₀ ≈ 51,900·Σv, unconditional via the T_s ladder); ANY μ₀ > 0 family → this bridge (live-route: q₀ = 2E/μ₀, conditional only on a μ₀-floor); both feed lonely_of_census with every hypothesis computable; (III) THE REDUCTION CAPSTONE — the residual in final form: with the assembly's closed branches (non-covering sieve, window ≤ 22, comparable, 91-dominant, detuned, coarse, common-residue), the cascade/gluing regimes (ratio ≥ 7/3 universal; block gluing), the dissociated branch (T_s ladder + THM-979), and this bridge, **LRC(14) reduces to: μ₀ > 0 for every trapped-core family** (covering, spread > 13, compressed, distinct, max ≥ 23, non-clusterable, chain-dense, NON-dissociated — i.e., carrying a small relation — and not census-certified at the supplied moduli) — the classical inf-L positivity on the thinnest family the program has ever stated, WITH the conversion guarantee: any μ₀-floor, however obtained (fee ledger, witness construction, covering rigidity), becomes a machine certificate through an explicit modulus. The nucleus is no longer entangled with certificate mechanics — it is pure measure positivity on a precisely-cut family
status: (I) proven (per-breakpoint counting, weight ≤ 1) + refereed three ways incl. the no-false-floor control; (II) assembly of proven pieces; (III) the reduction statement with every cited branch's file named — the honest final form of what remains
source: kind-pasteur-2026-07-17-S128 (cont.47; owner: finish the math, then all the formalization)
depends_on:
  - THM-979 (the B₅-route modulus), THM-950/census pipeline (the consumer)
  - THM-928/932/933 (cascade/gluing), the assembly branches (LRC14GrandAssembly)
related:
  - the fee ledger (death-star THM-938..945) and witness machinery (THM-936) — the two live routes to the remaining μ₀-floors
  - OPEN-Q inf-L (the classical nucleus, now in its final cut)
script: 04-computation/live_floor_bridge_referee_kps_S128c47.py -> 05-knowledge/results/live_floor_bridge_referee_kps_S128c47.out
---

# THM-980 — the live-floor bridge; the reduction capstone

## (I) The bridge

The good set is a finite union of intervals with ≤ E(V) endpoints; a modulus-q grid
misses at most E(V) of its mass-q points: liveCount ≥ qμ₀ − E. Three lines, weight 1.

## (II) The pair

| branch | route | modulus | conditional on |
|---|---|---|---|
| dissociated | B₅ (THM-979) | ≈ 51,900·Σv | nothing (T_s ladder) |
| μ₀ > 0 | live (this) | 2E/μ₀ | a μ₀-floor |

## (III) What remains — final form

μ₀ > 0 on trapped cores (the precise conjunction in the title). Every branch that could
be closed by certificate mechanics has been; the remaining question is pure measure
positivity, and any answer converts to a kernel certificate through an explicit modulus.

## Named next (the formalization phase, per the owner's sequencing)
- Lean: the weight-1 sampling lemma (interval point-counting — simpler than the B₅
  version; the q₀ arithmetic is norm_num) → wire to lonely_of_census.
- The μ₀-floors themselves: the fee ledger's trapped-core supply (death-star) and the
  witness constructions (opus THM-936 style) — the two live attack routes.
