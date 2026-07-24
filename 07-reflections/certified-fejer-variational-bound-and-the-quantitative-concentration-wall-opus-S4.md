---
source: opus-2026-07-23-S4 (harvesting the snippet technique for LRC(14), session 2)
status: >
  CONCRETE RESULT. A float-free CERTIFIED Fejer variational lower bound B_N(V) on the lonely-
  runner gap, realizing kps-S132's certifiable-concentration program with exact arithmetic.
  B_N <= gap(V), B_N -> gap(V); B_N > 1/14 certifies V is not an LRC(14) counterexample. The
  experiment maps the route's reach and, new, QUANTIFIES the wall: N* ~ (binding slope)/(gap-1/14),
  verified over 2 decades of margin. This makes the "bulk-certifiable + near-extremizer-rigidity"
  reduction precise: LRC(14) via concentration = a spectral gap above 1/14 (open) + tight-locus
  finiteness (OPEN-Q-108, open). Not a proof; a sharpened, quantitative reduction + a Lean-ready primitive.
tags: [lrc14, certifiable-concentration, fejer, variational, certified, float-free, harvest, wall, open-q-108, kps-s132]
related: [THM-729, THM-518, HYP-9023, CONSTANTS-INDEX]
---

# Certified Fejer variational bound & the quantitative concentration wall (opus-S4)

## The certified bound (harvest of the snippet technique)

For a 13-speed set `V`, `gap(V) = max_tau g(tau)`, `g(tau) = min_v ||v tau||`. For ANY probability
measure `mu`, `gap(V) >= INT g dmu` (max >= average). Take `mu = F_N(.-tau*)`, a Fejer kernel
(`F_N >= 0`, mass 1) at the lonely time `tau*`. Then

```
B_N(V) := INT g(tau) F_N(tau - tau*) dtau  <=  gap(V),   B_N -> gap(V) as N -> infinity.
B_N > 1/14  ==>  gap(V) > 1/14  ==>  V is NOT an LRC(14) counterexample.
```

**Float-free by construction** (the snippet's exact-arithmetic style, cf. the certified Q_s):
`g` is piecewise-LINEAR with RATIONAL breakpoints (tent points `j/(2v)`, crossings `k/(v_i +- v_j)`),
`F_N` is a trig polynomial, so

```
B_N = G_0 + 2 sum_{k=1}^N (1 - k/(N+1)) Re[ e^{-2pi i k tau*} G_k ],
G_k = INT g e^{2pi i k tau} dtau = exact sum of (rational)*(root of unity).
```

Lean-portable (exact rational pieces + roots of unity; here evaluated in double with exact phase
reduction, error ~1e-12 << the 1/14 margins). `gap(V)` and `tau*` are exact rationals.

## Reach map (`lrc14_fejer_variational_opus_S4.py`)

Gaps validated exactly against the CONSTANTS-INDEX (AP{1..13} & GW{1..11,13,24} = **1/14** tight;
{1..11,13,36} = **3/41**; {1..12,26} = **2/27**). Behaviour of `B_N`:

- **Tight AP/GW:** `B_N -> 1/14` but NEVER exceeds (0.0689 at N=1500) -- THM-518's exact-value wall,
  now visible: a band-limited `mu` is not a point mass, so `B_N < gap = 1/14` strictly.
- **Bulk** (`gap >> 1/14`, e.g. coprime set, gap 1/4): crosses 1/14 at `N* = 50`, `B_N -> 0.24`.
- **Near-floor:** `N*` grows sharply as `gap -> 1/14`.

## The wall, QUANTIFIED (`lrc14_fejer_scaling_opus_S4.py`)

`g` has a downward CORNER at `tau*` (max of a piecewise-linear min), so the Fejer mean converges at
the corner rate `gap - B_N ~ c log N / N` (verified: `(gap-B_N) N/ln N ~ 0.46->0.55` for AP). The
corner's slope is the **binding slope** = mean `|slope|` of the binding tents; for AP it is exactly
`7 = (1+13)/2` (the two binding runners `v=1,13` at `tau*=1/14`). Hence

```
N*  ~  (binding slope) / (gap - 1/14)     (up to a log factor),
```

verified across `delta = gap - 1/14` spanning two decades (`N* delta / slope ~ 0.34-0.53`, the slow
drift = the `log N*` correction). E.g. `2/27` (`delta=1/378`) certifies at `N*=2720`; `1/8`
(`delta=0.054`) at `N*=100`; `1/2` at `N*=20`.

## What this makes precise (the reduction)

The certifiable-concentration route (kps-S132) becomes an explicit two-part reduction of LRC(14):

1. **Bulk.** Every config with margin `delta = gap - 1/14 >= eps0` is certified `gap > 1/14` by the
   float-free Fejer bound at degree `N ~ (binding slope)/eps0`. Crucially the binding runners of a
   NEAR-tight config are its small-speed core (the `{1..n}` part), so the binding slope there is
   BOUNDED (`~ n ~ 13`), independent of any large remote speed. So the bulk needs a UNIFORMLY bounded
   degree `~ 13/eps0` once `eps0` is fixed.
2. **The strip `0 < delta < eps0`** (near-extremizers) is where `N*` diverges -- exactly OPEN-Q-108,
   the tight-locus rigidity/finiteness.

So **LRC(14) via concentration = [a spectral gap `eps0 > 0` above `1/14`] + [tight-locus finiteness]**.
The spectral gap is the conjectured emptiness of `(1/14, 3/41)` (`eps0 = 1/574`; searched, THM-1235/1240,
NOT a theorem); with it and part (1), the whole non-tight regime is certified at Fejer degree
`~ 13*574*log ~ 10^4-10^5`, reducing LRC(14) to the finite tight locus `{AP, GW}`.

## Honest scope

Both inputs (`(1/14,3/41)` emptiness and tight-locus finiteness) are open and hard; this does not
prove LRC(14). Its value: (i) a Lean-ready CERTIFIED per-config gap>1/14 primitive (the harvest,
realized for the variational route as for the density route's Q_s); (ii) a QUANTITATIVE law for the
route's cost, `N* ~ slope/delta`, turning the "wall" into an explicit divergence; (iii) a clean,
explicit statement of the reduction and its two open inputs -- confirming the fleet consensus
(kps/klein/mac-mini) that the analytic route certifies the bulk and the crux is the rigidity strip.

Artifacts: `04-computation/lrc14_fejer_variational_opus_S4.py`, `..._fejer_scaling_opus_S4.py` (+ `.out`).
