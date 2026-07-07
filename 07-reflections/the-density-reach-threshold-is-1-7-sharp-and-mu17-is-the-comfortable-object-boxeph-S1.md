---
source: boxeph-2026-07-07-S1 (HYP-4760)
status: validity audit + strategic clarification (not a proof). Resolves the
  1/7-vs-2/7 threshold ambiguity between THM-527 and the recent burst; confirms
  (exactly) that E[maxgap] is NOT AP-minimized; reconciles the conflicting
  origin-gap numbers; and argues the reverse-Markov E[maxgap] detour trades a
  COMFORTABLE tail (mu_{1/7}) for a RAZOR-THIN mean.
tags:
  - lonely-runner
  - LRC14
  - density-floor
  - max-gap
  - three-gap
  - validity-audit
---

# The density->reach threshold is 1/7 (sharp), and mu_{1/7} is the comfortable object

Owner (remote): come to a deep understanding of the LRC(14) work, how the proof
state has changed under corrections, assess validity, find what's been missed.
This reflection audits the **density-floor -> reach** link and the recent
**reverse-Markov / E[maxgap]** burst, cross-checking every number exactly.

## 1. The sharp threshold is 1/7, not 2/7 (THM-527's 2/7 was the robust cousin)

THM-527.A states the Vmax-ruler criterion with circular **max-gap > 2/7** on the
co-offset config `{frac(e_i x)}` (`e_i = Vmax - u_i`, so `0 in E`); the recent
burst (opus-S130+) uses **max-gap > 1/7** on the speed config `{frac(e x)}`. These
looked like different objects. They are the **same object at two clearance
levels**, and the correct (sharp) one is **1/7**:

- A ruler-period at slow `x` is GOOD iff `exists phi in (1/14,13/14)` with
  `||phi - e x|| >= 1/14` for all `e`. Each config point forbids a `1/7`-wide
  `phi`-band (two `1/14` half-bands); a witness `phi` **exists** iff some gap
  `> 1/7`. Demanding the *free `phi`-arc* itself exceed `1/7` (a robust margin)
  is the `> 2/7` condition. So **`1/7` = "a witness exists", `2/7` = "a witness
  with `1/7` slack"**.
- **Verified exactly (`lrc_maxgap_threshold_mechanism_audit`):** the ACTUAL
  good-period fraction of finite instances `S = {Vmax} u {Vmax-e}` converges to
  the **1/7**-density, not the 2/7-density. For the binding `APcoff {0..12}`:
  predicted `mu_{1/7} = 477/1078 = 0.4425`; measured `rho_K = 0.46 (Vmax=200),
  0.44 (Vmax>=300)`. The 2/7-density (`0.1794`) is a strict under-count.
- **Speed-config = co-offset-config exactly.** `{frac(e x): e in E}` and
  `{frac((V-e)x)}` are related by `x -> -x` (a measure-preserving reflection), so
  `mu_theta(speed) = mu_theta(co-offset)` for every family and every `theta`
  (checked exact k=8..13). The recent `mu_{1/7}(speed set)` IS THM-527's object.

**Consequence:** the recent burst is aimed at the *correct sharp object*; THM-527's
`2/7` was conservative by exactly one `1/7`-clearance. The load-bearing quantity is
the **1/7 good-density `G2(E) = mu_{1/7}(E)`**.

## 2. E[maxgap] is NOT AP-minimized (exact) -- correcting klein-S153

The reverse-Markov reduction (kps-S57) `mu_{1/7} >= (7/6)(E[maxgap]-1/7)` sends the
floor to a **mean** `E[maxgap] > 1/7`. kps-S58 and klein-S153 named the crux
"**AP-minimality of E[maxgap]**"; klein-S153 reported "AP is the minimizer, 48%
margin (descent converges to AP)". This is **false**, exactly:

> `E[maxgap](AP_13) = 93/440 = 0.211364`, but
> `E[maxgap](GW={1..11,13,24}) = 140413631/669278610 = 0.209798 < 93/440.`

GW is the second M-tight family (`AP[12->24]`); the one-swap minimizer of
`E[maxgap]` is exactly this `12->24` swap, sitting `0.0016` below the AP -- inside
klein's grid precision, which is why a random-restart local descent "converged to
AP" and missed it. (Independently confirmed same day by death-star-S1 and
opus-S133; opus found lower shapes `~0.205`, monad-explorer HYP-4787 found
`2*{1..11}u{11,13} = 0.19699`.) So: **the AP minimizes the max-gap TAIL
(`mu_{1/7}`, uniquely) but NOT the max-gap MEAN (`E[maxgap]`).** The two functionals
have different extremizers -- kps-S58's "shared AP-minimality core" does not hold.

The reduction itself survives (all these minima are `> 1/7`), but the clean
"prove `E[maxgap] >= E[maxgap](AP)`" strategy is unavailable.

## 3. The origin-gap numbers reconciled (klein 0.156 vs opus 0.134 vs kps 0.137)

`E[gap@0] = 2 E_x[min_i frac(e_i x)]` (speed config, `x->-x` symmetry). Exactly:
- `E[gap@0](AP_13) = 93/440 = E[maxgap](AP_13)` **exactly** -- so klein's "origin
  saturates the AP max gap (gap@0 = maxgap a.e.)" is CORRECT for the AP, and
  kps-S58's `0.137` was NOT the AP value.
- `inf_E E[gap@0]` is a **knife-edge at 1/7**: inhomogeneous APs `{a+dj}` bottom out
  at `0.1619` (`{17,23,..,89}`); opus's `{6,11,..}` gives `0.170` (not `0.134`);
  irregular spread minimizers reach `~0.147`. So the true single-anchor inf is
  `~0.147`, just above `1/7` and too thin to be a robust route. klein's `0.156` was
  optimistic, opus's `0.134` too low.
- The **2-anchor** floor `E[max(gap_0, gap_{1/2})]` has `inf ~ 0.187` (comfortable);
  `{j/8}` recovers `~0.210`. So klein's anchor-floor *idea* is sound, but needs
  `>= 2` anchors -- the single anchor is a knife-edge.

## 4. The finite-Vmax bridge is benign; the tail route is comfortable, the mean route is razor-thin

The density->reach reduction is exact only as `Vmax -> inf` (`rho_K -> mu_{1/7}`);
a fixed instance needs `rho_K > 0`. Because `mu_{1/7}` is dilation-invariant, any
compact binding family dilates to large `Vmax` with `mu_{1/7}`, `M` unchanged.
**Probe (`lrc_finite_vmax_bridge_probe`):** `V_0 <= 14` for every **bounded-spread**
binding shape (a good period already at the minimal ruler), and `rho_K` tracks the
floor by `Vmax ~ 100`. So for the bounded-spread regime -- exactly the extremal
regime THM-527.D reduces to -- the `O(#arcs/Vmax)` finite error is dominated by the
floor and `Vmax < V_0` is a tiny finite check.

**Honest caveat (credit monad-explorer HYP-4787).** This benignity is for
*bounded-spread* co-offsets. The genuinely open Part-A sub-case is the
**compressed all-big / spread~Vmax** family, where the co-offsets `e_i = Vmax-u_i`
are themselves `~Vmax` (fast, not frozen), so the slow-fast separation the ruler
rests on degrades and no `o(Vmax)` arc bound is yet written. So the Part-A gap is
benign where the extremizer lives, but not universally -- the correct factoring is
monad-explorer's `[G2 >= m_P uniform] + [Vmax >= V0(m_P)] + [finite check below V0]`,
with the arc bound outstanding for spread~Vmax shapes.

Now the decisive contrast, using monad-explorer HYP-4787's honest bar
`T* = 1/7 + (6/7) m_P ~ 0.19127` (the skeleton consumes `G2 >= m_P ~ 0.0565`,
not `G2 > 0`):

| route | object | needs | binding value | margin |
|---|---|---|---|---|
| **TAIL** | `mu_{1/7}(E) = G2` | `>= m_P ~ 0.0565` | AP `= 0.4425`, inf ledger `~0.32` | **~0.28 (huge)** |
| **MEAN** | `E[maxgap]` | `> T* ~ 0.19127` | min `~0.197` (and dropping) | **~0.006 (razor)** |

The reverse-Markov step is *lossy* (Markov's `7/6` throws away the tail shape), and
it trades the comfortable tail floor (`0.44` vs `0.057`) for a mean that must clear
`0.191` when its minimizer is already at `0.197` and every new adversary lowers it.

## 5. Strategic reading (what's been missed / where the aim drifted)

- The **load-bearing, comfortable object is `mu_{1/7}` (the sharp-1/7
  good-density)** -- AP-minimized (uniquely), dilation-invariant, floor `~0.44`,
  and it is exactly the finite good-period fraction. Its remaining Part-A gap
  (the `o(Vmax)` arc bound) is benign.
- The **reverse-Markov `E[maxgap]` detour, while elegant and "one order statistic",
  moved off the AP-extremal structure and onto a razor-thin mean**. It is the right
  tool only for a crude lower bound, not the primary crux. The fleet's energy on
  "`inf_E E[maxgap] > 1/7`" is chasing a threshold that is really `0.191`, with
  `~0.006` of room -- a moving target as adversaries improve.
- **Recommendation (independently reached by monad-explorer HYP-4787, same day):**
  re-anchor the crux on `mu_{1/7}(E) >= m_P` (the per-k tail floor (A')) for
  single-scale primitive families -- either via AP-minimality (opus's (A')/S130
  exact constants, the clean rigidity) or a *crude* uniform lower bound (only
  `0.057` is needed, so a loose three-gap / Riesz-product estimate suffices, no
  sharp extremal). This is strictly easier than a sharp mean bound. My distinctive
  additions to this shared conclusion: the **1/7-sharp/2/7-robust identification**
  (ties THM-527's ruler to the burst's object, and shows the recent aim is
  correct), the exact **`E[gap@0](AP)=93/440=E[maxgap](AP)`** identity, the
  **single-anchor knife-edge** reconciliation of the klein/opus/kps numbers, and the
  **`V_0<=14` bounded-spread benignity** of Part A.

## Ledger

- **New this session:** the 1/7-sharp-vs-2/7-robust identification (THM-527 <->
  burst), verified against actual finite good-period fractions; speed=co-offset
  reflection identity; exact `E[maxgap](GW) < E[maxgap](AP)`; exact
  `E[gap@0](AP)=93/440=E[maxgap](AP)`; the single-anchor knife-edge (`~0.147`) vs
  2-anchor comfort (`0.187`); `V_0 <= 14` finite-Vmax benignity; the
  tail-vs-mean margin table.
- **Corrects:** klein-S153 ("AP minimizes E[maxgap], 48% margin"; "inf E[gap@0]
  ~0.156"). **Confirms:** death-star-S1, opus-S133, monad-explorer HYP-4787
  (E[maxgap] not AP-min; T* bar). **Credits:** kps-S57 (reverse-Markov),
  opus-S130 (exact mu_{1/7} constants), THM-527 (the ruler reduction), mac-mini-S15
  (three-gap frame).
- **Not claimed:** any proof of `mu_{1/7} >= m_P` or LRC(14). Stays on the sup side
  (survives MISTAKE-117).
- **Files:** `04-computation/lrc_maxgap_threshold_mechanism_audit_boxeph_S1.py`,
  `lrc_Emaxgap_exact_apmin_boxeph_S1.py`, `lrc_gap0_anchorfloor_exact_boxeph_S1.py`,
  `lrc_gap0_inhomog_ap_exact_boxeph_S1.py`, `lrc_finite_vmax_bridge_probe_boxeph_S1.py`
  (+ `.out` in `05-knowledge/results/`).
- **Next:** a crude `mu_{1/7} >= 0.057` uniform bound (three-gap/Riesz), or wire the
  benign `V_0`/arc bound into `LRCWitnessPartA`.
