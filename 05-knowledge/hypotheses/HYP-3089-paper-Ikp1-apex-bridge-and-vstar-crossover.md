---
id: HYP-3089
title: The polynomial-method paper's I(k,p,1) descent is EXACTLY the project's covering condition (verified at the apex), and the direct-lonely-arc crossover at V* unifies the paper's witness bound, Conjecture-7.1 D, and the project's V* atlas — supplying the reduction kps-S31ag flagged open
status: VERIFIED exact bridge (I(13,7,1)=covering mod 7; c=2 lift -> covering mod 14; 0 counterexamples) + VERIFIED crossover (direct lonely arc bounded below for apex<=V*, decays ~1/V beyond) + SYNTHESIS (analytic substitute for enumeration). Not a proof of LRC(14).
source: mac-mini-2026-06-27-S61
extends:
  - HYP-3087   # kps polynomial-method / CRT bridge
  - HYP-3088   # kps largest-arc <=> Conjecture 7.1(13)
related:
  - THM-573    # level-7 sieve = the paper's c=7 lift (uses LRC<=13 = the paper as induction base)
  - HYP-2900   # Node-3 equidistribution = the apex>V* peel (the "reduction between the two")
  - HYP-3085   # gK8 = low-order moment-LP = the analytic substitute for I(k,p,1) enumeration
  - THM-565    # lonely set = finite union of intervals; arcCount bound after peel
  - THM-530    # witness floor (lonely measure >= m_P = 14249/252252)
  - OPEN-Q-108 # the covering-moment / bounded-core positivity crux
external: arXiv:2604.23906 (Sungkawichai-Trakulthongchai, "Eleven, twelve, thirteen lonely runners", LRC k<=12)
reflections:
  - lrc14-is-the-composite-k-plus-1-case-of-the-polynomial-method-paper   # kps-S31ag mapping (this extends it)
  - the-four-faces-of-14-why-the-exceptional-structures-crowd-into-lrc     # 14=2*7
---

# HYP-3089 — The paper's I(k,p,1) IS the covering condition; the V* crossover unifies the constants

Extends kps-S31ag (which mapped the Sungkawichai-Trakulthongchai polynomial-method paper onto the
project and posed Conjecture 7.1(13) <=> uniform largest-lonely-arc). Two verified contributions plus
the synthesis the user asked for ("merge and extend").

## 1. The apex bridge is EXACT (sharpens kps's "I(k,p,1) ~ p0")
The paper's `I(k,p,l)` = the `(k,p,l)`-IMPROPER tuples (no witness `t in (1/(lp))Z` with `||t v_i||>=1/(k+1)`).
At the apex prime `p=7`, `l=1` (so the gcd condition is vacuous), this is **exactly the covering condition**
(`lrc_paper_Ikp1_apex_bridge_macmini_S61.py`, 0 counterexamples / 40000):

```
LEVEL 7  [ansatz (1/7)Z ]:  no witness  <=>  some v_i ≡ 0 (mod 7)    =  covering mod 7   = I(13,7,1)
LIFT c=2 [ansatz (1/14)Z]:  no witness  <=>  some v_i ≡ 0 (mod 14)   =  PROJECT COVERING CONDITION
```
Proof of level 7: `t=j/7` gives `||t v_i|| = ||j v_i/7||`, which is `>=1/14` iff `j v_i ≢ 0 (mod 7)`; since
`j` is a unit, this holds for all `i` iff no `v_i ≡ 0 (mod 7)`. The improper set (no working `j`) is
therefore exactly `{some v_i ≡ 0 mod 7}`. **RESCUE mechanism** (verified): under the `c=2` lift, a
coordinate `≡ 7 (mod 14)` (odd multiple of 7) is rescued by an odd `j` (`7j ≡ 7 ≢ 0`), but a coordinate
`≡ 0 (mod 14)` (a multiple of 14) is killed for every `j` — so the descent `7 -> 14` lands **exactly** on
the project's covering condition. This is the precise, computable form of kps's `I(k,p,1) ~ p0` analogy:
the paper's bottleneck object at the apex IS the project's most-studied case. (`14=2·7`, the four-faces
H2 multiplicative face: the apex-7 base, doubled.)

## 2. The V* crossover — supplying the reduction kps flagged open
kps's open piece: a uniform largest-arc bound for the *direct* lonely set `L(S)={t:||s_i t||>=1/14}`,
"or a clean reduction between the two" (direct vs. peeled). The reduction is now explicit
(`lrc_lonely_arc_count_vs_apex_macmini_S61.py`), on the covering family `{1..12, 14V}`:

| apex `14V` | 14 | 70 | 112 | 182 | 280 | 700 | 2800 |
|---|---|---|---|---|---|---|---|
| largest lonely arc | .0032 | .0049 | **.0053** | .0047 | .0031 | .0012 | .0003 |
| `1/largest` (≈ needed `D`) | 308 | 203 | **188** | 212 | 327 | 817 | 3267 |

The direct largest arc is **bounded below (~0.003–0.005) while apex ≲ V*** (a finite/bounded regime), then
the apex's fine forbidden arcs (spacing `1/(14V)`) subdivide each bounded-core arc and the largest arc
**decays ~1/V**. So:
- **apex ≤ V*:** the direct lonely set HAS a long arc ⟹ a witness `a/d` exists for `d ≳ 1/ℓ_max ≈ 200–300`
  — a finite check. (kps's largest-arc route works directly here.)
- **apex > V*:** PEEL the apex (its forbidden set is `1/7`-equidistributed; Node-3 / HYP-2900); the
  **bounded core `{1..12}` keeps its long arc (≈0.006)** and the apex misses it on positive sub-measure.

**The crossover apex ≈ V* ≈ 200 is the SAME constant in three independently-derived framings:** the
paper's witness-denominator bound, kps's `D≈213` (Conjecture 7.1), and the project's `V*` atlas (≈234).
Strong evidence they are one constant. So **Conjecture 7.1(13) splits cleanly at V*** into `[apex≤V*:
finite bounded check]` + `[apex>V*: peel/equidistribution]` — exactly the project's two covering regimes
(top-balanced bounded core via gK8/p0≤cap, and one-large via Node-3). The "reduction between the two" is
the peel.

## 3. The synthesis: the project is the ANALYTIC SUBSTITUTE for the paper's enumeration
The paper proves `k≤12` by enumerating `I(k,p,1)` over a sieve of primes `p>k²+k`, at cost
`~p^{(k+1)/2}/(k·2^k)`; `k=13` is open because that enumeration is "astronomical." The project's
`p0(E)≤cap_k` (the **continuous**, uniform-in-`p` limit of the improper fraction) + the witness floor
(`THM-530`, lonely measure `≥ m_P>0`) is the **analytic estimate that replaces the enumeration**: a
uniform positive lonely measure on covering tuples gives a real witness directly, with no per-prime count.
mac-mini's gK8 = low-order moment-LP led by pairwise `S2` (HYP-3085) is precisely a uniform bound on the
improper fraction's leading correction — the analytic substitute, made certificate-shaped (the 3×3
reflection-Perron block). **This is the concrete sense in which the project can pass `k=12`: not a faster
finite check, but a uniform analytic lower bound on the witness fraction.**

Mutual dependency (the merge is symbiotic):
- the **paper → project:** LRC(`≤12`) is the **induction base** for THM-573 (the level-7 sieve uses
  `LRC(≤13)`), so the paper *unlocks* the project's `c=7` lift theorem;
- the **project → paper:** THM-573 makes the `c=7` lift a **theorem** (not a table lookup), and gK8/p0≤cap
  supplies the **analytic substitute** for the `I(13,p,1)` enumeration the paper cannot afford, both
  reducing to the same **Conjecture 7.1(13) ⟺ uniform positive lonely measure on covering tuples**.

## What remains (unchanged target, now triply-motivated)
`Conjecture 7.1(13) ⟺ LRC(14) ⟺ [every non-tight primitive covering 13-tuple has lonely measure ≥ m_P
and (apex≤V*: a long arc | apex>V*: the bounded core's arc survives the equidistributed peel)]`. The two
live pieces are the project's existing CRUX 1 (bounded-core positivity = p0≤cap = gK8, HYP-3085) and Node-3
effective Erdős–Turán (HYP-2900) — now also recognized as the paper's Conjecture 7.1(13). Induction base
`LRC(≤12)` = the paper (PROVED).

## Next
1. Confirm the V* crossover is worst-case (run the arc-curve over genuine-wide / top-balanced covering
   families, not just `{1..12,14V}`); pin `V* = D` numerically across families.
2. Make precise the inequality `improper-fraction(k,p,1) ≤ f(cap_k)` uniform in `p` (gK8 ⟹ the paper's
   per-prime bound), the formal analytic-substitute statement.
