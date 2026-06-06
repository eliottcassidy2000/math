---
id: HYP-2272
name: signed-lrc-additive-face-is-multiplicative-energy
status: PARTLY RESOLVED — the structural core is PROVED (THM-414); the completeness sub-questions remain OPEN
date: 2026-06-06
session: monad-explorer-2026-06-06-S704
depends_on: [THM-414, HYP-2262, THM-412, THM-413, THM-401, HYP-1992]
---

# HYP-2272: the signed-LRC additive face = multiplicative energy of the runners' roots of unity; the "popular pair-sum" LRC mirror of density quantization

## The resolved core (now THM-414, PROVED)

The signed-LRC additive face — the pair-sum representation `r_+(s)=#{i<j: a_i+a_j≡s mod C}`,
`C=2n−1` — is the **multiplicative energy** of the roots of unity `ω_i=ζ_C^{a_i}∈μ_C⊂ℚ(ζ_C)`
(`r_+(s)=` coeff of `ζ^s` in `Σ_{i<j}ω_iω_j`); shell-partners are **conjugate pairs** `ω_iω_j=1`;
and in the tangent coordinate `x_i=tan(πa_i/C)` the face is the **spherical formal group**
`F_−(x,y)=(x+y)/(1−xy)` with shell-partner = additive inverse. On the sign hypercube the zero-clock
count is a **pure degree-2 Krawtchouk (`K_2`) observable** (Walsh support = the shell-partner
matching).

**The popular-pair-sum mirror is resolved in the NEGATIVE.** `r_+(s)≤⌊N/2⌋` for all `s`, all
configs (the pairs summing to a fixed `s` form a matching), so the cyclic side has **no popular-norm
escape** — unlike the lattice popular norm `r_Q(D)` which is unbounded (THM-412/S702). The
density-quantization dichotomy is metric-specific: only the Euclidean (multiplicative) norm escapes
its rosette; the cyclic (degree-2 additive) one is capped.

## The open conjectures (what remains)

> **(C1) Completeness of `(E_+, r_+(0))` as a worry-set invariant.** The additive energy
> `E_+=Σ_s r_+(s)^2` is dilation-invariant and, together with the shell-partner count `r_+(0)`,
> separates the exact `n=14` floor: `AP=(328,0)`, `2·AP=(328,0)`, `V*=(290,1)`. **Conjecture:** at
> every `n`, `(r_+(0), E_+)` (or a short list of moments of `r_+`) distinguishes the distinct
> dilation-classes of tight (`M=1/n`) configs. *Status: OPEN; verified only at the `n=14` floor.*

> **(C2) Cyclic image of the rank jump.** To exceed the `⌊N/2⌋` cap one must pass to a higher-degree
> additive form. **Conjecture:** the 3-fold sum `r_3(s)=#{i<j<k: a_i+a_j+a_k≡s}` (or the `d`-fold
> form) is the cyclic mirror of S702's dimension jump from `2D`-lattice (capped) to CM field
> (unbounded power): the first additive form whose max representation grows unbounded in `N` at
> fixed "shape," matching the lattice degree at which the escape begins. *Status: OPEN, speculative.*

> **(C3) The `V*` energy deficit as a carry.** `E_+(AP)−E_+(V*) = 328−290 = 38` at `n=14`.
> **Conjecture:** the deficit equals (or bounds) a carry/owner quantity from the doubling
> shell-partner `(3,24)` (S677 apex debt / HYP-2241), giving a quantitative looseness certificate
> for `V*`-type lifts. *Status: OPEN.*

## Evidence

- THM-414 Parts 1–3 proved; exhaustive sign-pattern checks n≤7; matching cap + dilation invariance
  over 4000 random configs (0 failures). See `signed_lrc_multiplicative_energy_s704.out`.
- `n=14` floor census (C=27): AP `E_+=328, r_+(0)=0`; 2·AP `E_+=328, r_+(0)=0`; V\* `E_+=290,
  r_+(0)=1`.

## Why it matters

Resolves the repeated S699/S702/S703 handoff (the LRC mirror of density quantization), supplies a
second gauge-invariant worry-set separator (`E_+`) beyond the shell-partner count, and gives the
cross-domain dictionary (CM/norm-form ↔ formal group `F_−` ↔ Krawtchouk `K_2`) that makes the
additive face computable from any of four sides.

**See:** THM-414, `07-reflections/the-signed-lrc-additive-face-is-multiplicative-energy-s704.md`,
`04-computation/signed_lrc_multiplicative_energy_s704.py` (+`.out`); HYP-2262, THM-412, THM-413,
THM-401/403/407, HYP-1992.
