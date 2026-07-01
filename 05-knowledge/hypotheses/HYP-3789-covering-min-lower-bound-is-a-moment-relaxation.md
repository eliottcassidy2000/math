---
id: HYP-3789
title: THE COVERING-MIN LOWER BOUND IS A MOMENT RELAXATION (Lasserre/trig-moment), and its finite levels ARE the repo's existing threads -- level-1 = union bound / Fejer average (HYP-3785 spectral gap, fails), level-2 = pair correlations = the 2nd-moment floor (HYP-3571) = the S75 signed correction (HYP-3787), level-infinity = exact = the lazy-cut (HYP-3782). The extremal lonely set is a FEW ATOMS (measure zero: the construction's binding t*=n/Phi6 and its iota-image; the AP {1..13}'s six atoms = the units (Z/14)*), so the lonely measure's Toeplitz moment matrix is LOW-RANK PSD (rank = #atoms) and the covering-min lower bound is a FLAT-EXTENSION (Curto-Fialkow) certificate. MECHANISM (level-2, exact-on-grid): the small-speed dangers are POSITIVELY correlated (pile up redundantly near rationals, S_2 excess +0.60 over independence) so they waste coverage and leave t* lonely, while the far element n(n-1) is ~INDEPENDENT (I_{v,182} = (2r)^2 to 3 digits) and equidistributes -- the level-2 moment DIRECTLY sees S74's near/far split. Because the extremum is atomic, no finite Lasserre level is exact = the moment-hierarchy face of the S61 'finite-D certificates cannot close Step 3' rigidity and of HYP-3785's spectral gap.
status: CONFIRMED (computational, exact-on-grid + exact-Fraction atoms + exact Toeplitz rank). A UNIFYING REFRAME + one new mechanism (level-2 positive correlation of the near core = why loneliness survives) + one new invariant (the extremal lonely set is finitely atomic, rank = #atoms: 2 for construction, 6=phi(14) for AP). NOT a new proof of the covering-min; it identifies precisely what each existing lever is (a level of one hierarchy) and localizes the residual (>=7-huge-speed = a level-2 cross-harmonic moment estimate). Honest: the flat-extension certificate is rank-1 (a single witness point), so the moment view confirms the difficulty is the SEARCH over covering S, not the certificate degree.
source: mac-mini-2026-06-30-S76
related:
  - HYP-3784   # klein Delsarte dual near-tight at bounded V (= finite moment/LP level); weakens as V grows
  - HYP-3785   # klein SDP/Fejer spectral gap (averaged relaxation blind to the pointwise spike) = the finite-level gap
  - HYP-3571   # the 2nd-moment floor 1/(2 zeta(2)) = 3/pi^2 = level-2 of this hierarchy; extremal lonely set = units
  - HYP-3549   # extremal = AP {1..n-1}; lonely set = units mod n in phi/2 antipodal pairs (confirmed here: 6 atoms)
  - HYP-3787   # S75 signed correction = the pair-correlation (level-2 Fourier) moment; far element ~independent
  - HYP-3788   # S74 equidistribution on L_C; the near/far split is the level-2 correlation structure
  - HYP-3782   # lazy-cut = level-infinity (exact Positivstellensatz cuts)
  - HYP-2948   # Beurling-Selberg minorant = the dual extremal certificate of the LP relaxation
results:
  - 04-computation/moment_relaxation_covering_min_macmini_20260630.py
  - 05-knowledge/results/moment_relaxation_covering_min_macmini_20260630.out
---

# HYP-3789 -- the covering-min lower bound is a moment relaxation

The owner's seed -- "moment relaxation and covering-min" -- names the frame that unifies the repo's
covering-min lower-bound levers. The inner problem `M(S) = max_t min_v ||vt||` is a polynomial (here
trigonometric) optimization over the circle, so its lower bound has a **Lasserre / trigonometric-moment
hierarchy**, and the existing threads are exactly its levels.

## The hierarchy (inclusion-exclusion = Bonferroni = moment truncation)

    |L_S(r)| = integral_0^1 prod_{v in S} (1 - g_v(t)) dt,   g_v(t) = 1_{||vt|| < r}
             = sum_{A subseteq S} (-1)^|A| I_A,   I_A = meas( intersection_{v in A} danger_v )

Grouped by level `k`: `S_k = sum_{|A|=k} I_A = integral C(c(t),k) dt` where `c(t) = #{v: ||vt||<r}`; the
level-`m` truncation `T_m = sum_{k<=m} (-1)^k S_k` brackets `|L_S|` (Bonferroni). The levels:

| level | object | repo thread | status |
|------|--------|-------------|--------|
| 1 | `T_1 = 1 - |S|*2r` (union bound / Fejer average) | HYP-3785 spectral gap | **FAILS** (`= -0.97 < 0`) |
| 2 | `S_2 = sum pair-overlaps` (pair correlations) | 2nd-moment floor HYP-3571; S75 signed correction HYP-3787 | the mechanism (below) |
| ... | higher correlations | joint equidistribution HYP-3788 | Bonferroni oscillates |
| inf | exact `|L_S|` | lazy-cut HYP-3782 (Positivstellensatz) | exact |

Verified exact-on-grid (`G=2^22`, `r=0.99 M`, construction and AP): `|L_S|>0` (`0.0019`, `0.0026`); the
level-1 union bound `T_1 = -0.97 < 0` is **useless** (confirming HYP-3785: the averaged/PSD relaxation is
blind to the pointwise spike); the truncations `T_m` **oscillate** (`+1.4, -3.0, +5.4, -7.3, ...`) and
converge only at full level `m = |S|`. This IS "no finite level certifies."

## The extremal lonely set is a FEW ATOMS => flat-extension certificate

At `r = M` (exact Fraction), the lonely set (global maximizers of `G(t)`) is finite (measure zero):
- **construction `{1..12,182}`**: exactly **2 atoms** `{14/183, 169/183}`, iota-symmetric (`t <-> 1-t`),
  denominator `Phi6 = 183`.
- **AP `{1..13}`**: exactly **6 atoms** `{1,3,5,9,11,13}/14` = the **units `(Z/14)*`** (independent
  confirmation of HYP-3571/HYP-3549: `phi(14)=6` in `3` antipodal pairs).

So the lonely measure `mu` (uniform on the atoms) has a **Toeplitz moment matrix** `M_K=[y_{j-l}]`,
`y_k = (1/|A|) sum_a e^{-2pi i k a}`, that is **PSD** (Bochner, min-eig `~ -1e-14`) with **rank = #atoms**
(Caratheodory-Toeplitz): rank stabilizes at `2` from `K>=1` (construction) and at `6` from `K>=5` (AP).
Localization verified: every atom is lonely for every `v` (`||va|| >= M`). Hence the covering-min lower
bound is a **flat (finitely-atomic) moment problem** -- a Curto-Fialkow flat-extension certificate.

## The mechanism (the new content): level-2 positive correlation of the near core

The pair correlations `I_{vw}` (level-2 moment) split cleanly by arithmetic (exact-on-grid, `r=M`):
- **small-speed pairs** `(1,2),(1,3),(6,12),(4,8)`: `I_{vw} = M` (`+0.053` **above** independence `(2r)^2`)
  -- strongly **POSITIVELY correlated**; their danger arcs pile up redundantly near rationals.
- **far element** `182`: `I_{v,182} = (2r)^2` to 3 digits -- **~INDEPENDENT** of everything; it
  equidistributes (this is exactly S74/HYP-3788's near/far split, seen at level 2).

Full `S_2 = 2.43` vs independence `C(13,2)(2r)^2 = 1.83`, **excess `+0.60`**. So:

> **The covering-min lower bound survives because the small-speed dangers are positively correlated:
> they waste their coverage budget piling up redundantly near the rationals, leaving the binding time
> `t*` uncovered, while the far speed `n(n-1)` equidistributes (independent) and cannot finish the cover.**

This is the level-2 (pair-correlation) *mechanism* for "no beater," quantitative and exact-on-grid.

## Integration with the repo's moment/SDP machinery (S76 survey)

The primal object here (the inclusion-exclusion hierarchy on `|L_S|`) has a **dual** already built:

- **DUAL (proved)**: THM-534 (`the sector moment-LP dual certificate`) -- `meas(S_7) <= L_y(E)`, the magic
  function `g(t) >= 1_{t=0}`; the project's core LP, AP-extremal for bounded spread. klein's Delsarte LP
  (`covering_min_delsarte_lp_klein.py`, HYP-3784) is its speed-restricted instance.
- **DUAL SOS (verified)**: the `Z_7` cyclotomic **Fejer-Bochner SOS** (opus,
  `lrc_z7_cyclotomic_sos_floor_opus_20260629.py`): kernel `h` PSD with eigenvalues
  `rho_j = 1_safe(j/7) = (0,1,1,1,1,1,1)`, floor `c=1` on the six safe modes -- the *actual* SOS
  certificate of the level-2 floor, set-independent (apex-7 transitivity).

**Crucial caveat (resolves an apparent tension).** A naive Lasserre lift does NOT help:
- My own `lrc14_threadC_sos_lasserre_lift_macmini.py` found the **degree-2 Lasserre lift of the discrete
  mod-7 sector-hit distribution COLLAPSES to level-1** (CJJ Prop 1.2: missing-sector events aren't closed
  under linear combination).
- opus's `lrc_lasserre_sdp_Mofs_opus_20260701.py` (M(S) moment-SOS, exact on `{1,2,3}->1/4`) is the
  **wrong direction**: it upper-bounds `max_t min_v`, i.e. it certifies *beaters*, not the lower bound;
  and `M(S)` is computable exactly anyway.

The resolution sharpens the hierarchy picture: the **discrete** (mod-`p`) Lasserre collapses, but the
**continuous** 2nd moment -- the pair correlation of the danger *arcs* -- does NOT: it is the HYP-3571
floor `1/(2 zeta(2)) = 3/pi^2 = 0.304`. So the moment relaxation that carries content is the
**inclusion-exclusion on the continuous lonely measure** (this file), with dual = THM-534 LP + the `Z_7`
Fejer-Bochner SOS, NOT a generic sector-hit SDP.

**Loose inspiration (survey)**: the LRC magic function is the Cohn-Elkies / Chebyshev-equioscillation
extremal (`V_7(u)-2 = (u-2)(de Moivre cubic)^2`, kps reflection), the 1-D analog of Viazovska's
sphere-packing magic function -- the same "positive-definite minorant meets the extremal configuration"
duality, here cyclotomic (`E_2`/hexagonal) rather than modular.

## What it buys and does not buy

- **Unification**: the union bound (HYP-3785), the 2nd-moment floor (HYP-3571), the signed correction
  (HYP-3787), the Delsarte LP (klein HYP-3784), and the lazy-cut (HYP-3782) are one hierarchy, levels
  `1, 2, 2, 1-2(LP dual), inf`.
- **Residual localized**: the S75 `>=7`-huge-speed case is a **level-2 cross-harmonic moment estimate**
  (the off-diagonal Toeplitz entries `hat(1_{L_C})(jw_i - j'w_j)`).
- **New invariant**: the extremal lonely set is finitely atomic, rank `=` #atoms (`2`, `6=phi(14)`).
- **Honest limits**: the flat certificate is rank-1 (one witness point suffices), so the moment view
  *confirms* that the hard part is the `for all covering S` search, not the certificate -- the
  moment-hierarchy face of S61 ("finite-D certificates cannot close Step 3") and of HYP-3785's gap.
  cvxpy/SDP solver unavailable, so the full Lasserre SDP value is not computed here (klein's Fejer proxy
  HYP-3785 stands for the level-1 dual). NOT a new proof; a clarifying reframe + one mechanism + one
  invariant.
