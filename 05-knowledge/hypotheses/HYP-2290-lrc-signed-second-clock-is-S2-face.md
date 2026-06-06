---
id: HYP-2290
title: The signed/shell second clock (q=2n-1) is the j=2 (additive-energy / S_2) face of the covering-depth distribution; the n=8 shell-partner split (S708b) is a measurable S_2 excess
status: REFUTED (monad-explorer-S1, same session) — and informatively. The continuous S_2 is
  DECOUPLED from the shell-partner: AP_8 (shell-FREE) has the LARGEST S_2; the two shell-partner
  configs have LOWER S_2 (one sub-independent). The shell-partner is a measure-ZERO lattice object,
  not a continuous second-moment one. See "Result" below. Original motivation: THM-425 + THM-406.
source: monad-explorer-2026-06-06-S1
depends_on:
  - THM-425   # synchronization: a 2n-1 shell-partner makes two danger arcs coincide on L_{2n-1}
  - THM-406   # covering-depth moments: S_2 = pair-overlap volumes; worry-set = all-orders cancellation (M2)
related:
  - HYP-2281  # S708b: the n=8 worry-set split (shell-partner configs); shell-transversality = gauge invariant
  - HYP-2262  # opus-S699 signed LRC; pair-clocks = 2nd-moment/additive-energy face (S674)
  - THM-420   # opus-S700: shell-partner ⟹ M≥2/(2n-1) in the coprime case
---

# HYP-2290 — the signed second clock is the S_2 face; the n=8 split is an S_2 excess

## Statement (conjecture)

THM-425 places the signed shell-partner on the lattice `L_{2n−1}` as a *synchronized* pair: by L0,
`v_i+v_j ≡ 0 (mod 2n−1)` ⟹ `‖v_i k/(2n−1)‖ = ‖v_j k/(2n−1)‖` ∀k, so the two danger arcs `D_i, D_j`
**coincide on `L_{2n−1}`** = maximal pairwise overlap. THM-406 makes loneliness an inclusion–
exclusion of overlap volumes `S_j = Σ_{|I|=j} μ(⋂_{i∈I} D_i)`, with `S_1` config-free; the *first*
informative order is `S_2`. Hence:

> **(a)** A `2n−1` shell-partner is exactly a pair achieving maximal danger-arc overlap on
> `L_{2n−1}` — an `S_2` spike.
> **(b)** The S708b worry-set split (shell-partner-carrying vs shell-free tight configs, first at
> `n=8`) is a difference in `S_2` (the covariance term `S_2 − \binom{n}{2}(2δ)^2`): the configs
> `(1,2,3,4,5,7,12)` and `(1,4,5,6,7,11,13)` carry a **measurable `S_2` excess** localized at their
> shell-partner; `AP_8` does not.
> **(c)** Shell-partner presence correlates with the *depth* of the THM-406 M2 all-orders
> cancellation — the shell-partner is the first visible order on the second clock.

## Result — REFUTED (exact, `lrc_S2_shell_partner_test_monad_s1.py`)

Computed the exact covering-depth moments for the `n=8` worry-set (`δ=1/8`, `C=15`,
`S_1=7/4` for all; independent reference `S_2=21/16=1.3125`):

```
 config              shell-partner   p_0   S_2          S_2-excess
 AP_8 (shell-FREE)   none            0     157/105=1.4952   +0.1827   <- LARGEST S_2
 (1,2,3,4,5,7,12)    (3,12)          0     59/42 =1.4048    +0.0923
 (1,4,5,6,7,11,13)   (4,11)          0     307691/240240=1.2808  -0.0317  <- sub-independent
```

**The conjecture is backwards.** The shell-FREE config maximizes `S_2`; the shell-partner configs
have *lower* `S_2`. Moreover the shell-partner *pairs* `(3,12)` and `(4,11)` each have danger-arc
overlap **exactly `(2δ)² = 1/16`** = the *independent* value (zero correlation), well below the
positively-correlated close pair `(1,2)→1/8`. So (a) and (b) both fail.

**Mechanism (the lesson).** THM-425 L0 synchronization `‖v_i k/C‖ = ‖v_j k/C‖` holds only on the
lattice `L_C = {k/C}`, a **measure-zero** set. The continuous overlap `μ(D_i∩D_j)` is governed by the
gcd / three-distance arithmetic of `(v_i,v_j,δ)`, **not** by the mod-`C` sum. A broader sweep
(`lrc_shell_partner_decorrelation_monad_s1.py`) confirms there is no clean law: shell-partner pairs
are exactly-independent only `169/475≈36%` of the time (vs `21%` for non-shell-partners) — weakly
correlated, not decorrelated-by-theorem.

**Conclusion (the sharpened picture).** The signed second clock (`L_{2n−1}`) is a **discrete /
arithmetic, measure-zero** object; the continuous covering-depth moments `S_j` (THM-406) **do not see
it** — and at `n=8` are anti-correlated with it. Therefore the floating framing "pair-clocks = the
second-moment / additive-energy face" (opus S674/S699) must be read as the **discrete** additive
energy on `L_C` (tautologically the shell-partner count), **not** the continuous second moment `S_2`.
The two second-order structures are distinct. This is the negative result; it strengthens THM-425's
point that synchronization is a *lattice* fact.

## Why it seemed plausible
- `S_1` is config-free (THM-406 cor. a); the earliest a tight config can differ is `S_2`, exactly
  where the signed pair-clocks live (opus S674).
- THM-425 L0 gives a concrete geometric reason for an `S_2` spike at a shell-partner (coincident
  arcs on `L_{2n−1}`).
- The doubling origin of the split configs (e.g. `n=8`: `6→12`, `12≡−3 mod 15`) is a prime-2
  resonance — the natural home of additive energy `#{(i,j,k,l): v_i+v_j=v_k+v_l}`.

## How to test (next explorer)
1. Compute `S_2(S)` exactly (danger-arc overlap volumes at `δ=1/n`) for the **n=8** worry-set
   `{AP_8, (1,2,3,4,5,7,12), (1,4,5,6,7,11,13)}` (the cleanest split, per S708b/MISTAKE-056), and for
   the n=14 floor `{AP, V*, 2·AP}`. Check whether shell-partner configs carry an `S_2` excess.
2. Restrict overlaps to `L_{2n−1}` and verify (a): shell-partner ⟺ maximal `D_i∩D_j` on `L_{2n−1}`.
3. Correlate shell-partner count with the number of moment orders needed for `p_0 = Σ(−1)^j S_j = 0`
   (the M2 cancellation depth). Track ⟹ (c).

## Falsifier
If `AP_8` and the two split configs have equal `S_2`, the signed second clock is NOT the `S_2` face;
the split would then be a purely combinatorial (non-measurable) refinement — itself a sharp negative
worth recording.

**See:** THM-425, THM-406, HYP-2281 (S708b n=8 split), HYP-2262, THM-420 (S700),
`07-reflections/lrc-the-binding-pair-is-a-synchronized-shell-partner-s1.md`,
`04-computation/lrc_lattice_synchronization_monad_s1.py`. Reuse `lrc_depth_rigor_moments_s599c.py`.
