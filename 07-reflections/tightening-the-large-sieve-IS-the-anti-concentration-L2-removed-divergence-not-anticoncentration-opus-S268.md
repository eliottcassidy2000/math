---
source: opus-2026-07-11-S268
status: HONEST NEGATIVE on "close the 3.1x gap" + two structural advances + a correction of S267. Tightening the
  large-sieve energy bound does NOT close the gap, because the TIGHT large-sieve bound IS the anti-concentration
  (eps_v = O(1) uniformly). ADVANCES: (1) the gap is NOT uniform -- it is isolated to core=1 (the deep well =
  runner-1 lemma S265); core>=2 has a ~24x margin. (2) the scaling is eps_v=O(1), NOT sqrt(|G'|): max|eps_v|<=0.099
  stays flat as |G'|->0.2 while Bessel allows sqrt(6/49/|G'|)~0.78. CORRECTS S267: L2 removed the L1 DIVERGENCE
  (S266's alternating cancellation), NOT the anti-concentration, which reappears as the |G'|-scaling / the positive
  off-diagonal. Cleanest form of the irreducible core: eps_v = O(1)  <=>  <g(v.),1_G'> = O(|G'|).
tags:
  - lrc14
  - covering-min
  - L2-energy
  - large-sieve
  - anti-concentration
  - honest-negative
  - corrects-S267
---

# Tightening the large-sieve bound IS the anti-concentration: L² removed the divergence, not the anti-concentration

**opus-2026-07-11-S268.** Owner: tighten the large-sieve energy bound to close the 3.1× gap (S267). It does
not close — and finding out why is the value. The tight large-sieve bound *is* the anti-concentration. But the
attempt sharpens the structure decisively and corrects S267's over-optimism.

## Advance 1 — the 3.1× gap is not uniform; it is isolated to core=1

Splitting `core·Σ_v ε_v²` by the number of core (coprime-to-30030) speeds:

| core count | max `core·Σε²` | vs target 0.735 |
|---|---|---|
| **1** | **0.328** | the only case near threshold |
| 2 | 0.020 | 37× margin |
| 3–7 | ≤ 0.030 | ≥ 24× margin |

The worst case — `0.328` — is **core=1**, the deep well / near-AP family, where speed 1 is the *only* coprime
speed. There Cauchy–Schwarz is **exact** (a single term), so `core·Σε² = ε_1²` and the requirement is just
`ε_1 < 6/7`: this is precisely the **runner-1 lemma (S265)**, handled by measure ∪ equidistribution — *not* a
large-sieve object at all. For **core ≥ 2** the energy has a **~24× margin**. So the "gap" is not a uniform
3.1×; it is one hard family (the deep well) plus a wide-margin remainder.

## Advance 2 — the scaling is ε_v = O(1), not √|G'|

For core ≥ 2, `max|ε_v| ≤ 0.099` stays **flat** as `|G'|` shrinks to 0.2, while the Bessel/operator-norm bound
only gives `|ε_v| ≤ √(6/49 / |G'|)`, which **grows to 0.78**. So the truth is

> `⟨g(v·), 1_{G'}⟩ = ε_v·|G'| = O(|G'|)` — **linear** in `|G'|`, i.e. `ε_v = O(1)`,

whereas Bessel permits `√|G'|`-decay only, i.e. `ε_v = O(1/√|G'|)`. The Bessel bound is therefore loose by
exactly `~1/√|G'|` — which *is* the 3.1× at `|G'|≈0.27`, and worse for smaller `|G'|`. The same fact shows up
in the energy expansion as a **positive off-diagonal** (~64% of the total for core ≥ 2): the diagonal (= Bessel)
**underestimates** because `1_{G'}` is far from the frame operator's top eigenvector, and closing *that* gap is
the anti-concentration.

## Advance 3 (the correction) — L² removed the divergence, not the anti-concentration

S267 claimed the L² route "sidesteps the multi-linear cancellation." That is **half right**, and this session
pins down which half:

- **What L² removed:** the **L¹ divergence** — the harmonic tail `Σ|b_k| = ∞` and S266's alternating
  multi-order cancellation. The L² energy `Σε²` genuinely converges. Real progress.
- **What L² did NOT remove:** the **anti-concentration** itself. A *tight* bound on the convergent energy still
  requires `ε_v = O(1)` uniformly, i.e. `⟨g(v·),1_{G'}⟩ = O(|G'|)` — and that is exactly the anti-concentration
  statement. It reappears, no longer as a divergence, but as the **|G'|-scaling** (and the positive
  off-diagonal). Bessel cannot deliver it because the operator norm only sees `‖1_{G'}‖² = |G'|`, not that
  `1_{G'}` decorrelates from the core dilates.

So "tighten the large sieve" and "prove the anti-concentration" are **the same task**, not sequential ones.

## Net — the cleanest statement of the irreducible core

Across three sessions the same hard core has appeared in three costumes:

- **S266:** "multi-linear inverse theorem" — the L¹ / alternating-cancellation costume.
- **S267:** "tight large-sieve energy bound" — the L² costume (over-optimistic: L² converges but the tight
  bound is still anti-concentration).
- **S268 (this):** its cleanest, most quotable form — a **uniform L∞ bound on the density deviation**:

> **`ε_v = O(1)` uniformly over covering families and core speeds `v`**  ⟺  **`⟨g(v·), 1_{G'}⟩ = O(|G'|)`**,

realized most extremely at the **deep well (core=1)**, which is the **runner-1 lemma (S265)**. For **core ≥ 2**
the bound holds numerically with a ~24× margin (`ε_v ≤ 0.099`), but a rigorous proof of even that still needs
the same `ε_v = O(1)` input. This is a strictly sharper localization than S266/S267: the entire difficulty of
LRC(14) for covering families is a single uniform O(1) bound on `ε_v`, hardest at one explicit family. The
honest verdict on the owner's task: the gap does not close by tightening the large sieve, because the tight
large-sieve bound and the anti-concentration are one and the same.

→ opus-S267 (corrected — L² removes divergence, not anti-concentration), opus-S266 (same core, L¹ costume),
opus-S265 (core=1 = runner-1 lemma, the hard family), opus-S262 (`Cov ≤ 1/(3vv')` — the rigorous-but-loose
Bessel input). Files: `lrc14_L2_tightening_hits_anticoncentration_opus_S268.py` (+`.out`).
