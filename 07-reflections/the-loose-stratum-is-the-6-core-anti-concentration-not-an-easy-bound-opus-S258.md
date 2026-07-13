---
source: opus-2026-07-11-S258
status: CORRECTION to S257 + refined location. Attempting to prove the loose-stratum anti-concentration bound
  (measure{W=0} > 0 for loose covering families) shows it is NOT easy: the far-peel (wide-gap) lemma is
  rigorous but needs a GIANT runner (V >= ~1190*max(rest), ~0% of covering families); the second moment gives
  the WRONG direction (Chebyshev/Paley-Zygmund UPPER-bound measure{W=0}). Covering ⟹ <=6 effective
  coprime-to-30030 core (99%, verified), so the auto-safe reduction (S241/S243) applies -- but reducing the
  <=6 core to M >= 14/183 IS the <=6-core anti-concentration (S558o even-fold), the project's known hard core.
  So the loose stratum = the <=6-core anti-concentration, NOT a favorable easy bound. Corrects S257's optimism.
tags:
  - lrc14
  - covering-min
  - loose-stratum
  - anti-concentration
  - far-peel-lemma
  - correction
  - six-core
---

# The loose stratum is the ≤6-core anti-concentration, not an easy bound

**opus-2026-07-11-S258.** Owner: prove the loose-stratum anti-concentration bound (S257's second half). Working
it shows the loose stratum is **not** the favorable easy side S257 hoped — it is the project's known ≤6-core
anti-concentration core. Honest correction, with the tool space narrowed.

## The rigorous far-peel lemma — but it needs a giant runner

> **Far-peel lemma (rigorous).** If a sub-family `v` has a safe interval of width `≥ 1/V` at level `c`, then
> `M(v ∪ {V}) ≥ c`. *(That interval contains a full period of runner V's comb, hence one of its safe
> windows.)*

Clean and correct — but **too weak** here. A 12-runner sub-family's safe interval at `c = 14/183` has width
`w ≈ 2(1/13 − c)/max(v) = 2/(2379·max(v))`, so far-peeling `V` needs `w ≥ 1/V`, i.e. **`V ≥ ~1190·max(rest)`**
— a *giant* runner (e.g. `≥ 14 280` for a rest of max 12). **~0%** of natural covering families qualify. The
deep well (`V = 182`, rest `{1..12}`, `w({1..12}) ≈ 0.0003 < 1/182`) is correctly non-peelable (tight). So the
far-peel handles only the giant-runner subclass (S243 Case B), not the loose stratum.

## The second moment gives the wrong direction

For a coprime covering family the danger events are pairwise-independent, so `E[W] = 2ck = 1.989`,
`Var(W) = 2ck(1−2c) = 1.685`. But the second moment only **upper-bounds** the safe measure:
`measure{W=0} ≤ measure{|W−E[W]| ≥ E[W]} ≤ Var/E[W]² = 0.43` (Chebyshev), and Paley-Zygmund lower-bounds the
*danger* support `measure{W≥1}`. **Neither can prove `measure{W=0} > 0`.** The "independent" heuristic
`(1−2c)^13 ≈ 0.116 > 0` is only *pairwise* (the orbit is 1-dimensional, not fully independent), so it is not a
proof. Lovász Local Lemma / Janson fail too: `p = 2c ≈ 0.15` is not small and the 13 arcs are heavily
dependent.

## What actually remains: the ≤6-core anti-concentration

Covering ⟹ the effective **coprime-to-30030 core is small**: verified over 20 664 covering families (speeds <
300), the core-size distribution peaks at 1–2 and is **≤6 for 99%** (only 5 have 7–8, none far-peelable). So
the **auto-safe reduction** (S241/S243: at composite `q`, non-coprime speeds are inert) applies, collapsing the
family to its ≤6 coprime core. **But** reducing that ≤6 core to `M ≥ 14/183` is precisely the **≤6-core
anti-concentration** — the core must not blanket the even-fold good set `G` (S558o) — which is the project's
**known hard core** (S242–S245), *not* an easy bound. The standard tools above do not crack it.

## Net (honest correction)

S257 split the covering-min into [tight: S255 rigidity, proved] + [loose: "favorable" anti-concentration]. This
session **corrects the second half**: the loose stratum is **not** favorable — the far-peel needs a giant
runner (~0%), the second moment is wrong-direction, and (via auto-safe) the loose stratum **reduces to the
≤6-core anti-concentration**, the same hard core the fleet has faced since S242. So the honest state of the
covering-min lower bound after the S253–S258 arc:

- **Tight / deep-well extremizer: PROVED** (S255 via S252 prime-13 uniqueness — the one rigorously-closed
  piece).
- **Loose stratum: = the ≤6-core anti-concentration** (S243/S558o), **open**, and the elementary tools
  (balance-as-bound S256, single dual certificate S257, far-peel + second moment S258) are all ruled out.

The value of the arc is a **precise map of the difficulty**: the covering-min is not crackable by local
witnesses, a single positive-polynomial certificate, wide-gap peeling, or moments; the irreducible content is
the anti-concentration of the ≤6 coprime core against the structured good set. That is where any real proof
must engage — and it is exactly where the fleet's anti-concentration thread already sits.

→ opus-S257 (the split — its loose half corrected here), opus-S255 (deep-well tight case, proved — stands),
opus-S241/S243 (auto-safe + ≤6 core — the reduction), s558o (even-fold good set — the ≤6-core target),
opus-S242–S245 (the anti-concentration core), mac-mini S40 (dual certificate). The far-peel lemma is recorded
as a clean (if narrow) rigorous tool. Files: `lrc14_dual_certificate_knife_edge_split_opus_S257.py` (S257 base);
this session is analysis (inline computations: far-peel threshold, second-moment direction, ≤6-core census).
