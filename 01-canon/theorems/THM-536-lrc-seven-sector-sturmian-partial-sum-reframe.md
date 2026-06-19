---
id: THM-536
title: The Sturmian partial-sum-walk reframe of the LRC(14) seven-sector cover for the AP — meas(S7(consec_k)) = (1/7)·meas{θ∈[0,7): the partial sums S_e=⌊eθ⌋ mod 7 (e=0..k-1) of a slope-frac(θ) mechanical word cover Z/7}; PROVED vanishing meas(S7(consec_k))=0 for k≤6 (only k<7 residues visited); the PROVED pointwise subset-domination lemma (E⊆{0..N} ⟹ meas(S7(E))≤meas(S7(AP_{N+1}))); and three clean refutations (per-block IE-extremality, span/spread-monotonicity, and meas(S7) translation-invariance are all FALSE) that explain WHY AP-extremality of meas(S7) is irreducibly an aggregate phenomenon
status: MIXED. PROVED: the θ=7x Sturmian reparametrization (exact, verified k=3..14); meas(S7(consec_k))=0 for k≤6 (cardinality: |{S_e}|≤k<7); the pointwise subset-domination lemma (set inclusion at every x); meas(S7) scale-invariance (re-confirmed) and NON-translation-invariance (counterexample). VERIFIED (exact/adversarial): THM-534's deg-3 dual L_y closes k=9,10 with consec the maximizer and ZERO over-cap on widened boxes (maxE≤17, 24310 sets) — an independent adversarial confirmation of THM-534. REFUTED (exact): per-IE-block extremality of meas(B_M) (AP is neither the per-M max nor min; ~half the signed per-M differences are negative); span-monotonicity and elementary-spread-monotonicity of meas(S7) (both have increases); 7-term-AP ⟹ 31/210 (false for translated/non-0-based APs — only scale-orbits of {0..6} give 31/210). LRC(14) NOT proved; the live gap is unchanged ("consec maximizes meas(S7)"/L_y on k=8,9,10).
source: mac-mini-2026-06-18-S7
depends_on:
  - THM-532   # seven-sector relation-height split (meas(S7)=M7+corr; scale invariance)
  - THM-534   # the moment-LP dual certificate (this gives the independent adversarial check of L_y)
  - THM-535   # the subadditive cap-split (reduced the finite check to k=8,9,10 — exactly the rows this reframes)
  - HYP-2604  # the AP-frontier envelope (the meas(S7) extremality this reframes)
related:
  - THM-531   # AP-orbit invariance (mu_theta is translation+scale invariant — CONTRAST: meas(S7) is only scale-invariant)
  - HYP-2603  # codex seven-sector net-cap reduction
external: Lonely Runner Conjecture (first open case = 13 speeds). Sturmian / mechanical words; three-distance (Steinhaus) theorem; cutting sequences of a line.
---

# THM-536 — The Sturmian partial-sum-walk reframe of the seven-sector cover

This is **Angle B** (symbolic dynamics / cutting sequences) applied to the LRC(14) crux
`meas(S7(E)) ≤ cap_k`. It re-expresses the seven-sector cover of the **consecutive cluster**
(the conjectured extremiser, HYP-2604 / THM-535's last open rows) as a **combinatorics-on-words**
object, proves two clean sub-results, and cleanly **refutes** three plausible routes to the
AP-extremality — pinning down exactly why that extremality is hard.

## A. The reparametrization (PROVED, verified k=3..14)

The sector index `σ_e(x) = ⌊7ex⌋ mod 7 ∈ Z/7` is the **cutting sequence** of the slope-`7e`
line through the unit grid. Substitute `θ = 7x ∈ [0,7)` (measure-preserving up to the `1/7`
factor). Then for the AP `E = {0,1,…,k−1}`,

> `σ_e(x) = ⌊e·θ⌋ mod 7`,   so   `meas(S7(consec_k)) = (1/7)·meas{ θ∈[0,7) : {⌊e θ⌋ mod 7 : e=0..k−1} = Z/7 }`.

The residues `S_e := ⌊e θ⌋ mod 7` are the **partial sums** of the increment word
`d_e := ⌊(e+1)θ⌋ − ⌊e θ⌋ ∈ {⌊θ⌋, ⌊θ⌋+1}` — a **mechanical (Sturmian) word** of slope
`frac(θ)`, alphabet `{j, j+1}` on `θ∈[j,j+1)`. So the AP's seven-sector cover is exactly the
**Sturmian cover-time problem**: *for which θ do the length-k partial sums of the slope-frac(θ)
mechanical word, read mod 7, hit all seven residues?* (Verified: the Sturmian-cover engine
matches the breakpoint engine exactly for all `k=3..14`.)

Exact values (in `x`): `meas(S7(consec_k))` = `0` for `k≤6`, then `31/210, 481/1470,
2447/5880, 8899/17640, 3419/5880, 121103/194040, 14573/21560, 14109/20020` for `k=7..14`.

## B. Two PROVED sub-results

**B1 — Vanishing below k=7 (PROVED, trivial via the walk).** The walk visits at most `k`
distinct partial sums `{S_0,…,S_{k−1}}` (including `S_0=0`). For `k≤6` that is `<7` values,
so `Z/7` is never covered: **`meas(S7(consec_k)) = 0` for all `k≤6`.** (The first nonzero
value is `31/210` at `k=7`, where a full residue *permutation* is needed.) This is the exact
sector analogue of "you need at least 7 runners to surround the lonely one with a `1/7`-net."

**B2 — Pointwise subset domination (PROVED, set inclusion).** If `E ⊆ {0,1,…,N}` then at
**every** `x`, `{σ_e(x):e∈E} ⊆ {σ_e(x):e∈{0..N}}`, so `S7(E) ⊆ S7({0..N})` as *sets*. Hence

> `meas(S7(E)) ≤ meas(S7(AP_{N+1}))`   whenever `E ⊆ {0,…,N}`.

This rigorously certifies `meas(S7(E)) ≤ cap_k` for every primitive `E` of span `≤ N*(k)`,
where `N*(k)` is the largest `N` with `meas(S7(AP_{N+1})) ≤ cap_k`: `N*= 7,8,10,13,20,20` for
`k=8..13`. It is **sharp at the AP** (span `k−1`) but weak for small `k` (at `k=8`, `N*=7=k−1`
certifies only the AP itself); the residual is the larger-span shapes, which are empirically
far under cap (`≈0.20` vs `cap_8≈0.38`).

## C. Three REFUTATIONS — why AP-extremality is irreducibly aggregate

The natural "structural" proofs of `meas(S7(E)) ≤ meas(S7(consec_k))` all **fail**:

- **C1 — per-IE-block extremality is FALSE.** Writing `meas(S7) = Σ_{M⊆Z/7}(−1)^{|M|}meas(B_M)`
  (B_M = avoid the missed-residue set M), the AP is **neither** the IE-correct max nor min of
  the individual `meas(B_M)`: at `k=8`, AP is the per-M extremiser on `0/6` size-1, `0/15`
  size-2, `0/20` size-3 blocks; for the top competitors `~half` of the signed per-M differences
  `D_M(E)=(−1)^{|M|}(meas(B_M(AP))−meas(B_M(E)))` are **negative**. So AP-extremality survives
  **only after summing all 64 IE terms** — there is real sign cancellation; no term-by-term
  (Bonferroni-block) majorization exists.
- **C2 — span- and spread-monotonicity are FALSE.** `meas(S7)` is **not** monotone in the
  cluster span: e.g. `meas(S7({0..6,10}))=0.190 > meas(S7({0..6,9}))=0.188`, and at `k=12,13`
  there are genuine AP-beaters (HYP-2604). Elementary "spread the top element" moves both
  increase and decrease `meas(S7)`.
- **C3 — `meas(S7)` is NOT translation-invariant** (only scale-invariant). `meas(S7(E+t)) ≠
  meas(S7(E))` (e.g. `{0..7}→0.327` but `{1..8}→0.345`, `{2..9}→0.316`), because the `e=0`
  element pins residue `0` at `x=0`. **Correction to a common slip:** THM-531's "every AP has
  the same μ_θ" (translation+scale invariance) is a statement about `μ_θ` (gap functional),
  **not** about `meas(S7)`. A "7-term AP" certifies `meas(S7)=31/210` **only** when it is a
  scale-orbit of `{0,…,6}` (i.e. `{0,d,2d,…,6d}`); a *translated* AP like `{3,5,…,15}` does
  **not** (it gives `0.188`, not `0.148`). This kills the over-eager "contains-7AP ⟹ certified".

## D. Independent adversarial confirmation of THM-534 (VERIFIED, a by-product)

Using the **correct-degree** THM-534 duals (deg-4 for k=8, deg-3 for k=9,10), an adversarial
widening confirms THM-534: `g(t)≥1[t=0]` on `{0,…,6}` (the per-E Bonferroni bound is valid);
and `L_y(consec_k)` is the **maximum** over widened primitive boxes with **zero** over-cap —
`k=8` (maxE≤16, 11432 sets), `k=9` (maxE≤17, 12869 sets), `k=10` (maxE≤17, 24310 sets). A
deg-2 functional has an *apparent* k=10 beater (`{0,2,4,5,6,7,8,10,12,14}`, `L_y=0.6345>cap`),
but that is an artifact of the **wrong degree** — with THM-534's actual deg-3 dual the beater
disappears. THM-534's k=9,10 certificate is robust.

## Net

Angle B contributes a clean **word-theoretic coordinate** for the last open rows: the AP cover
is a Sturmian partial-sum cover-time. It yields two rigorous lemmas (the `k≤6` vanishing; the
pointwise subset-domination certificate) and, more importantly, **three sharp refutations** that
show the AP-extremality of `meas(S7)` is *not* reducible to any per-block, monotonicity, or
translation-symmetry argument — it is an irreducibly aggregate (full-IE-sum) rearrangement
inequality. This sharpens the target of THM-534/THM-535/HYP-2604 and rules out three blind
alleys. **LRC(14) is NOT proved.**
