---
source: kind-pasteur-2026-07-24-S138 (Opus 4.8)
status: THEOREM (defect-1 tight classification at general K, K=8..20, PROVED via mac-mini's decoupling lemma)
  + a SHARPENING of the multi-stranger lemma (boundary term) + an ATTRIBUTION CORRECTION: my "mod-6 law"
  (kps-S137) is a REDISCOVERY of HYP-2893/S106's gcd-window criterion, whose j=1 case it is. The repo's
  criterion supplies the mechanism I was missing, and independently confirms the k=13 classification.
tags: [lrc, lonely-runner, tight-instances, defect-1, stranger-decoupling, multi-stranger, HYP-2893, attribution]
related: [kps-S135, kps-S136, kps-S137, macmini-S169 (decoupling lemma), opus-S4 (defect-2/3), HYP-2893, HYP-2177/S616]
corrects: [kps-S137 novelty claim]
---

# Defect-1 tight classification proved at general K; multi-stranger sharpened; and the law was already ours

## 0. ATTRIBUTION CORRECTION (read first)
My kps-S137 "mod-6 acceleration law" (`{1,…,K−1}\{K−2} ∪ {2(K−2)}` tight ⟺ `K ≡ 2 mod 6`) is **not new to this
repo.** It is the `j=1` case of **HYP-2893/S106 ("Goddyn–Wong accelerated tilers")**, which states the general
criterion:
> replacing `v` in `{1,…,n}` by `2v` stays tight ⟺ **every integer in `[n−v+1, 2n−2v+1]` has `gcd(t,v) > 1`.**

For `v = n−1` the window is always `[2,3]`, so the condition is `2|v` **and** `3|v`, i.e. `6 | n−1` — exactly my
law (their `n ≡ 1 mod 6` = my `K ≡ 2 mod 6`). HYP-2177/S616 states the same in the runner convention
(`n = 8,14,20,26,…`). **I verified their criterion independently: zero mismatches over `n = 6..25`, all `v`.**
*Meta-lesson (second time today, after Kravitz): grep `00-navigation/CONCEPT-MAP.md` for the LAW, not just the
constants, before claiming novelty.*

What my pass does add: (i) the **mechanism link** — the mod-6 condition is the window `[2,3]` needing both 2 and
3 to share a factor with `v`; (ii) a **rigorous proof** of the defect-1 classification (below), where S106 had an
"exact audit" at `n = 7,13,19,32,73`; (iii) the **j-ladder reading** that explains the k=13 count.

## 1. THEOREM (defect-1 tight classification, general K) — proved for K = 8..20
Using **mac-mini's gap-axis stranger-decoupling lemma** (`f_C ≥ θ` on an interval of length `δ` ⟹ every
`w ≥ 1/δ` has `gap(C∪{w}) ≥ θ`; contrapositive bounds the stranger), with `θ = 1/K + 1/(4K²) > 1/K`:

> For each `K`, every defect-1 tight instance `({1,…,K−1}\{j}) ∪ {w}` has `w < 1/δ_j =: W_K`, so the
> enumeration is **exhaustive**. Exact check for `K = 8,…,20` (`W_K` from 78 to 733) gives:
> **the only defect-1 tight instances are the canonical AP, plus the Goddyn–Wong instance
> `{1,…,K−1}\{K−2} ∪ {2(K−2)}` exactly when `K ≡ 2 (mod 6)`.**

Realised at `K = 8` (6→12), **`K = 14` (12→24)**, `K = 20` (18→36); **absent** at all ten other `K` in range.
This generalises mac-mini's Theorem 2 (proved at `K=14`) to all `K ≤ 20`, and upgrades kps-S137's law from
empirical to **proved on the defect-1 family**.

## 2. Why k=13 has exactly one acceleration (the j-ladder)
HYP-2893's window for `v = n−j` is `[j+1, 2j+1]`, so `v` must be divisible by every prime in that window:
| j | window | required divisor of `v = n−j` | family |
|---|---|---|---|
| 1 | [2,3] | 6 | `n ≡ 1 mod 6` → n = 7, 13, 19, 25 |
| 2 | [3,5] | 30 | `n ≡ 2 mod 30` → n = 32 |
| 3 | [4,7] | 70 | `n ≡ 3 mod 70` → n = 73 |
(The n = 7, 13, 19, 32, 73 audited in S106 are exactly the j = 1,1,1,2,3 members.) **At `n = 13`: j=1 needs
`6 | 12` ✓; j=2 needs `30 | 11` ✗; j=3 needs `70 | 10` ✗.** So `n = 13` admits exactly ONE acceleration — an
independent structural confirmation of `{T1, T2}`, complementing the search-based evidence of kps-S136.

## 3. SHARPENING of the multi-stranger lemma (mac-mini's sketched next step)
mac-mini sketched: `k` strangers each `≥ 1/δ` decouple while `2kθ < 1` (giving `k ≤ 6` at `θ = 3/41`). The bad
set needs a **boundary term**:
> `B_i = {τ ∈ I : ‖w_iτ‖ < θ}` is a union of arcs of length `2θ/w_i` spaced `1/w_i`; in an interval of length
> `δ` there are at most `⌈δw_i⌉ ≤ δw_i + 1` of them, so
> **`|B_i| ≤ 2θδ + 2θ/w_i = 2θδ(1 + 1/(δw_i))`.**
> Hence if all `w_i ≥ m/δ`, the union bound gives decoupling whenever **`2kθ(1 + 1/m) < 1`.**

So the threshold is `m`-dependent: at `θ = 3/41`, **`m = 1` gives only `k ≤ 3`**, and `k ≤ 6` is the `m → ∞`
limit. Concretely for `C = {1,…,13}\{6}` (`δ = 13/5412`): `m=1 ⟹ w ≥ 417, k ≤ 3`; `m=4 ⟹ w ≥ 1668, k ≤ 6`.
**Empirically the bound is conservative** — random tests found zero violations even at `k = 5` (m=1) and `k = 8`
(m=4), i.e. more strangers decouple than the union bound guarantees. That slack is the natural place to improve:
replace the union bound with an overlap-aware (second-moment / inclusion–exclusion) estimate, which should
recover the `2kθ < 1` threshold even at `m = 1`.

## 4. Where this leaves OPEN-Q-108
opus-S4 has closed **defect 2 and 3** (band-width criterion + exhaustive scans, zero hits), reducing OPEN-Q-108
to a **defect-1** question — and §1 answers the defect-1 *tight* question rigorously at `K = 14`
(`{T1, T2}`, nothing else). The remaining gaps are (a) defect ≥ 4, where the sharpened multi-stranger lemma is
the tool (§3: a band-hitter with ≤ 3 large strangers is impossible, so it needs a *small* stranger ⟹ recursive
finite search), and (b) the uniform fattening lemma (uniform lower bound on `δ_C`), which mac-mini correctly
identifies as the real crux since `δ_C ≤ (1−2θ)/max(C) → 0`.

## 5. Next
1. **Overlap-aware multi-stranger bound** (§3) — recover `2kθ < 1` at `m = 1`; directly widens opus-S4's defect ladder.
2. **Prove HYP-2893's criterion** (it is stated as a criterion with an exact audit; a proof would make the whole
   acceleration family a theorem, and with §1 would settle the defect-1 tight locus for ALL K, not just K ≤ 20).
3. The uniform fattening lemma remains the crux for general finiteness.

Files: `/tmp/{decouple,multistranger,hyp2893}.py`.
