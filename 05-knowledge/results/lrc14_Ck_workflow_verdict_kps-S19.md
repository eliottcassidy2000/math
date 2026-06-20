# LRC(14) "one open constant" workflow — final adversarial verdict (kps-S19)

4-angle workflow (per-sector-koksma, per-scale-cluster, resonance-direct, finite-check) + independent
verifiers, on the target `C(k) := sup_{E',w} w|Δ_w| ≤ c·k`. Two verifiers died on API 529; the synthesis
re-derived their checks independently. EXACT-Fraction throughout.

## Verdict on C(k)
**`C(k) = +∞` at FIXED k — the literal target is DISPROVED (PROVED refutation, witnessed).**
- Witness 1 (per-sector): `E'={0,1,2,3}∪{M,M+1,M+2,M+3}, w=22M` (primitive, w=max): `w|Δ_w| = 1.18,2.05,3.46,7.94,12.42,25.44` at M=10,20,40,80,160,320 — `~0.08·M`, k=9 fixed.
- Witness 2 (resonance): #clusters=2 FIXED, scale M grown: `sup_w w|Δ_w| = 1.73,3.44,6.22,12.0` at M=30,60,120,240 — linear in M.
- Supersedes HYP-2645 (1.95), HYP-2655 (3.91), HYP-2653b/c — all artifacts of a CAPPED w-range.

## REFUTED: the per-scale-cluster "bounded C(k)" angle
`lrc14_ck_per-scale-cluster_kps-S18-wf.py` claims `C(k) ≤ 2.72·(k−2)` (single cluster saturates at K₁≈1.11
"independent of cluster size L"). **REFUTED.** Flaw: it conflates a cluster's internal *length L* with the
*scale separation M between clusters*. Witness 2 (2 clusters, growing M) directly violates it. The
decomposition steps (telescope, per-cluster attribution) are exact and fine; only the per-cluster O(1) bound
is false. Its peel thresholds (T_8=88, T_9=144, T_10=139) are VOID. The `finite-check-to-threshold` angle
imported that false constant, so its thresholds are also void — but its finite-check engine and direct-p0
facts survive.

## SURVIVING rigorous facts (the route closes ONLY via the direct/joint route)
| Link | Status |
|---|---|
| LRC(14)-S3 ⟺ p0(E) ≤ cap_k (primitive k-sets) | PROVED upstream |
| engine identity w·Δ_w = Σ_{|miss|=1}[G0(wb−s/7)−G0(wa−s/7)] (HYP-2653) | VERIFIED exact |
| per-sector telescope (THM-PSK-1) + σ-bound w|Δ_w| ≤ (6/7)σ(E') (THM-PSK-2) | PROVED (but σ-dependent, LOOSE at resonance — cannot close alone) |
| C(k) ≤ c·k | **DISPROVED** |
| p0(consec_k) < cap_k: 0.327/0.416/0.504 < 0.381/0.494/0.604 | VERIFIED exact |
| finite check span ≤ 14, k=8..12, consec argmax, 0 violations | **PROVED (kps-S19; stronger than the workflow's span-12)** |
| **wide ⟹ small p0 (span>B ⟹ p0 ≤ cap_k) = HYP-2675** | **CONJECTURE (sampled, margins ≥0.21; NOT proved)** |

## The SOLE residual = HYP-2675 (codex-s47), now triply-converged
`span(E)>14 ⟹ p0(E)=meas(S7(E)) ≤ cap_k` for primitive k-sets, k=8..12. Empirical wide max p0 (sampled):
k=8 ≤ 0.175, k=9 ≤ 0.19, k=10 ≤ 0.30 (margins ≥ 0.21/0.27/0.30). Span-monotone decay (k=8):
0.327(consec) → 0.175 → 0.126 → 0.074 for span bands [8,12]/[13,16]/[17,28]/[40,60].

## Proof direction (synthesis's key insight)
Prove HYP-2675 as a **coverage lemma**, bypassing Δ_w/C(k): the SAME multi-scale residue collapse that
*destroys* the Δ_w cancellation (making w|Δ_w| large) is exactly what *spreads sector coverage thin*
(making p0 small). A wide set's clusters decorrelate (Weyl/cross-scale independence), so p0 ≤ ∏_clusters
(bounded per-cluster coverage) < cap. Turn the span-monotone decay into an Erdős–Turán coverage bound on the
widest cluster ⟹ explicit B with p0 ≤ f(B) < cap_k. → HYP-2675, HYP-2653d, THM-PSK-1/2/3/4, OPEN-Q-108.
