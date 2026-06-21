# LRC(14)-S3: honest end-to-end rigor status (what is PROVED vs VERIFIED vs OPEN)

**kind-pasteur-2026-06-21** (rigor-audit thread of the formalization session). An exact,
link-by-link accounting of the kps-S23 ledger, distinguishing **PROVED** (a theorem with a checked
proof), **VERIFIED** (exhaustive/finite computation, not symbolic), and **OPEN** (a genuine gap).
The purpose is to stop "L7 closed" from being mis-read as "LRC(14) proved."

## Caps (re-verified exact): `cap_k = m_k = min meas(G_P)`
`cap_8=2243/5880, cap_9=1979/4004, cap_10=55/91, cap_11=66/91, cap_12=6/7, cap_13=1`.
The values `cap_11=25/91, cap_12=2243/5880` seen in some files are **WRONG** (below the plateau `P2`).

## Link-by-link

| Link | Statement | Honest status |
|---|---|---|
| L1 | S3 ⟺ `p0(E) ≤ cap_k` | **PROVED reduction MODULO HYP-2602.** Rests on `N(E)⊆S7(E)` (PROVED, HYP-2603) ⟹ `meas S7 ≤ cap ⟹ μ_{1/7} ≥ 1-cap_k`, and the union closure THM-530/531 (PROVED-exact *given* `μ_{1/7} ≥ thr_k`). The residual gap = **HYP-2602** "consec/AP minimises `μ_{1/7}`" — VERIFIED (k=8 792 shapes; k=9,10,11 bounded spread) but OPEN symbolically for unbounded spread. AP-orbit invariance THM-531-B is genuinely PROVED. |
| L2 | `k≤7` pigeonhole | **PROVED** unconditional. |
| L3 | finite span `≤14` | **VERIFIED** exhaustive (k=8 fully, 0 violations). |
| L4 | single-far collar | **PROVED reduction (THM-547) + UNRUN finite check** `14<w≤w*` (`w*=54/90/103`, k=8/9/10) — spot-verified, not exhaustively run. |
| L5 | single-block far cluster | **PROVED** (THM-557; coherent-block quotient + BV diagonal-freeze `|p0-D_m|≤7C(m,2)/M`). The strongest link. |
| L6 | separated multi-far (`f2/f1 ≥ ρ`) | **PROVED** (iterated comb chain, `V≤7σ`, ρ-cutoffs 1.59–2.15) — separated regime only. |
| L7 | balanced multi-far (`f2/f1 ∈ (1,2.15)`) | **Analytic tail PROVED**, now SHARP: `D_{p,q}=4f(‖p‖₇,‖q‖₇)/(7pq) ≤ 44/(7pq)` (HYP-2739, residue closed form; supersedes the elementary `D≤14/p`). `|R|≤D` is the clean `0≤g_B≤1` step. **Closure VERIFIED not symbolic**: finite atlas (0 viol k=8..12), worst-base = dense even AP (verified 400+ bases), finite-f1 (THM-546 reduction), r≥3 dominated. |
| S1 | non-covering case | **PROVED** unconditional (`M ≥ 1/q ≥ 1/14`). |
| S2/S3 | covering case | reduces via THM-525 to **OPEN-Q-108** `meas(G_C) ≥ c`; **GAP G2** (transversality: `G_C` can't concentrate in `O(1/w)`-nbhds of `{a/w}`) is **OPEN**, no general argument. |

## The deepest finding: two parallel sufficient stacks

There are **two** reduction stacks, and they are distinct:
1. **Sector-cover / L7 stack** (`p0 ≤ cap`): a VALID sufficient route to LRC via `N⊆S7` + `cap_k=m_k`.
   This is the one L7 closes. Its remaining symbolic gap is **HYP-2602**.
2. **Lonely-measure stack** (OPEN-Q-108 proper, `meas(G_C) ≥ c`): what THM-525 reaches; its gap is
   **G2**.

**`L7 closed` advances stack 1, NOT OPEN-Q-108 directly** — because (kps HYP-2732, same session proves it)
the lonely side's resonance correction is *conditionally* convergent (THM-504 `|T|≥3` wall, sign-
oscillating), whereas L7's `R(p/q)` is an *absolutely* convergent measure (residue atlas + `O(1/pq)`
tail). The two stacks share the bottleneck **HYP-2602** (since `cap_k=m_k` exactly, the net-cap route
gives `μ_{1/7} ≥ 1-cap_k` at zero slack, while the extremality route gives the positive union margin
`≥ 0.32`).

## What an honest end-to-end LRC(14) proof still needs

1. **HYP-2602 symbolic**: "consec/AP minimises `μ_{1/7}(E)` for ALL integer `E`, unbounded spread."
   THM-531-B already proves every AP attains `μ_{1/7}(consec)`; the missing piece is that no non-AP
   set drops below. **Lever (this session's recommendation):** the L7 torus-line discrepancy machinery
   (the residue closed form / Koksma comb) is the right tool for the multi-block tail of HYP-2602
   (decorrelating blocks add points, shrink gaps); the union slack `≥0.28` means a CRUDE rate suffices.
2. **Run the L4 finite check** (`14<w≤w*`, ~2.6e5 sets at k=9) — upgrades L4 to fully VERIFIED.
3. **L7 finite atlas + worst-base lemma made symbolic** (currently exhaustive-verified).
4. **GAP G2 or a bounded-speed reduction** (Tao-style) for the covering case — OR rely entirely on
   stack 1 (sector route), which bypasses G2 but inherits HYP-2602.

**Bottom line:** the sector route's last *analytic* mystery (L7) is resolved and now sharp+partly
Lean-formalized; LRC(14) is **not** proved. The single highest-leverage remaining target is the
symbolic **HYP-2602**, the shared bottleneck of both sufficient routes. -> kps-S23 ledger, HYP-2602/
2603/2732/2739, THM-525/530/531/546/547/557, OPEN-Q-108, GAP G2.
