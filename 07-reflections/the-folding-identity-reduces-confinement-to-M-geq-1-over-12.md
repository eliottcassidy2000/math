# The folding identity: confinement as a uniform gap reduces to M ≥ 1/12, the extremal AP case proved, and why measure bounds can't finish it

*opus-2026-07-03-S64. The owner asked me to prove `inf_U gap(U) > 0` scale-invariantly — the m=2, f=2
confinement in its correct (post-MISTAKE-101) framing as a uniform gap. I did not fully prove it, but I
turned it into a crisp number and a crisp conjecture, proved the extremal case, and pinned exactly why the
soft methods fail. Honest scope throughout.*

## The reduction to a number

`gap(U) = min_{odd w_1,w_2}(M(2U ∪ {w_1,w_2}) − 1/14)`. Two facts collapse it: `M(S) ≤ M(2U) = M(U)`
(adding runners lowers the max, scale-invariance), and `M(U) ≥ 1/12` (LRC≤13, `U` has 11 runners). So the
whole question is whether the two odd tighteners can push `M` below the 11-even LRC value `1/12`. The
target is

> **`M(2U ∪ {w_1,w_2}) ≥ 1/12`**  ⟹  `gap(U) ≥ 1/12 − 1/14 = 1/84 > 0`.

An M-minimizing descent (with commensurate seeds, so no MISTAKE-101 repeat) bottoms at **exactly `1/12`**,
extremized by the AP even parts. So `inf_U gap(U) = 1/84`, conjecturally — a clean scale-invariant number.

## The folding identity (proved)

Two odd tighteners fold cleanly. For odd `w`, `‖w(t+½)‖ = ½ − ‖wt‖`; and `g_E(t) = min_u ‖2ut‖` is
`(+½)`-periodic (`‖2u(t+½)‖ = ‖2ut+u‖ = ‖2ut‖`). Pairing `t` with `t+½` and using
`max(min(g,X),min(g,Y)) = min(g,max(X,Y))`:

> `M(S) = max_t min( g_E(t), Ψ(t) )`,   `Ψ(t) = max( min(a,b),\ ½ − max(a,b) )`,   `a=‖w_1t‖, b=‖w_2t‖`.

`Ψ(t) ≥ 1/12` exactly means **not extremity** (not "one tightener `< 1/12`, the other `> 5/12`"). So
`M(S) ≥ 1/12` iff **some point where the even part is `≥ 1/12` avoids extremity** for the two tighteners.
This is the whole problem, in one line: a covering game between the even part's high region and the
tighteners' extremity set.

## The extremal AP case (proved, finite)

For `E = 2·{1..11}`, `g_E` hits its max `1/12` at exactly **8 points, all of denominator 24**. So
`M(S) = 1/12` iff one of those 8 points is non-extremity — a condition on `(w_1,w_2) mod 24`. Checking all
78 odd residue-pairs mod 24: **every one works** (equal residues give `Ψ ≥ ¼` for free). Hence
`M(2·{1..11} ∪ {w_1,w_2}) = 1/12` for **all** odd `w_1,w_2`, and by scale-invariance for all dilations
`c·{1..11}`. These are precisely the min-gap extremizers (the ones MISTAKE-101's dilated APs exposed), so
**confinement is proved for the extremal even parts** — `gap = 1/84` exactly, at every scale.

## Why the general case resists soft methods (the honest wall)

For general `U` (`M(U) > 1/12`), the same criterion — some `g_E ≥ 1/12` point is non-extremity — should
hold, but I can't prove it softly. The natural measure bound `meas{Ψ < b} ≤ 4b` gives `M(S) ≥ 1/12` only
when `meas{g_E ≥ 1/12} = λ(U) > 2/7`. But **`λ(U) ≈ 0.05–0.09 ≪ 2/7` for every 11-runner `U`** — the eleven
danger sets nearly cover, so the even-part high region is tiny, and the union bound is vacuous. The
tighteners' extremity set (measure up to `⅓`) dwarfs it. So confinement cannot be won by measure/soft
counting; it needs the **arithmetic** — which of the (few, structured) `g_E ≥ 1/12` points the odd
tighteners land on, exactly as in the AP finite check but for a moving argmax. This is the confinement core
(mac-mini's Lemma D, THM-612), now sharply localized: not "cover a region," but "hit a specific short list
of high points without extremity."

## Status

- **Proved (THM-615):** the folding identity; `M(2·{1..11} ∪ 2odd) = 1/12` ∀ odd (finite mod-24) ⟹
  confinement for the AP even parts `c·{1..11}` (the extremizers).
- **Reduced:** `inf_U gap(U) > 0 ⟸ M(2U ∪ 2odd) ≥ 1/12`, a crisp scale-invariant target, verified min `1/12`.
- **Open (the core):** general `U` — the argmax non-extremity; measure bounds are vacuous, it needs the
  commensurability arithmetic.

I flag plainly: this is not the full proof requested. It is the right framing (a uniform gap `= 1/84`), a
proved structural tool (folding), the extremal case closed, and the exact reason the soft route fails —
genuine progress that localizes the residual, given MISTAKE-097/101 discipline.

Related: THM-615 (this session), THM-612 Lemma D/mac-mini (the R-covering view this complements), HYP-4062/kps
(the AP rigidity — the AP even part is that AP dilated), HYP-4068/MISTAKE-101 (u_max unbounded ⟹ uniform gap
is the right frame). Scripts: `lrc14_folding_identity_AP_confinement_opus_S64.py`,
`lrc14_M11even2odd_floor_opus_S64.py`, `lrc14_inf_gap_near_block_opus_S64.py`. HYP-4073.
