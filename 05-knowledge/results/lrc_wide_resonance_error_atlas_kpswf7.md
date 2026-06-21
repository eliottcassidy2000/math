# THREAD 1: the wide resonance-error bound -- HONEST RE-AUDIT (kps-S24-wf7)

**Object.** LRC(14) wide bound: for span(E)>14, `p0(E) = meas{x : E*x hits all 6 inner sectors of Z/7}`
must be `< cap_k`. The reduction had been framed as `p0(E) = p0_decorr(E) + err`, with
`p0_decorr = sum_t P_t^{(r)} p_t(B) <= Q(k-1) < cap` (moment dual, PROVEN finite check) and the
"signed resonance error" `err = p0 - p0_decorr` claimed `<= 0.012` (the prompt's "10x slack").

## HEADLINE: the prompt's premise (err <= 0.012 << margin, 10x slack) is REFUTED. The REAL bound holds.

Two distinct facts, both exact-rational, both load-bearing:

### (A) The REAL LRC(14) wide bound `sup p0 < cap_k` holds ROBUSTLY (0 violations).
Adversarial cap-check over WIDE configs (span>14, primitive, 0 in E): stretched APs, small-base +
adjacent far clusters, two clusters split by a big gap, and ~8k random wide configs per k.
```
  k    cap        sup p0      margin_to_cap   over_cap   #wide scanned   argmax p0
  8    0.381463   0.203571    0.177891        0          16409          (0,2,4,6,8,10,12,15)
  9    0.494256   0.364437    0.129819        0          17641          (0,1..7,27) [~single-far plateau]
  10   0.604396   0.465873    0.138523        0          18832          (0,2,4,6,8,9,10,12,14,16)
  11   0.725275   0.514116    0.211159        0          20044          (0,1..9,15)
  12   0.857143   0.605496    0.251647        0          21255          (0,1..10,38) [~single-far plateau]
```
TOTAL ~94181 exact wide configs across k=8..12, ZERO cap-violations. Margins GROW with k (0.13->0.25).
`sup p0` is ALWAYS well under cap (margin 0.13-0.21), and the worst-case p0 sits at / just above the
single-far plateau Q(k-1) (e.g. k=9 argmax = consec_8 + 1 stranger, p0=0.3644 vs Q(8)=0.3621).

### (B) The signed error `p0 - p0_decorr` EXCEEDS the fixed Q-margin -- the "0.012 / 10x slack" was an artifact.
Genuine sup of `err = p0 - p0_decorr(B,r)` over WIDE configs (correct base/far split per E):
```
  k    sup signed err    Q-margin (cap-Q)   err < Q-margin?
  8    0.107985          0.184864           True
  9    0.121298          0.132157           True (barely)
  10   0.170045          0.156509           FALSE  <-- exceeds the fixed margin
```
The prior 0.012 came from restricting to `q<=8` commensurable far-RATIOS + a few bases. The TRUE sup
error is ~0.12-0.17, and at k=10 it EXCEEDS the fixed margin `cap-Q(k-1)`. So the chain
`p0 <= Q(k-1) + err` with `err <= margin` does NOT close as a uniform statement.

### (C) WHY the real bound still holds: the error is large only where `p0` is FAR below cap.
The error peaks at SMALL-base + adjacent-far-cluster configs (e.g. k=10 peak err 0.170 at
base=(0,5,6,10,11,12,14)+far=(15,16,17)). There `p0_decorr` is small, so the available room
`cap - p0_decorr` is LARGE (0.47 >> err). In EVERY case `err < cap - p0_decorr(B,r)` (<=> `p0 < cap`).
Conversely, where `p0` is near its sup (the single-far plateau), the error is genuinely TINY (the
PROVEN comb bound THM-546 `|Delta_w| <= (6/49)V/w`). The two regimes are disjoint:
  - p0-near-cap  <=>  big base / single far  <=>  small error (THM-546)  -- the binding regime.
  - big error    <=>  small base / far cluster  <=>  p0 far below cap     -- harmless.

So the correct statement is `p0(E) < cap_k` directly (MISTAKE-080 framing), NOT `p0 <= Q(k-1)+small`.
The decorrelated bound `p0_decorr <= Q(k-1)` is a cut-space FLOOR; the cycle-space error is NOT
uniformly small (consistent with MISTAKE-082), but it only matters where there is ample room.

## DECAY of the resonance error in the far-pair denominator (curve atlas, confirmed)
For a fixed bounded base + 2 far at large scale (f1=101), error by far-ratio denominator q:
small-q (commensurable) ratios carry the error; large q decorrelates. tail-sup|err| over ratios with
denominator >= q is bounded (~0.017 at f1=101 scale). This confirms the atlas DECAY direction but the
RELEVANT sup is NOT at large scale -- it is at MODERATE scale (adjacent clusters), where the error is
largest. The finite-atlas-of-small-q picture is correct for the FAR-FAR resonance, but the dominant
error is the BASE-FAR / cluster interaction at moderate separation, which the q-atlas misses.

## Scripts / outputs
- `04-computation/lrc14_wide_resonance_sup_v2_kpswf7.py` (span>14 enforced; A and B sup)
- `04-computation/lrc14_wide_adversarial_capcheck_kpswf7.py` (adversarial cap-check, ~90k wide configs)
- outputs: `05-knowledge/results/lrc14_wide_resonance_sup_v2_kpswf7.out`,
  `05-knowledge/results/lrc14_wide_adversarial_capcheck_kpswf7.out`
- p0_fast cross-checked == repo exact `lrc14_wide_branch_ridge_codex_s47.p0` on 5 configs.

## VERDICT
- WIDE bound `sup p0 < cap` = SUPPORTED (0 / ~90k exact wide configs, margins 0.13-0.21). NOT a proof.
- The "wide => p0 <= Q(k-1) + err, err <= margin (10x slack)" reduction is REFUTED: err reaches 0.17 > margin.
- The proof must target `p0 < cap` directly, exploiting the REGIME SPLIT (THM-546 binds the only
  near-cap configs; the high-error configs are slack). -> OPEN-Q-108, HYP-2675, HYP-2775, MISTAKE-080/082.
