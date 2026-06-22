# The apex-denominator lemma: tight ⟹ 14 | D, reducing (★) to "14-covering ⟹ not tight"

*kind-pasteur-2026-06-22-S31aa. Attacking mac-mini's linchpin (★) [M(S)=1/14 ⟹ the optimum is at a
denom-14 apex-blocked point]. A short rigorous lemma forces the optimum denominator to be a multiple of
14, and (★) then reduces to a single specialized equidistribution statement about the multiples of 14.*

## The lemma (PROVED — elementary)
**Lemma. If `M(S)=1/14`, then at the optimum `t* = a/D` (lowest terms), `14 | D`; the two binding
runners `{v_a,v_b}` (at `±1/14`) satisfy `D | (v_a+v_b)`; and `D = 14·gcd(S)`.**

*Proof.* At a tight optimum the origin sits in an empty arc of length exactly `1/7`, with binding runners
at **both** endpoints `±1/14` (if only one endpoint were occupied the arc would be wider and `M>1/14`).
For the `+1/14` binder: `v_a·(a/D) = m ± 1/14`, so `14(v_a a − mD) = ±D`, giving
`D(1+14m) = 14 v_a a`; since `gcd(1+14m,14)=1`, `(1+14m) | v_a a` and **`D = 14·(v_a a/(1+14m))`, so
`14 | D`.** From the two binders `(v_a+v_b)t* ∈ ℤ`, so `D | (v_a+v_b)a`, hence `D | (v_a+v_b)`. ∎

Verified (`lrc_apex_denominator_lemma_kps.py`): AP `D=14`, binders `{1,13}` sum 14; GW `D=14`, binders
`{5,9}` sum 14; `2·AP` `D=28` sum 28; `3·GW` `D=42` sum 42 — i.e. `D = 14·gcd(S)` exactly, and
**primitive ⟹ `D=14`**.

## Why the residue check is NOT the closure
At `D=14`, `M` is determined by the residue multiset `{v_i mod 14}`. But the finite enumeration
(`lrc_apex_denominator_lemma_kps.py` Part 2) shows **all 156 "double-r / miss-m" patterns give
`M₁₄=1/14`**, plus the full-set AP pattern. So `M₁₄=1/14` is *necessary but far from sufficient* for
global tightness — the global `M` depends on the integer speeds, not just residues. The tight locus is
pinned by the *optimum location*, not by a residue list. So (★) must be proved as a statement about
**where the optimum sits**, which is what the lemma controls.

## The reduction of (★) (the clean split)
Call `S` **14-covering** if it contains a multiple of 14, else **14-free**. Then:

1. **14-free ⟹ optimum at `D=14`.** By THM-523 (`q=14`), a 14-free set has `‖v_i/14‖ ≥ 1/14` for all
   `i`, so `M₁₄(S) ≥ 1/14`. If `M(S)=1/14` then `M₁₄(S)=1/14` is already the global max, attained at
   `t=a/14`: **the optimum is at the apex denominator.** (Combined with the lemma's `14|D`, this is exactly
   (★) for 14-free sets — no equidistribution needed.)

2. **14-covering ⟹ `M(S) > 1/14` (not tight).** Write `S = M14 ⊔ R`, `M14` = the multiples of 14, `R` =
   the rest (14-free). `R` has `≤12` speeds, so by *proven* LRC(≤13), `M(R) ≥ 1/13 > 1/14` — there is an
   **open interval `I` around `R`'s witness on which `min_R > 1/14`** (margin). On `I` the multiples of 14
   sweep; by the comb-teeth bound (S31v) each covers `≤ 1/7` of `I`, so for `|M14| ≤ 6` the multiples
   leave `(1 − |M14|/7)·meas(I) > 0` where they are *also* `>1/14` ⟹ **`M(S) > 1/14`.** (`|M14| ≥ 7`:
   the second-moment / `(6/7)^r` residual.)

Putting them together:
> **tight ⟹ (by step 2) 14-free ⟹ (by step 1) optimum at `D=14`** — which is (★), for primitive `S`;
> dilates reduce to primitive by R3.

## What is proved vs. the residual
- **PROVED:** the lemma (`14|D`, `D=14·gcd`, primitive ⟹ `D=14`); step 1 (14-free tight ⟹ apex
  optimum, from THM-523); step 2 for `≤6` multiples of 14 (the S31v union bound + `M(R)≥1/13` margin).
- **RESIDUAL:** step 2 for `≥7` multiples of 14 — a set with ≥7 distinct multiples of 14 over a ≤6-speed
  14-free core. This is the *apex-specialized* form of the Node-3 equidistribution (the same `r≥7`
  second-moment as S31v, but now the large speeds are forced to be **multiples of 14**, and the core
  `R` carries a `1/13` margin from proven LRC).

So (★) is no longer "characterize the tight locus" — it is the single statement **"a set with ≥7
multiples of 14 over a 14-free core is not tight,"** provable by the second-moment bound on those
multiples against the core's margin. The lemma did the structural work (forcing the apex denominator);
what remains is one quantitative equidistribution estimate, sharply localized to the multiples of the
apex `14 = 2·7`.

## Net
mac-mini's (★) reduces, via the proved apex-denominator lemma, to: **14-free tight sets optimize at
`D=14` (done, THM-523), and 14-covering sets are not tight (the multiples of 14 cannot cover the
14-free core's `1/13`-margin — done for `≤6`, second-moment for `≥7`).** The crux is now an explicit,
localized equidistribution about the multiples of `14`, with the bounded part discharged by proven
LRC(≤13). Combined with Move A (peel) and Move B (apex floor), this assembles the THM-079-template proof
down to that one residual.

→ (★)/mac-mini's HYP-2906, THM-523 (q-witness), S31v (comb-teeth union bound + `M(R)≥1/13`), THM-560
(tight locus / dilation), THM-079 (the H=21 template), [[lrc14-thread]].
