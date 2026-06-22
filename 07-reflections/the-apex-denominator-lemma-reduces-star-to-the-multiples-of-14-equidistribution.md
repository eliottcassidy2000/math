# The apex-denominator lemma: tight implies an apex shell, not shell collapse

*kind-pasteur-2026-06-22-S31aa. Attacking mac-mini's linchpin (★) [M(S)=1/14 ⟹ the optimum is at a
denom-14 apex-blocked point]. A short rigorous lemma forces the optimum denominator to be a multiple of
14.  S120 corrected the strongest reading: the local arithmetic proves an apex shell `D=14h`, not
primitive shell collapse `h=1`; THM-571 later closes the `>=7` multiples-of-14 branch, so the live
strictness target is the `<=6` scale-separated / finite-core branch.

## The lemma (PROVED — elementary, corrected)
**Lemma. If `M(S)=1/14`, then at the optimum `t* = a/D` (lowest terms), `14 | D`; the two binding
runners `{v_a,v_b}` (at `±1/14`) satisfy `D | (v_a+v_b)`.**

*Proof.* At a tight optimum the origin sits in an empty arc of length exactly `1/7`, with binding runners
at **both** endpoints `±1/14` (if only one endpoint were occupied the arc would be wider and `M>1/14`).
For the `+1/14` binder: `v_a·(a/D) = m ± 1/14`, so `14(v_a a − mD) = ±D`, giving
`D(1+14m) = 14 v_a a`; since `gcd(1+14m,14)=1`, this gives **`14 | D`.** From the two binders
`(v_a+v_b)t* ∈ ℤ`, so `D | (v_a+v_b)a`, hence `D | (v_a+v_b)`. ∎

Verified (`lrc_apex_denominator_lemma_kps.py`): AP `D=14`, binders `{1,13}` sum 14; GW `D=14`, binders
`{5,9}` sum 14; `2·AP` `D=28` sum 28; `3·GW` `D=42` sum 42.  These examples are compatible with
`D = 14·gcd(S)`, but that equality is not proved by the local binder calculation alone.

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

2. **14-covering strictness is the residual.** Write `S = M14 ⊔ R`, `M14` = the multiples of 14, `R` =
   the rest (14-free). `R` has `≤12` speeds, so by *proven* LRC(≤13), `M(R) ≥ 1/13 > 1/14` — there is an
   **open interval `I` around `R`'s witness on which `min_R > 1/14`** (margin). On `I` the multiples of 14
   sweep as combs.  S31v is the relevant union-bound / component-budget mechanism for the `|M14|<=6`
   branch, but still needs the finite-core compression or uniform interval budget that makes the margin
   survive on every shell.  THM-571 closes the complementary `|M14|>=7` branch by descending to the
   `7`-phase and pigeonholing the lifts.

Putting them together:
> **tight ⟹ apex shell `D=14h`; 14-free tight rows already have a denominator-14 survivor; 14-covering
> rows must be made strict.**  Shell collapse `h=1` or covering strictness is the remaining theorem.

## What is proved vs. the residual
- **PROVED:** the local apex-shell arithmetic (`14|D`, `D|(v_a+v_b)`), formalized in
  `LRCApexShell.lean`; the 14-free denominator-14 survivor from THM-523/THM-569; and the
  `|M14|>=7` apex-majority branch by THM-571.
- **RESIDUAL:** the `|M14|<=6` scale-separated / finite-core strictness package, or an equivalent
  shell-collapse/state-lift theorem that rules out primitive tight rows on shells `D=14h` with `h>1`.

So (★) is no longer a raw residue-atlas problem.  The lemma did the local structural work (forcing the
apex shell); the remaining work is to cash out the shell height or prove that the few-multiple covering
branch cannot cover the smaller-core margin.

## Net
mac-mini's (★) reduces, via the corrected apex-shell lemma, to: **14-free tight sets have a
denominator-14 survivor (done, THM-523/THM-569), `|M14|>=7` covering rows are safe by THM-571, and the
live branch is `|M14|<=6` covering strictness / shell collapse.** Combined with Move A (peel) and Move B
(apex floor), this assembles the THM-079-template proof down to that one residual.

→ (★)/mac-mini's HYP-2906, THM-523 (q-witness), S31v (comb-teeth union bound + `M(R)≥1/13`), THM-560
(tight locus / dilation), THM-079 (the H=21 template), [[lrc14-thread]].
