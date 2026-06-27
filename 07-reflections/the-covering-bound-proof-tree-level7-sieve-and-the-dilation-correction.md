# The LRC(14) covering-bound proof tree: the level-7 sieve, and the dilation correction

*kind-pasteur-2026-06-27-S31af. The owner asked to integrate the incoming work (mac-mini-S59's
covering-bound REDIRECT, codex-S122 THM-571, kps Node-3/S31v) and past work into one comprehensive
picture of what remains, and let that synthesis suggest an improved argument. It did two things: a clean
**sharpening** (the level-7 lift sieve, THM-573) and a clean **correction** (covering sets are NOT
bounded away from 1/14 — the redirect's "strictly weaker" framing is wrong, though its logic stands).*

## 1. The reduction (mac-mini-S59) — valid, and confirmed
Split every 13-set `S` by the apex residue `mod 14`:
- **14-free** (no `v ≡ 0 mod 14`): `t = 1/14` gives `‖v/14‖ ≥ 1/14` for every speed. `M ≥ 1/14`. **Trivial** (THM-523). [Verified: 0 failures on thousands of 14-free sets.]
- **covering** (some `v ≡ 0 mod 14`): the only case needing work.

So **LRC(14) ⇔ [every covering 13-set has `M ≥ 1/14`]** — the covering bound. This reduction is correct
and the 14-free half is genuinely free. mac-mini's logical point stands: the **bound** `M ≥ 1/14` is what
LRC(14) needs, and that does not *logically* require the tight-locus census.

## 2. The sharpening (THM-573): the level-7 lift sieve
THM-571 closed `|M14| ≥ 7` (≥7 multiples of **14**) with a Case-1/Case-2 split. Its Case 2 is already a
level-7 sieve, fired only at `|H| ≥ 9`. The clean fact: **the level-7 sieve alone needs only `|H| ≥ 7`**,
where `H` = the multiples of **7**. One argument, no split, strictly larger domain (`|H| ≥ |M14|`):

> **THM-573.** If `≥ 7` of the 13 speeds are divisible by 7, then `M(S) > 1/14`.

Proof in one breath: `H = 7P`, take a `P`-safe phase `v*` (`‖pv*‖ ≥ 1/13` by LRC(≤13)); the seven lifts
`t_j=(v*+j)/7` keep all of `H` safe (`‖7p·t_j‖=‖pv*‖`); each 7-coprime speed sits among 7 points spaced
`1/7` so the length-`1/7` forbidden arc catches `≤ 1` of its lifts; with `13-|H| ≤ 6 < 7` coprime speeds,
a lift survives. (Verified: 1500/1500 fired; the `≤1`-forbidden crux confirmed over 20000 trials.)

**This relocates the residual** from "`≤ 6` multiples of 14" to "**`≤ 6` multiples of 7**" — strictly
smaller. The descent is single-step: `7` is forced (forbidden arc `1/7` = lift spacing `1/d` ⟹ `d=7`);
level 2 is useless; no other common prime is guaranteed. Below `|H| = 7` the arithmetic descent stops.

## 3. The correction: covering sets are NOT bounded away from 1/14
The tempting next step was a **margin**: "covering ⟹ `M ≥ 1/13`", which would make the covering bound
quantitatively easier than the census. **It is false.** Two refuters (exact):
- `2·{1,…,13} = {2,4,…,26}` is **covering** (it contains `14 = 2·7`) and **tight**: `M = 1/14` by
  dilation invariance `M(dS)=M(S)`. So a *dilated AP* is a covering set sitting exactly at `1/14`.
- The **primitive** `{1,…,12,182}` has `M = 14/183 ≈ 0.0765 < 1/13`. Here `182 = 14·13` **aliases to 0
  mod 13** (`14 ≡ 1 mod 13`, so `14·13 ≡ 13 ≡ 0`), landing on the optimum `t = k/13` of the LRC(13)-tight
  core `{1,…,12}` and dragging `M` below `1/13`. (S3 confirms: adding `14j` to `{1..12}` stays `≥ 1/13`
  for **all** `j` except `j ≡ 0 mod 13`.)

**Consequence for the redirect.** mac-mini's claim that "the census `{AP,GW}` is non-covering, hence off
the critical path, and the bound is *strictly weaker*" is half right and half wrong:
- ✔ the **reduction** is valid and the census is not *logically* needed for the bound;
- ✘ the covering bound is **tight** (achieved by dilated APs), so it is **not** "strictly weaker" with a
  margin. The hard tight cases (dilations of `AP`/`GW`) live **inside** the covering family, all with
  `|H|` small — e.g. `2·AP` has `|H| = 1`. They sit squarely in THM-573's residual.

The honest distillation: **WLOG primitive** (dilation invariance). Then 14-free primitive is trivial, and
the open content is *primitive covering with `≤ 6` multiples of 7* — which still includes near-tight rows
like `{1..12,182}`. No free margin; the difficulty was relocated, not dissolved.

## 4. The integrated proof tree (what is proved, what remains)
```
LRC(14)  ==WLOG primitive==>
 ├─ 14-free                         t = 1/14                      DONE (THM-523, trivial)
 └─ covering (|H| ≥ 1):
     ├─ |H| ≥ 7  (≥7 mult of 7)     level-7 lift sieve            DONE (THM-573 ⊇ THM-571)
     └─ |H| ≤ 6  (≤6 mult of 7)  ── RESIDUAL ── :
         ├─ a large speed present   scale peel / equidistribution  r≤6 RIGOROUS (S31v),
         │                          (Node-3, danger combs)         r≥7 needs resonance-defect bound
         └─ bounded core            finite check                   astronomical V*  (OPEN)
```
The two prior partitions — codex's *by apex count* and kps's *by scale* — are **orthogonal axes**; they
compose. THM-573 collapses the apex axis to `|H| ≤ 6`; the scale axis (Node-3) then attacks the rest, with
the genuine open piece the bounded core of `≥ 7` comparable speeds coprime to 7 (the `r ≥ 7`
second-moment / resonance-defect estimate, and the unreachable finite `V*`). This is the same hard core
the project has circled — now **precisely located**: `≤ 6` multiples of 7 over a bounded coprime core.

## 5. Net
- **Gained:** THM-573 (clean, verified, subsumes THM-570+571 and enlarges the closed domain); the residual
  is now "`≤ 6` multiples of 7", strictly smaller and stated by a single invariant.
- **Corrected:** the covering bound is *tight* (dilated AP at `1/14`); "covering ⟹ `M ≥ 1/13`" is
  refuted; the apparent `1/13` margin in the first scan was a search artifact (it omitted dilations).
- **Unchanged (honest):** LRC(14) is **not** proved. The bounded core (`≥ 7` speeds coprime to 7, no
  large speed to peel) remains open — the effective Erdős–Turán / finite-`V*` wall. The reduction and the
  sieve sharpen the bookkeeping and the perimeter; they do not breach that wall.

→ THM-573, THM-571, THM-568, THM-523, HYP-3084 (margin refuted), S31v (Node-3 union bound),
mac-mini-S59 (the redirect), [[lrc14-thread]].
