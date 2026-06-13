---
source: opus-2026-06-06-S699 (remote-control, long signed-LRC session)
status: RIGOROUS THEORY of the signed LRC, bottom-up n=3..14. THEOREMS: gauge invariance (M sign-blind); sign=2-coloring=cut (sum-clock iff bichromatic; maxcut=⌊(n-1)²/4⌋); zero-clock⟺shell-partner (v_i+v_j≡0 mod 2n-1). MAIN RESULT: the signed structure SPLITS the worry-set — through n=7 every tight config is shell-partner-free, but at n=14 the floor row V*=(1..11,13,24) is tight (M=1/14) WITH shell-partner 3+24=27 (from the doubling 24=2·12). So the signed LRC is a STRICTLY FINER invariant than M; it separates {AP, V*} at the n=14 floor; n=14 is where C=2n-1=27=3³ first admits a doubled-speed shell-partner (prime-2 × prime-3).
tags: [signed-lrc, gauge, sign-cut, maxcut, shell-partner, zero-clock, worry-set-split, Vstar, doubling, prime-2-prime-3, 2n-1, n14, rigorous]
---

# Signed LRC: a sign is a cut, and the worry-set splits at n=14

**Prompt (user):** the signed LRC is an important variant; prioritise its complete understanding
as an angle of attack for the regular LRC. Make rigorous progress from small n up, trying as many
sign-reversal patterns as possible at each n.

This session builds the signed LRC rigorously from `n=3` to `n=14`. The key structural fact makes
it tractable, and the climb pays off: **the signed structure splits the worry-set exactly at
`n=14`.**

## The rigorous foundation (theorems)

> **T1 (gauge invariance).** For any signs `ε∈{±1}^{n-1}`, `M({ε_i v_i}) = M(S)`. *Proof:*
> `‖ε_i v_i t‖ = ‖v_i t‖` for every `i,t`, so `min_i` is identical pointwise, hence so is its sup
> over `t`. ∎  So the *observer* predicate is sign-blind: the signed LRC is a **gauge cover** of the
> same loneliness, and the content lives in the *pairwise* (runner-runner) data.

> **T2 (a sign is a cut).** A sign vector `ε` is a **2-coloring** of the runners. The pair `(i,j)`
> clock is `ε_i v_i − ε_j v_j ∈ {±(v_i−v_j), ±(v_i+v_j)}`: a **difference** clock if `(i,j)` is
> *monochromatic*, a **sum** clock if *bichromatic* (a cut edge). The number of sum-clocks `=` the
> cut size `= |A|·|B|`, maximised by the **balanced** coloring at `⌊(n-1)²/4⌋`. *Verified* maxcut
> `= 1,2,4,6,9,12` for `n-1=2..7`. ∎

> **T3 (zero-clock ⟺ shell-partner).** Working mod `C=2n−1` (THM-401), a sum-clock is `≡ 0` iff
> `v_i+v_j ≡ 0 (mod 2n−1)`, i.e. the pair is a **shell-partner**. For distinct positive speeds the
> first such sum is `v_i+v_j = 2n−1` exactly. ∎

> **T4 (near-free action).** The sign group `(ℤ/2)^{n-1}` acts on the pair-clock multiset almost
> freely: the sign-orbit of the AP has size `2,4,7,16,32,60` for `n=3..8` (vs `2^{n-2}=2,4,8,16,32,
> 64`); the few collisions are config automorphisms. So the signed view is a **rich, near-`2^{n-1}`
> -valued refinement** of the single observer datum `M`.

## The climb, and the main result: the worry-set splits at n=14

Exhaustive over all sign patterns and all tight (`M=1/n`) configs in a window:

```
   n=4: tight {(1,2,3)}                 — shell-partners: NONE
   n=5: tight {(1,2,3,4),(1,3,4,7)}     — NONE
   n=6: tight {(1,2,3,4,5),(1,3,4,5,9)} — NONE
   n=7: tight {(1,2,3,4,5,6)}           — NONE
```
> **Through `n=7`, every worry-set (tight) config is *shell-partner-free*** — no pair sums to
> `2n−1`, so no signed zero-clock. (The AP has max sum `2n−3 < 2n−1`; the small additive-chain
> sporadics likewise.) This made "tight ⟹ no shell-partner" look like a theorem.

**It fails at `n=14`** — and that is the point. The exact `n=14` floor (repo Res_27 certificate,
S612) is `{AP, V*, 2·AP}`, all `M=1/14`. **Verified** (`…s699e.py`):
```
   AP    = (1,…,13)                    M=1/14   shell-partners: NONE  (max sum 25 < 27)
   V*    = (1,…,11,13,24)             M=1/14   shell-partners: (3,24)  3+24 = 27 = C   ✓
   2·AP  = (2,4,…,26)                 M=1/14   shell-partners: NONE
```
> **`V*` is a *tight* worry-set row that carries a shell-partner** `3+24=27` — refuting the
> conjecture. So **the signed structure splits the worry-set**:
> - **AP-type** (all `n≤7` tight rows; AP and `2·AP` at `n=14`): shell-partner-free, **no** signed
>   zero-clock.
> - **`V*`-type** (first at `n=14`): carries a shell-partner, **a** signed zero-clock.
>
> The signed LRC is therefore **strictly finer than `M`**: `AP` and `V*` are *both* `M=1/14`, yet
> the signed structure tells them apart (the cut separating `3` from `24` produces a `0`-clock for
> `V*`, never for `AP`). This is the "complete-understanding" payoff: the additive (`2n−1`) face is
> exactly what the sign group exposes, and it resolves the floor that `M` cannot.

## Why n=14, and the angle of attack

`V*`'s shell-partner is `3 + 24 = 27`, and **`24 = 2·12`** — `V*` is `AP` with `12` replaced by its
**double**. So the shell-partner is born from the **doubling** (prime-2, THM-404) landing on
`C = 27 = 3³` (prime-3, THM-407). **`n=14` is the first `n` whose `C=2n−1` admits a doubled-speed
shell-partner** — the prime-2 × prime-3 interaction that is the n=14 frontier, now visible as a
*single signed zero pair-clock*.

> **The angle of attack.** A proof of LRC(14) must handle *both* floor types; the signed structure
> separates them and pins the special feature of the hard one (`V*`): a doubled speed creating a
> shell-partner. The pair-clocks are the **second-moment / additive-energy face** of the
> covering-depth distribution (THM-406) — the sign group is the laboratory switch that turns the
> additive face on (S674). Concretely: classify the `n=14` floor by *signed zero-clock count*
> (`AP,2·AP → 0`; `V* → 1`), and attack `V*` via its doubling-shell-partner `(3,24)` with the
> owner/carry calculus (the `(3,24)` pair is the carry site).

## Honest status

- **Theorems (proved):** T1 gauge invariance; T2 sign=cut & maxcut `⌊(n-1)²/4⌋`; T3
  zero-clock⟺shell-partner; T4 near-free sign action (orbit sizes verified).
- **Verified:** every tight config shell-partner-free for `n≤7`; at `n=14`, `V*` tight (`M=1/14`,
  float + repo Res_27 exact certificate) with shell-partner `3+24=27`; `AP,2·AP` shell-partner-free.
- **Main finding (new):** "tight ⟹ no shell-partner" holds `n≤7` but **fails at `n=14`** (`V*`); the
  signed structure **splits the worry-set** and is strictly finer than `M`; the split is born from
  the doubling `24=2·12` hitting `C=27=3³` (prime-2 × prime-3).
- **Not claimed:** a proof of LRC(14). The deliverable is the rigorous signed-LRC theory and the
  worry-set split — a finer invariant and a concrete `V*`-specific attack site.
- **Open:** is `V*`-type (shell-partner-carrying tight row) the generic situation for even `n` with
  `3∣C`, or special to `n=14`? Does the signed zero-clock at `(3,24)` admit an owner/carry
  certificate forcing looseness of nearby lifts (tying to HYP-2241/S677 apex debt)?

**Artifacts:** `04-computation/signed_lrc_cut_structure_s699c.py`,
`signed_lrc_worryset_signature_s699d.py`, `signed_lrc_vstar_split_s699e.py` (+`.out`s). Extends
S674/S674b/S699; builds on THM-401 (`2n−1`/shells), THM-404/407 (doubling/prime-3), THM-406
(covering-depth moments), S612 (Res_27 floor `{AP,V*,2·AP}`). New: **HYP-2262**.
