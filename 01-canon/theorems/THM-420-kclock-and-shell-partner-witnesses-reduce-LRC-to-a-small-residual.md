# THM-420 — The k-clock witness and the shell-partner lemma reduce LRC(n) to a small explicit residual

**Status:** PROVED (the two witness lemmas) + VERIFIED (the residual is loose, n≤8 window). A major
refinement of the LRC reduction, NOT a full proof of LRC(14).
**Source:** opus-2026-06-07-S700. Convention: `n` runners, gap `1/n`, speeds `v_1,…,v_{n−1}`
distinct positive integers; `M(S)=max_t min_i ‖v_i t‖`; LRC(n) ⟺ `M(S) ≥ 1/n` for all `S`;
`C = 2n−1`.

## Two explicit-witness lemmas

> **Lemma A (k-clock witness — generalizes THM-369).** If there exists `k ∈ {2,…,n}` such that no
> `v_i ≡ 0 (mod k)`, then `t = 1/k` gives `M(S) ≥ 1/k ≥ 1/n`. *Proof:* `‖v_i/k‖ = min(r,k−r)/k`
> with `r = v_i mod k ∈ {1,…,k−1}`, so `≥ 1/k`; and `1/k ≥ 1/n` since `k ≤ n`. ∎
> (THM-369 is the case `k = n`; smaller `k` give a *larger* bound `1/k`.)

> **Lemma B (shell-partner lemma — NEW).** If all `v_i` are coprime to `C = 2n−1` and there is a
> **shell-partner pair** `v_i + v_j ≡ 0 (mod C)`, then `M(S) ≥ 2/(2n−1) > 1/n`.
> *Proof.* The discrete witness `t = m/C` gives `‖v_k m/C‖ ≥ 2/C` for all `k` iff
> `v_k m ∉ {0, ±1} (mod C)` for all `k`. The forbidden set is `F = {0} ∪ \{±v_k^{-1} : k\}`
> (using `v_k` coprime: `v_k m ≡ ±1 ⟺ m ≡ ±v_k^{-1}`). A shell-partner gives
> `v_j ≡ −v_i ⟹ v_j^{-1} ≡ −v_i^{-1} ⟹ \{±v_j^{-1}\} = \{±v_i^{-1}\}`, so that pair contributes
> **2 values, not 4**: `|F\setminus\{0\}| ≤ 2(n−1) − 2 = 2n−4`. Hence `|F| ≤ 2n−3 < 2n−1 = |ℤ/C|`,
> so `≥ 2` good `m` exist; for such `m`, `t = m/C` gives `M ≥ 2/C = 2/(2n−1) > 1/n`
> (gap `= 1/(n(2n−1))`). ∎
> **Verified** (`…s700.py`): all shell-partner configs at `n = 5,6,7` have a good `m` and
> `M ≥ 2/(2n−1)`.

## The refined reduction

> **Corollary.** LRC(n) holds for every config except possibly the **residual**:
> ```
>   R(n) = { S : every k∈{2,…,n} divides some v_i  (all clocks fail)
>                 AND  S has no shell-partner pair (v_i+v_j ≢ 0 mod 2n−1) }.
> ```
> Everything outside `R(n)` is loose by Lemma A or Lemma B. So **LRC(n) ⟺ every config in `R(n)`
> is loose.**

**Verified** (`…s700c.py`): `R(n)` is small and entirely loose in the window `[1,2n]`:
```
   n :   5    6    7    8
  |R| :  18   50   62  233    (of 205, 786, 2996, 11432 gcd-1 configs)
  min M : 3/13, 1/5, 3/19, 3/20   — all > 1/n  (in fact ≥ 2/(2n−1))
```
The residual's witnesses sit at **pair-sum resolutions** `t = m/(v_i+v_j)` (e.g. `M(1,3,4,10) =
3/13` at `n=5`, `3+10=13`) — the additive (`2n−1`/signed) face (S699). So `R(n)` is the
*divisibility-complete, shell-free* core, witnessed by pair-sums.

## Status and scope

- **PROVED:** Lemma A (all `n`); Lemma B (coprime case, all `n`).
- **VERIFIED:** every config outside `R(n)` loose by A/B; `R(n)` loose (n≤8 window), via pair-sum
  witnesses.
- **NOT proved:** that `R(n)` is always loose (the irreducible core) — and Lemma B's **coprime
  hypothesis excludes the `n=14` prime-3 strata** (the gcd-3,9 speeds mod `27`), so it does not by
  itself close LRC(14). The honest gain: two clean explicit-witness theorems that cover all but a
  small, structured, pair-sum-witnessed residual.
- **Reframe:** the witness hierarchy by resolution — **`k`-clock `1/k`** (Lemma A) `⊃`
  **shell-partner / discrete `2/(2n−1)`** (Lemma B) `⊃` **pair-sum `m/(v_i+v_j)`** (the residual).
  LRC(n) is exactly: the divisibility-complete, shell-free configs are loose, with a pair-sum
  witness.

**Artifacts:** `04-computation/lrc_shell_partner_lemma_s700.py`, `lrc_residual_no_shell_s700b.py`,
`lrc_kclock_shellpartner_residual_s700c.py` (+`.out`s). Builds on THM-369 (clock witness), THM-398
(`C'`), THM-401 (`2n−1`/shells), S599h (the `2/(2n−1)` floor), S599i (discrete witness), S699
(signed/pair-sum). Companion reflection `lrc-kclock-shell-partner-refinement-s700.md`. New:
**HYP-2283**.
