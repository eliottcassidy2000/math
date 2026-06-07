# THM-421 — LRC clock-witness leaks are localized to the prime-torsion of ℤ/n (CRT)

**Status:** PROVED (the CRT localization + the geometric margin) + VERIFIED (n=12,14,15,18,20,21,30,
33,35,45,50). Sharpens the k-clock witness THM-420/Lemma A by saying *where* the obstruction lives.
**Source:** opus-2026-06-06-S701, from the user's observation: "the n=14 half-turn leak at residue 7
is 1 mod 2 and 0 mod 7 — the 2-torsion projecting to zero in the 7-runner base; n=15's order-3 leaks
at 5,10 are the 3-torsion projecting to zero in the 5-runner base; the leakage is locked into the
algebraic torsion of the composite divisors."

## Setup

`n` runners, speeds `v_i`, full clock `t = 1/n`; runner `i` sits at `v_i/n`, position governed by
`v_i mod n`. A **clock leak** at prime `p ∣ n` = a runner that defeats every clock dividing the
`p`-cofactor (`≡ 0` there) yet survives in the `p`-direction. For `n = ∏_p p^{a_p}`, CRT gives
`ℤ/n ≅ ∏_p ℤ/p^{a_p}`; the **`p`-torsion subgroup** is `T_p = {x : p·x ≡ 0 (mod n)} = ⟨n/p⟩`, of
order `p`.

## The theorem

> **(1) CRT localization.** Every nonzero `x ∈ T_p` is `≡ 0` in the cofactor base
> `m_p = n/p^{a_p}` and nonzero in the `p`-base. In particular the generator `n/p` is the residue
> the user named: `n=14, p=2 ⟹ n/p = 7 ↔ (1 mod 2, 0 mod 7)`; `n=15, p=3 ⟹ n/p = 5 ↔ (2,0)` and
> `10 ↔ (1,0)` mod `(3,5)`.
>
> **(2) Geometric margin (all `n`).** At the full clock `t = 1/n`, a `p`-torsion runner `x ∈ T_p\{0}`
> sits at an **exact order-`p` rotation** `x/n = j/p` (`j = x/(n/p)`), so `‖x/n‖ ≥ 1/p`. The
> `p`-torsion leak is *intrinsically safe by margin `1/p`*.
>
> **(3) Squarefree vs prime-power dichotomy.** The canonical generator `n/p = p^{a_p−1} m_p`:
> - **`a_p = 1` (squarefree `p`):** `n/p = m_p` is coprime to `p`, so the socle *survives mod `p`* —
>   the leak is plugged by the **prime clock `t = 1/p`** (margin `1/p`). (n=14, 15, 21, 30, 33, 35.)
> - **`a_p ≥ 2` (prime-power `p`):** `n/p = p^{a_p−1}m_p ≡ 0 (mod p)`, so the socle *dies mod `p`*
>   too — the prime clock fails on it; control needs the **deeper `p`-adic clock** `t = 1/p^{a_p}`.
>   (n=12, 18, 20, 45, 50: socle `survive mod p = False`, verified.)
>
> **(4) The hard core.** The runners that leak at *every* prime simultaneously are
> `⋂_p (≡0 \text{ in the }p\text{-comp}) = \{x ≡ 0 (\bmod\, \mathrm{rad}\,n)\}` (multiples of the
> radical). For **squarefree `n` this is `{0}` only** — a runner is all-prime-dangerous iff
> `v ≡ 0 (mod n)`, recovering exactly the **THM-398 multiple-of-`n` residual**. For prime-power-
> containing `n` the hard core is larger (e.g. `{0,6}` at `n=12`).

## Proofs

**(1)** `n/p` divides `x`, and `n/p = p^{a_p−1}m_p` is a multiple of `m_p`, so `x ≡ 0 (mod m_p)`.
`x = j·(n/p)` with `1 ≤ j ≤ p−1`; in `ℤ/p^{a_p}` the element `n/p ≡ p^{a_p−1}m_p` is nonzero (`m_p`
invertible mod `p`), so `x ≢ 0` in the `p`-base. ∎
**(2)** `x/n = j(n/p)/n = j/p`, an order-`p` point of the circle; `‖j/p‖ = \min(j,p−j)/p ≥ 1/p`. ∎
**(3)** `gcd(n/p, p) = 1 ⟺ a_p = 1`; when `a_p = 1`, `(n/p)\bmod p = m_p \bmod p ≠ 0` so `t=1/p`
gives `‖x/p‖ ≥ 1/p`; when `a_p ≥ 2`, `p ∣ n/p` so `x ≡ 0 (mod p)` and `t=1/p` sends `x` to the
origin — the least controlling clock exceeds `p` (verified). ∎
**(4)** CRT: `x ≡ 0` in every `p`-component `⟺ x ≡ 0 (mod \prod_p p) = \mathrm{rad}\,n`. ∎

## Significance

- **Confirms the user's thesis:** the LRC clock-leaks are *completely locked into the algebraic
  torsion of the composite divisors* — not scattered across `ℤ/n` but supported on `⋃_p T_p`, each
  piece a single-prime-surviving CRT element.
- **The prime-power face is the hard one (3):** squarefree `n` plugs every prime-torsion leak with a
  prime clock; prime powers force descent down the `p`-adic tower. This is the *same* prime-power
  obstruction as THM-420's caveat — there the modulus `2n−1 = 27 = 3³` (for `n=14`) is a pure prime
  power with no coprime plug; here it is the `p`-adic depth of the torsion. **Both faces (the clock
  `ℤ/n` and the shell `ℤ/(2n−1)`) hide LRC's difficulty in prime-power torsion towers.**
- **Recovers THM-398 cleanly (4):** for squarefree `n`, "all clocks leak" forces `v ≡ 0 (mod n)`,
  the multiple-of-`n` residual — so the THM-398 residual is the CRT hard core, `= {0}` exactly when
  `n` is squarefree.

**Tournament reading (per the standing directive):** the order-`p` rotation clock `t = j/p` is the
`C_p` cyclotomic **round tournament** (THM-403); `T_p` is the lattice of cyclotomic sub-clocks. The
leak-localization says the worry-set's cyclotomic witnesses are organized by the prime-torsion
lattice of `ℤ/n`; the prime-power depth (3) = the height of the cyclotomic tower `C_p ⊂ C_{p²} ⊂ …`.

**Artifacts:** `04-computation/lrc_torsion_localization_s701.py` (+`.out`). Builds on THM-420
(k-clock witness / Lemma A), THM-398 (multiple-of-`n` residual), THM-403 (cyclotomic witnesses).
Companion reflection `lrc-clock-leaks-are-prime-torsion-s701.md`. New: **HYP-2285**.
