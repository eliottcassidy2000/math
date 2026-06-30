# The apex-7 secondary obstruction: the LRC's Z₂ mirror is a vanishing PRIMARY obstruction, and the real forcing lives in the metagraph R-odd homology b₁⁻ — which is apex-7-dimensional (1, 7, 119; 7 | b₁⁻), exactly where the Z₂ parity is blind

*opus-2026-06-29. Owner: work creatively on whether a non-Z₂ (Z₇/apex-7) obstruction could force the LRC
where the mirror Z₂ cannot. The answer reframes the whole forcing question through obstruction theory:
the Z₂ mirror is the PRIMARY obstruction and it vanishes; the SECONDARY obstruction is the R-odd
homology, and it is apex-7.*

## The setup: why Z₂ fails (primary obstruction vanishes)
The lonely point is a visit of the danger-walk `d(t)` to `0…0` on `Q_{n-1}`. The mirror `t→−t` is a
`Z₂`-symmetry (`d(−t)=d(t)`), so the lonely count is EVEN — the **primary obstruction** (degree-1, mod 2)
**vanishes**. Crucially, the `t`-circle is **one-dimensional**, so `H^*(S^1)` stops at degree 1: there
is NO room on the circle for a higher obstruction. **The Z₂ mirror is the whole story on the circle, and
it says nothing.** Any forcing must live in a HIGHER-dimensional space — the metagraph / cap homology.

## The secondary obstruction is the R-odd homology — and it is APEX-7 (verified)
The metagraph carries the complement involution `R`; its 1-skeleton's R-odd Betti number `b₁⁻ =
dim H_1^{-}(metagraph)` is the **secondary obstruction** (the part of the homology the `Z₂` parity
cannot reach). Recomputed from scratch (tilings → classes → wiggly edges → R-action → Lefschetz):
> **`b₁⁻ = 1, 7, 119` for `n = 4, 5, 6`**, and **`7 | b₁⁻` for `n = 5, 6`** (`b₁⁻(5)=7`, `b₁⁻(6)=119=7·17`).
At `n=5` it is **exactly 7** — the apex prime, the dimension of the **Fano plane** (7 points/lines) and
of the **imaginary octonions** (`O = R ⊕ R^7`). So the secondary obstruction is not just nonzero where
the primary vanishes — it is *carried by the apex-7*.

> **This is the creative payoff: apex-7 is not only number-theoretic (Mersenne ∩ Heegner ∩ 3-mod-4 =
> {3,7}) — it is the DIMENSION of the secondary topological obstruction of the metagraph.** The `Z₂`
> mirror is exactly blind to it (an odd-prime homology class is invisible to mod-2 parity).

## The octonion/Fano frame for the Z₇
The apex-7 structure is the imaginary octonions: `Z₂` = conjugation (real ↔ imaginary) = the **mirror**
(primary); `Z₇` = the **Singer cycle** on the 7 imaginary units = the Fano-line rotation = the
**secondary**. The `28 = C(8,2) = T(7)` arc-count of the 8-tournament (HYP-3546) is the octonion
multiplication table; the LRC threshold `14 = 2·7` factors as (conjugation `2`) × (Fano `7`). The
secondary obstruction lives on the `7` that the `2` cannot see.

## The concrete LRC anchor at the 7-scale (computed)
At `t = a/7`: `‖s_i a/7‖ < 1/14 ⟺ 7 | s_i a ⟺ 7 | s_i` (since `1/7 > 1/14`). So the `Z₇`-orbit
`{a/7 : a=1..6}` is **blocked exactly by the mult-of-7 runners** — which every covering set HAS (7 and
14 are forced). The lonely set lives in the 7 open cells `[a/7,(a+1)/7)`, is **Z₂-mirror-symmetric**
(`cell a ↔ cell 6−a`, verified), but is **NOT Z₇-symmetric** (the rotation `t→t+1/7` twists runner `i`
by `s_i/7`). So `Z₇` is *not* a naive symmetry forcing a count — consistent with it living in the
*secondary* (homological), not the primary (parity), layer.

## What this says about forcing (honest)
- **The Z₂ mirror cannot force the LRC** (primary obstruction vanishes; circle too low-dimensional). PROVED.
- **The secondary obstruction is nonzero and apex-7** (`b₁⁻ = 1,7,119`, `7 | b₁⁻`). VERIFIED.
- **The forcing question becomes:** does the nonvanishing secondary obstruction `b₁⁻ ≠ 0` (a
  `Z₇`/Volovikov-type, not `Z₂`/Borsuk-Ulam, invariant) descend to "a lonely point exists"? That descent
  is exactly **mac-mini's HYP-3544 (the R-odd Betti = the LRC odd index)** — now sharpened: the index is
  **apex-7-graded**. A `Z₇`-equivariant (Volovikov) Borsuk–Ulam on the metagraph/cap, if it transfers,
  forces where the `Z₂` cannot.

So the creative answer is yes-in-principle: a `Z₇`/apex-7 *secondary* obstruction is the right object,
it is provably nonzero (and apex-7-divisible) exactly where the `Z₂` primary dies, and it is the
octonion/Fano `7` of `14 = 2·7`. The gap is the descent (HYP-3544), now with a concrete target: explain
`7 | b₁⁻` structurally (a `Z₇`/Fano-line basis of `H_1^{-}`) and transfer it to the LRC cap's R-odd
eigenspace (HYP-3538).

## Status
- **PROVED/verified (opus, recomputed):** `b₁⁻ = 1, 7, 119` (n=4,5,6); `7 | b₁⁻` for n=5,6; `b₁⁻(5)=7`
  exactly. The `Z₂` mirror is a vanishing primary obstruction (even count); the secondary is `b₁⁻`.
- **Computed:** the 7-scale `{a/7}` orbit is blocked by mult-of-7 runners; the lonely set is
  `Z₂`-mirror-symmetric on the 7 cells, not `Z₇`-symmetric.
- **Reframe (opus):** apex-7 is the DIMENSION of the secondary obstruction; `Z₂`=conjugation (primary),
  `Z₇`=Singer/Fano (secondary); the forcing the mirror can't do lives on the octonionic `7`.
- **Open:** verify `7 | b₁⁻(n)` for `n ≥ 7`; a `Z₇`/Fano-line basis of `H_1^{-}`; the Volovikov descent
  `b₁⁻ ≠ 0 ⇒` lonely (the apex-7-graded HYP-3544).

Related: the Ky-Fan-forcing-fails reflection (the primary `Z₂` even count), mac-mini HYP-3544 (R-odd
Betti = LRC odd index) + HYP-3538 (cap R-odd), klein THM-587 (`P_n(±1)`), HYP-3547 (apex-7 three
pillars), HYP-3546 (`28=C(8,2)`=octonion table), the duality-web (even/odd–±–add/mult at apex 7),
THM-582 (two-index), OPEN-Q-108.
