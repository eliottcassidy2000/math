# Computing A038375 (max Hamiltonian paths) via the H-atom recursion + new tournament sequences

**Author:** opus, 2026-06-29 (H-atom recursion → extremal sequences).
**Artifacts:** `04-computation/maxH_circulant_opus_20260629.py`,
`04-computation/maxH_rotational_sequence_opus_20260629.py` (+ `.out`),
`condensation_monoid_opus_20260629.py`, `H_factorization_opus_20260629.py`.

## The technique: max-H lives on a single circulant prime

`H` is **multiplicative over the condensation** (`H(X⊕Y)=H(X)·H(Y)`, verified n≤6), so the
maximizer is a **single strongly-connected prime** — and (Paley-maximizer lore) the maximizer
is **circulant**. Therefore:

> **A038375(n) = max H over the `2^{(n-1)/2}` circulant tournaments** (odd n), not over all
> `2^{C(n,2)}` tournaments. Compression: n=13 → **64** circulants vs `~10^{23}` tournaments;
> n=19 → 512 vs `~10^{51}`.

**Validated:** circulant search reproduces A038375 exactly for n=3,5,7,9 (`3,15,189,3357`) and
gives `95095 = H(Paley₁₁)` at n=11 (the known max). **Extended (circulant-optimal; lower bounds
on A038375, conjecturally exact):**
| n | 11 | 13 | 15 |
|---|---|---|---|
| max-H (circulant) | 95095 | **3,711,175** | **198,464,295** |

## New structural finding: the maximizer is Paley OR rotational

The maximizing circulant's connection set:
- **prime n ≡ 3 (mod 4)** (7, 11): the **Paley/QR** tournament (out-set = non-residues). Known.
- **n = 13, 15**: the **rotational** tournament `R_n` (i beats the previous `(n-1)/2`); the QR
  circulant is strictly *sub*-maximal there (`13`: QR `3669497 < 3711175`).
- **n = 9**: a third circulant (`out={1,5,6,7}`), neither Paley nor rotational.

So "Paley maximizes H" is **family-specific** (prime ≡3 mod 4), not universal.

## New sequences

1. **Rotational-tournament Hamiltonian-path count `R_n`** (i beats next `(n-1)/2`), n=3,5,…,21:
   `3, 15, 175, 3267, 93027, 3711175, 198464295, 13689269499, 1184212824763, 125547534942879`.
   It equals A038375 at n=3,5,13,15 but is sub-maximal at n=7,9,11. **Honest negative:** `R_n`
   is *not* low-order P-recursive (order≤3, deg≤3, strict held-out validation fails) — the
   growing connection band defeats a fixed-size transfer matrix, so no cheap recurrence.
2. **Min-H over circulants** ("smallest rotational core"), n=3,5,…,15:
   `3, 15, 175, 3159, 92411, 3669497, 193215375`.
3. **H-atom spectrum** = the set of `H`-values of strongly-connected primes per size
   (n=5: `{9,11,13,15}`; n=6: `{15,17,19,23,25,27,29,31,33,37,41,43,45}` — a *proper* subset of
   the odds, missing 21,35,39). Since every tournament `H` is a product of these odd atoms, the
   achievable-`H` set is the multiplicative monoid they generate.
4. **#SC primes** `SC(n)=1,0,1,1,6,35,353` (the irreducibles; = #H-atoms).

## The Alon envelope (data point)

`H_max(n)·2^{n-1}/n!` = `2, 2, 2.4, 2.37, 2.44, 2.44, 2.49` (n=3..15) — slowly increasing,
consistent with Alon's `Θ(n^{3/2})` envelope (the `/n^{1.5}` column decreases monotonically).

## The catalog of extensive tournament variables (from `05-knowledge/variables/`)

`H(T)` (→A038375), `F(T,x)`/`F_k` (forward polynomial), `W(T,r)` (weighted HP), `SF`/`S(T)`
(signed HP permanent), `t_k` (k-cycle counts; `t_3=C(n,3)/4 ⟺ regular`), `M[a,b]` (transfer
matrix, `tr M=H`), `alpha_k` (disjoint odd-cycle tuples, `H=1+2Σα_k`), `c_k`/`D_k` (W-Fourier),
`mu(C)` (OCF cycle weight), `Ω(T)` (conflict graph, `H=I(Ω,2)`), `d(T)` (normalized
determinant, `E[d]=A000085`), good-cut/bucket/residue-rank invariants. Condensation-multiplicative
ones (`H`, and `F` up to shift) are exactly the ones the prime-factorization computes fast.

## Status
- **Verified:** A038375 via circulant = known values n≤11; `H` multiplicative over primes.
- **Extended (conjecturally exact):** A038375(13)=3711175, A038375(15)=198464295.
- **New:** the rotational `R_n` sequence, min-H-circulant, the H-atom spectrum, the
  Paley-vs-rotational maximizer dichotomy.
- **Negative (honest):** `R_n` not low-order holonomic; circulant values for n≥13 are lower
  bounds (the maximizer is only conjectured circulant for non-Paley n).
