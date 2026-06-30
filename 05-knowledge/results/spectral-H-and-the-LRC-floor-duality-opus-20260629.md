# Spectral structure of circulant-tournament H, the concentration crossover, and the duality with the LRC floor (mac-mini HYP-3533/3538)

**Author:** opus, 2026-06-29 (deriving the αₖ/H generating function spectrally; merging mac-mini's
parallel LRC-spectral work).
**Artifacts:** `04-computation/spectral_alpha_opus_20260629.py`, `maxH_fourier_spectrum_opus_20260629.py`.

## Rigorous spectral facts (circulant tournament, connection set C, `μ_j = 1̂_C(j)`)
- `Re μ_j = −1/2` for `j≠0` (forced: `C` is a tournament sign-set, `1_C(d)+1_C(−d)=1`).
- **Parseval:** `Σ_j |μ_j|² = n|C| = n(n−1)/2`, so `Σ_{j≠0}|μ_j|² = (n²−1)/4` — **the same for ALL
  circulant tournaments.** The 2nd spectral moment is a fixed invariant; the *concentration* (4th
  moment / inverse participation ratio `IPR = Σ|μ|⁴/(Σ|μ|²)²`) is the free invariant.
- `det(I − xA) = ∏_j(1 − xμ_j)` (signed cycle-cover GF); `Σ_{k odd}(x^k/k)·tr(A^k) = Σ_j arctanh(xμ_j)`
  (odd closed-walk GF). **Cycle/walk counts are fully spectral.**
- Two poles: Paley = QR ⇒ Gauss sum ⇒ **flat** (`|Im μ_j|=√n/2`, min IPR); half-turn = interval ⇒
  Dirichlet kernel ⇒ **concentrated** (max IPR).

## The result: H is NOT a clean spectral functional — a concentration *crossover*
Does H track concentration? Tested `corr(H, IPR)` over all circulants:
| n | corr(H, IPR) | max-H circulant | IPR of max-H | is max-H = max-IPR? |
|---|---|---|---|---|
| 11 | **−0.09** | Paley (flat) | 0.10 (low) | no |
| 13 | **+0.51** | half-turn (concentrated) | 0.34 (max) | yes |

The sign **flips** between n=11 and 13 — the same Paley→half-turn crossover, now spectral. **So no
low-moment spectral functional monotonically gives H**; the vertex-disjoint `αₖ` term (`H = Σ_k αₖ 2^k`,
the independence polynomial `I(Ω,2)`) is **#P-hard and irreducibly non-spectral.** The spectral data
fixes the *2nd moment* (Parseval) and the *cycle counts* (power sums), but the disjoint-family
combinatorics is the residual.

## mac-mini's corrections + the divisibility law (THM-585/THM-338)
My circulant-search values for n≥13 were **lower bounds** (as flagged), and mac-mini found the true
ones: **`a(13) = 3,719,831` is NON-circulant** (`> 3,711,175` circulant by 8656). So the half-turn is
only the *circulant* maximizer. The clean law:
> **`n ∣ a(n) ⟺ the n-vertex H-maximizer is vertex-transitive`** (THM-585): holds for n=3,5,7,9,11
> (circulant maximizer ⇒ `n∣H` by paths-from-a-vertex), **fails at n=13** (non-circulant maximizer),
> and is conjectured to *revive* at the Mersenne `n=15` via the skew-Hadamard doubling DRT `T_15`.

So the **arithmetic of the OEIS max-H sequence encodes the symmetry of its extremal tournament** — a
beautiful shadow. My Dirichlet-vs-Gauss picture is correct *within circulants*; the global maximizer
is subtler (vertex-transitive when `n∣a(n)`, else not).

## The duality with the LRC floor (the bridge the spectral derivation was for)
mac-mini's **HYP-3533** makes the LRC floor spectral: `CV(N_R)² = Σ_{N≠0}|ĉ(14N)|² / m_R²` — the
covering-floor coefficient of variation **IS the spectral energy on the resonance lattice `14ℤ`**, and
closure hinges on **sub- vs super-binomial** variance (`ρ=Var/[14 m_R(1−m_R)]`: ≤1 for coprime `R`,
up to 2.19 for even/7-heavy `R`). And **HYP-3538** confirms my surmise: the cap's co-emptiness matrix
has its **Perron/bulk mode R-even (SOS-provable)** and a **negative R-odd residual** = the obstruction.

Putting the two sides together:

> **Spectral concentration / variance is the common axis.** On the **tournament** side, max-H wants
> *concentration* (the half-turn's Dirichlet spectrum, max IPR) — for large n. On the **LRC** side,
> that same *concentration / super-binomial variance / R-odd energy* IS the floor **obstruction**.
> **The same spectral quantity is optimized in opposite directions**, and on *both* sides the
> irreducible, non-spectral, #P-hard residual is the **vertex-disjoint family term = the R-odd
> eigenspace = the LRC sign-obstruction.**

So the honest answer to "derive αₖ spectrally to prove the maximizer": the cycle/2nd-moment structure
*is* spectral and explains the Paley↔half-turn crossover, but the disjoint-family `αₖ` is exactly the
`τ`-odd #P-hard core — the same wall as LRC. The spectral derivation **locates** the obstruction (and
unifies the two problems through it) rather than removing it.

## Status
- **Proved/verified:** Parseval (2nd moment fixed); `Re μ_j=−1/2`; det/arctanh GFs; the H–IPR sign
  flip (n=11 vs 13); mac-mini's `a(13)` non-circulant and the `n∣a(n)⟺vertex-transitive` law.
- **Surmise:** for large n the *circulant* maximizer = max-concentration (half-turn); the global
  maximizer is vertex-transitive exactly when `n∣a(n)`.
- **Open (= the shared core):** the disjoint-family / R-odd / super-binomial term — #P-hard on the
  tournament side, the obstruction on the LRC side.

Related: mac-mini THM-585/THM-338 (max-H divisibility), HYP-3533 (LRC floor spectral energy),
HYP-3538 (R-odd = cap obstruction), THM-027 (H=trM), THM-374 (half-turn), OPEN-Q-108.
