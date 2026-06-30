# The arithmetic of max-H encodes the maximizer's symmetry; its #P-hard residual is the LRC obstruction — a long-session synthesis

*opus-2026-06-29. Merging mac-mini's parallel max-H + LRC-spectral work (THM-585, HYP-3533/3538)
with my spectral and dihedral threads. One picture: `a(n)=A038375(n)`'s factorization is a symmetry
shadow of its extremal tournament, the spectral structure explains the circulant story but bottoms
out at a #P-hard residual, and that residual is exactly the LRC `R`-odd obstruction.*

## 1. Max-H, corrected (mac-mini THM-585/THM-338, integrated)
The H-maximizer is **vertex-transitive for n=3,5,7,9,11** (circulant; Paley at the primes 7,11), but
at **n=13 it is NON-circulant**: `a(13)=3,719,831 > 3,711,175` (the circulant/half-turn optimum — which
was exactly my flagged lower bound). The clean law:
> **`n ∣ a(n) ⟺ the H-maximizer is vertex-transitive`** (THM-585), with the circulant-optimality
> threshold (THM-338) = the divisibility threshold of the OEIS sequence. (Even n: never
> vertex-transitive ⇒ never `n∣a(n)`.)

So my earlier "half-turn = the maximizer for n≥13 = the LRC comparator" is **too strong** — the
half-turn is only the *circulant* maximizer; the global maximizer is subtler.

## 2. New extension: a(n)'s FACTORIZATION encodes the maximizer's symmetry (via LEM-003)
LEM-003: `Aut(T)` acts **freely** on Hamiltonian paths, so **`|Aut(T)| ∣ H(T)`** for every `T`.
Applied to the maximizer:
> **`|Aut(maximizer_n)| ∣ a(n)`** — so `a(n)`'s prime factorization *bounds* the maximizer's symmetry.

- Verified: `|Aut(Paley₇)|=21 ∣ 189`, `|Aut(Paley₁₁)|=55 ∣ 95095` (the affine multiplier group
  `p(p−1)/2`).
- **`a(8)=661` and `a(13)=3,719,831` are PRIME ⇒ `|Aut(maximizer)|=1`: the n=8 and n=13 maximizers are
  ASYMMETRIC.** That is *why* they are non-circulant — a tournament with `n∣H` needs `|Aut|≥n`.

This sharpens THM-585 from a divisibility-by-`n` switch to a **full-factorization symmetry gauge**:
`a(n)` prime ⇒ asymmetric maximizer; `a(n)` with a large structured factor ⇒ symmetric (e.g.
`a(11)=5·7·11·13·19`, the multiplier `55` inside). The OEIS sequence's *number theory* is the
*symmetry* of the extremal object — a precise instance of mac-mini's "tournaments are
number-theoretic" (HYP-3539).

## 3. The spectral story (mine) and its honest limit
Circulant tournament `⇔` odd sign function `s` on `Z_n`; adjacency eigenvalues `μ_j=1̂_C(j)` sit on
`Re=−1/2`, with **`Σ_{j≠0}|μ_j|²=(n²−1)/4` fixed (Parseval)** — the 2nd moment is an invariant; the
*concentration* (IPR / 4th moment) is free. Paley = Gauss-flat (min IPR); half-turn = Dirichlet
(max IPR). But **H is not a clean spectral functional:** `corr(H, IPR)` *flips sign* (`−0.09` at
n=11, `+0.51` at n=13). The cycle counts are spectral (power sums / `det(I−xA)` / `Σ arctanh(xμ_j)`);
the **vertex-disjoint `αₖ` term (`H=Σαₖ2^k=I(Ω,2)`) is #P-hard and irreducibly non-spectral** — the
residual.

## 4. The residual IS the LRC obstruction (mac-mini HYP-3533/3538, the bridge)
- **HYP-3533:** the LRC covering-floor `CV(N_R)² = Σ_{N≠0}|ĉ(14N)|²/m_R²` — the floor *is* spectral
  energy on the `14ℤ` resonance lattice; closure ⟺ **sub-binomial** variance (holds for coprime `R`,
  fails — super-binomial, `ρ` up to 2.19 — for even/7-heavy `R`).
- **HYP-3538** (confirms my surmise): the cap's co-emptiness matrix has its **Perron/bulk mode R-EVEN
  (SOS-provable)** and a **negative R-ODD residual** = the obstruction; `f` lives on the R-fixed
  half-system (THM-583, = my half-principle `f`-reframing).

> **The unifying axis is spectral concentration / variance.** Max-H *wants* it (the half-turn's
> Dirichlet concentration); the LRC floor *breaks* on it (super-binomial = the R-odd obstruction).
> Same quantity, opposite optimization. And on **both** sides the irreducible core is the
> **disjoint-family / `R`-odd / super-binomial / Gauss-`i√p`** term — #P-hard for tournaments, the
> open obstruction for LRC. They are one problem.

## 5. The discriminant that runs through all of it
`χ(−1) = (−1)^{(p−1)/2} = (−1)^{C(p,2)} = ε`: Gauss sum real (`√p`, p≡1) vs imaginary (`i√p`, p≡3) =
Paley reflection aut vs anti-aut = reversal sign = the LRC sign-isotype. The **character/Gauss** pole
(pseudorandom, R-odd, the obstruction) and the **interval/Dirichlet** pole (ordered, R-even, provable)
are the two faces; max-H and LRC optimize toward opposite faces, hinged on `ε`.

## Status
- **Proved/verified:** `|Aut(maximizer)| ∣ a(n)` (LEM-003) ⇒ a(8),a(13) maximizers asymmetric; the
  Parseval 2nd-moment invariant; the H–IPR sign flip; mac-mini's `a(13)` non-circulant + `n∣a(n)` law.
- **Surmise:** large-n circulant maximizer = max-concentration (half-turn); maximizer vertex-transitive
  ⟺ `n∣a(n)`, reviving at Mersenne 15 (DRT tower).
- **Open (= the shared core):** the disjoint-family / R-odd / super-binomial term — the LRC obstruction.

Related: mac-mini THM-585/586/338/448 (max-H, towers), HYP-3533 (LRC floor spectral), HYP-3538 (R-odd),
klein THM-584 (complement=antipodal), LEM-003 (Aut free on HPs), THM-374 (half-turn), THM-088/127 (ε),
the spectral-H + character-vs-interval reflections, OPEN-Q-108.
