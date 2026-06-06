---
id: HYP-2321
title: Quadratic reciprocity = the master 2-adic-seam law; the Legendre symbol is the XNOR with an out-of-phase 0-middle (the Paley conference matrix)
status: OPEN (synthesis); Legendre-XNOR / Gauss-phase / conference-matrix / reciprocity all VERIFIED
source: claudebox-2026-06-03-S638
related:
  - THM-126  # Paley circulant maximizes H; THM-162 Paley flat spectrum
  - THM-135  # Paley local-max but NOT global H at p=19 (= metastable, S637 glass)
  - HYP-2316 # the polarized delta field / glass transition (the 2-adic seam)
  - HYP-2306 # parity-defect ladder (even/odd); HYP-2301 chi-bridge (CM fields)
  - HYP-2311 # character-ratio spectrum bounds H/chi_di (Gauss sum = the spectrum)
---

# HYP-2321 — reciprocity, XNOR, and the 2-adic seam

The user's "binary choice with an always-on middle out of phase with both" is, precisely, the
**quadratic character / Legendre symbol**, and quadratic reciprocity is the law that couples its
phase across primes — the same even/odd 2-adic seam running our whole program.

## The Legendre symbol IS the XNOR-with-out-of-phase-middle (VERIFIED)

- **XNOR = `±1` multiplication** (agreement). The Legendre symbol is a *multiplicative* `±1` character:
  `(ab/p)=(a/p)(b/p)` (verified, 0 violations) — a global XNOR on residue classes.
- **The always-on middle `0`**: `(0/p)=0`, present for the ramified value, **out of phase** with `±1`.
- **Matrix form = the Paley conference matrix** `C_{ij}=(j-i/p)`: `±1` off-diagonal (the binary
  choice), `0` diagonal (the always-on middle), `C Cᵀ = (p-1)I - J` (orthogonality = "out of phase").

## Reciprocity governs the PHASE = the 2-adic seam (VERIFIED)

The Gauss sum `g_p = Σ_a (a/p) ζ^a`:
- **`p ≡ 1 mod 4`:** `g_p = √p` (**real, in phase**); `−1` is a QR; the character is **symmetric**;
  `C` is **symmetric** = the Paley **graph**. The "even" side.
- **`p ≡ 3 mod 4`:** `g_p = i√p` (**imaginary, out of phase, arg 90°**); `−1` is a non-residue; the
  character is **skew/antisymmetric**; `C` is **skew** = the Paley **tournament**. The "odd" side.

So the **out-of-phase middle is `i√p`**, appearing exactly on the skew (tournament) `p≡3` side. The
switch is `(p-1)/2 mod 2 = [p≡3 mod 4]` — **the 2-adic seam**: the same even/odd that gives
`χ(C_n)=2/3` (HYP-2295), the alternating-group parity defect (HYP-2306), the LRC even-`n` hardness,
and the glass transition at even `n` (HYP-2316). Reciprocity's sign flip,
`(-1)^{((p-1)/2)((q-1)/2)}`, is the **AND of the two oddness bits** (`p≡3 ∧ q≡3`): two out-of-phase
characters interfere; otherwise they commute.

## Implications across the program

1. **Paley = the QR-trienement = the conference matrix = the spectral extremum.** Its character-ratio
   spectrum is the flat Gauss sum `±√p`, which maximizes `H` (THM-162) and (via Hoffman, HYP-2311)
   pins `χ_di = 2`. The "out-of-phase" imaginary spectrum (p≡3) is the skew tournament. Paley is where
   the character-ratio spectrum (S636), the trienement (S631), and reciprocity coincide.
2. **The 2-adic seam is reciprocity's parity.** real↔imaginary = symmetric↔skew = graph↔tournament =
   `p≡1`↔`p≡3` = even↔odd. Every "even is hard" result (glass transition, parity defect, LRC `n=14`)
   is this seam; reciprocity is its arithmetic.
3. **CM fields (the χ-bridge / grid disproof, HYP-2301/2230).** Split/inert/ramified of `p` in
   `ℚ(√d)` is the Legendre symbol `(d/p)` = reciprocity; the modulus-1 elements and unit-distance
   Cayley graphs are governed by which primes split. The out-of-phase `i` is `ℚ(i)` = the CM rotation.
4. **XNOR = the Walsh/character basis.** Walsh characters are products of `±1` bits = iterated XNOR;
   `H`'s even-degree Walsh structure (THM-163) and the polarized delta field (HYP-2316) live in this
   basis. XNOR is the atom; the quadratic character is its arithmetic lift; reciprocity is its law.
5. **Paley metastability (THM-135).** Paley is a local `H`-max but NOT global at `p=19` — a
   **metastable state** of the antiferromagnetic landscape (HYP-2316). So reciprocity's extremal
   object sits in the glass's rugged phase: the seam and the glass are one.

## The abstraction
**XNOR (`±1` agreement) is the atom; the quadratic character is its arithmetic incarnation with an
out-of-phase `0`-middle (the conference matrix / Paley trienement); quadratic reciprocity is the law
coupling the 2-adic-seam parities of two primes; and that seam — real/imaginary, symmetric/skew,
even/odd — is the single fault line under the LRC even-`n` frontier, the glass transition, the parity
defect, and the CM/character spectrum.**

## To do
1. Read the LRC `2n-1` shell through reciprocity: when `2n-1=p` prime, the multiplier dodge
   ((ℤ/p)*, S625) is the Legendre/QR structure; is the unramified-vs-ramified (n=19 vs n=14) the
   `(·/p)` split governed by reciprocity?
2. Paley as the H-landscape's metastable seam: chart Paley's basin (THM-135) vs the glass (HYP-2316)
   across `p≡3 mod 4`; does it freeze (go metastable) at a reciprocity-determined `p`?
3. Make "XNOR atom → quadratic character → reciprocity" a clean tower (Walsh degree 1 → quadratic
   character → the reciprocity 2-cocycle).
