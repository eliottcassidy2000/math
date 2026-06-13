# Reed–Muller on the Tiling Cube: the blue code, the hypotenuse defect, and the digit ladder

**Source:** mac-mini-2026-06-11-S1. Dispatch: RM(r,m)^⊥ = RM(m−r−1,m) — repetition↔SPC,
biorthogonal↔extended-Hamming, k = n/2 self-dual — as an analogy engine for tournaments.
Canon: THM-477, THM-478, HYP-2406..2408, T779. Builds on THM-474 (tilings = switching
classes), THM-163, Grinberg–Stanley 2023, kind-pasteur's THM-466 digit tower.

## The dictionary that emerged

| Reed–Muller world | tournament/tiling world |
|---|---|
| ambient cube F_2^{2^m} | tiling cube F_2^m, m = C(n−1,2) (= switching classes, THM-474) |
| RM(0,m) repetition code | Rédei: digit_0(H) ≡ 1 — H mod 2 IS the repetition codeword |
| RM(D,m), low-degree stratum | digit_1(H) = c_odd mod 2 has degree exactly D = 2⌊(n−1)/2⌋ (THM-478) |
| RM(r)^⊥ = RM(m−r−1): min-weight duals are flats | XOR of H mod 4 over every (D+1)-flat = 0 — an exact recurrence on H |
| k = n/2 self-dual codes | the blue code B_n = ker(1+σ): orbit-locally self-dual [2,1] repetition; defect = hypotenuse |
| Plotkin (u, u+v) doubling | Mode-B leg-swap: b(n) = b(n−2) + (n−2) |
| biorthogonal RM(1,m) = Sylvester/Hadamard | the skew tower (THM-447/451), Sylvester through order 8 |
| cut code ⊥ cycle code of K_n | switching (cut cosets) vs even graphs (cycle space) — the old GF(2) duality |

## Three findings

1. **The blue code (THM-477).** Grid-symmetric tilings are not just a symmetric
   subset — they are a LINEAR, DUAL-CONTAINING code: B_n^⊥ = im(1+σ) ⊆ B_n with
   self-dual defect exactly f = ⌊(n−1)/2⌋ = the hypotenuse. The user's instinct
   ("k = n/2 self-duality ≈ blue line pairs") is exact: on every σ-paired tile the
   blue constraint IS the unique self-dual code of length 2, and the entire failure
   of global self-duality is carried by the staircase's hypotenuse — the same side
   canon already credits with H = 1+2^d, the fiber fraction, and Walsh order-2.
   The glue group B_n/B_n^⊥ ≅ F_2^f, and the complement translation (blue LINES)
   acts on it by the all-ones glue vector.

2. **The digit ladder (THM-478, HYP-2406/2407).** The 2-adic digits of H, as
   Boolean functions on the tiling cube, climb the RM hierarchy from the bottom:
   digit_0 ∈ RM(0,m) (Rédei), digit_1 ∈ RM(D,m) with degree EXACTLY D = the
   THM-163 band limit (proved ≤ via the cycle-reversal pairing — each odd cycle's
   two orientations cancel their top term, the same cancellation that makes a
   cyclic triple's parity quadratic, i.e. the Babai–Cameron two-graph cocycle);
   digits ≥ 3 saturate to (near-)full degree m. **The 2-adic depth of H is graded
   by Reed–Muller level: arithmetic depth = polynomial degree, structured for two
   layers, then pseudorandom.** Combined with kind-pasteur's THM-466 (digits =
   odd-cycle collection counts): collection counts at depth 1 are low-degree,
   collection counts at depth ≥ 3 are not — the OCF has a shallow algebraic shell
   and a hard core, quantitatively.

3. **The recurrences asked for.** deg(digit_1) = D is EQUIVALENT to: the XOR of
   H mod 4 over every (D+1)-dimensional affine flat of the tiling cube vanishes
   (flat indicators = minimum-weight codewords of the dual RM code). At n=7 this
   says: on every 7-dim subcube of the 15-cube of tilings, H mod 4 at one corner is
   determined by the other 127. The blue recursion b(n) = b(n−2) + (n−2) is the
   (u,u+v) Plotkin step in Mode-B time.

## Where this points

- **σ-invariance is total**: every digit of H is grid-transpose invariant; none is
  complement-translation invariant (MISTAKE-033/034 live here). So the digit
  functions descend to the σ-orbit space, and their restrictions to B_n (the blue
  slice) are natural next objects.
- **The d_2 sequence 3, 6, 7, 11 (n = 4..7)** has no law yet — the first
  unexplained number sequence of this thread (HYP-2407).
- **Ax shadow (HYP-2408)**: degree-≤D level sets have forced 2-adic divisibility
  (Ax–Katz); forbidden H-values are joint digit-level-set emptiness — a possible
  second, coding-theoretic route to H ∉ {7,21}.
- The biorthogonal row of the dictionary (RM(1,m) = Sylvester = the skew tower,
  which exits the RM family at order 16, THM-451) suggests the dual question:
  what tournament object plays extended Hamming RM(m−2,m)? The dual of the skew
  tower code is unexamined.
