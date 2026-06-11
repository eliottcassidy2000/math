# THM-475 — The DRT flag construction attains the tournament Barba value 2(n−1)^((n−1)/2)

**Status:** PROVED (claudebox-2026-06-11-S1) + VERIFIED exact (n = 9, 13, 17, 25, 29).
**Provenance:** claudebox-2026-06-11-S1. Proves the lower-bound (construction) half of
HYP-2389 / OPEN-Q-058 (the tournament Barba problem, mac-mini-2026-06-10-S2).
**Companions:** THM-472 (ceiling — proves the upper context), THM-468 (d(T) invariant),
THM-447/448 (skew-Sylvester doubling / DRT Mersenne tower), HYP-2405 (even sibling).
**External anchors:** Tao optimization-constants repo C23a (smallest open Hadamard order 668);
the n ≡ 1 mod 4 problem is absent from the literature (Cameron Problem 104 thread, Armario's
skew papers, and Klanderman–Montee–Piotrowski–Rice–Shader LAA 2024 all treat det(S), which
vanishes for odd n, not det(I+S); OEIS negative for 1,2,4,16,32,160,512,…).

## Statement

Let m ≡ 3 (mod 4) and let S₀ be the skew-adjacency matrix of any doubly regular tournament
(DRT) on m vertices (S₀² = J − mI, S₀𝟙 = 0). Put n = m + 2 ≡ 1 (mod 4) and define the
**flag tournament** F = Flag(S₀) on n vertices: two new vertices u, v with

- u → x and v → x for every DRT vertex x (both apexes beat the whole DRT),
- u → v.

Then the skew matrix S of F satisfies

1. SSᵀ has spectrum { m with multiplicity m−1, 2m+1 with multiplicity 2, 0 once },
   i.e. char poly of S is x·(x² + n−2)^((n−3)/2)·(x² + 2n−3);
2. **det(I + S) = (1+m)^((m−1)/2)·(1+(2m+1)) = 2(n−1)^((n−1)/2)** — exactly the HYP-2389
   conjectured maximum, with exactly the HYP-2389 maximizer spectrum (flat base level n−2,
   one excited pair at 2n−3).

In particular B_t(n) := max det(I+S) ≥ 2(n−1)^((n−1)/2) for every n ≡ 1 (mod 4) such that a
DRT on n−2 exists. By THM-448 / Paley / GF(27) this is unconditional for
n−2 ∈ {q : q ≡ 3 mod 4 prime power} ∪ {2^k − 1 doubling-tower orders}; under the
skew-Hadamard conjecture it holds for ALL n ≡ 1 (mod 4).

## Proof

Write blocks in the order (DRT, u, v). The flag's skew matrix is

```
S = [ S₀   −𝟙   −𝟙 ]
    [ 𝟙ᵀ    0   +1 ]
    [ 𝟙ᵀ   −1    0 ]
```

**Gram blocks.** Using the two DRT identities S₀S₀ᵀ = −S₀² = mI − J and S₀𝟙 = 0:

- (SSᵀ)_{DRT,DRT} = S₀S₀ᵀ + 𝟙𝟙ᵀ + 𝟙𝟙ᵀ = mI + J;
- (SSᵀ)_{i,u} = (S₀𝟙)_i + 0 − 1 = −1; (SSᵀ)_{i,v} = (S₀𝟙)_i + 1 + 0 = +1;
- (SSᵀ)_{u,u} = (SSᵀ)_{v,v} = m + 1; (SSᵀ)_{u,v} = m + 0 + 0 = m.

**Spectrum.** For x ⊥ 𝟙 supported on the DRT block, SSᵀ(x,0,0) = (mx, −𝟙ᵀx, 𝟙ᵀx) = m·(x,0,0):
eigenvalue m = n−2 with multiplicity m−1. The complementary invariant 3-space
span{(𝟙,0,0), e_u, e_v} carries the matrix (columns = images)

```
B = [ 2m   −1    1  ]
    [ −m   m+1   m  ]
    [  m    m   m+1 ]
```

Row-reduce R₂ ← R₂ + R₃ in det(B − λI): the new row is (0, 2m+1−λ, 2m+1−λ), so

det(B − λI) = (2m+1−λ)·[(2m−λ)(1−λ) − 2m] = (2m+1−λ)·(λ² − (2m+1)λ) = −λ(λ−(2m+1))².

Eigenvalues {2m+1, 2m+1, 0} = {2n−3, 2n−3, 0}. Together: SSᵀ ~ {m^(m−1), (2m+1)², 0}. ∎(1)

det(I+S)² = det(I + SSᵀ) = (1+m)^(m−1)·(2m+2)²·1 = 4(n−1)^(n−1), and det(I+S) > 0
(THM-468: a sum of squares of odd Pfaffians), so det(I+S) = 2(n−1)^((n−1)/2). ∎(2)

## Verification (exact, fraction-free Bareiss over ℤ; script tournament_barba_flag_cbx1.py)

| n  | DRT(n−2) used                  | det(I+S) flag           | 2(n−1)^((n−1)/2)        | match |
|----|--------------------------------|-------------------------|-------------------------|-------|
| 9  | Paley QR₇                      | 8192                    | 2·8⁴                    | ✓     |
| 13 | Paley QR₁₁                     | 5971968                 | 2·12⁶                   | ✓     |
| 17 | doubling-tower DRT(15) (THM-448 core of skew-H(16)) | 8589934592 | 2·16⁸      | ✓     |
| 25 | Paley QR₂₃                     | 73040694872113152       | 2·24¹²                  | ✓     |
| 29 | GF(27)-Paley (F₃[x]/(x³−x−1))  | 364118239659885068288   | 2·28¹⁴                  | ✓     |

Cross-check against the mac-mini exhaustive censuses: at n = 9 the flag's S-char-poly
x⁹+36x⁷+462x⁵+2548x³+5145x = x(x²+7)³(x²+15) is EXACTLY the unique char poly shared by all
216 maximizer classes (hadamard_det_n9_macmini_s2.out); at n = 5 the construction degenerates
to C₃+flag and gives 32 = the exhaustive maximum. So at n = 5 and 9 the flag value is the TRUE
maximum; OPEN-Q-058 (upper bound for all n ≡ 1 mod 4) remains open.

## Remarks

1. **Why two apexes.** One borrowed eigenvalue needs a 2-dimensional carrier. The DRT supplies
   the flat level n−2 with its kernel 𝟙; the flag pair {u,v} couples to the DRT only through 𝟙,
   and the 3×3 block above converts {kernel + two apex modes} into {excited pair 2n−3, new
   kernel}. The integrality obstruction of THM-472 (no DRT at n ≡ 1 mod 4) materializes as
   exactly one excited pair — the minimal spectral defect, and the flag realizes it.
2. **Apex narrative.** The two flag vertices are literally apexes (dominant vertices); the
   tournament Barba maximizer is "a DRT watched by two stacked observers". The +1/observer/apex
   motif of the perspective-key arc (HYP-2130 …) recurs as the extremal mechanism.
3. **Doubling corollary (one line from THM-447's spectral law).** det(I + S_{D(T)}) =
   2^n·det(I+S_T)²: the skew-Sylvester double preserves the Hadamard ceiling exactly
   (skew-H ↦ skew-H, both sides n^(n/2) ↦ (2n)^n), but is strictly suboptimal at the other
   residues (e.g. n=5 max 32 ↦ 32768 < B_t-candidate 64000 at n=10): the maxdet ladder is NOT
   doubling-closed off the Hadamard residue.
4. **Tao C23b hook.** κ(I+S)² = (1+λ_max(SSᵀ))/(1+λ_min(SSᵀ)); the flag maximizer has
   κ = √(2n−2) (odd order forces σ_min = 1) — odd-order skew matrices are bad approximate
   Hadamard matrices; the conditioning frontier (17/92 exponent, AJM 2025 arXiv:2511.14653)
   lives at even orders, where THM-476's square law rules.
