# Independent audit of THM-4761, HYP-9370 and HYP-9369 §§3–5 (opus-2026-10-22-S54)

The auditor was an independent agent. It used its own code (scripts and outputs in this folder) and edited and committed nothing. The audit tool could not write this file, so the lane saved the report verbatim from the auditor's final message.

**Audited text:** HEAD `f949b0c234`. The three documents are unchanged since `62e239388d`, except HYP-9369 §2, which was out of scope. The lane rewrote the `deficit_*.txt` files during the audit (header plus self-converse rows only), and the auditor re-checked against the new files.

## Verdicts

| document | verdict | main issues |
|---|---|---|
| **THM-4761** | **CONFIRMED with one MEDIUM correction.** Every number and the asymptotics check out exactly. | Display (1) is missing a sign and is false as written. Score-scale wording (LOW). Prior-art additions (LOW–MEDIUM). |
| **HYP-9370** | **Data CONFIRMED; the text needs corrections.** No counterexample to (a) or (b), including new structured families. | The "min t₀" values at n = 12, 14, 16 are sample minima, and structured tournaments go lower (MEDIUM). The mechanism paragraph has false statements (MEDIUM). The analytic condition `a < r` is wrong (MEDIUM-LOW). Several LOW fixes. |
| **HYP-9369 §§3–5** | **CONFIRMED with corrections.** Square identity correct; claw identity independently re-verified to 5e-15; all deficit numbers reproduced. | The claw formula's `Q` is attributed to Session30's notation, which makes it false (MEDIUM). "n³X(−1,3) → 0" contradicts the cited data (LOW-MEDIUM). The |Aut| bound is Dixon's, not Moon's (LOW). |

## A. THM-4761

**(i) Chaos argument: CORRECTED (MEDIUM).**
- Orthogonality is fine: distinct edge sets give distinct Walsh characters.
- The display `N_L = Σ_F 2^p(F) χ_F` with `χ_F = ∏_{ij∈F, i<j} s_ij` is false. A directed path multiplies `±s_ij`, picking up `−1` for each arc that goes from a larger to a smaller label.
- The correct form is `N_L = Σ_F 2^p(F) ε_F χ_F` with `ε_F = (−1)^{#descending arcs}`. This is well defined because every path has even length.
- `chaos_check.out` §3: the unsigned formula fails on 15 of 16 tournaments with n = 3..6. For example, transitive T₃ has N₂ = −2, while the formula gives 6. The signed formula matches all 16.
- Since `ε_F² = 1`, statements (2)–(4) are unaffected.

**(ii) Exact second moment: CONFIRMED.**
- Brute force over all 2^C(n,2) tournaments gives `E κ² = 4/3, 4/3, 79/60, 58/45, 635/504` for n = 3..7. This equals formula (2), evaluated literally over multisets.
- Independent proof:
  - `E κ² = (1/n!) Σ_{τ∈S_n, no step a→a−1} 2^{#steps a→a+1}`.
  - Inclusion–exclusion over even intervals gives formula (2), term by term.
  - Brute force of this sum agrees for n ≤ 10.
- Equivalent closed form: `E N_L² = Σ_p 2^p (n)_{L+p} C(L/2−1, p−1)/p!`.

**(iii) Asymptotics: CONFIRMED.**
- A sympy series of the shapes gives `2/n − (14/3)n⁻³`, and the share is `1 − 1/n + (4/3)n⁻²`.
- `n·Var κ` = 1.871429, 1.952858, 1.985485, 1.994017, 1.997968.
- `n³(Var κ − 2/n)` → −4.666860 at n = 10⁵.
- `n²(share − 1 + 1/n)` → 1.333397.
- `n² E(κ−1−κ₂)²` → 1.99997.
- The n⁻⁴ coefficient of Var κ is −58/3.

**(iv) N₂ = n(n−1) − ‖S1‖²: CONFIRMED** algebraically and on 80 random tournaments.

**(v) Section 2: CONFIRMED in order of magnitude, with LOW wording fixes.**
- The table and the −0.70/−0.53 values reproduce from the lane's files, with H matching the audit DP on 1500/1500 and 800/800 samples.
- On 20,000 fresh samples at n = 12 and 3,000 at n = 14:

| quantity | n = 12 | n = 14 |
|---|---|---|
| corr(log κ, log ρ + T(1)) | 0.964 | 0.973 |
| corr(log κ, −W) | 0.962 | 0.971 |
| corr(log κ, log ρ) | −0.04 | −0.06 |
| mean of log(κ/ρ) − T(1) | +0.049 | +0.048 |

- At n = 12, corr(κ, 1+κ₂) = 0.9654 against √share = 0.9649.
- Wording: "`‖S1‖_∞ ≍ √n`" should be "typical `|(S1)_v| ≍ √n`"; the sup norm is ≍ √(n log n).
- Wording: the "second-order classes" added are X(1,0) + E₀[a]/n plus the lollipop, X(0,2) and the star. That is not Session30's second-order predictor.

PROVED is justified once the sign is fixed.

## B. HYP-9370

**Small exact cases: all confirmed** (`zeros_small.out`), every nearest zero negative real.
- C₃: P = 1 + (2/3)w, t₀ = √(3/2).
- Almost-regular 4: P = 1 + w/2, t₀ = √2.
- C₅: P = 1 + (4/5)w, t₀ = √5/2.
- Regular 7, all 2640 labelled, 3 classes: P₇ 1.1311, C₇(1,2,3) 1.1414, third class 1.1445.
- Almost-regular 6: min 1.2335, median 1.2735.
- Almost-regular 8: min 1.2222, median 1.2603, max 1.2950. The first-class coefficients match the text.

**Class-list completeness: CONFIRMED** (`census_check.out`).
- Brute-force canonical forms show the n = 4, 6, 8, 10 lists are pairwise distinct.
- Every |Aut| matches the file.
- The masses equal the audit's labelled counts: 24, 2640, 3,230,080, 48,251,508,480.
- Self-converse counts are 1, 3, 15, 229.

**Larger runs: all confirmed** (`zb_*.out`).
- n = 10, all classes: min 1.2467, all nearest zeros negative real.
- Regular 11 completions: min 1.2194. Regular 9: 15 classes, min 1.1790.
- Samples at 12, 14, 16 reproduce the stated min/median/max.
- DRT(23), recomputed from mac-mini's forward-step data: min 1.49384, median 1.49390.

**New structured families: no counterexample** (`zeros_structured.out`). Every nearest zero is negative real.
- All circulant regular tournaments on 13, 15 and 17 vertices: min 1.2763, 1.3100, 1.3370.
- Their vertex deletions (almost-regular 12, 14, 16): min 1.2789, 1.3095, 1.3357.
- Rotational R₁₉ = 1.3589; R₁₉ − v = 1.3576.
- R[K₂] for all regular 5- and 7-tournaments, C₃[C₅], C₅[C₃] and their deletions, C₃[C₃[K₂]] and C₉[K₂]: all ≥ 1.31.

**Finding (MEDIUM): the stated minima at n = 12, 14, 16 are not family minima.**
- The deleted rotational tournament R_{n+1} − v is the exact minimiser of the complete censuses at n = 6, 8, 10.
- At n = 12, 14, 16 it gives 1.2789, 1.3095, 1.3357. These are below the title and table "min" values 1.300, 1.349, 1.396, which are MCMC sample minima. The chain starts at exactly this tournament.
- Claim (c) survives with the corrected values: 1.2222 < 1.2467 < 1.2789 < 1.3095 < 1.3357 < 1.3576.
- "Regular orders sit below the neighbouring almost-regular ones" holds only for n ≤ 11. From q = 15 on, R_q − v lies below R_q.
- DRT(23) is not a family minimiser, so "on the trend" is loose.

**Random tournaments: zeros do enter the disc, but less often than stated.**
- Lane's file: 5/1500, min 0.9446.
- 20,000 random 12-tournaments: 21 (0.105%), min 0.8535.
- 3,000 random 14-tournaments: 1, min 0.9026.
- Every inside-disc zero is real (w > 0). These cases have very negative N₂ (W > n(n−1)), which gives a positive root of P_T.
- The kinds are 42.0% / 38.8% / 19.2%; the text says 43% / 38%.

**Mechanism paragraph: not sound (MEDIUM)** (`szego_check.out`).
- "P_T(w) with w > 0 is a positive sum" is false. Negative coefficients and positive real roots occur in 22/85 almost-regular 8 classes and 5379/13,333 almost-regular 10 classes, among others.
- "P_T(−t²) alternates, and its first sign change is the nearest zero" is not a valid argument.
- The Szegő model gives complex nearest zeros at even w-degree (n = 10, 14; arg w ≈ 2.66–2.80). Every exact P_T has a negative real nearest zero. So (b) is not explained by the model.
- The model numbers 1.351–1.716 are confirmed. The 0.95–0.79 ratios are measured in w, not |z|.

**"What this gives": partly unsound.**
- MEDIUM-LOW: "`C n⁻³ a^L` with `a < r`" must be "`C n⁻³ a^{−L}` with 1 < a < r".
- LOW: (a) guarantees an asymptotic ratio of 10/11 in general and 5/6 for almost-regular n ≥ 8. The text uses 1/1.22 and, in the title, 1/1.118; both are observed values, not guarantees.
- LOW: a zero on |z| = 1 does not force divergence at z = 1; only the open disc does.
- LOW: the spectrum provably cannot bind below √2 (Session30's bound gives 1/max|x_j| ≥ √(2n/(n−1))), so this can be stated rather than observed.
- LOW: the regular-table medians are taken over completions. Per-class medians are 1.1414 (n = 7) and 1.1876 (n = 9).

## C. HYP-9369 §4 and §3

**Square identity: CONFIRMED.** G = −Ŝ²(I+Ŝ²)⁻¹ is symmetric, so −(d+1)ᵀG(d+1) is exactly the dumbbell plus lollipop plus path-tree closed forms.

**Claw predicate: CONFIRMED.** With 4 vertices the edge set is forced, so `is_claw` picks exactly one core in each of (2,0), (1,1), (0,2), (−1,3). The audit's own predicate on Session30's core enumeration selects the same four cores.

**Claw identity: independently confirmed** (`claw_check.out`).
- The lane's brute-force self-test covers only L ≤ 6. The (1,1) and (2,0) claw cores first appear at L = 8 and L = 10, so they were not covered.
- All cores were recomputed with Session30's engine on nine tournaments, including uniform random and transitive ones.
- Claw max difference: 5.3e-15. Two-arm max difference: 1.1e-15. The identities appear to hold for all tournaments, not only almost-regular ones.

**Notation error (MEDIUM).**
- Y(u) holds with `Q̃ = Ŝ(I+Ŝ²)⁻¹`.
- The text says "Session30's notation (G, Q)", but Session30's Q is zŜ. With Q = Ŝ the identity fails by 2e-4 up to 0.70.
- "1e-17" is an absolute difference; relative agreement is about 1e-15.

**§3 checks.**
- RMSE 0.210 vs 0.048 is confirmed, as are the weights 13.0 and 52.0 and the φ₆ comparison.
- "n³X(−1,3) → 0" is wrong. `class_scaling.out` shows n³X ≈ −2.2; the intended statement is n²X(−1,3) → 0.
- The truncation offsets are +8.34, −8.92, +5.15, not "+5.1 to −8.9".

## D. HYP-9369 §5

**(i) Recomputation: CONFIRMED** (`deficit_check.out`; own DFS, blue test and |Aut| backtracking).

| tournament | (h, b) | Δ |
|---|---|---|
| C₅ | (3, 3) | 0 |
| C₃[K₂] | (15, 9) | 6 |
| P₇ | (9, 3) | 6 |
| C₃[C₃] | (39, 11) | 28 |

- Almost-regular 6 and regular 7: tables confirmed.
- Every row of `deficit_ar8.txt` and `deficit_reg9.txt` matches, with 0 mismatches in both file versions.
- n = 10 and 11 file statistics match (229 self-converse, min Δ 1592/1692, max b 193); spot-checked rows agree.

**(ii) Bound b ≤ 2^⌊n/2⌋⌊n/2⌋!: CONFIRMED for both parities.**
- An anti-involution has no fixed point for even n and one for odd n.
- The first half of a blue word comes from distinct orbits.
- The number of anti-involutions is at most |Aut|.

**(iii) Asymptotic argument: CONFIRMED.**
- Cuckler's lower bound plus the completion injection is valid for almost-regular T.
- The constant is 48 = 4·3·2·2.
- Attribution: |Aut| ≤ √3^{n−1} is Dixon (1967), Canad. Math. Bull. 10, 503–505. It is Exercise 5 of §28 in Moon's book; Moon's own theorem there gives only 2.5ⁿ.

**(iv) Unverified closure: correctly labelled.** The values 520,738 / 1,039,055 / 3,848,613 appear only as user-reported in Session30's notes. The onset-19 algebra re-derives exactly with sympy: same numerator, denominator p²u², gcd 1, residual −1 at k = 18 and 0 from 19 on, same g₅..g₁₆.

## E. Prior art

- **Janson (1995)**, RSA 6, doi:10.1002/rsa.3240060209. The abstract, seen only through a search snippet of the publisher listing (page 403), says the Hamilton-cycle count is asymptotically normal. Full text not read.
- **Wormald**, arXiv:2610.16245, says only that "routine calculations" give a variance of o(g(n)²).
- **Moon, *Topics on Tournaments*** (Gutenberg #42833, TeX source read):
  - §10 has Szele's mean.
  - §6, Theorem 6: the number of 3-cycles c₃ has variance (3/16)C(n,3) and is asymptotically normal (Kendall–Babington Smith 1940 conjectured it; Moran 1947 proved it).
  - There is no variance for spanning paths.
- **Why Moon's Theorem 6 matters:**
  - W = n(n²−1)/3 − 8c₃, verified on all samples. So κ₂ = 8c₃/(n(n−1)) − (n−2)/3.
  - THM-4761's Var κ₂ is therefore the classical variance.
  - Combined with E(κ−1−κ₂)² = O(n⁻²), this gives √n(κ−1) ⇒ N(0, 2) immediately.
  - The exact formula (2) and the −14/3 coefficient were not found anywhere.
- **Zeros of the forest polynomial:** no prior work found on zeros of linear-forest or path-cover polynomials of tournaments. The nearest objects are:
  - the Rédei–Berge function (arXiv:2307.05569, 2402.07606);
  - D'Antona–Munarini's cycle-path indicator polynomial (Adv. Appl. Math. 2000).
  - Szegő (1924) should be cited for the partial-sum zeros.
- **Decoration rule:** no identical statement found. It is analogous to tadpole dressing of an external source and to the articulation-point "topological reduction" in Mayer expansions. Present it as an analogy, not prior art.

## Required edits (ranked by severity)

1. **MEDIUM — THM-4761 (1).** Replace the unsigned display with "`N_L = Σ_F 2^(p(F)) ε_F χ_F`, `ε_F = (−1)^(#arcs traversed from larger to smaller label)`, well defined since paths have even length". Add "Since `ε_F² = 1`, `E N_L² = Σ_F 4^p(F)` is unchanged."
2. **MEDIUM — HYP-9369 §4.** Replace "With Session30's notation (`G`, `Q`)" with "With Session30's `G` and `Q̃ = Ŝ(I + Ŝ²)^(−1)` (Session30's `h = Q̃1`)", and write `(Q̃u)_w` in Y(u). Replace "better than `1e-17`" with "to machine precision (relative ~1e-15; an independent audit with Session30's engine on nine tournaments, including non-almost-regular ones, gives max |diff| 5e-15)".
3. **MEDIUM — HYP-9370 title, table and (c).** Replace "almost-regular 1.222, 1.247, 1.300, 1.349, 1.396 at n = 8..16" with "almost-regular 1.222, 1.247 at n = 8, 10 (complete) and ≤ 1.279, 1.310, 1.336 at n = 12, 14, 16, attained by R_(n+1) − v, the exact minimiser at n = 6, 8, 10". Relabel the n = 12–16 table entries as "sample min".
4. **MEDIUM — HYP-9370 Mechanism.** Delete "Along the real z-axis … positive sum. Along the imaginary axis … its first sign change is the nearest zero." Replace it with:
   - Negative high coefficients and positive real roots occur (22/85 classes at n = 8; 5379/13,333 at n = 10), but never as the nearest zero.
   - At even w-degree the truncated models have complex nearest zeros (arg w ≈ 2.66–2.80), so (b) is not explained by the Szegő picture.

   In the title, replace "The zeros are deformed Szego zeros" with "comparable in modulus (0.79–0.95 in w) to the truncated one-loop exponential".
5. **MEDIUM-LOW — HYP-9370.** Replace "`≤ C n^(−3) a^L` with `a < r`" with "`≤ C n^(−3) a^(−L)` with `1 < a < r`".
6. **LOW-MEDIUM — HYP-9369 §3.** Replace "`n³X(−1,3) → 0`" with "`n²X(−1,3) → 0` (`n³X(−1,3) ≈ −2.2`)". Replace "+5.1 to −8.9" with "+8.3 (T0), −8.9 (T1), +5.1 (T2)".
7. **LOW-MEDIUM — THM-4761 prior art.** Add:
   - κ₂ = 8c₃/(n(n−1)) − (n−2)/3, so Var κ₂ is Kendall–Babington Smith's (3/16)C(n,3) (Moon §6, Theorem 6).
   - Moran's CLT gives √n(κ−1) ⇒ N(0, 2).
   - Janson's (1995) abstract: the Hamilton-cycle count is asymptotically normal.

   Then revise the Janson sentence accordingly.
8. **LOW — HYP-9369 §5.** Replace "(Moon)" with "(Dixon 1967; Moon, Topics on Tournaments, §28 Ex. 5)".
9. **LOW — HYP-9370, random tournaments and circulants.**
   - Replace "0.33%" in the text and "frequency 0.003" in the status line with "≈0.1–0.3% (21/20,000 in the audit sample, min 0.854; also 1/3000 at n = 14); all such zeros real".
   - Replace "43% … 38%" with "42% … 39%".
   - Replace "Circulant families beyond n = 11 were not tested" with the circulant result (no zero in |z| ≤ 1.27 for n = 13–17).
10. **LOW — HYP-9370 consistency.**
    - Make the ratios consistent: 10/11 in general and 5/6 for almost-regular n ≥ 8, given (a).
    - Change "closed unit disc" to "open unit disc".
    - Label the regular-table medians as taken over completions.
    - Restrict "regular below almost-regular" to n ≤ 11.
11. **LOW — THM-4761 §2.** Replace "`‖S1‖_∞ ≍ √n`" with "typical `|(S1)_v| ≍ √n`". Name the classes actually added instead of "the almost-regular second-order classes".

## Main scripts and outputs (this folder)

- `forest_audit.c`, `canon_brute.c`, `ekappa2_brute.c`, `tourn.py`
- `chaos_check.out`, `census_check.out`, `zeros_small.out`, `zb_*.out`, `zeros_structured.out`, `rotational_check.out`, `szego_check.out`
- `drt23_check.out`, `rand_reading.out`, `claw_check.out`, `deficit_check.out`, `onset_check.out`, `misc_check.out`

Bulky intermediates (`*.rows`, `canon_ar10.out`) and compiled binaries are not committed.

## Sources

- [Janson 1995 (Wiley)](https://onlinelibrary.wiley.com/doi/abs/10.1002/rsa.3240060209)
- [Wormald, Tournaments with many Hamilton cycles](https://arxiv.org/html/2610.16245)
- [Moon, Topics on Tournaments (Gutenberg)](https://www.gutenberg.org/ebooks/42833)
- [Cuckler 2007 (CPC)](https://www.cambridge.org/core/journals/combinatorics-probability-and-computing/article/abs/hamiltonian-cycles-in-regular-tournaments/5CB605A35014A30F60F740D74199D3E7)
- [Grinberg–Stanley](https://arxiv.org/pdf/2307.05569)
- [Rédei–Berge Hopf algebra](https://arxiv.org/pdf/2402.07606)
- [Cycle-path indicator polynomial](https://www.sciencedirect.com/science/article/pii/S0196885800906894)
