# Audit B: number theory (THM-4566, THM-4567, the 2026-10-07 updates, and §4–5 of the note)

Auditor: independent adversarial audit, 2026-10-07. Repository read-only. Scripts and outputs are in this directory (listed at the end). openai/math results are taken as correct, per the owner directive.

## Verdict summary

| Object | Verdict |
|---|---|
| THM-4566 identification (genus theory, m ≡ 1 mod 4, F = Q(√p : p∣m)) | CONFIRMED |
| THM-4566 list of 19 m and completeness modulo #003 | CONFIRMED |
| THM-4566 χ table | CONFIRMED WITH CORRECTION. No stated inequality is false, but one sentence is false (no conjugation-fixed prime below 60, for m = 105 and 357), and the lower bounds for m = 177, 253, 345 and 357 are improvable using the repo's own tools: m = 177 is exactly 4 |
| THM-4567 (all five statements) | CONFIRMED |
| THM-4555 update (gap lemma, the 8207 core, the switch 53803/26901, the census 27/12/14) | CONFIRMED |
| THM-4512 alternative closure (Ellison) | CONFIRMED WITH CORRECTION (one constant) |
| THM-4558 pointer | CONFIRMED (should carry the m = 177 value) |
| Note §4–5 | Mostly CONFIRMED. Three wording corrections are needed (Borwein–Choi, Serre, Fano gadget) |

---

## 1. THM-4566

### 1a. Identification: CONFIRMED

* Gal(H/Q) = Cl(K) ⋊ ⟨c⟩, with c acting by inversion. H is abelian over Q iff exponent ≤ 2 iff H equals the genus field.
* **A missing step should be added.** If H = F(i) with F real and Galois over Q, then F = H ∩ R and ⟨c⟩ = Gal(H/F) is normal of order 2. So c is central, and the exponent is ≤ 2. The text only states the abelian ⇔ exponent ≤ 2 direction.
* i lies in the genus field iff −1 is, modulo squares, a product of the prime discriminants. The only way is the factor −4 itself, so d_K = −4m with m ≡ 1 mod 4. The genus field is then Q(i, √p : p∣m).
* F = H ∩ R is unique, and F determines m, so there are no duplicate planes.
* h(−4m) = 2^ω(m) = [genus field : K] for all 19 m (`idoneal_check.out`).

### 1b. List and completeness: CONFIRMED

* **Independent enumeration** (`exp2_audit.c`): every N ≡ 0, 3 (mod 4) with N ≤ 2.1·10¹¹, fundamental or not, was tested for a primitive reduced form with 0 < b < a < c.
  * Method: direct test to 10⁸; above that, a wheel sieve on split primes p ≤ 400 (with 4p² < N), and survivors re-checked in full.
  * Run time 46 s. Result: exactly **101** values, identical to the reader's list and to `b003171.txt`.
  * The 4∣N entries give exactly Euler's 65 idoneal numbers. Of these, exactly 19 are squarefree and ≡ 1 mod 4: the stated list.
* **Completeness chain, checked against the sources in `scratchpad/oai3/nt/`.** EKN Thm 4 assumes L(s, χ) ≠ 0 on [1 − 1/(4 log|D|), 1), and Tatuzawa's exceptional zero (EKN Lemma 2) lies above 1 − ε/4 > 7/8; #003 (Re s > 7/8) excludes both.
  * #003's 11/12 paper explicitly states the idoneal consequence (`qrh1112.txt`, lines 90–94).
  * The non-fundamental cases are covered by the conductor bound f ∣ 840, which I re-derived. My sieve also covers every f²D₀ ≤ 2.1·10¹¹.
* **Precision.** EKN's threshold is d₁₁ = 200560490130 ≈ 2.006·10¹¹, not "2·10¹¹". This is harmless, since the sieve covers the gap.
* The reader's 2.34·10¹⁴ sieve was not re-run; its code (`exp2_sieve.c`) is sound and it is not needed for completeness.

### 1c. Chromatic table: CONFIRMED WITH CORRECTION

**Decomposition groups (`lemmaR_audit.py`).** I recomputed c ∈ D and c ∈ I for every prime q < 200. Every stated upper bound and its prime is confirmed.

**κ values (SAT, `kappa_check.out`).**
* κ(3) = 3, κ(7) = 4 and κ(11) = 5 are exact.
* κ(19) ≤ 5: a 5-colouring was found, which confirms the m = 385 upper bound independently of decalion89.
* κ(2) = 4 is K₄.

**Rows confirmed.**
* **m = 1, 5, 13, 37, 85: χ = 2** (2 ramified, c ∈ I; no odd closed walk is forced in Q(√5)², Q(√13)², Q(√17)², Q(√37)²).
* **m = 21, 57, 93, 273: χ = 3.**
* **m = 133: χ = 3.**
  * Independent certificate: a **7-cycle** in Q(√7)², namely (0,0) → (√7/4, 3/4) → (√7/2, 3/2) → (√7/4, 9/4) → (0,3) → (0,2) → (0,1) → (0,0).
  * Upper bound: Lemma R at 3.
  * I could not check the "9-cycle" attribution to Madore. The value itself is right.
* **m = 33: χ = 4,** verified without Fischer. Lower bound: the Moser spindle, SAT-UNSAT for 3 colours. Upper bound: the prime over 2.
* **m = 165: χ = 5.** The upper bound (prime over 11) is verified. The lower bound (Heule's graphs in Z[ω₁, ω₃, ω₄]) is cited and was not re-verified here.
* **m = 385: 3–5,** confirmed.

**Corrections.** All lower bounds below come from THM-4558's generalized spindle (two rigid Eisenstein patches through 0 and v ∈ √−3·Z[ω] with |v|² = N = 3n, n Loeschian, the second patch rotated by ω_N ∈ Q(√−(4N−1))). The n-odd condition there is needed only for the upper bound. For each N I verified |v − ω_N v|² = 1 exactly and SAT-UNSAT for 3 colours (`spindle_check.out`).

| m | Stated | Correct | Evidence |
|---|---|---|---|
| 177 | 3–4 | **4** | N = 723 = 3·241, v = −31 − 14ω, ω = (1445 + 7√−59)/1446. So Q(√−3, √−59) ⊂ F(i). Upper bound: the prime over 2 |
| 345 | 3–5 | **4–5** | N = 144 = 12², ω = (287 + 5√−23)/288. Q(√−3, √−23) ⊂ F(i) |
| 357 | ≥ 3 | **≥ 4** | N = 3600 = 60², 4N − 1 = 7·17·11². Q(√−3, √−119) ⊂ F(i) |
| 253 | 2–4 | **3–4** | 11-cycle in Q(√11)²: six alternating steps (±√11/6, 5/6), then five steps (0, −1) (`oddcycle_more.out`). In general, √p ∈ F with p ≡ 3 (mod 4) gives a p-cycle, so for these 19 planes χ = 2 iff m has no prime factor ≡ 3 (mod 4) |

**A false sentence.** "No conjugation-fixed prime below 60" is false for two planes:
* m = 105: Frob₅₉ = c.
* m = 357: 47 and 59 are both c-fixed.

It is true for m = 1365, where the first c-fixed prime is 131. The reader's own `exp2_analysis.out` lists 59, and 47 and 59, for these planes. The script `oai3_20261007_hcf_planes.py` prints "None" both when no prime is c-stable and when κ is unknown. Its spindle search tries only n ∈ {3, 7, 9, 13}, which is why it misses 144, 723 and 3600.

**Consistency.** No spindle field fits m = 21, 57, 93 or 273 (`spindle_search.out`), as their χ = 3 requires.

---

## 2. THM-4567: CONFIRMED

All checks are in `jacobian2_audit.out` and `jacobian_collisions.out`.

1. **Mod 2.**
   * det JF = −2.
   * F₂y² + F₂³ + F₁²F₃ has all coefficients even.
   * Two further identities hold mod 2: x²(F₂ + y²F₃) ≡ F₃ and x⁶z² ≡ F₃² + x⁴y².
   * Some 2×2 minor of JF is nonzero mod 2.
   * Hence K² ⊂ F₂(F) and [K : F₂(F)] = 2^(3 − rank) = 2. So F mod 2 is dominant and purely inseparable of degree 2.
   * The theorem's "So" uses only the first identity. The reader's c5b supplies the rest (x ∈ F₂(F, y), z ∈ F₂(F, x, y), and degree ≠ 1 because det ≡ 0).
2. **2-adic fibres.**
   * Hensel is rigorous from level 2: since J(c)Z₂³ ⊇ 2Z₂³ (index |det| = 2), F(c + 2^kZ₂³) = F(c) + 2^kJ(c)Z₂³ bijectively for k ≥ 2. With k = 2 this shows N(w) depends only on w mod 8.
   * Exact computation over the 64 classes mod 4 gives N = 0, 1, 2, 3 on 336, 112, 48, 16 of the 512 classes. So μ = 11/32 and ∫N = 1/2.
   * Merge-partner law: P(j other preimages) = 2(j+1)·#/512 = 7/16, 3/8, 3/16.
   * Direct counts mod 2^k give 11/32 for k = 3..7, with fibre histogram {2:112, 4:48, 6:16} scaling by 8 per level. "Each preimage fills 2 classes mod 8" is correct.
3. **Conservation.**
   * Classes (0,1,0), (0,1,1), (1,0,1): exactly one in-disc partner and no other partner, on all 64 subclasses.
   * The other five classes: no in-disc partner.
   * Explicit curves in the mod-2 fibres pass through the five improper points (the x-axis; {y = xz, x²z = 1}; {y = xz, x²z² + z = 1}; {xy = 1, z = (1+x)/x³}).
4. **Collisions.**
   * The family identity holds symbolically for w = 2t + 1, and all coordinates are integral.
   * F(1,−1,5) = F(−1,2,8) = F(0,2,−16) = (0,2,0).
   * It is the plane family at s² = v² − 16u = 4.
   * The smallest collision by max-norm is F(2,−1,2) = F(0,−1,−4), together with its mirror (−x, −y, z).
5. **Local multiplicities.** These equal 1 trivially, since the map is étale over Q.

---

## 3. Updates

### THM-4555 update: CONFIRMED (`collatz_audit.out`)

* N_w is right: f_w(−1) = N_w/2^(A−1), checked on six words.
* The gap lemma proof is sound. Every partial sum preceding exponent B_(k+1) is a nonzero multiple of 2^(B_(k+1)), so 2^(B_(k+1)) ≤ Σ|c_j|·2^(B_k).
* f_(2,2,10,a)(−1) = f_(6,3,2,1,a+2)(−1) = 8207·2^(−13−a), symbolically in a.
* The switch: 53803 has exponents 1, 2, 2, 10, 4 and 26901 has exponents 6, 3, 2, 1, 6. Both reach U⁵ = 25, and step 5 is their first meeting.
* Census by my own code: 12 sporadic values for letters ≤ 8 and 27 for letters ≤ 22. Of the 15 new ones, 14 are in the (2,6,c) ~ (4,1,1,c+2) family and the remaining one is 8207.

### THM-4512 alternative closure: CONFIRMED WITH CORRECTION (`ellison_check.out`)

* **Citation.** It is verbatim in Waldschmidt's survey (`pillai_survey.txt`, lines 326–337): "for x ≥ 12 with x ≠ 13, 14, 16, 19, 27, and all y". The author is W. J. Ellison (Sém. Théorie des Nombres Bordeaux, Exp. 12).
* **Exceptions.** Checked with exact integers and a rigorous interval enclosure of e^(x/10): exactly {13, 14, 16, 19, 27} on [12, 20000]. A float scan finds none on (20000, 10⁶].
* **Wrong constant.** The bound is 2^(−j)e^(A/10) ≤ e^0.1·e^(−(ln 2 − log₂3/10)j), and ln 2 − log₂3/10 = 0.53465, which is less than 0.535.
  * So "< 1.11·e^(−0.535j)" does not follow. The ratio of the derived bound to it is 1.016 at j = 65 and 30.7 at j = 10⁴.
  * The conclusion survives: max over j ≥ 65 is 8.93·10⁻¹⁶ < 10⁻¹⁵.

### THM-4558 pointer: CONFIRMED

It should also mention the third 4-chromatic Hilbert-class-field plane, m = 177.

---

## 4. Note §4–5

| Claim | Verdict |
|---|---|
| #003 excludes Siegel zeros, then Tatuzawa and EKN | CONFIRMED (see 1b) |
| "#003's own 11/12 paper says so" | CONFIRMED |
| 101 discriminants of exponent ≤ 2, including −4n for the 65 idoneal n | CONFIRMED. Caution: separately, 65 of the 101 are fundamental (EKN: 9 + 56). This is a different 65-element subset; for example −3 is fundamental but not −4n, and −16 is −4·4 but not fundamental |
| Borwein–Choi "every positive integer except 18 is xy + yz + zx" | The mathematics is CONFIRMED: exceptions up to 3·10⁶ are exactly {1, 2, 4, 6, 10, 18, 22, 30, 42, 58, 70, 78, 102, 130, 190, 210, 330, 462}, all idoneal, and non-idoneal n are ab + bc + ca with 0 < a < b < c. The **wording** reads as "except the number 18" and must be fixed |
| EKN exponent-4 (203) and exponent-8 (778) lists complete | CONFIRMED. The counts are EKN Table 1. Thm 2's no-Siegel-zero hypothesis is the Lemma 3 / Chen zero, which lies above 7/8 |
| Class-group exponents → ∞ effectively | CONFIRMED with justification. The Re s > 7/8 strip gives a least split prime ≪ (log\|D\|)^8. If 𝔭^e is principal then p^e ≥ \|D\|/4, so E(D) ≫ log\|D\|/log log\|D\|. EKN note this is open unconditionally |
| CM labels 30, 42, 70, 105, 210 are all idoneal | CONFIRMED. The labels are in #004's own text (`h10/04-height.tex`, lines 113–172); h = 4, 4, 4, 8, 8 with class groups (Z/2)^r; there are 8 fixed points each, 40 in total. "Forced" needs both the Riemann–Hurwitz count (one branch value per label, so (Z/2)³ acts simply transitively on Fix(w_m)) and Shimura reciprocity (Pic acts freely and commutes with W). That is correctly flagged "modulo Shimura reciprocity" |
| Serre positivity "has content only for non-CM modules in dimension ≥ 4 over a ramified base" | True as a necessary condition, but **not sharp**. See below |
| Fano gadget = closed in-neighbourhood design of P₇ | CONFIRMED |
| "Every DRT of order ≡ 7 (mod 8) gives the same gadget" | Imprecise. What holds in general, proved here, is the parity structure, not the same gadget. See below |
| THM-4567 and THM-4555 summaries | CONFIRMED |

**Why the Serre threshold is not sharp.**
* If both modules are CM, Tor vanishes over any regular local ring: Tor_i(M,N) = Ext^(c−i)(M^∨, N), and grade(ann M^∨, N) = depth N = c.
* If dim R ≤ 4, a complementary pair of prime quotients has dimensions (≤ 1, ≥ d − 1), where both are CM (a 1-dimensional domain and a principal hypersurface), or (2, 2). In the (2, 2) case each complete local domain D has a CM normalization D̄, and Roberts/Gillet–Soulé vanishing gives χ(D, E) = χ(D̄, Ē) = ℓ(D̄ ⊗ Ē) > 0.
* Serre's argument also covers V[[x]] for a ramified complete DVR V (#193, introduction).
* So "dimension at least 4" should read "at least 5".

**The DRT gadget in general.** For a DRT of order n = 4t + 3, AAᵀ = AᵀA = (t+1)I + tJ, so the blocks N⁻[i] form a 2-(n, (n+1)/2, (n+1)/4) design. If n ≡ 7 (mod 8), blocks and intersections are even, so the F₂-span is self-orthogonal; D + Dᵀ = J + I has F₂-rank n − 1, so its dimension is exactly (n−1)/2.

---

## Required corrections (exact replacement wording)

1. **THM-4566 title:** "4 for m = 33 (Fischer)" → "4 for m = 33 (Fischer) and m = 177".
2. **THM-4566 table, row 177:** "| 177 | Q(√3, √59) | **4** | generalized spindle N = 723 = 3·241 (THM-4558), since Q(√−3, √−59) ⊂ F(i); prime over 2 |".
3. **Row 253:** "| 253 | Q(√11, √23) | 3–4 | 11-cycle in Q(√11)²; prime over 7 |".
4. **Row 345:** "| 345 | Q(√3, √5, √23) | 4–5 | generalized spindle N = 144, Q(√−3, √−23) ⊂ F(i); prime over 11 |".
5. **Row 357:** "| 357 | Q(√3, √7, √17) | ≥ 4 | generalized spindle N = 3600, Q(√−3, √−119) ⊂ F(i); its conjugation-fixed primes below 60 (47, 59) have unknown κ |".
6. **Row 105/1365,** replace "no conjugation-fixed prime below 60" with: "no usable conjugation-fixed prime (for 105 the only one below 60 is 59, where κ(59) ≥ 6; for 1365 none below 131)".
7. **THM-4512 update:** "N(w)/2^A < 1.11·e^(−0.535 j) < 10^(−15)" → "N(w)/2^A < 2^(−j)e^(A/10) ≤ e^(0.1)·e^(−0.5346 j) ≤ 8.93·10^(−16) < 10^(−15)".
8. **Note §4, Borwein–Choi:** → "every positive integer except the 18 numbers 1, 2, 4, 6, 10, 18, 22, 30, 42, 58, 70, 78, 102, 130, 190, 210, 330, 462 is xy + yz + zx with x, y, z ≥ 1".
9. **Note §4, chromatic bullet:** add "4 for m = 33 and 177; 5 for m = 165", and use the corrected ranges from items 2–6.
10. **Note §5, Serre:** → "Positivity was already known in equicharacteristic and over power series rings over any complete DVR (Serre), when both modules are CM, and for dim R ≤ 4 (small CM modules plus Roberts/Gillet–Soulé), as well as in the Skalit and KC–Soto Levins cases. #193's new content lies in ramified mixed characteristic with dim R ≥ 5 and a prime quotient of dimension ≥ 3."
11. **Note §5, Fano gadget:** "gives the same gadget" → "gives an even gadget of the same kind: a 2-(n, (n+1)/2, (n+1)/4) design with even blocks and intersections spanning a self-orthogonal F₂-code of dimension (n−1)/2".

## Optional improvements

* **THM-4566 Identification:** add the "F Galois ⇒ c central" step (1a).
* **THM-4566 completeness:** state the EKN threshold as d₁₁ = 200560490130.
* **THM-4566 row 133:** cite the explicit 7-cycle instead of, or alongside, Madore's 9-cycle.
* **THM-4566, statement of χ = 2:** "for these planes χ = 2 iff m has no prime factor ≡ 3 (mod 4)".
* **`oai3_20261007_hcf_planes.py`:** separate "no c-stable prime" from "κ unknown"; search all Loeschian n.
* **THM-4567 statement 1:** add the two extra mod-2 identities and the rank argument.
* **THM-4567, THM-1345 reference:** "THM-1345" is a known double ID. Cite `THM-1345-plane-family-section-radical-inverse-trace-module.md`.
* **THM-4567 status:** the partner law can be upgraded to PROVED (level-2 Hensel plus the exact 64-class enumeration).

## Scripts in this directory

`exp2_audit.c` (+ `exp2_sieve_2p1e11.txt/.err`, `exp2_direct_1e8.txt`), `idoneal_check.py`, `lemmaR_audit.py`, `spindle_check.py`, `spindle_search.py`, `oddcycle_lib.gp`, `oddcycle_sqrt7.gp`, `oddcycle_more.gp`, `kappa_check.py`, `jacobian2_audit.py`, `jacobian_collisions.py`, `collatz_audit.py`, `ellison_check.py`; each with its `.out`.
