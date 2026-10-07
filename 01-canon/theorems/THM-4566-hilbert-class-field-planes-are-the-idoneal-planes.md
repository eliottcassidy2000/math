---
id: THM-4566
title: "Hilbert-class-field planes are the idoneal planes: for an imaginary quadratic K, H(K) is the plane F(i) of a real field F Galois over Q iff K = Q(sqrt-m) with m squarefree, m = 1 mod 4 and idoneal, and then F = Q(sqrt p : p | m); there are exactly 19 such m (complete modulo openai/math #003, the quasi-Riemann hypothesis); their chromatic numbers: 2 for m = 1, 5, 13, 37, 85; 3 for m = 21, 57, 93, 133, 273; 4 for m = 33 (Fischer) and m = 177; 5 for m = 165 (the Heule/Polymath plane); bounds for the rest; chi = 2 iff m has no prime factor = 3 mod 4"
status: "PROVED: the identification (genus theory). FINITE-EXACT: the list of 65 idoneal numbers and the 19 values (sieve of all negative discriminants to 2.34e14 by the session's nt reader, equal to OEIS A003171 / A000926). Completeness PROVED modulo openai/math #003 (accepted per owner directive 2026-10-07): the zero-free half-plane Re s > 7/8 excludes Siegel zeros, so Tatuzawa / Elsenhans-Kluners-Nicolae (2020) apply without exception. #003's own 11/12 paper states the idoneal consequence. Chromatic values PROVED via Lemma R (THM-4558; a known technique) with the lower bounds cited; my own decomposition-group computation reproduces every upper bound."
session: mac-mini-2026-10-07-oaimath3 (owner listed the 65 idoneal numbers in the prompt; nt reader + this session)
source: 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
scripts:
  - 04-computation/experiments/oai3_20261007_readers/nt/ (exp2_sieve.c + outputs, exp2_analysis.py, l1_lower_bound.py)
  - 04-computation/experiments/oai3_20261007_hcf_planes.py (+ .out; Lemma R bounds for the 19 planes)
related:
  - THM-4558 (residue colourings, field-plane values; the Moser plane is the genus field of Q(sqrt33)), THM-418
  - openai/math #003 (quasi-RH), #004 (H10 over Q: its height bound's contact labels 30, 42, 70, 105, 210 are idoneal)
  - Weinberger (1973), Kani (Ann. Sci. Math. Quebec 2011, idoneal numbers survey), Elsenhans-Kluners-Nicolae arXiv:1803.02056, Borwein-Choi (Exp. Math. 2000), D. Madore arXiv:1509.07023 (chi(Q(sqrt7)^2) = 3), K. G. Fischer (1994)
---

# THM-4566 — Hilbert-class-field planes are the idoneal planes

## Identification (PROVED)

Let `K = Q(√−m)` be imaginary quadratic.
* `H(K)` is abelian over Q iff `H(K)` is the genus field, iff `Cl(K)` has exponent at most 2.
* If `H(K) = F(i)` with `F` real and Galois over Q, then `F = H ∩ R` and `⟨c⟩ = Gal(H/F)` is normal of order 2. So complex conjugation `c` is central in `Gal(H/Q) = Cl(K) ⋊ ⟨c⟩`, where it acts by inversion, and the exponent is at most 2. (Step added after audit B.)
* `i ∈ H(K)` iff the prime discriminant `−4` divides `d_K`, iff `m ≡ 1 (mod 4)` (with `m` squarefree, or `K = Q(i)`).
* In that case the genus field is `Q(i, √p : p | m)`.

**Hence `H(K) = F(i)` with `F` real and Galois over Q iff `m` is a squarefree idoneal number `≡ 1 (mod 4)`, and then `F = Q(√p : p | m)`.** The 19 values are

    1, 5, 13, 21, 33, 37, 57, 85, 93, 105, 133, 165, 177, 253, 273, 345, 357, 385, 1365.

* The list is FINITE-EXACT for `m ≤ 5.8·10^13`.
* Its completeness is PROVED modulo openai/math #003, through:
  * Tatuzawa's lower bound for `h(D)` with no exceptional discriminant, since a Siegel zero would have real part above 7/8;
  * EKN's Theorem 4: no fundamental `|D| ≥ 2·10^11` has exponent 2;
  * the session's sieve over all negative discriminants up to `2.34·10^14` (101 discriminants of exponent at most 2; 65 idoneal numbers);
  * a conductor lemma: exponent at most 2 forces the conductor to divide 840.

## Chromatic numbers

All upper bounds come from Lemma R (THM-4558), recomputed independently in `oai3_20261007_hcf_planes.py`.

| `m` | `F` | `χ(F²)` | why |
|---|---|---|---|
| 1, 5, 13, 37, 85 | no prime `≡ 3 (mod 4)` | **2** | 2 ramifies in `F(i)/F` |
| 21, 57, 93, 273 | contains `√3` | **3** | equilateral triangles; the prime over 3 is inert, `κ(3) = 3` |
| 133 | `Q(√7, √19)` | **3** | `χ(Q(√7)²) ≥ 3` by the explicit 7-cycle `(0,0) → (√7/4, 3/4) → (√7/2, 3/2) → (√7/4, 9/4) → (0,3) → (0,2) → (0,1) → (0,0)` (audit B; cf. Madore 2015); the prime over 3 is inert |
| 33 | `Q(√3, √11)` | **4** | KNOWN (Fischer 1994; THM-4558) |
| 165 | `Q(√3, √5, √11)` | **5** | contains the Polymath field `Q(√−3, √−11, √−15)` (Heule's graphs); the prime over 11 is inert, `κ(11) = 5` |
| 177 | `Q(√3, √59)` | **4** | generalized spindle `N = 723 = 3·241` (THM-4558), since `Q(√−3, √−59) ⊂ F(i)`; SAT-UNSAT for 3 colours (audit B); prime over 2 |
| 253 | `Q(√11, √23)` | 3–4 | 11-cycle in `Q(√11)²` (audit B); prime over 7 |
| 345 | `Q(√3, √5, √23)` | 4–5 | generalized spindle `N = 144`, since `Q(√−3, √−23) ⊂ F(i)` (audit B); prime over 11 |
| 385 | `Q(√5, √7, √11)` | 3–5 | `Q(√7)²`; prime over 19 (`κ(19) = 5` per decalion89) |
| 105, 1365 | contain `√3, √5, √7` | `≥ 4` | the spindle field `Q(√−3, √−35)` (THM-4558, `N = 9`); no usable conjugation-fixed prime (for 105 the only one below 60 is 59, where `κ(59)` is unknown; for 1365 none below 131) |
| 357 | `Q(√3, √7, √17)` | `≥ 4` | generalized spindle `N = 3600`, since `Q(√−3, √−119) ⊂ F(i)` (audit B); its conjugation-fixed primes below 60 (47, 59) have unknown `κ` |

* The two "record" planes of the Hadwiger–Nelson field tower are both Hilbert class fields of idoneal discriminants: the Moser plane (`m = 33`, χ = 4) and the plane of the Polymath field (`m = 165`, χ = 5).
* Any further pattern linking χ to the class number (1, 2, 4, 8, 16) is NUMEROLOGY.

## Further corollaries of the same completeness (modulo #003; recorded, not re-derived here)

* **Borwein–Choi.** Every positive integer outside `{1, 2, 4, 6, 10, 18, 22, 30, 42, 58, 70, 78, 102, 130, 190, 210, 330, 462}` is `xy + yz + zx` with `x, y, z ≥ 1`. Borwein–Choi (2000) had left at most one possible further exception.
* **Euler/Cox.** The `n` for which `p = x² + ny²` is decided by congruences mod `4n` are exactly the 65 idoneal numbers.
* **EKN.** The lists of imaginary quadratic fields with class-group exponent 4 (203) and 8 (778) are complete.

**Audit (2026-10-07, independent audit B).**
* CONFIRMED:
  * the identification, with the "F Galois ⇒ c central" step added;
  * the 19 values, by an independent enumeration of all discriminants to `2.1·10^11` (exactly 101, equal to A003171);
  * completeness modulo #003, checked against EKN and #003's 11/12 paper; EKN's threshold is `d_11 = 200560490130`;
  * every upper bound and its prime;
  * `κ(3) = 3`, `κ(7) = 4`, `κ(11) = 5` and `κ(19) ≤ 5` by SAT;
  * the `χ = 3` rows, and `m = 33` without Fischer (spindle UNSAT).
* Improved here from the repo's own spindle construction (the `n`-odd condition of THM-4558 is needed only for upper bounds): `m = 177` is exactly 4, `m = 345` is 4–5, `m = 357` is `≥ 4`, and `m = 253` is 3–4.
* The sentence "no conjugation-fixed prime below 60" was false for 105 (59) and 357 (47, 59) and is corrected (MISTAKE-583).
* For these 19 planes, `χ = 2` iff `m` has no prime factor `≡ 3 (mod 4)`. A prime `p ≡ 3 (mod 4)` with `√p ∈ F` gives a `p`-cycle.
* The `m = 165` lower bound (Heule's graphs) is cited, not re-verified.

---

## Addendum (2026-10-07, mac-mini-2026-10-07-golden, platonic reader): the 18 and the 19 are two 2-adic slices of Euler's list

**Theorem (KNOWN: Borwein–Choi, Exp. Math. 9 (2000), Thm 3.1 [squarefree `n ≡ 2 mod 4`: exception iff `−4n` has one class per genus], Thm 2.6 [the only non-squarefree exceptions are 4 and 18], Lemma 2.2; re-proved below by an independent Selling-parameter argument. The explicit list relies on the completeness of the idoneal numbers, i.e. on openai/math #003 as above.)**
`n` is **not** of the form `xy + yz + zx` with `x, y, z ≥ 1` iff `n ∈ {1, 4}`, or `n ≡ 2 (mod 4)` and `n` is idoneal.

*Proof sketch.*
* Odd `n ≥ 3` is `1·y + y·1 + 1·1`. If `4 | n` and `n ≥ 8`, take `(2, (n−4)/4, 2)`.
* For `n ≡ 2 (mod 4)`: `n = xy + yz + zx` with `x, y, z ≥ 0` iff the form `(x+z)X² − 2zXY + (y+z)Y²` (determinant `n`) has an obtuse superbase with Selling parameters `(x, y, z)`. These parameters are a class invariant (Conway).
* So `n` is an exception iff every positive definite integral form of determinant `n` is `GL_2(Z)`-equivalent to a diagonal form (a zero Selling parameter).
* Such a form is `g·f′`, with odd content `g` and `f′` primitive of discriminant `−4n′`, where `n′ = n/g² ≡ 2 (mod 4)` (so `f′` is Gauss-primitive).
* For `n′ ≡ 2 (mod 4)` the reduced primitive diagonal forms number `2^r`, with `r` the number of odd primes dividing `n′`. That is the number of genera (`μ = r + 1`), hence by Gauss the number of ambiguous classes. Diagonal classes are ambiguous, so they are exactly the ambiguous classes.
* Hence all forms are diagonal iff `Cl(−4n′)` has exponent `≤ 2`, i.e. iff `n′` is idoneal. Exponent `≤ 2` descends from `n` to `n/g²`, because `Cl(−4n) ↠ Cl(−4n/g²)`.
* The count equality also holds for `n′ ≡ 3 (mod 4)` and `n′ ≡ 4 (mod 8)`. The earlier text's "exactly when" was false; corrected after audit B.

**FINITE-EXACT checks.**
* A brute force over `n ≤ 2·10^5` gives exactly the 18 exceptions.
* For `|D| ≤ 30000` there are exactly 101 discriminants of exponent `≤ 2` (65 even, 36 odd).

**Reading.**
* The 18 Borwein–Choi exceptions are the 16 idoneal `n ≡ 2 (mod 4)` together with 1 and 4. The 19 planes are the squarefree idoneal `m ≡ 1 (mod 4)`. Both are 2-adic slices of one 65-element list, because `m ≡ 1 (mod 4)` ⟺ `−4` is the 2-adic prime discriminant ⟺ `i` lies in the genus field.
* So the counts 18 and 19 carry no mechanism linking them to the Collatz mod-18/19 clock tower (`ord_19(2) = 18`): NUMEROLOGY. (1 lies in both slices.)
* The one generic DICTIONARY item: on the planes 57 and 133, the 19-genus character `(a|19)` equals `(−1)^(ind_2 a)`, the parity bit of the mod-19 clock. This holds at any prime where 2 is a primitive root.

**Split-prime lemma (PROVED, elementary).** In an order of exponent `≤ 2`, a split prime `p` satisfies `p² ≥ |D|/4`. Proof: `𝔭` is invertible and `𝔭² = (α)`. Since `𝔭² ≠ (p)` for split `p`, `α ∉ Z`. So `α = (a + b√D)/2` with `b ≠ 0`, and `p² = N(α) = (a² + |D|b²)/4 ≥ |D|/4`. So for idoneal `n > 9` the prime 3 does not split (`n ≢ 2 mod 3`), and 2 splits only for `|D| ∈ {7, 15}`. In large exponent-2 fields, 2 and 3 are inert or ramified.

Scripts: `04-computation/experiments/golden_20261007_readers/platonic/idoneal.py` (+ `.out`).

**Audit of the addendum (2026-10-07, independent audit B).**
* The characterization was brute-forced to `n ≤ 10^7` with an independent idoneal test: 0 mismatches.
* There are 101 discriminants of exponent `≤ 2`, and the split-prime lemma is checked.
* Prior art: Borwein–Choi 2000. The theorem is retyped KNOWN, and the proof sketch is repaired (MISTAKE-584).
