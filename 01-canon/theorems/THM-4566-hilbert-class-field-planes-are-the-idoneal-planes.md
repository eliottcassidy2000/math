---
id: THM-4566
title: "Hilbert-class-field planes are the idoneal planes: for an imaginary quadratic K, H(K) is the plane F(i) of a real field F Galois over Q iff K = Q(sqrt-m) with m squarefree, m = 1 mod 4 and idoneal, and then F = Q(sqrt p : p | m); there are exactly 19 such m (complete modulo openai/math #003, the quasi-Riemann hypothesis); their chromatic numbers: 2 for m = 1, 5, 13, 37, 85; 3 for m = 21, 57, 93, 133, 273; 4 for m = 33 (Fischer); 5 for m = 165 (the Heule/Polymath plane); bounds for the rest"
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
| 133 | `Q(√7, √19)` | **3** | `χ(Q(√7)²) = 3` (Madore 2015, a 9-cycle); the prime over 3 is inert |
| 33 | `Q(√3, √11)` | **4** | KNOWN (Fischer 1994; THM-4558) |
| 165 | `Q(√3, √5, √11)` | **5** | contains the Polymath field `Q(√−3, √−11, √−15)` (Heule's graphs); the prime over 11 is inert, `κ(11) = 5` |
| 177 | `Q(√3, √59)` | 3–4 | triangle; prime over 2 |
| 253 | `Q(√11, √23)` | 2–4 | prime over 7 |
| 345 | `Q(√3, √5, √23)` | 3–5 | prime over 11 |
| 385 | `Q(√5, √7, √11)` | 3–5 | `Q(√7)²`; prime over 19 (`κ(19) = 5` per decalion89) |
| 105, 1365 | contain `√3, √5, √7` | `≥ 4` | the spindle field `Q(√−3, √−35)` (THM-4558, `N = 9`); no conjugation-fixed prime below 60 |
| 357 | `Q(√3, √7, √17)` | `≥ 3` | no conjugation-fixed prime below 60 |

* The two "record" planes of the Hadwiger–Nelson field tower are both Hilbert class fields of idoneal discriminants: the Moser plane (`m = 33`, χ = 4) and the plane of the Polymath field (`m = 165`, χ = 5).
* Any further pattern linking χ to the class number (1, 2, 4, 8, 16) is NUMEROLOGY.

## Further corollaries of the same completeness (modulo #003; recorded, not re-derived here)

* **Borwein–Choi.** Every positive integer outside `{1, 2, 4, 6, 10, 18, 22, 30, 42, 58, 70, 78, 102, 130, 190, 210, 330, 462}` is `xy + yz + zx` with `x, y, z ≥ 1`. Borwein–Choi (2000) had left at most one possible further exception.
* **Euler/Cox.** The `n` for which `p = x² + ny²` is decided by congruences mod `4n` are exactly the 65 idoneal numbers.
* **EKN.** The lists of imaginary quadratic fields with class-group exponent 4 (203) and 8 (778) are complete.
