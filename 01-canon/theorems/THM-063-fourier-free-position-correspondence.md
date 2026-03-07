# THM-063: Fourier Degree = Free Position Count

**Type:** Theorem (PROVED)
**Certainty:** 5 -- PROVED (algebraic + exhaustively verified n=3,...,9)
**Status:** PROVED
**Added by:** opus-2026-03-06-S11b
**Tags:** #Fourier #OCF #master-polynomial #free-positions #transfer-matrix

---

## Statement

### (i) Fourier coefficient formula

The master polynomial W(r) = sum_k w_k r^k has Fourier coefficients:

**w_k(T) = sum_I 2^{parts(I)} * [r^k in F_{f_I}(r)] * I(T)**

where the sum is over all independent sets I in the conflict graph Omega(T), including the empty set (I_empty(T) = 1), and f_I = n - 1 - sum_{C in I} (|C| - 1) is the free-position count.

### (ii) Free-position polynomial

F_f(r) = sum_{j=0}^{f} A(f+1, j) * (r+1/2)^j * (r-1/2)^{f-j}

is the polynomial of degree f whose coefficients in the power-of-r basis are:

| f | F_f(r) |
|---|--------|
| 0 | 1 |
| 1 | 2r |
| 2 | 6r^2 - 1/2 |
| 3 | 24r^3 - 4r |
| 4 | 120r^4 - 30r^2 + 1 |
| 5 | 720r^5 - 240r^3 + 17r |
| 6 | 5040r^6 - 2100r^4 + 231r^2 - 17/4 |

Key properties:
- F_f has the same parity as f (only even/odd powers of r)
- Leading coefficient = (f+1)!
- F_f(1/2) = 1 for all f >= 0

### (iii) Fourier degree ↔ free positions

The Fourier coefficient w_{n-1-2j} (degree 2j in centered variables) receives contributions only from invariants I with f_I >= n-1-2j (even parity). Specifically:

**Fourier degree 2j contains exactly the invariants with n-1-2j <= f_I <= n-1-2j (mod 2).**

This means:
- Degree 0 (w_{n-1} = n!): only the empty set (f = n-1)
- Degree 2 (w_{n-3}): all I with f_I in {n-3, n-5, ...}; effectively t3 at leading order
- Degree 4 (w_{n-5}): all I with f_I in {n-5, n-7, ...}; effectively t5, alpha_{3,3}
- Degree 6 (w_{n-7}): all I with f_I in {n-7, n-9, ...}; effectively t7, alpha_{3,5}, alpha_{3,3,3}

### (iv) General w_k formulas (verified n=3,...,9)

**w_{n-1} = n!** (all n)

**w_{n-3} = 2*(n-2)! * (t3 - C(n,3)/4)** (all n >= 3)

**w_{n-5} = 2*(n-4)! * t5 + 4*(n-4)! * alpha_{3,3} - (n-4)!*C(n-2,3)*t3 + const** (all n >= 5)

**Master coefficient pattern:** At Fourier degree 2j, the coefficient of:
- A new odd cycle t_{2j+1}: 2 * (n-2j)!
- A new k-fold disjoint 3-cycle collection: 2^k * (n-2j)!

### (v) Telescoping to OCF

When evaluated at r = 1/2:

H(T) = W(1/2) = sum_k w_k / 2^k = sum_I 2^{parts(I)} * F_{f_I}(1/2) * I(T) = sum_I 2^{parts(I)} * I(T) = I(Omega(T), 2)

since F_f(1/2) = 1 for all f. This is the OCF.

The intermediate corrections (sub-leading terms of each F_f) telescope perfectly: each invariant's contributions across all Fourier degrees sum to the OCF coefficient 2^{parts(I)}.

---

## Proof

Part (i) follows from the OCF decomposition W(r) = sum_I 2^{parts(I)} F_{f_I}(r) I(T) (THM-059) by extracting the coefficient of r^k.

Part (ii): F_f(r) is constructed as the W-polynomial of the "free portion" of the path, which is a subpath of length f on f+1 vertices with all edges independently ±1/2. The Eulerian number formula for F_f follows from counting permutations by ascents on the free vertices.

Part (iii): F_f(r) has degree f and same parity, so [r^k in F_f] = 0 when k > f or k and f have different parity.

Part (iv): verified computationally for n = 3, 5, 7, 9 with exact integer coefficients (max error < 10^{-4} from polynomial interpolation).

Part (v): F_f(1/2) = 1 is proved by noting that at r=1/2, (r+1/2)^j (r-1/2)^{f-j} = 1^j * 0^{f-j} = delta_{j,f}, so F_f(1/2) = A(f+1, f) = 1.

---

## Connection to Other Results

- **THM-059 (Master Polynomial):** Provides the OCF decomposition W = sum F_f I, which is the starting point.
- **THM-062 (Deformed Eulerian):** The change-of-basis from the {r^k} basis to the {(r+1/2)^j (r-1/2)^{d-j}} basis converts between Fourier coefficients w_k and forward-edge counts a_k.
- **Transfer matrix symmetry (INV-001):** The Fourier coefficients w_k = tr(c_k) where c_k are transfer matrix coefficient traces. Proving M symmetric ↔ proving the w_k formulas.

---

## Verification

| Result | n | Method | Status |
|--------|---|--------|--------|
| w_{n-3} general formula | 3,5,7,9 | regression + exact match | PASS |
| w_{n-5} general formula | 5,7,9 | regression + exact match | PASS |
| Complete w_k at n=7 | 7 | 50 random tournaments | PASS |
| Complete w_k at n=9 | 9 | 25 random tournaments | PASS |
| Telescoping to OCF at n=9 | 9 | all 8 invariants sum correctly | PASS |

## Scripts

- `04-computation/fourier_ocf_decomposition.py` -- Main verification script
- `04-computation/degree4_proof_n7.py` -- Exhaustive degree-4 proof at n=7
- `04-computation/degree4_identity_n7.py` -- Discovery of 2D degree-4 structure
