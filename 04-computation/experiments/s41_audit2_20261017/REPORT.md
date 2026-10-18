# S41 second audit: THM-4668 and the product-direction censuses

Independent referee (general-purpose subagent), 2026-10-17. Audited HEAD cd99bb7a46; THM-4668 is unchanged since
7af0b145fe, so the line numbers hold. The referee could not write this file itself; it was saved from the referee's
final message by the session. A session checkpoint (cd99bb7a46) swept the referee's files in mid-run:
b5_forcing4_check.out went in empty (the full output is in the working tree), and b6_identities.py/.out and
b7_session_proof_check_rerun.out were untracked. All are committed together with this report.

**Verdict.** Two errors and one overstatement, all in wording and citations. The proof is correct and complete, and
PROVED stands.

## Errors

1. **Census size (THM line 8).** The e-census covers `dA ∈ D(A) ∪ {0}`, `d1, d2 ∈ D(C) ∪ {0}`, and
   `d3 ∈ (D(C) ∪ D(L)) \ {0}`. That is `1603 x 83 x 83 x 130 = 1,435,598,710 ≈ 1.44 x 10^9`, recounted exactly. The
   figure 9.06e8 belongs to the v^2 census. Fix:
   - write the count as 1.44e9;
   - describe the method as a three-prime sieve followed by exact checks of the 26 survivors;
   - list all four coefficient slots with their zeros.
2. **Wrong lemma cited (THM lines 17, 38, 84).** These cite "G5" of joint_contact_field for the units of
   `K + uK + vK`.
   - G5 (lines 95-106 there) is the three-contact theorem for the spindle contacts. It does not cover `d = 1 + zeta5`.
   - The classification actually used (axes plus `q(d, 1, 1)`) is lemma G1 of marked_contacts_root_current_20261017e.md
     section 5 (lines 175-203), with `e = 1`.
   - Its hypotheses for the golden contact are certified in joint_contact_field section 5 (lines 187-199): `|d|^2 ≠ 1`,
     `0 < N < 4`, and `s` not in `K(u)`.
   - The same mis-citation sits in the docstrings of golden_forcing_search4_20261017.py (line 8) and
     golden_forcing_search_20261017.py (line 4).
   - The content is right: the referee's census of the whole source difference box finds exactly G1's list.

## Overstatements

3. **Colouring extension without its hypothesis.** Three places drop the corollary's condition `k >= χ(L)`: the THM title
   (line 3), the note (lines 126-127), and HYP-9289 (line 36).
   - The corollary itself requires it. For six colours, L's unit graph must be 6-colourable, which holds for
     clock_duality's 13-cycle but is not known for arbitrary finite L.
   - The claim also assumes the ports lie in one e-layer.
   - Fix: "every k-colouring with k >= χ(L) extends, ports in one e-layer".
4. **Scope (minor).** HYP-9289 line 57 ("ue and ve add no cross edges") and note lines 127-128 omit the scope. The
   result is FINITE-EXACT only for the coefficient sets D(A), D(C), D(C) ∪ D(L).

## Minor

5. **λ derivation (THM line 76).** After eliminating κ' with the u^2-equation, the u-equation is quadratic in λ:
   `d̄ λ^2 + (2 - 2|d|^2) λ + (d̄ m1 - d - |d|^2 m2) = 0`.
   - The printed formula follows only after substituting `λ^2 = -d̄`. The logic holds, since all three equations are
     necessary.
   - Equivalently, the resultant `Res_λ = q1^2 (λ^2 + d̄)` is nonzero, which the referee computed exactly.
   - The proof-check docstring ("[u] gives lam linearly") has the same wording issue; its code substitutes correctly.
6. **Credit.** The argument "N not in K" and the c ≠ 0 pattern are clock_duality's (G1's proof,
   marked_contacts...e.md lines 201-203). The colouring c + c_L mod k is their formula (3), credited only in the Reading.
7. **The two families (note lines 143-144).** Both are said to "follow from |d + u + v| = 1". That holds for uv. The e^2
   identity uses only `e^2 = -Z/Z̄` and `ū = 1/u`. golden_product_relations_20261017.py identifies the families only in
   floating point; the referee confirmed them exactly in K.
8. **Witnesses.** golden_forcing_search4 records only SAT/UNSAT. The referee's replay (section C) closes the gap.
9. **Wording.** "With distinct addresses" is automatic (1, u, v, e are K-independent), and the "independent" sieve is
   independent of the proof only.

## Verified as correct

**A. The proof, line by line.**
- The four norm components and the coordinates check out symbolically.
- `h` is fixed by conjugation: at the intended embedding `N = 2.7091185` and `R = 1.2909 > 0`. Also `D = -0.4765 < 0`,
  so `s -> -s`.
- Degree 192: D, R and DR are non-squares in F. The referee reproduced clock_duality's two places at p = 31, and over
  320 degree-one places D, R and DR are non-residues at 164, 140 and 152 of them.
- Field degrees: `[K:Q] = 16` and `[F:K] = 3` (the cubic has no root mod 181, 751 or 1291, which split completely in K),
  so `[F:Q] = 48`.
- Case c ≠ 0 is correct, and N is not in K.
- Case c = 0:
  - the derivation is correct;
  - the minimal polynomial `u^3 = -((20 + 11j)/19) u^2 + ((20 - 11j)/19) u + 1` was recomputed from
    `h^3 - 7h^2 + 11h - 1` (`h = w + 3`); it is irreducible over `Q(sqrt(-11))`, with residual 2.7e-51 at 50 digits;
  - the three coefficient equations and the λ formula check out.
- Nonvanishing: `λ^2 + d̄ ≠ 0`. Its absolute norm is `11^2 · 1639885141 · 1091860726096141 / (2^16 · 19^16)`, and
  `|λ^2 + d̄| >= 0.2647` at all 8 embeddings. The session's printed numerator matches exactly, and its proof-check script
  reproduces byte for byte.
- No case is missing: b = 0 forces a = 0, and κ = 0 is covered. The corollary follows given `k >= χ(L)`.

**B. Censuses, by two methods independent of the session's.**
- A float census of the full box at tolerance 1e-8 identified every candidate exactly in K:
  - e: 26 pure units;
  - ue: 26;
  - ve: 26, plus two near misses that are 3.306e-9 off 1 at 60 digits;
  - uv: 26 + 10;
  - e^2: 26 + 10.
- A modular sieve at primes 5281, 5581, 5641 agrees.
- The families are exact in K: uv is `q(0, 1, 1, d̄)` and e^2 is `q(d, 1, 0, d̄)`.
- The identities are symbolic: `u + v + d̄uv = uv·conj(d + u + v)`, `d + u + d̄e^2 = Zū/Z̄`, and
  `e^2 = Z^2 D/R = -Z/Z̄ ∈ F`.

**C. golden_forcing_search4.**
- The box census including d3 = 0 (9.17e8 differences) finds exactly 112 unit vectors:
  - axis units 26 / 22 / 22 / 22;
  - v-contact 10;
  - e^2-contact 10.
- The graph has 54,000 vertices and 389,980 edges, and coincides with the script's.
- A re-solve with a different encoding gives 0 of 570 blocked, with both shared and independent coefficient colourings.
  All 1,140 colourings were checked on every edge and port. So "still additively free" is justified for (A, C, C, C).

**D. Typing and credit.**
- PROVED is justified: the case analysis is complete and the computer step is exact and was reproduced.
- Novelty is fair. The whole-space statement is not in clock_duality's notes; they had only the finite source
  `G ∪ (P + eL)` and named this as the target. The only uncredited item is in finding 6.

Scripts: audit2_field.py, b1_proof_algebra.py, b2_degree192.py, b2b_counts.py, b3_census_numeric.py,
b3b_nearmiss_highprec.py, b4_census_modular.py, b5_forcing4_check.py, b6_identities.py,
b7_session_proof_check_rerun.out. Each runs in about 70 s or less.
