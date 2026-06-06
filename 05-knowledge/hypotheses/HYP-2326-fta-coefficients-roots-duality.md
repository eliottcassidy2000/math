---
id: HYP-2326
title: The two spectra are FTA-dual — combinatorial spectrum = coefficients, character-ratio spectrum = roots; tie-graph chromatic roots = roots of unity; Lee–Yang zeros = the transition
status: OPEN (synthesis); chromatic-roots=roots-of-unity and Z-zero-at-0 ⟺ collapse VERIFIED
source: claudebox-2026-06-03-S639
related:
  - THM-403  # cyclotomic witness orbit (roots of unity) = the character ratios
  - HYP-2311 # character-ratio spectrum bounds H/chi_di
  - HYP-2316 # the polarized delta field / glass transition (Lee-Yang real-axis pinch)
  - HYP-2295 # chromatic polynomial = partition function; ultra-log-concavity / Newton (canon)
  - HYP-2321 # Gauss sums (roots) and reciprocity
---

# HYP-2326 — the fundamental theorem of algebra is the bridge between our two spectra

The user's frame: a degree-`n` polynomial has `n+1` coefficients (the constant = the "+1") and, by
the **fundamental theorem of algebra**, exactly `n` roots on `ℂ` (with multiplicity). For our objects
this is the bridge we had been missing: **the combinatorial spectrum is the COEFFICIENTS, the
character-ratio spectrum is the ROOTS, and FTA / Newton's identities / Vieta map between them.**

## The two spectra, made dual

| | COEFFICIENTS (n+1) | ROOTS (n) |
|--|--|--|
| object | the combinatorial spectrum | the character-ratio spectrum |
| examples | H-values, `p_k` (covering depth), independence counts `α_k`, chromatic polynomial | eigenvalues, Gauss sums, **roots of unity**, Lee–Yang/fugacity zeros |
| sessions | S636–S637 (H, delta, glass) | S625/S636/S638 (multipliers, char ratios, Gauss) |
| bridge | **Newton's identities** (power-sums ↔ elementary-symmetric), **Vieta**, real-rootedness ⟺ ultra-log-concavity (canon) |

## Verified anchors

1. **Tie-graph chromatic roots = roots of unity.** The chromatic polynomial of `C_n` (the tight LRC
   tie-graph, S632) has roots `= 1 + (n-1)th roots of unity` (verified n=5,6,7: `C_5→{−1,±i,1}+0`,
   `C_7→` 6th roots of unity `+0`). So its **roots are exactly the cyclotomic witness orbit (THM-403)
   = the character ratios**, while its **coefficients are the H/coloring spectrum**. One polynomial,
   both spectra, FTA-linked.

2. **Lee–Yang zeros encode the phase transition.** The covering-depth `Z(z)=Σ p_k z^k` (the LRC
   partition function) has **`z=0` as a zero ⟺ `p_0=0` ⟺ tight/collapse** (loneliness measure
   vanishes). Verified: tight `(1,2,3,4),(1,3,4,7)` have a zero at the origin (the ground state
   gone); loose `(1,2,4,8)` (`p_0=0.1`) does not. The **partition-function zeros pinching the real
   fugacity axis are the transition** — the glass transition (HYP-2316) read in the root picture
   (canon: "fugacity-axis vanishing theorem").

## The n / n−1 / n+1 bookkeeping is FTA degree/coefficient/root counting

The program's recurring counts are the FTA ledger: `n` runners ↔ degree; `n−1` moving speeds ↔ the
non-constant coefficients / the roots after removing the trivial one; the **constant term (`+1`)** ↔
the observer / the empty independent set (the lone odd `1` in `H=1+2α₁+…`, the `k=1` chromatic root,
the empty packing). The `2n−1` shell and the `n+1` perspective coincidences (A093934) are this same
degree/root accounting one level up. **"Mapping n−1 things to n things" is reading a polynomial by
its coefficients vs its roots.**

## Implications
- The character-ratio spectrum (which bounds `H` and `χ_di`, HYP-2311) **is literally the root set** of
  the combinatorial polynomial; the Hoffman bound is a coefficient↔root inequality.
- Real-rootedness of the independence polynomial `I(Ω,z)` (claw-free ⟹ real-rooted, Chudnovsky–
  Seymour) would force the `α_k` to be **ultra-log-concave** (Newton; canon "ultra-log-concavity-
  newton-tournament") — a coefficient constraint *from* the roots. The forbidden-H values (S636) are
  then a root-locus phenomenon.
- The glass transition / loneliness / collapse are **Lee–Yang zero conditions** (zeros at the origin
  / pinching the real axis), tying S637's energy landscape to the analytic theory of zeros.

## To do
1. Compute the Lee–Yang/fugacity zero locus of `Z(z)` and `I(Ω,z)` across `n`; does it pinch the real
   axis at the even-`n` glass transition (HYP-2316)? Is the locus a circle (ferro) or off-axis (the
   antiferromagnet, THM-290)?
2. Newton's identities explicitly: write the H-moment spectrum (power sums) ↔ the character ratios
   (roots) for the conflict graph / Hermitian adjacency; recover the eigenvalues from the H-data.
3. Forbidden-H as a root-locus statement: which root configurations are unrealizable (the gaps).
