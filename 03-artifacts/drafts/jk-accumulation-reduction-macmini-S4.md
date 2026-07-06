# The J-K accumulation reduction of (G) — exact statements, normalization, and the finiteness program

**mac-mini-2026-07-06-S4 (HYP-4262).** Full-text extraction of Jain–Kravitz
(arXiv:2411.12684, "Relative Lonely Runner spectra") + the Giri–Kravitz
no-accumulation-from-below theorem, translated to the fleet's normalization;
the (G) tail program made precise. Species test:
`05-knowledge/results/lrc_double_lift_species_macmini_S4.out`.

## The dictionary

J-K work with D(T) = L∞-distance from a subtorus T ⊆ (ℝ/ℤ)ⁿ to the center
(1/2, …, 1/2). For speeds v₁,…,vₙ (nonzero integers):
**M(W) = 1/2 − D(⟨(v₁,…,vₙ)⟩_ℝ)**. Our objects at n = 12:

- the 12-spectrum = 1/2 − 𝒮₁(12), where 𝒮₁(n) = D-values of PROPER 1-dim
  subtori (no zero coordinate). Repeated-|coordinate| subtori reduce to
  fewer distinct speeds (values ≥ 1/12 in M — above our window), so inside
  the gap window the spectrum IS the distinct-12-set spectrum.
- our gap M ∈ (1/13, 2/25) ⟺ D ∈ (21/50, 11/26). Directions flip:
  "accumulation from below" in D = "from above" in M.

## The three citable inputs (exact statements as extracted)

1. **[G-K] acc⁻(𝒮₁(n)) = ∅** (Giri–Kravitz, cited as [8] "Giri and the
   second author"): no accumulation from below in D ⟹ **in M: no sequence
   of distinct spectrum values decreases to a limit.**
2. **[J-K eq. (4)] acc(𝒮₁(n)) ⊆ 𝒮₂(n)**: accumulation points are D-values
   of proper 2-dim subtori. With **[eq. (5)]**: near d ∈ acc(𝒮₁(n)),
   𝒮₁(n) ∩ (d, d+ε) = (⋃_{j=1}^t 𝒮₁(U_j)) ∩ (d, d+ε) for FINITELY many
   2-dim U_j with D(U_j) = d.
3. **[J-K Thm 1.1]** Each relative spectrum 𝒮₁(U) has finite symmetric
   difference with a finite union ⋃ᵢ (D(U) + 1/Prog(αᵢ, βᵢ)),
   Prog(α,β) = {αs + β : s ∈ ℤ≥0}; the αᵢ, βᵢ are algorithmically
   computable ("quite explicit... an algorithm, albeit not so simple");
   the finite symmetric difference "is itself a finite calculation."
   Their Section 3 parametrizes candidate 2-subtori by SATURATED rank-2
   lattices: u, v ∈ ℤⁿ with ⟨u,v⟩_ℝ ∩ ℤⁿ = ⟨u,v⟩_ℤ.

## The two unconditional freebies

**(F1) The bottom endpoint is safe.** Accumulation at 1/13 from inside the
gap would be M_n ↓ 1/13 with M_n > 1/13 — a decreasing sequence of distinct
values, impossible by [G-K]. No stability-rate analysis is needed for the
1/13 end of the tail. (klein's (U) — uniqueness of the attainer AT 1/13 —
remains a separate, still-needed input to hdich.)

**(F2) Product 2-tori are safe.** U = T_A × T_B (blocks A ⊔ B, independent
times) has M(U) = min(M(A), M(B)) ≥ 1/12 > 2/25 by the LRC(≤ 12) citation
(each block ≤ 11 speeds). One line. Dangerous 2-tori are COUPLED — the
lift-family limits (v_i = r_i + N·ℓ_i, N → ∞ along a ray).

## The reduction

Suppose (G) fails infinitely: infinitely many distinct spectrum values in
(1/13, 2/25). They accumulate at some μ ∈ [1/13, 2/25]; μ ≠ 1/13 by (F1);
by [eq. (4)] μ = M(U) for a proper 2-subtorus, coupled by (F2), and by
[eq. (5)] the gap values near μ lie in finitely many relative spectra
𝒮₁(U_j), each a reciprocal progression climbing to μ from below (up to
finite Δ). Hence:

> **(G) ⟺ (A) no coupled proper 2-subtorus of (ℝ/ℤ)¹² has
> M(U) ∈ (1/13, 2/25], PLUS (C) a FINITE census of 1-dim gap values.**
>
> More precisely: (A) kills all accumulation in the half-open gap; then
> spectrum ∩ (1/13, 2/25) is finite outright (bounded witness denominator
> q), and (C) each residual value is a Farey cell (c, q) killed by the
> anchored machinery (HYP-4232/4242/4252: unit pairs, witness determinism,
> k-stratification, gcd-ladder).
> If instead some U has M(U) = 2/25 exactly, insert (B): compute 𝒒₁(U) by
> J-K's finite procedure and check its gap side is empty/finite-Δ.

## The species test (new data, this session)

The natural candidate for a 2-torus AT 2/25 is the double-lift direction of
the attainer: W_s = {1,2,3,5,7,8,9,10,11,12} ∪ {4+13s, 6+13s}. EXACT
computation, s = 1..30: M_1 = 2/25 (the attainer), then the species
RETREATS: M_s ≥ 2/19 for every s ≥ 2, dominant value 2/17 (19 of 30),
dips only to [1/9, 2/17] at CRT-structured exceptional s (the t = 4/17
witness breaks at s ≡ 15,16,0 mod 17, handled by neighboring witnesses).
**Zero in-gap values.** So the double-lift 2-torus sits near 2/17, NOT at
2/25: the attainer is an isolated small-s exception of its own species —
exactly J-K's finite-symmetric-difference phenomenology — and no
progression climbs through the gap toward 2/25 along this species.
Step (B) is plausibly VACUOUS: the program is (A) + (C).

## What (A) actually asks (the lift-floor program, reframed)

Enumerate coupled saturated rank-2 lattices ⟨u, v⟩ (12 nonzero-projection
coordinates, up to the symmetries); U's in the window need
M(U) = max_{(t,s)} min_i ||u_i t + v_i s|| ≤ 2/25 — a TWO-parameter
no-2/25 covering system. The fleet's tools apply verbatim one level up:
dist-0 pinning at 2..12 (11 divisor demands on the PAIRS (u_i, v_i):
at (t,s) = (j/q', j'/q') grids), the cluster-gcd ladder on the u-side,
and lift_floor_beta_ladder ≥ 2/25 is exactly the first proved stratum
(single-coordinate coupling ℓ = e_i along the THM-621 ladder rays).
The HYP-4109-series strata are the remaining named surface.

## Honest caveats (for canon discipline)

- J-K (2411.12684) and G-K ([8]) are PREPRINTS about the spectrum at
  general n — NOT covered by the owner's LRC(≤13) take-as-true policy
  (which is about ≤13-runner bound statements). Citing (F1)/[eq. 4]/[eq. 5]
  in the final assembly needs either the standard preprint caveat or
  self-contained re-proofs (Chabauty compactness of the space of closed
  subgroups + upper structure of limits of 1-subtori; the G-K
  no-below-accumulation is the deeper input).
- The exact hypotheses of eq. (4)/(5) (properness, the ε-uniformity, the
  𝒮₃ refinement at points of acc(𝒮₁) ∩ 𝒮₃ — Corollary 1.2 excludes
  𝒮₃(n) points) must be re-verified against the full PDF when this enters
  the DAG; the window's 𝒮₃ values are plausibly ∅ (3-block products are
  ≥ 1/11 in M; coupled rank-3 = double-lift towers, expected ≥ 2/25 by
  the same floor program) — UNVERIFIED, flagged.
- The Giri–Kravitz full citation must be pulled from the PDF bibliography.

## Status

Reduction: STATED, citation-conditional (preprint caveats above).
(F2) + species test: verified this session. (A)-classification and the
q-bound extraction: the named next steps; (A) is the lift-floor program's
new coordinates, (C) is the anchored-census machinery already running.

-> HYP-4252 (uniform cell lemma: the (C) machinery), HYP-4242/4232 (the
3/38 instance), HYP-4212 (lift_floor_beta_ladder = (A)'s first stratum),
HYP-4167 (kps: the J-K lead), THM-621/622, OPEN-Q: the (A) enumeration;
the self-contained re-proof of [G-K]/[eq. 4] if preprint caveat refused.
