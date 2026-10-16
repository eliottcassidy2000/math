# S38 independent audit: THM-4659, THM-4660, the dense census, the pre-registration and the S38 note

Independent adversarial referee, 2026-10-16 (subagent of opus-2026-10-16-S38; resumed once after the session restart).

**Scope.** Worktree `/tmp/math-wt-explicit-modal-38`. THM-4659, THM-4660 and the note were checked as committed in
`c3bfa60c5c` (byte-identical to the working copies). The pre-registration was checked at `dc233fcf45`. Following the
coordinator's message, two items are deferred: the Q1-Q6 evaluation (Run B had not finished) and the irrational-grid fit.
No existing file was edited; no file with dense42b, irr38, dense43 or 2p43 in its name was opened.

**Summary.**
- The mathematics of THM-4659 and THM-4660 is correct, apart from minor fixes.
- The main problem is credit. The central content of both theorems, and the zero ρ*, already exist in the repository
  (E1, E2).
- The aliasing claims go too far (E3).
- The residue numbers in THM-4660's remark are wrong (E5).
- The per-block invariants of the census do not cover the moments and energy sums that later results depend on (E6).
- The 2^40 dense census passed every check made. The pre-registration is clean and reproduces exactly.

## ERRORS

**E1. MAJOR (credit and novelty; the mathematics is correct).** THM-4659 statements 2-3, the constant C and the zero ρ*
were already in the repository:
- **Codex note** `collatz_signed_resonance_pair_mass_20261013b.md` (2026-10-13): Eq. (9), `exp(-w) + exp(αw + 2πiθ) = 2`,
  is THM-4659's Φ(w, ε) = 0 with θ = ε/2, and the root is "continued from w(0)=0". It states the rule "q odd, p even".
  Its root table has the 332.5 line: frequency 332.50007587 and decay 5.47564541e-7 per octave. That decay equals
  ln 2 · (1 − Re ρ*), so this is ρ* = 0.99999921 + 3014.0202i. The 26.5 line's decay 0.00124127 equals
  ln 2 · (1 − Re s_0).
- **THM-4629 statement 3** (already PROVED): it gives w = 2πiθ_f/(2−θ) − 4π²θ_f²/(2−θ)³ + O(θ_f³). With θ_f = ε/2 and
  w₄₆₂₉ = ln 2 · w₄₆₅₉, this is exactly Im w = πε/κ and Re w = −Cε². It lists the parity-matched pairs (53,84),
  (159,252), (559,886), (665,1054), and contains the converse argument ("w → 0 forces both terms to 1").
- **HYP-9230** has the same coefficient as a Gaussian variance, 2/(2−θ)³.
- **HYP-9260 statement 4** has the pole distance (same constant). Its PROVED line bound |1/D(1+it)| ≤ θ²/(4θ_f²),
  combined with a bound on |D′|, gives statement 4 directly.
- **renewal_shapes_oct11's certificate** (`collatz_pole_certificate_20261013c.py`, lines 28-30) already finds s_0 by
  homotopy in the defect ε/2 of the same equation.

THM-4659 credits THM-4629 only for the selection rule and lists HYP-9260 only under `related`. The note cites none of
these and presents ρ* as an S38 finding ("the near-critical zero 1054/665"). Genuinely new: q-uniform neighbourhoods
with error orders O(ε⁴) and O(ε³); a converse that starts from Re ρ alone; statement 1; statement 4 as a corollary;
statement 5.

**E2. MODERATE (credit).** The identity F·D = G in THM-4660, and the conditional meromorphic continuation with residue
G(ρ)/D′(ρ), are Codex's P4 ((9a)-(9d) of the same note: identical E and R) and renewal_shapes_oct11's Lemma B1 /
Theorem B3. THM-4660 cites only "the identity behind Theorem B3", in `depends_on`. New: hypothesis (W), Lemmas A-C, and
part 2 (the Riesz-mean explicit formula).

**E3. MODERATE (aliasing; THM-4659 title and statement 5, note §3, pre-registration notes).**
- "Alias pairwise on every log-uniform grid with an integer number of points per octave" is false. A partner at
  q′ = 2Pk − q exists only when 2Pkθ is close to an even integer, because ε + ε′ = 2Pkθ − (p + p′). Counterexample:
  P = 1000; 2000θ = 3169.925 is 0.075 from an even integer, and none of the 23 family members has a partner. 1024 is a
  resonant density: 1024θ = 1623.0016.
- "Alias exactly" is wrong. All k = 1 pairs share the detuning (2048θ − 3246) · ln2/(2κ) = 0.00386 cycles per octave,
  25° of phase over [2^22, 2^40) and 2.8° over [2^40, 2^42), on top of different growth rates.
- "Fits that resolve the family need a non-integer number of points per octave" is false. A fit of M2 to a synthetic
  signal built from all 295 boundary residues:

  | Mode | Input | 1024/octave (observed) | 1000/octave | 1000√2/octave |
  |---|---:|---:|---:|---:|
  | ρ* (3014.0) | 0.00999 | 0.00625 (0.00712) | 0.01024 | 0.01000 |
  | 4455.5 | 0.01239 | 0.138 (0.151) | 0.01209 | 0.01284 |

  The integer density 1000 recovers the residues as well as the irrational one. Statement 5 is correct as a
  conditional identity once both zeros are family members.

**E4. MODERATE.** Note §7 and open thread 3 say "the unsmoothed sum needs |G| = o(|s|) on average along horizontal
lines". This is not enough: with |G| = O(|s|) the vertical Perron segment contributes x^(1/2+δ)·T; the formula also
needs ∫|G(σ′+it)|dt/|t| = O(T^ε) on the vertical line and horizontal bounds of order T^(1/2+δ) at T ≍ x; and
Σ r_ρ x^ρ/ρ is not absolutely convergent if |r_ρ| ≍ 1. THM-4660's own remark is correct but incomplete.

**E5. MODERATE.** THM-4660's remark says the estimates are "bounded (by about 5)" and "about 0.005-0.02 on the Mercator
family". The maximum |r̂| is 11.05, at 0.501924 + 16556.222i. On the family the range is 0.0052-0.184, and 9 of 23
members fall outside [0.005, 0.02] (s_0 0.184, 720.74 0.070, 5307.3 0.057, 8561.6 0.047, 2533.5 0.037, and four more).
"The numerics beat [O(|ρ|)] by far" is not supported: at 2^40 under the (2,0) budget only 11 of 295 zeros have |r̂|
larger than the radius; the median radius/|r̂| is about 47,000 for Re ≤ 0.6.

**E6. MODERATE (reliability).** Single-bit faults injected into `block_invariants`: detected for class counts, the
register total, Σt per sub-block and the H increment (20 of 20 each); never detected for Σ(n−lo)t (which feeds the
residues and the Arb certificates) and hs1/hs2 (prefix sums and energies used by Q3-Q5) (0 of 20); the 16 GB parity
table and the grid pass are not covered. This machine produced two bit flips in one run. Run B is a single run, so it
should not be typed FINITE-EXACT until it is replicated. For [0, 2^40) the question is closed by the checks in the
verdicts table.

**E7. MINOR (typing).** The 295-zero count comes from a floating-point argument principle: NUMERICAL, not FINITE-EXACT.
The count itself is right.

**E8. MINOR (numbers).** "199.17 at 1054/665" should be 199.1638, below C = 199.1647 (the ε⁴ coefficient is +206120.8;
the 199.17 comes from rounding 1 − Re to 7.9e-7). The near-1/2 zero at 1527.9 has Re = 0.500806, not "0.5006". The
almost-period median is 0.00061, not 0.0012. The pre-registration says "ten largest low modes" but lists five.

**E9. MINOR.** "This verifies the S37 2^40 census on [2^37, 2^40)" overstates it: the rerun is the same code, so it is a
replication, not an independent check.

**E10. MINOR.** The good-height step in Lemma C ignores zeros just outside [1/4, 3] × [m, m+1]; zeros in [1/4, 3] can
be as close as 0.0018 in height. Fix: count over the enlarged rectangle (Jensen ratio 11.34 at radius 1.61) and halve η_K.

**E11. MINOR.** "|G(s)| ≤ K′|s| on Re s ≥ 1/2 + δ" holds only on bounded strips, because of the factor (3/2)^σ. The
proof only uses σ ≤ 2.

**E12. MINOR.** THM-4660 cites "the two pre-registered tests" as supporting evidence, but the S38 test was still pending.

**E13. MINOR.** The ρ* certificate "PENDING at 2^42 and 2^43": at 2^42 the radius is 0.0130 > |r̂| = 0.00999, so it
cannot certify there; at 2^43 the radius is 0.0092, a margin of about 8%.

**E14. MINOR (reproducibility and disclosure).** `dense40.{npz,json,_state.json}` and `cert_rhostar_2p40.json` are not
committed. `Q2_bound` and `Q6_bound` were added to the prediction JSON by hand; the committed script does not write
them. Run B was relaunched at 14:20:09 with the post-registration script edit; the evaluation grid is provably
unchanged, but this should be disclosed.

**E15. MINOR.** `--check-SH` only compares the state file's own scalars with the command-line values; it does not
recompute them from the counts and does not abort. Recomputed by the referee: they match.

## IMPRECISIONS
- **I1.** The ε-series of w(ε) has radius of convergence 0.0175005; the branch point is at ε = −0.0175005i, where the
  zero merges with the one coming from Re s = 2. w = 10.92036iε − 199.16466ε² − 5757.14iε³ + 206120.8ε⁴ + …. The 8/5
  member lies about 4× beyond the radius; continuing along the real ε axis still lands on 0.693918294 + 22.253194573i.
- **I2.** The "±240.25" sidebands sit around 2·Im ρ*, i.e. (2·1054 ± 84, 2·665 ± 53).
- **I3.** "At most one zero per unit height" is true for 1/2 < Re < 1 (minimum gap 44.5) but not for Lemma C's
  rectangles, which can hold 2.
- **I4.** |D′(ρ)| depends only on Re ρ. Its minimum over 1/2 < Re < 1 is κ/2 = 0.143841 (approached as Re → 1); over
  0 ≤ Re ≤ 3 the infimum is 0.1308 (as Re → 2⁺).
- **I5.** The 7.6% mismatch at 720.74 is also aliasing (partner 8561.6; synthetic fit 0.0656). The drop in M2's
  holdout error from 0.0118 to 0.00167 is absorption of aliased power by the near-Nyquist zeros, not resolution of the
  high modes; the synthetic test reproduces the drop (0.0110 to 0.00122).
- **I6.** M2 was picked from 12 candidates on the holdout itself; the resulting bias is negligible (0.001670-0.001677).
- **I7.** In §4, the tail bound at σ = 0.85 first exceeds the numerator at t ≈ 528 on the grid, not 200 (it does so at
  96% of grid t ≥ 200).
- **I8.** Further missing credit: clock_duality_oct15's W1/W2 use the same function (W1's γ₁(σ) = f(σ), contracting
  exactly on (1,2); W2 uses f′(1) = −κ/2). The boundary estimator is renewal_shapes_oct11's formula (4), not an "S37
  certificate formula". The note does not cite THM-4632 or HYP-9260.
- **I9.** HYP-9283, HYP-9284 and HYP-9285 had no S38 edits; HYP-9285 still showed the S37 state.

## VERDICTS

| Claim | Verdict | Evidence |
|---|---|---|
| THM-4659 st.1 (no zeros for 1 ≤ Re ≤ 2) | CONFIRMED | proof checked; torus enumeration finds none |
| Identity D(1+iπq/L+w) = Φ(w,ε) | CONFIRMED | agrees to 2.3e-57 |
| St.2 (IFT and expansion) | CONFIRMED (credit, E1) | a₁ = iπ/κ and a₂ = −C to 1e-41 |
| St.3 (converse) | CONFIRMED | all 23 zeros with Re > 0.98 below 20000 are ρ_(p,q), p even, q odd; distance ≤ 5e-12 |
| St.4 (zero-free region) | CONFIRMED | also follows from HYP-9260 |
| St.5 and title (aliasing) | CONFIRMED WITH FIXES | E3 |
| Quoted zeros | CONFIRMED | 199.17 → 199.164 |
| F·D = G | CONFIRMED | exact cut-off identity to 1.4e-12 at 24 (s, N) pairs |
| Lemmas A, B | CONFIRMED | |
| Lemma C | CONFIRMED WITH FIXES | E10 |
| THM-4660 part 1 | CONFIRMED WITH FIXES | E11 |
| THM-4660 part 2 | CONFIRMED | arcs, separation and Riesz kernel checked |
| (W) weaker than HYP-9283 budgets | CONFIRMED | |
| Unsmoothed sum | the note's characterization REFUTED | E4 |
| THM-4660 remark (residue sizes) | REFUTED | E5 |
| Resume logic | CONFIRMED | fresh = split run at 2^24→2^26 = own census |
| Agreement with the S37 64-point census | CONFIRMED | 1,920 points; moments array and octave sums identical |
| Independent census below 2^34 | CONFIRMED | all 24,489 grid points, 131,072 moment rows, 14 octave sums; S(2^34) = 24953018, H(2^34) = 26869 |
| Above 2^34 | CONFIRMED | 96 random sub-blocks; 18 dense-only points; 20,008 random n in [2^40, 2^42) plus 836 record-trajectory points; S38 kernel on 1,796 values passing 2^64; 0 mismatches |
| Reliability story | CONFIRMED WITH FIXES | E6 |
| Pre-registration precedes Run B | CONFIRMED | commit 11:43:32/35, push 11:43:41, first launch 11:43:48, relaunch 14:20:09; no Run B data file existed; pre-registration files unchanged |
| Holdout table and M1/M2 predictions | CONFIRMED | exact reproduction |
| Q1-Q6 evaluation | NOT CHECKED | deferred |
| 295 zeros | CONFIRMED (NUMERICAL) | independent lattice-point count |
| Residues reproduced | CONFIRMED | to 8.6e-11; meaningful only for Re ≳ 0.7 |
| Stability of r̂ at ρ* | CONFIRMED | 0.009969-0.009995 over N = 2^32 to 2^40 |
| Boundary vs least-squares residues | CONFIRMED | reproduced exactly |
| Aliasing explanation | CONFIRMED WITH FIXES | E3 |
| Arb radius scaling | CONFIRMED | ∝ N^(−1/2): 0.0261, 0.0130, 0.0092, 0.0065 at 2^40, 2^42, 2^43, 2^44 |
| Irrational-grid fit | NOT CHECKED | deferred; the synthetic test supports the method |
| S38 note | CONFIRMED WITH FIXES | E1-E9, E12-E14 |

## Scripts (this folder; each with a matching .out and .json)
a1_thm4659_family.py (identity, Taylor coefficients, branch point, quoted zeros, the 23 family members, alias pairs for
P = 1024 and P = 1000, convergents); a2_zeros_torus_count.py (independent lattice-point enumeration of the zeros, Lemma
A/C numbers, near-1/2 zeros, almost-period; + a2_zeros_half_one.npy); a3_thm4660_identity_constants.py (exact cut-off
identity, Lemma C constants, Riesz kernel); a4_independent_census34.py (own census to 2^34, random sub-blocks, random n;
+ a4_census34_grid.npz); a4b_s38_kernel_resume_invariants.py (the S38 kernel on values passing 2^64, fault injection,
resume test); a4c_dense_only_points.py; a5_census_replication.py; a6_boundary_residues.py; a7_prereg_and_alias.py
(pre-registration reproduction and the synthetic aliasing test); a8_misc_checks.py.
