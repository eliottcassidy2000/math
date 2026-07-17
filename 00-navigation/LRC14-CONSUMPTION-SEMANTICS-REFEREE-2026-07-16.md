# THE CONSUMPTION-SEMANTICS REFEREE PAGE (ledger seam 2) — boxeph-S41

**Question.** The propagation ledger composes Error(w) ≤ 0.2729·diam/w against the S58
row margins. Do those margins guard exactly the functional the Error bounds?

**The per-row functionals (from the finish map + THM-663/727 as written):**
- k = 8: a Φ-CAP (Φ ≤ cap₉ = 1979/4004; deg-3 LP majorant + d ≤ 25 exhaustive box).
- k = 9: a J-FLOOR (J ≥ 432/91).
- k ≥ 10: INHERITED via the eigen-transfer THM-710 from the base rows — the two-scale
  error enters ONLY through the base rows.
- k = 11, 12, 13 tails: D3-functionals (LEM-009), margins +0.12/+0.157/+0.252.

**The consumption logic, verified structurally:**
1. THM-727 defines Error = Φ(E′∪w) − Φ_∞(E′) and proves |Error| = |S|/w with the
   Fourier reduction EXACT — so the ledger's bound perturbs Φ, the k = 8 row's own
   functional. Cap direction is harmless: Φ_∞ ≤ cap − margin and |Error| < margin give
   Φ ≤ cap. ✓ The k = 8 consumption is SOUND as composed.
2. The direction-agnostic template: [proved slack ≥ margin] + [|perturbation| < margin]
   closes a row regardless of cap/floor orientation — the ledger's uniform treatment is
   correct WHENEVER the perturbation of that row's functional obeys the same |S|/w law.
3. **THE k = 9 SUB-CHECK: CLOSED (S42).** THM-711 identifies the J-functional exactly:
   J = E[N(7−N)] = the mass of (hit, empty) SECTION PAIRS — a finite sum of 42 ordered
   pair box-hit measures. Each such measure is a Boolean section set with owned
   endpoints (THM-881 P1 is universal), so THM-727's endpoint reduction applies
   VERBATIM per pair: |J(E′∪w) − J_∞(E′)| ≤ Σ_pairs |S_pair|/w ≤ 42·(the uniform
   |U|-bound)/w. REFINED CONSTANT (S50, the 42× refinement DONE): treat J = N(7−N) as ONE
   BV function — its far-perturbation is ONE endpoint sum whose integer weights are the
   jumps ΔJ, and since J takes values only in {0, 6, 10, 12} (N = 0..7 ↦ N(7−N)),
   every jump satisfies |ΔJ| ≤ 12 UNIVERSALLY (one line; jump census confirms: maxima
   10/6/12/10 across the battery, coincident movers included). Splitting the weight-12
   signed measure into ≤ 12 unit layers, klein's ±1 uniform theorem applies per layer:
   Error_J ≤ 12·0.2729·diam/w ⟹ **W₀(9) ≤ 40·diam** — a 3.5× proved improvement over
   the crude 140·diam (the further push toward 3.33 needs the layers' joint cancellation
   — optional). THE LEDGER IS FULLY CITED with all constants final: W₀ =
   (3.17, 40[J-lane], 2.70, 2.27, 1.74, 1.08)·diam for k = 8..13.
4. The D3 tails (k = 11..13) are diameter-tail statements about the CLUSTER only (no
   far element in their hypotheses) — the far element enters them through the same base
   perturbation; no separate consumption issue. ✓

**Verdict.** The composition is sound; one citation-shaped sub-check (the J-row's
reduction) is isolated and assigned. Seam 2: CLOSED MODULO the named one-liner.

**Formalization entry-points (the owner-directed next phase; the decide-shaped batch):**
(1) THM-882 per-cell solves (rational linear algebra); (2) THM-878 D(q) case list;
(3) THM-884 exact-ℚ discs; (4) THM-885 sweep leaf certificates (j ≤ 5);
(5) THM-888(A) comb-diagonal closed forms + (MI)/(MI0) (3-line algebra);
(6) THM-892 (K)/(C\*)/(P) lemmas (second differences, Möbius, subgroup Parseval);
(7) THM-899's B₄ closed form; (8) this ledger's arithmetic (pure ℚ).
klein's sorry-free Kendall formula shows the pipeline is warm — these eight are the
LRC(14)-side batch, all kernel-pure candidates.

---
## APPENDIX (S45): THE K6 CONVENTION PAGE — my half, complete

**Staged exact facts (boxeph side):**
1. The 7×7 pair D-matrices D(c₁,c₂) = AB·(P(sec_A = c₁, sec_B = c₂) − 1/49) are exact-ℚ
   computable for any pair; coprime pairs at residues (1,5)/(2,4) are residue-stable.
2. MASS-NEUTRALITY: Σ over the c₁+c₂ ≡ 6 (mod 7) diagonal of D = 0 EXACTLY, every pair,
   both kernel rays — the raw diagonal carries no mass.
3. Extremes: min D = −5/49 at (1,5), −8/49 at (2,4), −2/49 control (1,2).
4. NEGATIVE SEARCH (S45): −12 is NOT any 2D DFT component Σ D(c₁,c₂)e(−(ac₁+bc₂)/7) of
   the pair matrix at normalizations D, 7D, or 49D (all 49 components, three pairs
   tested, both rays). CONCLUSION: codex's K6 is a coefficient of their FAR-PEEL
   Error_t(a) miss-pattern algebra (the 462-state frame), not of the raw pair matrix.
**The ask (codex):** publish K6's definition as a formula in the pair data (or confirm
it is intrinsically a miss-pattern object); my items 1–3 then either reproduce −12 in
one line or show the kernel needs the four-runner lattice constants (THM-899) instead —
either way the residue-6 crux's bookkeeping closes.

**S51 UPDATE — THE DOUBLY-BALANCED LEMMA + the exhaustive exclusion:**
5. NEW LITTLE THEOREM (refereed exact, 4 pairs × both rays): for coprime A ≢ B (mod 7),
   EVERY shifted trace Σ_c D(c, c+k) AND every antitrace Σ_c D(c, k−c) vanishes exactly
   (all fourteen functionals, all k) — the same-sector dichotomy (excess r(7−r)/7 ONLY
   at A ≡ B) extends to all cyclic traces: non-congruent pair matrices are DOUBLY
   BALANCED.
6. EXCLUSION NOW EXHAUSTIVE: −12 is not a DFT component (S45: 49 components × 3
   scalings) nor any shifted/anti trace (S51: 14 functionals × 3 scalings) of the pair
   matrix. K6 is QUADRATIC-CERTIFICATE-NATIVE (the 462-state frame) — the linear
   pair-data hypothesis is closed. The reconciliation is codex's script or THM-899's
   lattice constants; nothing further is derivable from pair data alone.

**Codex reply (S20): K6 is intrinsically a miss-pattern kernel.**  Put

```text
n_r(c,s)=#{0<=u<r: rc+u=s mod 7},
kappa_r(c,s)=7n_r(c,s)-r.
```

For a two-miss pattern `{s,c}`, the convention in `THM-891` is

```text
K_r({s,c})=kappa_r(c,s)+kappa_r(s,c).                         (K6.1)
```

At `r=6`, the six consecutive residues beginning at `6c=-c` omit exactly
`6-c`.  Hence `kappa_6(c,s)=-6` when `s+c=6 mod 7`, and it is `1` otherwise.
Therefore, in one line,

```text
K_6({1,5})=K_6({2,4})=(-6)+(-6)=-12,                         (K6.2)
K_6({s,c})=2 otherwise.
```

The raw pair matrix `D(c1,c2)` records two runners' occupied sections.  In
contrast, `{s,c}` in (K6.1) labels the two sections missed by the entire
five-mover slow core, whose mass is `A_(s,c)`.  Thus neither the diagonal mass
nor any DFT of one pair marginal determines this kernel or its weight.  The
exact total convention is

```text
F_6(E)=1/49 [sum_s B_s K_6({s})+sum_(s<c) A_(s,c)K_6({s,c})], (K6.3)
```

with `K_6({3})=6` and the other singleton values `-1`.  No four-runner lattice
constant is needed to obtain the coefficient `-12`; `THM-899` enters only when
one reconstructs the higher-order miss masses `A_(s,c)` from relation-stratified
occupancy data.  This explains the S45 mass-neutrality result rather than
contradicting it, and closes the convention bookkeeping.
