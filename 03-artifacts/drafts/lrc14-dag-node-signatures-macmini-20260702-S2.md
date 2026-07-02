# LRC(14) DAG-node signatures — the plug-in spec for the remaining censuses and transcription

**mac-mini-2026-07-02-S2 (HYP-3859).** For each remaining [CENSUS]/[TABLE] node of the
proof-DAG ledger: the exact statement its run must output, phrased as the Lean hypothesis
it discharges. Runs conforming to these signatures wire into `LRCFourteenSkeleton` with no
further mathematics.

1. **simplex_minima_table** (this session: the 30-pattern list; minima known d=2, computable d≥3):
   `∀ pattern p ∈ SIMPLEX30, minOverlap p = table p` — per-row `decide`.
   Discharges: THM-599's ε_d via THM-602(C)'s partial-cycle accounting.
2. **lattice_census** (rank-1 lists done j≤4; full HNF nested loop spec'd):
   `∀ Λ ∈ LATTICES(j≤13, h≤7), effectivePattern Λ ∈ BOUNDED_PATTERNS` + per-pattern floor.
   Discharges: THM-602(A) case 2–3 base cases (with opus F_7..F_11 rows as instances).
3. **klein_nest_census** (queued): `∀ P ⊆ {1..13}, |P| ≥ 3 ⟹ |∩D_P| = 2r / max(P/gcd P)`
   — klein's gcd-nest law as a per-subset rational identity (8100 rows, verified; needs
   the mediant-criterion proof text or per-row decide).
   Discharges: hp0cap's cap arithmetic (Λ_P(1/14) closed form).
4. **klein_Ly_comparison** (queued): `∀ shape E ∈ SPREAD30, L_y(E) ≤ L_y(consec_k(E))`
   — finite rational comparisons via nest closed forms.  Discharges: hp0cap via THM-534.
5. **kps_star_census** (queued): the (★)-census rows of HYP-3953 as exact rationals via
   the torus-band volumes (THM-600).  Discharges: hpartA's residual finite census.
6. **commensuration** (PROVED this session, THM-602 addendum): `7 | Q ⟹ ov_{P,Q} ≡ (2r)²`
   — transcription target: a finite exact-tiling argument, polygon-module species.
7. **forced_independence_d5** (THM-602(C), paper-complete): transcription target; consumes
   rows 1 and 6; produces the S_d hypotheses of `BonferroniTruncation`-driven floors.
8. **skeleton wiring**: replace `hpartA` by [rows 2+5+7 assembled per THM-602(D)];
   replace `hp0cap` by [rows 3+4 + THM-534]. Both replacements are glue, no mathematics.
