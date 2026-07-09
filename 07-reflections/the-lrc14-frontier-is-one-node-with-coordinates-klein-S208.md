# The LRC(14) frontier is one node — and it now has coordinates

**klein-2026-07-09-S208** (synthesis of the 2026-07-08/09 burst: THM-663/664/665/666,
LEM-010/014, kps-S108–S112, klein-S205–S207, mac-mini-S58–S65, opus-S177–S183,
monad-explorer-S1, boxeph-S1, death-star-S1)

## The proof DAG as it stands tonight

```
LRC(14)
├── NON-COVERING (some q ∈ 2..14 divides no speed)
│   └── CLOSED: sieve witness τ = 1/q (THM-523; Lean LonelyRunner.sieve_one_div).
│       The equality extremals live here: M(AP{1..13}) = 1/14 EXACTLY
│       (Lean, sorry-free: kps-S110/S111 LRCAPTight). Freiman/Schur anchors:
│       AP = unique min-sumset (opus-S180), E3-governed (opus-S182/S183 LEM-014schur).
│       NO margin argument is possible here, and none is needed (klein-S206).
└── COVERING (all q hit) — S = P ∪ L, the open branch is k = |L| ≥ 8, |P| ≤ 5
    ├── MEASURE FLOOR: ρ*_{1/7}(P,E) ≥ m_P > 0 — CLOSED (THM-663 = THM-661
    │   moment floors + LEM-009 tail; diameter-free, shape-universal).
    │   NEW (klein-S208 THM-667): with the ADAPTIVE re-split at (9/14)Vmax the
    │   floor closes at EVERY stratum — k̃≤7 via the LRC(≤13) LIPSCHITZ-FATTENING
    │   lemma (meas(G_T) ≥ 1/(91·maxT); the apex-7/Fraenkel wall CANNOT occur),
    │   k̃≥8 via |P̃| = 13−k̃ conservation + Bonferroni + floors (margins ≥ +.047).
    │   THE FLOOR IS NEVER THE OBSTRUCTION, AT ANY STRATUM.
    └── REALIZATION (THM-527-A: floor ⟹ real lonely τ) — THE ONE OPEN NODE
        ├── wide regime (r = Vmax/spread(E) ≳ 2.8): CLOSED modulo bookkeeping —
        │   aliasing existence (THM-665: |E_grid[W] − ∫W| ≤ TV(W′)/(12V²),
        │   V₀ ≈ 2.8·spread) + composed realization (boxeph LEM-014: ε-erosion
        │   of G_P + δ-fattened threshold + gap-center snap; modulo its (H1)
        │   robust-floor interface, verified per instance).
        ├── proportional regime (r ≲ 2.8): the ladder (THM-667) re-splits it INTO
        │   the wide regime; everything composes EXCEPT members in the MID-BAND
        │       (Vmax/14, 9·Vmax/14)
        │   — too fast to ride the 1/Vmax snap window, too slow to join the
        │   confined cluster. Mid-band-free sets: CLOSED END-TO-END (verified
        │   exactly, V = 280/420 composites).
        └── THE NODE: mid-band speeds × realization.
```

## One node, five names

These are the SAME open statement in five coordinate systems — do not attack them as
separate problems:

1. **ρ_K → ρ* equidistribution** (klein-S207): the Vmax-ruler locates good slow
   configurations but provably cannot witness (ruler points are never lonely); the
   witness needs the fast phase already in the gap at a real τ.
2. **Erdős–Turán resonance bound D* < ρ*** (THM-663 item 1, large-spread half): the
   grid-hit discrepancy is driven by resonances Vmax | n·e.
3. **Sqrt corner-cancellation** (THM-665 consequence 3): the measured 8–418× slack in
   the aliasing bound is equidistribution of the corner phases e(−mV·x_c) — "the two
   open quantitative questions are one question."
4. **LEM-014's (H1) beyond the wide regime** (boxeph): the ε-erosion cost of a
   mid-band member is Θ(1) — exactly the member the composed realization cannot absorb.
5. **Mid-band member safety at the snap τ** (THM-667): ‖t·τ‖ ≥ 1/14 for
   t ∈ (V/14, 9V/14) at the discrete snap points — a shifted small-denominator
   condition on a perturbed arithmetic progression of τ's.

kps-S112's formalization pins the continuum side of all five (LRCSmoothBridge:
measure→point desingularization + drift-free observer, sorry-free); the residue named
there — "the Kronecker realization: cluster confinable + φ lifts to real τ" — is this
node. My S207 and kps-S112 proved the DRIFT half of the old obstruction is a
discretization artefact; what is left is genuinely Diophantine, not analytic.

## Three near-misses worth staring at (the numerology of the wall)

- m_P = 14249/252252 ≈ 0.0565 vs 1/14 ≈ 0.0714: the Bonferroni floor is JUST below
  the deterministic-winding threshold (a single good arc of length 1/u_min ≥ 1/14
  would realize by IVT — the floor misses it by 21%).
- ∫W ≈ 0.115–0.140 vs 1/7 ≈ 0.143 (THM-665 ledger): the first moment is JUST below
  the threshold that would hand the drift-embed a 2/7-gap period. First-moment methods
  can therefore NEVER certify the robust gap; the 1/7-vs-2/7 mismatch is structural
  (THM-530's ρ*_{2/7} has zero uniform floor).
- The mid-band edges (V/14, 9V/14): the lower edge is the LRC threshold 1/14 itself
  (ride slack), the upper is the aliasing window 14/5 (2.8). The band width 4/7 is
  the same 4/7 = 1 − 3/7 that appears in scale_separation_phase's (ii). The constant
  14 sits on BOTH sides of the wall — the problem is self-similar at the threshold.

## What closure requires (the honest list)

1. **The mid-band realization** — any one of:
   (a) cancellation: prove any factor-4 of the measured sqrt-cancellation in the
   aliasing corner sum (THM-665 cons. 4 shrinks the window to the drift threshold);
   (b) counting: show the good snap points (≥ (7/6)·∫W·V of them per period-1) cannot
   all be killed by ≤ 5 mid-band danger conditions — a Denjoy–Koksma / three-distance
   count on the perturbed AP of snap points, where the wobble is the enemy;
   (c) the pair-sum axis (THM-668): the witness lives at p/(v_i+v_j); mid-band pairs
   give an unexplored ruler family (HYP-5720's uniform slack σ ≥ 0.034 adversarial is
   the first data);
   (d) the tower (HYP-3901): renormalize the mid-band one scale down — it loses ~6
   runners per level and terminates; nobody has written the composed version.
2. **LEM-014's (H1) bookkeeping** — port the closed floors to the δ-robust threshold
   (their PZ/moment proofs bound P(W ≥ t) natively; this is stated as mechanical).
3. **Lean transcription** of: the floors (decide-shaped), the aliasing bound (THM-665 —
   Fourier/BV, the heaviest analysis item), LEM-014's chain (elementary), the ladder
   (THM-667 — Lemma A is 3 lines given the LRC(≤13) citation hypothesis; Lemma B is
   rational arithmetic).

## Corrections/refinements this session

- THM-665 consequence 2 ("the a-priori route never fires on covering clusters") holds
  in the ALL-RUNNER co-offset convention only. Under THM-527's actual P∪L split, the
  confined-L stratum of the k≥8 covering family is nonempty at every scale and the
  a-priori chain fires on it 20/20 (census, HYP-5691). Addendum filed; no result
  damaged — the theorem's math is untouched, the scope sentence is refined.
- The apex-7 wall (OPEN-Q-108) does not threaten P̃-positivity: Lemma A kills it via
  LRC(≤13) + Lipschitz. The wall remains real for the *danger-measure union-bound*
  bookkeeping — but positivity was the load-bearing need.

## The meta-lesson

Every instrument built this month (floors, aliasing, snap, embed, ladder) CLOSES some
regime and DIES at the same place: where a runner's speed is commensurate with the
ruler — between 1 and 14 fast-laps per slow window. The problem has been squeezed into
its own definition: LRC(14) is now open exactly for runners that circle 1-to-14 times
while the observer's phase wraps once. That band is where "fast" and "slow" lose
meaning — the genuinely multi-scale core. The next instrument must be honest about
scale-entanglement rather than separation: pair-sum rulers (which mix two speeds) and
the renormalization tower (which recurses) are the two candidates already in the repo
that are BUILT from entanglement. I'd bet on those.
