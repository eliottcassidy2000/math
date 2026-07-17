---
id: THM-955
title: THE CLUSTER GAP LEMMA — whenever its explicit numerator is positive, any k ≤ 6 speeds leave a safe subinterval of width ≥ [(1−k/7)L − k/(7m)]/(1 + k + L·Σx) (exact periodic discrepancy + internal-tooth pigeonhole); the a-priori width floor that feeds norm_glue's base certificates analytically and is the continuous face of the formalization picture's item 3 (a-priori liveCount floors)
status: PROVED (proof repaired in-file: the original raw-tooth count did not prove the displayed constants, and unconditional existence on very short windows was false; exact period-one discrepancy and fully internal tooth counts prove the stated conditional bound) + verified on all 291 positive-bound rows in the original 400-row battery and 100 targeted positive c=6 rows (zero violations; cluster_gap_verify_opus_S336); Lean L1 is production-green in LRCClusterGapBrick (kernel-pure sorted-gap pigeonhole with the necessary positivity hypothesis); L2 periodic-comb enumeration/discrepancy and L3 arcSafe assembly remain
source: opus-2026-07-17-S338; proof correction codex-2026-07-17-S62
depends_on: [THM-928/932 (the gluing consumers), LRCLacunaryNest.window_tail_glue (the formal consumer of the base certificate this lemma produces)]
scripts: 04-computation/cluster_gap_verify_opus_S336.py -> 05-knowledge/results/cluster_gap_verify_opus_S336.out
---

# THM-955 — the cluster gap lemma

**Lemma.** Let B be k ≤ 6 distinct positive integer speeds, m = min B,
and [a, b] any window, L = b − a. With λ = 1/14 (teeth ‖wt‖ < λ), assume

```text
A := (1-k/7)L-k/(7m) > 0.                              (1)
```

Then:

> the window contains a closed subinterval [a′, b′] with ‖w t‖ ≥ 1/14 for
> every w ∈ B and every t ∈ [a′, b′], of width
> b′ − a′ ≥ A / (1 + k + L·Σ_{w∈B} w).                (2)

*Proof.* Put `δ=2λ=1/7`.  First fix `w` and scale the window by `w`; its
length is `s=wL`.  Write `s=q+r`, where `q=floor(s)` and `0≤r<1`.  The danger
indicator is one-periodic and has mass `δ` in every full period.  On the
remaining circle arc of length `r`, its intersection with the danger arc has
length at most `min(r,δ)`.  Hence

```text
|[a,b] ∩ {t: ||wt||<1/14}|
 ≤ L/7 + [min(r,1/7)-r/7]/w
 ≤ L/7 + 6/(49w)
 < L/7 + 1/(7w).                                     (3)
```

The maximum `6/49` occurs at `r=1/7`.  Summing (3) and using `w≥m` shows
that the common safe set has measure strictly greater than `A`, so (1) also
guarantees that it is nonempty.

It remains to count its components without losing the sharp denominator.
A danger tooth for `w` can split the window only if both of its endpoints are
strictly internal.  Its centre `n/w` then satisfies

```text
wa+1/14 < n < wb-1/14.
```

This open interval has length `wL-1/7`, so it contains at most
`max(0,wL+6/7) < wL+1` integers.  Teeth meeting an endpoint of the original
window can trim its leftmost or rightmost safe component but cannot create a
new one.  Consequently the number `C` of common safe components satisfies

```text
C ≤ 1 + Σ_w (# fully internal w-teeth)
  < 1+k+L Σ_w w.                                      (4)
```

The common safe set is a finite union of closed intervals.  By (3)--(4), its
widest component has length at least its total measure divided by `C`, and
therefore at least the right-hand side of (2). ∎

The positivity hypothesis in (1) is necessary for an every-window existence
statement: for example, a sufficiently short window contained in the danger
tooth of speed one contains no safe point.  All block-cascade applications
below impose an entry-length inequality that makes `A>0`.

**Why it matters.** (1) It converts `window_tail_glue`'s base-certificate
hypothesis from per-instance decidable DATA into an ANALYTIC guarantee:
any ≤ 6-speed base block, no structure assumed, yields a certified window;
a 7/3-tail above it then glues (norm_glue) — so ≤6-block + gap + lacunary
towers are lonely with NO search. (2) It is the continuous face of the
formalization picture's open item 3: an a-priori floor on surviving
structure inside any window — width floors and liveCount floors are the
same pigeonhole seen from the two sides of the discrete/continuous bridge.
(3) k ≤ 6 is the union-bound horizon (k/7 < 1); beyond 6 the covering
program's machinery is genuinely needed — the lemma marks the exact
boundary where elementary methods end.

**Lean route.** Three layers.  (L1) is now production-green in
`TournamentH7.LRCClusterGapBrick`: `sorted_gap_pigeonhole` proves that
removing `N` open rational intervals of total clipped length at most `B` from
`[a,b]`, under `0<b-a-B`, leaves a closed avoiding interval of width at least
`(b-a-B)/(N+1)`.  The same module kernel-checks a counterexample to the
hypothesis-free form.  The live pieces are (L2) finite tooth enumeration and
the exact period-one discrepancy for each comb, and (L3) assembly into
`arcSafe` with the internal-tooth component count.  The scope-corrected draft
contains those two remaining `sorry`s.  Consumer:
`window_tail_glue`/`norm_glue` immediately after L2--L3.
