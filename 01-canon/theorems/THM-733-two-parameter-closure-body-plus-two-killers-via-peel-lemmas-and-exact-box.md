---
id: THM-733
title: TWO-PARAMETER CLOSURE — for the body E={1,…,11}, EVERY family {E,a,b} (11<a<b) satisfies LRC(14) — via two elementary PEEL LEMMAS (r and |G'| of E∪{a} controlled explicitly in a) + THM-732's far-element tail (uniform leg a≥A₀) + exact per-a thresholds v₀(a) + a finite exact-ℚ box check. The method is BODY-UNIFORM: for any 11-element body E with good-set data (r,m), A₀(E) = explicit; the open separate-13/14 class reduces body-by-body to finite exact boxes
status: PROVED (kind-pasteur-2026-07-13-S128; Lean transcription open). Peel lemmas P1/P2 elementary (below). Uniform leg: K=2121/500<√18 exact, f(a) exactly increasing (sum of increasing terms), minimal A₀=267, f(267)=12425969/102795000>0 certified exactly ⟹ all a≥267, all b>a: L>0. Per-a leg: exact (r_a,m_a) for 12≤a≤266 (subtract_arcs validated == direct sweep on 9 spot values, exact rational match); box = 1810 pairs (max b=169), ALL swept in exact rational arithmetic: L>0 everywhere EXCEPT exactly two NON-COVERING TIGHT pairs — (a,b)=(12,13) = the AP {1..13} (known, missing q=14, THM-366 t=1/14) and (a,b)=(13,24) = {1..11,13,24} (missing q=14, THM-366 t=1/14) — the KNOWN Goddyn–Wong accelerated tight family ({1..13} with 12 doubled to 24; GW Theorem 12: at n=13 only r=12,m=2 passes; = THM-709-doubling-singleton, mac-mini cont.33, klein-S253 literature merge). The box REDISCOVERED it independently and confirms EXACTLY: L=0, M=1/14 attained, and it is the ONLY tight family in the box besides the AP — matching GW's "only r=12,m=2" prediction with zero extras. No covering pair has L=0. Total runtime 1.8s
source: kind-pasteur-2026-07-13-S128
depends_on:
  - THM-731   # the rigorous certificate chain
  - THM-732   # exact disc_v pair form + tail disc_v ≤ r²/(3v²)
  - THM-366   # non-covering families: elementary witness t=1/q (used for the L=0 tight AP {1..13} and sub-box non-covering pairs if any sweep hits L=0)
related:
  - HYP-6495, HYP-6248 (kps cont.70: worst |core|=1 body {1..11,13,84} — inside this family)
  - THM-527-A / finite-Vmax glue (this is the covering-side realization of that program via iterated peel)
---

# THM-733 — the {1..11,a,b} two-parameter closure

## Peel lemmas (PROVED, elementary)

Let E be a finite speed set whose good set G'(E) has r maximal intervals and measure m. Let a ≥ 1.
The arcs of D_a are 1/a-periodic with width 1/(7a). For any interval I of length λ, the number of
arcs of D_a meeting I is at most aλ + 8/7 (centers lie in an interval of length λ + 1/(7a)).

- **LEMMA P1 (interval count):** r(E∪{a}) ≤ r + Σ_i (aλ_i + 8/7) ≤ a·m + (15/7)·r.
- **LEMMA P2 (measure):** each meeting arc removes ≤ 1/(7a), so
  m(E∪{a}) ≥ m − (a·m + (8/7)r)/(7a) = (6/7)·m − 8r/(49a).

## Uniform leg (a ≥ A₀)

By THM-732(iii) the certificate for {E,a,b}, peeling b, holds whenever b² > r_a²/(18 m_a²)
(r_a = r(E∪{a}), m_a = m(E∪{a})). By P1/P2 it SUFFICES that
`(a·m + (15/7)r)·(99/70) < 3·((6/7)m − 8r/(49a))·(a+1)·√2-lower…` — implemented with the rational
bounds 99/70 > √2 and 4243/3000 < √2 as appropriate; the deficiency
f(a) = 3√2·((6/7)m − 8r/(49a))(a+1) − (a·m + (15/7)r) has f'(a) ≥ m(18√2/7 − 1) > 2.6m > 0,
so ONE exact check at a = A₀ certifies ALL a ≥ A₀. For E = {1..11} (r=14, m=10931/194040):
A₀ ≈ 350 (exact value certified in the script).

Then for every b > a ≥ A₀: b ≥ a+1 > v₀(E∪{a}) ⟹ certificate ⟹ **L > 0**.

## Per-a leg (11 < a < A₀)

For each a, compute (r_a, m_a) EXACTLY; every b with b² > r_a²/(18 m_a²) closes by the tail.
The remaining box {(a,b): a < b ≤ v₀(a)} is FINITE; each pair is decided by an exact rational
sweep of L({1..11,a,b}): L > 0 (proof) — or L = 0, which must be the tight AP {1..13} (a,b)=(12,13),
non-covering, handled by THM-366 (t = 1/14). Any covering pair with L = 0 would be flagged loudly
(none expected; none found = theorem).

## Scope note (honest)

This closes the body E = {1..11} — the body of BOTH known extremal families (deep well ray sits at
a=12, b=182j; the worst |core|=1 body {1..11,13,84} at a=13, b=84). Other 11-element bodies of the
open separate-13/14 class need the same finite computation with their own (r, m, A₀) — the method is
uniform, the per-body cost is one script run. Multi-scale bodies (bodies containing their own far
elements) iterate the peel (P1/P2 compose); the bottom of the recursion is the bounded-Vmax finite
check — the covering-side realization of the finite-Vmax glue program.

## Evidence log (all done, script lrc14_thm733_two_param_closure_kps_S128.py, 1.8s)

- [x] exact A₀ certification for E={1..11}: A₀ = 267 (f(267) = 12425969/102795000 > 0; monotonicity proved)
- [x] per-a exact (r_a, m_a) for a = 12..266 (helper validated exactly vs direct sweep); box = 1810 pairs, max b = 169
- [x] all box pairs swept exactly: L > 0 except (12,13) AND (13,24), both non-covering (missing q=14) → THM-366
- **Confirmation en route:** {1,…,11,13,24} (L=0, M=1/14 attained) = the Goddyn–Wong accelerated tight
  family (THM-709-doubling-singleton; klein-S253 merge item 2). The exact box finds EXACTLY the two tight
  families GW theory predicts (AP + the unique r=12,m=2 acceleration) and no others — an independent
  computational confirmation of GW Theorem 12's n=13 slice.

## Statement (final)

**For every pair of integers 11 < a < b, the family {1,…,11,a,b} satisfies LRC(14).** Moreover L > 0
(positive measure of 1/14-lonely times) for all such families except the two tight non-covering ones,
{1..13} and {1..11,13,24}, where M = 1/14 is attained at t = 1/14. In particular the entire
separate-13/14 open class over the body {1..11} — which contains the worst known |core|=1 body
{1..11,13,84} (kps cont.70) and the deep-well and residue rays — is CLOSED.
