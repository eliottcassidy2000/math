---
id: HYP-3796
title: KACZMARZ + CHRISTOFFEL/RKHS + BLASCHKE-DYNAMICS -- three more facets of the LRC group-action frame (S79). (A) The witness search IS a CONSTRUCTIVE FEASIBILITY problem: alternating projections (POCS/Kaczmarz) onto the safe sets S_v={t:||vt||>=r} builds a lonely time (verified: found in 18/40 random starts for the construction), and its conditioning = the resonance = the singular-series crux. (B) The Christoffel/Christoffel-Darboux reproducing kernel is the moment-computable local-density (= the RKHS face of flat-extension S76 / orbit-count S79); it connects to the Beurling-Selberg / Delsarte extremal (HYP-2948/THM-534); the naive observer-lens is subtle (peaks at near-atom, not lonely, times). (C) From arXiv:2604.16750 (generalized Blaschke products, Julia sets, rotation numbers, Arnold tongues): the runner is a degree-1 CIRCLE MAP, loneliness = avoiding all Arnold tongues (mode-locking), and the construction's witness t*=[0;n-1,n] has BOUNDED partial quotients = Diophantine / Herman-rigid = the dynamical face of the deep-well isolation (S77). +6 new loop-functions.
status: CONFIRMED (Kaczmarz witness-search verified constructive; Diophantine CF picture exact; Blaschke/dynamics = a documented merge) with one CAVEAT (the Christoffel observer-lens numeric is direction-subtle, flagged). A MERGE + exploration, NOT a proof: all three are facets of the same resonance/singular-series crux (S79 orbit-count Ramanujan error), each giving a tool (constructive algorithm / moment-extremal / dynamical vocabulary) but none bypassing it. Proof-relevant takeaways: the Kaczmarz witness-finder is a real algorithm (engineering) and its conditioning IS the resonance; the Diophantine rigidity of t*=[0;n-1,n] is why the construction is isolated.
source: mac-mini-2026-06-30-S80
related:
  - HYP-3797   # S79 loop-function dictionary + group action (this adds dynamics + feasibility + RKHS facets)
  - HYP-3792   # S77 deep-well isolation (= the Diophantine/Herman rigidity of t*=[0;n-1,n])
  - HYP-3789   # S76 flat-extension (= the CD-kernel rank / Christoffel)
  - HYP-3793   # kps-S7 Verblunsky/Ramanujan (the CD kernel lives on the Toeplitz moments)
  - HYP-2948   # Beurling-Selberg minorant (= the Christoffel/Fejer extremal)
  - HYP-3782   # lazy-cut (= the exact/LP version of the Kaczmarz feasibility search)
results:
  - 04-computation/kaczmarz_christoffel_blaschke_lrc_macmini_20260701.py
  - 05-knowledge/results/kaczmarz_christoffel_blaschke_lrc_macmini_20260701.out
references:
  - "Carranza Trejo & Moreno Rocha, arXiv:2604.16750 -- Puzzle Pieces, Bi-accessibility, and Connectivity of the Julia Set for Generalized Blaschke Products (rotation numbers, Arnold tongues, circle-map dynamics)"
  - "Kwapien-Mycielski -- the Kaczmarz algorithm and reproducing kernels (stationary sequences / POCS)"
---

# HYP-3796 -- Kaczmarz + Christoffel/RKHS + Blaschke-dynamics: three facets of the group-action frame

Owner's seed: merge arXiv:2604.16750 (generalized Blaschke products, Julia sets) + prior Kaczmarz +
Christoffel/reproducing-kernel work; explore new loop-functions; push for proofs. The paper is about circle-
map DYNAMICS of Blaschke products -- exactly the dynamics of the Blaschke/Verblunsky loop-functions from S79.

## (A) Kaczmarz / alternating projections = the witness search, made constructive
The safe set `S_v = {t : ||v t|| >= r}` (a union of `v` arcs); the lonely set `L = intersect_v S_v`. The
method of alternating projections (POCS, non-convex) cycles projecting `t` onto each `S_v`; a fixed point in
`L` is a lonely time -- a **constructive existence certificate**. Verified: POCS finds a lonely time in
**18/40** random starts for the construction (fewer for the tight AP/GW, whose `L` at `r=1/14` is measure
zero). The **convergence conditioning = the transversality/angle between the safe sets = the resonance**:
maximally-resonant (near-construction) configs are the hardest, matching the deep well (S77). So Kaczmarz is
(i) a real **algorithm** (an engineering deliverable: a POCS LRC witness-finder), and (ii) yet another face
of the singular-series crux (its rate is the resonance). Its exact/LP form is the lazy-cut (HYP-3782); its
finite decision-procedure scope is the MSS ceiling (S78).

## (B) Christoffel / Christoffel-Darboux reproducing kernel
`lambda_N(z) = min{ int|P|^2 dmu : P(z)=1, deg<=N } = 1/K_N(z,z)` is the **moment-computable local density**
of the runner-measure -- the **RKHS face** of the flat-extension (S76) / orbit-count-=-#atoms (S79): the CD
kernel `K_N(z,w) = sum_j phi_j(z) conj(phi_j(w))` has rank `= #atoms`. It is the natural home of the
**Beurling-Selberg / Fejer extremal** (HYP-2948) and the Delsarte moment-LP (THM-534). CAVEAT: the naive
"observer-lens" `max_t K_N(1,1)` peaks where a runner is *closest* to the observer (an atom at `z=1`), not
at the lonely time -- the direction is subtle; a correct loneliness detector needs the extremal polynomial
pinned off the support, i.e. the Beurling-Selberg minorant. (opus-S13's lens needs this care.)

## (C) Blaschke dynamics (arXiv:2604.16750) -- 6 new loop-functions + Diophantine rigidity
The runner `t -> t + v t` is a **degree-1 circle map**; its **rotation number** is the phase; **mode-locking
= rational rotation number = the danger/resonance**; the danger arcs are **Arnold-tongue** slices; the
**Blaschke product** is the disk-dynamics generator whose `|alpha| -> 1` degeneration is the atomic/lonely
limit (the Verblunsky recursion). New dictionary entries (24-29): rotation number, Arnold tongue, Blaschke
product, devil's staircase (`rho` vs parameter = the resonance Cantor set = klein-S67's `13Z` lattice),
Poincare return map (= the three-gap map, S73), Denjoy/Herman conjugacy.

**Diophantine rigidity (the proof-relevant merge).** The construction's witness `t* = n/Phi6 = [0; n-1, n]`
has **bounded, balanced partial quotients** `(n-1, n)` -- it is the *badly-approximable / Diophantine* point
at its denominator scale, hence (Denjoy-Herman) the **rigid** rotation number farthest from mode-locking.
Restructured beaters bind at shallow small-PQ witnesses (`[0,6,3]`, `[0,3,5]`); huge patches walk the
Stern-Brocot ray `[0; n-1, nk]` upward (last PQ grows `28, 42, ...`, `M` increasing). So the deep-well
isolation (S77) IS the **Herman-rigidity of the Diophantine witness `[0; n-1, n]`**: the construction sits at
the one arithmetically-rigid node, and the covering obligation forces the `(n-1)`-multiple onto that lattice.

## Synthesis + honest proof status
Kaczmarz (feasibility), Christoffel (RKHS extremal), and Blaschke (dynamics) are three more coordinates on
the S79 group-action object: the resonance is the Arnold-tongue mode-locking, the extremal certificate is the
Beurling-Selberg/Christoffel minorant, the search is a convergent alternating projection, and the rigidity is
Diophantine. None bypasses the crux (bound the resonant Ramanujan error, S79/THM-501), but each contributes a
tool and the dynamics vocabulary geometrizes *why* the construction is special (a rigid Diophantine rotation
number). NOT a proof; a merge that widens the toolkit and sharpens the rigidity picture.
