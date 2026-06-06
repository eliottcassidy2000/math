---
id: HYP-2341
title: The LRC frontier (unramified n) is the ±-transversal / quasi-random core; C′(n) [2n-1 prime] modulo the transversals
status: OPEN (the residual); the non-transversal part is PROVED (THM-420)
source: claudebox-2026-06-06-S642
related:
  - THM-420  # non-transversal dodge (proved): non-transversal ⟹ M≥2/(2n-1)
  - HYP-2336 # Rado/random tournament (the transversals = the random-like configs)
  - HYP-2321 # Paley = QR transversal; HYP-2280 LRC(19); HYP-2197 constructive route
---

# HYP-2341 — the transversal core is the LRC frontier

THM-420 proves: for `2n-1` prime, every multiple-of-n config that is **not** a ±-transversal mod
`(2n-1)` is loose (`M ≥ 2/(2n-1)`), by a one-line shell dodge — and that is ~100% of configs. The
entire residual is the **±-transversals**: residues mod `2n-1` hitting every ±-pair exactly once.

## The reduction (a major refinement of the constructive route)

> **C′(n) for `2n-1` prime ⟺ every ±-transversal multiple-of-n config is loose.**
> (LRC(n) for the unramified family `n = 15,19,21,22,24,…`, all open `>13`, reduces to this rare core.)

The transversals are the **maximally-spread / spectrally-flat / quasi-random** residue patterns —
the QR set (Paley, HYP-2321) is one; the random tournament's residues are transversal-like (HYP-2336).
So the elementary dodge clears every config *except the quasi-random ones*, and the LRC frontier is
exactly the quasi-random core — the same locus as the character-ratio/Gauss-sum spectrum (S638) and
the Rado-tournament limit (S641). **The hard part of LRC, for unramified `n`, is the random-like
configurations** — poetic and precise.

## Status of the residual (verified loose, not proved uniformly)

Sampled/generated ±-transversal multiple-of-n configs (n=15,19) are **all loose**, via:
- lower 1-clocks `t=1/m`, `m ≤ n-1` (freed when the config has no multiple of `m`): e.g. AP-with-top-
  bumped `{1,…,n-2,n}` is a transversal with `M = 1/(n-1)` (S630);
- larger margins: sampled residuals have `M = 1/4, 1/6, …` ≫ `1/n` (`lrc_transversal_M_s642.out`);
- the measure dodge B′ (THM-398) for the all-short cases.
No counterexample found.

## To finish C′(n) [2n-1 prime] (→ LRC for the unramified family)
1. Prove every ±-transversal multiple-of-n config has a **free lower clock or a long safe component**
   (B′). Candidate: a transversal mod `2n-1` is CRT-independent of residues mod `m≤n-1`, but a
   multiple-of-n transversal *blocking all 1-clocks* would need a multiple of every `m≤n-1` — show
   this forces either a collision (contradicting transversality) or a dominant runner (B′).
2. The residue-profile DP (HYP-2256) restricted to the transversal core — now a finite, *small* check
   (transversals are rare), feasible per `n` (certify n=15, then n=19).
3. Generalize the dodge to `2n-1 = p^k` (ramified, n=14): the ±-pair structure on `(ℤ/p^k)^*` and the
   inner shells (the 3-adic tower) — the non-transversal dodge with the ramified correction.
