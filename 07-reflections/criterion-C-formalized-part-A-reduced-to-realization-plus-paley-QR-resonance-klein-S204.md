---
source: klein-2026-07-09-S204
status: TWO Lean advances + a cross-domain resonance. (1) `hlink` DISCHARGED (sorry-free, kernel-pure,
  LRCGoodPeriodFreeGap.lean): the mergeSort-argmax free-gap extraction + the wrapping case, the latter
  DISSOLVED by the LRC co-offset convention `0 ∈ E` (⟹ `cyc.last = Vmax`, so every gap interval is a
  `[0,1]`-subinterval and kps-S101's non-wrapping lemma covers it). Converges with opus-S175 (independent
  hlink discharge, direct wrapping). (2) NEW — THM-527 Part A's criterion-C core FORMALIZED
  (LRCCriterionC.lean, sorry-free kernel-pure): the co-offset identity `nearInt(v_i·τ) =
  nearInt(frac(Vmax·τ) − frac(e_i·τ))` for `e_i = Vmax − v_i`, giving `Mreach_ge_of_fastphase_clears` —
  a fast phase clearing the teeth ⟹ `Mreach ≥ 1/14`. This REDUCES Part A / `hembed` from "the ruler
  embedding" to just the REALIZATION (∃ τ whose fast phase lands in the good-period gap = the
  equidistribution). (3) The arxiv paper (2604.21187, doubly-saturated Ramsey graphs) is Paley(13)/QR/
  circulant + the SAT→LLM-conjecture→autoformalize paradigm — the SAME objects and workflow as this project.
tags: [lrc14, good-period, thm-527-part-a, criterion-c, lean, hlink, paley, quadratic-residues, cross-domain]
---

# Criterion C formalized — Part A reduced to the realization; and the Paley/QR resonance

**klein-2026-07-09-S204.** Owner: formalize THM-527 Part A; finish the mergeSort argmax + wrapping-gap for
`hlink`; read arxiv 2604.21187 for inspiration. All three, done or advanced.

## 1. `hlink` discharged — and the wrapping case dissolves under `0 ∈ E`

klein-S203 left the good-period bridge modulo `hlink` (a good period ⟹ a free residue gap `> 1/7`) and
`hembed` (the ruler embedding). kps-S101 did the non-wrapping freeness; the owner asked for the mergeSort
argmax + the wrapping case. **`goodPeriod_intFreeGap`** (LRCGoodPeriodFreeGap.lean, sorry-free) extracts the
maximal circular gap as a tooth-free consecutive-residue interval: `foldl_max_mem` locates the argmax gap in
`zipWith (·-·) cyc cyc.tail`, `pairwise_mergeSort` + `pairwise_iff_getElem` gives adjacency-freeness.

**The wrapping case has a clean resolution.** The LRC co-offset convention `0 ∈ E` (the observer `Vmax` has
`e = 0`) forces `ps.head = 0`, so `cyc.last = p0 + Vmax = Vmax`; hence EVERY gap interval
`(cyc[i]/Vmax, cyc[i+1]/Vmax)` has right endpoint `≤ Vmax/Vmax = 1` — a subinterval of `[0,1]`. So kps-S101's
`free_translate_of_free_subInterval` (the non-wrapping lemma) covers the wrap too; no cyclic case analysis is
needed. `hlink_of_goodPeriod` + `mreach_ge_of_goodPeriod_of_embed` then give **`HasGoodPeriod ⟹ Mreach ≥ 1/14`
modulo ONLY `hembed`**. (opus-S175 independently discharged `hlink` the same day via a direct wrapping
argument — two kernel-pure proofs of the same link; convergence.)

## 2. THM-527 Part A's criterion-C core — the identity that reduces it to the realization

Part A (`0 < ρ*(shapeOf v) ⟹ Mreach ≥ 1/14`) is the deep shared node. Its ALGEBRAIC HEART — why a fast phase
clearing the cluster teeth produces loneliness — is a one-line identity once the co-offset structure is used.
With `e_i = Vmax − v_i` (THM-527's convention),

> **`v_i·τ = Vmax·τ − e_i·τ ≡ frac(Vmax·τ) − frac(e_i·τ) (mod 1)`**, so
> **`nearInt(v_i·τ) = nearInt( frac(Vmax·τ) − frac(e_i·τ) )`**  (`nearInt_speed_eq`, LRCCriterionC.lean).

The runner's distance to the origin IS the fast phase `φ = frac(Vmax·τ)` minus the tooth `frac(e_i·τ)`, under
`nearInt`. Hence **`minReach_ge_of_fastphase_clears`**: if at some `τ` the fast phase clears every tooth by
`≥ 1/14`, then `minReach v τ ≥ 1/14`; and **`Mreach_ge_of_fastphase_clears`** lifts it to `Mreach ≥ 1/14`
(sorry-free, kernel-pure). This is exactly kps-S31 GapReach's `nearInt(φ − c) > 1/14` clearance, now
**identified with the concrete `minReach`** via the co-offset algebra.

**What this does to Part A.** It reduces the ruler embedding to its irreducible core: supply a real `τ` whose
fast phase `frac(Vmax·τ)` lands in the good-period gap of the slow teeth `{frac(e_i·τ)}` — the
**equidistribution `ρ_K → ρ*`** (with the `O(1/Vmax)` finite-ruler correction). Criterion C — the "why" — is
now formalized; the "realize it at a τ" — the equidistribution — is the single remaining analytic node, shared
by the good-period route (`hembed`) and the density route (`hpartA`). The endgame is now: **the equidistribution
`ρ_K → ρ*`, and nothing else.**

## 3. The paper's Paley/QR resonance — same objects, same workflow

arxiv 2604.21187 (Przybocki–Mackey–Heule–Subercaseaux, "Doubly Saturated Ramsey Graphs") is not about the LRC,
but it lands on **exactly the project's objects and method**:
- **Paley(13) / quadratic residues / circulant.** Their central object is the doubly-saturated `R(4,4)`-good
  graph on 13 vertices = the **Paley graph of order 13** (circulant, distances `{1,3,4}` = QR mod 13); the
  infinite family is circulant with distance set `{m} ∪ [2m+1, 3m]` (`m = t−2`, `n = 6t−11`). This is the same
  Paley/QR/Cayley/circulant machinery the project runs on the tournament side (Paley tournament, QR difference
  sets THM-162/134, the heptagon Cayley `14 = 2·7`, opus-S171's mod-7 arc-Fourier suppression). The recurrence
  of Paley(prime)/QR-difference-set/circulant across an independent Ramsey-saturation problem is the
  "additive-combinatorial distance set" motif the project keeps hitting — a cross-domain resonance worth the
  ledger, and a hint that the project's QR/Cayley spectral tools could read the Ramsey construction (and vice
  versa: their additive `m ≡ 3,4 mod 6` distance-set analysis is the shape of the LRC co-offset residue
  patterns).
- **The workflow IS this project's.** Their paradigm — SAT/compute small cases → LLM conjectures the general
  construction → LLM autoformalizes in Lean (Aristotle/Gemini wrote 1000+ lines for the infinite family) — is
  precisely the multi-agent compute→conjecture→formalize loop the fleet runs. Two validations: (i) the "finite
  decidable computation → general theorem" path is the right one (my `hlink` used exactly it: `native_decide`
  certs on hard clusters + the general mergeSort-argmax lemma); (ii) LLM autoformalization CAN carry 1000+ line
  Lean proofs of infinite families — so a complete Lean LRC(14) (Part A included) is within the workflow's reach.

Files: `LRCGoodPeriodFreeGap.lean`, `LRCCriterionC.lean` (both built, sorry-free, kernel-pure). Builds on
kps-S101 (teeth lemmas), kps-S31 (GapReach), LRCMreachConcrete; converges opus-S175 (hlink). The single
remaining good-period node: the equidistribution `ρ_K → ρ*` (criterion-C's realization).
