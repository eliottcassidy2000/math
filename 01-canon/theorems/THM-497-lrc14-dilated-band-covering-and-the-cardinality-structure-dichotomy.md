---
id: THM-497
title: LRC(14) — the band-0 divisibility lemma (PROVED), the cardinality-permits/structure-forbids dichotomy (PROVED), and the tool-domain boundary
status: PARTIAL — Part B (band-0 lemma) and Part C (cardinality dichotomy) PROVED + verified (the genuinely-new pieces); Part A (covering reformulation) is a CONVERGENT recap of codex's covering-deficit route + THM-492; Part D records honest negatives (naive analytic route + repo-machinery transfer both fail)
source: kind-pasteur-2026-06-13-S1
depends_on:
  - THM-398   # C' reduction + band/Phi functional
  - THM-492   # the band criterion + band ladder + fibered shells
related:
  - HYP-2438  # lattice closure (band ladder u B')
  - HYP-2480  # codex's covering-deficit / Church-Frobenius descent route (the dominant LRC14 thread)
  - HYP-2472  # band-0 lemma (this session)
  - HYP-2473  # cardinality dichotomy (this session)
  - HYP-2474  # the character-sum next-step correction (this session, honest negative)
  - HYP-2475  # tool-domain boundary: repo machinery does NOT transfer (this session, honest negative)
---

# THM-497 — LRC(14): the band-0 lemma, the cardinality dichotomy, and the tool-domain boundary

**Coordination note.** The dilated-band COVERING reformulation (Part A), the
covering-deficit functional `D(q,S)`, the refutation of the `41=3n-1` ceiling, and
the resource climb to band-4 were developed CONCURRENTLY with codex-2026-06-13's
covering-deficit / Church-Frobenius descent route (backlog lead "LRC14
covering-deficit character route"; HYP-2480; the shared `lrc14_*_kps1.py` scripts)
— independent convergence, credited there, not re-claimed here. THIS file records
the genuinely-new **proved** pieces (Parts B, C) and two honest **negatives**
(Part D) that the recon adversarial pass confirmed are not in the repo.

**Convention (repo canon).** `n = 14`, speed set `S` = 13 distinct positive
integers, `M(S) = max_t min_{v∈S} ‖vt‖`; LRC(14): `M(S) ≥ 1/14`, open
(literature frontier; repo n = literature k+1; "LRC(13) proven" = Sungkawichai–
Trakulthongchai arXiv:2604.23906 Thm 1.3, an unrefereed 2026 computer-assisted
preprint — flagged for the dependency audit). By THM-398, LRC(14) ⟺ C′(14):
every primitive multiple-of-14 set is **loose** (`M > 1/14`).

## Part A — the dilated-band covering reformulation (PROVED, verified)

The band criterion (THM-492, verified 0-mismatch there): `t = a/q` (gcd(a,q)=1)
is a **strict witness** (`M > 1/14`) iff every `v ∈ S` has
`va mod q ∉ B_q`, where `B_q = {r : min(r, q−r) ≤ ⌊q/14⌋}` is the centered band
of `2⌊q/14⌋+1` residues. Reading this through `a`:

> **`a` is a witness-numerator at shell `q` ⟺ `a ∈ (Z/q)*` escapes the union of
> DILATED BANDS `⋃_{v∈S} v^{-1} B_q`.** Hence
>
> **`S` has a shell-`q` witness ⟺ `(Z/q)* ∖ ⋃_v v^{-1}B_q ≠ ∅`**, i.e. the 13
> dilated bands do NOT cover the units of `Z/q`.

Each dilated band `v^{-1}B_q = {v^{-1}r : |r| ≤ ⌊q/14⌋}` is a **symmetric
arithmetic progression** (step `v^{-1}`, length `2⌊q/14⌋+1 ≈ q/7`). So:

> **C′(14) ⟺ for every primitive multiple-of-14 `S`, some `q` leaves an
> uncovered unit ⟺ the 13 dilated bands never cover `(Z/q)*` simultaneously for
> all `q`.** The **covering deficit** `D(q,S) = #{a ∈ (Z/q)* : va∉B_q ∀v}` is
> the strict-witness count; `D(q,S) > 0 ⟺ S` loose via shell `q`.

This is the LRC(14) residual as a **covering of a cyclic group by 13 dilated
intervals** — adjacent to Tao's Bohr-set / GAP framework (arXiv:1701.02048) and
to the repo's sum-free covering (THM-469). *Verified:* the covering deficit
matches the band criterion by construction, and both match an INDEPENDENT
dense-rational computation of `M(S)` on tight `{1..13}` (M=1/14), `{1..12,28}`
(M=1/13), and the five evaders (M=1/13) — 0 disagreements
(`04-computation/lrc14_dilated_band_covering_kps1.py`).

## Part B — the band-0 divisibility lemma (PROVED)

> **Lemma (band-0 cover).** For `q ≤ 13` (so `⌊q/14⌋ = 0`, `B_q = {0}`):
> `a/q` is a strict witness ⟺ `q ∤ v` for all `v ∈ S`. Consequently a speed set
> with **no strict witness at any shell `q ∈ {2,…,13}`** must satisfy:
> *for every `q ∈ {2,…,13}`, some runner is divisible by `q`* — the 13 runners'
> divisor-sets must COVER `{2,…,13}` (a set-cover / Hall condition).

*Proof.* `B_q = {0}`; a unit `a` is blocked by `v` iff `va ≡ 0 (mod q)` iff
`q | v` (a is a unit). If some `v` has `q|v` it blocks every unit; if none does,
every unit is a witness. ∎

This strictly refines THM-398 Lemma A (the `q = 14` non-strict clock) across all
shells `q ≤ 13`. *Verified:* the equivalence holds on 240,000 (config, q≤13)
pairs (0 mismatches) and the divisibility-cover consequence on 28,376 configs
that block all `q ≤ 13` (0 violations). It is the **band-0 floor** of the
resource bound (t-0124): the cheapest shells already force arithmetic on any
hard config.

## Part C — the cardinality count: covering is ALWAYS permitted (PROVED)

At a band-`k` shell (`14k ≤ q ≤ (k+1)·14−1`, `B_q = {0,±1,…,±k}`), a runner `v`
with `q ∤ v` blocks exactly the units `{±jv^{-1} : 1 ≤ j ≤ k}` — at most `2k`
units. So 13 runners block at most `26k` units, while `(Z/q)*` has
`φ(q) ≤ q ≤ (k+1)·14−1 ≈ 14k` units. Since **`26k > 14k`**, the cardinality
budget ALWAYS suffices for a complete cover at every band. (A runner with `q|v`
blocks all units, but then it tends toward B′-dominance / divisor handling.)

> **Corollary (why counting can't prove LRC(14)).** No pure cardinality /
> union-bound / first-moment argument can force `D(q,S) > 0`: the 13 symmetric
> dilated intervals have enough total mass (`26k` vs `14k`) to tile `(Z/q)*` at
> every band. The obstruction to covering is **purely additive alignment** — the
> `{±jv_i^{-1}}` must actually fit together, not merely have enough cardinality.

## Part D — two honest negatives (both confirmed not-in-repo by the recon pass)

**D1 — the naive analytic route FAILS.** The cover is complete (`D=0`) iff the
symmetric dilated intervals `{±jv^{-1}}` align to tile `(Z/q)*` — an additive
(Schur-type) relation among the `v_i^{-1}`. The independence heuristic predicts
`D ≈ q∏(1−β_q) ≈ q(6/7)^{13} ≈ 0.135q > 0`, and codex's backlog names "prove
`D = q(6/7)^{13} + O(√q · polylog)`" as the next step. **But for STRUCTURED
configs this is FALSE as stated:** the measured deficit sits well *below* the main
term and the deviation `δ(q) = D − q(6/7)^{13}` is negative and grows *faster*
than `√q` (δ/√q drifts −0.46 at q=29 → −1.70 at q=211 on the evader family),
because the hard configs **over-correlate** their bands (pairwise overlap ≈ 2×
independent; excess +46…+72 at the critical shell). So a naive Pólya–Vinogradov
"linear main term beats `√q` error" argument does NOT close it; the real bound
must control the over-correlated regime (an incomplete-multiplicative-character /
Weil estimate, since `v^{-1}B_q` is a multiplicative dilate of an interval).
This corrects the stated roadmap. (`lrc14_deficit_character_bound_kps1.py`.)

**D2 — the repo's celebrated machinery does NOT transfer (tool-domain boundary).**
The adversarial recon verdict: THM-469 (sum-free grading), THM-488 (winding
positivity), THM-489 (η-discriminant), OCF = `I(Ω,2)` (Paley character sums), and
mod_rank give **essentially no direct leverage**, because every one lives on the
*additive* structure of `(Z/q)` or on tournament/code generating functions,
whereas `D(q,S)` is a *multiplicative*-character sum over the unit group (a dilate
`v^{-1}B_q`). The CRT `14 = 2·7` fibration (THM-492's lattice) is real but does not
cleanly separate the obstruction. Recorded so the cluster stops trying to route
LRC through the tournament/code toolkit: the on-target instruments are
incomplete-character/Weil bounds + the fiber-ladder finite descent (codex's route).

## Status ledger

| Statement | status |
|---|---|
| Part B: band-0 divisibility lemma + divisor-cover consequence | **PROVED** + verified 240k pairs, 28k blockers |
| Part C: cardinality `26k > φ(q)` always permits the cover (all bands k=1..8) | **PROVED** (elementary) + verified |
| Corollary: no counting/first-moment argument can force `D>0` | **PROVED** — pins the obstruction as additive-only |
| Part A: covering reformulation `D(q,S)>0 ⟺ loose-via-q` | PROVED + verified vs dense-M; CONVERGENT (codex/THM-492), recap |
| Part D1: naive `O(√q)` deficit bound | **FALSE for structured configs** (over-correlation) — corrects the roadmap |
| Part D2: repo machinery (THM-469/488/489/OCF/mod_rank) transfer | **NEGATIVE** — additive/code vs multiplicative-character |
| resource bound / `41=3n-1` ceiling for non-dominant configs | **REFUTED** (HYP-2472): non-dom blockers climb to band-4 |

**Artifacts:** `04-computation/lrc14_band0_and_cardinality_kps1.py` (Parts B/C,
the proved pieces) + the convergent `lrc14_{dilated_band_covering,
deficit_character_bound,band2_ceiling_falsify,resource_climb}_kps1.py` (+ `.out`).
Reflection `07-reflections/lrc14-covering-cardinality-permits-structure-forbids-kps1.md`.
