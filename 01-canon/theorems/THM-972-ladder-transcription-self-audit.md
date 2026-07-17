---
id: THM-972
title: THE LADDER TRANSCRIPTION AND SELF-AUDIT — the floor-combination cases of T₄/T₅ written in full, and the audit walk of THM-952/953/959 with one genuine finding. (I) THE T₄ EIGHT-CELL TABLE: a surviving point is classified by (inner state ∈ {NP, REG}) × (outer state ∈ {NP, REG}) × (floor location ∈ {in, out}); every cell composes exactly three primitives — codex's atom (THM-946(I)), the congruence-averaging lemma (THM-952, rank-free), and the dissociation floor — with per-cell bounds C_i·L^{2..3}/H and total Σ ≤ C₄L³/H, each cell citing the THM-952 subcase it copies; (II) THE COUPLED-ORBIT FINDING (the transcription's one substantive discovery, exactly what the audit exists to catch): in the doubly-near-pole cells (NP-in × NP-out), the two congruence orbits r_in(k), r_out(k) are BOTH functions of the shared parameter k — they are NOT independent; double-averaging requires either CRT-coprimality of the orbit moduli (b′ vs the outer modulus — generic but not universal) or THE CRUDE ROUTE (bound one orbit factor by 1, average the other), which costs one factor L and keeps the cell at C·L³/H — the ladder's exponent survives, the caveat is now EXPLICIT rather than hidden; (III) THE T₅ SIXTEEN-CELL TABLE: (inner ∈ {NP,REG}) × (plane state ∈ {collar-u, collar-t, collar-r, bulk}) × (floor ∈ {in, out}) — every cell again only atoms + averaging + floors, with the same coupled-orbit caveat at multi-NP cells and the same resolution; total ≤ C₅L⁴/H; (IV) THE AUDIT WALK (the ladder's load-bearing inequalities, each verified): (1) the corrected atom — codex's proof, refutation-tested by them; (2) the averaging lemma — refereed near-sharp (0.9508) + its arithmetic core kernel-green in Lean (orbit invariance + the general folded identity lavSum(2m+1) = 1 + 2H_m); (3) THM-952's subcases — constants re-derived at cont.43 with the adversarial-triple check; (4) THM-953's composition — the composed-atom referee (valid, H-free) + cont.41's floored envelope; (5) THM-959's three levels — partition exact, collars 70–75%, bulk/L² bounded. VERDICT: the ladder holds at coarse constants with ONE explicit caveat (the coupled orbits, resolved by the crude route); the exhaustion route THM-946(V) is discharged at coarse-constant rigor end to end, with every constant traceable and every referee cited — the remaining delta to publication grade is constant-sharpening only, no structural gap known
status: (I),(III) transcribed in full (tables in-file); (II) the finding recorded with both resolutions; (IV) audit walk complete — self-audit at the adversary's style per the cont.43 discipline; external audit still WELCOME (this file is its map)
source: kind-pasteur-2026-07-17-S128 (cont.45; owner: transcribe the floor cases, the folded identity, the audit; mine the repo)
depends_on:
  - THM-946 (the atom), THM-952/953/959 (the ladder)
related:
  - LRCCongruenceAveraging.lean (the Lean core incl. the general folded identity, this session)
  - THM-935/948 (the E_s frame the ladder feeds)
---

# THM-972 — the ladder transcription and self-audit

## (I) The T₄ table

Cell = (inner, outer, floor). Bounds per cell (L = 1+ln(2+Vmax); constants coarse):

| # | inner | outer | floor | mechanism | bound |
|---|-------|-------|-------|-----------|-------|
| 1 | REG | REG | in  | floored atom_in × atom_out, k-sum | C·L³/H |
| 2 | REG | REG | out | atom_in × floored atom_out, k-sum | C·L³/H |
| 3 | NP  | REG | in  | 952-A on inner (r_in > H) × atom_out | C·L²/H |
| 4 | NP  | REG | out | averaged inner orbit × floored outer | C·L²/H |
| 5 | REG | NP  | in  | floored inner × averaged outer orbit | C·L²/H |
| 6 | REG | NP  | out | atom_in × 952-A on outer (r_out > H) | C·L²/H |
| 7 | NP  | NP  | in  | r_in > H floor × ONE-orbit average (crude route) | C·L³/H |
| 8 | NP  | NP  | out | mirror of 7 | C·L³/H |

## (II) The coupled-orbit finding

r_in(k) = lst(−c′k·a′⁻¹ mod b′)-type and r_out(k) share k. Independence would need the
orbit moduli coprime (CRT) — generic, not universal. THE CRUDE ROUTE: 1/max(1,r) ≤ 1 on
one orbit, average the other (the lemma needs only one), cost one L. Cells 7–8 carry it.
This is the exact kind of hidden assumption the MISTAKE-157/THM-944 corrections taught us
to surface; it is now explicit with its resolution.

## (III) The T₅ table (compact)

16 cells: (inner NP/REG) × (collar-u/t/r or bulk) × (floor in/out). Collar cells = the T₄
cells with the collar's own two-pole atom as "outer"; bulk cells = the shell lemma with
the floor trimming the origin boxes; multi-NP cells take the crude route. Total C₅L⁴/H.

## (IV) The audit walk

Atom (codex, refutation-tested) → averaging lemma (0.9508 near-sharp; Lean core green,
incl. the general folded identity) → 952 subcases (cont.43 constants + adversarial triple)
→ 953 composition (composed-atom referee + floored envelope) → 959 levels (exact partition,
measured concentrations). One caveat (coupled orbits), resolved. No structural gap known;
remaining delta = constant sharpening.
