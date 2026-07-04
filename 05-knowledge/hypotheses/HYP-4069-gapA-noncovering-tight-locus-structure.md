---
id: HYP-4069
title: GAP-A structure — the non-covering tight locus is {AP, GW}, reduced to a finite check via forced-residues + lift-rigidity + coverer-bound; GW is the UNIQUE non-AP single-swap tight family (only q=12 admits a second coverer, 24)
status: PARTIAL (forced-residues RIGOROUS; single-swap uniqueness VERIFIED exact; the two contingent pieces — odd/even lift-rigidity and the coverer magnitude bound — have strong evidence + the single-swap proof, not a full proof). GAP-A = tight-locus finiteness for non-covering families = OPEN in the literature (Perarnau–Serra).
source: mac-mini-2026-07-03-S35
related:
  - THM-612   # the tight-locus tower; S34 reframing split the open core into GAP-A (this) + GAP-B (covering-min)
  - THM-523   # q-witness lemma — gives "covers {2..13}" and the deletion-hiding mechanism
  - HYP-2913  # ±units-covering (gives the units ⊆ residues) + three-gap
  - HYP-2561  # tight-locus finiteness (inf L=1/1260, locus={AP,GW})
results:
  - 04-computation/gapA_structure_macmini_20260703.py
  - 05-knowledge/results/gapA_structure_macmini_20260703.out
  - 04-computation/tight_locus_enumerate_macmini_20260703.py
external: LRC(14); Goddyn–Wong tight instances Integers 6 (2006) #A38; tight-locus finiteness open (Perarnau–Serra arXiv:2409.20160).
---

# HYP-4069 — GAP-A: the non-covering tight locus is {AP, GW}

**Object.** GAP-A (from THM-612 S34): classify the *primitive tight* families that **miss 14** (non-covering).
By S34 these are automatically at `q*=14` (tight at `t=1/14`), cover `{2,…,13}`, phases on the 14th-root grid.

## 1. Forced residues (RIGOROUS): odds `{1,3,5,7,9,11,13}` are all present
- **units `{1,3,5,9,11,13} ⊆ R'`**: `M=1/14` at every `t=u/14` (`u` a unit) ⟹ `u·R'` covers `±1` ⟹ `R'`
  covers `±u^{-1}`; over all units, `R' ⊇` the ±units (HYP-2913).
- **`7 ∈ R'`**: covering `q=7` needs `7|v`; then `v mod 14 ∈ {0,7}`, and miss-14 rules out `0`, so `=7`.
- For odd `q∈{3,5,9,11,13}` (coprime to 14) the covering multiple can carry ANY residue, so no further
  residue is forced; `q=2` forces only "some even residue". **Net rigorous floor: all 7 odd residues + ≥1 even.**
Verified: AP residues `{1..13}`, GW residues `{1..11,13}` — both contain all odds, miss 14.

## 2. GW is the UNIQUE non-AP single-swap tight family (VERIFIED exact)
Replacing AP runner `k` by a multiple `j·k` (`j≥2`, which blocks the spot `t=1/k` that deleting `k` opens):
**only `k=12`** yields a tight family, and only `X=24`. Every other `k∈{2,…,11,13}` and every multiple
gives `M>1/14` (loose). So the AP/GW dichotomy is exactly **"how to cover `q=12`": 12 (→AP) or 24 (→GW)**.
`q=12` is special because `24=2·12` vacates the **non-unit** residue 12 and doubles residue 10, preserving
unit-tightness; other swaps either vacate structure or open a new spot.

## 3. The deletion-hiding + monotone-loosening mechanism (the speed-bound engine)
For `k≥7` (`2k>13`), `{1,…,13}∖{k}` hides at `t=1/k` (min `1/k>1/14`, q-witness). A tight replacement
must block `t=1/k` (a multiple of `k`) AND open no new spot. For the `q=12` axis, `{1..11,13,X}` (`12|X`):
`M=1/14` iff `X∈{12,24}`; for `X≥36`, `M(X)` **increases monotonically** (`3/41, 4/53, 5/65, 7/89, …`
`→ ~0.081`), the family hiding at a spot where **`X` itself binds** (`||Xt||=M`). So large coverers are
loose ⟹ **the coverer is magnitude-bounded** (here `≤24`). Likewise every residue-preserving **lift loosens**
(odd `r→r+14`: all `M>1/14`, verified; even lifts: S32) ⟹ **runners sit at minimal lift**.

## 4. The reduction (GAP-A ⟹ finite check ⟹ {AP,GW})
Granting §3's two mechanisms (minimal lifts + bounded coverers), a non-covering tight family is
**7 fixed odd runners `{1,3,5,7,9,11,13}` + 6 even coverers of `{2,4,6,8,10,12}` with bounded speed** — a
FINITE set of candidates, exhaustively `{2,4,6,8,10,12}` (AP) or `{2,4,6,8,10,24}` (GW). Exact-`M`
enumeration (all single/double swaps+lifts + `120k` random to speed 60/80) realizes **exactly `{AP,GW}`,
max speed 24, `g≤2`**.

## What is rigorous vs. open
- **RIGOROUS**: §1 (odds forced); §2 (single-swap uniqueness of GW, exact); "covers {2..13}", "miss-14 ⟹ q*=14".
- **STRONG EVIDENCE, not proved**: the two magnitude bounds of §3 (lift-rigidity: minimal lifts; coverer-bound:
  coverers `≤` explicit `B`) — each holds on the single-swap axis (proved) and over all searches, but the
  general (multi-swap, arbitrary lift) bound is the residual = the finiteness itself. Closing §3 in general
  closes GAP-A. This is the LRC(14) non-covering extremal-uniqueness — open in the literature.

NET: GAP-A is reduced to **two magnitude bounds** (minimal lifts + bounded coverers) on the 14th-root grid,
both rigorous on the single-swap axis and strongly evidenced generally; granting them, GAP-A = {AP, GW}.
