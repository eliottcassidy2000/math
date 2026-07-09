---
source: opus-2026-07-08-S164
status: the LRC(14) capstone lemma (good period j*=O(k), THM-527-A) characterized structurally.
  (1) The AP is NOT the j*-maximizer (random/near-AP clusters exceed it at k=8,9,10). (2) j* is
  GOVERNED by the longest-AP L (dilation-invariant axis): j* is small (<=~7) for small L, and large
  (~k, up to ceil(7(k-1)/6)) only for near-APs. (3) THE DICHOTOMY (490k configs, 0 non-AP
  exceptions): the ONLY hard-region clusters with NO good period are the exact full complete-residue
  APs {0..k-1} at Vmax=k, k prime (11,13) = exactly the tight M=1/k LRC cases, cited via LRC<=13.
  So the good-period lemma is [j* = O(k)] OR [tight AP], no genuine counterexample.
tags:
  - lrc14
  - good-period
  - THM-527
  - longest-AP
  - dilation-invariance
  - tight-locus
  - capstone
---

# The good-period dichotomy: the exception is exactly the tight AP

**opus-2026-07-08-S164.** Owner: work the remaining math creatively.  mac-mini-S59 reduced ALL of
LRC(14) to ONE lemma — the good period `j* = O(k)` (THM-527-A): for a `k`-cluster
`E = {0=e_0<..<e_{k-1}}` and ruler `Vmax` (hard region `spread >= 6Vmax/7`),
`j*(E,Vmax) = min{ j>=1 : maxgap{frac(e_i j/Vmax)} > 1/7 }` should be `O(k)`.  The AP case is proved
(Dirichlet, `j* <= ceil(7(k-1)/6)`) and there are `0` counterexamples / 90k.  This note characterizes
the capstone structurally, on the DILATION-INVARIANT longest-AP axis (`j*(cE,cVmax)=j*(E,Vmax)` — the
same gauge as the density floor, opus-S155–S163).

## (1) The AP is NOT the j*-maximizer

Natural guess: the AP (most rigid) needs the largest `j*`.  FALSE.  At `k=8,9,10` random/near-AP
clusters reach `j* = 7` while the AP reaches only `2–3`.  The AP is actually EASY (it dilates to a
clustered set at small `j`).  The `j*`-maximizers are NEAR-APs / dilated APs (e.g. `2·{0..10}` at
`k=11`, `j*=11`).  So "AP maximizes j*" is the wrong reduction.

## (2) j* is governed by the longest AP `L`

Stratifying max `j*` by the longest AP `L(E)` (dilation-invariant):

| k | small L (`<=5`) | large L (near-AP) | bound `ceil(7(k-1)/6)` |
|---|-----------------|-------------------|------------------------|
| 8 | `<= 2` | `7` (at L=7) | 9 |
| 11 | `<= 5` | `11` (L>=6) | 12 |
| 13 | `<= 7` | `13` (L>=6) | 14 |

`j*` is SMALL (`<= ~7`) for small longest-AP and SPIKES to `~k` only when `L` is large (near-AP).
This is exactly kps-S90's interlock (longest-AP `<= 8` pigeonhole vs `>= 9` density floor), here made
quantitative for `j*` directly.  The universal bound `j* <= ceil(7(k-1)/6) = O(k)` holds in all data.

## (3) THE DICHOTOMY: no-good-period ⟺ tight AP (the safety check)

The worry: does a good period ALWAYS exist in the hard region?  NO — but the exceptions are exactly
the tight cases.  Census over **~490,000** `(E, Vmax)` configs, `k = 8..13`:

> **The ONLY clusters with NO good period are the exact full complete-residue APs `{0,1,…,k-1}` at
> `Vmax = k`, and only for PRIME `k` (`k=11, 13`).  ZERO non-AP exceptions.**

Why: for `k = Vmax` prime, every `j ∈ {1,…,k-1}` is coprime to `k`, so `{e_i j mod k}` is a
permutation of ALL residues `{0,…,k-1}` — evenly spaced, `maxgap = 1/k < 1/7`.  No good period.  For
composite `k`, a `j` sharing a factor with `Vmax` collapses to a sub-lattice with a `>1/7` gap; for
`k < Vmax` (or any dilated/defected AP), a good period exists.  These `k = Vmax` prime full APs are
precisely the **tight `M = 1/k` LRC instances** (the AP `{1,…,k}`), handled by the LRC(`<=13`)
citation / the tight-locus analysis (`= {AP, GW}`), NOT the good-period route.

**Consequence.**  The good-period lemma is a clean DICHOTOMY:
> **every hard-region cluster has a good period at `j* = O(k)`  OR  it is a tight complete-residue AP
> (cited).**  There is no hidden counterexample — the `0`-counterexamples mac-mini observed is
> explained: the only exceptions are the tight APs, which are excluded from this leg by construction.

## The proof route this exposes

The capstone `j* = O(k)` (for non-tight clusters) splits cleanly on the longest-AP axis:
- **small `L`** (dissociated cluster): `j*` small — the arc-count pigeonhole (kps-S90) is
  non-vacuous because the resonance count is bounded by the AP cap; a good period appears at
  `j* <= O(1/theta) = O(7)`.
- **large `L`** (near-AP, `L` close to `k`): the cluster is `AP + o(k)` defects; the good period comes
  from the embedded AP's Dirichlet structure (mac-mini's `j* <= ceil(7(k-1)/6)`), and the defect only
  HELPS (it breaks the exact-AP full-residue rigidity that is the *sole* no-good-period case).
- **exact AP** (`L = k`, full residues, prime `Vmax`): the tight `M=1/k` case, cited.

So the residual reduces to: "small-`L` ⟹ `j* <= O(7)`" (a bounded-resonance pigeonhole) + "near-AP
(`L` large, not exact) ⟹ `j* <= ceil(7(k-1)/6)` via the embedded AP."  Both are `O(k)`; the tight AP
is cited.  This is the dilation-invariant unification of mac-mini's AP bound and kps-S90's interlock.

## Ledger / next

- FOUND: AP is not the `j*`-maximizer; `j*` is longest-AP-governed; the no-good-period set = exactly
  the tight complete-residue APs (`k=Vmax` prime), 0 non-AP exceptions over 490k configs.
- REFRAMED: the capstone `j*=O(k)` as [small-`L` pigeonhole `O(7)`] + [near-AP embedded-AP `O(k)`] +
  [exact AP cited], on the dilation-invariant longest-AP axis.
- NEXT: prove "small-`L` ⟹ `j* <= O(7)`" (bounded-resonance pigeonhole, the genuinely new quantitative
  piece) and "near-AP ⟹ embedded-AP Dirichlet" rigorously (mac-mini owns the AP bound; kps owns the
  interlock — this hands them the exact exception set + the axis).
- Files: `lrc14_jstar_AP_extremal_opus_S164.py`, `lrc14_jstar_extremal_structure_opus_S164.py`,
  `lrc14_nogoodperiod_dichotomy_opus_S164.py` (+outs).
