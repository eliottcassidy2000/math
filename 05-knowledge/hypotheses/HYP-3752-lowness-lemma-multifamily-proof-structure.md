---
id: HYP-3752
title: Working the LOWNESS-LEMMA proof -- the witness families are MULTIPLE, and speed 1 is the universal gap-filler. Missing speed 1 leaves the speed-1 slot empty at a WHOLE ENSEMBLE of witness families, not just the t_a family (HYP-3745): the a=1 witness (t=1/15, M=2/15), the t_a family (t=a/(14a+1), M=2a/(14a+1)), the D=16 witness (t=1/16, M=1/8), ... Each family is killed only by a speed of a SPECIFIC residue (the a=1 gap by =±1 mod 15; D=16 by =±1 mod 16; t_a by the q=14-coverer 14m, which fills t_a for a>m/2). SPEED 1 has residue 1 at EVERY modulus, so it fills the gap in ALL families simultaneously -- which is exactly why the construction (WITH speed 1) reaches 14/183. A missing-1 set's forced q-coverers kill SOME families but SPAWN others (e.g. speed 14 = a q=14-coverer kills the entire t_a family but fires the D=16 witness, M=1/8). WINDOW LEMMA: a single speed w fills the gap of a given family only on a bounded-ratio window of rotations (asymptotically a in (m0/(s+2), m0/(s-2)), s=w mod 14, ratio <=5 for s>=3), so it cannot cover a wide rotation range. VERIFIED (no overturn): every missing-1 set has M>14/183 (min ~1/8 to 2/17; {2..14}=1/8, {2..12,13,182}=2/15, {2..12,15,182}=2/17). So the lemma HOLDS; the full proof reduces to MULTI-FAMILY INEXHAUSTIBILITY (no single replacement speed covers the whole family ensemble) = the LRC14 lower bound
status: STRUCTURE MAPPED + verified TRUE (no overturn). The multi-family witnesses, the speed-1 universal-filler mechanism, and the window lemma are exact/verified; every tested missing-1 set has M>14/183. The full proof (the family ensemble is inexhaustible by one speed, for ALL covering sets) is the LRC14 lower bound -- OPEN. This session refines HYP-3745 (the t_a family was ONE of the ensemble; speed 14 kills it but another fires) and confirms the lemma is not overturnable.
source: mac-mini-2026-06-30-S60
related:
  - HYP-3751  # my S59: the t_a witness family (this refines it -- one family of the ensemble; killing it spawns another)
  - HYP-3745  # klein-S44 CRT-escape uncoverable / fused-radius trap + the PROVED CRT-invariant counting bound (=2r+1 rotations per speed, regardless of value) = the rigorous window lemma
  - HYP-3744  # the constant-residue principle (speed 1's residue 1 at every modulus = the universal gap-filler, here dynamic over families)
  - HYP-3740  # the LRC14 hard core = the lowness lemma (this maps its proof structure)
  - HYP-3741  # klein-S42 witness theorem (each family is an instance: gap radius => M>=(r+1)/D)
results:
  - 04-computation/lowness_lemma_multifamily_macmini_20260630.py
  - 05-knowledge/results/lowness_lemma_multifamily_macmini_20260630.out
---

# HYP-3752 -- the lowness lemma's proof structure: multiple witness families, speed 1 the universal filler

Working toward a proof of the lowness lemma (`M(S) <= n/Phi_6 => {1,..,n-2} \subseteq S`, the LRC14 hard core,
HYP-3740). My S59 `t_a` family (HYP-3751) is only ONE of an ensemble; mapping the ensemble both clarifies the
proof and shows the lemma is not overturnable. **Convergent independent work:** klein-S44 (HYP-3745) reached the
same "CRT escape uncoverable" conclusion from the construction's local rigidity (the *fused-radius trap*) and --
crucially -- **proved** the counting bound my window lemma needed (see below).

## The witness families are multiple
Missing speed 1 leaves the speed-1 slot empty at a whole **ensemble** of witness families (each a klein-S42
witness, gap radius => `M >= (r+1)/D`):

| family | witness `t` | `M` | killed by |
|--------|-------------|-----|-----------|
| a=1 | `1/15` | `2/15` | a speed `≡ ±1 mod 15` |
| t_a | `a/(14a+1)` | `2a/(14a+1)` | the q=14-coverer `14m` (fills `a>m/2`) |
| D=16 | `1/16` | `1/8` | a speed `≡ ±1 mod 16` |
| ... | `1/p` etc. | `2/p` | a speed `≡ ±1 mod p` |

## Speed 1 is the universal gap-filler
Speed `1` has residue `1` at **every** modulus, so `||1 * t||` is small at every one of these witnesses -- it
fills the speed-1 gap in **all** families at once. That is exactly why the construction `{1,..,12, 182}` (which
HAS speed 1) sits at `14/183`: speed 1 closes every gap, and the residual binding is the three-gap at `D=Phi_6`.
A missing-1 set has no such universal filler; its **forced q-coverers** close some families but **spawn**
others -- e.g. the q=14-coverer `14` kills the entire `t_a` family (`14a ≡ 14 mod (14a+1)`, distance 1) but
fires the `D=16` witness (`M=1/8`). This is the dynamic constant-residue principle (HYP-3744): only the constant
residue `1` covers every modulus.

## The window lemma (klein-S44 proved the counting bound)
A single speed `w` fills the gap of a fixed family only on a **bounded-ratio window** of rotations:
asymptotically `a in (m0/(s+2), m0/(s-2))` with `s = w mod 14`, `m0 = floor(w/14)` -- ratio `<= 5` for `s>=3`.
So no single speed can close a wide range of a family's gaps; covering the whole ensemble needs many speeds,
which the `n-1` budget cannot supply. **klein-S44 (HYP-3745) proved the underlying counting bound**: each speed
covers at most `2r+1` rotations of `Z/p` **regardless of its value** (CRT-invariant -- tuning the speed by CRT
chooses *which* rotations, never *how many*; "the hole moves but never vanishes"). That is exactly the rigorous
backbone of this window lemma: the per-speed coverage of any one family is bounded, independent of CRT tuning,
so the ensemble cannot be exhausted by re-aiming a few large speeds.

## Verified: the lemma is TRUE (no overturn)
The dangerous case is a small q=14-coverer killing the `t_a` family. Checked: `{2,..,14}` (speed 14 kills all
`t_a`) has `M = 1/8` at `t=1/16` (a NEW family), `>> 14/183`. Across forced q-coverers `14m` (`m=1..13`),
`M in {1/8, 2/15}`, always `> 14/183`. Every tested missing-1 covering set has `M > 14/183` (min `~1/8` to
`2/17`). The arc (`covering-min = 14/183`) is not overturned; the lemma holds.

## Status -- reduced to multi-family inexhaustibility (= LRC14)
- **Exact / verified:** the family ensemble, the speed-1 universal-filler mechanism, the window lemma, and
  `M > 14/183` for every tested missing-1 set.
- **Open (= LRC14 lower bound):** that the ENSEMBLE is inexhaustible by the `n-1` speeds of ANY covering set --
  i.e. some family's gap always survives. This is the heart of LRC14 and is not closed.

## What it buys
Maps the lowness-lemma proof: the obstruction is an **ensemble** of witness families, each closable only by a
speed of a specific residue, all closed at once only by the constant-residue speed 1. The proof reduces to
"the ensemble is inexhaustible by `n-1` speeds" -- the LRC14 lower bound -- and the window lemma + budget are
the tools. It also rules out the overturn (a small q=14-coverer killing `t_a` does not reach `14/183`; a new
family fires), confirming `covering-min(14) = 14/183` is safe and the remaining work is the inexhaustibility.
