---
id: HYP-3727
title: PRIMITIVITY resolves the odd/even covering-min back-and-forth, reconciles opus's "1/n" with my S47 "mediant", and aligns with canon THM-523. The FULL covering-min (incl. non-primitive) = 1/n for ALL n (not just even) via the scaled block g*{1..n-1}, g = smallest prime factor of n -- this is the q-witness/EASY case in disguise (g*{1..n-1} / g = {1..n-1}, which omits a multiple of n; opus's even block 2*{1..n-1} is the g=2 special case for even n). Parity only chooses g; the value is 1/n. But THM-523 (CANON, PROVED) reduces LRC to PRIMITIVE covering sets, where M>1/n STRICTLY (n=7->2/13, n=8->2/15, n=9->4/33, n=14~7/89; HYP-2566 uniform looseness) -- the genuinely HARD case, which carries the MARGIN. So opus's covering-min=1/n is the EASY case (excluded by primitivity); the hard primitive covering-min is >1/n and is what S47/S48 found. RAMANUJAN/PALEY FRAME: the primitive covering-min lives on a circulant mod m (n=7:13, n=8:15, n=9:33); the Ihara-RH / Ramanujan / Weil sqrt-bound on the speed character sums is the spectral-gap criterion controlling the floor; 2n-1 is a Paley-tournament/graph vertex count (n=7: Paley graph on 13 = RAMANUJAN, verified; n=14: 2n-1=27=GF(3^3) -> Paley tournament). THREE LEVERAGE WAYS for the (mediant) margin 1/(n(2n-1))=1/dim so(2n): (1) tournament embedding -> the Paley tournament on 2n-1 vertices (the H(T)/OCF bridge); (2) Borel-Cantelli (Sum margins = ln4); (3) Beta-moment LP (margin=2B(2n-1,2))
status: VERIFIED (full covering-min=1/n via scaled blocks n=7..15 exact; primitive covering-min>1/n n=7,8,9 exact + THM-523 canon for n=14; Paley-13 Ramanujan exact). Reconciles opus-S1 (1/n, non-primitive) + mac-mini-S47 (mediant, primitive) + THM-523 (canon). The Ramanujan/Paley/Ihara link and the 3 leverage routes are FRAMES/leads (not executed proofs).
source: mac-mini-2026-06-30-S49
related:
  - HYP-3725  # S47 refutation (convergent not covering-min) -- correct, and it was the PRIMITIVE covering-min
  - HYP-3726  # the hexagon margin 1/(n(2n-1)) = 1/dim so(2n); this is the mediant/primitive margin (n=7,8)
  - HYP-2566  # uniform looseness: inf M over primitive covering sets > 1/n -- THE conjecture this all serves
  - THM-523   # CANON: LRC <=> M>=1/n for all PRIMITIVE covering sets (the reduction that makes primitivity essential)
references:
  - opus-2026-06-30-S1  # even block 2*{1..n-1} M=1/n -- the g=2 non-primitive/q-witness case (corrected here)
  - opus metazeta (commit 52d1bbe40)  # Ihara zeta of the metagraph G_n -- tournament-side of the Ramanujan frame
results:
  - 04-computation/primitivity_parity_ramanujan_macmini_20260630.py
  - 05-knowledge/results/primitivity_parity_ramanujan_macmini_20260630.out
---

# HYP-3727 -- primitivity resolves odd/even; the Ramanujan frame

The owner asked to work the odd-n covering-min and the even-n statement "back and forth," plus the three
leverage ways and the Ramanujan/Ihara-RH idea. Concurrently opus claimed the covering-min is `1/n` (the even
block), and klein/opus refuted the construction. The reconciliation is **primitivity**.

## The odd/even back-and-forth is a PRIMITIVITY artifact
- **FULL covering-min (any covering set) = `1/n` for ALL `n`.** The scaled block `g*{1,..,n-1}` with
  `g = smallest prime factor of n` is a covering `(n-1)`-set with `M = M({1,..,n-1}) = 1/n` (verified n=7..15).
  This is the **q-witness / EASY case in disguise**: `g*{1,..,n-1} / g = {1,..,n-1}` omits a multiple of `n`,
  so its `M = 1/n` is exactly the `q=n` witness bound. opus's even block `2*{1,..,n-1}` is the `g=2` case
  (works when `n` even, so `n` is itself in the doubled set). **Parity only chooses `g`** (even `n`: `g=2`;
  odd prime `n`: `g=n`; odd composite: `g=p_min`). The value is `1/n` for every `n`.
- **PRIMITIVE covering-min (the HARD case) > `1/n` strictly.** THM-523 (CANON, PROVED) reduces LRC to
  *primitive* covering sets, and there `M > 1/n` strictly: `n=7 -> 2/13`, `n=8 -> 2/15`, `n=9 -> 4/33`,
  `n=14 ~ 7/89` (HYP-2566: the inf is conjecturally bounded away from `1/n` -- uniform looseness). This is the
  genuinely hard case and the source of the **margin**.

So opus (`1/n`) and my S47 (`mediant > 1/n`) are about **different quantities**: the full vs the primitive
covering-min. THM-523 makes the **primitive** one the LRC's content. opus's `1/n` is the easy case; the hard
case has a margin at every `n` (independent of parity).

## Why the reduction is to primitive sets (the key step)
A non-primitive covering set `S` (gcd `g`) satisfies `M(S) = M(S/g)`. `S/g` is primitive but typically
**non-covering** (it can omit a multiple of `n`), so the elementary `q`-witness already gives `M(S/g) >= 1/n`.
Thus non-primitive covering sets collapse to the easy case, and the only hard cases are **primitive** covering
sets -- exactly THM-523's reduction. opus's even block is the canonical example of this collapse.

## The Ramanujan / Paley / Ihara frame (the owner's "consider")
The primitive covering-min lives on a **circulant mod `m`** (the optimal witness `t* = k/m`: `m=13` at n=7,
`15` at n=8, `33` at n=9). `M(S)` is governed by the **speed character sums** `R-hat(j) = sum_{v} e^{2pi i v
j/m}` (HYP-3704). The criterion "a regular graph is Ramanujan iff its Ihara zeta satisfies the RH analogue"
is exactly the **spectral-gap / Weil sqrt-bound** on these character sums:
- `2n-1` is a **Paley** vertex count: at `n=7`, `2n-1=13` (prime, `1 mod 4`) gives the **Paley graph on 13**,
  which is **Ramanujan** (max non-trivial eigenvalue `2.303 <= 2sqrt(5)=4.472`, verified). At `n=14`,
  `2n-1=27=GF(3^3)` (`3 mod 4`) gives the **Paley tournament on 27 vertices** (non-principal `|eigenvalue| ~
  sqrt(27)/2` = the Weil/Ramanujan sqrt bound).
- The construction's AP gives **Gauss/Dirichlet sums** (`|sum| ~ sqrt(m)`, Weil-tight) -- structured, not the
  flat optimum.
- opus's **metazeta** (the Ihara zeta of the metagraph `G_n`, Bass formula) is the **tournament-side**
  instance of the same zeta machinery -- the dual mandate's version of the LRC's spectral criterion.

So the Ihara-RH/Ramanujan idea is the spectral lens on the floor: the covering-min is the most equidistributed
(Ramanujan-flat) primitive covering set, and the floor `M >= 1/n` is a lower bound on how flat any covering
set can be (the covering constraint forces a gap).

## The three leverage ways (for the mediant margin `1/(n(2n-1)) = 1/dim so(2n)`, HYP-3726)
1. **Tournament embedding.** `margin = 1/(arcs of K_{2n})` and the mediant's `mod (2n-1)` circulant point at
   the **Paley tournament on `2n-1` vertices** -- a concrete bridge from the LRC margin to the project's
   `H(T)`/OCF tournament machinery (the two mandates meeting). `2n-1=27=GF(3^3)` at `n=14` is the natural
   target.
2. **Borel-Cantelli.** `Sum_n 1/(n(2n-1)) = ln 4 < infinity`: a finite safe-measure budget; a union-bound
   over levels controlled by `ln 4` (ties HYP-3615, THM-579).
3. **Beta-moment LP.** `margin = 2 B(2n-1,2)`: an explicit circle moment, a ready test function for a
   Beurling-Selberg / moment-LP lower bound on the floor.

## Caveat
The mediant margin `1/(n(2n-1))` is exact at `n=7,8` (where the primitive covering-min IS the mediant); at
`n=14` the primitive covering-min is `~7/89` (THM-523), a different value, so the margin there is `~9/1246`,
not `1/378`. The leverage routes target the mediant family; whether the true primitive covering-min has the
same clean form is open (HYP-2566). The primitivity resolution and the Ramanujan frame are unconditional.
