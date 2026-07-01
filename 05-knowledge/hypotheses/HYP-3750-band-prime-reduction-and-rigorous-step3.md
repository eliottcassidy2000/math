---
id: HYP-3750
title: TWO rigorous results toward the LRC14 lower bound. (A, MINOR THEOREM PROVEN) THE BAND-PRIME REDUCTION -- LRC14 <=> M(S)>=1/14 for every covering 13-set that is a +-transversal-OR-multiple modulo EACH of the band primes {17,19,23}; every other covering set has M>=2/23>1/14 (rigorous corollary of the klein-S42 radius-1 witness theorem). The band primes are EXACTLY the primes in (n, 2n*183/14 ~ 26.1] that exceed n=14; primes p<=13 are satisfied FREE (covering => a multiple of p). (B, STEP 3 MADE RIGOROUS) klein-S45's hand-wavy budget step ("the n-1 budget leaves ~1 slot, forcing a CRT band-coverer", verified by exhaustive search) is replaced by: (i) the (T_p) reformulation [rigorous]; (ii) PATCH IDENTIFICATION -- a missing core speed k breaks the pair {k,p-k} mod p for p in (12+k,25], whose only patch is a speed =+-k mod p that is LARGE (>= p-k >= 13) [rigorous]; (iii) the CARDINALITY FLOOR -- no multiple of 23 => >=11 speeds in distinct pairs mod 23 (each speed lies in exactly ONE pair), so |S|>=11 with <=2 spare, TIGHT at the construction (residues 1..12,21 mod 23) [rigorous]. THEN the decisive negative result: a SINGLE huge CRT speed (w=0 mod 182, w=+-k mod 17*19*23) satisfies ALL large-speed obligations at once, so the cardinality lower bound on |S| is EXACTLY 11 and never reaches 14 -- COUNTING CANNOT CLOSE STEP 3. The entire residual is therefore Step 4 (does that CRT patch dig an M-hole) = HYP-3745 (perturbation-proved) / multi-family inexhaustibility (HYP-3749) = the LRC14 lower bound. NET: removes the exhaustive search from Steps 1-3 (now a clean reformulation) and PROVES the difficulty is residue-structural, not budgetary -- localizing 100% of LRC14 to one statement about the +-k CRT patch.
status: TASK A (reduction) PROVED -- rigorous corollary of the witness theorem; verified (243 not-(T_p) covering sets, 0 violations; construction is a triple band-transversal). TASK B: the (T_p) reformulation, patch identification, and cardinality floor are RIGOROUS; the tightness (cardinality min = 11, one CRT speed serves all large obligations) is PROVED by explicit CRT construction -- so the residual is provably NOT cardinality but Step 4. Does NOT prove LRC14; it REDUCES it (cleanly) and localizes the residual.
source: mac-mini-2026-06-30-S61
depends_on:
  - HYP-3741   # klein-S42 radius-1 witness theorem (PROVED): no runner in {-1,0,1} mod p at some rotation => M>=2/p
related:
  - HYP-3747   # klein-S45 the full lowness lemma 4-step chain (this makes its Step 3 rigorous + localizes the residual to Step 4)
  - HYP-3745   # klein-S44 CRT escape uncoverable / fused-radius trap (= Step 4, the now-sole residual; perturbation-proved)
  - HYP-3749   # mac-mini-S60 multi-family inexhaustibility (= Step 4 for arbitrary patch = the residual)
  - HYP-3740   # the lowness lemma (mac-mini, exhaustive-search verification this partially de-randomizes)
results:
  - 04-computation/lrc14_band_reduction_step3_rigorous_macmini_20260630.py
  - 05-knowledge/results/lrc14_band_reduction_step3_rigorous_macmini_20260630.out
---

# HYP-3750 -- the band-prime reduction (minor theorem) + Step 3 made rigorous

Two rigorous deliverables toward the LRC14 lower bound, both flowing from one clean corollary of the
klein-S42 radius-1 witness theorem (PROVED, HYP-3741).

## The (T_p) reformulation (rigorous foundation)
**Witness theorem (klein-S42, PROVED).** If some rotation `a in (Z/p)*` has no speed `s in S` with
`s a in {-1,0,1} mod p`, then `M(S) >= 2/p` (witnessed by `t=a/p`).

**Corollary (T_p).** For a prime `p`: `M(S) < 2/p  =>  S has a multiple of p  OR  {+-s mod p : s in S}`
covers every nonzero pair of `(Z/p)*` (a **±-transversal**). *Proof.* If some unit `u` has neither
`p | s` nor `s ≡ ±u` for any `s`, take `a = u^{-1}`: then no `s` has `sa in {-1,0,1}`, so `M(S) >= 2/p`. `□`

This is elementary and rigorous (no search). `M(S) <= 14/183 < 2/p` for every `p <= 23`, so it applies to
all of `p in {2,...,23}`.

## TASK A -- THE BAND-PRIME REDUCTION (a minor related statement, proved)
The **band primes** are exactly the primes `p` with `2/p > 14/183` (`<=> p <= 26.14`) that exceed `n=14`:
`{17, 19, 23}`. For primes `p <= 13`, a covering set automatically has a multiple of `p`, so `(T_p)` is
**free**. Hence:

> **Reduction.** A covering 13-set `S` that is NOT a `±`-transversal-or-multiple modulo some
> `p in {17,19,23}` has `M(S) >= 2/p >= 2/23 = 0.087 > 1/14`. Therefore
> **LRC14 `<=>` `M(S) >= 1/14` for every covering 13-set that is a `±`-transversal-or-multiple modulo
> EACH of `17, 19, 23`.**

Every other covering set already satisfies LRC14 with room to spare. The hard subclass is non-empty and
contains the extremal candidate: the construction `{1,...,12,182}` is a triple band-transversal (residues
`182 ≡ 12, 11, 21 (mod 17,19,23)`; no multiples). Verified: 243 covering sets failing some `(T_p)`, all had
`M >= 2/p`.

## TASK B -- klein-S45 STEP 3 MADE RIGOROUS (removing the exhaustive search)
klein-S45's Step 3 read: *"the small core-minus-k (`n-3`) + killers (`~2`) exhaust the `n-1` budget, leaving
`~1` slot, forcing a single CRT band-coverer,"* certified for `n=14` by mac-mini's exhaustive search
(HYP-3740). Replace it by three rigorous facts plus one decisive negative result.

1. **(T_p) reformulation** (above) -- rigorous, no search.
2. **Patch identification** (rigorous). If core speed `k` is missing, the pair `{k, p-k} mod p` is broken for
   every `p in (12+k, 25]`. Its only patch is a speed `≡ ±k mod p`; the smallest non-missing representative is
   `p-k >= 13`, so **every patch is a LARGE speed `>= 13`** (or a multiple of `p`). (Table matches klein-S45
   Step 1 exactly: `k=1,3 -> {17,19,23}`, `k=6 -> {19,23}`, `k=10 -> {23}`, `k=12 -> {}`.)
3. **Cardinality floor** (rigorous). If `S` has no multiple of 23, then `±`-transversality mod 23 requires all
   11 pairs hit, and **each speed lies in exactly one pair mod 23**, so `|S| >= 11` -- at most **2 spare**.
   TIGHT at the construction (residues `1..12, 21 mod 23`: 11 pairs, exactly 2 redundant).

4. **The decisive negative result -- counting CANNOT finish.** For missing `k <= 4` (no multiples of the band
   primes) the large-speed obligations are `O1: 13|L`, `O2: 14|L`, `O3..O5: L ≡ ±k mod 23, 19, 17`. A **single**
   CRT speed `w ≡ 0 mod 182`, `w ≡ ±k mod 17·19·23` satisfies **all of them simultaneously** (explicit
   `w` printed for `k=1..4`). So the minimum number of speeds meeting `O1..O5` is **1**, and the cardinality
   lower bound on `|S|` stays at the mod-23 floor `11` -- it **never reaches 14**. Therefore **Step 3 cannot be
   closed by any counting/LP/budget argument** (this is why an exhaustive search was used).

## No finite-witness certificate either (a second negative result)
One might hope to close the near-construction (`core + <= 2` extra speeds) by a **finite** set of small-modulus
witnesses, since `min_v ||v a/D||` depends only on residues mod `D`. It fails: the bare core `{1..12}\{k}` has
30-60 witnesses at `D <= 30`, but a greedy hitting-set shows **just 2 speeds suffice to kill all of them** (for
every `k`). So no `D <= 30` certificate proves the lowness lemma; the binding witness of a well-chosen
completion sits at `D > 30` ("the hole moves but never vanishes," klein-S44) and the required `D` is
**unbounded**. This confirms the residual is irreducibly about **unbounded moduli** -- exactly the multi-family
inexhaustibility (HYP-3749) -- not any finite check.

## What this buys -- the residual is localized
The entire remaining difficulty is whether the forced `±k` CRT patch digs an `M`-hole above `14/183` -- i.e.
**Step 4** (klein-S45), `=` HYP-3745 (CRT-escape uncoverable, perturbation-proved) `=` multi-family
inexhaustibility (HYP-3749) for an arbitrary patch. Steps 1-3 are now a clean reformulation with **no search**;
the difficulty is **residue-structural, not budgetary**, and 100% of LRC14 sits in one statement:

> *The `±k` CRT band-patch forced by a missing core speed satisfies `M(S) >= 2/(2n-3) > n/Phi_6`.*

This both gives a genuine **minor theorem** (the band reduction) and tells the team precisely where to spend
effort: **not** on the budget (provably hopeless) but on the patch's hole (Step 4 for general `S`).

## Honest scope
Task A is a complete rigorous reduction (corollary of the proved witness theorem). In Task B, items 1-3 are
rigorous and item 4 is a proof (explicit CRT witness) that cardinality is insufficient. None of this proves
LRC14; it removes the exhaustive search from Steps 1-3 and pins the sole residual to Step 4. The construction
being the *consecutive* transversal (not just *some* transversal) remains the spread-regime gap inside Step 4.
