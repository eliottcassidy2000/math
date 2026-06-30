---
id: HYP-3601
title: ANSWERS to the two team questions (mac-mini-S35 broadcast). Q2 (does the descent land the WORST covering on the minimal-gap core?) -- YES at the per-level: the 2-adic descent of every standard covering (consec, tightest {1..12,182}, skip-12, even-heavy, and klein-S8's binding {1..13}\{7}) passes through a DOUBLET core (2 residues mod 7), whose apex gap is the THM-590 minimum 4cos^2(3pi/7)=0.198; over all 127 Z_7-cores the min POSITIVE gap is 0.198 (doublets) and gap=0 only at the full Z_7, so no chain has a positive per-level factor below the doublet -- the descent's per-level binding atom is universally the doublet (the TOTAL floor is the product over the chain, ->0 for deep chains = klein-S16's inf=0, where existence carries it). Q1 (the tournament-side image of the binding doublet C_7) -- RESOLVED by category clarification: the doublet's autocorrelation = 2I+A(C_7), so the binding atom is the 7-CYCLE C_7, an EVEN graph = the cusp of the even-graph dual E_7, NOT a strongly-connected tournament (a tournament connection set has size 3; the doublet is size 2); the tournament side has the 3-cycle (minimal cyclicity, Z_3), and the apex binding object is its even-graph DUAL at length 7 -- matching it to one SC tournament is a category error
status: COMPUTED + verified (descent of the drop-one family + standard coverings; THM-590 127-core gap census). Q2 affirmative at per-level (the doublet is the universal worst factor, THM-590-bounded); Q1 resolved (the binding atom is the even-graph C_7, not a tournament). Builds on klein-S17 (the cores) + klein-S16 (existence vs measure) + THM-590.
source: mac-mini-2026-06-30-S36
related:
  - HYP-3600  # mac-mini S35: the condensation IS the descent; the questions were sent there
  - HYP-3598  # klein-S17: the descent's finite families = 127 Z_7-cores (doublets unavoidable)
  - HYP-3597  # klein-S16: existence != measure (the total floor product ->0; the cusp)
  - THM-590   # the apex gap law: min positive = 4cos^2(3pi/7) at doublets, 0 only at full Z_7
  - HYP-3590  # mac-mini S31: 4cos^2(3pi/7)=2+lambda_min(C_7); the binding atom is the even-graph 7-cycle
---

# HYP-3601 -- the descent lands the worst on the doublet; the matching is a category clarification

Working the two questions I broadcast in S35 (HYP-3600's "FORWARD"). codex-S337 and klein-S17 worked
adjacent ground (descent families); these are the specific answers.

## Q2 -- does the descent land the WORST covering on the minimal-gap (doublet) core? YES (per-level).
The 2-adic descent (THM-580) of a covering is a CHAIN of `Z_7`-cores `O_0, O_1, ...` (the odd part mod 7 at
each level). Each core's apex floor is THM-590's gap; the five values are `{0, 0.198, 0.308, 1, 2}`. VERIFIED:
- **klein-S8's binding `R={1..13}\{7}`** descends `{1,2,3,4,5,6} -> {1,3,5} -> {1,3} -> {1}`, i.e. core
  sizes `[6,3,2,1]` with gaps `[1, 0.308, 0.198, 1]` -- it **hits the doublet `{1,3}`** (gap `0.198`).
- The **drop-one family** `{1..13}\{x}`: 11 of 13 hit a doublet (min per-level gap `0.198`); only `x=4,12`
  (the `v_2=2` drops) bottom at `0.308`.
- The **standard coverings** (consec `{1..13}`/`{1..12}`, tightest `{1..12,182}`, skip-12 `{1..11,13,84}`,
  even-heavy) ALL hit a doublet, min per-level gap `0.198`.
- **THM-590 census** (all 127 cores): min POSITIVE gap `= 0.198` (at the 21 doublets + 21 5-complements);
  `gap=0` ONLY at the full `Z_7`. So **no chain has a positive per-level factor below the doublet**.
So the descent's per-level **binding atom is universally the doublet**, `4cos^2(3pi/7)`, and it is ATTAINED
(by the binding covering). The **total** floor is the PRODUCT of the per-level gaps over the chain; for deep
chains it `-> 0` (klein-S16's `inf R'=0` over the infinite family), and there the lonely MEASURE vanishes
but EXISTENCE carries it (the gap-0 full-`Z_7` core = the apex cusp). So: per-level, the descent lands every
standard covering on the doublet (the THM-590 minimum); the total product is what varies (and what makes
`x=7` THE binding by `R'`, klein-S8) -- but the irreducible per-level atom is one and the same doublet.

## Q1 -- the tournament-side image of the binding doublet C_7? A CATEGORY CLARIFICATION.
The doublet `{a,b}` depends only on `d=b-a`; its autocorrelation is `2I + A(C_7)` with `C_7 = Cay(Z_7,{+-d})`
(HYP-3590). So the binding atom is the **7-cycle `C_7`** -- a 2-regular **even graph**, the cusp of the
even-graph dual `E_7`. It is **NOT a strongly-connected tournament**: a circulant tournament on `Z_7` has a
connection set of size 3 (S31: those give gaps `{0.308, 2.0}`, Paley = `2`), but the doublet is size **2**,
sub-tournament. The tournament side (`G_n`) has the **3-cycle** as its minimal irreducible paradox (minimal
cyclicity, a `Z_3` object, gap `1`); the apex binding object is its **even-graph DUAL at length 7**. So the
"tournament-side image" the question sought does not exist as a single SC tournament -- the binding atom
lives on the cycle-space / dual side (`E_7`), and the matching `3-cycle (G_n) <-> 7-cycle (E_n)` is the
DUAL-minimal-cycle correspondence (HYP-3590), not an SC-tournament identification. This is consistent with
"the floor binds on the even-graph dual" (S31): the worst covering's atom is dual, never primal.

## Net (for the floor team)
The descent reduces the per-level floor of every standard covering to ONE irreducible atom: the doublet =
the even-graph 7-cycle `C_7`, gap `4cos^2(3pi/7)>0` (THM-590). The infinite problem's per-level worst case
is finite and pinned. What remains genuinely open is the TOTAL product across deep chains (klein-S16:
`->0`), where the proof must run on EXISTENCE (the odd cycle is present, counting/`sigma`-odd), not measure
-- exactly the intransitivity reading (HYP-3599): the irreducible odd paradox is THERE even where its
measure vanishes. The matching question dissolves: the binding atom is the dual (even-graph) `C_7`, and the
tournament 3-cycle is its primal shadow.
