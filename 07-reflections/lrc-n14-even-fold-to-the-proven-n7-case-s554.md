---
source: opus-2026-06-01-S554 (remote-control)
status: rigorous reduction (half the configs) + strong computation; NOT a full proof of LRC@14
tags: [LRC, n14, n7, even-fold, parity, LRC7-lever, doubling, antipodal, S520, S554]
---

# LRC@14 via the proven n=7 case: the even-fold lever ("the 7 impossibility")

**Prompt (user):** consider the 7-impossibility the key to the proof of LRC n=14.

**Reading of the hint.** `n=14 = 2·7` and **LRC(7) is a theorem** (Barajas–Serra) —
the "impossibility" of a 7-runner counterexample. There is an exact bridge from
n=14 to n=7 through the factor 2, and it is the natural lever.

## 0. The exact bridge (doubling)

An even speed `v = 2u` satisfies `||v t|| = ||u·(2t)||`. With `fold(S) = {v/2 :
v∈S even}` and `s = 2t`, the even runners' collar is the n=7-collar of the fold:
`g_even(t) = min_{u∈fold} ||u s|| = g_fold(2t)`. Hence, for the max-collar
`M(S)=max_t min_i ||v_i t||`,

> **(*)  `M14(S) ≤ M(fold(S))`.**   (verified exactly, 0 violations / 224 sets)

## 1. LRC(7) always engages — `min(e,o) ≤ 6`

There are 13 runners (an **odd** count), so one parity class has `≤ 6` members:
`min(#even, #odd) ≤ 6`. Thus LRC(7) is *always* applicable to the minority
parity. Two regimes (≈51% / 49% of random primitive 13-sets):

- **e ≤ 6 (even-minority): clean fold.** By LRC(7), `M(fold) ≥ 1/7`, so the even
  runners are `≥ 1/7 (≥ 1/14)` safe on `E_good = {t : g_fold(2t) ≥ 1/7}`, a set of
  **positive measure** (the doubling-preimage of the n=7 safe set). The even half
  is *fully settled by the proven case*.
- **o ≤ 6 (odd-minority):** LRC(7) gives an interval where the ≤6 odd runners are
  `≥1/7`; but the even half then has `≥7` speeds and dodging it needs LRC(8+)
  (open). Not closed here.

## 2. The residual is the ODD interleaving (antipodal preimages)

For e ≤ 6, all that remains is: find `t∈E_good` with the (≥7) odd runners also
`≥ 1/14`. The two doubling-preimages of an even-good `s`, namely `s/2` and
`(s+1)/2`, place every **odd** runner at **antipodal** positions (differ by 1/2),
so each odd runner is unsafe (`<1/14`) at **at most one** preimage. Therefore a
single even-good `s` yields a full n=14 witness **unless the odds "split"** (some
unsafe at `s/2`, others at `(s+1)/2`). So:

> **LRC(14)[e ≤ 6]  ⟸  some even-good `s` has no odd-split.**

The difficulty is neither parity alone — the odd part `{1,3,…,13}` has collar
`1/2`, the even fold has collar `≥1/7` — but their **coupling**: the even-optimal
`s` is exactly where the odds tend to split, pinning `M14` at `M(fold)/2`.

## 3. The reduction works (computation)

`lrc_n14_even_fold_to_n7_s554.py`:
- **(*) holds exactly**, 0/224 violations.
- **Fold reduction succeeds 127/127** on the e≤6 configs: an even-good `s` whose
  better preimage is a full n=14 witness was found in *every* case (random,
  structured, AP, V*). So the e≤6 case of LRC(14) is, empirically, reduced to the
  proven LRC(7) plus a no-odd-split choice.
- Every e≥7 sample is lonely (0 would-be counterexamples), consistent with LRC.

## 4. The tight configs through the lens

The two known n=14 tight configs both have `e = 6` and the **same** 7 odd runners
`{1,3,5,7,9,11,13}` (there are exactly 7 odd numbers ≤ 13):

| config | even fold | M(fold) at n=7 | M14 |
|--------|-----------|----------------|-----|
| AP `{1..13}` | `{1,2,3,4,5,6}` = the n=7 AP | **1/7 (n=7-tight)** | 1/14 |
| V* `{1..11,13,24}` | `{1,2,3,4,5,12}` | **2/13 (n=7-loose)** | 1/14 |

So the n=14 AP **literally contains the n=7 AP** as its halved even part — the
regular 14-gon's even sub-runners are the regular 7-gon. But tightness at 14 does
**not** require the fold to be n=7-tight (V*'s fold is loose at 2/13): the tightness
is carried by the odd coupling, with the even fold only required to stay `≥1/7` at
the witness. This matches S520's "the menu collapses at n=7": past n=7 the even
fold has slack, and the n=14 extremal structure is governed by the odd interleaving
sitting on top of an n=7-safe even base.

## 5. Honest status

**Not a proof of LRC@14.** What is rigorous: (*), the `min(e,o)≤6` dichotomy, the
LRC(7) even-slack for e≤6, and the antipodal "each odd unsafe at ≤1 preimage".
What is open: (a) prove some even-good `s` avoids the odd-split (closing e≤6 to a
genuine theorem conditional only on LRC(7)); (b) the e≥7 (o≤6) branch, which needs
to dodge ≥7 even runners — LRC(8+) territory.

**Why this is the right lever.** It puts the *proven* n=7 case to work on a full
half of the problem and isolates the entire residual difficulty into one crisp
statement about the odd runners' antipodal split — a much smaller target than
LRC(14) itself, and exactly the place where `2·7` structure lives.

**Handoff:** prove no-odd-split over `E_good` for e≤6 (try: `|E_good|` lower bound
vs. odd-split measure; or use that odds are the residue-odd class mod 14). For
e≥7, seek an odd-side fold or a Davenport-style covering of the ≥7 even runners.

**Artifacts:** `04-computation/lrc_n14_even_fold_to_n7_s554.py` (+`.out`).
Builds on: proven LRC(7); S520 (n=7 menu collapse); HYP-2055 / V* (S553);
HYP-2052 (S551). New: **HYP-2056**.
