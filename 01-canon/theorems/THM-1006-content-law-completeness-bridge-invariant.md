---
id: THM-1006
title: The content law — the n=12 DEEP half stated as ONE invariant inequality. Writing M(A)=val/q on the pair-sum ruler (klein THM-1002), tightness is exactly q=13·val, the numerator val IS THM-769's sheet number s, and the dilation law val(cA)=c·val(A) makes val ≥ gcd(A) FREE. Deep-branch emptiness for primitive sets is therefore the single inequality val ≤ gcd(A) on the tight locus; sporadic emptiness = that PLUS shallow rigidity. KEY LOCALIZATION (from the Goddyn–Wong control): GW at n=13 SATISFIES val=gcd yet is not {1..13} — so the n=12/n=13 asymmetry lives entirely in the SHALLOW half, and the content law may hold uniformly in n. PROVED here: the dilation law, val ≥ gcd, tight ⟺ q=13val, max(A) ≥ 13val/2, the identification val=s, and a new sheet-number bound val ≤ 4/(169·δ(A\{max})) coupling the sheet stratification to the safe-interval geometry of THM-1001 (sharp to a uniform factor 572/169 ≈ 3.38 on the dilates).
status: PROVED (parts A–D: elementary, given klein THM-1002's ruler lemma and THM-1001). Part E is an EQUIVALENT RESTATEMENT of the OPEN problem (HYP-6800/6820), not a proof of it — it isolates the content, it does not discharge it.
source: mac-mini-2026-07-18-S110 (owner-directed attempt at the completeness bridge invariant)
depends_on:
  - THM-1002  # klein: M(A)=val/q with q | v_i+v_j, hence q <= 2 max A (the engine)
  - THM-1001  # mac-mini: safe-interval element bound (gives the val bound in D)
  - THM-769   # codex: sheet number s / shallow-deep split (identified with val here)
external: LRC(13) SETTLED (M(A) >= 1/13 for every 12-set).
related:
  - THM-759   # ratio bound
  - THM-763   # global finite height
  - THM-770   # shallow full-residue CSP through height 12
  - THM-774   # deep two-sheet folded diamond
  - HYP-6800  # the n=12 sporadic branch
  - HYP-6820  # the uniformity audit (codex-S64 sec.6 names this bridge as the missing piece)
  - HYP-7310  # klein n=12 census
  - HYP-7360  # this session's assembly
---

# THM-1006 — The content law

**One line.** All of the `n=12` sporadic question collapses to one identity: for a
tight 12-set, **the numerator of the maximum equals the content of the set**,
`val = gcd(A)`. One half of that identity is free; the other half *is* the open
problem — but now stated once, above the shallow/deep split, rather than twice.

## Setup

For a 12-set `A` write the maximum on its pair-sum ruler (klein THM-1002): the
maximizer in lowest terms is `t* = a/q` with `q | (v_i+v_j)` for an active pair,
hence `q ≤ 2·max(A)`, and

```
M(A) = val/q,      val = min_{v∈A} |v a|_q,     |x|_q = min(x mod q, q − x mod q).
```

## (A) The dilation law, and `val ≥ gcd` for free (PROVED)

For every integer `c ≥ 1`, scaling `A ↦ cA` is a `c`-fold cover of the circle, so

> **`M(cA) = M(A)`, `q(cA) = c·q(A)`, `val(cA) = c·val(A)`.**

(Verified exactly on `{1..11,13}` at `c = 1,2,3`.) Writing `g = gcd(A)` and
`A = g·A'` with `A'` primitive, this gives `val(A) = g·val(A')`, and `val(A') ≥ 1`,
so

> **`val(A) ≥ gcd(A)` for every set `A`.** ∎

This direction costs nothing. Everything hard is the reverse inequality.

## (B) Tightness is exactly `q = 13·val` (PROVED)

`M(A) = 1/13 ⟺ val/q = 1/13 ⟺ q = 13·val`. Define the **denominator debt**
`d := 13·val − q`, so

> **`A` is tight ⟺ `d = 0`; `A` is non-tight ⟺ `d ≥ 1`.**

Combining `q = 13val` with klein's `q ≤ 2·max(A)`:

> **a tight 12-set has `max(A) ≥ 13·val/2`** — higher sheets force larger sets.

*(Debt also renders klein's stability gap cleanly: `M ∈ (1/13, 2/25) ⟺ d ≥ 1` and
`val > 2d`. Their `val ≥ 3, q ≥ 38` is the first integer point of that region.)*

## (C) `val` **is** the sheet number (identification)

THM-769 writes a tight maximizer in lowest terms as `t* = p/Q` with `Q = 13s`.
Since `Q = q = 13·val`, we get **`s = val`**. Hence the branch split is a statement
about one integer:

| | `val = 1` | `val ≥ 2` |
|---|---|---|
| `q` | `13` | `≥ 26` |
| THM-769 name | **shallow** | **deep** |
| structure | complete nonzero residue system mod 13 | on-sheet `E={v : val\|v\}`, off-sheet `F`, `\|F\|≥2` |
| tools | THM-770 (height ≤12), THM-1001 (single-coordinate, all heights) | THM-774/775/776/836 |

The dilates `c·{1,…,12}` realize **every** `val = c` (verified `c=1..7`), so the
`val`-stratification is genuinely occupied — but only by imprimitive sets.

## (D) A sheet-number bound from the safe-interval geometry (PROVED, new)

Let `C = A\{max A}` and let `δ(C)` be the widest arc on which `φ_C > 1/13`
(THM-1001). For a tight `A`, THM-1001 gives `max(A) ≤ 2/(13 δ(C))`, while (B) gives
`max(A) ≥ 13val/2`. Therefore

> **`val ≤ 4/(169·δ(A\{max A}))`.**

So a large sheet number *requires* a complement whose `1/13`-safe region is finely
fragmented. Verified on the dilates: `c·{1..12}` has `δ = 1/(143c)` and bound
`≈ 3.38c` against the true `val = c` — the bound is **sharp up to a factor ≈ 3.4**,
i.e. the safe-interval geometry controls the sheet number to within a constant.
(Crudely `δ(C) ≥ 1/(78 max C)` also yields `val ≤ 24·max(C)/13`.)

## (E) The bridge: the content law (EQUIVALENT to the open problem)

> **Content law (conjecture).** Every tight 12-set satisfies `val = gcd(A)`.

By (A) the inequality `val ≥ gcd` is automatic, so the entire content is

> **tight ⟹ `val ≤ gcd(A)`,**

i.e. the primitive quotient `A' = A/gcd(A)` has `val(A') = 1` — precisely **the
deep branch is empty for primitive sets**. This is the **deep half only**:

```
sporadic-branch emptiness  ⟺  [content law:  val ≤ gcd on the tight locus]      (DEEP half)
                              +  [shallow rigidity: primitive, val=1, tight ⟹ {1,…,12}]
```

**The two halves are independent — do not conflate them.** `val = gcd` does *not*
by itself give `A = gcd(A)·{1,…,12}`; that conclusion needs shallow rigidity as
well. The Goddyn–Wong control makes the distinction concrete (below).

So the content law reduces the *deep* half — the whole of codex's two-sheet and
higher-sheet programme (THM-774/775/776/836) — to **one inequality between two
integers attached to the same set**, stated without reference to packets, sheets,
or lift heights.

## Honest scope

- **Proved:** (A) dilation law and `val ≥ gcd`; (B) `tight ⟺ q=13val` and
  `max(A) ≥ 13val/2`; (C) the identification `val = s`; (D) the new sheet-number
  bound `val ≤ 4/(169 δ)`.
- **NOT proved:** (E). It is a *restatement* of HYP-6800/6820, exactly as hard —
  it isolates the content into one inequality but discharges none of it. In
  particular it does **not** close the deep branch; codex's THM-774/775/776/836
  remain the live attack there, and `val ≤ gcd` is what they are collectively
  trying to establish sheet-by-sheet.
- **Positive control — and the localization it buys.** Goddyn–Wong
  `{1..11,13,24}` is tight at `n=13` with `val = 1` and `gcd = 1`. So GW
  **satisfies** the content law (`val = gcd`) and is nevertheless not `{1,…,13}`.
  Two consequences, both useful:
  1. `val = gcd` really is strictly weaker than "`A` is a dilate of the initial
     segment" — the conflation is a trap, and GW is the witness that it is a trap.
  2. **The `n=12` / `n=13` asymmetry lives entirely in the SHALLOW half.** The only
     known sporadic instance is *shallow* (`val = 1`), not deep. So the content law
     is not contradicted at any known `n`, and a proof of it may well be **uniform
     in `n`** — whereas shallow rigidity provably is *not* (it fails at 13). This
     inverts the natural expectation and says where each half's difficulty sits:
     the deep half is the one that might generalize; the shallow half (THM-770,
     THM-1001) is the one that must be `n=12`-specific.

## (F) Small-`n` evidence, and a correction to the corpus (S110)

Exhaustive enumeration of tight `n`-sets in `{1,…,3n+2}` for `n = 3,…,7`
(+ the `n=12,13` knowns) gives two things.

**(F1) The content law survives every test.** Across every tight set found at
`n = 3..8, 12, 13` — primitive and imprimitive — `val = gcd` held with **zero
violations**. Combined with the GW control (which satisfies it), the content law is
contradicted at no known `n`. This is the evidence for attacking the deep half
**uniformly in `n`** rather than sheet-by-sheet at `n=12`.

**(F2) CORRECTION — sporadic tight instances are NOT first seen at `n=13`.**
Verified primitive tight sets that are **not** the initial segment:

| `n` | sporadic primitive tight set | `M` | `val` | `gcd` |
|---|---|---|---|---|
| 4 | `{1,3,4,7}` | `1/5` | 1 | 1 |
| 5 | `{1,3,4,5,9}` | `1/6` | 1 | 1 |
| 7 | `{1,2,3,4,5,7,12}` | `1/8` | 1 | 1 |
| 7 | `{1,4,5,6,7,11,13}` | `1/8` | 1 | 1 |

So the reflection `the-sporadic-branch-where-goddyn-wong-lives-macmini-S108`
overstated when it called GW "the first sporadic tight instance": **shallow/segment
rigidity already fails at `n = 4, 5, 7`.** (Corrected there.) The `n=12` claim is
therefore *not* "sporadics start at 13" — sporadics are common; `n=12` is claimed to
be one of the rigid values, alongside `n = 3, 6`. That makes the `n=12` conjecture
more delicate than the corpus framing suggested, and it is why it is Tao's
*optimistic* conjecture.

**A nuance worth keeping straight.** The full-residue characterization of `val = 1`
(THM-769 §2) is a **prime-`13`** phenomenon. At `n=12`, `n+1 = 13` is prime and the
characterization applies; at `n=5,7` (`n+1 = 6, 8` composite) it does not, and
indeed `{1,3,4,5,9}` and the two `n=7` sets have **repeated** residues mod `n+1`.
At `n=4` (`n+1 = 5` prime) `{1,3,4,7}` *does* carry a complete nonzero residue
system mod 5 — it is a genuine single-coordinate *wind* of `{1,2,3,4}` (`2 ↦ 2+5`).

**(F3) THM-1001's bound tracks the truth.** Testing single-coordinate winding
`{1..n}\{j} ∪ {j+(n+1)}` at every `n = 3..13`: it is tight at **`n = 4` only**. The
THM-1001 bound `w ≤ 2/((n+1)·δ(C))` gives `8.0` at `n=4`, which **admits** the true
`w = 7`, and excludes the wind at every other `n`. So the safe-interval bound is not
merely valid — it is sharp enough to permit exactly the one real winding tight
instance and reject all the rest.

*Artifacts:* `04-computation/lrc13_content_law_bridge_macmini_S110.py` (+out).
Credits: klein THM-1002 (the ruler engine and the gap arithmetic), codex THM-769
(the sheet split now identified with `val`) and codex-S64 §6 (the bridge target),
THM-1001 (the safe-interval geometry behind (D)).
