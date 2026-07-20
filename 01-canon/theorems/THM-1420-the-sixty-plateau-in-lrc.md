---
id: THM-1420
title: "THE SIXTY IN LRC IS A PLATEAU, AND IT ENDS EXACTLY AT THE TARGET. (I) 60 does appear in the n=12 AP-uniqueness problem, twice, and both are INITIAL-SEGMENT LCMs: the distance-spectrum numerators of the extremal AP are exactly {1,…,⌊n/2⌋}, giving lcm(1..6) = 60 at n = 13; and the covering-system modulus lcm(2..Q) equals 60 at Q = 5,6. Same attractor THM-1415 quantified, not a new phenomenon. (II) THE PARAMETRISED FORM IS THE USEFUL OUTPUT: the constant attached to LRC(n) is lcm(1..⌊n/2⌋), which is 60 for the whole band n = 10,11,12,13 and jumps to 420 at n = 14 — so 60 is a PLATEAU constant that expires exactly at the real target, and the covering moduli make the same jump at the same place. Anything built on '60' does not survive the move to LRC(14). (III) THE ONE SIXTY WITH INDEPENDENT CONTENT PROVABLY DOES NOT ENTER: ord₁₃(2) = 12 says 2 is a primitive root mod 13, but every unit j/n is an equally good witness for the AP optimum — verified exhaustively at n = 11, 13, 17 — so no generator is distinguished and the primitive-root status of 2 has no role in the extremal structure. (IV) AND THE CRUX INTERVAL HAS NO ROOM FOR ONE: 2/25 is exactly the MEDIANT of 1/13 and 1/12, which are Farey neighbours (determinant 1), and the runner-up gap 1/156 has 156 = 12·13 — denominators 12, 13, 25, 156, no 60"
status: >
  (I),(II),(IV) VERIFIED-EXACT — the spectrum identity {min(k,n−k)} = {1,…,⌊n/2⌋} is
  asserted in the script and checked by assertion at every n in range; the mediant and
  determinant identities are exact rational arithmetic.
  (III) PROVED-BY-EXHAUSTION over all units at n = 11, 13, 17: the witness set for the
  AP optimum is the FULL unit group, so it cannot see a generator.
  This is a deflation with a parametrised replacement, not a new result.  It advances
  no open problem; it closes a line of inquiry and supplies the constant that should be
  used instead.
source: kind-pasteur-2026-07-20-S128c111 (owner: see how the three sixties relate to the LRC n=12 AP uniqueness)
depends_on:
  - THM-1415    # the three-sixties deflation and the lcm-attractor statistic
related: [THM-1400]
script: 04-computation/sixty_in_lrc_kps_S128c111.py (+ .out)
---

# THM-1420 — the sixty is a plateau constant, and the plateau ends at n = 14

## I. Where 60 actually is

Two places, both initial-segment lcms.

**The distance spectrum.** For LRC(n) the extremal family is the AP `{1,…,n−1}` with
witness `t = 1/n`, and the distances are `‖k/n‖ = min(k, n−k)/n`. So the numerators are
**exactly** `{1, …, ⌊n/2⌋}` — checked by assertion at every `n` in range — and the natural
constant is `lcm(1..⌊n/2⌋)`. At `n = 13` that is `lcm(1..6) = 60`.

**The covering moduli.** The covering-system reduction tests families mod `lcm(2..Q)`, and
`lcm(2..5) = lcm(2..6) = 60`.

Both are instances of the mechanism THM-1415 measured: `60 = lcm(1..6)` is the smallest
highly-composite target, so lcms of small integers pile onto it. Nothing new is happening.

## II. The plateau, and where it ends

| n | 10 | 11 | 12 | **13** | **14** | 15 | 16 |
|---|---|---|---|---|---|---|---|
| `lcm(1..⌊n/2⌋)` | 60 | 60 | 60 | **60** | **420** | 420 | 840 |
| `lcm(2..Q)`, `Q=⌊n/2⌋` | 60 | 60 | 60 | **60** | **420** | 420 | 840 |

This is the substantive point. **60 is not the constant of the 12-speed problem; it is the
constant of a four-wide band `n = 10..13` that happens to contain it** — which is why it
feels ubiquitous around this material. And it **expires exactly at `n = 14`**, the actual
target, where both the spectrum lcm and the covering modulus jump together to `420`.

So any structure built on the number 60 does not survive the step from LRC(13) to LRC(14).
If a constant of this shape is wanted for the real problem it is `420`, and `lcm(1..⌊n/2⌋)`
is the parametrised form to carry.

## III. The one sixty with independent content does not enter

`ord₁₃(2) = 12` — 2 is a primitive root mod 13, and the number 12 is the speed count. This
is *not* an lcm coincidence, so it is the only one of the three sixties with a chance of
independent content. It fails on a direct test.

The powers of 2 mod 13 are all of `{1,…,12}`, i.e. the AP itself — but so is any full
residue system, so that alone distinguishes nothing. The real question is whether the
**witness set** sees a generator. It does not:

| n | witnesses `t = j/n` attaining the AP optimum | `ord_n(2)` |
|---|---|---|
| 11 | **all 10 units** | 10 |
| 13 | **all 12 units** | 12 |
| 17 | **all 16 units** | 8 (2 not primitive) |

Every unit is an equally good witness at every `n`, including `n = 17` where 2 is *not*
primitive. So the primitive-root status of 2 has no role in the AP's extremal structure,
and `ord₁₃(2) = 12` is a fact about 13 that LRC never consults.

## IV. The crux interval has no room for a 60

The gap is `(1/13, 2/25)`, and

> `mediant(1/13, 1/12) = 2/25` exactly, and `|1·12 − 1·13| = 1` — Farey neighbours.

The runner-up gap in canon is `1/156` with `156 = 12·13`. Denominators throughout: `12, 13,
25, 156`. The interval is pure Farey structure on `1/13` and its neighbour `1/12`; there is
no slot for a 60 and none appears.

## Verdict

60 appears in the 12-speed problem twice, both times as an initial-segment lcm — a fourth
and fifth instance of the cheap mechanism, not a bridge to the Pisano period, the order of
2 mod 1001, or `|A₅|`. The one sixty that is not an lcm provably does not enter. And the
constant that *is* structurally attached to LRC(n) is `lcm(1..⌊n/2⌋)`, which equals 60 only
on the band `n = 10..13` and becomes 420 at the target.

## Named next

- If the covering-modulus route is worked further, the relevant jump is at `Q = 7`
  (`lcm(2..7) = 420`), not at 5 or 6. Any bound tuned to 60 should be re-derived at 420
  before being applied to LRC(14).
- The plateau structure itself is worth one line in the proof map: `lcm(1..⌊n/2⌋)` is
  constant on `n ∈ {2k, 2k+1}` and jumps only at even `n` with a new prime power, so the
  LRC constants move in steps, and `n = 14` is a step.
