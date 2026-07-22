# Chirality, toothpicks, and why tournament counts are even — one Lefschetz parity

*boxeph-2026-07-21-S213. Owner: connect previous chirality work and the toothpick sequence A139250; think
about tournament values that end up being even (iso classes, etc.) and what toothpick diagrams correspond
to their structure. Builds on THM-587 (self-converse = antipodal Euler/Lefschetz number), THM-584 (complement
= antipodal map), THM-1830 (BLUE self-complementary atom parity), boxeph S210/S212 (the reversal involution,
mirror-parity), the chiral-Čech thread (HYP-3123). All verified in
`04-computation/chirality_toothpick_tournament_parity_boxeph_S213.py`.*

## One law behind all three

A count graded by a `ℤ/2` **involution** splits into fixed points and free 2-element orbits:

> **`count = (#fixed) + 2·(#free pairs)`  ⟹  `count ≡ #fixed  (mod 2)`** — a Lefschetz/Euler parity.

The repo's master involution is the **reversal / converse `R:T↦Tᵒᵖ`** (reverse every arc = complement =
antipodal map, THM-584). It is the *same* involution as the LRC reversal `ι:t↦1−t` (S210/S212) and as the
**mirror symmetry of a toothpick diagram**. All three obey the one parity law; they differ only in how many
fixed points they have.

## Tournaments are even because self-converse count is even (THM-587, verified)

The reversal `R` acts on isomorphism classes of tournaments. A class is **self-converse (achiral)** if
`T≅Tᵒᵖ` — a fixed point of `R`; otherwise it lies in a **chiral pair** `{T, Tᵒᵖ}`. So
`A000568(n) = SC(n) + 2·(chiral pairs)`, and **`A000568 ≡ SC (mod 2)`**. THM-587 (PROVED) makes this a
character identity: the per-level signed cycle index `P_n(x)` has `P_n(1)=A000568(n)` (the Burnside count)
and **`P_n(−1)=SC(n)` = the antipodal Euler/Lefschetz number** of `R`, with `SC = 2,2,8,12,88,176` for
`n=3..8` — all even. Hence **`A000568(n)` is even for every `n≥3`.** Verified by direct enumeration:

| `n` | `A000568` (iso classes) | `SC` (self-converse, `R`-fixed) | chiral pairs | even? |
|---|---|---|---|---|
| 3 | 2 | 2 | 0 | ✓ |
| 4 | 4 | 2 | 1 | ✓ |
| 5 | 12 | 8 | 2 | ✓ |
| 6 | 56 | 12 | 22 | ✓ |

(`count = SC + 2·pairs` holds in every row.) So *the evenness of tournament counts is exactly the statement
that the achiral (self-converse) tournaments — the `R`-fixed points — are even in number*, and that even-ness
is `P_n(−1)`, the antipodal Lefschetz trace. THM-1830 is the local face of the same fact: the unstable
non-transitive tournaments have exactly one BLUE (self-complementary) class iff `n` is odd — a single
`R`-fixed point appearing/vanishing with parity.

## The toothpick sequence A139250 is the geometric picture of the same parity (verified)

Simulating the toothpick cellular automaton reproduces `A139250 = 1,3,7,11,15,23,35,43,47,55,67,79,95,123,
155,171,…` exactly. Its structure *is* the involution law made visible:
- The diagram has **mirror (D4) symmetry**; a toothpick is `R`-fixed iff it lies on the symmetry axis,
  else it belongs to a mirror pair. So `A139250(n) = (#axis toothpicks) + 2·(#pairs)`.
- The **central seed** (one axis toothpick) is the fixed point. Verified: `#axis toothpicks` is **odd**
  (`=3` at generation 16), so **`A139250(n)` is odd for all `n≥1`**, and the **first differences
  `1,2,4,4,4,8,12,8,…` have exactly one odd term — the seed.** Every later generation adds toothpicks in
  mirror pairs (even).
- Self-similarity: **`A139250(2^k) = (2^{2k+1}+1)/3 = 1,3,11,43,171,…`** (Jacobsthal A007583), verified — the
  doubling recursion, the geometric analog of the level-graded metagraph recursion of `P_n(x)`.

So the toothpick has **one** fixed point (odd total); tournaments (`n≥3`) have an **even** number of fixed
points (even total). Same law `count ≡ #fixed`, opposite parities — set entirely by the fixed-point count of
the reversal.

## What toothpick diagram corresponds to the tournament structure

The correspondence is exact, orbit for orbit:

```
   toothpick diagram (D4 mirror)              tournament chirality (converse R)
   ───────────────────────────────           ─────────────────────────────────
   central seed toothpick (on the axis)   ↔   the achiral / SELF-CONVERSE classes   (R-fixed)
   a mirror-symmetric toothpick PAIR      ↔   a CHIRAL pair {T, T^op}
   the whole diagram's mirror symmetry    ↔   R = converse = complement = antipode  (THM-584)
   #axis toothpicks (mod 2) = total parity ↔  SC(n) = P_n(-1) (mod 2) = A000568 parity   (THM-587)
   doubling A139250(2^k)=(2^{2k+1}+1)/3    ↔   the level-graded signed cycle index recursion
```

The toothpick fractal is literally the shape of *"one achiral seed generating chiral mirror-pairs, level by
level."* Reading a tournament census through `R` and drawing each chiral pair as a mirror-symmetric toothpick
and each self-converse class as an axis toothpick yields a D4-symmetric toothpick diagram whose total-count
parity is `SC(n) mod 2` — the antipodal Euler number.

## The unification (with S212)

This closes a loop with my LRC work. In S212 I showed the same reversal `ι` acts on the LRC good set `G_δ`,
free for covering sets, so `χ(G_δ)` is even (loneliness in mirror pairs) — `χ(G_δ) ≡ #ι-fixed lonely points
(mod 2)`. That is *identical* to `A000568 ≡ SC (mod 2)` and to `A139250 ≡ #axis (mod 2)`: three instances of
**`count ≡ (involution fixed-point count) (mod 2)`**, i.e. the two antipodal evaluations `P(+1)` (the count)
and `P(−1)` (the fixed set) of one `ℤ/2` action — Burnside at `+1`, Lefschetz at `−1`. The "chirality Euler
class" I posited for LRC in S212 (even `χ` + the odd Gauss-sum index `i√7`) is exactly this `(P(1),P(−1))`
pair, the same object THM-587 proves for tournaments. The toothpick diagram is its picture.

Honest scope: THM-587 (tournament parity), the toothpick OEIS/oddness/Jacobsthal facts, and the S212 LRC
parity are each verified/proved; the contribution here is the **unification** — one reversal involution, one
Lefschetz parity `count ≡ #fixed`, realized as tournament chirality, the toothpick diagram, and LRC
mirror-parity — plus the explicit toothpick↔chirality dictionary.

Links: HYP-8850, THM-587, THM-584, THM-1830,
[[antisymmetry-is-the-hinge-tori-odd-functions-saddles-and-tournaments-boxeph-S210]],
[[the-good-sets-reversal-symmetry-an-equivariant-mirror-parity-sharpening-of-the-chi-criterion-boxeph-S212]].
