---
source: opus-2026-06-06-S702 (poke Steering Task 1 + 2)
status: ANSWERS poke Task 1 — the binding shell-partner q=a+b (HYP-2296) and the prime-torsion leak
  (THM-421/427) are ONE object: the antipodal involution σ:x↦−x. (A) synchronization (THM-425) =
  σ-invariance of ‖·‖; shell-partners = σ-orbits. (B) self-partners = σ-fixed = 2-torsion {0,q/2}
  (the half-turn, n=14 r=7). (C, PROVED) the signed floor is NEVER the half-turn — a half-turn
  relative speed gives ‖·‖∈{0,½}, never the small binding value — so poke's 2-torsion leak is
  structurally EXCLUDED from setting the floor. (D, PROVED) σ_2=id ⟹ the antipodal/shell-partner
  content is an ODD-PRIME phenomenon (mod 2 always trivial self-partner; genuine orbits mod odd p).
  (E) 2n−1 and the block q=4n−5 are always ODD ⟹ shell antipodally free; clock n even ⟹ half-turn
  n/2 leak. THM-430, HYP-2297. Task 2: the realized q's (19,27,42,29,20,8,25) tend to share a prime
  with the CLOCK n (gcd(q,n)>1), not the shell (gcd(q,2n−1)=1) — observed, not proved.
tags: [lonely-runner, LRC, antipodal-involution, sigma, shell-partner, torsion-leak, half-turn,
  2-torsion, synchronization, CRT, odd-prime, signed-floor, max-cut, poke-task, THM-430, HYP-2297]
---

# The antipodal involution σ:x↦−x unifies the shell-partner floor and the torsion leak

**Prompt (poke, via `comms/POKE-COORDINATION.md` Steering Task 1):** "If the minimal gap is realized
at a synchronized shell-partner `a+b=q`, how does this `q`-denominator relate to the torsion subgroup
of the divisor projections mod 2 and mod 7? Investigate if the shell-partner binding is constrained
by the same fiber alignments observed in the n=14 and n=15 leaks."

This is a bridge question between two lanes that the repo had so far kept **complementary**: the
torsion-leak (n-grid, THM-421/427, my S701) and the signed shell-partner floor (mod `q`,
THM-425/429/HYP-2296, the monad lane). THM-428 already showed the two *moduli* `n` and `2n−1` are
CRT-coprime (orthogonal towers). Poke asks about the *third* modulus — the **actual optimal witness
denominator `q = a+b`** — which is neither `n` nor `2n−1` (the realized values are 19, 27, 42, 29,
20, 8, 25). The answer is that all three are the **same antipodal involution** read on different
moduli.

## The one object: σ

`‖x‖ = ‖−x‖`. So every loneliness modulus carries the involution `σ_N : x↦−x` on `ℤ/N`, with
`Fix(σ_N) = {x:2x≡0} = T_2^{(N)} = {0, N/2}` (the 2-torsion). Then:

1. **A shell-partner `{a,b}` with `a+b≡0 mod q` is a σ_q-orbit** (`b≡−a`), and THM-425's
   synchronization `‖a k/q‖=‖b k/q‖` is *literally* the σ-invariance of `‖·‖`. The monad lane's
   "binding pair = synchronized shell-partner" becomes "the floor is an **antipodal orbit**."
2. **The self-partners** (a pair collapsing to one residue, `a≡b`) are the **σ-fixed points = the
   2-torsion** `{0, q/2}`. Poke's n=14 half-turn `r=7` is exactly `Fix(σ_{14})\{0} = {14/2}` — the
   clock's 2-torsion, the same residue THM-427 flags as maximal cell-leak.

## The sharp finding: the floor never sits on the half-turn (THM-430 C)

The minimizer census (5 published + 7 searched) is unanimous: **every binding pair is a genuine
2-orbit, none a self-partner.** And this is provable, not just observed: a half-turn relative speed
`w=q/2` gives `‖w k/q‖=‖k/2‖∈{0,½}` — zero (kills loneliness) or one-half (the maximum). It can
never be the *minimising* speed of a positive floor `M=k/q<½`. So the 2-torsion leak poke points at
is **structurally excluded from setting the signed floor.**

This reconciles a tension. THM-427 C2 says the half-turn (`g=n/2`) is a *maximal* cell-leak — the
"most dangerous" defect in the counting model. Yet it never binds the signed floor. The involution
resolves it: **a σ-fixed point leaks maximally in the count but is a max (not a min) of `‖·‖`, so it
can index danger without ever being the witness's binding constraint.** Fixed points are loud,
orbits are binding.

## The fiber answer to "mod 2 vs mod 7" (THM-430 D)

σ commutes with CRT, so a shell-partner splits into per-prime σ-orbits. But **`σ_2 = id`** (`−1≡1
mod 2`): on the 2-fiber every pair is a trivial self-partner. The genuine antipodal content lives in
the **odd-prime fibers.** Verified on `n=7`'s `{19,23}`, `q=42=2·3·7`: self/fixed mod 2 `(1,1)`,
genuine orbit mod 3 `(1,2)` and mod 7 `(5,2)`. So poke's literal question — relate `q` to the
torsion mod 2 and mod 7 — answers itself: **mod 2 is always the degenerate fiber; the shell-partner
binding is an odd-prime (mod 7, mod 3, …) alignment.** The "fiber alignment" is real, and it is an
odd-prime phenomenon.

## Why the shell modulus is odd, structurally (THM-430 E)

Because the antipodal action is trivial at 2, the *interesting* shell-partner structure wants an odd
modulus — and indeed both canonical shell moduli are forced odd: `2n−1` (THM-401) and the
consecutive-block witness `q=4n−5` (HYP-2296 B). On an odd modulus `σ` has **no** nonzero fixed
point, so the face is *antipodally free*: `(N−1)/2` clean 2-orbits, no half-turn. The clock `n`,
free to be even, is the only face that carries the half-turn. **This is the antipodal explanation of
THM-428's parity asymmetry:** n=14's hardness is the *odd* prime-cube `3³` shell tower, never the `2`
in its clock, because `2` is antipodally inert. The deep lane (S708/S710 homometry at `C=3³`) is
odd-prime by necessity.

## Task 2 (the q-denominators and their relation to n)

Poke Task 2 asked to report the denominators `q` realizing the infima and their relation to `n`.
From the census:

| n | config (a min.) | floor M | binding `{a,b}` | `q=a+b` | factor q | gcd(q,n) | gcd(q,2n−1) |
|---|---|---|---|---|---|---|---|
| 6 | (2,3,4,6,8) | 3/19 | {9,10} | 19 | 19 | 1 | 1 |
| 6 | (4,5,8,10,15) | 4/27 | {4,23} | 27 | 3³ | 3 | 1 |
| 7 | (2,4,7,10,11,12) | 5/42 | {19,23} | 42 | 2·3·7 | 7 | 1 |
| 6 | (5,7,8,9,11) | 3/29 | {11,18} | 29 | 29 | 1 | 1 |

**Observed (not proved):** `gcd(q, 2n−1)=1` in every case, while `gcd(q,n)>1` is common — the
witness denominator aligns with the **clock** primes, not the shell `2n−1`. And `q=4n−5` for the
consecutive-block family is the shell modulus `2N−1` of the **doubled** system `N=2n−2`. Both are
recorded in HYP-2297 for a larger census; neither is yet a theorem. The minimizers all carry an
**irreducible small-speed cluster** (`{2,3,4}`, `{3,4,5}`) — the `r_min≥n` obstruction of THM-429
Cor 2 — whose torsion alignment is the natural next probe (does the forced cluster sit in a low-order
fiber?).

## Honest status

- **PROVED:** σ-orbit recast of synchronization (A); self-partner = 2-torsion (B); floor never a
  self-partner (C); σ_2=id ⟹ odd-prime antipodal content (D); odd-shell/even-clock face split (E).
- **VERIFIED:** 12 minimizers (n=5,6,7; searched + published) all genuine σ-orbits; face table n=5..15.
- **OBSERVED, open (HYP-2297):** `gcd(q,2n−1)=1` & `gcd(q,n)>1` clock-alignment of the witness
  denominator; the small-cluster torsion alignment (Task 2 frontier).
- **Not resolved:** any open LRC case. This is a *structural unification* that localizes and explains
  (it reconciles the half-turn's maximal-leak vs never-binding roles, and the parity of the shell),
  not a new bound.
- **Concurrency credit:** THM-425/426 (monad-S1), THM-427/428/429 + HYP-2296 (monad-S3) are the
  substrate; this adds the antipodal involution as their common root and answers poke Task 1/2.

**Artifacts:** `04-computation/lrc_antipodal_shellpartner_torsion_s702.py` (+`.out`). Theorem:
**THM-430**. New: **HYP-2297**. Logged to `comms/CLUSTER-FEED.md` per poke protocol.
