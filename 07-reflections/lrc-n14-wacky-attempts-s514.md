---
source: oracle-2026-06-01-S514
status: exploratory (wacky proof attempts at LRC n=14; honest partial results)
tags:
  - lonely-runner
  - n14
  - anti-concentration
  - bounded-ansatz
  - parity-split
  - wacky
---

# Wacky Attempts at LRC for n = 14

Three deliberately bold angles at the first open case (13 speeds, threshold
1/14), each *screened computationally* (`lrc_n14_wacky_attempts_s514.py`). None
proves it; one is genuinely proof-shaped. Honest throughout.

## Idea A — anti-concentration (the second moment)

Let `N(t) = #{i : ‖v_i t‖ < 1/14}` (runners crowding the observer). Loneliness =
`N(t)=0`. Each bad set has measure `2/14`, so `E[N] = 13·2/14 = 13/7 ≈ 1.857`. If
the speeds behaved "independently," `N ~ Poisson(13/7)` and
`P(N=0) ≈ e^{-13/7} ≈ 0.156 > 0` — loneliness for free. The conjecture's teeth are
exactly the *resonant* sets where pairwise overlaps make `N` concentrate away
from 0:

```
set                         P(N=0) = lonely measure   (Poisson ref 0.156)
initial 1..13 (tight)       0.00000   ← only exact resonance reaches 0
drop13 add14                0.02390
drop8  add14                0.04519
primes-ish + 14             0.09869
super 1..12, 360360         0.02923
```

So the picture is an **anti-concentration wall**: generic sets sit near the
Poisson value; arithmetic structure pushes `P(N=0)` toward 0, but *only the exact
initial segment touches 0* (and there only on the boundary). A proof would be a
lower bound on `P(N=0)` (or a guaranteed boundary witness) surviving maximal
concentration. The 2nd moment alone won't do it — this is the known reason
density/Fourier methods stall — but it frames the target sharply.

## Idea B — bounded ansatz (the proof-shaped one)

Does every 13-speed set have a lonely time at a *small-denominator* rational
`t = j/(14 s)`? Screened across structured and random sets, including **speeds up
to 1000**:

```
witness cofactor s (Q = 14 s):  always s ∈ {1, 2, 4, 5, 7}
worst s, speeds<40:   4      worst s, speeds<300:  5
worst s, speeds<100:  4      worst s, speeds<1000: 5
```

Every set tested has a lonely witness at `t = j/(14 s)` with `s ≤ 7` — and,
strikingly, **the cofactor does not grow with speed magnitude** (still `≤5` at
speed 1000). The witnesses live at denominators `14·{1,2,4,5,7}` — small
multiples of 14 echoing the `2·7` structure (`{1,2,4}` the 2-part, `{5,7}`
small odds). This is exactly the **bounded-ansatz** shape the recent proofs
(Rosenfeld; Sungkawichai–Trakulthongchai) exploit: combine a minimal-counter-
example speed bound with a finite check over ansatz times `a/(n s)`, `s ≤ s_0`.

**Honest caveat.** Random and simply-structured sets are *easy* — they have
loose witnesses. The genuine hard regime is the adversarial large-speed
structured sets near the counterexample bound, which random sampling cannot
reach. So this is strong evidence that the ansatz cofactor is small, not a proof
that it is bounded adversarially. Still, it is the most concrete proof route: the
open lemma is "the n=14 ansatz cofactor is bounded by a small `s_0`."

> **Wacky conjecture (n=14 ansatz).** Every primitive 13-speed set has a lonely
> time `t = j/(14 s)` for some `s ≤ 7` (cofactor supported on `{1,2,4,5,7}`).

## Idea C — parity split (n = 14 = 2·7) into two proved halves

Split speeds into odd `O` and even `E = 2W`. At time `t`, odd runners see
`‖o t‖`, even runners see `‖(2w) t‖ = ‖w (2t)‖`. Since `|O| + |E| = 13`, the
smaller side has `≤ 6` speeds — a **proved** LRC@≤7 instance (Barajas–Serra). So
both halves *individually* have loneliness witnesses at threshold `1/7 ≥ 1/14`:
there is `τ_O` lonely for `O` and `τ_W` lonely for `W` (at the halved time scale).

The reduction is clean; the obstruction is **coupling**: the odds want a good
`t`, the evens want a good `2t`, and these must be the *same* `t`. So LRC@14 sits
on top of two solved subproblems joined by the doubling map `t ↔ 2t` — the very
2-adic/supersingular seam flagged in the rapidity work (S19). A CRT alignment of
`τ_O` and `τ_W` modulo the 2-part vs 7-part would close it; making that alignment
unconditional is the missing step.

## Verdict

- **A** names the enemy (concentration) but the 2nd moment can't defeat it.
- **B** is the real route and the data is encouraging (small, speed-stable
  cofactor `s ≤ 7`); the open lemma is an adversarial bound on `s`.
- **C** reduces 13 speeds to two *proved* ≤7-speed halves, leaving only the
  doubling-coupling — the cleanest conceptual decomposition.

A through-line: all three localize the difficulty at the **2·7 arithmetic of 14**
(parity split, the `14·{2-part, odd}` ansatz denominators, the resonance that
concentrates `N`). None settles n=14, but B + C together sketch a genuine
program: *bound the ansatz cofactor by CRT-aligning the parity halves.*

## Next
1. Test the ansatz conjecture against adversarial structured large-speed sets
   (not random): does `s ≤ 7` survive sets engineered to be near-tight at large
   scale?
2. Make the parity-split CRT coupling explicit for the case `gcd` structure is
   benign; attempt LRC@14 for sets where the even half is `2·(initial segment)`.
3. Connect the ansatz denominators `14·{1,2,4,5,7}` to the sieve (THM-369) /
   marked-source target (HYP-1987).

## Artifacts
```
04-computation/lrc_n14_wacky_attempts_s514.py
05-knowledge/results/lrc_n14_wacky_attempts_s514.out
```
