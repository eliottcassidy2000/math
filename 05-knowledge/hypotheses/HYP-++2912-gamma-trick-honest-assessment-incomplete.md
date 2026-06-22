---
id: HYP-+2912
title: HONEST assessment of the gamma-trick for LRC(14) -- a real partial advance (closes |A|>=7 coprime + r<=6), but does NOT complete the proof; the bulk covering case + the prime-tower descent gaps remain
status: ASSESSMENT (verified); the gamma-trick is incomplete -- LRC(14) remains open
source: mac-mini-2026-06-22-S52 (owner asked to complete LRC(14) via the gamma-trick; disciplined honest result)
related:
  - the-gamma-trick-closes-the-14-covering-residual-by-apex-periodicity  # kps S31ad
  - THM-568   # covering => M>1/14 reduction
  - HYP-2909  # the (star) crux + census
---

# HYP-+2912: the gamma-trick is a real advance but does NOT complete LRC(14)

kps's gamma-trick (S31ad) is genuinely clever: a multiple of the modulus is apex-periodic
(||14m(t+1/14)|| = ||14m t||), so the multiples of 14 decouple to a fine phase and the rest is
pigeonholed on 14 points. It CLOSES real cases. But verifying its actual reach:

## What it proves (verified)
- **r = |multiples of 14| >= 7 with R coprime to 14:** the 14-point pigeonhole (each coprime runner
  marks <=2 of 14 bad points, |R| <= 6 => <=12 < 14) gives a good point. REAL, correct.
- **r <= 6:** the union bound (S31v) -- with the condition meas(safe(R)) > r/7.

## What it does NOT prove (verified, the residual)
- The 14-level pigeonhole bound is 2|Codd| + 4|Ceven| + 7|B| (by gcd(s,14): coprime <=2, even <=4,
  odd-mult-of-7 <=7). For 4000 RANDOM covering 13-sets, this bound was < 14 for **0 of 566** -- the
  simple pigeonhole closes essentially NONE of the bulk (random covering sets have few multiples of 14,
  so |R| is large and the bound exceeds 14).
- The residual (R not coprime to 14: even runners, multiples of 7) is left to a prime-tower descent
  14 -> 7 -> 2 -> 1. But the descent has a COUNTING GAP: at the p-point level (p = 7 or 2), the
  coprime-to-p runners mark <=2 each, and there are up to ~13 of them -- far exceeding p points
  (2*13 >> 7 or 2). The pigeonhole does not close at the small-prime levels.

## The covering case IS true, just unproven by the trick
Searching 1059 covering sets, **min M = 1/11 ~ 0.0909 > 1/14 ~ 0.0714** (margin 0.0195): covering sets
are lonely with COMFORTABLE room (0 counterexamples). So THM-568 (covering => M>1/14) holds, but proving
it UNIFORMLY over the infinitely many (unbounded) covering sets still needs either the full descent
(gaps above) or effective equidistribution (Node 3, the analytic core) -- neither complete.

## Honest conclusion
The gamma-trick closes important special cases and is the right *mechanism* (apex-periodicity decoupling).
It does NOT complete LRC(14): the bulk covering case is comfortable but unproven by the trick's
pigeonhole, and the prime-tower descent has a genuine small-prime counting gap. LRC(14) -- open for 13
runners in the literature -- remains open. I will not record a completion that the verification does not
support.