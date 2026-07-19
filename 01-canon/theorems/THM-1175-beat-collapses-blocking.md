---
id: THM-1175
title: THE BEAT CONSTRAINT COLLAPSES BLOCKING CAPACITY — PROVED — at a family's own beat frequency q = vᵢ + vⱼ the two speeds kill IDENTICAL numerator sets, because W_q is symmetric (r ∈ W_q ⟺ −r ∈ W_q) and vⱼ ≡ −vᵢ (mod q); the same holds at difference beats q = |vᵢ − vⱼ|, where vⱼ ≡ vᵢ directly. Verified 0/200 and 0/174 discrepancies. So 13 speeds cannot act as 13 independent blockers at their own beats — measured collapse to as few as 4–8 distinct kill-sets — which directly weakens THM-1110's construction, whose blocking assumed FREE residues. The unified beat certificate (classical sieve at q ≤ 14, counting theorem at q > 14) fired on every family tested, including all three tight families at their own q = 14 = 1 + 13, and survived adversarial hill-climbing at margin ≥ 1
status: the pairing/redundancy lemma is PROVED (elementary: symmetry of W_q plus vⱼ ≡ ±vᵢ) and verified exactly. The blocker-collapse figures and the certificate's survival are MEASURED, over 30-family hill-climbs — and this is the same class of evidence that misled me in MISTAKE-152/154/156/157, so the certificate is a CONJECTURE, not a theorem. Note also it is strictly stronger than LRC(14) (a union bound is sufficient, not necessary), so it could be false while LRC(14) holds
source: opus-2026-07-17-S384 (owner: work whether the beat constraint breaks the blocking construction)
depends_on: [THM-1170 (the optimum sits at beat frequencies), THM-1110 (the blocking construction this weakens), THM-1100 (the band condition), MISTAKE-172]
scripts: 04-computation/beat_blocking_opus_S384.py, beat_certificate_unified_opus_S384.py -> 05-knowledge/results/
---

# THM-1175 — why blocking is harder at a beat

## The pairing lemma (proved)

The forbidden window W_q = {r : min(r, q−r)·14 < q} is **symmetric**:
r ∈ W_q ⟺ −r ∈ W_q. At a sum beat q = vᵢ + vⱼ we have vⱼ ≡ −vᵢ (mod q), so
for every numerator p,

> vⱼp ≡ −vᵢp (mod q),  hence  vᵢp ∈ W_q ⟺ vⱼp ∈ W_q.

**The two speeds kill exactly the same numerators.** At a difference beat
q = |vᵢ − vⱼ| the same conclusion holds even more directly, since vⱼ ≡ vᵢ
(mod q). Verified with **0 discrepancies in 200 sum beats and 174 difference
beats**.

## Why this matters against THM-1110

THM-1110 established that 13 speeds can block *any* single modulus, because
each kills k_q numerators and 13·k_q ≥ φ(q) always. That construction
assumed the residues were **free**. At a family's own beat frequency they
are not: the defining pair collapses to one blocker, and in practice far
more coincidences occur — measured **distinct kill-sets as low as 4–8** out
of 13 speeds.

So THM-1110's blocking argument **does not apply verbatim** at beat
frequencies, which is exactly the question THM-1170 raised. The blocking
capacity of a family against its own beats is strictly less than against an
arbitrary modulus.

## The unified beat certificate (conjecture)

For a beat frequency q:
- **q ≤ 14**: W_q = {0}, so a speed blocks only if q | v. If q divides no
  speed, every numerator survives — the classical sieve.
- **q > 14**: if (#distinct kill-sets)·k_q < φ(q), some numerator survives.

Tested on all three tight families and under adversarial hill-climbing:

| family | certifying beat | kind | margin |
|---|---|---|---|
| {1,…,13} | **q = 14 = 1+13** | sieve | 6 |
| {1,…,11,13,24} | **q = 14 = 1+13** | sieve | 6 |
| 2·{1,…,13} | q = 28 | count | 10 |
| {1,3,…,25} | q = 2 | sieve | 1 |
| adversarial worst | q = 2 | sieve | **1** |

The hunt bottomed at margin 1, never at 0 or below.

## Status — deliberately guarded

The pairing lemma is proved. **The certificate is not.** Two cautions I
want on record:

1. The supporting evidence is adversarial hill-climbing, which is precisely
   the evidence type that misled me four times (MISTAKE-152, 154, 156, 157) —
   each time a sampled extremum looked like a bound and was not.
2. The certificate is **strictly stronger than LRC(14)**: a union bound is
   sufficient for loneliness, not necessary. So it could fail on some family
   while LRC(14) still holds there. Proving it would prove LRC(14);
   refuting it would prove nothing.

What is genuinely new and solid is the collapse: **at its own beats, a
family's speeds are provably not independent blockers.** That is the first
structural fact this programme has found that weakens the blocking side
rather than strengthening the covering side.
