---
id: THM-1190
title: THE NATURAL COUNTERCONSTRUCTION AGAINST THE BEAT CERTIFICATE IS SELF-DEFEATING — the certificate fires at a beat q iff (distinct blockers)·k_q < φ(q), and distinct blockers ≤ |{±v mod q}| (verified 0/150 violations, strictly fewer in 53 — so the count is a CONSERVATIVE proxy). Since φ(q)/k_q is small exactly for highly composite q, the natural attack is to force every beat composite: put all speeds in one class r mod m, making differences multiples of m. But then sums are 2r + multiples of m, and making THOSE composite too requires 2r ≡ 0 (mod m), i.e. r = 0 or m/2 — both forcing gcd(V) > 1, a NON-PRIMITIVE family that dilation invariance reduces away (verified: gcd = 210, 105, 120, 60, 60, 30 in the six cases). A targeted sweep of 188 primitive fixed-class families never drove the margin below 1
status: the conservative bound distinct ≤ |{±v mod q}| is verified 0/150 and follows from the S384 pairing lemma; the self-defeat argument is exact arithmetic. The certificate itself remains a CONJECTURE — this rules out one natural construction class, not all counterexamples. A claim in this session that distinct EQUALS |{±v mod q}| was WRONG (53/120 mismatches) and is corrected here
source: opus-2026-07-17-S387 (owner: work whether the beat constraint breaks the blocking construction — second pass)
depends_on: [THM-1175 (the pairing lemma, S384), THM-1170 (beat frequencies), THM-1110 (the blocking construction), THM-1050 (dilation invariance, which does the work here)]
scripts: 04-computation/beat_targeted_opus_S387.py, beat_targeted_corrected_opus_S387.py -> 05-knowledge/results/
---

# THM-1190 — why the certificate resists the obvious attack

THM-1175 proved that at a beat frequency the defining pair collapses to one
blocker, so THM-1110's construction — which assumed **free** residues —
does not apply. It left open whether a cleverer family could still block
all its own beats. This file attacks that directly.

## The conservative count

W_q is symmetric, so ±v ≡ ±v′ (mod q) ⟹ identical kill-sets. Hence

> **#distinct blockers ≤ |{±v mod q : v ∈ V}|**

Verified with **0 violations in 150** trials, and **strictly fewer in 53** —
non-coprime speeds collapse further still. So the ± class count is a *safe,
conservative* proxy: the true blocking capacity is often lower, which only
helps the certificate.

**Correction.** Earlier in this session I claimed these were *equal*. That
was wrong — 53/120 mismatches — and the verification I ran first caught it.
Only the ⟹ direction holds.

## The attack, and why it defeats itself

The certificate fires at q when (blockers)·k_q < φ(q), so it is hardest
where φ(q)/k_q is small — i.e. for **highly composite** q (small φ). The
natural counterconstruction is therefore to force every beat composite.

Put all speeds in a fixed class r mod m. Differences become multiples of m
— highly composite, exactly as wanted. But then **sums are 2r + multiples
of m**, and making those composite too requires

> 2r ≡ 0 (mod m),  i.e. **r = 0 or r = m/2**

and both force a common factor. Verified:

| m | r | gcd(V) |
|---|---|---|
| 210 | 0 | 210 |
| 210 | 105 | 105 |
| 120 | 0 | 120 |
| 120 | 60 | 60 |
| 60 | 0 | 60 |
| 60 | 30 | 30 |

Every case is **non-primitive**, and dilation invariance (THM-1050) reduces
it to a smaller family. The construction cannot be completed while staying
primitive: killing the differences is free, killing the sums as well costs
primitivity.

A targeted sweep of **188 primitive fixed-class families** confirmed it —
the certificate never failed, with lowest margin **1**, falling back to
q = 2 (all-odd, classical sieve).

## Answering the question

**Does the beat constraint break the blocking construction?** Yes, in a
precise sense:

1. THM-1110's construction assumed free residues; at a family's own beats
   the pairing lemma forces collapse, so it does not apply (THM-1175).
2. The natural repair — engineering composite beats — is **self-defeating**,
   because suppressing sums and differences simultaneously forces
   non-primitivity.

**What is still open:** the certificate in general. This rules out one
natural construction class, not all counterexamples, and the certificate
remains strictly stronger than LRC(14) (a union bound is sufficient, not
necessary), so it could fail somewhere LRC(14) holds. I am not claiming it.
