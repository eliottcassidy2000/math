---
id: THM-1050
title: NO SPEED THRESHOLD EXISTS — BONF5 is DILATION-INVARIANT (every S_k satisfies μ(∩ D_{kvᵢ}) = μ(∩ D_{vᵢ}); verified 0/120 failures), so any BONF5-failing family dilates to arbitrarily large speeds while still failing; this REFUTES the S362 hypothesis that the residual might be bounded by a min-speed floor V₀ and hence closable by finite census — no such V₀ can exist, and the small-speed signature observed in S362 is a property of the sampling, not a bound; the surviving refined question is whether the PRIMITIVE failures are finite in number, and the S362 failures are all verified primitive (gcd 1), so they are genuine residual members
status: dilation invariance verified 0/120; the refutation is a two-line consequence; the S362 failures checked primitive; the refined primitive-finiteness question is OPEN
source: opus-2026-07-17-S363 (owner: work the last open regime)
depends_on: [THM-1045 (the S362 taxonomy and its now-refuted framing), THM-1025 (which uses the same dilation invariance constructively), THM-1026 (the ledger)]
scripts: 04-computation/residual_bounded_opus_S363.py -> 05-knowledge/results/residual_bounded_opus_S363.out
---

# THM-1050 — no speed threshold, and what that leaves

**Dilation invariance.** Every S_k is dilation-invariant: the set
∩ D_{k·vᵢ} is the k-fold compressed-and-repeated copy of ∩ D_{vᵢ}, so
their measures over a unit window agree. Verified on 120 random
(subset, dilation) pairs, zero failures. Hence

> **BONF5(k·V) = BONF5(V)** for every dilation factor k.

**The refutation.** In S362 I observed that BONF5 failures concentrate at
small speeds and proposed: *if BONF5 is uniformly positive above some
min-speed floor V₀, the residual is bounded, hence a finite census.*
Dilation invariance kills this outright. Take any failing family V with
minimum speed m; for any candidate floor V₀ choose k with km > V₀. Then
k·V has minimum speed km > V₀ and **BONF5(k·V) = BONF5(V) ≤ 0**. So no
threshold V₀ exists, and the residual is **unbounded as a set of
families**. The small-speed signature in S362 was a property of how I
sampled (random families in [m, 13m] are typically primitive), not a
bound.

**What survives.** The residual is closed under dilation, so the honest
question is whether it is bounded **up to dilation** — i.e. whether only
finitely many PRIMITIVE families fail BONF5. That reformulation is not
vacuous: it is exactly the structure THM-1025 exploits, where reduction
to primitive pairs left a 15-element finite table. The S362 failures are
all verified primitive (gcd = 1 for every one of
{32,60,95,144}, {49,55,59,72}, {63,88,107,186}, {71,72,136,158},
{84,91,167,200}, {89,118,219,322}), so they are genuine residual members
rather than artefacts of a common factor — the reformulated question has
real content and is OPEN.

**Method note.** The invariance that makes THM-1025 work — reduce to
primitive, handle a finite table — is the same invariance that forbids a
speed threshold. A tool and an obstruction with one source; worth
remembering before proposing any bound phrased in absolute speed.
