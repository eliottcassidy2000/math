---
id: THM-1030
title: THE KILLER-COLLAPSE IS THE CLASSICAL SIEVE — the owner's modulus collapse (q = ⌈v₁/m⌉ makes every residue small, so a clustered 13-family behaves at modulus q like a family of speeds ≤ 12) is REAL and vivid: residues drop from ~10⁴ to a median of 292; but the ratio-13 band condition e_max ≤ 13·e_min carries almost no information for q ≤ 27 (residues automatically lie in [1, ⌊q/2⌋] ⊆ [1,13]), and it is NOT sufficient — when the residue pool has ≤ 13 values, 13 speeds must realise ALL of them by pigeonhole, and such a family has NO lonely p at that modulus (checked q = 15,18,21,24,27); the sharp form of the mechanism is exactly the classical sieve at q ≤ 14 with witness t = 1/q, which the repo already holds kernel-pure
status: verified exactly; includes a SELF-CORRECTION — an earlier run of this script reported "100%" for adaptive-q, but it measured only the band condition and not loneliness; the corrected figure is 111/120 = 92%, with the 8% failures being precisely the pigeonhole-forced families
source: opus-2026-07-17-S359 (owner: consider the killer-collapse idea)
depends_on: [LonelyRunner.lonely_of_no_multiple + counterexample_needs_all_divisors (the classical sieve, already kernel-pure), THM-1026 (the five-slot ledger this does NOT fill)]
scripts: 04-computation/modulus_collapse_opus_S359.py -> 05-knowledge/results/modulus_collapse_opus_S359.out
---

# THM-1030 — the killer-collapse, and what it is

**The collapse is real.** With q = ⌈v₁/m⌉ we get v₁ = qm − s, 0 ≤ s < m,
so v₁ ≡ −s (mod q); for a cluster v_i = v₁ + d_i every residue is
d_i − s, small. Measured on 267 clustered 13-families: residues fall from
speeds of order 10⁴ to a median e_max of **292**. The picture the owner
described — "the family behaves at modulus q like a family of speeds
≤ 12" — is accurate.

**The band condition is nearly vacuous, and not sufficient.** For q ≤ 27
the residues (distance to the nearest multiple) lie in [1, ⌊q/2⌋] ⊆
[1, 13], so e_max ≤ 13·e_min holds almost automatically. Worse, when the
residue pool has ≤ 13 values, 13 speeds must realise **all** of them by
pigeonhole — and a family realising every residue has no lonely p at all:

| q | residues present | lonely p |
|---|---|---|
| 15 | 1..7 | none |
| 18 | 1..9 | none |
| 21 | 1..10 | none |
| 24 | 1..12 | none |
| 27 | 1..13 | none |

**Where the mechanism is sharp.** At q ≤ 14 the band requirement
"residue ≥ q/14" reduces to "residue ≥ 1", i.e. simply q ∤ v_i, and the
witness is t = 1/q. That is exactly `LonelyRunner.lonely_of_no_multiple`,
with `counterexample_needs_all_divisors` as its contrapositive: a
counterexample must have, for every 2 ≤ q ≤ 14, some speed divisible by q.
Both are already kernel-pure in the corpus. **The collapse does not extend
the modulus range past 14** — the q > 14 extension fails on exactly the
pigeonhole families above.

**Self-correction (recorded deliberately).** An earlier run of this
script reported that adaptive q works "100% of the time". It measured
only whether some q satisfies the band condition — not whether a lonely
p exists. Corrected, the figure is **111/120 = 92%**, and the 8% failures
are precisely the pigeonhole-forced families. The 100% was measuring the
wrong predicate.

**Value.** The collapse is a genuinely illuminating reformulation: it
explains *why* the classical sieve is the natural first tool for
clustered families (their residues really do collapse), and it identifies
the dense core concretely as the stratum that blocks every small modulus
— "needs all divisors". It does not add a new theorem to the five-slot
ledger (THM-1026).
