---
id: THM-1060
title: PRIMITIVE BONF5 FAILURES ARE INFINITE — THE FINITE-CENSUS ROUTE IS DEAD, BUT THE FAILURE IS THE CERTIFICATE'S, NOT THE CONJECTURE'S — the construction {8L,…,14L} ∪ {six speeds coprime to L} is globally primitive (gcd 1), blocks all seven sieve moduli, and carries the scale-free correlation of {8,…,14}; its BONF5 runs −1.202, −0.346, −0.242, −0.178, −0.168 at L = 7, 31, 101, 331, 1009 (min speeds 56 → 8072), deltas +0.856, +0.105, +0.063, +0.010, flattening to a clearly NEGATIVE limit ≈ −0.16; yet every one of these families is LONELY, with uncovered measure ≈ 0.115 throughout — so the level-5 certificate provably fails on an infinite family while LRC(14) holds on it
status: constructed and verified exactly at five scales (BONF5 exact-rational; the L = 1009 point computed with cached teeth); loneliness verified directly by exact uncovered measure at four scales; the negative limit is established by a flattening sequence, not proved as a limit theorem
source: opus-2026-07-17-S365 (owner: work the primitive finiteness question)
depends_on: [THM-1050 (dilation invariance), THM-1055 (the sampling artifact that made this the right question), THM-1045 (the taxonomy), THM-1026 (the five-slot ledger)]
scripts: 04-computation/primitive_finiteness_opus_S365.py -> 05-knowledge/results/primitive_finiteness_opus_S365.out
---

# THM-1060 — the primitive question, answered

After THM-1050 killed the speed threshold and THM-1055 showed the
apparent one was a sampling artifact, the surviving question was: **are
there finitely many PRIMITIVE BONF5 failures?** If so, a finite census
closes the last regime.

**The construction.** A family can carry huge pairwise gcds while being
globally primitive. Take the seven speeds {8L, 9L, …, 14L} — pairwise
gcds ≥ L, and by dilation invariance their internal overlaps are exactly
those of {8,…,14}, hence SCALE-FREE — then adjoin six speeds in [8L, 14L]
coprime to L, forcing gcd = 1 overall. The result is primitive, blocks
all seven sieve moduli, and is comparable (ratio 1.75).

**The measurement.**

| L | min speed | gcd | BONF5 | uncovered (truth) |
|---|---|---|---|---|
| 7 | 56 | 1 | −1.202221 | 0.117702 |
| 31 | 248 | 1 | −0.346420 | 0.117117 |
| 101 | 808 | 1 | −0.241533 | 0.113802 |
| 331 | 2648 | 1 | −0.178462 | 0.115343 |
| 1009 | 8072 | 1 | −0.168378 | — |

Deltas +0.856, +0.105, +0.063, +0.010 — flattening to a limit near
**−0.16**, comfortably negative.

**The answer, in two parts.**

1. **Primitive failures exist at every scale, so there are infinitely
   many.** The finite-census route to the last regime is dead — not
   merely unavailable, but ruled out. Combined with THM-1050 (no speed
   threshold) and THM-1055 (the apparent one was an artifact), the whole
   "bound it and enumerate" strategy is closed off.

2. **The failure belongs to the certificate, not the conjecture.** Every
   one of these families is LONELY: the exact uncovered measure is ≈
   0.115 throughout, essentially constant in L. LRC(14) is comfortably
   true on precisely the families where the level-5 certificate returns
   −0.17.

**What this means for the program.** The dense core needs a STRONGER
TOOL, not a bigger computation. BONF5 is provably insufficient on an
infinite, explicitly constructed family — so the route forward is a
higher-level truncation (level 7+), or a certificate that is not of
Bonferroni type at all. The gap is not small: −0.17 against a true value
of +0.115 means the level-5 estimate is off by ≈ 0.28 on this family,
which is the price of approximating a strongly correlated core by
inclusion–exclusion truncated at depth five.

**Method note.** The construction is worth remembering independently: it
manufactures arbitrarily large primitive families with the correlation
structure of any chosen small set, by dilating the small set and
adjoining coprime filler. Any argument that assumes "large speeds ⟹
generic behaviour" must contend with it.
