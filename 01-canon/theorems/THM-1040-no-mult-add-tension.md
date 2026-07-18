---
id: THM-1040
title: THE MULTIPLICATIVE/ADDITIVE OVER-DETERMINATION DOES NOT EXIST — a REFUTATION of my own S360 proposal: blocking all seven sieve moduli {8,…,14} and having minimal additive doubling are INDEPENDENT, not conflicting (median doubling 6.92 among blockers vs 6.92 among non-blockers — identical), and explicit arithmetic progressions with doubling 1.92 block all seven (a = 665 d = 53; a = 1927 d = 137; a = 259 d = 49), so the two constraints are simultaneously satisfiable and no over-determination is available; the genuine positive salvaged is that the seven-moduli sieve alone kills ~89% of random comparable 13-families
status: refuted (my own S360 conjecture), verified on 6000 comparable families and 20000 AP trials; the 89% sieve kill rate is measured on 4000 random comparable families
source: opus-2026-07-17-S361 (owner: work the multiplicative/additive over-determination tension)
depends_on: [THM-1035 (the seven-moduli sieve, kernel-pure — the multiplicative side), THM-1026 (the ledger), the repo's Freiman/AP modules (the additive side)]
scripts: 04-computation/mult_add_tension_opus_S361.py -> 05-knowledge/results/mult_add_tension_opus_S361.out
---

# THM-1040 — the tension does not exist

In S360 I proposed that the dense core is over-determined: it must satisfy
a purely MULTIPLICATIVE constraint (contain multiples of all seven moduli
8…14, else the sieve produces a lonely time) and simultaneously carry
ADDITIVE structure (small doubling, which is what makes BONF5 fail). I
suggested that tension was "where I would look next." **It is not there.**

**The two constraints are statistically independent.** Over 6000 random
comparable 13-families:

| class | n | median doubling \|A+A\|/\|A\| |
|---|---|---|
| blocks all seven moduli | 732 | 6.92 |
| does not | 5268 | 6.92 |

Identical. Blocking exerts no pressure on doubling whatsoever.

**And they are simultaneously satisfiable at the additive extreme.**
Arithmetic progressions — doubling 1.92, the minimum — block all seven
moduli 30.6% of the time. Explicit witnesses: (a,d) = (665,53),
(1927,137), (259,49). The reason is structural rather than accidental:
for an AP a + kd with gcd(d,q) = 1, the residues a + kd run over all of
ℤ/q as k ranges over 13 consecutive values whenever q ≤ 13, so *every*
such modulus is automatically hit; only q = 14 can be missed, and only
for one residue class of k.

**Why the idea was wrong.** I was treating "multiplicative" and
"additive" as competing demands on a scarce resource. They are not: an AP
is a translate of a dilate, and dilation is exactly the operation that
distributes multiples of every small modulus through the family. Additive
structure *helps* satisfy the multiplicative constraint rather than
fighting it. The Freiman analogy (arXiv:2408.00183) suggested a smallness
hypothesis forcing rigidity — but here the smallness hypothesis and the
divisibility demand point the same way.

**What survives, and it is worth having.** The seven-moduli sieve is far
stronger on comparable families than I had appreciated: only **10.8%** of
4000 random comparable 13-families block all seven, so the kernel-pure
sieve of THM-1035 disposes of roughly **89%** of the comparable regime
outright, before any measure-theoretic machinery is invoked. The dense
core is the ~11% that survive, and among those the additively structured
ones (APs) are already covered by the corpus's AP theorems.
