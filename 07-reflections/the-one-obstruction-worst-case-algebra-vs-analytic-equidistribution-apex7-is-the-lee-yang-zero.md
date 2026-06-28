# The one obstruction: config-blind worst-case algebra (the apex-7 Lee-Yang zero at λ=1/p) vs config-specific analytic equidistribution

*mac-mini-2026-06-27-S74. The owner asked for a long, more-abstract, synthesizing session. Pulling the
project's recent threads together — my Lee-Yang/Joukowski/Carathéodory/Perron, kps's covariance/associator/
Bonferroni, codex's Asano/equidistribution — they collapse to ONE statement about WHERE the difficulty lives
and WHY every algebraic method stalls at the same place. Integrates [[lee-yang-zeros-explain-bonferroni-failure-kps-S254]],
[[the-moment-problem-merge-caratheodory-toeplitz-griffiths-gks-verblunsky-read-the-lrc-ferromagnetic-phase]],
[[the-joukowski-de-moivre-bridge-lrc-circle-is-the-tournament-real-rooted-class]], HYP-3160/3210.*

## The universal even/odd split (one obstruction, many faces)
Every "hard half" the project has isolated is the **same object**, and every "easy half" is its complement:

| EVEN / easy (provable) | ODD / hard (the obstruction) | level |
|---|---|---|
| κ₂ pairwise covariance | **κ₃ 3-way = the associator** | cumulants (kps S31ai) |
| even Bonferroni terms | **odd inclusion-exclusion terms** | incl-excl (kps S254) |
| Q-block (apex comb, ≤6 events, Lee-Yang) | **R-block (14-free, ≥7 events, interior zero)** | covering (kps S254) |
| cut space (score, abelian) | **cycle space (OCF, H)** | GF(2) (CLAUDE.md) |
| real-rooted / hard-core (repulsive) | **circle / Lee-Yang (ferromagnetic)** | partition fn |
| sum-of-squares (SOS) | **non-SOS (Motzkin/AM-GM-of-3)** | real alg geom |
| commutative / associative | **non-associative (octonion/Fano apex-7)** | algebra |
| degree ≤ 4 (Galois ⊆ S₄, solvable) | **the odd Worpitzky/Eulerian residue** | resolvent (S70/S72) |

The left column is degree-2 / pairwise / abelian / SOS — and it is provable by classical machinery (Perron-
Frobenius, Carathéodory-Toeplitz λ_min, Hermite-Biehler). The right column is **3-way / non-associative / non-
SOS / odd-parity / cycle / apex-7** — one irreducible object wearing eight masks.

## The apex prime 7 IS a Lee-Yang zero at λ = 1/p (NEW, exact)
Model the loneliness as `P(M=0)`, `M = #danger events` (each an arc of measure `p=1/7`). For the
**equidistributed** apex comb, `M ~ Binomial(k, 1/7)`, partition function `G_M(z)=(6/7 + z/7)^k`, whose only
zero is `z = −6`, i.e. in the Lee-Yang variable `λ=1−z`:
```
   the Lee-Yang zero sits EXACTLY at  λ = 7 = 1/p = the apex prime,
   zero-free margin |λ|−1 = 6, growing like k  ->  the certificate.
```
(`lrc_one_obstruction_algebraic_vs_analytic_macmini_S74.py`: min|λ-zero| ≈ 7.000 for all k, 0 inside the
disk.) **The apex 7 is not a metaphor here — it is literally the partition function's zero, at one over the
danger probability.** This is the cleanest appearance yet of the project's 7.

## The meta-theorem (TESTED): algebra is worst-case, the resolution is analytic
| `k` | comonotone `P(M=0)` (worst case) | equidist `P(M=0)` | apex-comb `P(M=0)` |
|---|---|---|---|
| 6 | 1/7 | 0.397 | 0.457 |
| **7** | **0 (saturates)** | 0.340 | 0.335 |
| 10 | 0 | 0.214 | 0.138 |
- **Comonotone (the algebraic certificate's adversary)** saturates at `k=7` (7×1/7 = one full circle) and its
  partition function has **interior Lee-Yang zeros at k=4,5,6** — the Bonferroni failure, the Asano block, the
  non-SOS gap, ALL at once.
- The **equidistributed / actual apex comb stays positive for all k** (zero parked at λ=7).
- Therefore **every config-blind algebraic certificate** (single-variable Lee-Yang zero-freeness, Bonferroni,
  Asano contraction, the moment-LP, SOS) **must fail at k=7**, because it is forced to certify the worst case,
  and the worst case genuinely has no lonely point. The actual config survives **only** by its arithmetic
  (the integer speeds equidistribute) — a **config-specific analytic** fact (Weyl/Erdős-Turán), not an
  algebraic identity.

> **The one obstruction, stated once:** the LRC(14) hard core is the **odd / non-associative / non-SOS /
> coherent (comonotone) worst case**, where the danger partition function pulls a Lee-Yang zero **into** the
> unit disk; no zero-free / SOS / Bonferroni / Asano certificate can cross it; and the **only** thing that
> does is the **equidistribution of the apex comb** (zero pinned at λ=7), which is analysis, not algebra.

## The ferromagnetic twist (why the AP is both extremal and worst-case)
The comonotone config = the **maximally coherent = ferromagnetic = the AP/consec** (S73c/d: consec is the FM
ordered phase, positive covariance, the comonotone end). So the **coverage-maximizer (AP) is the same as the
algebraic worst case** (interior zero). The obstruction *is* the ferromagnetic config — and what rescues it is
that an actual AP has **distinct integer speeds**, whose fine-scale equidistribution (the apex comb's λ=7
clearance) defeats the coherence. **The proof is the fight between coherence (which the AP maximizes) and
arithmetic equidistribution (which the integers enforce) — and arithmetic wins by the margin `1/p − 1 = 6`.**

## Honest status
- **TESTED:** the apex-7 = the equidistributed comb's Lee-Yang zero at λ=1/p (exact); comonotone saturates at
  k=7 with interior zeros at k=4–6; equidistributed/apex-comb stay positive. So the algebraic/analytic split is
  real and located precisely at k=7.
- **SYNTHESIS:** the universal even/odd table — one obstruction (odd/non-associative/non-SOS/coherent/cycle/
  apex-7) with eight faces; the resolution is irreducibly analytic (equidistribution), which is exactly where
  codex's R'≥0.642 spectrum thread and the Erdős-Turán program live.
- **NOT a proof.** This is a meta-diagnosis: it explains WHY the algebraic attacks (Lee-Yang, Asano, SOS,
  moment-LP, Bonferroni) all stall at the same threshold, and certifies that the finishing mathematics must be
  the apex comb's equidistribution. LRC(14) remains open. Value: stop seeking an algebraic certificate for the
  R-block; invest the remaining effort in the effective equidistribution (Erdős-Turán) of the apex comb.

Related: HYP-3221 (this), HYP-3160/3210 (the node), HYP-3202 (moment merge), kps S254 (Bonferroni=Lee-Yang),
codex R'≥0.642 (equidistribution floor), THM-546 (the comb), OPEN-Q-108.
