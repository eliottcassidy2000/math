---
source: opus-2026-07-07-S143
status: two theorems PROVED (|S|=2 linearization; n=12 parity-forced mirror breaking) +
  the mu=M bleeding-edge hunt (companion .out); R2 positioning
tags:
  - lonely-runner
  - circular-chromatic
  - linearization
  - crossing-parity
  - mirror-breaking
  - motzkin
---

# Two small theorems: |S|=2 linearization, and the n=12 parity force

**opus-2026-07-07-S143.** Owner: work the queued handoffs + the LRC bleeding edge. Three of
my S142 handoffs were claimed by the fleet within the hour (mac-mini-S50: Q-vs-H
**anti-correlation** with min Q = Z(n), the n=8 mirror-break mechanism, the two-circle seed;
klein-S168: the c3-parity law; kps-S70: **(A′) reduced to R2** — spread-AP global
PA₂-minimality — the current bleeding edge). This session lands the two that were still
open, plus the μ = M hunt.

## Theorem 1 (|S| = 2 linearization: the bottom rung of GRAPH-LRC is closed)

> For coprime `a < b`: `χ_c(Cay(ℤ, ±{a,b})) = 1/M({a,b})`.

*Proof.* Upper bound: the witness coloring (S141 L1) gives `χ_c ≤ 1/M` always.
Lower bound, two cases.
- **a, b both odd:** every edge changes vertex parity, so the graph is bipartite and
  `χ_c = 2`; and `M = 1/2` (t = 1/2 puts both odd speeds at distance 1/2). Equality.
- **a + b odd:** `M({a,b}) = ⌊(a+b)/2⌋/(a+b)` (verified as part of the run; the witness is
  `t = ⌊(a+b)/2⌋/(a+b)`-adjacent, classical two-distance form), so `1/M = 2(a+b)/(a+b−1)`.
  The explicit closed walk `x, x+a, x+2a, …, x+ba, x+ba−b, …, x` (b steps of +a, then a
  steps of −b) is a **simple cycle of odd length a+b** (distinctness: `ia ≡ (a−j)b` forces
  `b | i` by coprimality — verified for all coprime pairs ≤ 40). An odd cycle `C_{2k+1}`
  forces `χ_c ≥ (2k+1)/k = (a+b)/((a+b−1)/2) = 1/M`. Equality. ∎

So on two-generator distance graphs the lonely-runner quantity is *exactly* the circular
chromatic number, certified by a single odd cycle. The generalization pressure is now
precise: for |S| = 13 the pairwise odd cycles only give `χ_c ≥ 2(a+b)/(a+b−1) → 2`, and
cliques in `G_S` cap near 6 — **no local (bounded-radius) subgraph can certify χ_c ≥ 14**,
so GRAPH-14's content is genuinely global/fractional. The linearization gap is a
long-range-order question, which is exactly the shape of the density-floor rigidity (R2).

## Theorem 2 (the n = 12 mirror breaking is parity-forced)

From the S142 affine parity law: `Q(x) ≡ C(n,4) + Σ_{(l−1)(n−1−l) odd} x_e (mod 2)`.
σ-fixed (grid-symmetric) chords satisfy `a + b = n + 1`, so their length `l = b − a` is
**odd** (n even) — a σ-fixed chord is never parity-carrying — and the parity-carrying
(even-length) chords come in σ-pairs with equal bits under a gridsym assignment. Hence

> **every grid-symmetric page assignment has `Q ≡ C(n,4) (mod 2)`;** in particular, if
> `Z(n) ≢ C(n,4) (mod 2)`, no grid-symmetric assignment attains the 2-page optimum.

The criterion fires at **n = 12** (Z = 150 even, C(12,4) = 495 odd) — the n=12 anomaly of
the S142 table is now a theorem — and at n = 20 in the extended table. Working the
arithmetic gives the full characterization:

> **Theorem 2′.** Parity-forced mirror breaking occurs exactly at **n ≡ 4 (mod 8)**.
> *Proof sketch.* Odd n: no parity-carrying chords, criterion never fires. Even n = 2m:
> `C(n,4)` is odd ⟺ bit 2 of n is set (Lucas) ⟺ n ≡ 4,5,6,7 (mod 8), so even-n oddness
> means n ≡ 4, 6 (mod 8); and `Z(2m) = m(m−1)²(m−2)/4` is even for m ≡ 2 (mod 4)
> (m = 4t+2 gives the factor (4t+2)(4t)/4 = t·(8t+4)/2… even) while for n ≡ 6 (mod 8)
> both Z and C(n,4) are odd. Net: mismatch ⟺ n ≡ 4 (mod 8). ∎ (Verified n ≤ 36:
> fires at 4, 12, 20, 28, 36 only.)

n = 8 (0/16 σ-fixed optima, S142 census) is parity-consistent, so the criterion is
sufficient not necessary; mac-mini-S50's pairing quantization supplies the n=8 mechanism.
Two mechanisms, one phenomenon: the SC/NS mirror breaking of the crossing landscape — now
with an exact arithmetic locus (n ≡ 4 mod 8) for the parity-forced half.

## Engine correction (honest, caught in-session by the verification design)

The |S|=2 closed-form check exposed a bug in **my** `M_exact` (S141/S143): the candidate
denominators scanned the pair sums/differences/doubles *themselves* instead of their
**divisors**, missing e.g. q = 2 — so all-odd sets were under-reported (M(3,7) is 1/2 via
t = 1/2, not 5/14). Fixed with divisor closure; **every previously published value
revalidated unchanged** (the S141 table and the 13-family witnesses contain no all-odd
sets; GW = 1/14 and deep well = 14/183 re-verified). The math of Theorem 1's both-odd
branch was always correct — the mistake was the engine's, and it is the repo's recurring
divisor-closure pattern ("witness denominators DIVIDE pair sums") mis-transcribed; logged
here so the next engine author doesn't repeat it.

## The μ = M hunt (bleeding edge of the collapse conjecture)

Companion run: exact periodic Motzkin (transfer DP, periods ≤ 320) vs exact M over every
primitive 3- and 4-element set with max ≤ 14 — hunting a single μ > M instance, which would
separate the fractional rung of the ladder from LRC. Result in the .out; every set tested so
far sits at μ = M. (The witness-independent-set direction μ ≥ M is L1; the hunt probes the
converse. Cite-checks for Haralambis/Cantor–Gordon still queued for a web session.)

## Positioning the bleeding edge (for the record)

kps-S70: (A′) is reduced to **R2 — the spread inhomogeneous AP is the global
PA₂-minimizer** — which survived exact 1-move local-minimality testing at k=8 but is not
proved; kps classifies it with μ-AP-minimality as the σ-even (measure) core (kps-S67's
grading). Theorem 1's moral applies: R2, like GRAPH-14, is a long-range-order statement
that no local certificate can settle; the tools with a chance are the global ones the fleet
has been sharpening — the Farey-cell exactness (roof machinery on two anchors), the CE/PZ
moment assemblies, and now possibly the LP/duality face of the homomorphism ladder (a dual
feasible solution for the circular-clique LP *is* a global certificate).
