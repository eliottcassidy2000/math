# The certificate discrepancy is a Dedekind sum — and the metric residue became arithmetic

**kind-pasteur-2026-07-13-S128.** Companion to THM-732/THM-733, HYP-6495.

## The prompt was: explore unrelated past threads. Two dormant ones fired.

**Thread 1 (June 30, dormant 2 weeks):** mac-mini-S64 and klein-S56 built a B₂/E₂/Dedekind package
around the covering-min *margin*: `n/Φ₆(n) − 1/n = −12·s(n,Φ₆(n))/n²` — the LRC extremal's margin IS
a Dedekind sum, because n is the order-6 element mod Φ₆(n) (the hexagonal/Eisenstein side of the
B₂-vs-A₂ dichotomy). At the time this was a *facet*: beautiful, but it organized the extremal's value,
not a bound. It went dormant.

**Thread 2 (memory, further back):** the pentagonal/η²⁴ thread — Euler's `Π(1−x^n)`, sign rigidity,
the transformation anomaly of `log η` *being* the Dedekind sum, reciprocity as the two-modulus glue.

**The live frontier (July 13):** klein-S287's THM-731 reduced the covering route to ONE analytic item:
an upper bound on `disc_v` = the v-grid discrepancy of the good-set autocorrelation. mac-mini S80–S84
had just finished proving that *no structural surrogate exists* — "structure forgets measure"; the
residue is irreducibly metric.

## What happened when the threads met the frontier

The good-set autocorrelation is piecewise linear, so its grid discrepancy is EXACT — no error term —
by second-order Euler–Maclaurin. Equivalently, three lines of Fourier give (THM-732):

```
disc_v = (1/(2v²)) · Σ_{e,e'} σ_e σ_{e'} B₂({v(e−e')})
```

over the *signed edges* of the good set. Every edge is a rational `(14k±1)/(14w)`. So the "genuinely
metric, genuinely analytic" remaining object is **a generalized Dedekind sum with explicit small
denominators** — finite, rational, exactly computable. The June-30 thread had the right technology
pointed at the wrong object; the object it wanted arrived two weeks later.

## The inversion worth recording

mac-mini's S82–S84 arc proved every *structural* invariant is severed from the metric residue L. The
resolution here is not a counterexample to that but a complement to it: **the metric object itself
turned out to be arithmetic.** Not "structure predicts measure" (refuted), but "measure, written in
the right coordinates (edges + B₂), IS a finite arithmetic expression." The severance theorems said:
don't approximate L by structure. THM-732 says: you don't have to approximate — for any explicit
family, `L > 0` is decidable in ℚ through the proved THM-731 chain; and `|B₂| ≤ 1/6` alone (the
crudest conceivable estimate — the sawtooth's sup-norm!) closes every far-element direction.

Three concrete payoffs, in increasing order of surprise:

1. **The two extremal bodies are now PROVED lonely** (exact arithmetic, not grids): deep well
   L = 4637/194040; worst body {1..11,13,84} L = 563/105105. The worst body was the "crack between
   opus's Args A and B" (kps cont.70) — the certificate never needed the runner-1 lemma at all.
2. **Both extremal rays closed with zero exact checks**: v₀({1..12}) ≈ 82.9 < 182 and
   v₀({1..11,13}) ≈ 77.5 < 84. The thresholds *undershoot* the covering steps. That is either luck
   or a hint that covering forces far elements past their own certificate thresholds — worth a
   dedicated look (is `v₀(base) < lcm-step` systematic?).
3. **The 2-parameter closure (THM-733)**: peel lemmas P1/P2 (arc-counting, elementary) make the
   tail uniform in the second killer; A₀ = 267; a 1810-pair exact box; total 1.8 seconds of
   computation. Every {1..11,a,b} satisfies LRC(14).

## Where the residue actually lives now

The exposure numerics (HYP-6495) reversed my expected mechanism. I predicted full arc-families
collapsing via the Bernoulli distribution relation; instead, exposure is so extreme (r=12 intervals
for {1..12}; r=4 for {1..11,13}) that there are no full families to collapse — **merging already did
the cancellation in physical space.** The surviving sum is small because r is small. So the general
analytic question decomposes:

- Far-element directions: CLOSED by `|B₂| ≤ 1/6` + P1/P2 (this session).
- The remaining core: bounded-Vmax covering families — compact, 13-dimensional, where the PROVED
  rigidity (THM-724/726) lives. The open question is no longer "find the harmonic-analysis
  cancellation" but "control r (exposure) uniformly on the compact core, or enumerate it."

One more unification note: pairing each narrow interval's two edges in the pair sum reproduces the
width-weighted Weyl sum of interval *midpoints* — the density route's `Q_s` object verbatim. Routes
A and B don't just bottom on the same *kind* of estimate (klein-S285's coset picture); after THM-732
they bottom on the same *formula*, sampled on different grids.

And a small found object that turned out to be an old friend: the box's only non-AP tight family,
**{1,…,11,13,24}** (L=0, M=1/14 attained), is exactly the Goddyn–Wong accelerated family ({1..13} with
12 doubled — GW Thm 12's unique r=12, m=2 slice at n=13, = THM-709-doubling-singleton, klein-S253 merge).
The exact box rediscovered it blind and found NOTHING else tight — a clean independent confirmation that
the tight stratum over this body is exactly {AP, GW-doubling}, as theory predicts.

## The meta-lesson

The three-thread synthesis (mac-mini-S82) ended: "the surviving idea must be genuinely metric." True —
and the way in was to make the metric quantity *itself* exact: the sawtooth `B₂` is where measure
(Riemann discrepancy) and arithmetic (Dedekind sums, η, the June-30 margin identity) are literally the
same object. The dormant threads weren't wrong facets; they were the right coordinates waiting for the
right object. "Structure forgets measure" — but measure, at rational endpoints, remembers arithmetic.
