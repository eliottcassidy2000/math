# Exact Collatz transfers from tournament recursion

Session106, 2026-10-29. **PROVED scoped elementary results / independently
checked finite computations / OPEN universal convergence.** No literature-wide
novelty, Collatz solution, or new large verification record is claimed.

The recent tournament work yields useful Collatz tools when we transfer its
**representation rules**, rather than its numerical spectrum. Extend on the
exact carrier, retain the original source and native clocks, then compress only
what later tests cannot distinguish.

## Main results

1. **Past source-payment tests compress to two bounds.** In an exact valuation
   cylinder write the immutable source as n=r+Mq and the endpoint as y=t+2Aq.
   Exact extension is a congruence split of q. Every conjunction of affine
   source/endpoint guards leaves one integer interval, so the full history of
   those inequalities is represented by its two extrema. The arithmetic lift
   is inherited; interval-filtered source survival is the new assembly.
   [Native transfer](recurrence/SOURCE_INTERVAL_TRANSFER.md) and
   [general guarded macros](recurrence/GUARDED_MACRO_INTERVALS.md).
   For native first-descent tests this sharpens to a single upper cutoff,
   the minimum of the exact prefix carry thresholds.
   [Upper-cutoff theorem and CST boundary](recurrence/NATIVE_UPPER_CUTOFF.md).

2. **A narrow reverse forest, with clocks resolving its ambiguity.** For the
   five strictly growing ternary returns, restricted to shifted coordinates
   m=(n+1)/2 in 81Z, any endpoint has at most **three** ancestors at a fixed
   completed-return depth. The bound is sharp. Those ancestors have different
   native clock pairs (odd steps, divisions), so endpoint plus both clocks
   identifies at most one source and word. The entire predecessor forest needs
   O(log m) inverse-node checks, not exponential word enumeration. This is
   scoped to this guarded language and stable layer. It does not bound a
   forward orbit or repay an original source.
   [Proof and exact decoder](reverse/REVERSE_FOREST.md).

3. **Actual infinite-family descent certificates.** A replayable atlas of 751
   exact words proves descent below the original source on infinite arithmetic
   progressions covering 29,864 of 32,768 odd residue classes modulo65,536.
   The remaining 2,904 classes are unresolved at that precision, not purported
   exceptions. The old guarded compiler is credited; this is its independently
   checked finite realization. The nonnegative residue-polynomial rank route
   is ruled out for active native steps, while unbounded carry/clock ranks
   remain possible.
   [Atlas and state barriers](descent/SOURCE_PAID_ATLAS_AND_STATE_BARRIER.md).

4. **The right spectral target is strong stability, not a uniform gap.**
   The killed native Collatz operator on ordinary ℓ¹ has norm 1 at every power
   and its complex spectrum is the entire closed unit disk, regardless of whether all orbits terminate. An explicitly terminating inverse ray already has that spectrum; its finite cutoffs are all nilpotent. A mixed
   refresh estimate transfers only if it beats the exact no-refresh survival
   cost. A positive discounted incoming-flow supersolution gives an auditable
   source-atomic certificate, but constructing one with the needed mass bound
   remains open.
   [Spectral boundary and Abel certificate](spectral/NATIVE_SPECTRAL_BOUNDARY.md).

5. **Exact future responses still require growing linear memory.** Explicit
   terminating sources exist at every exact hitting clock. Their survival
   response matrix is triangular with determinant 1, forcing dimension at least N
   for all-source linear observation through N clocks. This does not rule out
   a compact nonlinear integer rank, nor a short scalar counting recurrence.
   [Source-response rank](observability/SOURCE_RESPONSE_RANK.md).

## A concrete sharp reverse example

The shifted endpoint 97,551,864 has three two-return ancestors:

| Shifted source m | Growing-return word | Native clocks J,D |
|---:|:---:|:---:|
|43,356,384|00|2,2|
|60,901,146|02|8,12|
|85,545,639|22|14,22|

All are positive multiples of 81 and pass their actual forward guards.
The exact clock relation 5J−3D=2k recovers return count k, while defect D−J
separates the at-most-three histories. This is a rigorous instance where
already-retained arithmetic data recovers joint history information.

## What this says to the convergence effort

The viable next object is a **growing guarded source carrier with short
checkable certificates**, not a fixed finite-state surrogate or a randomized
refresh chain. Interval convexity genuinely erases redundant old inequalities;
clock uniqueness genuinely erases inverse branching. Neither erases unpaid
source debt. The new inverse theorem also rules out distinct same-clock
collision switches staying wholly inside the stable 81 layer and five-growing
language: a successful switch must change clocks or leave that interface.

Three obligations remain separate:

* escaping the language of indefinitely growing auxiliary returns;
* controlling a complete excursion's original-source weighted mass;
* proving its endpoint falls below the original source, or supplying another
  genuinely well-founded terminal certificate.

The recent lower-74 bound concerns only the first language at a finite horizon.
We inspected and credited it, but did not repeat its expensive full search or
promote it to an asymptotic statement. The five-cycle 13→33→83→13 at multiplier 5,
denominator 101→35/175 reclosure, 27-layer collision, and source-mark overshoot
remain explicit hostile controls.

## Verification and provenance

The arithmetic uses Python's standard library, exact integers and fractions,
and explicit exceptions that remain active under optimized Python.

* [Transfer audit](audit_transfer/REPORT.md): independent cylinder rejection,
  interval guards, substochastic kernels, Abel identities, and source-rank checks.
* [Reverse/descent audit](audit_reverse/REPORT.md): literal native queue census,
  independent affine inverse, all 751 symbolic atlas certificates and all 32,768
  residue statuses. No producer implementation is imported by either auditor.
* [Pinned incoming sources](history/source_pins.json), with clearly labelled
  [atomic](history/ATOMIC_INCOMING.md) and [reverse-decoder](history/SEVEN_INCOMING.md)
  historical snapshots.

External context was checked at [Tao's paper, version7](https://arxiv.org/abs/1909.03562v7):
its almost-all logarithmic-density theorem is not an all-source atomic result.
Our notes make no claim to strengthen that theorem.

Run `python3 verify_packet.py` from this directory. It copies the packet to an
isolated temporary directory, checks every declared source hash, runs all eight
checkers in normal and optimized modes, compares receipts, and regenerates the
atlas. See `verification.json` for the exact completed scope.
