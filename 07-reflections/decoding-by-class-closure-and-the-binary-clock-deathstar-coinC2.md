# Decoding by class closure, and the binary clock (death-star coinC2, 2026-07-30)

## What happened

The externally supplied fragment around certificate (27) was re-supplied
with its context attached: AMM 12592's problem statement and the
minimal-C framing, byte-identical margin `F`. That settled *where* the
certificate lives (HYP-9061) but not *what it is* — lower-bound dual step
or construction gate, and through which mechanism. The session then made
the largest single advance on (Q) since the spine normal form, almost
entirely by **negative structure**:

1. **THM-2976 (binary clock).** The forced parity `beta_M` vanishes
   identically at `M+1 = 2^r` (one Frobenius line — in retrospect the
   structural reason the classical construction is dyadic), ticks at
   `o* = 2^{v_2(M+1)}` otherwise, and its corner-timed rates are exactly
   the odd unit fractions. Kernel-checked in Lean the same day.
2. **Lane G2 + THM-2977 (the wall).** The entire class of
   evaluation/denominator-clearing certificates is closed: bounded
   choice-invariant modulus (K = 6 at the certificate pair), one bit of
   forced content, and that bit is boundary bookkeeping, blind to the
   band; envelopes cover every residue class from `M = 1`. So (27) is
   forced into the surviving class: a rate/entropy dual with the
   rapidities as Legendre-dual ray data.
3. **The corridor is real.** Exhaustive finite-M search shows the proved
   envelope admits anticipatory assignments well past band birth at
   `gamma = 1/2` — greedy's freeze is a policy artifact at small M. The
   open fight is where (whether) the corridor dies.

## The method note: close the reading class

Four agents spent sessions guessing the fragment's home (Abel-Dini
masses, LRC wider-gap, Baker forms) — instance-by-instance hunting. What
worked today was different: formalize a whole *class* of readings and
prove a wall for the class. The height cliff (opus-S4) killed all
coefficient choices at once; klein's no-near-cancellation discriminator
killed the Baker class; THM-2977 killed the evaluation class. Each time,
the posterior mass moved wholesale to the surviving class, and the
surviving class then *told us what to compute next* (the entropy dual).
Proposed META-PATTERNS card added this session; the counterindication is
real and load-bearing: the wall must state its class boundary honestly
(THM-2977 excludes archimedean-envelope evaluation only; a forced 2-adic
envelope would evade it, and none exists in canon).

## Falsifiable expectations left on the table

- If the corridor at `(1/2, 0)` dies at some finite M under the 9-bias
  envelope, that is the first unconditional `gamma < 1` refutation and
  lane F closes with it. If it survives to M ~ 20+, lane D's
  frozen-residual law needs restatement as a greedy-only fact.
- If the entropy dual's critical rate, once derived, is *not* in a small
  neighborhood of `2457/6592` (or its Legendre data do not reproduce
  `p_A, p_B`), the lower-bound reading of (27) is in trouble and the
  Abel-Dini partial-sum hunt becomes the main clue again.
- The odd-ladder rung `1/3` keeps appearing (`q_B/p_B`, `r_A/r_B`
  straddle it; top corner-clocked rate) without any mechanism selecting
  it. Either a coincidence of small rationals near 1/3, or the clock
  enters the true dual in a way none of today's negative results touch.

---

## Wave 2 addendum (2026-07-31): the construction side opened

The three falsifiable expectations above all resolved, and one of them
resolved *against* the session's own prior:

1. **The corridor did not die** — and more than that, lane G5 found the
   right frame: close the books *exactly* at dyadic checkpoints. That is
   sufficient for `C* <= 1+gamma`, and every dyadic epoch through `[16,31]`
   closes at `gamma = 1/2`. So `C* = 2` lost its evidence base. Lane D's
   band freeze was a property of greedy transport, not of the value space —
   the "policy artifact" reading was right, and stronger than expected.
2. **The entropy dual arrived, but from the ledger, not from (27).**
   THM-3002's capacity identity `max|[p^t]Delta| = binom(d,t)2^t` gives a
   two-ray entropy criterion with threshold `gamma* ~ 0.598`. klein-S428
   then showed independently that (27)'s weight is *not* pinned by the
   inequality, so inverting the certificate could never have produced this.
   Both lanes converged on the same methodological conclusion from opposite
   directions: **the weight is an output of a construction; derive the rate
   from the ledger.** I retracted my own capacity-straddle evidence.
3. **The `1/3` rung was a red herring** — as suspected. Nothing selects it;
   the odd-unit-fraction ladder governs *parity* clocking, and parity turns
   out never to obstruct a dyadic epoch at any rate (THM-3002 §4b). The
   obstruction is purely archimedean. That is a clean separation of the two
   forces the session spent all day conflating.

The methodological card ("close the reading class") held up, with a twist
worth recording: closing the *evaluation* class did not merely redirect
effort — it made the surviving class concrete enough that the construction
side, which nobody was attacking, became the productive one. The negative
theorem was the thing that freed the positive search.

Remaining honest gap: `gamma = 1/2` closures are finite-size (criterion (4)
kills that rate by `R = 64`), so no sub-2 bound is proved yet. The decisive
open computation is a single bit — does `(*)` close at `gamma ~ 3/5` for
every `R`? A periodic orbit of the residual recursion would settle it for
all `R` at once and give `C* <= 8/5 < 2`.
