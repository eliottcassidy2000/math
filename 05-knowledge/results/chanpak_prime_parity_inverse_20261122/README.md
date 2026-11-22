# Chan–Pak: visible prime parity and absolute inverse tests

**2026-11-22 · PROVED in the scopes below. Chan–Pak OPEN.**
This extends the [positive affine sieve and literal +2 charts](../chanpak_affine_sieve_metagraph_20261122/README.md).
The aim is to keep three distinct predicates separate: a local phase,
a physical history, and an exact ordinary count.

## Main results

- [Small-prime factor parity on the actual words](SMALL-PRIME-PARITY.md).
  The full multiplicities of factors up to a small constant times word
  length give balanced signs. Two strong outputs and their scaled gap
  have joint sign asymptotics. A separate sieve makes both actual
  outputs and the scaled gap bounded-almost-prime simultaneously.
  The gap in this full-state construction is not fixed.
- [Repair the missing 411/413 phases](DYADIC-ATLAS.md).
  One new legal letter preserves the affine invariant and the exact
  backward history decoder. It fills infinitely many **absolute counts
  outside the earlier atlas**, in its two omitted classes modulo 512.
  The expanded middle bank gives infinitely many in every odd class
  modulo $2^{18}$ but has two empty classes modulo $2^{19}$. Separate
  literal one-arc witnesses attain the small **absolute** pair 411/413.
- [The three-fraction equation and word-height bounds](FRUIT-INVERSE.md).
  On the prior orientation quotient, a dyadic guard plus Hensel lifting,
  an Archimedean height bound and rational reconstruction give an exact
  inverse without a Mordell–Weil discrete logarithm. The same
  local-versus-height discipline keeps the physical tournament
  selectors honest; the equation's coordinates are not Hamilton counts.

## Outside sources and limits

[Allahyar–Yianni's 180-gap preprint](https://github.com/sahibzada-allahyar/math/blob/main/preprint/prime-gaps-180/prime_gaps_at_most_180.pdf)
sharpens pair detection by charging the forbidden extreme pair. Its
bound of 180 is conditional; in particular Lean does not discharge
its unestablished Type III hypothesis A4. Our local calculations keep
the joint endpoint forms and isolate their exceptional primes. They
do not import a prime distribution estimate along tournament words.

The [prior fruit quotient](../session97_overnight_arithmetic_carriers_20261029/fractions/README.md)
and [Bremner–Macleod](https://real-j.mtak.hu/22313/1/AnnMathInf_43_2014.pdf)
supply the actual cubic and its geometric context. The new dyadic
test and reconstruction extend it without identifying local feasibility
with existence of an ordinary tournament at a specified count.

The full parity of the factors **above** the controlled cutoff and
an exact selector for every admissible integer are still open.
