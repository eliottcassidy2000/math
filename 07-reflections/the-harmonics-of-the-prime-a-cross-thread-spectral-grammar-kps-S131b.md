# The harmonics of the prime: a cross-thread spectral grammar

**kind-pasteur-2026-07-26-S131b.** Provenance note, not truth source.

The owner asked to hear how the harmonics of the prime show up in our
number-theory questions. After today's work the answer has a precise
shape, and it is the same shape in at least five threads: **every
question that has yielded decomposes into a smooth carrier plus
discrete lines indexed by primes, one exactly computable local factor
per prime; the hard open problems are exactly the ones that require
certifying that a specific line does not cancel.**

## 1. The literal instance: twin gaps (THM-2447, today)

The distribution of consecutive twin-rank gaps factors, measurably,
as (continuous availability spectrum) x (discrete prime lines): the
p-th line has amplitude exactly `(p-2)/(p-4)` on `p | g` and
`(p-3)/(p-4)` on `p | 9g^2-1` (local lemma, proved; measured within
2.4% for all twelve amplitude checks, p <= 19). Classically this is
the Hardy-Littlewood singular series, whose Ramanujan expansion
`S(g) = sum_q (mu(q)/phi(q))^2 c_q(g)` literally writes the
arithmetic correction as a Fourier series over the primes' additive
characters -- `c_q` is the q-th harmonic. Two lessons with reach:

- **You cannot see a prime's line through the carrier.** The naive
  class-mean estimate of the p=19 line was off by 23% -- not noise,
  but the availability decay confounded with class placement. After
  detrending on the generic class, every line snapped to its exact
  local value. Any repo thread comparing arithmetic classes of a
  decaying statistic must detrend first or it will "refute" true
  local laws.
- **The owner's "which previous terms combine" is spectral.** The
  favoured partners (30, 42, 60, 180) are the gap values whose prime
  lines constructively interfere -- smoothness-ranked, not
  size-ranked (HYP-9025).

## 2. LRC(14): the harmonics of 7 and 13, and the cost of nonvanishing

The rank-eleven endgame lives entirely on the joint spectrum of the
two primes of `91 = 7 x 13`. The stalk `Z/91` is the CRT product of
the `F_7` and `F_13` character groups; THM-2436's decisive object is
a **mixed unit character mod 91** -- nontrivial on both factors, a
genuinely two-dimensional harmonic -- and its proof mechanism (the
punctured stalk cannot stay vertical) is a statement that the defect
cannot be a pure `F_13`-harmonic with trivial `F_7` component.
klein's fresh MSG-2390 lands the same grammar one step deeper: a
same-parent `C_7 x C_13` cell of mass `>= delta/338` whose **tensor
Fourier coefficient `Jhat_l * Ahat_k` is nonzero for every
`l in F_7` and `k != 0`** -- a certified nonvanishing bi-spectrum.
Meanwhile the analytic lane (THM-2051's Fejer--BV alternative) is
the archimedean carrier: Fejer kernels handle the continuous
spectrum, and the whole residual program is the discrete-line
bookkeeping the carrier cannot absorb. The frontier's live item --
owner-conditioned **no-cancellation** for the noncirculant graft
(THM-2445, codex) -- is in this language exactly a
non-vanishing-of-a-harmonic obligation. That is why it is hard: the
easy direction of spectral arguments is bounding lines from above;
the open residues all need lines bounded away from zero.

## 3. Tournaments: the prime 3 owns the rigid spectrum

Today's THM-2444 refutation put the metagraph in the same grammar.
The pure-blue census's rigid stratum has `H = |Aut|` values
`1, 3, 9, 9, 27` (n=9) and `1, 9, 9, 9` (n=10) -- **powers of three
only**, consistent with THM-643 C1's `3^floor((n-2)/2)` cap on
`H_sym`. The prime 3 enters as the odd cycle `C_3`, the unique
smallest non-transitive atom, and (conjecturally, tower probe in
flight) the rigid classes are its iterated lexicographic
compositions -- the tournament world's 3-adic tower. The lone
nonrigid pure-blue class `(H,|Aut|,tc) = (15,5,3)` is the Paley
circulant `T_5`: the quadratic-residue character mod 5, i.e. the
prime 5's order-2 harmonic, appearing at every odd n and refusing
even n. Paley tournaments at `p = 3 mod 4` were already the repo's
canonical spectral extremizers (THM-256, Paley cycle advantage); the
pure-blue inventory now shows the same two primes (3 by towers, 5 by
characters) partitioning a census we had mistaken for an interleave
formula. The involution-parity reflection (S131) is the `p = 2`
instance of the same statement.

## 4. The older threads, re-read spectrally

- **Sum-free/Erdos 592:** the classical extremal dichotomy is
  carrier-vs-line: interval-type sum-free sets (archimedean/middle
  segment) against congruence-type (odd numbers = the 2-line;
  `+-2 mod 5` = the 5-line). THM-469's dichotomy sits on exactly
  that fork.
- **Pentagonal/eta (THM-487/488/489):** `eta^24`'s coefficients are
  the Hecke-coherent package of ALL prime harmonics at once; the
  Euler-sign rigidity thread is about how much of that coherence
  survives in a combinatorial code -- the `[72,36,16]` obstruction
  being combinatorial says the modular harmonics do not by
  themselves force the code.
- **Rosetta/Faulhaber:** the hypotenuse rungs `2, 3, 5, 17` are the
  2-adic harmonic tower (Fermat primes `2^{2^k}+1`), and the Moser
  break at row 6 is where the polynomial carrier stops absorbing the
  combinatorial correction.

## 5. The working law (candidate META-PATTERN, not yet carded)

**Split every arithmetic statistic into carrier x prime lines; prove
local line amplitudes exactly; detrend before comparing classes; and
budget the real difficulty at the no-cancellation steps.** Evidence
from three threads today (twin gaps, LRC tensor cells, rigid
towers); counterindication to record before carding: the black
self-line count at n=8 (THM-849) broke a parity-style prediction --
involution/spectral grammars say nothing when the acting group has
fixed-point-free pieces of unknown size. A second counterindication:
naive class comparison FAILED at p >= 13 until detrended -- the
grammar misleads exactly when the carrier is left in place.
