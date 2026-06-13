---
source: oracle-2026-06-01-S546
status: synthesis + verified computation (odd/even cycles <-> odd/even numbers <-> Goldbach/Lemoine; doubled primes = the resonant diagonal)
tags: [doubled-primes, odd-cycles, burnside, goldbach, lemoine, parity, A000568, lrc, doubling, resonance]
---

# Doubling is pairing: what the doubled primes are, across cycles, numbers, and LRC

**Prompt (user):** understand odd/even cycles, odd/even numbers, and primes
(Goldbach: even = 2 primes; odd = prime + a *doubled* prime), and figure out the
importance of the doubled primes.

The pieces line up into one principle: **doubling (×2) is the operation that
PAIRS, and the doubled primes `2q` are the resonant diagonal where structure
concentrates** — in additive number theory, in the tournament's odd-cycle Burnside,
and in LRC. Here is the dictionary, all verified.

## 1. Parity = the parity of the number of ODD parts

```
even n = p + q        Goldbach          two odd primes        (an EVEN # of odd parts)
odd  n = p + 2q       Lemoine/Levy      prime + DOUBLED prime  (THREE odd parts: p, q, q)
doubled prime 2q = q + q = the DIAGONAL of Goldbach (the even number that is twice a prime)
```

A number's parity is the parity of how many odd primes it takes: even needs an even
count (2), odd needs an odd count (3, with two equal: `p + q + q = p + 2q`). The
**doubled prime is the repeated pair `(q,q)`** that turns the 2-prime even story
into the 3-prime odd story. (Verified: Lemoine reps `p+2q` for odd `n`; the doubled
primes `2q = 4,6,10,14,22,26,...`.)

## 2. Tournaments live on ODD cycles; even cycles vanish

A000568 (tournament iso classes) is the odd-cycle Burnside (verified exactly,
n=3..8):

```
A000568(n) = (1/n!) sum_{ all-odd cycle types lambda } |class| * 2^{ e(lambda) },
   e(lambda) = sum_i (l_i - 1)/2  +  sum_{i<j} gcd(l_i, l_j),
   and Fix(sigma) = 0 the moment sigma has ONE even cycle (an even cycle reverses an
   antipodal edge -> no sigma-invariant orientation).
```

So **an even cycle is a `0` — the additive analogue of "even is not a single prime."
Only odd cycles carry tournament structure**, exactly as only odd primes build the
Goldbach/Lemoine sums. The number of odd cycles plays the role of the number of odd
primes; its parity is the parity of `n`.

## 3. The doubled prime = the equal odd-cycle pair = the maximal cross-term

Here is the importance, made precise. The Burnside exponent's between-cycle term is
`sum_{i<j} gcd(l_i, l_j)`. A **pair of EQUAL odd cycles `(q,q)` contributes
`gcd(q,q) = q`** — the *maximal* possible pairwise term — whereas distinct coprime
cycles contribute only `1`. Verified:

```
 n=6 : (3,3) -> cross gcd = 3   |  (5,1),(3,1,1,1) distinct -> 1 each
 n=10: (5,5) -> cross gcd = 5   |  (7,3),(9,1)              -> 1
 n=14: (7,7) -> cross gcd = 7   |  (11,3),(13,1)            -> 1
```

So the **equal odd-cycle pair `(q,q)` — the cycle form of the doubled prime `2q = q+q`
— is where the tournament's fixed-point symmetry CONCENTRATES.** A doubled prime is
not two random primes; it is the *resonant diagonal* `q = q`, and resonance (equal
lengths -> `gcd = q`, not `1`) is exactly what makes the symmetry count blow up.
This is the same diagonal as Goldbach's `q+q`.

## 4. LRC: the doubled primes `n = 2q` are the rank-one (cleanest) cases

And it closes onto our own thread. From S542, the LRC channel rank is `omega(n/2)`.
For `n = 2q` (a doubled prime), `n/2 = q` is prime, so **rank `= 1`: a single p-adic
tower, the cleanest hard LRC instance**, with apex (source-sink) speed `= n/2 = q`.
Verified rank-one for `n = 6, 10, 14, 22, 26`. So `n = 14 = 2*7` — the canonical hard
LRC target — is a *doubled prime*, and that is *why* it is the simplest non-trivial
case (one tower, apex 7).

## The answer: the importance of the doubled primes

A doubled prime `2q` is the **resonant diagonal** that recurs in every layer:

| layer | the doubled prime `2q` is... |
|---|---|
| additive | the **diagonal of Goldbach** (`q+q`); the even step in Lemoine `odd = p + 2q` |
| tournaments | the **equal odd-cycle pair `(q,q)`** — maximal Burnside cross-term `gcd = q` |
| twins (S521) | the wheel `6 = 2*3` (the doubled prime `3`); twin centers on `6k` |
| LRC (S542/S530/S541) | `n = 2q` = **rank-one** tree, apex `= q = n/2`, the cleanest case |

The unifying fact: **doubling is pairing.** `q -> 2q` turns a single odd prime into a
*resonant pair* `(q,q)`, and pairs of equal things are precisely where `gcd` is
maximal, where Goldbach is diagonal, where the LRC tower is rank-one. The lone `2`
(the only even prime) is the **parity-fixer**: it is what converts the even
2-object story into the odd 3-object story (`p + q + q`), and what halves `n` to the
apex. So the doubled primes are not a side species — they are the **diagonal / the
self-paired / the resonant locus**, the place all three subjects (additive primes,
odd-cycle symmetry, LRC tightness) put their extremal structure.

## Open (→ HYP)
- Is the AP-tightness of LRC at `n = 2q` (rank-one) governed by the `(q,q)` cycle
  type's dominance in the Burnside (the apex speed `q` = the repeated cycle length)?
  i.e. does the doubled-prime cycle type `(q,q)` control the loneliest configuration
  at `n = 2q`?
- Lemoine `odd = p + 2q` as the "one odd cycle + one doubled cycle" decomposition of
  an odd tournament's symmetry: does the count of Lemoine reps track an odd-`n`
  tournament invariant the way Goldbach reps might track an even-`n` one?

## Anchor
`04-computation/doubled_primes_odd_cycles_goldbach_s546.py` (+ `.out`): A000568 via
odd-cycle Burnside (n=3..8); equal odd-cycle pairs `(q,q)` give maximal `gcd=q`;
Goldbach/Lemoine + doubled primes; `n=2q` rank-one LRC. Builds on the repo's odd-cycle
Burnside / OCF, S521 (twin wheel), S542 (rank), S530/S541 (apex).
