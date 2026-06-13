---
source: oracle-2026-06-01-S546o
status: synthesis + computation (doubled primes as the parity hinge across numbers, cycles, and LRC channels)
tags:
  - goldbach
  - lemoine
  - doubled-primes
  - parity
  - odd-even-cycles
  - burnside
  - lonely-runner
  - channels
---

# Doubled Primes as the Parity Hinge: Numbers, Cycles, and LRC Channels

The factor of `2` — the unique even prime — is one operator wearing three masks:
the **doubled prime** `2q` (additive), the **even cycle / antipodal involution**
(combinatorial), and the **halving `n* = n/2`** (LRC). Understanding the doubled
prime means understanding this operator, and it pays off directly on the LRC
frontier.

## The parity dictionary (additive)

```
EVEN  E = p + q        (Goldbach):  odd + odd = even.
ODD   O = p + 2q       (Lemoine/Levy):  odd + EVEN = odd.
```

Verified for `E = 4..28`, `O = 7..29` (`lrc_doubled_primes_parity_cycles_s546.py`):
every even has a Goldbach pair, every odd a Lemoine pair `p + 2q`. The decisive point
is parity: an odd number is `odd + even`, and the **doubled prime `2q` is the unique
prime-structured *even* summand**. So:

> **The doubled primes are the PARITY-COMPLETION of the prime additive basis.** Even
> numbers are reachable by primes alone (Goldbach); odd numbers require a doubled
> prime (Lemoine) to supply the parity-flipping even term. `2q = (the even prime) ×
> (an odd prime)` is exactly where the prime basis "crosses the parity divide."

## The same operator in cycles (Burnside / A000568)

A000568 (tournament iso-classes: `2,4,12,56,456,6880`) is computed by Burnside, and
the cycle parity is decisive (verified): **a permutation with any EVEN-length cycle
fixes `0` tournaments; an ALL-ODD cycle type fixes `2^{pairs}`.** So

> **A000568(n) is an *all-odd-cycle* sum; the DOUBLING (an even cycle = an
> antipodal-paired structure) annihilates the fixed count.**

This is the cycle face of the same dictionary: **odd cycles are the free, atomic
contributors (like the odd primes); even cycles are the "doubled" structures that
kill** — because an even cycle carries the antipodal involution (the `×2`), forcing
ties that admit no fixed orientation. Odd = atomic/free; even = doubled/antipodal.

## The payoff: doubled primes control LRC channel cleanliness

In LRC the character modulus is `n* = n/2` (even `n`) or `n` (odd `n`) — the same
halving. The channel structure (S533/S534) is cleanest when `n*` is **prime** (clean
QR/Paley channels) and messy when `n*` is a prime power (filtered, S534) or composite
(CRT, S17). The table (verified):

```
 n   n*   n* type          n structure          clean prime channels?
  3   3   PRIME            odd prime              YES
  6   3   PRIME            2*3  (DOUBLED PRIME)   YES
 10   5   PRIME            2*5  (DOUBLED PRIME)   YES
 14   7   PRIME            2*7  (DOUBLED PRIME)   YES   <<< open
 16   8   prime-power 2^3  2*8  (2*2^3)           no    <<< open
 18   9   prime-power 3^2  2*9  (2*3^2)           no    <<< open
 22  11   PRIME            2*11 (DOUBLED PRIME)   YES
```

> **Clean prime channel modulus `n*` occurs iff `n` is an ODD PRIME or a DOUBLED
> PRIME `n = 2p`.** The doubled-prime dimensions `n = 4,6,10,14,22,…` are the **even
> shadows of the odd primes** (`n* = p`); they inherit the clean prime-character
> (QR/Paley) channel structure. This is the LRC importance of doubled primes.

And it sharpens the frontier: among the open even cases,
- **`n = 14 = 2·7` is a doubled prime** → `n* = 7` prime → **clean channels** → the
  *easiest* open even case;
- **`n = 16 = 2·8 = 2·2^3`** and **`n = 18 = 2·9 = 2·3^2`** are `2 ×` (prime power) →
  `n*` a prime power → **filtered channels** (the S534 prime-power mess).

So the doubled-prime structure of `n` predicts *which* open even case is tractable:
`14` (doubled prime) is qualitatively cleaner than `16, 18` (twice a prime power).

## The unifying operator

> **`×2` (the even prime) = the antipodal/doubling involution = the universal parity
> operator.** Its three faces: the **doubled prime `2q`** (parity-completion of the
> additive basis, Lemoine); the **even cycle** (`= 2 ×` odd via antipodal pairing;
> kills Burnside Fix, so A000568 is an all-odd sum); and the **halving `n* = n/2`**
> (so `n = 2p` doubled-prime dimensions inherit `p`'s clean prime channels). In every
> register, **odd = atomic/free/prime, even = doubled/antipodal**, and the doubled
> prime `2q` / dimension `n = 2p` is the bridge that carries odd-prime structure
> across the parity divide.

## Verdict / next
- The importance of doubled primes, answered three ways: (additive) the parity
  hinge that completes the prime basis for odd numbers; (cyclic) the antipodal
  doubling whose even cycles annihilate the A000568 Burnside sum; (LRC) the
  dimensions `n = 2p` with clean prime channel modulus `n* = p`.
- Frontier payoff: `n = 14` (doubled prime, clean) is the easiest open even case;
  `n = 16, 18` (twice a prime power) carry the filtered prime-power channels.
- Concrete next: (1) exploit `n = 14`'s prime `n* = 7` (Paley/QR channels, S535) for
  a sharper attack than `16/18`; (2) the even-cycle-kills-Burnside fact as the cyclic
  Lemoine analog — is there a "doubled-cycle" decomposition mirroring `O = p + 2q`?;
  (3) twin-doubled-primes / the `2q` structure vs the twin-Goldbach exceptions (S516).

## Artifacts
```
04-computation/lrc_doubled_primes_parity_cycles_s546.py
05-knowledge/results/lrc_doubled_primes_parity_cycles_s546.out
```
Related: S533/S534 (channels, n* = n/2, prime-power filtration), S535 (QR/Paley),
S529 (cyclotomic), S516 (twin-Goldbach 35 exceptions), S17 ("if 15 were prime"),
CLAUDE.md (Burnside: even cycles -> Fix 0, all-odd -> 2^{pairs}).
