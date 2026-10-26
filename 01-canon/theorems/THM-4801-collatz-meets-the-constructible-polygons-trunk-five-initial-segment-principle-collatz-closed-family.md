---
id: THM-4801
title: "Collatz meets the constructible polygons (n with phi(n) a power of 2, OEIS A003401: 3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24, 30, 32, 34, ...). (1) Trunk five, unconditional: the Collatz trunk number (4^k - 1)/3 (the odd n with 3n+1 a power of 2) is constructible iff k is one of 1, 2, 4, 8, 16; then it equals F_1 F_2 ... F_(j-1) for k = 2^(j-1), and the Collatz step on it is Goldbach's product identity F_0 F_1 ... F_(j-1) + 2 = F_j shifted by one, because the multiplier 3 is F_0. (2) Initial-segment principle: 2^k - 1 is constructible iff k is one of 1, 2, 4, 8, 16, 32 (unconditional, because F_5 is composite), while 2^k + 1 is constructible iff it is a Fermat prime, so 'k is one of 1, 2, 4, 8, 16' is there equivalent to the completeness of the five Fermat primes. Families built on initial runs of Fermat numbers are finite for free; families that test them one at a time carry the open problem. (3) Collatz-closed family: below 2^(2^33) (everywhere, if the five Fermat primes are complete), the constructible n whose whole Collatz orbit stays constructible are exactly 2^a n with n one of 1, 3, 5, 85, 21845, 1431655765, which are also exactly the odd constructible n with 3n+1 constructible."
status: >
  PROVED (elementary; proofs below). It uses Gauss-Wantzel (CITED), Euler's factor 641 of F_5 (CITED), the pairwise
  coprimality of Fermat numbers, and, in (2) only, Zsigmondy's theorem for 2^(2k) - 1 (CITED; its exception k = 3 is
  treated directly).
  INDEPENDENTLY AUDITED (05-knowledge/results/s58_fermat_pepin_collatz_20261026/AUDIT.md, MISTAKE-673): all
  statements and finite checks confirmed with independent code (totient test). Corrected: a wrongly invoked
  squarefreeness in (1), replaced by coprimality; the Zsigmondy exception k = 3 in (2). The identity 2^(2^j) - 1 = F_0 ... F_(j-1) is classical (THM-4749).
  The observations (1)-(3) are new as far as the repository and the searches made here show; the trunk link and the
  Collatz-closed family are not stated in klein's FIVE_LADDER, the S10 Gilbreath/Fermat/Platonic note or the doubling-
  tower fives inventory.
  FINITE-EXACT: collatz_constructible_20261026.py.
  - (1) for k <= 4096.
  - (2) for k <= 4096.
  - (3) on all 32 odd constructible numbers below 2^(2^33).
  CONTEXT (KNOWN, not claimed): Somer-Krizek, Czech. Math. J. 54 (2004).
  - Thm 4.1: the squaring map on units mod n has no cycle of length > 1 iff phi(n) is a power of 2, i.e. iff the
    n-gon is constructible.
  - Cor 4.3: F_m is composite iff that digraph has a longer cycle.
  - Thm 4.4: the components have nu_2(lambda(n)) + 1 levels.
  So "Collatz for squaring" holds exactly on the constructible orders, with depths 1, 2, 4, 8, 16 at the five Fermat
  primes.
session: opus-2026-10-26-S58
source: 05-knowledge/results/s58_fermat_pepin_collatz_20261026/README.md
depends_on:
  - 01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md
related:
  - 01-canon/theorems/THM-4793-the-fermat-line-under-collatz-valuation-two-descent-exact-merge-of-consecutive-lines-and-pepin-as-an-affine-cycle-type.md
  - 01-canon/theorems/THM-4800-pepin-as-the-valuation-two-collatz-branch-trivial-cycle-half-turn-fermat-line-arc-mersenne-mirror.md
  - 05-knowledge/results/fermat_five_pepin_collatz_20261025/FIVE_LADDER.md
  - 05-knowledge/results/gilbreath_fermat_platonic_20260926.md
scripts:
  - 04-computation/fermat_pepin_s58/collatz_constructible_20261026.py (+ collatz_constructible.json)
---

# THM-4801: Collatz meets the constructible polygons

**Constructible orders.** By Gauss–Wantzel, the regular `n`-gon (`n >= 3`) is constructible iff
`n = 2^a * (product of distinct Fermat primes)`, i.e. iff `phi(n)` is a power of 2. Every odd constructible
`n < 2^(2^33)` divides `2^32 - 1 = 3 * 5 * 17 * 257 * 65537`: its prime factors are Fermat primes below `2^(2^33)`,
hence among `F_0, ..., F_32`, and `F_5, ..., F_32` are composite.

**(1) The trunk five.** The odd numbers whose Collatz successor `3n+1` is a power of 2 are the trunk numbers
`t_k = (4^k - 1)/3`.
- **Claim:** `t_k` is constructible iff `k` is one of `1, 2, 4, 8, 16`.
- **`k` must be a power of 2.** Suppose `k` has an odd prime factor `r`. Then `2^r - 1` divides `4^k - 1`. Every prime
  factor `q` of `2^r - 1` has `ord_q(2) = r`, since `r` is prime and `q` does not divide `2 - 1`; Zsigmondy is not
  needed. Then `r | q - 1`, so `q` is not a Fermat prime, and `q != 3` because `ord_3(2) = 2`. So `q | t_k`, and
  `t_k` is not constructible.
- **The power-of-2 case.** For `k = 2^(j-1)` we have `4^k - 1 = 2^(2^j) - 1 = F_0 F_1 ... F_(j-1)`, so
  `t_k = F_1 ... F_(j-1)`.
  - If `F_1, ..., F_(j-1)` are all prime, `t_k` is a product of distinct Fermat primes and so constructible.
  - If some `F_i` among them is composite, its prime factors are not Fermat primes, because Fermat numbers are
    pairwise coprime and every Fermat prime is some `F_k`. Then `t_k` is not constructible.
  - Squarefreeness of Fermat numbers, which is open beyond small indices, is not used.
  - Since `F_5` is composite, `t_k` is constructible iff `j <= 5`.
- So the constructible trunk numbers are `1, 5, 85, 21845, 1431655765`, with `3t + 1 = 4, 16, 256, 65536, 2^32`,
  i.e. `F_j - 1` for `j = 1..5`.
- Because Collatz multiplies by `3 = F_0`, the step `t -> 3t + 1` on these numbers *is* the classical identity
  `F_0 F_1 ... F_(j-1) = F_j - 2`, read as `3 (F_1 ... F_(j-1)) + 1 = F_j - 1`.
- The five indices `k = 1, 2, 4, 8, 16` are the owner's second floor `{1, 2, 4, 8, 16}`.

**(2) The initial-segment principle.**
- **`2^k - 1`.** It is constructible iff `k = 2^j` with `F_0, ..., F_(j-1)` all prime, by the same order argument
  and coprimality. Hence iff `k` is one of `1, 2, 4, 8, 16, 32`, and this is unconditional.
- **`2^k + 1`.** Let `k = 2^a r` with `r > 1` odd.
  - For `k != 3`, Zsigmondy gives `2^(2k) - 1` a primitive prime divisor `q`. It divides `2^k + 1` (not `2^k - 1`),
    and its order `2k` is not a power of 2, so `q` is not a Fermat prime.
  - Zsigmondy's one exception, `2^6 - 1`, is the case `k = 3`, and `2^3 + 1 = 9 = 3^2` is not squarefree.
  - So `2^k + 1` is constructible iff it is a Fermat prime, and "`k` is one of `1, 2, 4, 8, 16`" here is
    *equivalent* to the completeness of `{3, 5, 17, 257, 65537}`.
- **The principle.** Any family that needs a whole initial run of Fermat numbers to be prime is unconditionally
  finite, because the run breaks at `F_5`. This covers the trunk, `2^k - 1`, the products `2^(2^j) - 1`, and the
  Sierpiński row `2^j - 1`. A family that tests Fermat numbers one at a time carries the open problem: `2^k + 1`, the
  prime rows `2^i` of the binary Sierpiński triangle, the Paley and rigid-circulant families, and Pépin.

**(3) The Collatz-closed constructible family.**
- **Claim:** for `n < 2^(2^33)`, every term of the Collatz orbit of `n` is constructible iff `n = 2^a m` with `m` one
  of `1, 3, 5, 85, 21845, 1431655765`. Equivalently, `m` is an odd constructible number with `3m + 1` constructible.
- **Halving stays inside.** If `n` is even and constructible, so is `n/2`, because `phi(n/2)` is `phi(n)` or
  `phi(n)/2`. So an orbit leaves the family only at an odd step `m -> 3m + 1`.
- **The finite check.** Testing the 32 odd constructible `m` (all `< 2^32`, with `3m + 1 < 2^(2^33)`) leaves exactly
  six with `3m + 1` constructible. Their orbits are
  - `3 -> 10 -> 5 -> 16 -> ... -> 1`;
  - `t -> 2^(2^j) -> ... -> 1` for the five trunk numbers.
  All terms are constructible.
- **Without the bound.** A further exception would need a Fermat prime `F_m` with `m >= 33`.
- **Geometric reading.** Collatz on regular polygons sends an `n`-gon to an `n/2`-gon or a `(3n+1)`-gon. It never
  leaves the ruler-and-compass world exactly for these six odd families. Five of them are the trunk five; the sixth,
  the triangle, enters the trunk at the pentagon.
