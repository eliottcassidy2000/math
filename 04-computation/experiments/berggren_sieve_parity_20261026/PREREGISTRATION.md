# Pre-registered predictions (written and committed before the runs that test them)

mac-mini-2026-10-25-pyth, 2026-10-26. Filed under HYP-9418 (3)-(4). Each prediction names its computation; nothing
listed here had been computed when this file was committed.

## (3) Twin sieve at the critical level: the Mertens ratio

Let `pi2(X) = #{p <= X : p, p + 2 both prime}` and

    rho(X) = (X/2) * prod_{3 <= q <= sqrt(X)} (1 - 2/q) / pi2(X),

the naive "independent residues" density of the exact level-sqrt(X) sieve on the D ray (odd x, both x and x + 2 free
of odd prime factors <= sqrt(X)), divided by the true count. Model: Hardy-Littlewood `pi2 ~ 2 C2 Li2(X)` and Mertens
`prod ~ 4 C2 e^(-2 gamma) / log^2 sqrt(X)` give `rho(X) ~ 4 e^(-2 gamma) (X / log^2 X) / Li2(X)`, tending to
`4 e^(-2 gamma) = 1.2609...` from below.

Predictions (computed by our own sieve, `twin_sieve_ratio.c`):

- `rho(10^8)` in `[1.100, 1.130]`;
- `rho(10^9)` in `[1.120, 1.145]`;
- `rho(10^10)` in `[1.130, 1.160]`;
- strictly increasing across the three, and all three below 1.2609.

## (4) Branch-character height sum

`S(X) = sum over primitive triples with hypotenuse c <= X of (-1)^(#A)`. Observed through `10^9` (parity_height_sum.out,
parity_height_remainder.out): `S(X) - (#root-ray triples) = R(X)` oscillates in sign with `|R| <= 0.21 sqrt(X)`.
Predictions for the not-yet-computed heights `2*10^9` and `5*10^9` (`parity_height_big.c`):

- `S(X)/sqrt(X)` lies in `[1.00, 1.42]` (that is, `(1/2 + 1/sqrt 2) +- 0.21`) at both heights;
- `S(X) > 0` at both heights (even branch parity stays ahead).

## (5) Ramanujan exceptions of the Berggren graphs G_q (added before the run)

Observed for primes 5 <= q < 1000 (`berggren_ramanujan_scan.out`): G_q is Ramanujan except at q = 331, 569, 719, 809,
887 (5 of 165 primes), with largest ratio max|mu|/(2 sqrt 3) = 1.0133. Model: G_q behaves like a random 4-regular
graph on about q/4 vertices, whose second eigenvalue fluctuates at scale N^(-2/3) around 2 sqrt 3. Predictions for
the primes 1000 < q < 4000 (`berggren_ramanujan_scan.py` with a lower bound), not yet computed:

- the fraction of non-Ramanujan G_q lies in `[0.03, 0.35]`;
- the largest ratio max|mu|/(2 sqrt 3) over the range is below 1.03;
- every G_q in the range has max nontrivial |mu| <= 3.8 (a uniform spectral gap, as Selberg's theorem guarantees in
  some form).

## Erratum (added after the runs; the predictions above are unchanged)

The independent audit (audit_thm4808/REPORT.md §6) confirmed the timing:
- entries (3) and (4) were committed at `2fe36c8088`, before the producing binaries were compiled;
- entry (5) was committed at `a0a2bbd805`, before its output existed.

It also found two labelling slips:
- **Section numbering.** Entries (3), (4) and (5) here correspond to HYP-9418 (3), (2) and (4), not to "(3)-(4)" as
  the header says.
- **Prime count.** "5 of 165 primes" in entry (5) should read 166 primes (`5 <= q < 1000`).

Two caveats on scope:
- Entry (3) predicts quantities already fixed by published data (OEIS A007508) together with an elementary product,
  so it checks arithmetic rather than a model.
- Entry (4)'s `+-0.21` band was drawn from sampled heights. The full-resolution maximum of `|E|/sqrt X` up to
  `10^9` is 0.368.
