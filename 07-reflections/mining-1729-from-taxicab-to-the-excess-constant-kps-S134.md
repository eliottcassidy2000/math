# Mining 1729: from taxicab coincidences to the excess constant

**kind-pasteur-2026-07-26-S134.** Provenance note, not truth source.

The owner asked to mine connections from the appearance of
`tc(T_11) = 1729` in the pure-blue census (THM-2454 data note) and
then to generate directions of the same flavour. The mine had three
levels; only the deepest one produced durable mathematics, and the
record of WHICH levels died is part of the value.

## Level 1: the coincidences (recorded, then released)

- `1729 = 7 * 13 * 19 = 19 * 91`: it contains `91 = 7 * 13`, the
  LRC(14) common stalk modulus, and its factors are
  `(6k+1)(12k+1)(18k+1)` at `k = 1` -- Chernick's universal-form
  Carmichael factors, i.e. the `+1` companions of the first three
  twin ranks `1, 2, 3` (A002822). Three of this month's threads in
  one integer.
- `H(T_11) = 95095 = 5 * 7 * 11 * 13 * 19` -- five consecutive
  primes with `17` skipped; also `95 * 1001` with `1001 = 7*11*13`.
- Verdict: charming, not load-bearing. The Carmichael property was
  KILLED at the next data point (`tc(QR_19)` fails Korselt,
  HYP-9028), which is exactly why the check was run before anything
  was built on it.

## Level 2: the structural residue of the coincidence

What survives scrutiny in `tc(T_11) = 1729` is not its factorization
but its SIZE: `tc` grows explosively along the rotational/Paley
family (`1, 3, 9, 1729, ...`) while blue multiplicity crawls
(`1, 3, 3, 37`), so the purity ratio collapses and `T_5` is the last
pure point (THM-2454 data note). The right question was never "why
is this number famous" but "what is this number's growth law".

## Level 3: the growth law is an excess constant (HYP-9028)

Normalizing `H` by the random-tournament mean `n!/2^{n-1}` turns the
family data into a clean monotone sequence climbing to `~2.52` at
`n = 17..19` in BOTH circulant families, consistent with

```text
H ~ e * n!/2^{n-1} * (1 - alpha/n),   alpha ~ 1.25.
```

Two mechanism readings: the Paley expansion
`H = 2^{-(p-1)} sum_sigma prod (1 + chi(d_i))` whose empty term IS
the random mean and whose corrections are Weil character sums (the
prime-harmonics grammar, third appearance); and the permanent-style
`e` from exact-regularity conditioning. Bonus proved lemma: the
rotation group acts freely on Hamiltonian paths, so `n | H` across
both families.

## The meta-move, for reuse

The 1729 prompt is a template for a repeatable research move:
**when a famous constant appears in your data, (1) enumerate its
coincidences fast and cheaply, (2) kill them with the next data
point if they are killable, (3) keep only the normalization that
made the number comparable across the family -- the surviving
object is usually a growth constant, and growth constants have
mechanisms.** Here the survivor is the Hamiltonian excess constant;
the graveyard holds Carmichael-tc, taxicab structure, and the
missing 17. Compare the S131b detrending lesson (THM-2447): in both
cases the discovery step was dividing out the carrier so the
structure could be seen.
