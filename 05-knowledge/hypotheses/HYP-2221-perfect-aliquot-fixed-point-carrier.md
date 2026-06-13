# HYP-2221: Perfect Numbers Are Aliquot/Divisor Carrier Fixed Points

**Status:** OPEN synthesis / exact finite scout.

## Claim

Complementing HYP-2220's triangular/Vieta aliquot carrier, perfect numbers are
fixed points of the aliquot map

```text
s(n) = sigma(n) - n.
```

Equivalently, they are rows where the divisor-lattice pushforward satisfies

```text
sigma(n) = 2n
```

or where the abundancy carrier is exactly

```text
A(n) = sigma(n)/n = 2.
```

The important side channel is the prime-power product

```text
sigma(n)/n = product_{p^a || n} (1 + 1/p + ... + 1/p^a).
```

So the scalar fixed equation is only the visible quotient.  The divisor
factorization, parity, local prime-power weights, and aliquot graph position
are the retained carrier data.

## Evidence

S645 scans `n <= 100000` and finds aliquot fixed points

```text
6, 28, 496, 8128.
```

These are the even Euclid-Euler section:

```text
n = 2^(p-1) * (2^p - 1)
```

with `2^p-1` prime.  The fixed equation is then a product identity:

```text
sigma(2^(p-1)) * sigma(2^p-1)
= (2^p-1) * 2^p
= 2n.
```

In the same finite window, the aliquot graph also has `13` length-2 cycles
and one length-5 sociable cycle.  This makes the fixed-point reading concrete:
perfect numbers are length-1 sociable cycles, amicable pairs are length-2
cycles, and longer sociable rows are higher-period carrier loops.

## Odd-Perfect Reading

No odd fixed point appears in the window.  The closest odd near-fixed row by
`|sigma(n)-2n|/n` is

```text
n = 32445 = 3^2 * 5 * 7 * 103,
sigma(n)-2n = 6.
```

This is a useful warning: being very close to `A(n)=2` is not enough.  The
missing odd-perfect proof is a side-channel proof about prime-power weights,
congruence obligations, and the aliquot graph, not a scalar continuity
argument around `2`.

## Transfer

This fits the repo carrier program:

- pi/e trace-norm carriers: two shadows reconstruct hidden pair data.
- Goldbach/Lemoine pair carriers: `(E,O)` reconstructs an ordered prime pair.
- LRC/unit-distance carriers: scalar quotients need owner/carry/spine/bulk
  side channels.
- perfect-number carriers: `sigma(n)=2n` needs divisor-lattice and
  prime-power product side channels.

## Next Tests

- Build a local prime-power obstruction ledger for `A(n)=2`, especially odd
  near-fixed rows.
- Treat aliquot cycles as a directed graph and measure which side channels
  separate fixed points from amicable and sociable cycles.
- Search for scalar twins: rows with nearby or equal defect but different
  divisor-lattice structure.
- Add a tournament whose vertices are divisor atoms, prime-power factors,
  aliquot arrows, defect shells, and proof obligations.

**See also:** HYP-2220, `04-computation/aliquot_fixed_point_carrier_s645.py`,
`05-knowledge/results/aliquot_fixed_point_carrier_s645.out`,
`07-reflections/perfect-aliquot-fixed-point-carrier-s645.md`, HYP-2216,
HYP-2215, HYP-2211, HYP-2208.
