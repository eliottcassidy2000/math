# HYP-2229: Basel-Type Pi Identities Are a Period-Carrier Duality

**Status:** OPEN method hypothesis with exact classical identities and finite
carrier audit (codex-2026-06-05-S653).

## Claim

The infinite Basel family is best understood as a three-face carrier:

```text
disjoint elementary packets <-> power-sum moments <-> Bernoulli/p-adic sieve.
```

The classical first identity

```text
zeta(2) = pi^2 / 6
```

is not isolated.  It is the first member of the even zeta family

```text
zeta(2k) = c_k * pi^(2k),
```

and it is also the first member of the cleaner elementary family

```text
zeta({2}^m) = pi^(2m)/(2m+1)!.
```

The sine product carries both:

```text
sin(pi*x)/(pi*x) = product_{n>=1} (1 - x^2/n^2).
```

Its coefficients are elementary symmetric sums of `{1/n^2}`.  Taking the
logarithmic derivative turns those elementary/disjoint packets into the power
sums `zeta(2k)`.  Newton identities are the exact algebraic bridge.

## S653 Evidence

`04-computation/basel_pi_family_carrier_s653.py` computes the first twelve
even-zeta coefficients in two ways:

1. Euler's Bernoulli formula

   ```text
   zeta(2k) =
     (-1)^(k+1) * B_(2k) * (2*pi)^(2k) / (2*(2k)!).
   ```

2. Newton identities from the elementary coefficients

   ```text
   e_m = zeta({2}^m)/pi^(2m) = 1/(2m+1)!.
   ```

The two computations agree for `k=1..12`, giving:

```text
1/6, 1/90, 1/945, 1/9450, 1/93555,
691/638512875, ...
```

The disjoint-packet family begins:

```text
1/6, 1/120, 1/5040, 1/362880, 1/39916800, ...
```

The same script rechecks the repo's von Staudt chain:

```text
6 -> 42 -> 1806 -> 1806.
```

At weight `6`, von Staudt-Clausen selects primes `{2,3,7}`, matching the older
`B_6=1/42` Hurwitz/tournament-prime note.

## Repo Interpretation

The useful analogy is not "pi equals a tournament invariant."  The useful
analogy is that both Basel and OCF package a hard scalar as an evaluation of a
retained disjoint-packet object.

For Basel:

```text
elementary packets of {1/n^2} -> sine product -> pi period -> zeta moments.
```

For tournaments:

```text
disjoint odd-cycle collections -> I(Omega(T),2) -> H(T).
```

The elementary face is the side channel.  The power-sum face is the scalar
moment.  Bernoulli numbers add rational descent and p-adic denominator data via
von Staudt-Clausen.

## Consequence

Basel-type identities are a good model for "scalar plus retained packet"
proofs:

- if one only remembers `zeta(2k)`, the disjoint coefficient structure is
  hidden;
- if one remembers the sine product, the whole infinite family becomes one
  product identity;
- if one remembers Bernoulli denominators, the p-adic obstruction primes are
  visible.

This suggests a route for future tournament and LRC work: when a scalar looks
mysteriously rational, look for the product/elementary object whose
log-derivative or moment map produced it.

New nearby example: the triangular power-tower anchor problem can be rewritten
in midpoint form as

```text
F_p(c,n) = c^p - 2 * sum_{r odd} binom(p,r) c^(p-r) S_r(n),
```

so the entire asymptotic correction is carried by odd Faulhaber/Bernoulli
moments. This is the finite discrete analogue of the same “retain the moment
carrier before trusting the scalar” rule. See
`03-artifacts/drafts/triangular-power-anchor-asymptotics.md`.

## Assumption Challenge

The vertices of the Tournament Analysis need not be numbers `k` or zeta values.
S653 uses proof lenses as vertices: sine-product elementary packets,
Newton/log-derivative power sums, Bernoulli/von-Staudt sieve, Fourier boundary,
Dirichlet/modular periods, p-adic zeta audit, and raw numeric PSLQ.

The tournament over these proof lenses is transitive and ranks the sine-product
elementary packet carrier first, then Newton/log-derivative, then
Bernoulli/von-Staudt.  This is the same methodological bias as recent repo
work: keep the side channel before trusting the scalar.

**See:** `04-computation/basel_pi_family_carrier_s653.py`;
`05-knowledge/results/basel_pi_family_carrier_s653.out`;
`07-reflections/basel-pi-family-carrier-s653.md`; HYP-2228, HYP-2227,
HYP-2214, HYP-2211, HYP-1627, T252, THM-002.
