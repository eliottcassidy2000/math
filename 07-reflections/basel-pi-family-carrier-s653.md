# Basel Pi Family Carrier S653

I read "basil" as Basel: the Basel problem and its infinite family of
pi-identities.  The obvious family is

```text
zeta(2k) = rational * pi^(2k).
```

But the better carrier is one level lower:

```text
zeta({2}^m) = pi^(2m)/(2m+1)!.
```

That second formula is the elementary/disjoint face.  It comes straight from
Euler's sine product

```text
sin(pi*x)/(pi*x) = product_{n>=1} (1 - x^2/n^2).
```

The coefficient of `x^(2m)` is the elementary symmetric sum over distinct
indices:

```text
sum_{n1<...<nm} 1/(n1^2 ... nm^2).
```

Then Newton identities turn elementary symmetric sums into power sums, so the
ordinary even zeta values are the moment/log-derivative face of the same object.
That is the part that feels most alive for the repo.  It is the analytic cousin
of the OCF identity:

```text
H(T) = I(Omega(T),2) = sum alpha_j 2^j.
```

In OCF, `alpha_j` counts disjoint odd-cycle packets.  In Basel, the coefficient
`1/(2m+1)!` counts disjoint reciprocal-square packets after the pi-period
normalization.  The scalar moment is not primary.  The packet object is primary.

The Bernoulli side is a different retained channel.  Euler's formula

```text
zeta(2k) =
  (-1)^(k+1) B_(2k) (2*pi)^(2k)/(2(2k)!)
```

is rational descent plus p-adic memory.  Von Staudt-Clausen says the denominator
of `B_(2k)` is the product of primes `p` with `(p-1)|2k`.  That makes the old
repo note `B_6=1/42` feel less like numerology: weight `6` is the first level
where the primes `{2,3,7}` cohabit the Bernoulli denominator, and the chain

```text
6 -> 42 -> 1806 -> 1806
```

is a denominator-dynamics version of "retain the local obstruction primes until
the carrier closes."

The useful slogan:

```text
Basel is not one miracle sum.
Basel is period = product packets = moments = Bernoulli local data.
```

That also clarifies how to use the family without overclaiming.  It does not
prove anything new about `pi+e`, and it does not turn pi into a tournament
count.  It gives a proof architecture: when a scalar identity is too clean,
look for the product object and the local denominator side channel it forgot.

Next target: build a "Basel-OCF dictionary" with three columns:

```text
sin product packets      <-> disjoint odd-cycle packets
log derivative zeta      <-> H / scalar spectra
Bernoulli denominators   <-> local obstruction product ledgers
```

If that dictionary is any good, it should suggest new companion sequences:
elementary packet counts, power-sum moments, and local-prime denominator
closures for the same finite carrier.
