# B_6 = 1/42: Why Bernoulli Knows About Hurwitz

**CORRECTED -- 2026-08-25.** The original reflection made three false
promotions: it assigned the Hamilton quaternion algebra finite ramification
`{2,3,7}`, turned a finite negative `7`-adic valuation into a pole, and called
`B_6` the unique Bernoulli number of denominator `42`. The exact
von Staudt--Clausen calculation survives; the quaternion and pole bridge does
not. See MISTAKE-506 and the finite-exact companion
[`bernoulli_42_and_padic_pole_correction_20260825.py`](../04-computation/bernoulli_42_and_padic_pole_correction_20260825.py).

**Session:** kind-pasteur-2026-03-16-S116n

## The surviving fact

The sixth Bernoulli number is B_6 = 1/42.

## Why

The von Staudt-Clausen theorem (1840): denom(B_{2k}) = prod{p prime : (p-1) | 2k}.

For 2k = 6:
- Divisors of 6: {1, 2, 3, 6}
- Adding 1: {2, 3, 4, 7}
- Primes among these: {2, 3, 7}
- Product: 42

The set `{2,3,7}` here is exactly the von Staudt divisor set for index six.
It is not the finite ramification set of the Hamilton quaternion algebra.

## What the quaternion comparison loses

The Hamilton quaternion algebra `(-1,-1)_Q` ramifies at the real place and at
the finite prime `2`; its reduced finite discriminant is `2`. A maximal
Hurwitz order does not acquire finite ramification at `3` or `7`, and its
arithmetic does not turn `2*3*7` into a quaternion discriminant. Thus the old
claim that von Staudt--Clausen and the Hurwitz order detect one obstruction is
refuted. At most, the word "Hurwitz" and the number `42` formed a mnemonic.

## The p-adic correction

- `zeta(6)=pi^6/945`, and the standard Bernoulli formula carries the factor `7`.
- `zeta(-5)=-1/252=-1/(2^2*3^2*7)`.
- `B_42` has `43` in its denominator because `43-1=42` divides the index.

For the standard Kubota--Leopoldt interpolation convention,

```text
zeta_7(-5)=-(1-7^5)B_6/6=2801/42,
v_7(zeta_7(-5))=-1.
```

This is a finite element of `Q_7` with a factor `7` in its denominator, not a
pole. The p-adic zeta function's exceptional pole is at `s=1`, not at every
class `s=0 mod 6`.

## Denominator 42 is not unique to B_6

The old uniqueness claim already fails at `B_114`. The divisors `d` of `114`
for which `d+1` is prime give only `2,3,7`, so

```text
denom(B_114)=42.
```

An exact scan over every even index through `1000` finds

```text
6, 114, 186, 258, 354, 402, 426, 474,
582, 654, 762, 834, 894, 942, 978.
```

This is a finite census, not an infinitude claim. The genuine characterization
is simply

```text
denom(B_n)=42
iff {p prime: p-1 divides n}={2,3,7}.                      (1)
```

## Remaining numerical patterns

- B_2 = 1/6: denom = 2*3 = the smallest Hurwitz primes
- B_4 = -1/30: denom = 2*3*5 = the Platonic system
- B_6 = 1/42: denom = 2*3*7 = the Hurwitz system
- B_12 = -691/2730: denom = 2*3*5*7*13 = Platonic AND Hurwitz together

Only the first sentence of the intended sieve picture is mathematical: each
even index selects exactly the primes `p` for which `p-1` divides that index.
The labels "Platonic" and "Hurwitz" in the list above are historical
mnemonics, not proved common mechanisms.

The equality `weight(Delta)=phi(42)=12` is numerically true. No map in the
reflection identifies the modular discriminant with a "Cayley chromatic
number," so that concluding identification is withdrawn.
