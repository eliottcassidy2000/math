---
id: THM-514
title: The A000568 unit-tail peel isolates a two-parameter core kernel, and odd-prime tails add only divisibility-count statistics
status: PROVED (unit-tail factorization, support count, generating-function factorization, odd-prime tail lemma). VERIFIED computationally through n=100 by a000568_core_tail_ladder_codex.py.
source: codex-2026-06-15-S12
depends_on:
  - THM-305   # Burnside formula for A000568 is used as input
related:
  - THM-512   # Möbius sieve / peel-the-universal-tail viewpoint
  - OPEN-Q-103
  - reflection: the-burnside-unit-tail-is-soluble-and-the-next-obstruction-is-the-prime-tail-ladder-codex
---

# THM-514 — peel the 1-tail, isolate the hard core

Write the Davis/Burnside formula for the tournament count

```text
a(n) = A000568(n) = sum_{lambda odd partition of n} 2^e(lambda) / z(lambda),
```

where for odd cycle lengths `lambda = (l_1,...,l_k)`:

```text
e(lambda) = sum_i (l_i-1)/2 + sum_{i<j} gcd(l_i,l_j),
z(lambda) = prod_j (j^{m_j} m_j!).
```

The user's observation is exact: the `1`-part tail is completely soluble.

## Statement

Every odd partition `lambda` of `n` splits uniquely as

```text
lambda = mu U 1^r,
```

where `mu` has odd parts at least `3`. Let

```text
m = |mu|,   t = ell(mu),   r = n-m.
```

Then:

```text
e(lambda) = e(mu) + C(r,2) + r t,
z(lambda) = z(mu) r!.
```

Hence if

```text
B[m,t] = sum_{mu: |mu|=m, ell(mu)=t, odd parts >=3} 2^e(mu)/z(mu),
```

then

```text
a(n) = sum_{m,t} B[m,t] * 2^(C(n-m,2) + (n-m)t) / (n-m)!.
```

Equivalently, after clearing the core denominator with

```text
K[m,t] = m! B[m,t] in Z_{>=0},
```

we get the exact integer form

```text
n! a(n) = sum_{m,t} C(n,m) K[m,t] 2^(C(n-m,2) + (n-m)t).
```

Moreover

```text
K[m,t] = sum_{mu: |mu|=m, ell(mu)=t, odd parts >=3}
         (m!/z(mu)) 2^e(mu),
```

so `K[m,t]` is literally a weighted count of labelled permutations on `m`
letters whose odd cycle lengths are all at least `3`.

The active support of `(m,t)` is exactly

```text
(m,t) = (0,0) or [ t>=1, m>=3t, m==t (mod 2) ].
```

So the number of active outer states is

```text
S(n) = 1 + sum_{t=1..floor(n/3)} ( floor((n-3t)/2) + 1 ).
```

In particular:

```text
p_odd(100) = 444793  ->  S(100) = 834.
```

That is the exact outer collapse reported by the user.

## Proof

For `lambda = mu U 1^r`, the contribution of the `1`-cycles is completely explicit.

1. A `1`-cycle contributes no inside-pair orbit:

   ```text
   (1-1)/2 = 0.
   ```

2. Between two `1`-cycles, `gcd(1,1)=1`, giving

   ```text
   C(r,2).
   ```

3. Between a `1`-cycle and any part of `mu`, the gcd is also `1`, and there are
   `r * ell(mu) = r t` such pairs.

Adding these gives

```text
e(lambda) = e(mu) + C(r,2) + r t.
```

For the centralizer factor, `1`-cycles contribute exactly `r!`, so

```text
z(lambda) = z(mu) r!.
```

Substituting into the Burnside sum and grouping by `(m,t)` yields the displayed
formula for `a(n)`. Multiplying by `m!` gives the integer kernel `K[m,t]`, and
then multiplying the identity by `n!` gives the binomial form.

The integrality of `K[m,t]` is not accidental: `m!/z(mu)` is the number of
permutations in `S_m` with cycle type `mu`, so `K[m,t]` is an integer weighted
sum over labelled cycle types.

The support criterion is immediate: a partition of `m` into `t` odd parts at
least `3` exists iff `m>=3t` and `m==t (mod 2)`; conversely every such pair is
realized by

```text
(m-3(t-1), 3, ..., 3).
```

Counting the allowed `m` for each `t` gives the state-count formula. QED.

## Corollary — the universal unit-tail function

Define

```text
F(x) = sum_{r>=0} 2^C(r,2) x^r / r!.
```

Then the ordinary generating function of `a(n)` factors as

```text
A(x) = sum_{n>=0} a(n) x^n = sum_{m,t} B[m,t] x^m F(2^t x).
```

The tail function is universal, independent of the hard core. It also satisfies

```text
F'(x) = F(2x),
```

by termwise differentiation.

So the entire difficulty of A000568 is isolated in the finite-core coefficients
`B[m,t]`; the `1`-tail itself is a closed-form universal background.

## Odd-prime tail lemma

The same peel works one prime at a time.

Let `p` be an odd prime and split a partition with no `1`-parts as

```text
mu = nu U p^s,
```

where `nu` has no part equal to `p`. Write

```text
u = ell(nu),
c_p(nu) = #{parts of nu divisible by p}.
```

Then

```text
e(mu) = e(nu) + s(p-1)/2 + p C(s,2) + s( u + (p-1)c_p(nu) ),
z(mu) = z(nu) p^s s!.
```

Proof: each `p`-cycle contributes `(p-1)/2` internally; two `p`-cycles interact
by `gcd(p,p)=p`; and for any odd part `L` of `nu`,

```text
gcd(p,L) = 1 + (p-1)[p|L].
```

So the next exact refinement after the `1`-tail is not a whole new partition:
it asks only for one extra statistic, the count of remaining parts divisible by
`p`.

For `p=3`, this yields an exact `3`-tail reconstruction of the core kernel from
a `3`-free kernel indexed by `(m,t,c_3)`. The stored computation verifies that
reconstruction exactly through mass `100`.

Numerically at mass `100`, the residual partition universe with odd parts at
least `5` has size `7551`, while the exact `3`-free state space carries only
`2049` active `(m,t,c_3)` states. So the prime-tail ladder is still genuinely
compressive after the `1`-tail peel; it just no longer closes on `(m,t)` alone.
