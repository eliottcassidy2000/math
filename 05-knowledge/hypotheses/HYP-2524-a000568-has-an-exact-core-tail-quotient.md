# HYP-2524: A000568 has an exact core-tail quotient

**Status:** CONFIRMED exact identity through `n=20` with closed-form outer
quotient (`a000568_core_tail_quotient_codex.py`, codex-2026-06-15).

## Claim

Every odd partition contributing to A000568 has a canonical decomposition

```text
lambda = mu union 1^r,
```

where `mu` is the core consisting of odd parts `>= 3`, and `r` is the number of
ones. If `m=|mu|` and `t=#parts(mu)`, then the Burnside weight factors exactly:

```text
e(lambda) = e(mu) + C(r,2) + r*t,
z(lambda) = z(mu) * r!.
```

Therefore

```text
a(n) = Σ_{m<=n} Σ_t B[m,t] * 2^( C(n-m,2) + (n-m)t ) / (n-m)!,
```

where the **core kernel**

```text
B[m,t] = Σ_{mu core, |mu|=m, #parts(mu)=t} 2^e(mu) / z(mu)
```

depends only on odd parts `>=3`.

This is an exact second quotient on odd partitions: the outer sum is indexed only
by `(m,t)`, not by the full odd partition.

## Exact evidence

`04-computation/a000568_core_tail_quotient_codex.py` verifies the formula against
known exact values through `n=20`, including:

```text
a(14)=28401423719122304
a(15)=31021002160355166848
a(20)=645022068557873570931850526424042500096.
```

All checks through `n=20` pass exactly.

## Compression

The outer quotient is much smaller than the full odd-part partition space:

```text
n=20: odd partitions 64       -> active (m,t) states 34    (1.88x)
n=50: odd partitions 3658     -> active (m,t) states 209   (17.50x)
n=80: odd partitions 77312    -> active (m,t) states 534   (144.78x)
n=100: odd partitions 444793  -> active (m,t) states 834   (533.32x)
```

So the full odd-part Burnside sum has an exact `O(n^2)` outer address.

## Interpretation

This identifies the true locus of complexity:

- the **tail of ones** is now solved exactly;
- the remaining hardness is the **core kernel** `B[m,t]`.

So the partition problem has split cleanly into:

```text
outer quotient:  (m,t) endpoint recursion        -- exact and small
inner kernel:    core odd partitions >= 3        -- still hard
```

This is precisely the sort of "endpoint recursion" the repo has been circling:
the long tail of ones is an address variable, not deep combinatorial content.

## Consequences

1. A000568 does have a genuine second quotient on odd partitions.
2. Any future speedup should target `B[m,t]`, not the full odd-part sum.
3. The divisor-profile work of HYP-2523 belongs **inside** the core kernel, not
   outside it.
4. The natural next question is whether `B[m,t]` itself has a further quotient:
   divisor-profile strata, prime-support types, or a low-rank kernel recurrence.

## Assumption challenge

The obvious compression candidate was "compress all odd partitions uniformly."
The exact computation says the cleaner split is asymmetric:

- `1` is a special endpoint species and should be peeled off first;
- only after removing the `1`-tail does it make sense to hunt for a deeper
  quotient on the true core.
