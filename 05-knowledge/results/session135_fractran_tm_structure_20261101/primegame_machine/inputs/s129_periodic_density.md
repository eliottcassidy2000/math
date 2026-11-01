# Uniform cycle parameters and exact density of periodic points

**PROVED consumer of the independently audited [complete cycle classification](../partition/README.md); FINITE-EXACT integer-level controls.** Periodic points are states already on a cycle. This is not a density of eventual basins.

## 1. Two coprimality formulas classify every cycle

Let

`D2=3*7*13*17*19*23*29=58815393`,
`D6=3*5*7*17*29=51765`.

Every least-period-two cycle is uniquely

`11M <-> 13M`, with `gcd(M,D2)=1`.

Every least-period-six cycle is uniquely the six displayed base states multiplied by X:

`(4301,17765,16302,9867,8671,23023)*X`,

where `gcd(X,D6)=1` and `38` does not divide X.

These are equivalent to the exponent restrictions in the parent theorem. In the two-family, arbitrary powers of2,5,11 and spectator primes are allowed, and the other seven program primes are excluded from M. In the six-family arbitrary powers of2,11,13,19,23 and spectator primes are allowed, but2 and19 cannot both divide X. No other prime may divide X. Here spectator means any prime outside the ten program primes; it is not a bounded list.

The unique state immediately before instruction G in a six-cycle is4301X, proving uniqueness of that parameter. The two-cycle has one state with no13 and one with exactly one13, proving uniqueness there. Within each orbit the base states are distinct. Distinct phases of different allowed scalings cannot overlap: a common state would have the same deterministic successor orbit, then the same unique parameter and phase. The two families cannot overlap because their least periods differ.

## 2. Exact finite counts and natural densities

For a finite prime set P write

`C_P(Y)=sum_{d|product(P)} mu(d)*floor(Y/d)`.

For P6={3,5,7,17,29}, define

`C6(Y)=sum_{d|D6} mu(d)*(floor(Y/d)-floor(Y/(38d)))`.

The subtraction is valid since38 is coprime to D6. If A2(B),A6(B) count periodic states <=B of least period2,6, then the exact formulas are

`A2(B)=C_P2(floor(B/11))+C_P2(floor(B/13))`,
`A6(B)=sum_{b in (4301,17765,16302,9867,8671,23023)} C6(floor(B/b))`.

Unique phase decomposition proves these formulas without multiplicity corrections. Passing to the limit gives scale densities

`rho2=product_{p|D2}(1-1/p)=1216512/2800733`,
`rho6=(1-1/38)*product_{p|D6}(1-1/p)=18944/46835`.

Thus the exact natural densities of periodic **points** are

`delta2=rho2*(1/11+1/13)=2654208/36409529`,
`delta6=rho6*sum_b(1/b)=37396469504/151504041213525`,
`delta_total=144064009693952/1969552535775825`.

These are elementary finite inclusion-exclusion formulas, not probabilistic assumptions about independent trajectories. The floor errors are uniformly bounded by the finite numbers of terms. No claim about the density or decidability of eventual basins follows.

## 3. A fixed-cost periodic-point decision procedure

No factorization of the input is needed. For each of the eight base phases b, test b|N; if it holds, compute N/b and apply the corresponding fixed-modulus gcd test (and38-divisibility test in the six-family). At most eight divisions by fixed constants and fixed-modulus gcd tests decide periodic-point membership and return its orbit parameter and least period. Equivalently F^6(N)=N tests membership using six program steps, because both allowed least periods divide6. This decides whether N is already periodic, not eventual basin membership. The periodic-point set is a finite union of residue classes: each scaled coprimality condition is periodic modulo its fixed scale modulus times b.

## 4. Exact independent checks and asymptotic boundary

The standalone checker compares every integer1..100000 against actual first-eligible execution. It verifies that F^6(N)=N iff N lies in exactly one displayed phase family, and checks the least period directly. It finds7288 period-two points and28 period-six points; the exact inclusion-exclusion counts agree. Normal and optimized outputs are identical, with exception-based checks active under optimization. The unbounded classification is supplied by the parent proof, not inferred from this range.

Because the last FRACTRAN denominator is1, every positive integer has an infinite trajectory. Each trajectory either eventually enters one of the classified cycles or tends to infinity. This last dichotomy is elementary for any deterministic map on positive integers: if any state repeats the trajectory is eventually periodic; otherwise each bounded set is visited only finitely often. It does not decide the alternative for an arbitrary input or assert monotonic growth.
