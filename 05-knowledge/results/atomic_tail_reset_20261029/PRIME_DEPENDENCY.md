# A genuinely well-founded arithmetic control: the squarefree prime selector

**Status: PROVED elementary coordinate convergence and uniqueness;
FINITE-EXACT uniform bounded-prime test. Global finite-time stabilization
of every integer orbit is not proved here.** No literature-priority claim.

## 1. Recover the corrected operation before changing it

The repaired [March Bernoulli fixed-point note](../../../07-reflections/bernoulli-fixed-point.md)
and incoming [October29 archaeology](../twelve_twenty_two_trace_20261029/ARCHAEOLOGY.md)
use the arithmetic selector

    S(n)=product of primes p with p-1 dividing n.

We use this explicit definition; no theorem about Bernoulli numbers is
needed in the proofs below. Every S(n) is squarefree. If n is the product
of a finite prime set P, its next prime set is

    F(P)={p prime: p-1 is squarefree and Supp(p-1) subset P}.    (1)

Equation (1) also defines an operation on arbitrary, possibly infinite,
prime sets. It is monotone, and every dependency prime is smaller than
the prime being decided. In particular 2 always belongs to F(P).

The correct fixed set is P*={2,3,7,43}, with product 1806. To check it,
enumerate every divisor of 1806; adding one gives

    2,3,4,7,8,15,22,43,44,87,130,259,302,603,904,1807.

Only 2,3,7,43 are prime. Checking just 1807 would be insufficient.
The older proposed fixed set {2,3} fails because it selects 7. The
incoming packet already repaired that error; it is not a new correction.

## 2. A uniform expanding lower bound for the discrepancy labels

Let P_t=F^t(P_0), and let b_t be the least prime of `P_t symmetric_difference P*`,
or infinity when the sets agree. Then

    b_1>=3,
    b_(t+1)>=2b_t+1  for t>=1,                               (2)
    b_t>=2^(t+1)-1   for t>=1.                               (3)

The infinity case in (2) means equality with the fixed set persists.

**Proof.** Both P_t and P* contain 2 for t>=1. If a prime p differs in
the next sets, its squarefree dependency test (1) succeeds on one set
and fails on the other. Thus some dependency q divides p-1 and differs
in the preceding sets. That q cannot be 2. Since p is odd, 2q divides
p-1, so `p>=2q+1>=2b_t+1`. The base bound and induction prove (3).

This is the exact property still missing from the Collatz minimum
problem: a discrepant large label is forced to point to a smaller
discrepant label at the preceding time. The proof retains the same
coordinate predicate at both times, rather than invoking an averaged
probability of contraction.

Consequences:

1. **Unique fixed prime set.** Any fixed set contains 2; if it differed
   from P*, its least discrepancy would violate (2). This includes
   infinite prime sets.
2. **Every fixed prime coordinate stabilizes.** All coordinates p<=H
   agree with P* whenever `2^(t+1)-1>H`, uniformly over P_0.
3. **No other periodic prime-set orbit.** Coordinate convergence makes
   each periodic coordinate constant, so the periodic orbit is P*.
4. **Unique positive integer fixed point, and no other integer cycles.**
   An integer fixed point or cycle is squarefree after the first step,
   so the prime-set conclusions apply. The fixed integer is 1806.

For any summable nonnegative weights a_p attached to the **current prime
labels**, the discrepancy weight is at most
`sum_{p>=2^(t+1)-1} a_p`, hence tends to zero. This does not prove that
the discrepancy set becomes empty at one common finite time. An
unbounded sequence of new prime labels is not excluded by the argument.
We neither assert such an orbit exists nor call its existence an
established open problem in the literature.

## 3. Exact finite universal control and what it does not certify

[prime_dependency.py](prime_dependency.py) computes every dependency
for primes through 1,000,000. Higher primes cannot affect those
coordinates. Monotonicity sandwiches **every** initial prime bank
between the empty and full banks in this window. Their successive sizes
are:

|t|lower size|upper size|least possible discrepancy|
|---:|---:|---:|---:|
|0|0|78498|2|
|1|1|29397|3|
|2|2|5745|7|
|3|3|935|43|
|4|4|148|283|
|5|4|30|1319|
|6|4|11|55399|
|7|4|4|none in the window|

Thus seven updates settle every prime coordinate through one million
for every initial prime set, including initial sets containing arbitrarily
large primes. This is a finite-coordinate universal statement, not a
seven-step stabilization theorem for the full integer.

## 4. The transfer to Collatz is a target, not a substitution

For the Collatz E_k sets, a legal predecessor relation

    n in E_(k+r) -> y in E_k,   y+1 <= (n+1)/b

would give `s_(k+r)+1 >= b(s_k+1)`, the analog of (2). The
[seven-return candidate](ATOMIC.md) specifies a finite-tested scale.
THM-4555 supplies some equal-clock arithmetic deletions, but preservation
of the growing-return guard is still missing for general sources.

The prime-set result also exposes a forbidden shortcut: its small
weighted discrepancy can come from moving to larger **current labels**.
In Collatz the weight belongs permanently to the original source; a
persisting exceptional integer cannot become weightless merely because
its orbit grows. Similarly, moving a geometric port or changing the
Keller class can make an observation small without solving the original
boundary or polynomial problem. The source-retaining analog is the
predecessor map and repayment functional, not the bare convergence of
a weight on moving coordinates.
