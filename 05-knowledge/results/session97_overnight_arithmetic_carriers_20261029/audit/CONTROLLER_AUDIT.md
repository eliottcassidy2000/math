# Independent audit: cheap ternary integer exclusion

Status: **PROVED** (elementary argument), plus **FINITE-EXACT** hostile checks.
Independent lane; no producer imports. Source repository pin:
`5d6eb91dd76eee86b91e04d856daa32c32b3d61d`.

## Precisely quantified statement

For any positive odd integer source n0, iterate the native accelerated map
n'=(3n+1)/2^a, where a=v2(3n+1). Start Q=1. When 3 divides a,
replace Q by Q+2; otherwise replace Q by Q-1. Stop at Q=0.
This process terminates for every positive odd source. This is an auxiliary
ternary controller, **not** a descent theorem for the native Collatz map.

## Independent proof

At an active state put P=(n+1)(3/2)^Q. If a=1, then
(n'+1)/(n+1)=3/2, so the leaf operation preserves P exactly.
For every other leaf, a>=2, and

    (n'+1)/(n+1) = (3n+1+2^a)/(2^a(n+1)) < 3/2.

Thus no leaf increases P. For a branch, a belongs to {3,6,9,...}.
The six positive odd sources below 13 have valuations 2,1,4,1,2,1;
none is a branch. Hence n>=13, and

    (n'+1)/(n+1) <= (3n+9)/(8(n+1)) <= 3/7.

The first inequality follows from 2^a>=8; the second is equivalent to n>=13.
Multiplying by (3/2)^2 gives P'/P<=27/28. Equality is attained
at n=13, a=3 (with output n'=5).

After B branches, P <= (3/2)(n0+1)(27/28)^B. Every active state
has P>=3. If infinitely many steps were possible, then infinitely many
branches would be required: after the last branch Q could only decrease.
But infinitely many strict uniform contractions contradict the lower bound.
This proves termination without any probabilistic assumption.

Every branch output is active. Consequently the total number B of branches
satisfies the exact integer certificate

    2 * 28^B <= (n0+1) * 27^B.

Equivalently B<=floor(log((n0+1)/2)/log(28/27)). There are 2B+1
leaves, so the total number of accelerated source steps is N=3B+1.
This bounds accelerated steps, not the ordinary Collatz clock N+sum(a).

## Hostiles and scope

- n=5 has a=4: the assertion 'a>=3 implies n>=13' is false. The
  *multiple-of-three valuation guard* is indispensable.
- n0=3 terminates immediately at 5: stopping does not imply descent.
- n0=1 terminates immediately at 1: a native cycle is not a controller cycle.
- No statement here changes or solves the original odd-valuation ternary rule.
- No random-source weight, input valuation, source state, or carry is changed.
  The potential is a deterministic witness on the exact native trajectory.

`check_controller_independent.py` tests all 50,001 positive odd seeds up to
100001, exact branch ratios, integer branch certificates, equality conditions,
and the full-tree identity. Normal and `python -O` outputs agree. The longest
observed trace is 31 steps (10 branches), largest attaining source 68813.
Finite checks support but do not replace the preceding all-source proof.
