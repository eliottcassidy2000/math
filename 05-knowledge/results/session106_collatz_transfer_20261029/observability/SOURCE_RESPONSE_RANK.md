# All-source survival has unbounded linear memory

**PROVED elementary statement; FINITE-EXACT independent trajectory controls.**
This is a source-faithful observer lower bound, not a lower bound for deciding
Collatz convergence, computing one scalar average, or constructing nonlinear
integer formulas. The rank-factorization argument is classical linear algebra;
we claim an explicit application and witnesses, not literature priority.

## 1. A native source at every exact stopping clock

Let T(n)=(3n+1)/2^v2(3n+1) on positive odd integers, stopped at 1. Define
τ(n) as its first hitting time of 1, allowing infinity. There is an explicit
strictly increasing sequence x_j with x_0=1 and T(x_(j+1))=x_j, hence τ(x_j)=j.

Given odd y>0 with 3 not dividing y, try x=(2^a y−1)/3:

* for y=1 use a=4, giving x=5;
* for y=1 mod3 and y>1 try a=2,4;
* for y=2 mod3 try a=3,5.

Both candidates in each pair are integral and odd. Their exponents differ
by two. They cannot both have numerator zero mod9: if 2^a y=1 mod9,
then 2^(a+2)y=4 mod9. Thus at least one candidate is not divisible by3.
It is larger than y (the only equality possibility is y=1,a=2, explicitly
excluded). Since y is odd, the valuation of 3x+1 is exactly a, not merely
at least a. Choose the first admissible candidate. Every step satisfies
x<=(32/3)y, so x_j<=(32/3)^j. The first values are

    1, 5, 13, 17, 181, 241, 1285, 6853, 9137, 24365, ...

No convergence assumption is used; these sources have supplied native
finite trajectories by construction.

## 2. Exact rank obstruction

Let f_i(n)=1[τ(n)>i]. For any N>=1 consider rows i=0,...,N−1 and
columns x_j, j=1,...,N. The response matrix is

    M_(i,j)=1[j>i].

It is upper triangular with all diagonal entries one (row i matches
column j=i+1). Its rank is N over every field.

Suppose a representation assigns each source a vector e(n) in a fixed
d-dimensional vector space and supplies linear observers l_i with
l_i(e(n))=f_i(n) for every positive odd source and all i<N. Then M factors
through that space, so d>=N. In particular no fixed-dimensional linear
state, even with time-dependent output rows, preserves all native survival
questions for every source and every clock. A time-homogeneous recurrence
f_i(n)=l A^i e(n) is a special case and is equally excluded.

The result persists under multiplying each source column by any nonzero
weight. Keeping strictly positive original-source atoms does not make
these independent response channels disappear. Aggregating all columns
into one scalar series can erase the obstruction; we do NOT infer a
nonrationality result for that scalar generating function.

A finite deterministic observer with d states similarly cannot encode
sources with more than d distinct exact stopping times if its outputs
are determined by state. Unbounded integer coordinates or growing state
spaces are entirely compatible with the theorem.

## 3. What transfers from the tournament work

Session105's exact path-word attachment and earlier weighted-language
minimization both distinguish a small scalar counting observer from a
representation preserving every later named test. The map here is exact:

* objects: marked positive Collatz sources;
* continuation: one more native odd step;
* observations: still active after i steps;
* representation: linear source vector and continuation/output operators;
* necessary information: independent hitting-clock response channels;
* decisive test: the explicit triangular matrix above.

Thus a compact A000568 counting recurrence does not automatically offer
a bounded-dimensional all-source Collatz recurrence. The useful design
choice is a growing guarded carrier, or a decision certificate that proves
only the required inequality and does not promise full future reconstruction.
This matches the inherited theorem that exact affine carry plus clocks
reconstructs an entire valuation word: a handful of unbounded integers
can be algebraically compact without having bounded information capacity.

The theorem does not say that termination needs infinitely many stored
numbers. A single suitable unbounded integer rank would suffice. It also
does not prohibit a bounded-size formula defining that integer rank.

## 4. Separate uniformity obstruction

For each k, source n=2^(k+1)−1 has its first k exact valuations equal to1:

    T^i(n)=3^i 2^(k+1−i)−1  (0<=i<=k).

It survives all k steps. This is a changing-source witness, not one
infinite surviving orbit. It explains why a uniform stopping clock and
ordinary ℓ¹ operator-norm contraction fail; the accompanying spectral
note proves that distinction. By contrast, section2 uses sources which
are ALL known to terminate and still forces unbounded observer rank.

## 5. Provenance, literature and checks

Inherited internal sources actually read:

* Session105 `metagraph/WEIGHTED_NATIVE_RECURSION.md`: faithful attachment
  versus endpoint-class quotient, and native proposal versus refresh clock.
* Session101 `continuations/CONTINUATION_STATE_AND_SYMMETRY.md`: exact carry
  reconstruction, and denominator101→35/175 continuation hostile.
* Session104 `descent/POSITIVE_ATOMS_AND_INTEGER_DESCENT.md`: source atoms,
  exact rank obligations and inherited expanding-word controls.

For external scale, [Tao, arXiv:1909.03562v7](https://arxiv.org/abs/1909.03562v7)
was checked on 2026-10-29. Its theorem concerns almost all sources in
logarithmic density, reaching every prescribed bound tending to infinity.
That theorem is not all-source atomic extinction. Its abstract identifies
first-passage transport and 3-adic random walks as proof ingredients; no
claim is made here that the present finite matrices improve that theorem.
The inverse-tree perspective is classical; see
[Lagarias's survey](https://www.cecm.sfu.ca/organics/papers/lagarias/).

Run `python3 check_observability.py` and `python3 -O check_observability.py`.
The checker builds80 exact inverse levels, independently replays their
native orbits, computes rational matrix ranks through40, and checks40
long-survival Mersenne controls. All checks use explicit exceptions, not
assertions removed by optimized Python. The all-order claims follow from
the proofs, not finite extrapolation.
