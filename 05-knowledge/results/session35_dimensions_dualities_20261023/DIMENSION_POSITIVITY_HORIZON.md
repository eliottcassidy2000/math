# Dimension arithmetic, signed Burnside, and a precise positivity horizon

Session35, 2026-10-23. **PROVED elementary application of classical symmetric-group theory**, not a novelty claim for the underlying hook-content identity. Finite controls are separately labelled. This note does not identify an intermediate Lie-algebra label with a fractional vector-space dimension.

## Inheritance and question

Session34's Reynolds proof of Burnside counts invariant vectors by averaging permutation traces. Its warning was that irreducible representations are additive atoms, not tensor-primes. The older `05-knowledge/results/cartan_dimension_ladder_s95.md` distinguishes symmetric, antisymmetric, and central dimensions: a scalar dimension loses its decomposition. `05-knowledge/results/collatz_representation_positivity_20261005.md` distinguishes positivity in each finite quotient from pointwise positivity surviving refinement. These mechanisms suggest asking exactly how long a formal dimension can imitate a positive tensor trace.

The board is dimension / trace / tensor power / antisymmetrization / positivity / refinement. The hostile is t=7/2: every permutation-trace test through fourth tensor order is positive, but the fifth antisymmetric trace is negative.

## 1. Signed Burnside interpolates exterior dimensions

Let t be real. On S_k define f_t(sigma)=t^{c(sigma)}, with c the number of cycles, and the matrix

    G_k(t)[sigma,tau] = t^{c(sigma^{-1}tau)}.

For t=d a nonnegative integer, f_t is the character of the action permuting tensor factors of (C^d)^{tensor k}. Thus G_k(d) is positive semidefinite. The signed average is

    (1/k!) sum_sigma sign(sigma)t^{c(sigma)} = binomial(t,k).

Proof: insert the kth letter into a permutation on k-1 letters. Making a new singleton multiplies the signed cycle polynomial by t; inserting it in one of k-1 existing positions changes sign and does not change the cycle number. The recurrence is P_k(t)=(t-k+1)P_{k-1}(t). Dividing by k! gives the identity. At integer d this is the trace of the alternating projector and counts dim exterior^k(C^d).

This is Burnside with a sign character, not a count of ordinary set orbits when t is not an integer.

## 2. Exact finite-depth classification

**Theorem.** For k>=1 and t>=0,

    G_k(t) is positive semidefinite
    iff t is a nonnegative integer or t>=k-1.

Proof of sufficiency: at integers use the actual tensor representation. For general t>=k-1, convolution by f_t is central. On the S_k irreducible indexed by a partition lambda of k its scalar is

    product over boxes (i,j) of lambda of (t+j-i).

To verify this identity without assuming interpolation as an axiom, at every integer d>=k use Schur-Weyl decomposition: the scalar is k! s_lambda(1^d)/dim(Specht_lambda), and the hook-content and hook-length formulas give the displayed product. Both sides are polynomials in d, so equality holds for all t. Every content j-i>=1-k, proving nonnegativity.

For necessity, take t in (m,m+1) with m<=k-2. The principal submatrix indexed by the embedded S_{m+2} is t^{k-m-2}G_{m+2}(t). Its sign vector has eigenvalue (m+2)! binomial(t,m+2)<0. Hence G_k(t) cannot be positive semidefinite. t=0 was already handled.

**Corollary: positivity horizon.** For noninteger t>=0 the first failing tensor order is floor(t)+2. For all tensor orders simultaneously, positive permutation traces of this form force t to be a nonnegative integer. Negative t already fails order one.

In particular t=7/2 passes every fourth-order permutation-trace check but fails at fifth order with alternating trace -7/256. The formal dimension t=57/2 first fails at order30. An ordinary dimension57 instead has first vanishing exterior power58 and never a negative one.

## 3. What +2 and doubling do, and what they do not do

For actual vector spaces, V -> V direct_sum C^2 and V -> V tensor C^2 implement dimension maps t->t+2 and t->2t. They obey B A=A^2 B on dimensions, not commutation. For a noninteger t whose doubled value is still noninteger, the first-negative horizon changes from floor(t)+2 to floor(2t)+2. Adding two moves the horizon exactly two steps. If t is half-integral, doubling lands on an integer and eliminates this particular scalar obstruction.

This does NOT repair the original nonpositive category: the initial object's negative projector remains if that object remains in the theory. Nor is it an equivalence to the tournament constructions merely because sizes agree. A functor and its effect on the distinguished operators would be needed.

## 4. Rationality is not the positivity boundary

The criterion depends on crossing an integer, not on whether t is rational, algebraic irrational, or transcendental. All t in (1,2), including sqrt(2) and pi/2, have the same first failing order3. These real values have different arithmetic types but identical positivity horizon.

Deligne interpolation supplies genuine symmetric tensor categories with complex dimension parameter, without asserting a positive Hilbert-space realization. See Etingof--Gelaki--Nikshych--Ostrik, *Tensor Categories*, sections9.11--9.12, especially Theorem9.12.10 and the loop-evaluation construction:
https://math.mit.edu/~etingof/egnobookfinal.pdf

Our result is a statement about this particular positive permutation-trace family. It does not prohibit noninteger quantum dimensions in braided unitary categories: their braiding need not satisfy the symmetric-group relations, and the trace of their braid operators need not be t^{number of cycles}. It does not identify categorical dimension, central charge, Lie-algebra dimension, or conformal weight.

## 5. Connection contract and verification

Map: symmetric tensor-factor permutations -> cycle trace -> signed projector trace. Preserved: character, positivity at each fixed tensor order, exterior dimension. Lost: all multiplication/bracket/associator data and an actual underlying Hilbert realization when t is formal. Sidecar: the complete family of tensor-order Gram tests, not just scalar dimension or a finite initial segment. Decisive hostile: t=7/2 at order5.

`check_dimension_positivity.py` uses exact rational arithmetic, all partitions through order12 on a quarter-integer grid, direct signed permutation sums through order7, and explicit higher-order witnesses. The output `dimension_positivity.json` records460 checks. Explicit exceptions remain active under python -O; normal and optimized outputs match. This finite check is not the proof of the all-order theorem above.

The potentially reusable principle is a **positivity horizon**: interpolation may preserve all low-order observations while an explicit higher-order projector detects impossibility. This is analogous to the earlier refinement issue, but not a transferred Collatz theorem or an obstruction to every noninteger category.
