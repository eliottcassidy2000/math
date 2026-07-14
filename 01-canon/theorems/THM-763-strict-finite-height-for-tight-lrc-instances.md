---
id: THM-763
title: Strict finite height for tight lonely-runner instances
status: PROVED, conditional only on the lower-dimensional LRC citation
source: codex-2026-07-14-S3
depends_on:
  - LRCUpTo13
  - THM-759
related:
  - THM-668-pair-sum-ruler-witness-structure
  - HYP-6775
  - HYP-6800
  - HYP-6820
external:
  - Malikiosis--Santos--Schymura, Forum of Mathematics Sigma 13 (2025), Theorem A and Section 3.2
  - Giri--Kravitz, Mathematical Proceedings of the Cambridge Philosophical Society 180 (2026), Section 7
---

# THM-763 — Strict finite height for tight lonely-runner instances

## Statement

For `n >= 2` and a positive integer tuple `v = (v_1,...,v_n)`, write

```text
M(v) = max_t min_i ||v_i t||.
```

Assume the Lonely Runner Conjecture for `n-1` nonzero speeds.  If `v` is
primitive (`gcd(v_1,...,v_n)=1`) and

```text
                         n-1
sum_i v_i > binom(n+1,2)   ,
```

then the conclusion is strict:

```text
M(v) > 1/(n+1).
```

Consequently, every primitive tight `n`-speed tuple satisfies

```text
                         n-1
sum_i v_i <= binom(n+1,2)   .                 (1)
```

For the twelve-speed problem this gives the explicit uniform bound

```text
sum_i a_i <= 78^11 = 650190514836423555072.  (2)
```

If `1 <= a_1 < ... < a_12`, then the other eleven entries have sum at
least `1+...+11=66`, so

```text
a_12 <= 78^11 - 66 = 650190514836423555006.  (3)
```

Thus the primitive tight twelve-speed locus, including the putative
sporadic branch, is finite.

## Proof

The argument is a strict-interior refinement of the proof of Theorem A in
Malikiosis--Santos--Schymura (MSS).  It is useful to spell out the strictness,
because the theorem as normally quoted only asserts `M(v) >= 1/(n+1)`.

Set

```text
d = n-1,             C = binom(n+1,2) = n(n+1)/2.
```

Let `Z = Z_v` be the `d`-dimensional lonely-runner zonotope associated to
the primitive tuple `v`, and let `x` be its centre.  The standard volume
identity is

```text
vol(Z) = sum_i v_i.
```

The difference body is `Z-Z = 2(Z-x)`.  If `sum_i v_i > C^d`, then

```text
vol((Z-Z)/C) = 2^d vol(Z)/C^d > 2^d.
```

Minkowski's first theorem, used contrapositively with the strict volume
inequality, gives a nonzero lattice point in the interior of `(Z-Z)/C`.
Equivalently, for the first successive minimum,

```text
lambda_1(Z-Z) < 1/C.                                  (4)
```

Now use the symmetric body from MSS Section 3.2,

```text
K = ((n-1)/(n+1)) (Z-x).
```

Since `K-K = ((n-1)/(n+1))(Z-Z)` and `K` is centrally symmetric,

```text
lambda_1(K)
  = 2(n+1)/(n-1) lambda_1(Z-Z)
  < 4/(n(n-1)).                                        (5)
```

Choose a primitive lattice vector `w` attaining `lambda_1(K)` and project
along `w`.  MSS Theorem 3.1 shows that the projected zonotope contains a
lonely-runner zonotope of dimension `n-2` with the same centre.  The assumed
LRC for `n-1` speeds therefore supplies the lower-dimensional point `b` used
in their proof.  If `T` is the projection, `T(x)=y`, and `T(w)=0`, then

```text
b in ((n-2)/n)(T(Z)-y) intersect (y + Z^(n-2)).
```

Choose a lift

```text
a in ((n-2)/n)(Z-x),       T(a)=b.
```

This is also where the affine lattice on the fiber enters.  Because
`b-y` belongs to `Z^(n-2)=T(Z^(n-1))`, the line

```text
a + Rw = T^(-1)(b)
```

contains infinitely many points of the coset `x+Z^(n-1)`.  More precisely,
for a suitable `z in Z^(n-1)` its intersection with that coset is

```text
(a+Rw) intersect (x+Z^(n-1)) = (x+z) + Zw,
```

since `w` is primitive.  Consecutive points therefore have lattice spacing
one.  This is the paragraph immediately before MSS Lemma 3.6; primitivity of
`w` alone would not imply that an arbitrary parallel line meets the required
coset.

There are now the same two cases as in MSS Lemma 3.6.

- If `a in Rw`, then the fiber is `Rw` and its chord through `K` has lattice
  length `2/lambda_1(K)`.  By (5),

  ```text
  2/lambda_1(K) > n(n-1)/2 >= 1.                      (6)
  ```

- If `a notin Rw`, the triangle construction in the proof of MSS Lemma 3.6
  puts the segment

  ```text
  [a - 2w/(n(n-1)lambda_1(K)),
   a + 2w/(n(n-1)lambda_1(K))]
  ```

  inside `K`.  Its lattice length is

  ```text
  4/(n(n-1)lambda_1(K)) > 1.                          (7)
  ```

In either case, a segment of lattice length strictly greater than one on this
fiber contains a point of `x+Z^(n-1)` in its relative interior.  That point is
in `int(K)`.  In the first case this follows because every non-endpoint of a
central chord of the full-dimensional symmetric body `K` is interior.  In the
second case, first note that

```text
(n-2)/n < (n-1)/(n+1),
```

so `a` is already an interior point of `K`; every non-endpoint of either half
of the displayed segment is then a proper convex combination of `a` and an
endpoint in `K`.

Finally, the lonely-runner-zonotope equivalence sends an intersection with
`K` to a `1/(n+1)` witness.  An intersection with `int(K)` is stronger: it
lies in `c(Z-x)` for some `0 <= c < (n-1)/(n+1)`.  By the parameterized
zonotope equivalence in MSS Remark 1.6, writing `c=1-2 eta` gives
`eta>1/(n+1)` and hence a time at which every speed has distance at least
`eta`.  Therefore `M(v)>1/(n+1)`, as claimed.  Taking the contrapositive
proves (1), and `n=12` gives (2)--(3).  ∎

## Exact finite decision for the twelve-speed sporadic branch

The height theorem makes the equality problem genuinely finite, though not
remotely feasible by a naive enumeration.  The pair-sum ruler theorem makes
each individual test exact.  For a candidate `A={a_1,...,a_12}` set

```text
q_ij = a_i + a_j.
```

Then

```text
M(A) = max_{i <= j} max_{0 < p < q_ij}
       min_k dist(p a_k mod q_ij, {0,q_ij}) / q_ij.   (8)
```

Thus (2), primitivity, distinctness, and the finite integer calculation (8)
form a uniform decision procedure for `M(A)=1/13`.  At a tight maximum below
`1/2`, two distinct active runners may be chosen, so `i<j`.  If their ruler is
`q=a_i+a_j` and the minimum residue distance is `m`, tightness says

```text
m/q = 1/13,
```

and hence

```text
13 divides q = a_i+a_j.                               (9)
```

For the max-peel `P=A\{a_12}`, the sporadic condition is `M(P)>1/12`.
THM-759 adds the independent branch-specific restriction

```text
a_12 <= a_11/(13 M(P)-1) < 12 a_11.                  (10)
```

Equations (2), (8), (9), and (10) are therefore a quantifier-honest finite
reduction of uniform sporadic-branch emptiness.

## What this does not prove

This theorem does **not** prove that the sporadic branch is empty.  The bound
`78^11` is about `6.5 * 10^20`; the existing exact censuses to heights `16`,
`24`, or `48`, and the bounded residue-lift searches, do not exhaust it.  In
particular:

- the MSS strictness argument supplies finite height, not AP rigidity;
- divisibility by `13` of one active pair sum does not force the twelve
  integer speeds to be `{1,...,12}`;
- the pair-sum formula is an exact verifier after a candidate is supplied,
  not a practical enumeration of all candidates under (2).

So the honest advance is from an apparently unbounded equality problem to an
explicit finite one.  A structural compression theorem is still required to
turn that finite decision into a proof of branch emptiness.

## Sources and scope

- R. D. Malikiosis, F. Santos, and M. Schymura,
  *Linearly-exponential checking is enough for the lonely runner conjecture
  and some of its variants*, Forum of Mathematics, Sigma **13** (2025),
  e164, DOI `10.1017/fms.2025.10107`: equation (1.3), Remark 1.6,
  Minkowski's theorem as Theorem 2.4, Theorem 3.1, equation (3.7), and the
  paragraph immediately before and proof of Lemma 3.6.
- V. Giri and N. Kravitz, *The structure of Lonely Runner spectra*,
  Math. Proc. Cambridge Philos. Soc. **180** (2026), 343--361, DOI
  `10.1017/S0305004125101497`, Section 7.  Their
  independent spectrum argument also states that lower-dimensional LRC makes
  the tight locus finite; the proof above gives the sharper explicit
  sum-of-speeds bound needed here.
- `01-canon/theorems/THM-759-tight-instance-ratio-bound.md` supplies (10).
- `01-canon/theorems/THM-668-pair-sum-ruler-witness-structure.md` supplies
  (8), and (9) is its immediate equality specialization.
