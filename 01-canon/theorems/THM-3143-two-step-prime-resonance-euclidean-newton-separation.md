---
id: THM-3143
title: "Two-step prime resonance Euclidean-Newton separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  odd prime p, no exact quadratic has three zero factorial moments beginning
  at r=p.  At resonance d=p+2, the two original p-adic polygons share the
  Frobenius slope 2/p.  One exact Euclidean cancellation replaces the second
  polynomial by a remainder with slopes 2/(p+1) and 2/(p-1), disjoint from
  2/p.  Thus the overlap is removable and the resonant pair is coprime.
audit: >
  An independent hostile audit rederived both original polygons and their
  normalized Frobenius faces, the exact leading-ratio cancellation, every
  remainder divisibility threshold, the midpoint residue, the closed leading
  identity and valuation, p=3, the common-factor implication, and the final
  three-composite residual scope.  Fresh normal, optimized, and stored
  transcripts and LF-normalized hashes agree exactly.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3138-left-factorial-resonance-newton-slope-separation
script: 04-computation/factorial_two_step_prime_resonance_euclidean_newton_thm3143.py
output: 05-knowledge/results/factorial_two_step_prime_resonance_euclidean_newton_thm3143.out
script_sha256: 4386065bc5a7308c3e77bab0c672f1101ac95eb155a7860e9586586bf51a467c
output_sha256: 5744b86567ba6d96bdd759e00de31e78c63fbc4fb4e1f492f6c0a9ccbed82fe4
hash_basis: LF-normalized bytes
---

# THM-3143 -- two-step prime resonance Euclidean-Newton separation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(t^k)=k!,                    q(t)=a+bt+ct^2,                 (1)
```

with `abc!=0`.  For every odd prime `p`, put

```text
r=p,                         d=r+2=p+2.                      (2)
```

Then the three consecutive moments

```text
L(q^r), L(q^(r+1)), L(q^(r+2))                               (3)
```

cannot all vanish.

The direct Newton polygons overlap at this resonance.  The proof removes
their common initial face by one exact Euclidean step and then separates the
replacement polygons.

## 1. Resonant integral pair

THM-3124 says that `(3)` forces `b/a=-1/d`.  After dividing by `a`, writing
`u=c/a`, and making the invertible change `v=d u`, define

```text
A_n(v)=d^n L((1-t/d+u t^2)^n)
      =L((d-t+v t^2)^n) in Z[v].                             (4)
```

It is enough to prove that

```text
A=A_p,                         B=A_(p+1)                     (5)
```

are coprime.  If `a_(n,j)=[v^j]A_n`, then

```text
a_(n,j)=binom(n,j) sum_(ell=0)^(n-j)
  binom(n-j,ell)d^(n-j-ell)(-1)^ell(2j+ell)!.                (6)
```

Throughout, `d=p+2==2 (mod p)`.

## 2. The genuine Frobenius collision

For `A=A_p`, the constant coefficient is `2 (mod p)`.  Every coefficient
with `0<j<p` has valuation at least one from `binom(p,j)`; when
`j>=(p+1)/2`, its factorial supplies a second `p`.  Finally

```text
[v^p]A=(2p)!,                         v_p((2p)!)=2.           (7)
```

All intermediate points lie strictly above the endpoint line, so

```text
NP_p(A):(0,0)--(p,2),                    sole slope 2/p.      (8)
```

The polygon of `B` has vertices

```text
(0,0), (1,0), (p+1,2),                  slopes 0 and 2/p.    (9)
```

This is not an accidental numerical overlap.  Let `pi^p=p` in a ramified
extension and scale `v=pi^(-2)z`.  On their slope-`2/p` faces, the normalized
endpoint residues are

```text
A: 2+2z^p=2(1+z)^p,
B: 4z+4z^(p+1)=4z(1+z)^p.                                (10)
```

Thus both faces carry the same multiplicity-`p` residue root `z=-1`.
Plain slope comparison has genuinely lost the next-order datum.

## 3. Exact Euclidean cancellation

The leading-coefficient ratio is the integer

```text
C=(2p+2)(2p+1).                                             (11)
```

Set

```text
R(v)=B(v)-C v A(v).                                         (12)
```

The degree-`p+1` term cancels exactly.  Equation `(10)` predicts that the
whole common face cancels; the coefficient valuations show what replaces it.

First, `R(0)=B(0)==2 (mod p)`.  Every positive-degree coefficient of `R` is
divisible by `p`.  For `j=1`, this is the explicit congruence

```text
[v]B ==4 == C A(0)                              (mod p);     (13)
```

for `2<=j<p`, both terms in `(12)` have a binomial factor `p`, and the case
`j=p` is even more divisible.

Put

```text
m=(p-1)/2,                         k=m+1=(p+1)/2.             (14)
```

At index `m` of `A`, only the `ell=0` term survives after division by `p`:

```text
a_(p,m)/p == (-1)^m 2^(m+1)/m !=0                 (mod p).   (15)
```

At index `k` of `B`, the outer binomial and the factorial together supply
at least `p^2`.  Since `C==2 (mod p)`, `(12)` and `(15)` give

```text
[v^k]R/p == (-1)^(m+1) 2^(m+2)/m !=0              (mod p).  (16)
```

Hence the coefficient at `k` has valuation exactly one.

For every `k<j<p`, both terms of `[v^j]R` have valuation at least two: the
relevant binomial contributes one `p` and the factorial contributes another.
At the leading surviving index, direct two-term simplification gives the
closed identity

```text
[v^p]R
 =(1-p^2)(2p)!-C p(3-p)(2p-2)!
 =-4p(p+1)(p+2)(2p-2)!,                                  (17)

v_p([v^p]R)=2.                                             (18)
```

Therefore every intermediate point is strictly above one of the two endpoint
segments, and

```text
NP_p(R):(0,0)--(k,1)--(p,2),
slopes             1/k       1/(p-k)
                 =2/(p+1),  =2/(p-1).                       (19)
```

Neither slope in `(19)` equals `2/p`.

## 4. Coprimality and the factorial window

If `A` and `B` had a nonconstant common factor over `Q`, it would also divide
the exact combination `R=B-CvA`.  Base change to `Q_p` preserves a
nonconstant common factor, which would contribute a common lower Newton
slope to `A` and `R`.  Equations `(8)` and `(19)` forbid this.  Thus `A` and
`B` are coprime over `Q`, so the resonant pair has no common complex root.
This contradicts `(3)` and proves the theorem.  QED.

The constants of both `A` and `R` are units, so zero roots or infinite slopes
do not create an exception.

## 5. Exact controls

The companion reconstructs `(6)` over the integers for
`p=3,5,7,11,101,211`, checks the complete polygons of `A`, `B`, and `R`, and
verifies the exact cancellation and leading identity `(17)`.  An independent
repeated-product bivariate constructor checks both original polynomials for
`p<=11`.  It also verifies the endpoint residues in `(10)` and scans the
closed nonzero residues `(16),(18)` over all `1,228` odd primes below
`10,000`.  The finite scan is a control; the proof is symbolic for every odd
prime.

Run

```text
python3 04-computation/factorial_two_step_prime_resonance_euclidean_newton_thm3143.py
python3 -O 04-computation/factorial_two_step_prime_resonance_euclidean_newton_thm3143.py
```

and compare byte-for-byte with the declared output.

## 6. Connection contract and scope

The source object is THM-3124's resonant pair `(5)`.  Direct passage to the
`p`-adic associated graded destroys the distinction between the two
Frobenius faces `(10)`.  The repair map is the exact Euclidean syzygy
`(A,B)->(A,B-CvA)`, which preserves the common-factor predicate while
exposing the next valuation layer.  The load-bearing sidecars are the
midpoint residue `(16)` and exact leading identity `(17)`.

Combining THM-3131, THM-3138, this theorem, and THM-3124, any still-open bad
exact-quadratic window must satisfy

```text
r>=201,     d=r+2 composite,     d-1 composite,     d-2 composite. (20)
```

Thus its resonance lies after a prime gap of length at least three.  The
mechanism suggests an iterated Euclidean-Newton resolution for longer gaps,
but no such general induction is asserted here.

This is an exact `{0,1,2}` / `SFC(1)` factorial-moment theorem.  It does not
settle arbitrary-support `SFC(3)`, `GMC(2)`, `NC(2)`, or `LRC(14)`.

## 7. Independent hostile audit

The audit independently reconstructed the polygons of `A` and `B`, including
the common normalized faces `(10)`, before checking the Euclidean repair.  It
rederived `R(0)`, the cancellation at degree one, all coefficient valuation
thresholds, the nonzero midpoint residue `(16)`, and the exact leading
identity `(17)`.  It separately checked `p=3`, the common-factor/base-change
argument, and the combined `d,d-1,d-2` residual scope.  Fresh normal and
optimized runs byte-match the stored transcript and declared hashes.

**End of proof.**
