---
id: THM-2412
title: "Delta exponential and central Gregory--Newton layer split"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Differentiation on powers and forward
  difference on falling factorials are conjugate lowering operators under
  the Stirling transform. Their eigenfunctions are e^(lambda x) and
  (1+lambda)^n; hence 2^n, not an arbitrary power, is the eigenvalue-one
  unit for Delta. After the displayed one-index shift, the two supplied
  Pascal-half sequences give the lower weak count and, by complementing
  subsets, its upper strict complement; the central binomial layer is their
  exact sidecar, and
  Catalan convolution carries weak half to 4^n to strict half. This is an
  operator and labelled-operation theorem, not a causal derivation of the
  number two from tournaments.
source: codex-2026-07-26-delta-exponential
depends_on: []
related:
  - THM-361-product-sum-defect-normal-form
  - THM-362-natural-operation-graph-shadows
  - THM-710-factorial-moment-eigen-transfer
script: 04-computation/delta_exponential_central_newton_split_thm2412.py
output: 05-knowledge/results/delta_exponential_central_newton_split_thm2412.out
script_sha256: 85c3d1229480c6024a7dd06da495c4c78c58da20a08b5991bb0adc928443d080
output_sha256: 781aa07051597263073dfe7d2adc26680ed8b42c8725df96e43230908399dc67
hash_basis: working-tree bytes (LF)
---

# THM-2412 -- delta exponentials and the central Newton split

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

There are three exact objects here:

```text
continuous coordinate:     powers x^k, operator D;
discrete coordinate:       falling factorials x^(underline k), operator Delta;
binary layer coordinate:   subsets, graded by cardinality.
```

The first two are conjugate through the Stirling transform. The third gives
the coefficients of the discrete exponential. Keeping the middle layer
rather than quotienting it away identifies the two apparently unrelated
sequences in the prompt.

## 1. Maclaurin and Gregory--Newton are the same lowering scheme

Put

```text
D f(x)=f'(x),                 Delta f(x)=f(x+1)-f(x),

x^(underline k)=x(x-1)...(x-k+1).
```

Then

```text
D x^k = k x^(k-1),

Delta x^(underline k) = k x^(underline (k-1)).       (1)
```

The second identity follows by factoring:

```text
(x+1)^(underline k)-x^(underline k)

=x^(underline (k-1))((x+1)-(x-k+1))

=k x^(underline (k-1)).
```

Consequently every polynomial `f` of degree at most `d` has the two exact
expansions

```text
f(x)=sum_(k=0)^d f^(k)(0) x^k/k!,                   (2)

f(x)=sum_(k=0)^d Delta^k f(0) x^(underline k)/k!
    =sum_(k=0)^d Delta^k f(0) binom(x,k).            (3)
```

Equation (2) is Maclaurin in the derivative-lowering basis. Equation (3) is
Gregory--Newton in the difference-lowering basis. The proof of (3) is the
same triangular interpolation proof as (2): both sides have the same value
at zero, and applying `Delta` reduces the degree and shifts the coefficient
list by one.

If `S(j,k)` and `s(k,j)` are respectively Stirling numbers of the second kind
and signed first kind, then

```text
x^j = sum_k S(j,k)x^(underline k),

x^(underline k)=sum_j s(k,j)x^j.                    (4)
```

Thus the coefficient dictionaries are

```text
Delta^k f(0)/k!
 =sum_(j>=k) S(j,k) f^(j)(0)/j!,                    (5)

f^(j)(0)/j!
 =sum_(k>=j) s(k,j) Delta^k f(0)/k!.                (6)
```

No limiting argument is present for polynomials.

The dictionary has an exact lattice-spacing parameter. For `h!=0`, put

```text
D_h f(x)=(f(x+h)-f(x))/h,

x^(underline k,h)=product_(j=0)^(k-1)(x-jh),

U_h(x^k)=x^(underline k,h).                         (6a)
```

The images of `U_h` are monic of degree `k`, so `U_h` is a linear
automorphism of the polynomial vector space (not an algebra homomorphism),
and direct cancellation gives

```text
D_h x^(underline k,h)=k x^(underline (k-1),h),

D_h U_h=U_h D.                                      (6b)
```

Thus the two calculi are literally intertwined, not merely asymptotic.

Let `E f(x)=f(x+1)`. On polynomials the operator series terminate, so

```text
E=exp(D),             Delta=exp(D)-I,

D=log(I+Delta).                                      (7)
```

This is the exact continuous/discrete bridge. For analytic functions, (7)
requires the usual convergence/domain hypotheses; the polynomial identity
does not grant an unrestricted Newton series for every entire function.

## 2. Exponential eigenfunctions and the special role of two

For `lambda in C` and `n in Z_(>=0)`,

```text
D exp(lambda x)=lambda exp(lambda x),                (8)

Delta (1+lambda)^n=lambda(1+lambda)^n.               (9)
```

At `lambda=-1`, use the standard lattice convention `0^0=1`; equivalently
check the `n=0` endpoint separately.

The second line is immediate:

```text
(1+lambda)^(n+1)-(1+lambda)^n
=lambda(1+lambda)^n.
```

Therefore the eigenvalue-one units are

```text
D e^x=e^x,                    Delta 2^n=2^n.         (10)
```

The `2` in (10) is exactly `1+lambda` at `lambda=1`: one copy from the
identity part of the shift and one from `Delta`. It is not an additional
number-theoretic assumption.

For every nonnegative integer `n`, Gregory--Newton is finite and gives

```text
(1+lambda)^n
 =sum_(k=0)^n lambda^k binom(n,k)
 =sum_(k=0)^n lambda^k n^(underline k)/k!.           (11)
```

In particular,

```text
2^n=sum_(k=0)^n n^(underline k)/k!,                  (12)

4^n=sum_(k=0)^n 3^k binom(n,k).                      (13)
```

More generally, the unit-eigenvalue solution of `D_h f=f`, `f(0)=1`, on
`h Z_(>=0)` is

```text
f(Nh)=(1+h)^N.                                      (13a)
```

At a fixed endpoint `X=Nh`,

```text
(1+lambda X/N)^N -> exp(lambda X),

X^(underline k,X/N) -> X^k                         (13b)
```

as `N->infinity` for every fixed `k`. The operator, its basic basis, and
its eigenfunction therefore pass to the continuous limit together.

For generalized `x`, the binomial series

```text
(1+z)^x=sum_(k>=0) binom(x,k)z^k
```

is analytic for `|z|<1` after choosing the compatible logarithm. The value
`z=1` in (12) is deliberately asserted only at nonnegative integers, where
the series terminates.

## 3. The central Pascal sidecar

Define, for `n>=0`,

```text
A_n=sum_(k=0)^n binom(2n,k),

B_n=sum_(k=0)^n binom(2n+2,k).                       (14)
```

Their initial values are

```text
A: 1,3,11,42,163,638,...,

B: 1,5,22,93,386,... .                              (15)
```

The reflection `k -> 2n-k` pairs every noncentral layer of the row of
Pascal's triangle. Hence

```text
A_n=(4^n+binom(2n,n))/2,                             (16)

B_n=(4^(n+1)-binom(2n+2,n+1))/2.                    (17)
```

Equivalently, for every `n>=0`,

```text
A_(n+1)+B_n=4^(n+1),                                (18)

A_(n+1)-B_n=binom(2n+2,n+1).                        (19)
```

Thus `A_(n+1)` counts the lower weak family `|S|<=n+1`, while `B_n`
counts the lower strict family `|S|<=n`. Those two literal families overlap.
Complementation bijects the latter with the upper strict family
`|S|>=n+2`, which is the actual set complement of the lower weak family.
It is in this precise, count-preserving sense that (18) is a complementary
split. The central binomial coefficient is the exact tie layer. Forgetting
which side and which lower/upper realization was chosen preserves (18) and
destroys (19).

This supplies a precise tournament analogy. A labelled tournament on `m`
vertices is a binary choice on each of the

```text
T_m=binom(m,2)
```

unordered pairs, so there are `2^(T_m)` tournaments. Fixing a reference
orientation identifies tournaments with subsets of those `T_m` pairs.
The power two comes from the binary orientation coordinate; the triangular
number counts how many such coordinates there are. This is a combinatorial
realization of (11), not a proof of the analytic eigenfunction law.

It also isolates the tie guardrail: allowing a third, unoriented state would
replace the local binary alphabet by a ternary one. A weak/strict comparison
cannot silently delete the central layer and still be called an equivalence.

There is a second exact tournament realization. Given a reference
orientation with signs `a_(ij) in {+1,-1}`, vertex switching by
`x_i in {+1,-1}` sends

```text
a_(ij) -> a_(ij)x_i x_j,

Q_a(x)=sum_(i<j)a_(ij)x_i x_j.                      (19a)
```

Thus `Q_a(x)` is the signed edge bias of the switched representative and
`Q_a(-x)=Q_a(x)` records the global-sign gauge. This identifies the
plus-minus quadratic form as switching energy; it does not prove an
asymptotic min--max limit, because cycle products retain information that
the scalar bias destroys.

## 4. Catalan convolution is the self-similar ladder

Let

```text
C_n=binom(2n,n)/(n+1),             P_n=4^n
```

and write `*` for ordinary sequence convolution. Then

```text
C*A=P,                       C*P=B,

B=C*C*A.                                             (20)
```

Proof: set `s=sqrt(1-4z)`. The four generating functions are

```text
C(z)=2/(1+s),

A(z)=(1/s^2+1/s)/2=(1+s)/(2s^2),

P(z)=1/s^2,

B(z)=(1/s^2-1/s)/(2z)=C(z)P(z).                     (21)
```

Multiplying the first two functions gives `P(z)`, and multiplying by
`C(z)` again gives `B(z)`. This is an all-coefficient identity, not a
finite-prefix pattern.

The same boundary mechanism has the pointwise recurrence form. Put

```text
M_n=(4^n-binom(2n,n))/2,                    M_0=0.
```

Then `M_(n+1)=B_n`, and for `n>=1`,

```text
A_n=4A_(n-1)-C_(n-1),

M_n=4M_(n-1)+C_(n-1).                              (21a)
```

The Catalan number is therefore the exact one-step leakage through the
reflection-fixed middle layer: it leaves the weak half and enters the
strict half with the opposite sign.

The ladder

```text
weak Pascal half --Catalan convolution--> full binary cube
                 --Catalan convolution--> strict Pascal half
```

is self-similar but asymmetric: the middle layer in (19) records the
direction of travel.

## 5. A labelled summand--multiplicand collision

Equations (18) and the geometric recursion for `P_n` meet at one target:

```text
A_(n+1)+B_n = 4^(n+1) = 4*4^n.                      (22)
```

In the labelled operation-cospan language of THM-362, (22) has

```text
additive parents:       A_(n+1), B_n;
multiplicative parents: 4, 4^n;
common target:          4^(n+1);
lost quotient datum:    which operation and which Pascal half;
needed sidecar:         binom(2n+2,n+1).
```

For `n>=2`, the multiplicative parents are distinct nonunits, so neither a
unit convention nor a diagonal convention creates the collision.

The unlabelled additive shadow alone is only the order relation `x<z`, and
the multiplicative shadow alone is divisibility. They do not retain (22).

## 6. Equality and failure boundaries

1. Ordinary powers are not the exact difference-lowering basis:

   ```text
   Delta x^2=2x+1,               not 2x.
   ```

2. The eigenvalue-one discrete exponential is specifically `2^n`.
   For example, `Delta 3^n=2*3^n`.
3. The sixth term in (15) is `638`. The five-term prefix
   `1,3,11,42,163` does not determine a sequence; a continuation by `639`
   belongs to a different object and fails (16).
4. The central binomial layer is not lower-order bookkeeping: it is exactly
   the difference between the two halves.
5. Catalan convolution is a statement about the labelled coefficient
   sequences. It does not turn an LRC toothpick word or an arbitrary
   tournament quotient into the same object.
6. Switching energy in (19a) is a scalar quotient. It does not retain the
   switching-invariant cycle products.

## 7. Exact companion

The dependency-free companion:

- verifies `228` rational falling-factorial lowering identities;
- verifies both Stirling coefficient transforms on a nontrivial rational
  polynomial;
- checks `E=exp(D)` and `D=log(I+Delta)` at `22` rational points;
- verifies `231` finite exponential identities and the `2^n`, `4^n`
  specializations;
- verifies the tournament binary-coordinate count through `12` vertices;
- checks (16)--(20) through `n=50`; and
- retains monomial, nonunit-eigenvalue, tie-deletion, and sixth-term hostiles.

Run

```bash
python3 04-computation/delta_exponential_central_newton_split_thm2412.py
python3 -O 04-computation/delta_exponential_central_newton_split_thm2412.py
```

Both outputs must byte-match the stored transcript. Every executable check
raises explicitly under optimized Python.

## 8. Independent hostile audit

An independent audit rederived the lowering identities, both Stirling
transforms, the scaled-`h` intertwiner, the terminating exponential, the
Pascal-half split, and both Catalan convolutions. It separately checked the
reflection-fixed tie layer and the tournament arc-slot count, including the
non-homomorphism, ordinary-monomial, nonunit-eigenvalue, and sixth-term
hostiles. Normal and optimized transcripts byte-match the stored output, and
the recorded hashes match the audited files. No claim about nonintegral
Newton convergence or switching asymptotics is used.
