---
id: THM-2412
title: "Delta exponential and central Gregory--Newton layer split"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The scaled
  falling-factorial umbral map U_h(x^k)=
  x(x-h)...(x-(k-1)h) exactly intertwines the continuous derivative
  with the forward difference quotient. Its unit-eigenvalue lattice
  solution is (1+h)^(x/h), so h=1 gives 2^N while the fixed-endpoint
  h->0 limit gives e^x. On the central Bernoulli-triangle cut,
  A032443=(4^n+C(2n,n))/2 and shifted A000346=
  (4^n-C(2n,n))/2 partition 4^n; the central binomial coefficient is
  exactly the reflection-fixed tie layer and Catalan numbers are the
  signed one-step leakage. For tournaments the exponent is the number
  C(v,2) of independent arc slots, not the number of vertices and not a
  causal source of the base two. No nonintegral Newton-series
  extension is used, and no plus-minus quadratic-form limit is claimed.
source: codex-2026-07-26-delta-exponential
depends_on: []
related:
  - THM-438-paley-cluster-integrals-are-catalan
  - THM-710-factorial-moment-eigen-transfer
  - THM-1415-switching-is-the-canonical-star-quotient
  - THM-1430-graph-switching-is-exactly-E-n
  - THM-1470-even-tournaments-are-the-tournament-two-graph-theorem
external: >
  OEIS A032443 and A000346; MathOverflow question 413935,
  "Min max of a quadratic form of plus-minus ones".
script: 04-computation/delta_exponential_central_newton_split_thm2412.py
output: 05-knowledge/results/delta_exponential_central_newton_split_thm2412.out
script_sha256: b4c4decc432ca3b4948b28e03f6f0753ec7e69186bdd2285d2c1792495c7b088
output_sha256: 7995b067b71bc152a44393c2cd5efa596c2f5f6f3a16f5d0b8ea78fcaaa3b34a
hash_basis: working-tree bytes (LF)
cite_by_filename: true
---

# THM-2412 -- delta exponentials and the central Newton split

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The Maclaurin and Gregory--Newton expansions are not merely analogous
lists of coefficients. They are the same lowering-operator calculus in
two polynomial bases:

```text
continuous:  D(x^k)=k x^(k-1),
discrete:    D_h(x^(falling k,h))=k x^(falling k-1,h).               (1)
```

The associated eigenfunctions are `e^(lambda x)` and
`(1+h lambda)^(x/h)`. The familiar base `2` is the specialization
`h=lambda=1`, while `e` is the fixed-endpoint limit as the lattice
spacing tends to zero.

The three displayed sequences in the motivating prompt form a second
exact packet. They are the total and the two central half-spaces of a
Boolean cube:

```text
1,4,16,64,256,...          =4^n,
1,3,11,42,163,...          =(4^n+C(2n,n))/2,
1,5,22,93,386,...          =(4^(n+1)-C(2n+2,n+1))/2.                 (2)
```

The apparent startup asymmetry is precisely a fixed boundary layer.

## 1. The scaled Gregory--Newton intertwiner

Fix `h!=0` in a characteristic-zero field and define

```text
D_h f(x)=(f(x+h)-f(x))/h,                                            (3)

x^(falling k,h)
  =product_(j=0)^(k-1)(x-jh),

x^(falling 0,h)=1.                                                   (4)
```

Direct cancellation gives

```text
D_h x^(falling k,h)=k x^(falling k-1,h).                             (5)
```

For `k>=1`, one proof factors the numerator (with the remaining product
empty when `k=1`):

```text
(x+h)^(falling k,h)-x^(falling k,h)
 =[(x+h)-(x-(k-1)h)]
   x(x-h)...(x-(k-2)h)
 =kh x^(falling k-1,h).                                              (6)
```

Define the degree-preserving linear map

```text
U_h:K[x]->K[x],
U_h(x^k)=x^(falling k,h).                                            (7)
```

Every image in (7) is monic of degree `k`, so `U_h` is an isomorphism.
Equations (1) and (5) are the operator identity

```text
D_h U_h=U_h D.                                                       (8)
```

Thus powers and scaled falling factorials are not competing
approximations. They are the basic sequences of the derivative and
forward-difference delta operators.

## 2. The discrete exponential and the continuous limit

On the lattice `h Z_(>=0)`, solve

```text
D_h f=lambda f,                         f(0)=1.                       (9)
```

The recurrence is

```text
f((N+1)h)=(1+h lambda)f(Nh),                                         (10)
```

so its unique solution is

```text
f(Nh)=(1+h lambda)^N.                                                (11)
```

The terminating Gregory--Newton expansion is exactly

```text
(1+h lambda)^N
 =sum_(k=0)^N lambda^k (Nh)^(falling k,h)/k!
 =sum_(k=0)^N C(N,k)(h lambda)^k.                                   (12)
```

At `h=lambda=1`,

```text
2^N
 =sum_(k=0)^N N^(falling k)/k!
 =sum_(k=0)^N C(N,k).                                                (13)
```

More generally, setting `h=1,lambda=q-1` gives

```text
q^N=sum_(k=0)^N (q-1)^k N^(falling k)/k!.                           (14)
```

So base two is not forced by discreteness alone. It is
`1+h lambda` at unit step and unit eigenvalue.

For fixed `X,lambda in C` with `X!=0`, take `h=X/N` and let the positive
integer `N` tend to infinity in (11)--(12). Then

```text
(1+lambda X/N)^N
 =sum_(k=0)^N lambda^k X^(falling k,X/N)/k!
 -> e^(lambda X)
 =sum_(k>=0) lambda^k X^k/k!.                                       (15)
```

For each fixed `k`, the scaled falling factorial in (15) tends to
`X^k`, and the left side is the standard exponential limit. Equation
(15) is the exact discrete/continuous bridge: the operator, basis, and
eigenfunction converge together. The case `X=0` is the separate
constant identity `1=e^0`.

## 3. Bernoulli's triangle and the two central sequences

Define the Bernoulli-triangle entry

```text
B(r,s)=sum_(k=0)^s C(r,k),              0<=s<=r,
B(r,-1)=0.                                                           (16)
```

This is the partial-sum triangle of binomial coefficients, not the
Bernoulli-number sequence. Its bulk rule and terminal boundary are

```text
B(r,s)=B(r-1,s)+B(r-1,s-1),             0<=s<r,
B(r,r)=2^r=2B(r-1,r-1),                 r>=1.                        (17)
```

For `n>=0`, put

```text
P_n=B(2n,n)=sum_(k=0)^n C(2n,k),
M_n=B(2n,n-1)=sum_(k=0)^(n-1) C(2n,k),                              (18)

M_0=0.
```

Reflection `k<->2n-k` has one fixed layer, `k=n`. Therefore

```text
P_n=(4^n+C(2n,n))/2,
M_n=(4^n-C(2n,n))/2,                                                 (19)

P_n+M_n=4^n,
P_n-M_n=C(2n,n).                                                     (20)
```

Here

```text
P_n=A032443(n)=1,3,11,42,163,638,...,

M_n=0,1,5,22,93,386,1586,...,
M_(n+1)=A000346(n).                                                  (21)
```

The indexing shift in (21) explains why the prompt's third row begins
with `1` rather than `0`.

Equivalently, among all subsets of a `2n`-set:

```text
P_n counts |S|<=n,
M_n counts |S|<n,
C(2n,n) counts the balanced equator.                                 (22)
```

The ordinary generating functions, as formal power series (and
analytically for `|z|<1/4`), are

```text
sum_(n>=0) P_n z^n
 =1/2[(1-4z)^(-1)+(1-4z)^(-1/2)],

sum_(n>=0) M_n z^n
 =1/2[(1-4z)^(-1)-(1-4z)^(-1/2)].                                  (23)
```

Thus the two half-space sequences are the symmetric and antisymmetric
central projections of the same Boolean growth law.

## 4. Catalan numbers are the boundary leakage

Let

```text
Cat_n=C(2n,n)/(n+1).                                                  (24)
```

Using

```text
C(2n,n)=(4-2/n)C(2n-2,n-1),                                        (25)
```

in (19) gives, for `n>=1`,

```text
P_n=4P_(n-1)-Cat_(n-1),
M_n=4M_(n-1)+Cat_(n-1).                                             (26)
```

The same Catalan packet crosses the reflection boundary with opposite
sign. This is the exact mechanism behind the empirical recurrences

```text
3=4*1-1,       11=4*3-1,       42=4*11-2,       163=4*42-5,

5=4*1+1,       22=4*5+2,       93=4*22+5,       386=4*93+14.         (27)
```

The "irregularity" is therefore localized: `Cat_(n-1)` is exactly the
signed difference between fourfold growth of the old half-space and
the new central boundary. Catalan numbers separately count the
standard one-sided Dyck/ballot paths; no identification with all paths
touching this particular Boolean boundary is needed here.

## 5. Tournament arc slots and the source of the two

A labelled tournament on `v` vertices has

```text
E=C(v,2)                                                            (28)
```

unordered arc slots. Fix a reference orientation. Reversing any subset
of the `E` slots gives a unique tournament, and the Hamming shell with
exactly `k` reversals has size

```text
C(E,k)=E^(falling k)/k!.                                             (29)
```

Consequently

```text
# labelled tournaments=2^E=sum_(k=0)^E C(E,k).                       (30)
```

This gives a faithful tournament realization of (13): the factor two
is the pair of orientations available independently in each arc slot,
while the triangular number `C(v,2)` counts the slots. It would be
wrong to write `2^v` or to say that tournaments cause the algebraic
base `2`; both are realizations of the same Boolean choice product.

## 6. The plus-minus quadratic form is switching energy

For a fixed reference order, let `a_ij in {+1,-1}` encode the
orientation/sign of every edge of `K_v`. For `x_i in {+1,-1}`, define

```text
Q_a(x)=sum_(i<j) a_ij x_i x_j.                                      (31)
```

Vertex switching by `x` sends

```text
a_ij -> a_ij x_i x_j.                                                (32)
```

Thus `Q_a(x)` is exactly the total signed edge bias of the switched
representative. Moreover

```text
Q_a(-x)=Q_a(x),                                                       (33)
```

so the global sign is a gauge. This is the precise tournament/two-graph
reading of the linked MathOverflow min--max expression. The central
binomial term in (19) is the balanced sign equator when the number of
vertex-sign slots is even.

This reformulation does **not** prove that the normalized min--max limit
exists or determine its value. Cycle products, not the scalar total
bias alone, are the switching invariants that must be retained.

## 7. Scope and convergence guard

Equations (12)--(14) terminate because `N` is a nonnegative integer.
They do not license the unqualified identity

```text
sum_(k>=0) x^(falling k)/k!=2^x                                    (34)
```

for arbitrary complex `x`. The generalized series
`sum_k binom(x,k)z^k` is ordinary inside `|z|<1`. At `z=1`, a
nonintegral `x` gives absolute convergence for `Re(x)>0` and
conditional convergence for `-1<Re(x)<=0`; outside that range one
needs a separately declared summation convention such as an Abel
value. None of these nonterminating extensions is used here.

Likewise, (31)--(33) are an exact representation of the open
plus-minus quadratic-form problem, not evidence for its asymptotic
limit. The theorem establishes the operator dictionary, central layer
split, and tournament Hamming-shell interpretation only.

## 8. Exact companion

Run:

```text
python3 04-computation/delta_exponential_central_newton_split_thm2412.py
python3 -O 04-computation/delta_exponential_central_newton_split_thm2412.py
```

The standard-library rational companion checks:

- `39` scaled falling-basis identities;
- `33` mixed-sign umbral intertwiners;
- `117` terminating exponential identities;
- the central split and Catalan recurrences through `n=12`;
- tournament Hamming shells through `v=9`; and
- the hostile failures `D_1(x^2)!=2x` and
  `#T_4=2^6!=2^4`.

Both modes reproduce:

```text
05-knowledge/results/delta_exponential_central_newton_split_thm2412.out
```
