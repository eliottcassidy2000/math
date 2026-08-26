---
id: THM-4206
title: "Rule 30 characteristic-address contrast-deck entropy decomposition"
status: >
  PROVED universal binary left-permutive triangular contrast-deck theorem +
  PROVED exact entropy/covariance consequences + FINITE-EXACT Rule-30 hostile
  and slope-one correlation bank + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  One representative per distinct left characteristic address is a fresh
  conditional Haar bit. Repeated-address outputs are recovered by a directed
  deck of Boolean contrasts depending only on larger-address representatives.
  This gives H(Y|N)=d and H(Y)=d+H(D|R). Pair covariance is block diagonal by
  address but is not sufficient: a three-cell Rule-30 hostile has two uniform
  cross-address pairs and a nonzero third Walsh character. Haar input is
  load-bearing; every named-seed Rule-30 prize remains OPEN.
source: open-frontiers-incoming-20260826b / Rule-30 incoming-signal lane
audit: >
  PASS. The universal result is proved by the extreme-pivot normal form and a
  unit-diagonal triangular inverse. The primary exact path uses packed truth
  tables, exhausts the minimal three-cell hostile, verifies its conditional
  deck and entropy ledger, and evaluates all 22,369,620 slope-one tail words
  through depth 12. A no-import scalar audit independently enumerates the
  hostile and all 87,380 tail words through depth 8. Normal and optimized
  streams agree. No finite sign pattern is extrapolated.
depends_on:
  - THM-4085-rule30-left-characteristic-independence-and-haar-extremes
related:
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-4047-rule30-left-front-affine-monodromy-clock
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
  - THM-4065-rule30-temporal-cylinder-transfer-infinite-zero-columns-and-next-zero
script: 04-computation/rule30_characteristic_address_contrast_deck_thm4206.py
output: 05-knowledge/results/rule30_characteristic_address_contrast_deck_thm4206.out
script_sha256: 084b7fdc65a1dcc712e9575aa6aa48f554b390c1f8fa06a74be76930debad7df
output_sha256: 43480dd835d16e88024f4d8b8cee588e8b0ca63d2d486b4de4cb608a7c4d965d
independent_audit_script: 04-computation/rule30_characteristic_address_contrast_deck_thm4206_independent_audit.py
independent_audit_output: 05-knowledge/results/rule30_characteristic_address_contrast_deck_thm4206_independent_audit.out
independent_audit_script_sha256: ee75595ec29d8e5e5ff8d1bcd6bf6ecd3c9f33b593758bcb165387818b748819
independent_audit_output_sha256: 664dc1044bdbec1317e99678bb709ef4b20e75744a422754f81945dccb2e962e
hash_basis: raw LF bytes
---

# THM-4206 -- the characteristic-address contrast deck

**PROVED in the universal binary left-permutive scope below; FINITE-EXACT for
the displayed Rule-30 banks.** THM-4085 proves that outputs at distinct left
characteristic addresses are iid under Bernoulli Haar input. Repeated
addresses do not merely subtract rank. They carry an ordered deck of Boolean
contrasts, and that deck gives an exact entropy decomposition.

The address order is intrinsic and directed. Pair covariance is symmetric and
loses higher Walsh characters, so it neither supplies a tournament nor
reconstructs the contrast deck.

## 1. Universal extreme-pivot normal form

Let

```text
F:{0,1}^Z -> {0,1}^Z
```

be a radius-one binary cellular automaton with local rule `f(l,c,r)` that is
left permutive. For a spacetime cell `v=(t,j)`, put

```text
lambda(v)=j-t,                    rho(v)=j+t,
Y_v=(F^t(X))_j.                                           (1)
```

The dependency interval of `Y_v` is `[lambda(v),rho(v)]`. Extreme-input
permutivity and the two-element alphabet give the exact Boolean form

```text
Y_v=x_(lambda(v)) xor
    g_v(x_(lambda(v)+1),...,x_(rho(v))).                    (2)
```

Indeed, after fixing every input except the extreme left coordinate, the
iterate is a permutation of one bit, hence either that bit or its complement.
Equation `(2)` records the complement choice as `g_v`.

## 2. The triangular characteristic contrast deck

Let `V=(v_i)_(i=1)^m` be a finite family of cells. List its distinct addresses
in decreasing order:

```text
alpha_1>alpha_2>...>alpha_d.                              (3)
```

Choose one representative output `R_k` from the class with address `alpha_k`.
For every other output `Y_i` in that class, define the contrast

```text
D_i=Y_i xor R_k.                                          (4)
```

Let `P=(x_(alpha_1),...,x_(alpha_d))` be the pivot vector. Let `N` contain
every other initial coordinate in the union of the dependency intervals of
`V`.

### Theorem 2.1 (triangular deck)

For every fixed value of `N`:

1. the map `P |-> R=(R_1,...,R_d)` is a bijection of `F_2^d`;
2. each contrast in class `k` is a Boolean function only of
   `(N,R_1,...,R_(k-1))`.

Consequently, under Bernoulli Haar input, `R` is uniform on `F_2^d` and is
independent of `N`. Conditional on `N` and the preceding representatives, all
outputs in class `k` are one fresh fair bit `R_k` plus fixed contrast bits.

### Proof

By `(2)`, the chosen representative in class `k` has the form

```text
R_k=x_(alpha_k) xor
    h_k(N,x_(alpha_1),...,x_(alpha_(k-1))).                (5)
```

Only larger-address pivots can occur to the right of `alpha_k`. Thus `(5)` is
triangular with diagonal coefficient one. Solve successively for
`x_(alpha_1),...,x_(alpha_d)` to obtain the bijection.

In `(4)`, the two copies of `x_(alpha_k)` cancel. The remaining initial
coordinates are strictly to the right of `alpha_k`; among the pivots these are
only `x_(alpha_1),...,x_(alpha_(k-1))`. Substitution through the triangular
inverse expresses the contrast using only `(N,R_1,...,R_(k-1))`.

For each fixed `N`, every `R` has exactly one pivot preimage. Hence its
conditional probability is `2^(-d)`, independent of `N`. This proves all
claims. `QED`

## 3. Exact entropy and covariance consequences

All entropies below are in bits. The coordinate change

```text
Y_V <-> (R,D)                                             (6)
```

is lossless: representatives are retained, and every other output is its
representative XOR its contrast. Since `D` is deterministic from `(R,N)`,
Theorem 2.1 gives

```text
H(Y_V | N)=d,
H(Y_V)=d+H(D | R).                                       (7)
```

For fixed `N`, the conditional output law is uniform on exactly `2^d` words.
Therefore every unconditional atom has mass at most `2^(-d)`, so

```text
H_infinity(Y_V)>=d.                                      (8)
```

For output signs `S_i=(-1)^(Y_i)`, every `S_i` is centered. If two cells have
distinct addresses, THM-4085, Theorem 1.1, makes them independent, hence

```text
Cov(S_i,S_j)=0.                                          (9)
```

If they have the same address, their shared pivot cancels and

```text
Cov(S_i,S_j)=E[S_i S_j]
              =E[(-1)^(Y_i xor Y_j)].                  (10)
```

Thus the sign-covariance matrix is block diagonal on the address partition.
This is only the degree-two shadow of the deck. Lower-address contrasts may
depend on higher-address representatives, so higher Walsh characters can
couple different covariance blocks.

## 4. The minimal Rule-30 covariance hostile

For Rule 30, use three Haar initial bits and the cells

```text
A=a_0(1)=x_1,                         address 1,
B=a_0(0)=x_0,                         address 0,
C=a_1(1)=x_0+x_1+x_2+x_1*x_2,        address 0.          (11)
```

In lexicographic `(A,B,C)` order, exhaustive counts over the eight input
triples are

```text
(000,001,010,011,100,101,110,111)
    = (1,1,1,1,0,2,2,0).                                (12)
```

The pairs `(A,B)` and `(A,C)` are both uniform with counts `(2,2,2,2)`. The
same-address pair `(B,C)` has THM-4085's counts `(1,3,3,1)`. Hence

```text
Cov((-1)^A,(-1)^B)=Cov((-1)^A,(-1)^C)=0,
Cov((-1)^B,(-1)^C)=-1/2,                                (13)
```

but the cross-block third character is

```text
E[(-1)^(A+B+C)]=1/2.                                    (14)
```

The triangular deck explains the missing dependence exactly. Choose

```text
R=(A,B),       N=x_2,       D=C xor B=A or x_2.          (15)
```

For each fixed `x_2`, the output support has four equiprobable words. Moreover,

```text
H(A,B,C | x_2)=2,
H(A,B,C)=5/2=2+H(D | A,B),
H(D | A,B)=1/2.                                         (16)
```

The covariance graph makes `A` isolated from `{B,C}`, yet `(14)` couples it
to that block. Therefore covariance, an oriented covariance tournament, or
the bare address partition is not sufficient. The lossless object is `(R,D)`
together with the law of the nonpivot sidecar `N`. Because covariance is
symmetric and has zeros, there is no intrinsic tournament orientation here;
the native directed object is the acyclic larger-address-to-smaller-address
dependency of Theorem 2.1.

## 5. Exact rational-ray entropy decomposition

Use THM-4085's rational observer

```text
j_t=b+floor(pt/q),       Y_t=(F^t(X))_(j_t),
0<=p<=q.                                                   (17)
```

Among `t=0,...,M-1`, the number of distinct addresses is

```text
d_M=M-floor(p(M-1)/q).                                   (18)
```

Choose one time over each address and form its contrast deck. Equations
`(7)--(8)` sharpen the earlier entropy floor to the exact identity

```text
H(Y_0,...,Y_(M-1))
  =d_M+H(D_M | R_M),
H_infinity(Y_0,...,Y_(M-1))>=d_M.                       (19)
```

When `p<=0` in the physical light cone, all addresses are distinct and the
deck is empty, recovering the iid trace. At the right-characteristic boundary
`p=q`, every address equals `b`. Taking `R=Y_0=x_b` and

```text
D_t=Y_t xor x_b                                          (20)
```

gives

```text
H(Y_0,...,Y_(M-1))=1+H(D_1,...,D_(M-1)).                (21)
```

Here the contrast vector depends only on the strict right tail and is
independent of `x_b`. Thus the entire repeated-address entropy question has
been isolated in its native contrast process rather than compressed to a
count or covariance matrix.

## 6. Finite-exact slope-one correlation bank

By translation take `b=0`. Put

```text
c_t=E[(-1)^(a_t(t)+x_0)]=E[(-1)^(D_t)].                 (22)
```

The primary exact truth-table path exhausts all `2^(2t)` strict-right-tail
words for each `1<=t<=12` and obtains

```text
t     1      2      3      4       5        6
c_t  -1/2   1/4   -1/4    5/32   -5/64    77/1024

t     7          8        9             10
c_t  -141/2048   39/512  -3273/65536    2785/131072

t     11               12
c_t  -21759/1048576    27905/2097152.                  (23)
```

A no-import scalar evolution independently reproduces `(23)` through `t=8`.
The signs alternate through this finite bank. No sign theorem, recurrence,
nonvanishing statement, decay bound, or asymptotic extrapolation is claimed.

## 7. Preservation, loss, and prize firewall

| source -> target | preserves | destroys / does not control | required sidecar |
|---|---|---|---|
| cells -> address partition | conditional innovation count `d` | repeated-address response | Boolean contrast deck |
| outputs -> sign covariance | every two-point character | `(14)` and higher Walsh law | full ordered deck and `N` law |
| deck -> address DAG | allowed dependency direction | Boolean function values | labelled contrast functions |
| Haar row -> named single seed | local Rule-30 equation only | fresh independent pivots and entropy | a deterministic temporal-transfer theorem |

Bernoulli Haar input is load-bearing in `(7)--(23)`. The distinguished
single-seed orbit is one deterministic row, not a typical Haar sample. None of
the iid representatives, Shannon identities, covariance statements, or finite
slope-one biases transfers to that orbit without a new theorem. Therefore
center non-eventual-periodicity, asymptotic balance, fixed-seed query lower
bounds, and all named Rule-30 prizes remain **OPEN**.

## 8. Exact verification

Run from the repository root:

```bash
python3 04-computation/rule30_characteristic_address_contrast_deck_thm4206.py
python3 -O 04-computation/rule30_characteristic_address_contrast_deck_thm4206.py
python3 04-computation/rule30_characteristic_address_contrast_deck_thm4206_independent_audit.py
python3 -O 04-computation/rule30_characteristic_address_contrast_deck_thm4206_independent_audit.py
```

The primary path represents complete truth tables as packed integers. The
independent path loops over scalar assignments and evolves ordinary rows. Both
use exact integer arithmetic, contain no random sampling, and agree on every
overlapping hostile, deck, and correlation value. The universal claims are
the proofs in Sections 1--3; the computations audit the Rule-30 applications
and do not establish an infinite slope-one law.

**QED.**
