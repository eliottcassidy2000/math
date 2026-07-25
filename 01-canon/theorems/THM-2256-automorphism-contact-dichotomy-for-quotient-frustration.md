---
id: THM-2256
title: "Automorphism-contact dichotomy for quotient frustration"
status: >
  PROVED + VERIFIED-EXACT. For every fixed finite quotient tournament R,
  the minimum THM-2249 forced energy outside the scaled automorphism
  transports has exactly one of two growth scales. It is uniformly bounded
  if some automorphism layer has zero two-way interaction with a
  nonautomorphism layer, and is Theta(N), with universal lower bound N-1,
  otherwise. In the bounded class the contact witness gives an explicit
  constant ceiling. Every transitive quotient has exact floor one for all
  block sizes; the directed triangle lies in the linear class and retains
  THM-2249's sharper exact floor. This concerns the forced quotient envelope,
  not the residual internal-block response.
source: codex-2026-07-25-quotient-frustration-scale
depends_on:
  - THM-2249-directed-triangle-forced-quotient-frustration
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity
script: 04-computation/tournament_automorphism_contact_dichotomy_thm2256.py
output: 05-knowledge/results/tournament_automorphism_contact_dichotomy_thm2256.out
script_sha256: 27ef0d7a52950b9022df14691de13c0dfbffb903826d0436a31be87c8c80fdcc
output_sha256: 68b0ea6961fc1d7434adaab21d89b852e8edd7383cf9c95a9cac9680b9a9ee91
hash_basis: working-tree bytes (LF)
---

# THM-2256 -- quotient frustration is either bounded or linear

Let `R` be a tournament on a finite label set `C`, with `|C|>=2`. For
permutations `sigma,tau in Sym(C)`, use THM-2249's ordered layer
interaction

```text
f_R(sigma,tau)
 =#{(i,k):i->_R k and tau(k)->_R sigma(i)}.          (1)
```

Put

```text
d_R(sigma)=f_R(sigma,sigma),

w_R(sigma,tau)=f_R(sigma,tau)+f_R(tau,sigma).        (2)
```

The diagonal term is the number of quotient arcs reversed by `sigma`.
Consequently

```text
d_R(sigma)=0 iff sigma in Aut(R).                    (3)
```

THM-2249 also proves the less obvious off-diagonal rigidity

```text
sigma,tau in Aut(R), sigma!=tau
 implies w_R(sigma,tau)>=1.                          (4)
```

For an integer `N>=1`, let `T_R(N)` be the nonnegative integer matrices
with every row and column sum equal to `N`. Define the nonfree forced-energy
floor

```text
phi_R(N)
 =min{F_R(X):
       X in T_R(N),
       X!=N P_sigma for every sigma in Aut(R)},       (5)
```

where `F_R` is THM-2249's intrinsic forced-pair energy

```text
F_R(X)=sum_(i->_R k) sum_(l->_R j) X_ij X_kl.        (6)
```

The set in (5) is nonempty because `|C|>=2` and there is a permutation
which is not an automorphism.

## 1. The contact predicate

Call `(sigma,tau)` a **zero contact** when

```text
sigma in Aut(R),
tau notin Aut(R),
w_R(sigma,tau)=0.                                   (7)
```

This is an intrinsic finite predicate on the quotient. The theorem is the
following exact scale dichotomy.

> **Automorphism-contact dichotomy.**
>
> - If a zero contact exists, then for every `N>=1`,
>
>   ```text
>   1<=phi_R(N)<=
>      min_{(sigma,tau) satisfying (7)} d_R(tau).     (8)
>   ```
>
>   In particular, `phi_R(N)=Theta_R(1)`.
>
> - If no zero contact exists, then there is a constant `C_R` such that
>   for every `N>=1`,
>
>   ```text
>   N-1<=phi_R(N)<=C_R N.                            (9)
>   ```
>
>   In particular, `phi_R(N)=Theta_R(N)`.

Thus no intermediate unbounded sublinear scale and no quadratic minimum can
occur outside the automorphism axes.

## 2. Hall layers turn the problem into a positive quadratic form

Take any `X in T_R(N)`. Repeated Hall matching gives a layer
decomposition

```text
X=sum_(s=1)^N P_(sigma_s).                           (10)
```

Let `n_sigma` be the multiplicity of `sigma`. Expanding (6), exactly as in
THM-2249, gives

```text
F_R(X)
 =sum_sigma d_R(sigma)n_sigma^2
  +sum_(sigma<tau) w_R(sigma,tau)n_sigma n_tau.      (11)
```

Every coefficient in (11) is a nonnegative integer. Although `X` can have
several Hall decompositions, every resulting right side equals the same
intrinsic number `F_R(X)`; no uniqueness is asserted or needed.

THM-2249's zero theorem says that `F_R(X)=0` exactly on the excluded
matrices `N P_sigma`, `sigma in Aut(R)`. Therefore every value minimized in
(5) is a positive integer:

```text
phi_R(N)>=1.                                         (12)
```

## 3. A zero contact gives a constant defect

Suppose `(sigma,tau)` satisfies (7), and take

```text
X_N=(N-1)P_sigma+P_tau.                              (13)
```

For `N=1`, this is the nonautomorphism layer `P_tau`. For `N>=2`, one row
contains positive entries in the two distinct columns `sigma(i)` and
`tau(i)`, so `X_N` is not a scaled permutation matrix. It is therefore
admissible in (5).

Equations (2), (7), and (11) give

```text
F_R(X_N)=d_R(tau),                                   (14)
```

independently of `N`. Minimize (14) over all zero contacts and combine with
(12) to obtain (8).

The mechanism is literal: one local nonautomorphism defect is invisible to
all `N-1` automorphism layers, so only its diagonal reversal count remains.

## 4. Without a zero contact every defect sees the bulk

Now suppose (7) has no solution and fix a Hall decomposition of a nonfree
transport. Let

```text
k=sum_(tau notin Aut(R)) n_tau                       (15)
```

be the number of nonautomorphism layers.

If `k>=1`, (3) gives

```text
sum_(tau notin Aut(R)) d_R(tau)n_tau^2
 >=sum_(tau notin Aut(R)) n_tau^2
 >=k.                                                (16)
```

The absence of zero contacts makes every automorphism--nonautomorphism
coefficient at least one. Their total contribution to (11) is therefore at
least

```text
k(N-k).                                              (17)
```

Equations (16)--(17) imply

```text
F_R(X)>=k(N-k+1)>=N,             1<=k<=N.            (18)
```

If `k=0`, all layers are automorphisms. Because `X` is not free, at least
two distinct automorphisms occur. Equation (4) then gives

```text
F_R(X)
 >=sum_(sigma<tau)n_sigma n_tau
 >=N-1.                                               (19)
```

The last minimum is attained by multiplicities `(N-1,1)`. Equations
(18)--(19) prove the lower bound in (9).

For the upper bound, if `Aut(R)` contains two distinct elements, use
`N-1` copies of one and one copy of the other; (11) is a fixed multiple of
`N-1`. If `Aut(R)` is trivial, fix any nonautomorphism `tau` and combine it
with `N-1` identity layers. Its energy is

```text
w_R(id,tau)(N-1)+d_R(tau).                           (20)
```

Increasing a fixed constant handles `N=1` in the first case. This supplies
some finite `C_R` in (9) and completes the dichotomy.

## 5. Sharp boundary models

Let `T_q` be the transitive tournament on `q>=2` ordered labels.
Its automorphism group is trivial. If `tau` swaps one adjacent pair, then
exactly that arc is reversed, while neither cross interaction with the
identity contains a forced pair:

```text
d_(T_q)(tau)=1,
w_(T_q)(id,tau)=0.                                   (21)
```

Thus (8) and positivity give the exact all-size law

```text
phi_(T_q)(N)=1                         for every N>=1. (22)
```

This explains why a transitive quotient can have a strong whole-block
prefix theorem while its bare forced-pair envelope has no macroscopic
floor: the internal-block response is the missing sidecar.

For the directed triangle, THM-2249's six-type table has no zero contact.
It lies in the linear class and has the sharper exact law

```text
phi_(C_3)(1)=3,
phi_(C_3)(N)=3(N-1)                  for 2<=N<=4,
phi_(C_3)(N)=2N+1                    for N>=5.        (23)
```

Thus (22) and (23) realize both branches at the smallest nontrivial
quotient scales.

## 6. What transfers and what does not

The theorem classifies `F_R`, a canonically forced subset of full reversal
cost. Since THM-2249 proves

```text
F_R(X)<=G_R(X),                                       (24)
```

the linear branch supplies a genuine linear lower tax for every nonfree
transport. The bounded branch says only that this particular quotient
envelope cannot supply such a tax. It does not say the full response is
bounded: internal orientations, pins, and context cuts can add a larger
residual, as THM-2242 does for transitive quotients.

The typed connection is

```text
source:
  Hall layers of a pinned equal-block quotient transport;

target:
  the asymptotic scale of its forced quotient-pair energy;

map:
  replace the transport by the multiplicity vector of permutation layers
  and retain the complete nonnegative interaction matrix (d_R,w_R);

preserved:
  every quotient-forced reversal, automorphism sectors, integrality, and
  the exact zero set;

destroyed:
  internal-block reversal cost, pins inside blocks, continuation context,
  and the identity of individual vertices within one layer;

needed sidecar in the bounded class:
  an internal response, pin, or context observable that detects a
  zero-contact defect.                                      (25)
```

## 7. Independent exact audit

The companion checks every `1,098` labelled tournament of orders two
through five. Its `124,468` diagonal checks recover automorphisms exactly,
and `774` distinct-automorphism pairs verify (4). The exact class census
`(bounded,linear)` is

```text
q=2: (2,0),
q=3: (6,2),
q=4: (48,16),
q=5: (680,344).                                      (26)
```

It also minimizes (11) over all `1,312,036` Hall-layer multisets for
orders two through four and `1<=N<=4`, checking the appropriate side of
(8)--(9), and verifies the adjacent-transposition witness (21) through
order five. Normal and optimized runs reproduce the stored transcript
exactly.

QED.
