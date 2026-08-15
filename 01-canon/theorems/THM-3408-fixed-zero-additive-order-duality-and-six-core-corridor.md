---
id: THM-3408
title: "Fixed-zero additive-order duality and the six-core corridor"
status: >
  PROVED structural theorem + COMPUTER-ASSISTED PROVED finite arithmetic
  lemma + FINITE-EXACT q<=20000 census + INDEPENDENTLY AUDITED.  At fixed
  source centre zero, an owner of quotient order m covers the exact fraction
  alpha(m/gcd(m,q/n)) of the additive-order-n sheet stratum.  This gives an
  exact fractional-cover dual, and every family descends without loss to the
  lcm of its quotient orders.  Outside the five positive dilation bases, a
  cover by at most six owners must contain an owner of quotient order in
  {22,33,44,46,50,102}.  Exact primal/dual rational certificates exclude all
  base-free q from 15 through 20000; this bounded census is not an all-q
  classification.  No mobile-centre, arbitrary-cochain, physical-time, or
  LRC(14) conclusion is claimed.
source: root-2608-crouzeix-puzzle-2026-08-15
audit: self-contained stratum/fibre, lcm, weak-dual, arithmetic-cutoff, and prime-loss proofs; exact q21/q22/q102 hostiles; 15985 exact rational primal/dual games; independent fibre/cutoff/DFS/prime-loss/LP-type audit clean; normal and optimized outputs byte-identical
depends_on:
  - THM-3398-general-finite-mode-sheet-cover-cochain
  - THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight
related:
  - THM-3402-atomized-sheet-covers-and-constructive-cochain-locus
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
script: 04-computation/lrc_fixed_zero_additive_order_duality_thm3408.py
output: 05-knowledge/results/lrc_fixed_zero_additive_order_duality_thm3408.out
script_sha256: d95d459078a8cf3be5a61c34a8b6b46aaa20e3269726f35b8a4cfc94b28d5be6
output_sha256: c2d55d60d646ea77891a44333cf7390534b7824cc085e0b441affc1bff570605
semantic_sha256: ab976baf418e2610aee1558eb249a160998426812fa96dddbc83b5a4a528541e
hash_basis: LF-normalized bytes
---

# THM-3408 -- fixed-zero additive-order duality and the six-core corridor

**PROVED structural theorem + COMPUTER-ASSISTED PROVED finite arithmetic
lemma + FINITE-EXACT `q<=20000` census + INDEPENDENTLY AUDITED.**

## 1. Statement

Fix `q>=2`.  A transverse fixed-zero owner is a nonzero speed residue
`u mod q`.  Put

```text
m=q/gcd(u,q)>1,              u=(q/m)a,  a in (Z/mZ)^*,              (1)
D_(q;m,a)={ell in Z/qZ:14 min(a ell mod m,-a ell mod m)<m}.          (2)
```

Thus `m` is the owner's quotient order and `(2)` is exactly

```text
D_(q,u)(0)={ell:||u ell/q||<1/14}.                                  (3)
```

For every divisor `n|q`, let

```text
S_n(q)={ell in Z/qZ:additive_order(ell)=n}.                          (4)
```

For `n>1` and `m>1` define

```text
e_q(m,n)=m/gcd(m,q/n),                                               (5)
K(e)=floor((e-1)/14),
alpha(1)=1,
alpha(e)=2 #{1<=r<=K(e):gcd(r,e)=1}/phi(e)       (e>1).              (6)
```

Then the following hold.

### A. Exact additive-order density

For every mode `a in (Z/mZ)^*`,

```text
|D_(q;m,a) intersect S_n(q)| / |S_n(q)| = alpha(e_q(m,n)).           (7)
```

The fraction is independent of `a`.  This forgets the location of the
covered points inside the stratum, not their exact number.

### B. Exact lcm descent

For modes `(m_i,a_i)`, put `M=lcm_i(m_i)`.  Since every `m_i|q`, also
`M|q`, and reduction `rho:Z/qZ -> Z/MZ` satisfies

```text
D_(q;m_i,a_i)=rho^(-1)(D_(M;m_i,a_i)).                               (8)
```

Consequently

```text
union_i D_(q;m_i,a_i)=Z/qZ
iff union_i D_(M;m_i,a_i)=Z/MZ.                                     (9)
```

Thus every fixed-zero counterexample to a proposed quotient-order
classification descends, with all modes and the owner count preserved, to
the primitive modulus `M`.

### C. Exact fractional-cover dual

Let `lambda=(lambda_n)_(n|q,n>1)` be nonnegative rational weights with
`sum_n lambda_n=1`, and put

```text
beta_lambda(m)=sum_(n|q,n>1) lambda_n alpha(e_q(m,n)).                (10)
```

If `k` fixed-zero owners cover every sheet, then

```text
sum_(i=1)^k beta_lambda(m_i)>=1.                                     (11)
```

In particular, `max_(m|q,m>1) beta_lambda(m)<1/k` excludes a cover by at
most `k` owners.  Equivalently, the exact fractional game is

```text
tau(q)=min_lambda max_m beta_lambda(m)
      =max_y min_n sum_m y_m alpha(e_q(m,n)),                         (12)
```

where both simplices range over nontrivial divisor orders of `q`.

### D. The global six-core corridor

Call an integer **base-free** when it is divisible by none of

```text
A={15,16,18,20,24}.                                                  (13)
```

The exact exceptional sets are

```text
B={21,22,26,28,33,35,39,42,44,46,50,52,56,57,63,65,70,
   74,76,77,78,84,91,99,102,104,114,117,130,143,156,186,190},        (14)
H={22,33,44,46,50,102},                                              (15)
E=B\H.                                                              (16)
```

For every base-free `m>=2`,

```text
alpha(m)>4/25 iff m in B,
m notin B implies alpha(m)<=7/44<1/6,
alpha(m)=1/6 for m in E,
alpha(m)>1/6 for m in H.                                            (17)
```

If a base-free modulus `q` has a fixed-zero cover by at most six transverse
owners, then at least one owner has quotient order

```text
m_i in H={22,33,44,46,50,102}.                                      (18)
```

Hence `q`, and after `(9)` its primitive lcm modulus, is divisible by at
least one member of `H`.  This is a necessary corridor, not a construction
and not an all-`q` nonexistence theorem.

### E. Bounded exact census

For every base-free `15<=q<=20000`, exact rational primal and dual
certificates give

```text
tau(q)<=4/25<1/6.                                                    (19)
```

Thus no such `q` has a six-owner fixed-zero cover.  Among the `15,985`
tested moduli,

```text
tau(q)=4/25 iff 102|q,                                               (20)
```

which occurs `78` times; all other exact game values are strictly smaller.
Statements `(19)--(20)` are **FINITE-EXACT only**.  The suggested global
continuation `tau(q)<=4/25`, with equality exactly on base-free multiples of
`102`, remains open.

## 2. Inheritance and connection contract

[THM-3401](THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight.md)
proves the exact fixed-zero ranks through `q=28` and leaves open whether a
rank-at-most-six cover exists exactly on multiples of the five moduli in
`(13)`.  [THM-3398](THM-3398-general-finite-mode-sheet-cover-cochain.md)
supplies the typed quotient blocks.  The present theorem replaces literal
sheet subsets by their exact additive-order incidence matrix and pushes the
open reverse implication into the six-order corridor `(18)`.

The canonical hostile is `q=22`: its unit-stratum density is already
`1/5`, but two independent kernel strata force rank seven.  The corrected
near miss is the assertion that a density bound reconstructs a cover; it
does not retain within-stratum alignment.  The least-used sidecar is the
additive order of the sheet after projection to each owner quotient.

| field | exact connection |
|---|---|
| source | fixed-zero THM-3398 quotient modes `(m,a)` |
| target | the divisor-by-divisor matrix `alpha(e_q(m,n))` |
| map | send an owner to quotient order `m` and a sheet to additive order `n` |
| preserved | exact incidence density on every `S_n`, cover obstruction, lcm pullback |
| destroyed | mode orientation `a`, pointwise alignment, intersections, and mobile centre |
| required sidecar | literal mode blocks when one wants a construction rather than an obstruction |
| cheapest decisive tests | `q=21` prime-sheet loss, `q=22` kernel obligations, `q=102` sharp dual |

There is no intrinsic binary orientation here, so tournament analysis would
be cosmetic.  The relevant carrier is the rectangular owner-order versus
sheet-order payoff matrix.

## 3. Proof of the density formula

Every element of `S_n(q)` has a unique form

```text
ell=(q/n)x,                    x in (Z/nZ)^*.                         (21)
```

Let

```text
h=gcd(m,q/n),        m=he,        q/n=hb,        gcd(b,e)=1.          (22)
```

The image of `ell` modulo `m` is `h(bx mod e)` and has additive order
`e=e_q(m,n)`.  In particular `e|n`.  Reduction

```text
(Z/nZ)^* -> (Z/eZ)^*                                             (23)
```

is surjective.  Indeed a chosen unit modulo `e` can be lifted while avoiding
the single forbidden residue modulo each prime dividing `n` but not `e`;
the Chinese remainder theorem makes those choices simultaneous.  Therefore
every unit modulo `e` has exactly `phi(n)/phi(e)` preimages.  Multiplication
by the units `a` and `b` only permutes these fibres.

For a unit residue `r mod e`, cyclic distance scales by `h`:

```text
min(hr,he-hr)=h min(r,e-r).                                         (24)
```

When `e=1`, the image is zero and every sheet in the stratum is dangerous.
When `e>1`, the dangerous units are exactly

```text
{+-r:1<=r<=floor((e-1)/14), gcd(r,e)=1}.                             (25)
```

The two signs in `(25)` are disjoint.  Equal fibres in `(23)` now give
`(7)`.

For later use, the unit sheet stratum `n=q` gives the especially simple
identity

```text
|D_(q;m,a) intersect S_q(q)|/phi(q)=alpha(m).                        (26)
```

## 4. Proof of lcm descent and the fractional dual

The predicate in `(2)` depends only on `ell mod m_i`.  Because `m_i|M|q`,
the reductions to `m_i` factor through the surjection `rho:Z/qZ->Z/MZ`.
This proves `(8)`.  Taking unions and using surjectivity proves `(9)`.  At
modulus `M`, the mode is represented by speed `(M/m_i)a_i`, whose quotient
order is exactly `m_i`; no realizability information was lost.

Give every point of `S_n(q)` mass `lambda_n/phi(n)`.  This is a probability
measure on the nonzero sheets.  Formula `(7)` says that an owner of order `m`
has mass exactly `beta_lambda(m)`.  A literal cover of every sheet covers the
support of this measure, so the union bound gives `(11)`.  Linear-programming
duality on the finite rational payoff matrix gives `(12)`; only the weak
direction `(11)` is needed for every nonexistence result below.

## 5. Complete arithmetic cutoff

Put

```text
C_m(K)=#{1<=r<=K:gcd(r,m)=1},       omega=omega(m).                  (27)
```

Inclusion-exclusion gives

```text
C_m(K)=sum_(d|rad(m)) mu(d) floor(K/d)
       < K phi(m)/m + 2^omega.                                      (28)
```

Since `K=K(m)<m/14`,

```text
C_m(K)<phi(m)/14+2^omega.                                           (29)
```

If `alpha(m)>4/25`, then `C_m(K)>2phi(m)/25`; comparison with `(29)`
forces

```text
phi(m)/2^omega < 350/3.                                             (30)
```

Likewise `alpha(m)>7/44` forces

```text
phi(m)/2^omega < 616/5.                                             (31)
```

These are finite searches without an assumed bound on `m`.  If
`m=prod p^a`, then

```text
phi(m)/2^omega=prod_(p^a||m) p^(a-1)(p-1)/2.                        (32)
```

An increasing-prime recursion stops a branch as soon as the product in
`(32)` reaches the target; before that, the next prime is bounded directly
by `(p-1)/2`.  The exact companion exhausts every factorization once.  It
finds

```text
threshold       candidates including 1       largest candidate
350/3                       1710                    30030
616/5                       1829                    39270.             (33)
```

Exact inclusion-exclusion on these complete candidate sets gives `(14)` and
`(17)`.  Outside `B`, equality `alpha=7/44` occurs only at
`m=115,184,276`; this equality list is not needed below.

The six strict exceptions have values

```text
alpha(22)=alpha(33)=alpha(44)=alpha(50)=1/5,
alpha(46)=2/11,                 alpha(102)=3/16.                     (34)
```

Every member of `E` has value exactly `1/6`.

## 6. Proof of the six-core corridor

Suppose that a base-free `q` is covered by `k<=6` owners and that no owner
order lies in `H`.  Every divisor of `q` is base-free.  Apply `(11)` to the
unit stratum, or equivalently sum `(26)`:

```text
1<=sum_(i=1)^k alpha(m_i).                                          (35)
```

By `(17)`, each term is at most `1/6`, with equality outside `H` only on
`E`.  Hence `(35)` forces `k=6` and every `m_i in E`.

The complete prime-loss table is

```text
p= 7: 21,28,35,42,56,63,70,77,84,91
p=11: 77,99,143
p=13: 26,39,52,65,78,91,104,117,130,143,156
p=19: 57,76,114,190
p=31: 186
p=37: 74.                                                           (36)
```

It covers `E`, and in every displayed incidence

```text
2<=m/p<=14.                                                         (37)
```

Choose a displayed prime `p` from any one owner.  Then `p|q` and
`n=q/p>1`.  For every `E`-owner,

```text
e_q(m,n)=m                 if p does not divide m,
e_q(m,n)=m/p in {2,...,14} if p divides m.                           (38)
```

Thus nondivisible owners retain density `1/6`, every divisible owner drops
to `alpha(m/p)=0`, and the chosen owner drops strictly.  The density sum on
`S_(q/p)(q)` is at most `5/6`, contradicting `(11)`.  This proves `(18)`.

The forward half of the five-base picture is already exact: THM-3401 gives
covers of ranks `6,5,5,6,6` at `15,16,18,20,24`, respectively, and `(8)`
pulls each cover to every multiple.  The reverse half remains open beyond
the bounded census and the corridor above.

## 7. Three decisive hostiles

### 7.1 `q=21`: equality on units still loses prime sheets

At `q=21`, six order-`21` modes saturate the unit-stratum bound because
`alpha(21)=1/6`.  Their union misses exactly

```text
{3,6,7,9,12,14,15,18}=S_7(21) union S_3(21).                         (39)
```

The order-`3` stratum forces an owner of quotient order `7`, and the
order-`7` stratum forces one of quotient order `3`.  Hence the exact rank is
eight, with speed witness

```text
(1,2,3,4,5,7,8,10).                                                 (40)
```

This is the cheapest hostile to treating the unit stratum as sufficient.

### 7.2 `q=22`: strict density does not pay kernel obligations

At `q=22`, an order-`22` owner covers exactly `alpha(22)=1/5` of the ten
primitive sheets, while quotient orders `2` and `11` cover none of them.
Thus five order-`22` modes are necessary.  The stratum `S_2(22)` separately
forces an order-`11` owner, and `S_11(22)` forces an order-`2` owner.  The
exact rank is therefore seven, attained by

```text
(1,2,3,5,7,9,11).                                                    (41)
```

So membership in `H` is only a necessary corridor condition.

### 7.3 `q=102`: the sharp fractional hostile

Give the four sheet strata `(6,34,51,102)` weights

```text
(1,4,4,16)/25.                                                       (42)
```

Every owner order dividing `102` then has mass at most `4/25`; equality
occurs at owner orders `(2,3,17,102)`.  Conversely give those four owner
orders weights

```text
(7/150,7/150,4/25,56/75).                                           (43)
```

Every sheet stratum has expected coverage at least `4/25`; equality occurs
at `(6,34,51,102)`.  Hence

```text
tau(102)=4/25,
6 tau(102)=24/25<1,                                                  (44)
```

so six owners cannot cover.  The gap from the six-owner threshold is exactly
`1/6-4/25=1/150`.  This is sharp for the fractional method while remaining a
noncover.

## 8. Finite exact scan and its stopping boundary

For every base-free `15<=q<=20000`, the companion constructs an exact
rational saddle point for the divisor payoff matrix in `(12)`.  It checks
nonnegativity and normalization of both distributions, every owner upper
inequality, every sheet lower inequality, and equality of the two rational
values.  The universe contains `15,985` moduli; the largest matrix has
dimension `47`.  There are `59,210` nonzero sheet weights and `44,447`
nonzero owner weights across the stored certificate digest.

This proves `(19)--(20)` only in the displayed finite universe.  No
asymptotic compactness argument is supplied.  In particular, the six-order
corridor

```text
{22,33,44,46,50,102}                                                (45)
```

is still live globally, and a future base-free six-cover would have to use
one of those exact quotient orders before lcm descent.

There are two further scope losses.

1. Density forgets the alignment and intersections of modes.  It can refute
   a cover but cannot reconstruct one.
2. The source centre is literally zero.  THM-3405's `q=16`, centre `1/32`
   half-twist shows why an arbitrary common centre cannot be silently moved
   into this gauge.  Nonzero cochains and physical continuous time are still
   farther away.

Therefore the theorem neither proves the open reverse implication for all
`q` nor yields an LRC(14) decrement.

## 9. Exact companion

The standard-library companion uses strict integer danger inequalities,
inclusion-exclusion, an unbounded-by-design prime-factor recursion, and an
exact `Fraction` two-phase simplex.  It records

```text
32,493 literal stratum-density checks,
 2,954 equal-fibre checks,
 7,140 literal mode-pullback checks,
12,042 lcm-family checks,
15,985 exact primal/dual fractional games.                           (46)
```

It also exhausts the `q=21` and `q=22` literal set covers and freezes both
`q=102` saddle distributions.  There are no floating literals, random
choices, external solvers, network calls, or assertion-dependent gates.

Reproduce with

```text
python3 04-computation/lrc_fixed_zero_additive_order_duality_thm3408.py
python3 -O 04-computation/lrc_fixed_zero_additive_order_duality_thm3408.py
```

Artifact and semantic hashes are pinned in the frontmatter.

**QED.**
