---
id: THM-3429
title: "Prime-fibre activity descent for mixed-order half-twist seven-covers"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Transverse target-free
  odd joint-period cap-seven covers use only primes
  {3,5,13,17,29}; outside divisors 13,29,51 only the towers 17^a and
  5*17^a remain.  Q=51 forces an affine lift cocycle.  This is a reduction,
  not a classification; arbitrary common time and LRC(14) remain open.
source: codex2 base descent plus rank7-reflection fixed-fibre strengthening, 2026-08-15
audit: independent statement-universe, prime-power pullback, strict capacity, fixed-fibre parity, activity budgets, support/tower reconstruction, Q51 lift hostile, standalone p7 closure, dual companions, hashes, normal/-O, security, and documentation CLEAN
depends_on:
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3421-prime-half-twist-rank-seven-classification
related:
  - THM-3425-half-twist-rank-six-primitive-breaker-profile-closure
  - THM-3428-rough-maximal-order-half-twist-rank-seven-exclusion
script: 04-computation/lrc_prime_fibre_activity_descent_thm3429.py
output: 05-knowledge/results/lrc_prime_fibre_activity_descent_thm3429.out
script_sha256: 2de78875810a95eaa1ccdddfd4a91cac4db5ea32bcbe384a635f1045d5a692c7
output_sha256: 75b6d84866f38c826b4e26e1374f7d45aa233704de47e37dc3499bf5b0689f99
semantic_sha256: c46b24a98d32f83e1933276372e8bc1688ac9d3c3efa54c154ec19fcd1094c2d
independent_script: 04-computation/lrc_rank7_prime_fibre_activity_descent_thm3429.py
independent_output: 05-knowledge/results/lrc_rank7_prime_fibre_activity_descent_thm3429.out
independent_script_sha256: f89f31256715615babe7ffec792363f9149134a5d02cf20db906ed45c210f1a8
independent_output_sha256: ba5196cfa2f8ec82c2cd4fb0320181375334c85525fdb8a5b9e891f2f613c2cf
independent_semantic_sha256: c3b8596233d20f55e5c0bf0419c5f3f627cbf371e3004a5f40ab8a582a09b9f8
hash_basis: LF-normalized bytes
---

# THM-3429 -- prime-fibre activity descent for mixed-order half-twist seven-covers

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

For odd `Q`, a positive transverse representative `1<=r<=Q-1`, and a sheet
`ell` modulo `Q`, put

```text
B_(Q,r)={ell: ||r(2ell+1)/(2Q)||<1/14},
m_Q(r)=Q/gcd(Q,r).                                      (1)
```

Let

```text
A_6={8,9,10,11,12,15,23,25}.                           (2)
```

Assume that no member of `A_6` divides `Q`, and that **distinct** positive
transverse representatives `R=(r_1,...,r_s)`, `s<=7`, satisfy

```text
1<=r_i<=Q-1.                                             (3a)
```

This is the sign-normalized transverse universe `Q not|r_i`; it excludes the
universal residue `r=0` and the empty residue `r=Q`.  Repeated or
sign-equivalent masks can be deleted before this normalization.  Assume
further that

```text
union_i B_(Q,r_i)=Z/QZ,
L(R)=lcm_i m_Q(r_i)=Q.                                  (3)
```

For a prime `p|Q`, call owner `i` **p-active** when

```text
v_p(m_Q(r_i))=v_p(Q),
```

equivalently `p` does not divide `r_i`; otherwise call it **p-inactive**.
Then:

1. if at least one owner is `p`-inactive, then

   ```text
   p in {3,5,17,29};                                    (4)
   ```

2. if every owner is `p`-active, then

   ```text
   p in {13,29}.                                        (5)
   ```

More generally, `(5)` holds whenever there is no **even** `p`-inactive
residue.  In particular, each of the mixed primes `3,5,17` must divide an
even selected residue.

Consequently

```text
prime_support(Q) subset {3,5,13,17,29},                 (6)
Q/m_Q(r_i)=gcd(Q,r_i) is {3,5,17,29}-smooth for every i. (7)
```

Since `(2)` forces `v_3(Q),v_5(Q)<=1` and forbids `15|Q`, a composite `Q`
with no divisor among `13,29,51` is necessarily

```text
Q=17^a       (a>=2),
or Q=5*17^a  (a>=1).                                    (8)
```

Formula `(8)` is only a reduction of the remaining odd composite lane.  The
literal atoms at `13,29,51`, together with global divisor ancestry, already
give rank-seven covers on their target-free multiples; this does **not** say
that the primitive joint-period certificate in `(3)` descends from that
divisor.  The theorem does not exclude either tower in `(8)` and does not
assert completeness of the candidate rank-seven divisor antichain.

## 2. Connection contract

THM-3416 closes literal half-twist cap six on every modulus.  THM-3421 closes
literal full-order half-twist cap seven on primes.  The present theorem uses
one prime fibre to decide whether a lower-order arm descends or instead spends
all six remaining owners on that fibre.

| field | exact connection |
|---|---|
| source | transverse joint-period odd half-twist covers by at most seven literal blocks |
| target | active/inactive owner profiles over one prime coordinate |
| map | `Z/QZ -> Z/(Q/p)Z` together with one fibre of size `p` |
| preserved | strict endpoints, literal union, owner count, and the full `p`-adic order bit |
| destroyed | the integer lift of an active coefficient modulo `2p` |
| required sidecar | the affine lift character in Section 6 |
| positive boundary | the mixed `(3,17,17,51,51,51,51)` atom at `Q=51` |
| cheapest hostiles | `p=7` for strict endpoints and `Q=51,p=17` for recursive descent |

The corrected near miss is a support-only prime projection: away from the
reflection-fixed fibre, congruent active coefficients modulo `2p` need not
induce the same local block.  The least-used coordinate is their relative
integer lift.

## 3. Exact inactive descent

Fix `p|Q` and write

```text
Q=pN.                                                   (9)
```

The reduction map `pi:Z/QZ -> Z/NZ` has fibres

```text
F_y={y+Nt:t in Z/pZ}.                                  (10)
```

If `r=pu` is `p`-inactive, direct division of the strict phase word gives

```text
B_(Q,pu)=pi^(-1)(B_(N,u)).                              (11)
```

This remains true when `p^2|Q`; no coprimality between `p` and `N` is used.
Condition `L(R)=Q` says that at least one owner is `p`-active.  Hence, if an
inactive owner exists, there are at most six inactive owners.

If their descended masks in `(11)` covered all of `Z/NZ`, THM-3416 would
force a member of `A_6` to divide `N`, hence to divide `Q`, contradicting
`(2)`.  Choose a base sheet `y` missed by every inactive mask.  On `F_y`, the
active blocks alone must cover all `p` points.

## 4. The strict fibre capacity

For an active residue `r`, substitution of `(10)` into `(1)` gives the phase
set

```text
r(2(y+Nt)+1)/(2Q)
 = r(2y+1)/(2pN) + rt/p                 (mod 1).        (12)
```

Because `p` does not divide `r`, the second term runs through a translate of
the complete `p`-grid.  The danger region is one open circular arc of length
`1/7`.  If `n` grid points lie in such an arc, the first and last have
distance at least `(n-1)/p` and strictly less than `1/7`.  Therefore

```text
n<=ceil(p/7).                                           (13)
```

The strict inequality is load-bearing at `p=7`: one active block meets a
seven-fibre in at most one point, not two.

There are at most six active owners when an inactive owner exists, so `(13)`
gives

```text
p<=6 ceil(p/7).                                         (14)
```

Writing `p=7k+s` solves `(14)` exactly.  Its prime solutions are

```text
{2,3,5,11,17,23,29}.                                   (15)
```

Oddness removes `2`, while target-freeness `(2)` removes `11,23`.  This proves
`(4)`.  Notice in particular that `7` is not in `(15)`.

The broader activity invoice from the same fibre argument remains useful.
Let

```text
a_p=#{i:i is p-active},             c_p=ceil(p/7).       (15a)
```

Then every prime divisor satisfies

```text
a_p c_p>=p,
a_p>=ceil(p/c_p).                                      (15b)
```

Indeed, if `a_p c_p<p`, the active owners cannot fill a fibre missed by the
inactive owners, so the inactive owners descend to a cover of `Z/NZ`.
Joint period gives `a_p>=1`, hence at most six inactive owners, contradicting
THM-3416 and target-freeness.  This argument includes repeated prime factors.

A target-free cap-seven cover has exactly seven owners.  On the prime support
eventually allowed by `(4)--(5)`, the inactive-owner defect budgets are

```text
|D_3|<=4,  |D_5|<=2,  |D_13|=0,  |D_17|<=1,  |D_29|<=1, (15c)
```

where `D_p={i:p|r_i}`.  Thus

```text
#{i:gcd(Q,r_i)=1}
 >=max(0,7-sum_(p|Q)|D_p|).                            (15d)
```

Before the fixed-fibre gate is applied, the equality case at `p=7` would
force all seven owners active and exactly one hit by each owner on every
seven-point fibre.  Hence it would be a global exact partition.  Section 5
rules this case out already; the separate parity companion gives an
independent contradiction by forcing all seven residues even and then
colliding them at the common reflection-fixed sheet.

## 5. The all-active prime fibre and parity gate

Suppose first that every owner is `p`-active.  Use the reflection-fixed base sheet

```text
y_0=(N-1)/2.                                            (16)
```

Then `2y_0+1=N`, so `(12)` becomes

```text
r(2t+1)/(2p).                                           (17)
```

Thus restriction to `F_(y_0)` is a literal full-order half-twist block on the
prime modulus `p`.  The restrictions of the at most seven global blocks cover
that fibre.  Repeated local masks are harmless and may be discarded.  By
THM-3421,

```text
p in {11,13,23,29}.                                    (18)
```

Target-freeness removes `11,23`, proving `(5)`.

The same argument needs only the absence of an **even** inactive residue.  If
`r=pu` is inactive, its restriction to the fixed fibre is constant, and

```text
||r(2(y_0+Nt)+1)/(2Q)||=||u/2||.                       (19)
```

It covers the whole fibre exactly when `u`, equivalently `r`, is even, and
misses it exactly when `u` is odd.  Thus if every inactive residue is odd,
the active restrictions still cover the fixed fibre and `(18)` applies.
Intersecting `(18)` with the mixed list `(4)` leaves only `p=29`; hence a
mixed coordinate `p=3,5,17` necessarily has an even inactive owner.

If a prime divides the coindex `Q/m_Q(r_i)=gcd(Q,r_i)`, owner `i` is inactive
at that prime, so `(4)` proves `(7)`.  Combining `(4)` and `(5)` proves `(6)`,
and the elementary divisibility deductions after `(7)` prove `(8)`.

## 6. Why this is not yet a finite recursive compiler

For an active residue write, relative to one prime fibre,

```text
r=bar(r)+2pk,        0<=bar(r)<2p,        k mod N.      (20)
```

Equation `(12)` refines to

```text
bar(r)t/p + bar(r)(2y+1)/(2pN)
             + k(2y+1)/N                    (mod 1).   (21)
```

The last term is an affine character of the base sheet.  It vanishes modulo
one on the fixed fibre `(16)`, but not on a general missed fibre.  Therefore
the order multiset, collision polynomial, joint period, and coefficients
reduced modulo `2p` do not determine the recursive fibre blocks.  A faithful
factor-tree state must retain the relative lift vector `k mod N`, or a proved
quotient of its affine characters.

The `Q=51` atom is the sharp stopping hostile:

```text
R=(1,11,12,18,23,34,35),
orders=(51,51,17,17,51,3,51).                          (22)
```

Take `p=17,N=3`.  The unique inactive block is `r=34`; it descends to the
order-three mask and covers only the fixed base sheet `y=1`.  On each missed
base sheet `y=0,2`, the six active restrictions all have size three and cover
the seventeen-point fibre, with total overlap one.  In particular `r=1` and
`r=35` have the same reduction modulo `34`, but their lift indices `k=0,1`
make their missed-fibre blocks disjoint.  They coincide only on the fixed
fibre already covered by `r=34`.

Thus `(22)` is not a counterexample to the descent theorem; it is the reason
the theorem stops at finite prime support rather than claiming that a mixed
cover descends to a smaller cover.

## 7. Exact companion and non-consequences

Run

```bash
python3 04-computation/lrc_prime_fibre_activity_descent_thm3429.py
python3 -O 04-computation/lrc_prime_fibre_activity_descent_thm3429.py
```

The standard-library companion pins THM-3416 and THM-3421; solves `(14)` over
all primes below a rigorous cutoff; checks the exact projection and fibre
bound on a finite hostile bank including prime powers; reconstructs every
fibre in `(22)`; and verifies the strict `p=7` boundary.  Normal and optimized
outputs are intended to be byte-identical.

The independently written `lrc_rank7_prime_fibre_activity_descent_thm3429.py`
is retained as a broader finite audit of `(11)--(15)`: it checks every
transverse residue on every prime fibre of all odd composite `Q<=315`.  Its
`p=7` equality row is a sharp capacity boundary for seven active owners, not a
surviving mixed-prime lane; Section 5's fixed-fibre restriction is the
additional obstruction.

The independent parity companion
`odd_target_free_p7_half_twist_parity_obstruction.py` freezes the strict
`p=7` equality boundary by a different all-`Q` proof.

This theorem classifies neither tower in `(8)`.  It supplies no arbitrary-time
cover, physical runner row, rank-seven antichain completeness, or LRC(14)
decrement.  **QED.**
