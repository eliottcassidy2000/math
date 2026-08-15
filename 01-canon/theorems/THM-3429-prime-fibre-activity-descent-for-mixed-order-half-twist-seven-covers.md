---
id: THM-3429
title: "Prime-fibre activity descent for mixed-order half-twist seven-covers"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For an odd composite
  joint-period-Q half-twist cover by at most seven blocks, if Q/p has no
  six-block cover then the number a_p of owners active on the p-fibres obeys
  a_p*ceil(p/7)>=p.  Consequently every target-free primitive mixed-order
  seven-cover has a prime factor in {3,5,7,17,29}; in particular none exists
  when spf(Q)>29.  The exceptional-prime activity sets have exact Boolean
  defect budgets, and p=7 forces a fibrewise exact partition.  The five
  remaining small-prime lanes, fixed zero, arbitrary common time, and LRC(14)
  remain open.
source: independent prime-fibre proposal plus codex2 proof reconstruction and hostile audit, 2026-08-15
audit: direct fibre derivation including repeated prime factors; joint-period valuation audit; THM-3416 descent and THM-3428 full-order closure; Q=51 positive and Q=39 nonprimitive hostile; normal/optimized/stored exact replay
depends_on:
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3428-rough-maximal-order-half-twist-rank-seven-exclusion
related:
  - THM-3425-half-twist-rank-six-primitive-breaker-profile-closure
  - THM-3426-rough-composite-odd-interval-collision-and-dyadic-clique-law
script: 04-computation/lrc_rank7_prime_fibre_activity_descent_thm3429.py
output: 05-knowledge/results/lrc_rank7_prime_fibre_activity_descent_thm3429.out
script_sha256: 5ea4a9830f77a984aa7db83124ada94a498ed8f7687164533919545cd207f23f
output_sha256: 2ad7b43fd9afb14539438f2a17d4f458a7af6dbf3bd3deb947f8740099d7c316
semantic_sha256: d65463e762b3ba39e6cbbf13daeba73f9b689f9051ec460ac7d6ec27eaeece85
hash_basis: LF-normalized bytes
---

# THM-3429 -- prime-fibre activity descent for mixed-order half-twist seven-covers

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

For odd `Q` and `r` modulo `2Q`, put

```text
B_(Q,r)={ell in Z/QZ: ||r(2ell+1)/(2Q)||<1/14},          (1)
m_Q(r)=Q/gcd(Q,r).                                       (2)
```

Let `R` be a family of at most seven residues whose blocks cover `Z/QZ`, and
assume its joint quotient period is primitive:

```text
lcm_(r in R) m_Q(r)=Q.                                   (3)
```

Fix a prime `p|Q`, put `M=Q/p`, and assume that `Z/MZ` has no half-twist
cover by at most six blocks.  Define

```text
a_p=#{r in R:p does not divide r},
c_p=ceil(p/7).                                           (4)
```

Then

```text
a_p c_p>=p,
a_p>=A(p):=ceil(p/ceil(p/7)).                            (5)
```

This includes repeated prime factors: no squarefreeness hypothesis is used.
If a base sheet not covered by the `p`-divisible owners is covered at equality
`a_p c_p=p`, then the active blocks partition its `p`-point fibre and each
active block meets that fibre in exactly `c_p` points.

Call `Q` **rank-six target-free** when no member of

```text
T={8,9,10,11,12,15,23,25}                               (6)
```

divides `Q`.  By THM-3416, this is equivalent to the absence of a literal
half-twist cover by at most six blocks.  If `Q` is odd, composite, target-free,
and admits a cover satisfying `(3)`, then

```text
some p in {3,5,7,17,29} divides Q.                       (7)
```

In particular,

```text
Q odd composite, spf(Q)>29
  => no joint-period-Q half-twist cover by at most seven blocks.  (8)
```

Equation `(7)` reduces the arbitrary mixed-order composite frontier to five
small-prime lanes; it does not close those lanes.

## 2. Exact fibre dichotomy

Fix a base sheet `ell_0 modulo M`.  Its fibre under
`Z/QZ -> Z/MZ` is

```text
ell_j=ell_0+jM,                  0<=j<p.                  (9)
```

If `p|r`, write `r=pr_0`.  Directly from `(1)`,

```text
||r(2ell_j+1)/(2Q)||
 =||r_0(2ell_0+1)/(2M)||,                                (10)
```

because changing `j` adds the integer `jr_0`.  Thus `B_(Q,r)` is the full
pullback of `B_(M,r_0)`.

If `p` does not divide `r`, the `p` phases on `(9)` are

```text
x_0+jr/p modulo 1.                                      (11)
```

Multiplication by `r` permutes `Z/pZ`, so `(11)` is a translated regular
`p`-gon.  The danger arc has open circular length `1/7`.  If it contains `h`
grid points, the circular distance from the first to the last is
`(h-1)/p<1/7`; hence

```text
h<=ceil(p/7)=c_p.                                        (12)
```

This is a strict-interval count.  It does not cancel `r`, assume that `M` is
coprime to `p`, or replace nonunits by a quotient set.

## 3. Descent proof and the joint-period gate

Suppose `a_p c_p<p`.  If the `p`-divisible owners missed some base sheet,
then `(12)` would let all active owners cover fewer than `p` points of its
fibre, contradicting the cover.  Therefore the `p`-divisible owners alone
cover every base sheet via `(10)`.

Condition `(3)` forces `a_p>=1`: if every `r` were divisible by `p`, every
order `m_Q(r)` would lose one power of `p`, so their least common multiple
could not contain the full `p`-part of `Q`.  The descended cover consequently
uses at most

```text
|R|-a_p<=6                                               (13)
```

blocks, contrary to the hypothesis on `M`.  This proves `(5)`.  If equality
holds in the fibre-capacity invoice, coverage forces every inequality in
`(12)` to be an equality and forbids overlaps, proving the partition clause.

The joint-period hypothesis is load bearing.  At `Q=39`, scaling the order-13
atom gives

```text
R=(3,6,9,15,21,27,33).                                  (14)
```

These seven blocks cover and `39` is target-free, but all residues are
divisible by `3` and their joint period is only `13`.  Thus `(14)` is a
literal pullback, not a primitive mixed-order counterexample to `(5)`.

## 4. Five-prime reduction

If `Q` is target-free, every divisor `M=Q/p` is target-free as well.  THM-3416
therefore supplies the local rank-six hypothesis used in `(5)` for every
prime divisor.  The activity floors are

| `p` | `3` | `5` | `7` | `11` | `13` | `17` | `19` | `23` | `29` | `p>=31` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `c_p` | `1` | `1` | `1` | `2` | `2` | `3` | `3` | `4` | `5` | `ceil(p/7)` |
| `A(p)` | `3` | `5` | `7` | `6` | `7` | `6` | `7` | `6` | `6` | `7` |

The last entry follows directly at `p=31`; for `p>=37`,

```text
6 ceil(p/7)<=6(p+6)/7<p.                                 (15)
```

Target-freeness excludes prime factors `11` and `23`.  If none of the five
primes in `(7)` divided `Q`, the table would give `a_p=7` for every `p|Q`.
Thus there are exactly seven owners and every one is coprime to every prime
dividing `Q`; all seven blocks have full quotient order `Q`.  Since then
`spf(Q)>7`, THM-3428 says that a cover is possible only for
`Q in {11,13,23,29}`.  None is composite and target-free.  This contradiction
proves `(7)`, and `(8)` follows immediately.

For bookkeeping, an odd target-free modulus has no factors `11,23`, has
`v_3(Q)<=1`, `v_5(Q)<=1`, and cannot contain both `3` and `5`.  Equation `(7)`
adds that at least one of `3,5,7,17,29` occurs.  These are restrictions on the
modulus, not a classification of positive covers.

## 5. Boolean activity and defect budgets

A target-free cover in the theorem has exactly seven owners.  For each prime
factor define the intrinsic subsets

```text
A_p={i:p does not divide r_i},       D_p={1,...,7}\A_p.  (16)
```

The active set records which owners retain the full `p`-primary quotient
period; its complement is the descent defect.  On the five surviving lanes,
`(5)` gives

```text
|D_3|<=4,  |D_5|<=2,  |D_7|=0,  |D_17|<=1,  |D_29|<=1. (17)
```

Every other allowed prime has zero defect.  Consequently the number of
full-order owners obeys the Boolean union bound

```text
#{i:gcd(Q,r_i)=1}
 =|intersection_(p|Q) A_p|
 >=max(0,7-sum_(p|Q)|D_p|).                              (18)
```

Useful pairwise floors are

```text
3*17 or 3*29: at least 2 full-order owners;
5*17 or 5*29: at least 4 full-order owners;
17*29:        at least 5 full-order owners.              (19)
```

At `p=7`, `(17)` says every owner is active and `(12)` gives at most one hit
per seven-point fibre.  Coverage therefore makes the seven blocks an exact
partition on every fibre.  In particular the global multiplicity polynomial
is `Q X`: OR and XOR agree, and the overlap defect is zero.  This equality
case is a stronger prospective handle on the seven-lane than scalar density.

The sets `(16)` are a Boolean incidence carrier, not a tournament: there is
no intrinsic orientation between two owners.  Pairwise collision or overlap
data would be an additional sidecar.

## 6. Sharp positive and stopping boundary

The target-free atom at `Q=51=3*17`,

```text
R=(1,11,12,18,23,34,35),                                (20)
```

is a genuine positive boundary.  It covers, has order histogram

```text
3^1 17^2 51^4,                                          (21)
```

and has activities

```text
a_3=5,                  a_17=6.                          (22)
```

Thus the `17`-activity floor is attained exactly, while four owners are full
order.  On each of the two base fibres missed by the unique `17`-divisible
owner, every one of the six active owners hits three points; the resulting
total mass `18` pays one overlap above the `17` sheets.  This explains why the
prime-fibre theorem reduces but cannot exclude the `3*17` lane.

The exact companion directly checks `(10)--(12)` for every residue on every
prime fibre of all `93` odd composite moduli through `315`, including repeated
prime factors.  It audits `2,142,888` pullback cells and `11,819,480` active
cells, verifies the threshold table for every odd prime through `997`, and
reconstructs `(14)` and `(20)--(22)`.  Normal and optimized executions match
the LF-normalized stored output exactly.

No completeness result for the five lanes, fixed-zero cover, zero-mode
cochain classification at rank seven, arbitrary common time, physical runner
cover, current, or LRC(14) follows.  **QED.**
