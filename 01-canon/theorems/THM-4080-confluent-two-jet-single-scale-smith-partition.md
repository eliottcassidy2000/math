---
id: THM-4080
title: "Confluent two-jet single-scale Smith partition"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. Let p be prime and
  let 1<=s<p integer nodes lie in one p-adic ball with every pairwise
  difference of exact valuation e>=1. The first two Hasse jets on degree
  below 2s have p-Smith exponent multiset
  e*{0,0,1,...,s-1,s+1,...,2s-1}. Consequently, for two-jet observation on
  n consecutive nodes and n<=p(p-1), the complete p-primary partition is the
  multiset union of these profiles over the residue classes modulo p. This
  extends the earlier p<n<=2p pair band to a quadratic range. The bound s<p
  is load-bearing: at s=p the derivative-row matroid loses one rank and the
  displayed profile already fails for p=2 and p=3. Larger clusters and the
  complete global Smith form remain open. The independently audited
  2026-09-05 extension linked below closes the complete-residue s=p boundary
  with a one-unit saturation correction and extends consecutive n to p^2.
source: codex-frontier-synthesis-creative-20260825e / integral-observer wildcard
audit: >
  PASS. The primary DVR path checks 132 single-scale profiles across
  p=3,5,7,11, e=1,2,3, two unit-coordinate systems and translations; it also
  checks 34 direct consecutive-node matrices and 178 quadratic-range formula
  rows in 1,068,749 exact gates. The independent path reconstructs six Smith
  forms from every determinantal divisor, audits 136 residual matroids and
  1,896 unit witness minors through p=19, and checks 948 quadratic-range
  rows in 13,437,566 exact gates. Both normal and optimized outputs
  byte-match. The s=p hostiles and their derivative-rank failure are checked
  by both paths.
depends_on:
  - THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall
related:
  - THM-4000-centered-base-split-cubic-observer-and-tripotent-crt-atlas
  - THM-4064-rule30-cyclotomic-kernel-character-and-c60-alias-obstruction
script: 04-computation/confluent_twojet_single_scale_smith_thm4080.py
output: 05-knowledge/results/confluent_twojet_single_scale_smith_thm4080.out
script_sha256: 10c4a24dd28912cac1b167c89abb7b5cc26d389c7ee5400beb34695068cba733
output_sha256: 70f3bd1d5d8532b17e75a7d31c5b1be99fe81a4fc51022a9540e8960fa2b05c2
independent_audit_script: 04-computation/confluent_twojet_single_scale_smith_thm4080_independent_audit.py
independent_audit_output: 05-knowledge/results/confluent_twojet_single_scale_smith_thm4080_independent_audit.out
independent_audit_script_sha256: 02f270aadfd1c96fed7f3f9631a07d55dbb7dc108b5be177663442dea142229b
independent_audit_output_sha256: 7ade11c7d4338fef797e01ddd2bf83155ba8afafda6f88c60ddde9d46ef25eed
hash_basis: raw LF bytes
---

# THM-4080 -- a complete two-jet Smith layer inside one p-adic scale

**PROVED EXTENSION — 2026-09-05.** At `s=p`, the complete exponent list is
`0,0,e,...,(p-2)e,(p-1)e+1,(p+1)e,...,(2p-2)e,(2p-1)e-1`, with empty ranges
omitted. Thus the consecutive-node partition is now known through `n<=p^2`.
The old formula still fails at `s=p`; saturation of its missing derivative
line costs exactly one `p`, not one scale `e`.
[Full proof, exact controls and independent audit](../../05-knowledge/results/synthesis_20260905_wildcard_smith_boundary.md).
Larger multiscale clusters and higher jets remain open.

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** THM-4010
determines the zero layer, the first positive layer, and a two-node cluster.
The later layers become transparent while a residue cluster occupies only one
new p-adic scale: each minor pays for its column degrees and receives one
unit of scale back for every derivative row. The residual confluent
evaluation matroid says exactly when that lower tariff is attained.

## 1. Statement

Fix a prime `p`, integers

```text
e>=1,             1<=s<p,
```

and nodes

```text
x_i=a+p^e u_i,       0<=i<s,                         (1)
```

where the `u_i` are pairwise distinct modulo `p`. Equivalently,

```text
v_p(x_i-x_j)=e             for every i!=j.            (2)
```

On `Z[X]_<2s`, form the square two-jet Hasse observer

```text
J(P)=(P(x_i),D^[1]P(x_i))_(0<=i<s).                  (3)
```

Let

```text
alpha_1<=...<=alpha_(2s)                              (4)
```

be the p-adic valuations of its Smith invariant factors. Then their multiset
is exactly

```text
{alpha_i}
 =e*{0,0,1,2,...,s-1,s+1,s+2,...,2s-1}.             (5)
```

Empty ranges are omitted, so `s=1` gives `(0,0)` and `s=2` gives
`(0,0,e,3e)`. The latter recovers THM-4010's odd-prime two-node law at a gap
of exact valuation `e`.

## 2. Weighted minors

Integer translation acts on coefficients by a binomial unitriangular matrix,
so take `a=0`. In monomial degree `q`, the row of Hasse order `r in {0,1}`
has entry

```text
binom(q,r)(p^e u_i)^(q-r).                            (6)
```

Consider a nonzero `h` by `h` minor. Let its column degrees be

```text
q_1<...<q_h
```

and let `d` of its rows have order one. Multiplying those derivative rows by
`p^e` and extracting `p^(e q_j)` from column `q_j` gives

```text
v_p(minor)
 =e*(sum_j q_j-d)+v_p(residual minor),                (7)
```

where the residual rows are the value and derivative evaluations

```text
u_i^q,                 q u_i^(q-1)                    (8)
```

over `F_p`.

If `h<=s` and not every selected row is a derivative, then `d<=h-1` and

```text
sum_j q_j-d >= h(h-1)/2-(h-1)
                 =(h-1)(h-2)/2.                      (9)
```

If every row is a derivative, column zero vanishes; a nonzero minor must use
degrees at least `1,...,h`, which gives at least the same lower bound (and a
strictly larger one for `h>=2`). If
`h>s`, at most `s` derivative rows exist, so

```text
sum_j q_j-d >= h(h-1)/2-s.                           (10)
```

Thus, writing `Delta_h` for the `h`th determinantal divisor and
`v_p(Delta_0)=0`,

```text
v_p(Delta_h) >= e*(h-1)(h-2)/2,          1<=h<=s,
v_p(Delta_h) >= e*(h(h-1)/2-s),          s<h<=2s.     (11)
```

## 3. The residual matroid attains every lower bound

For `h<=s`, use degrees `0,...,h-1`, one value row at `u_0`, and derivative
rows at `u_0,...,u_(h-2)`. Translate the residual coordinates so `u_0=0`.
After the value row selects the constant column, the residual determinant is

```text
(h-1)! product_(0<=i<j<h-1)(u_j-u_i),                (12)
```

a unit modulo `p` because `h<=s<p`. This attains the first line of `(11)`.

Now let `s<h<=2s` and again use degrees `0,...,h-1`. The `s` derivative rows
have rank `s`: columns `1,...,s` give `s!` times a Vandermonde unit. The full
collection of value and derivative rows has rank `h` on polynomials of degree
below `h`. Indeed, a polynomial killed by every value and derivative row has
all `u_i` as double roots, hence is divisible by

```text
product_i (X-u_i)^2                                  (13)
```

of degree `2s`; a polynomial of degree below `h<=2s` must therefore vanish.
Linear algebra now selects `h-s` value rows extending the derivative rows to
rank `h`. Their residual minor is a unit, attaining the second line of `(11)`.

Taking consecutive differences of the now-equalities in `(11)` gives

```text
0,0,e,2e,...,(s-1)e,(s+1)e,...,(2s-1)e,              (14)
```

which proves `(5)`. The missing exponent `se` is the exact signature of the
one value row needed before all derivative rows become available.

## 4. Consecutive-node quadratic range

Let the nodes be `0,1,...,n-1` and suppose

```text
1<=n<=p(p-1).                                        (15)
```

Write `n=qp+r`, `0<=r<p`. The residue classes modulo `p` contain `r`
clusters of size `q+1` and `p-r` clusters of size `q`, with empty clusters
omitted. Condition `(15)` makes every nonempty size strictly less than `p`.

THM-4010's p-local CRT decomposition splits the complete exponent partition
as the multiset union of the residue-cluster partitions: distinct residue
factors have unit resultants. Each nonempty cluster is

```text
c,p+c,2p+c,...,(s_c-1)p+c,                           (16)
```

so `(5)` with `e=1` gives the closed formula

```text
Alpha_p(n)=multiset_union_(c mod p)
 {0,0,1,...,s_c-1,s_c+1,...,2s_c-1}.                 (17)
```

For odd `p` and `p<n<=2p`, `(17)` specializes to THM-4010's pair band

```text
0^(2p), 1^(n-p), 3^(n-p).                            (18)
```

For `p>=5`, `(17)` continues far beyond that band, through a quadratic number
of consecutive nodes.

## 5. Sharp boundary, connection ledger, and scope

The condition `s<p` is not cosmetic. At `s=p`, the derivative witness uses
the factor `p!`; over `F_p` the derivative rows on degrees below `2p-1` have
rank `p-1`, not `p`. The first two exact hostiles are

```text
p=2: actual (0,0,2,2),       false extension (0,0,1,3),
p=3: actual (0,0,1,3,4,4),   false extension (0,0,1,2,4,5). (19)
```

The source is the integral two-jet lattice; the target is its p-primary
weighted residual matroid. Scaling and translation preserve determinantal
valuations. Passing to residue clusters forgets interactions inside a cluster
at the next p-adic scale; the exact-size sidecar `s<p` is what makes the
residual matroid uniform. The cheapest decisive hostile is `s=p`, where that
sidecar fails and the derivative rank drops.

This theorem determines every positive p-layer in the stated single-scale
and consecutive quadratic scopes. It does not give the partition for a
cluster of size at least `p`, the full Smith form for arbitrary `m,k`, or any
automatic Rule-30, Jacobian, factorial, or LRC consequence.

## 6. Exact audits

Reproduce the primary DVR path with

```bash
python3 04-computation/confluent_twojet_single_scale_smith_thm4080.py
python3 -O 04-computation/confluent_twojet_single_scale_smith_thm4080.py
```

It checks 132 single-scale matrices, 34 direct consecutive matrices, 178
formula rows, three sharp hostiles, translations, and two residual coordinate
systems in `1,068,749` exact gates.

The independent path uses all minors rather than DVR pivoting in its small
load-bearing cases, then reconstructs the proof witnesses over finite fields:

```bash
python3 04-computation/confluent_twojet_single_scale_smith_thm4080_independent_audit.py
python3 -O 04-computation/confluent_twojet_single_scale_smith_thm4080_independent_audit.py
```

It computes every determinantal divisor in six cases, checks 136 residual
matroids, 1,896 witness minors and 948 quadratic-range rows, and independently
locates the boundary rank drop in `13,437,566` exact gates. Normal and
optimized outputs byte-match their frozen files.
