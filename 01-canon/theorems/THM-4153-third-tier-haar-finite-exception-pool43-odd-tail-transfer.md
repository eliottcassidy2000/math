---
id: THM-4153
title: "Third-tier Haar/finite-exception pool-43 universal odd-tail LRC(14) transfer"
status: >
  PROVED THIRD-TIER HAAR + FINITE EXCEPTIONAL-RATIO/SCALE TRANSFER +
  VERIFIED-EXACT + INDEPENDENT INTERVAL AUDIT; LRC(14) OPEN. Haar measure
  at least 4/91 excludes every primitive odd tail ratio except an explicit
  eleven-ratio set. One positive safe component reduces those ratios to
  finite odd-scale banks. An explicit 43-label pool closes its full bank at
  one clock and supplies 5,752,004,349 eleven-body families, 5,751,751,249
  outside the current THM-4148 min/max width gate. Arbitrary bodies,
  parity-class entry, and LRC(14) remain OPEN.
source: codex-lrc14-planar-jc-breakthrough-20260825
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4152-second-tier-haar-finite-exception-pool40-odd-tail-transfer
related:
  - THM-4142-common-safe-arc-clock-pool-universal-odd-tail-lrc14-completion
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
script: 04-computation/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153.py
output: 05-knowledge/results/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153.out
independent_audit_script: 04-computation/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153_independent_audit.out
script_sha256: eb5b7d85b260686b372c3a92f352157ad33573efb1c0819ca5830174babd78ab
output_sha256: 29b72543b5ce18574ac03b498504c1a5e9e1c13c5cfb769a013b024ce64feaff
semantic_sha256: 8fc29d9da1111f0b2cf5e13c6cd6dbca6e68072acd7f4b95f5455f7e412245f7
independent_audit_script_sha256: cc65e3c37d3f63fc7cc5af843f666cd93f9911a544ef529e2259899d4535e6b9
independent_audit_output_sha256: 8f43757f05e7967c7d643bdd16b9abf4dffc3bf5d4c97426ee63559c9f5ea7ab
independent_semantic_sha256: fec0ab982776064029be96988978f783e365914e2615e03cc5bd0e1a17f3a998
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic proves the product-161 cutoff, evaluates
  all 117 primitive odd ratios below it, matches the Bernoulli formula to
  strict wall geometry on 1,053 ratios through q=101, reconstructs all
  exceptional widths and scale banks, and rebuilds 4,776 pool walls, 46
  safe components, the clock ledger, and family counts. Normal, optimized,
  and hash-seeded streams byte-match.
independent_audit: >
  ACCEPT. A no-import implementation constructs danger arcs and their
  doubled quotient image, intersects closed safe teeth successively, and
  independently reproduces the threshold universe, exception table, pool
  geometry, clock gaps, and hereditary counts. Normal, optimized, and
  hash-seeded streams byte-match.
---

# THM-4153 -- third-tier Haar plus finite exceptional ratios

**PROVED THIRD-TIER HAAR + FINITE EXCEPTIONAL-RATIO/SCALE TRANSFER +
VERIFIED-EXACT + INDEPENDENT INTERVAL AUDIT; LRC(14) REMAINS OPEN.**

THM-2061 gives the folded two-tail obstruction. THM-4150 proves the strict
`4/63` Haar localization, and THM-4152 lowers the working level to `4/77`
by retaining one component and one exceptional ratio. Here the scalar level
is lowered again, while all discarded labels are restored before clocks are
used.

## 1. Reusable criterion

For a nonempty finite positive body `H`, put

```text
G_H={y in R/Z:min_(h in H)||hy||>=1/14}.                (1)
```

For coprime odd `0<p<q`, let `C_(p,q)` be THM-4136's open quotient set on
which both half-lifts are tail-bad. Define the labelled exception set

```text
E={(1,9),(1,11),(3,11),(1,23),(5,11),(1,37),
   (1,25),(3,25),(3,23),(1,51),(5,9)}.                 (2)
```

Its exact measures and largest open-component widths are:

| primitive ratio | `mu(C_(p,q))` | `beta_(p,q)` |
|---:|---:|---:|
| `(1,9)` | `4/63` | `2/63` |
| `(1,11),(3,11)` | `4/77` | `2/77` |
| `(1,23)` | `8/161` | `2/161` |
| `(5,11)` | `18/385` | `9/385` |
| `(1,37)` | `12/259` | `2/259` |
| `(1,25),(3,25)` | `8/175` | `2/175` |
| `(3,23)` | `22/483` | `2/161` |
| `(1,51)` | `16/357` | `2/357` |
| `(5,9)` | `2/45` | `1/45` |

Suppose `mu(G_H)>=4/91`. Positive measure gives a positive-length component;
let `L` be the largest component length. For `(p,q) in E` set

```text
T_H(p,q)={t positive odd:tL<beta_(p,q)}.                (3)
```

Assume that for every `(p,q) in E` and `t in T_H(p,q)` there are
`y_(p,q,t) in G_H` and `epsilon_(p,q,t) in {0,1}` such that, for
`x=(y_(p,q,t)+epsilon_(p,q,t))/2`,

```text
||ptx||>=1/14,                         ||qtx||>=1/14.    (4)
```

> **Third-tier transfer.** These hypotheses imply that every distinct
> positive odd pair `a,b` admits `x in R/Z` with
>
> ```text
> min_(v in 2H union {a,b})||vx||>=1/14.                (5)
> ```

## 2. Why only the eleven ratios survive

THM-4150's exact Bernoulli formula is

```text
mu(C_(p,q))
 =2/49+2[B_2(u_-)-B_2(u_+)]/(pq),                      (6)
u_-={1/2+(q-p)/14},             u_+={1/2+(q+p)/14},
B_2(u)=u^2-u+1/6.
```

Since the oscillation of `B_2` is `1/4`,

```text
mu(C_(p,q))<=2/49+1/(2pq).                              (7)
```

Now

```text
4/91-2/49=2/637,              1/(2*161)<2/637.          (8)
```

Thus every odd product `pq>=161` is strictly below `4/91`. The complete
primitive odd universe `pq<160` has `117` elements. Exact substitution in
`(6)` gives

```text
mu(C_(p,q))>4/91  iff  (p,q) in E,                      (9)
```

with equality exactly at

```text
(1,13),(3,13),(5,13),(1,39),(1,65).                    (10)
```

Strict-wall enumeration yields the component maxima in Section 1. The
independent audit obtains them instead by intersecting physical danger arcs
and doubling the resulting interval union.

## 3. Proof

Write distinct odd tails, after interchange, as

```text
a=pt,             b=qt,             t=gcd(a,b),
0<p<q,            gcd(p,q)=1.                            (11)
```

If `(5)` failed, every `y in G_H` would have both physical half-lifts
tail-bad, hence

```text
G_H subset m_t^(-1)(C_(p,q)),             m_t(y)=ty.   (12)
```

Multiplication by `t` preserves Haar measure. If `(p,q) notin E`, equations
`(9),(10),(12)` force equality of measures at worst. But `G_H` is compact
and nonempty, whereas the right side is open and proper because
`0 notin C_(p,q)`. Equality is impossible: a nonempty open difference has
positive measure, while an empty difference would be a nontrivial clopen
subset of the circle. This also closes the five equality ratios.

Take `(p,q) in E` and a closed component `J subset G_H` of length `L`. If
`tL>=1`, its circular image is the whole circle. If `tL<1` but the image
wraps, it contains zero. Otherwise it is a compact interval of length `tL`.
Containment in `C_(p,q)` would then force the strict inequality

```text
tL<beta_(p,q).                                          (13)
```

Thus width closes every scale outside `T_H(p,q)`, including equality. For a
scale inside the bank, `(4)` keeps both tails safe and `y_(p,q,t) in G_H`
keeps the doubled body safe. This contradicts failure and proves `(5)`.
**QED.**

## 4. Explicit pool and one-clock closure

Put

```text
P={1,2,4,5,8,10,16,17,19,20,23,25,29,31,32,34,38,40,
   41,43,47,50,51,53,58,62,64,67,69,71,73,75,76,80,
   82,86,89,93,95,111,141,159,285}.                    (14)
```

This is THM-4152's pool plus `111,159,285`. Its exact safe-set arrangement
has `4,776` walls and `46` positive-length components, with

```text
mu(G_P)
 =10080921463555906580413/211196778145191767531400
 =4/91
  +10368105800402918384569/2745558115887492977908200
 >4/91.                                                 (15)
```

The longest reflected pair is

```text
[911/3990,307/1330],      [1023/1330,3079/3990],
L=1/399.                                                 (16)
```

The residual odd-scale banks are:

| ratio | scales |
|---:|:---|
| `(1,9)` | `1,3,5,7,9,11` |
| `(1,11),(3,11)` | `1,3,5,7,9` |
| `(1,23),(1,37),(1,25),(3,25),(3,23)` | `1,3` |
| `(5,11)` | `1,3,5,7,9` |
| `(1,51)` | `1` |
| `(5,9)` | `1,3,5,7` |

Every row is closed by the single upper lift

```text
y=57/742 in G_P,                x=(y+1)/2=799/1484.     (17)
```

The doubled-body clearance is exactly `1/14`, uniquely owned by speed `53`.
The smallest tail gap across the entire bank is

```text
113/1484=1/14+1/212>1/14.                               (18)
```

Therefore `2P union {a,b}` is safe for every distinct positive odd pair.
Deleting body labels preserves the certificate.

## 5. Counts, boundary, and loss ledger

The pool supplies

```text
binom(43,11)=5,752,004,349                              (19)
```

eleven-body families. Since the 40-label pool is literally contained in
`P`, the increment is `3,440,202,909`. Inside precisely the universe of
eleven-subsets of `P`, the THM-4148 inequality

```text
27(13 min H-max H)>=4 min H max H                       (20)
```

holds for `253,100` bodies. Hence `5,751,751,249` members of this stated
universe fail that width gate; no claim about every other certificate is
made.

The `4/91` boundary is inclusive because compact-to-proper-open containment
closes equality. Below it, the five ratios in `(10)` become genuine
exceptions. The present clock is hostile there: for `(5,13)` it fails at
every residual scale `1,3,5,7`. A fourth tier therefore needs another clock
or a different pool/clock pair.

```text
source:       complete closed body-safe set G_H
target:       every doubled-body row with two distinct odd tails
map:          Haar mass -> ratio filter; component width -> scale filter
preserved:    total mass, strict/open type, labelled ratio, width, clocks
destroyed:    component addresses and owners at the scalar stage
sidecar:      one addressed component plus E-indexed clocks
positive:     P closes all E-banks with one clock
hostile:      that clock fails the next (5,13) bank
next test:    two-clock covering of the five equality-ratio banks.        (21)
```

No pool maximality or threshold optimality is claimed. The theorem does not
handle mixed/even tails, physical entry into the parity seam, arbitrary
bodies, or LRC(14).

## 6. Exact replay

```text
python3 -B 04-computation/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153.py
python3 -B -O 04-computation/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153.py
PYTHONHASHSEED=271828 python3 -B 04-computation/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153.py

python3 -B 04-computation/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153_independent_audit.py
python3 -B -O 04-computation/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153_independent_audit.py
PYTHONHASHSEED=314159 python3 -B 04-computation/lrc14_third_tier_haar_finite_exception_pool43_transfer_thm4153_independent_audit.py
```

All streams byte-match their frozen outputs. **QED.**
