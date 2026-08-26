---
id: THM-4152
title: "Second-tier Haar and finite-exception transfer for universal odd-tail completion"
status: >
  PROVED SECOND-TIER HAAR + FINITE EXCEPTIONAL-SCALE TRANSFER +
  VERIFIED-EXACT + INDEPENDENT IMPLEMENTATION AUDIT; LRC(14) OPEN. If the
  complete body-safe set has Haar measure at least 4/77, every primitive odd
  tail ratio except (1,9) is automatically excluded; one largest safe
  component reduces (t,9t) to finitely many odd scales, which any exact clock
  bank may close. The explicit 40-label pool needs only t=1,3,5,7 in this
  Haar/clock proof and has 2,311,801,440 eleven-body subfamilies. Exactly
  2,311,548,340 fail THM-4148's stated min/max width gate, but THM-4154
  proves that the concrete pool-family safety was already inherited from the
  small-denominator/divisor sieve. Arbitrary bodies, physical entry, and
  LRC(14) remain open.
source: codex-lrc-multiwindow-probe-20260825
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
related:
  - THM-4142-common-safe-arc-clock-pool-universal-odd-tail-lrc14-completion
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  - THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction
script: 04-computation/lrc14_second_tier_haar_finite_exception_pool40_transfer.py
output: 05-knowledge/results/lrc14_second_tier_haar_finite_exception_pool40_transfer.out
independent_audit_script: 04-computation/lrc14_second_tier_haar_finite_exception_pool40_transfer_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_second_tier_haar_finite_exception_pool40_transfer_independent_audit.out
script_sha256: f1c59b205c258c6127c595d04161bb49633ee3a59977e930953225a635de074b
output_sha256: 7f4c77741f74488e4c413056ea548351d62adaa123ac9f81831d8aebca63d9b3
semantic_sha256: f13989721e4bddfbfea0b48abdaba9201e06daed15e7cbf7ad139084ce00168c
independent_audit_script_sha256: a504d3486eadc40adbb55124a97ea5875ddeab2c18831a7fdb7ea2332bb3f20f
independent_audit_output_sha256: f9ccfec7191cd9d78928b455344f5951e1c0d65eebc3ef9224ba6480c5365d05
independent_semantic_sha256: ce2f02310458b21dbc247590ff969065f64008ef4d93241dab052118b09553ee
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic derives the product-45 cutoff, evaluates
  the complete 26-ratio residue universe, and matches the Bernoulli formula
  to independently rebuilt strict wall cells on all 1,053 primitive ratios
  through q=101. It reconstructs 3,744 body walls, 44 components, the exact
  measure and width, the four-scale clock, and all family counts. Normal,
  optimized, and hash-seeded streams byte-match the frozen output.
independent_audit: >
  ACCEPT. A no-import direct periodic-tooth implementation checks 2,350
  primitive ratios through q=151, intersects closed safe teeth successively,
  reproduces both longest components and every clock gap, and recomputes the
  hereditary and width comparison counts. Normal, optimized, and hash-seeded
  streams byte-match its frozen output.
---

# THM-4152 -- second-tier Haar plus finite exceptional scales

**PROVED SECOND-TIER HAAR + FINITE EXCEPTIONAL-SCALE TRANSFER +
VERIFIED-EXACT + INDEPENDENT IMPLEMENTATION AUDIT; LRC(14) REMAINS OPEN.**

THM-2061 already gives the exact folded strict obstruction and the sharp
`4/63` measure cap. THM-4150 supplies a Fourier--Bernoulli dual derivation
and upgrades equality by compact-closed versus proper-open containment. The
next operation is deliberately two-stage: lower the Haar threshold until
one primitive ratio survives, then restore one component address and a
finite clock bank only for that ratio.

## 1. Reusable two-stage criterion

Let `H` be a nonempty finite set of positive integers. Put

```text
G_H={y in R/Z:min_(h in H)||hy||>=1/14},                 (1)
```

and write `mu` for normalized Haar measure. Suppose

```text
mu(G_H)>=4/77.                                           (2)
```

Because `(2)` is positive, `G_H` has a positive-length component. Let `L`
be the largest component length and define the finite set

```text
T_H={t positive odd:tL<2/63}.                            (3)
```

Assume that for every `t in T_H` there are `y_t in G_H` and
`epsilon_t in {0,1}` for which the physical lift

```text
x_t=(y_t+epsilon_t)/2                                   (4)
```

satisfies

```text
||t x_t||>=1/14,                    ||9t x_t||>=1/14.    (5)
```

> **Two-stage theorem.** Under `(2)--(5)`, for every pair of distinct
> positive odd integers `a,b`, there is an `x in R/Z` such that
>
> ```text
> min_(v in 2H union {a,b})||vx||>=1/14.                 (6)
> ```

Thus Haar measure closes all but one primitive ratio, one component makes
the remaining scale bank finite, and clocks are required only on that bank.

## 2. The sharp second cross-comb level

For coprime odd `0<p<q`, let `C_(p,q)` be THM-4136's open quotient set on
which both physical lifts are tail-bad. In THM-4150 notation,

```text
mu(C_(p,q))
 =2/49+2[B_2(u_-)-B_2(u_+)]/(pq),                       (7)
u_-={1/2+(q-p)/14},                 u_+={1/2+(q+p)/14},
B_2(u)=u^2-u+1/6.                                       (8)
```

The oscillation of `B_2` is `1/4`, so

```text
mu(C_(p,q))<=2/49+1/(2pq).                              (9)
```

Already at `pq=45`,

```text
2/49+1/90 < 4/77,                                      (10)
```

so only the complete coprime odd universe `pq<45` remains. Exact
substitution in `(7)` gives:

| `q` | admissible `p : mu(C_(p,q))` |
|---:|:---|
| `3` | `1:0` |
| `5` | `1:0, 3:2/105` |
| `7` | `1:2/49, 3:2/49, 5:2/49` |
| `9` | `1:4/63` |
| `11` | `1:4/77, 3:4/77` |
| `13` | `1:4/91, 3:4/91` |
| `15,17,19` | `1:4/105, 1:4/119, 1:4/133` |
| `21,23,25` | `1:2/49, 1:8/161, 1:8/175` |
| `27,29,31,33` | `1:8/189, 1:8/203, 1:8/217, 1:8/231` |
| `35,37,39` | `1:2/49, 1:12/259, 1:4/91` |
| `41,43` | `1:12/287, 1:12/301` |

These are all `26` pairs. Consequently

```text
mu(C_(p,q))<=4/77                    if (p,q)!=(1,9),    (11)
```

with equality exactly at `(1,11)` and `(3,11)`. The only larger value is

```text
mu(C_(1,9))=4/63.                                       (12)
```

The normal and independent exact implementations reproduce `(7)` from
strict wall cells and direct periodic-tooth intersections, respectively.

## 3. Proof of the two-stage theorem

Fix distinct odd tails and write, after interchange,

```text
a=pt,             b=qt,             t=gcd(a,b),
0<p<q,            gcd(p,q)=1.                              (13)
```

If `(6)` fails, then every `y in G_H` has both half-lifts tail-bad, so

```text
G_H subset m_t^(-1)(C_(p,q)),              m_t(y)=ty.    (14)
```

Multiplication by `t` preserves normalized Haar measure. If
`(p,q)!=(1,9)`, equations `(2)`, `(11)`, and `(14)` force equality of the
two measures. But the left side is compact and nonempty, while the right
side is a proper open set because `0 notin C_(p,q)`. As in THM-4150, the
open difference would either have positive measure or the common set would
be a nontrivial clopen subset of the circle. Both are impossible. Hence
every primitive ratio other than `(1,9)` is closed.

It remains to take `(p,q)=(1,9)`. Let `J` be a closed component of `G_H` of
length `L`. THM-4136 proves that every component of `C_(1,9)` has length at
most `2/63`. If `tL>=1`, the circular image `tJ` is the whole circle. If
`tL<1` but the image wraps, it contains zero. Otherwise it is a closed
interval of length `tL`; containment in one open component of `C_(1,9)`
would force

```text
tL<2/63.                                                 (15)
```

Thus every odd `t notin T_H` is closed by component width, including the
equality boundary. For `t in T_H`, `(4)--(5)` give a full-row-safe lift
because `y_t in G_H` already preserves every doubled-body speed. This
contradicts failure in every remaining case and proves `(6)`. **QED.**

## 4. Exact geometry of the 40-label pool; safety inherited

Put

```text
P={1,2,4,5,8,10,16,17,19,20,23,25,29,31,32,34,38,40,
   41,43,47,50,51,53,58,62,64,67,69,71,73,75,76,80,
   82,86,89,93,95,141}.                                 (16)
```

The THM-4150 subpool has no multiple of `6`, and the seven added labels have
residues

```text
67:1, 82:4, 86:2, 89:5, 93:3, 95:5, 141:3 mod 6.
```

Thus `P intersect 6Z=empty`. THM-4154's fixed phase `x=1/12` already closes
`2P union {a,b}` for every odd pair, with body clearance at least `1/6` and
tail clearance at least `1/12`. The following Haar geometry and finite bank
remain an exact independent realization of the reusable criterion, but they
do not create previously uncovered pool families.

Its complete exact safe-set arrangement has `3,744` walls and `44`
positive-length components, with

```text
mu(G_P)
 =23518182747542658511201/441420293060220631236800
 =4/77
  +6459842759986025773611/4855623223662426943604800
 >4/77.                                                 (17)
```

This lies strictly below `4/63`, so THM-4150's first-tier scalar gate does
not apply. The two longest components are the reflected pair

```text
J=[299/658,321/700],           1-J=[379/700,359/658],
L=137/32900.                                             (18)
```

Exact comparison gives

```text
7L<2/63,                      9L-2/63=1697/296100>0,    (19)
```

and hence

```text
T_P={1,3,5,7}.                                           (20)
```

One body-safe endpoint and its upper physical lift close the whole bank:

```text
y=57/742,                    x=(y+1)/2=799/1484.         (21)
```

The doubled-body clearance at `x` is exactly `1/14`, uniquely owned by body
speed `53`. The exact tail gaps are:

| `t` | `||tx||` | `||9tx||` |
|---:|---:|---:|
| `1` | `685/1484` | `229/1484` |
| `3` | `571/1484` | `687/1484` |
| `5` | `457/1484` | `339/1484` |
| `7` | `49/212` | `17/212` |

Every entry is strictly above `1/14`. Sections 1--3 therefore prove that
`2P union {a,b}` is `1/14`-safe for every distinct positive odd pair.
Deleting body speeds preserves safety, so every eleven-element `H subset P`
gives a thirteen-speed universal odd-tail family.

The label `141` is one exact successful extension of the 39-label
provenance pool. No maximality of `P` is claimed.

## 5. Family counts and comparison universes

The complete hereditary family count is

```text
binom(40,11)=2,311,801,440.                              (22)
```

This count is exact. Its safety significance is inherited: every one of
these bodies already has THM-4154's common `x=1/12` certificate.

The THM-4150 pool is a strict subset of `P`; its `193,536,720` bodies are
therefore inherited literally. The number of additional eleven-subsets is

```text
2,311,801,440-193,536,720=2,118,264,720.                 (23)
```

For the comparison with THM-4148's stated gate, the universe is **exactly all
eleven-subsets of `P`**. Fixing the minimum and maximum at positions `i<j`
leaves `binom(j-i-1,9)` choices. Summing over the pairs satisfying the
stated unrestricted width inequality

```text
27(13 min H-max H)>=4 min H max H                       (24)
```

gives `253,100`. Thus exactly

```text
2,311,801,440-253,100=2,311,548,340                     (25)
```

members of this stated universe fail THM-4148's min/max gate. This does not
claim they evade every other certificate.

## 6. Loss ledger, hostiles, and scope

The first stage intentionally scalarizes all component locations; the
second restores only the largest-component length and finitely many clocks.

```text
source:       complete closed body-safe set G_H
target:       every doubled-body row with two distinct odd tails
map:          Haar measure -> primitive-ratio filter; one component -> t-bank
preserved:    total safe mass, strict/open tail type, largest width, clocks
destroyed:    all component addresses and owners at the Haar stage
sidecar:      one addressed component and exact clocks for tL<2/63
hostile 1:    mu(G_P)<4/63, so the first-tier scalar test misses P
hostile 2:    7L<2/63, so width alone misses the last finite scale
decisive test: second cross-comb level, then the finite (1,9)-scale bank. (26)
```

Condition `(2)` by itself does not close `(1,9)`: the finite clock sidecar is
load-bearing for the **abstract** two-stage criterion. It is not load-bearing
for this particular pool, whose safety is inherited from THM-4154. The
theorem does not give a uniform `4/77` lower bound for arbitrary bodies,
claim pool maximality, handle mixed/even tails, or provide physical entry
into the `11+2` parity seam. LRC(14) remains open.

## 7. Exact replay

Replay the two disjoint implementations with

```text
python3 -B 04-computation/lrc14_second_tier_haar_finite_exception_pool40_transfer.py
python3 -B -O 04-computation/lrc14_second_tier_haar_finite_exception_pool40_transfer.py
PYTHONHASHSEED=314159 python3 -B 04-computation/lrc14_second_tier_haar_finite_exception_pool40_transfer.py

python3 -B 04-computation/lrc14_second_tier_haar_finite_exception_pool40_transfer_independent_audit.py
python3 -B -O 04-computation/lrc14_second_tier_haar_finite_exception_pool40_transfer_independent_audit.py
PYTHONHASHSEED=1618033 python3 -B 04-computation/lrc14_second_tier_haar_finite_exception_pool40_transfer_independent_audit.py
```

All streams byte-match their frozen outputs. **QED.**
