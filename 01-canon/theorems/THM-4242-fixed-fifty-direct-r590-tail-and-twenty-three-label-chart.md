---
id: THM-4242
title: "Fixed-fifty direct r590 tail, twenty-three-label chart, and gcd-fibre tariff"
status: >
  PROVED RELATIVE TO THM-4150/4156/4170/4191/4211/4229/4234/4240 +
  FINITE-EXACT + INDEPENDENTLY AUDITED; LRC(14) OPEN. Every nine-body in
  the thirty-label pool is safe after adjoining 50 and every r>=590. The
  selected five-petal set {10,15,20,30,264} has every triple, quadruple,
  and quintuple finite head exhausted, proving chi_50>=23. A separate
  gcd-weighted fibre tariff is proved for every modulus h>=2; its prime p<=7
  corollary and p=5 census are only sufficient certificates. Arbitrary pair
  entry and LRC(14) remain OPEN.
source: codex-lrc14-breakthrough-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
  - THM-4229-fixed-fifty-nineteen-label-petal-haar-charts
  - THM-4234-fixed-fifty-twenty-label-pair-haar-charts
  - THM-4240-fixed-fifty-twenty-two-label-four-petal-haar-chart
related:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4238-forty-small-label-uniform-r590-haar-tail-closure
atlas_primary_script: 04-computation/lrc14_fixed_fifty_four_petal_haar_chart_thm4240.cpp
atlas_independent_script: 04-computation/lrc14_fixed_fifty_high_petal_batched_independent_audit_thm4242.cpp
bridge_primary_script: 04-computation/lrc14_fixed_fifty_r590_bridge_primary_thm4242.cpp
bridge_independent_script: 04-computation/lrc14_fixed_fifty_r590_bridge_independent_audit_thm4242.cpp
chart_primary_script: 04-computation/lrc14_fixed_fifty_five_petal_chart_primary_thm4242.cpp
chart_independent_script: 04-computation/lrc14_fixed_fifty_five_petal_chart_independent_audit_thm4242.cpp
fibre_primary_script: 04-computation/lrc14_five_fiber_certificate_census_thm4242.cpp
fibre_independent_script: 04-computation/lrc14_five_fiber_certificate_census_independent_audit_thm4242.cpp
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. Independent grouped-event, dual core-first, complement-scatter,
  literal midpoint, and cached divided-system census implementations reproduce
  every consequence-bearing count. The r590 bridge agrees on all 24 exceptional
  bodies and checks every one of their 258 literal rows by two independent
  integrators. The selected chart agrees on all 16 higher faces and all
  88,007,232 literal comparisons. No computation uses floating point or sampling.
---

# THM-4242 -- fixed-fifty direct r590 tail, twenty-three-label chart, and gcd-fibre tariff

**PROVED RELATIVE TO THM-4150/4156/4170/4191/4211/4229/4234/4240 +
FINITE-EXACT + INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

For a finite positive set `A`, write

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14},
alpha=4/63.                                             (1)
```

Retain the split pool

```text
C={8,16,40,42,80,84,85,88,95,
   120,126,143,145,168,193,240,252,286},

O={10,15,20,30,60,63,132,170,176,190,264,290},
P=C union O.                                            (2)
```

Then the following statements hold.

> **Direct fixed-fifty tail.** For every `B in binom(P,9)` and every integer
> `r>=590`,
>
> ```text
> mu(G_(B union {50,r}))>=alpha.                        (3)
> ```
>
> **Universal twenty-three-label chart.** Put
>
> ```text
> F={10,15,20,30,264},       S=C union F.               (4)
> ```
>
> For every positive integer `r notin P union {50}` and every
> `K in binom(S,9)`,
>
> ```text
> mu(G_(K union {50,r}))>=alpha.                        (5)
> ```

Consequently the fixed-fifty chart number defined in THM-4211 satisfies

```text
chi_50>=23.                                             (6)
```

The closest proved mechanism is the fixed-fifty component-discrepancy atlas
of THM-4234 and THM-4240. The canonical hostile remains THM-4207's warning
that individually safe newcomer marginals do not compose. The corrected near
miss is MISTAKE-520: zero-newcomer faces are inherited from full-pool heredity,
not from an unrelated cardinality certificate. The least-used sidecar here is
the component count of the unprojected circular safe set.

The live concept board was

```text
fixed-fifty hypergraph | finite heads | component discrepancy
split core/petal | gcd fibres | owner/tail boundary.    (7)
```

## 2. The complete fixed-fifty body atlas

For `B in binom(P,9)`, put

```text
V_B=G_(B union {50}),
mu(V_B)=m_B/D,
D=91,205,797,082,400,
c_B=# positive-length circular components of V_B.      (8)
```

THM-4170 gives, for every positive integer `r`,

```text
mu(V_B intersect G_r)
 >=(6/7)(m_B/D)-6c_B/(49r).                            (9)
```

Define

```text
delta_B=54m_B-4D,
kappa(B)=ceil(54c_B D/(7delta_B))                       (10)
```

when `delta_B>0`. Equations `(9)--(10)` imply `(3)` for every
`r>=kappa(B)`. This is a sufficient analytic cutoff, not necessarily the
first literal safe integer.

Let `k=|B intersect O|`. Exact integration gives the complete atlas:

| `k` | bodies | strict bodies | largest `kappa` | source |
|---:|---:|---:|---:|:---|
| 0 | 48,620 | 48,620 | 448 | THM-4211 |
| 1 | 525,096 | 525,096 | 502 | THM-4229 |
| 2 | 2,100,384 | 2,100,384 | 563 | THM-4234 |
| 3 | 4,084,080 | 4,084,080 | 589 | THM-4234 |
| 4 | 4,241,160 | 4,241,160 | 608 | THM-4240 |
| 5 | 2,423,520 | 2,423,520 | 626 | new exact atlas |
| 6 | 753,984 | 753,984 | 596 | new exact atlas |
| 7 | 121,176 | 121,176 | 517 | new exact atlas |
| 8 | 8,910 | 8,910 | 433 | new exact atlas |
| 9 | 220 | 220 | 336 | new exact atlas |

The layer count is lossless:

```text
sum_(k=0)^9 binom(12,k)binom(18,9-k)
 =14,307,150=binom(30,9).                              (11)
```

Thus every fixed-fifty body is strict. The analytic worst body is

```text
B_0={170,176,190,193,240,252,264,286,290},
m_(B_0)=18,257,430,123,108,
c_(B_0)=552,
delta_(B_0)=621,078,038,318,232,
kappa(B_0)=626.                                        (12)
```

## 3. Closing the analytic gap down to 590

The layers `k<=3` already have `kappa<=589`, while the layers `k>=7` have
`kappa<=517`. Only `k=4,5,6` can miss the desired boundary. Among their

```text
4,241,160+2,423,520+753,984=7,418,664                 (13)
```

bodies, exactly twenty-four have `kappa(B)>590`:

| layer | exceptional bodies | rows `590<=r<kappa(B)` |
|---:|---:|---:|
| 4 | 10 | 93 |
| 5 | 13 | 159 |
| 6 | 1 | 6 |
| **total** | **24** | **258** |

Exact literal integration of all 258 rows gives

```text
failures=0, equalities=0.                              (14)
```

The closest row is

```text
B={20,170,176,190,193,240,252,286,290}, r=592,

mu(G_(B union {50,592}))
 =523,724,171,982,358 / 3,374,614,492,048,800,

mu(G_(B union {50,592}))-4/63
 =51,577,155,731,993 / 562,435,748,674,800>0.          (15)
```

The primary and independent programs agree on the body, mass, component
count, surplus, and cutoff of every one of the twenty-four exceptions. The
independent program then integrates every literal row twice: once by grouped
endpoint toggles and once by fresh midpoint cells. All `258` paired values
agree. Together with `(9)` outside the finite gaps, this proves `(3)`.

This lowers the direct fixed-`50` boundary from the raw analytic value `626`
to the same `590` boundary used by THM-4238's genuine small-first-outsider
theorem. It does not address the remaining first-owner rays `51<=q<=1289`.

## 4. A universal five-petal chart

For a nonempty `T subset F` with `|T|=j`, a nine-body using exactly `T`
has the form

```text
L union T,       L in binom(C,9-j).                    (16)
```

The layers `j=0,1,2` are already closed by THM-4211, THM-4229, and THM-4234.
The new computation exhausts all

```text
binom(5,3)+binom(5,4)+binom(5,5)=10+5+1=16            (17)
```

triple, quadruple, and quintuple faces. Every limiting body is strict; every
permitted outsider below its analytic cutoff is then integrated literally:

| face size | faces | limiting profiles | literal rows | literal checks |
|---:|---:|---:|---:|---:|
| 3 | 10 | 185,640 | 3,884 | 72,102,576 |
| 4 | 5 | 42,840 | 1,742 | 14,925,456 |
| 5 | 1 | 3,060 | 320 | 979,200 |
| **total** | **16** | **231,540** | **5,946** | **88,007,232** |

There are zero literal failures and zero equalities. The referee directly
replays two limiting extremizers and one closest literal body for each face,
for `32+16` body-local controls, and agrees on all consequence-bearing fields.

The internal layer identity is

```text
sum_(j=0)^5 binom(5,j)binom(18,9-j)
 =48,620+218,790+318,240+185,640+42,840+3,060
 =817,190=binom(23,9).                                 (18)
```

For any permitted `r`, every eleven-subset of `S union {50,r}` is safe:
THM-4156 handles subsets containing neither newcomer, THM-4191 handles those
containing exactly one, and `(5)` handles those containing both. The disjoint
Pascal count is

```text
binom(23,11)+2binom(23,10)+binom(23,9)
 =1,352,078+2,288,132+817,190
 =4,457,400=binom(25,11).                              (19)
```

Equations `(16)--(19)` prove `(5)--(6)`. THM-4150 supplies universal distinct
odd-tail completions after common positive scaling; it does not supply entry
of an arbitrary LRC(14) speed row.

## 5. Gcd-fibre tariff

The finite chart exposed a reusable operation that is independent of the
fixed-fifty atlas.

> **Gcd-fibre lemma.** Let `h>=2` be an integer, and let `D,E` be finite sets
> of positive integers such that every `d in D` is divisible by `h`. For
> `e in E`, define
>
> ```text
> g_e=gcd(e,h),
> w_h(e)=g_e ceil(h/(7g_e)).
> ```
>
> Then
>
> ```text
> mu(G_(D union E))
>  >=((h-sum_(e in E)w_h(e))_+/h) mu(G_(D/h)),         (20)
> ```
>
> where `D/h={d/h:d in D}` and `G_empty=R/Z`.

**Proof.** Disintegrate Haar measure along the fibres of multiplication by
`h`: write the `h` preimages of `y` as `x_j=(y+j)/h`, `0<=j<h`. For `d=hd'`,

```text
||d x_j||=||d'y||,
```

so all constraints from `D` are fibre-constant. For fixed `e`, the residues
`ej mod h` form `h/g_e` equally spaced phases, each with multiplicity `g_e`.
An open circular arc of length `1/7` contains at most
`ceil((h/g_e)/7)` points of this grid. Thus the unsafe condition
`||e x_j||<1/14` removes at most `w_h(e)` fibre indices. The union bound leaves
at least `(h-sum_e w_h(e))_+` indices over every `y in G_(D/h)`. Conditional
uniformity of Haar measure on the `h`-point fibres proves `(20)`. When
`7` divides `h/g_e`, openness at both endpoints prevents an extra grid point;
this boundary convention is load-bearing. **QED.**

For prime `p<=7` and every `e` not divisible by `p`, one has `g_e=1` and
`w_p(e)=1`, so `(20)` specializes to

```text
mu(G_(D union E))>=((p-|E|)_+/p)mu(G_(D/p)).           (21)
```

The fibre weight itself is sharp. If `D` is empty, `E={e}`, and
`7` divides `h/g_e`, then both sides of `(20)` equal `6/7`; for example take
`h=14,e=2`. For a hostile to the naive unweighted extension, take `h=8,e=2`:
the false one-phase coefficient is `7/8`, while `mu(G_{2})=6/7`; `(20)` uses
the required gcd multiplicity. A tariff is informative only when its total
weight is less than `h`.

As a diagnostic, the `p=5` lemma was applied to the eleven higher faces of
`F` that contain `264`. It certifies exactly

```text
12,446,627 of 59,630,016                              (22)
```

literal comparisons. The remaining `47,183,389` comparisons are merely
**uncertified by this lemma**, not failures. The weakest certificate surplus is
`53/9,317,700` (the program's scaled gap is `315` times this), with four
nondivisible constraints and divided system
`{2,3,10,17,24,29,54}`. The closest literal outsiders left to the direct
engine are the common-factor-`5` resonances `65`, `70`, and `110`. Thus the fibre view
explains a real fraction of the finite chart but is not equivalent to it.

## 6. Independent audit architecture

The primary atlas constructs the exact common wall grid, evaluates midpoint
safety on open cells, records core failure masks, and applies subset-zeta
transforms over petal masks. The higher-petal referee instead groups
simultaneous wall events, fixes the core first, and uses complement scatter.
Across `k=5,...,9`, the two programs agree on all

```text
792+924+792+495+220=3,223                             (23)
```

labelled atlas lanes; all displayed masses, components, surpluses, cutoffs,
and extremizers agree, modulo harmless tied traversal where applicable.

The bridge primary is a small body-local literal integrator driven by the
analytic atlas. The bridge referee independently reconstructs the exceptional
set by dual core-first/O-superset zeta and checks literal rows by both grouped
events and midpoint cells. Neither engine imports the other's exception list.

The chart primary batches all sixteen higher faces on one exact grid. Its
referee replaces midpoint cells by grouped endpoint toggles and directly
replays every minimum, cutoff, and closest row. The fibre primary and referee use separate cached divided
systems and reproduce every lane count in `(22)`.

## 7. Reproduction, hashes, and scope firewall

From the repository root, representative commands are:

```bash
g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_r590_bridge_primary_thm4242.cpp \
  -o /tmp/lrc14-thm4242-r590-primary
/tmp/lrc14-thm4242-r590-primary

g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_r590_bridge_independent_audit_thm4242.cpp \
  -o /tmp/lrc14-thm4242-r590-independent
/tmp/lrc14-thm4242-r590-independent

g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_five_petal_chart_primary_thm4242.cpp \
  -o /tmp/lrc14-thm4242-chart-primary
/tmp/lrc14-thm4242-chart-primary --selected-five-petal

g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_five_petal_chart_independent_audit_thm4242.cpp \
  -o /tmp/lrc14-thm4242-chart-independent
/tmp/lrc14-thm4242-chart-independent
```

The higher-petal atlas is reproduced by running the THM-4240 primary with
`--all-petal-base k` and the new independent program with argument `k`, for
each `k=5,...,9`. The two fibre programs take no arguments. Compare each
stream with the correspondingly named frozen output in `05-knowledge/results`.

| artifact | SHA-256 |
|:---|:---|
| higher-petal independent source | `60e565bec0654c41fb97d46fd79652720f3917b18350a787c9dcf589b917722c` |
| r590 bridge primary source | `93cad190f7bec7fe1e439829a530dec752989a392df6b6d5ee5e2e88aebc6d5e` |
| r590 bridge independent source | `b0ed4b93396390b5d94bb6fd2d93ad486c90a43fab56fa1d2b2fa4cec87f3980` |
| r590 bridge primary output | `6a88fa40107d38c610366048dfb96a02f823dab246b8dd1f097d8723bfa50ff5` |
| r590 bridge independent output | `f567996fe0e62520b176ac768e657b80574ac77b5c717bbd591d2962d03dc869` |
| chart primary source | `7a837c03703f259180f774d84cfe754d2e91cd7af1ca3ee9d101740ad6fb8d05` |
| chart independent source | `0c53f7856d9f34e1622df6eb5b72ba013b7edd3931007894e8ffb813ae7b9a0d` |
| chart primary output | `14d4f971c7cf42ab22ef3ce709526549adaf5488d6520fdc5bef98db4e899f74` |
| chart independent output | `d620e1ab9aba91eaa6dd7da8e78a9be7ae2ee8624c5acc7857ddd34f87d733aa` |
| fibre primary source | `f207c002bbd7b5e6ea30c8fada0a272eb41709368f006722642d2ad1291c71c7` |
| fibre independent source | `76492bbcce17868f24d93da7772aa7e5b5572e9e712f2682db893c7cb868292e` |
| fibre primary output | `268c483f7cdf35d03e5660d840332954cf5a5460774ef0cbaff751d18aa98ce1` |
| fibre independent output | `6d976c65a7b5ef7d4f95c9d2b91f185d16fee97ef4f5e711d35b8f02036494e8` |

The primary/independent output hash pairs for `k=5,...,9` are respectively

```text
k=5  b07c339a605549506043f964633f337e20053d353833d6fa54078d4082e81e58
     d802adb16483d26777f86ede8027e5ccb9fce890a3b75113a4c64712ed6501a5
k=6  e33fc8a557d28d79a3a7294918aa17bb7c08a9e0d572ed20fb01e9650711c041
     6672d99a45568d1c499360d53df15088cda5bd6b78da516c7cc4b16973c67ab7
k=7  954a2187132c172c8d579fbbd7f0c5df9ce7e2c201b1c1f52be4f613438b11da
     51ba4dcbcc0a7bdb141eceb40fb5955a0e4d595145358a89f9652d45179543f3
k=8  27334408c9cf0d52ea58d046736e146e2166e4abb82ff562d9e359470b5548ed
     1aa6a3bcb90f8de93b95aaff9e6f5835386fa42e35b10d104d3f36feaac5b051
k=9  b3efeaae3347503e8b3cb53cc7398cf4acd036bef195bea05bfa8e1559037ab9
     d9c237e0628f8641aad7bede5cdec059e8e29c9b9b2308eddfe0526b015b3675
```

The pairs differ because line formats and traversal orders are deliberately
independent; all normalized numerical fields agree. No program uses floating
point or sampling.

This theorem does **not** prove that `590` is the minimal direct tail, that
`chi_50=23`, replacement of the fixed center `50`, arbitrary pair entry,
physical entry of an arbitrary speed row, or LRC(14). **QED.**
