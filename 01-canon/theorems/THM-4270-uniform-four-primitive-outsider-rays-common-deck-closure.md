---
id: THM-4270
title: "Uniform four-primitive-outsider-ray common-deck closure"
status: >
  PROVED RELATIVE TO THM-4231/4256 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  For the fixed thirty-label pool and every labelled nine-body, the rays
  3:5, 7:8, 8:9, and 11:12 are Haar-safe at every strict-above-pool scale.
  One frozen 8,192-repair carrier supplies four ratio-specific common decks;
  each deck covers all 14,307,150 bodies throughout its finite bridge, while
  THM-4231 supplies the cofinite tail. In the exact earlier proof order after
  THM-4261/4262, the four-ray component contains 152 residual edges and leaves
  180,470; this is an order-relative component ledger, not the current
  post-THM-4266/4267 contribution. No other ratio, smaller scale, physical
  entry, or LRC(14) follows.
source: codex-creative-frontiers-20260827
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4256-uniform-two-three-outsider-ray-endpoint-cocycle-closure
  - THM-4261-semantic-endpoint-band-prefix-union-lift
  - THM-4262-uniform-three-four-outsider-ray-common-deck-closure
common_header: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_common.hpp
primary_script: 04-computation/lrc14_four_primitive_rays_common_deck_bridge_thm4270.cpp
primary_output: 05-knowledge/results/lrc14_four_primitive_rays_common_deck_bridge_thm4270.out
independent_script: 04-computation/lrc14_four_primitive_rays_literal_wall_audit_thm4270.cpp
independent_output: 05-knowledge/results/lrc14_four_primitive_rays_literal_wall_audit_thm4270.out
postprocess_script: 04-computation/lrc14_four_primitive_rays_residual_postprocess_thm4270.py
postprocess_output: 05-knowledge/results/lrc14_four_primitive_rays_residual_postprocess_thm4270.out
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary uses exact primitive-period prefix integrals.
  The clean-room path instead builds literal safe intervals on a fresh common
  denominator at every scale and streams the candidate order through a
  priority queue. O2/O3 and assertion-enabled replays agree after newline
  normalization. The paths agree on every common deck and body scan; their
  displayed raw least margins use different scale-dependent integer
  normalizations and are deliberately not compared as numerical invariants.
---

# THM-4270 -- uniform four-primitive-outsider-ray common-deck closure

**PROVED RELATIVE TO THM-4231/4256 + FINITE-EXACT + INDEPENDENTLY
AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

Put

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.
```

For every `B in binom(P,9)`, the following four statements hold:

```text
mu(G_(B union {3g,5g}))   >=4/63 for every integer g>=97;
mu(G_(B union {7g,8g}))   >=4/63 for every integer g>=42;
mu(G_(B union {8g,9g}))   >=4/63 for every integer g>=37;
mu(G_(B union {11g,12g})) >=4/63 for every integer g>=27.       (1)
```

Each lower endpoint is exactly the first scale with `ug>max(P)=290`. It is a
proof-domain boundary, not an asserted failure or optimum below that scale.
The finite bridge and inherited tail are

| ratio `u:v` | bridge audited here | first THM-4231 tail scale |
|---|---:|---:|
| `3:5` | `97..153` | `154` |
| `7:8` | `42..96` | `97` |
| `8:9` | `37..85` | `86` |
| `11:12` | `27..64` | `65` |

Indeed, at and beyond the last column `vg>=770`; THM-4231 closes every
distinct outsider pair whose maximum is at least `770`.

## 2. Primitive observable with the scale coordinate retained

For a primitive ratio `u<v`, let

```text
A_(u,v)=G_{u,v},       N_(u,v)=14*lcm(u,v).
```

The exact primitive wall arrangements used here have

| ratio | `N_(u,v)` | safe ticks in one period | components |
|---|---:|---:|---:|
| `3:5` | 210 | 156 | 7 |
| `7:8` | 784 | 576 | 12 |
| `8:9` | 1008 | 742 | 14 |
| `11:12` | 1848 | 1360 | 20 |

Under circle multiplication `m_g(x)=gx`,

```text
G_{ug,vg}=m_g^(-1)(A_(u,v)).                            (2)
```

Let

```text
D=lcm(14s:s in P)=18,241,159,416,480.
```

For `R in binom(P,8)`, write `U_R=G_(P\R)` as nonwrapping lifts
`[a_i/D,b_i/D]`. Define the exact periodic prefix integral

```text
I_(u,v)(z)=N_(u,v)D integral_0^(z/D) 1_(A_(u,v))(t)dt. (3)
```

Then

```text
mu(U_R intersect G_{ug,vg})
 =sum_i(I_(u,v)(gb_i)-I_(u,v)(ga_i))/(N_(u,v)Dg).       (4)
```

Thus `R` is active at scale `g` exactly when

```text
63 sum_i(I_(u,v)(gb_i)-I_(u,v)(ga_i))
 -4N_(u,v)Dg >=0.                                      (5)
```

Passing only to the primitive ratio destroys `g`, the endpoint phase, and
repair activation. Equation `(5)` is the necessary scale sidecar.

## 3. Exact common-deck bridges

Order all `binom(30,8)=5,852,925` repair masks by

```text
(SplitMix64(mask xor 0x4245422842334245),mask).         (6)
```

The first 8,192 candidates have FNV `60148ca1fc61dbcb`. For each ratio,
retain a candidate only when `(5)` is nonnegative at every bridge scale in
Section 1. The four ratio-specific common decks are:

| ratio | deck count / FNV | least primary numerator margin | scale | repair |
|---|---|---:|---:|---|
| `3:5` | `4178 / 7a9034824db27f98` | `10690982966160` | 124 | `{8,10,20,63,85,145,264,290}` |
| `7:8` | `4011 / 7ecafda9de695f6d` | `6029824356864` | 88 | `{10,63,80,120,132,168,193,290}` |
| `8:9` | `3046 / 35a89d2659e92442` | `146850483971520` | 40 | `{16,42,60,85,120,168,193,264}` |
| `11:12` | `4073 / 5c48ba1405aa1453` | `108871444725120` | 35 | `{8,42,63,80,85,126,132,176}` |

All displayed margins are strictly positive. Exhaustive scans of all
`binom(30,9)=14,307,150` labelled bodies give, in the same row order,

```text
ordered repair/body checks = 461733464, 468020899, 492537890, 464221590;
maximum checked prefixes    =      1835,      2905,      2334,      1495;
failures                    =         0,         0,         0,         0. (7)
```

For a disjoint body and repair, `B subset P\R`, hence

```text
G_((P\R) union {ug,vg}) subset G_(B union {ug,vg}).    (8)
```

An active disjoint repair therefore lower-bounds the target in the correct
direction. Equations `(5)--(8)` prove `(1)` on every finite bridge, and
THM-4231 supplies every cofinite tail.

## 4. Independent literal-wall audit

The clean-room checker never evaluates `(3)`. At each bridge scale it forms
the literal safe intervals for speeds `ug` and `vg` on

```text
L_g=lcm(D,14ug,14vg),                                  (9)
```

intersects the two sorted interval lists, and sweeps their joint arrangement
against all 7,133 fixed-pool atoms. It independently obtains the first 8,192
candidates by a bounded priority queue rather than sorting the full universe.

The load-bearing matrix controls are

| ratio | nonnegative cells | matrix FNV | literal-mass FNV | common count/FNV |
|---|---:|---|---|---|
| `3:5` | 414601 | `0a26e5c64b22308c` | `1740e07378870f0e` | `4178 / 7a9034824db27f98` |
| `7:8` | 386745 | `098e632eb8aca044` | `179003dca73e0a6e` | `4011 / 7ecafda9de695f6d` |
| `8:9` | 345371 | `893df3f01101f38b` | `79d8e1d07608a602` | `3046 / 35a89d2659e92442` |
| `11:12` | 261078 | `f135d4414de36aa6` | `0412791788a972c8` | `4073 / 5c48ba1405aa1453` |

It repeats all four complete body scans and reproduces the check counts and
zero failures in `(7)`. Its printed raw least margins need not select the same
row as the primary: `L_g` changes with `g`, so those integers use
scale-dependent normalizations. The invariant comparisons are every sign,
the four common-deck hashes, and the exact body scans; all agree.

## 5. Exact earlier-order component ledger

Start with the exact residual after THM-4261 and THM-4262:

```text
count=180622,
FNV=0cef4e2887c8f24e,
SHA256=fa1c5672b0f2cd2490413e9b69a4720
       bf1dc4eef8aee694c1c73d390aba58e11.               (10)
```

Semantic selection of the four proved bridge rays finds respectively

```text
32, 46, 43, 31                                        (11)
```

still-open edges. The four lists are pairwise disjoint. Their 152-edge union
has

```text
FNV=80ed45e7e179d8ff,
SHA256=46b1c987fe337d879a721efab99da8265
       d126d12798788295b0829f6bd5741fd.                 (12)
```

Deleting exactly that union leaves

```text
count=180470,
FNV=d9afdc10f8e2aa88,
SHA256=e2ba42307e4c628ea9ef517a858456e4
       dd64e0b7a034fd6dd74ff707a8838f3f,
maximum endpoint=732.                                 (13)
```

The postprocessor asserts every per-ray count/FNV/SHA, the disjoint union,
and the updated ledger. The theorem closes all bridge scales in `(1)`; only
152 were new in the declared post-THM-4261/4262 proof order because those two
earlier theorems already owned the others. Later THM-4266/4267 close a larger
region, so `(10)--(13)` are retained component controls rather than a claim
that all 152 are new on the current proof surface.

## 6. Reproduction and scope

```bash
g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Werror -pthread \
  04-computation/lrc14_four_primitive_rays_common_deck_bridge_thm4270.cpp \
  -o /tmp/thm4270-primary
/tmp/thm4270-primary

g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Werror -pthread \
  04-computation/lrc14_four_primitive_rays_literal_wall_audit_thm4270.cpp \
  -o /tmp/thm4270-literal
/tmp/thm4270-literal

python3 -B 04-computation/lrc14_four_primitive_rays_residual_postprocess_thm4270.py
python3 -B -O 04-computation/lrc14_four_primitive_rays_residual_postprocess_thm4270.py
```

This theorem proves only the fixed pool, the four displayed primitive ratios,
and their strict-above-pool scales. It proves no nearby ratio, below-boundary
scale, physical entry of an arbitrary row into this pool, or LRC(14).
**QED.**
