---
id: THM-4270
title: "Uniform four-primitive-outsider-ray common-deck closure"
status: >
  PROVED RELATIVE TO THM-4231/4256 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  For the fixed thirty-label pool and every labelled nine-body, the rays
  3:5, 7:8, 8:9, and 11:12 are Haar-safe at every strict-above-pool scale.
  One frozen 8,192-repair carrier supplies four ratio-specific common decks;
  each deck covers all 14,307,150 bodies throughout its finite bridge, while
  THM-4231 supplies the cofinite tail. After exact overlap accounting against
  THM-4266/4267/4269, the current contribution is 146 edges and the residual
  has 177,323 edges with unique top edge (520,688). The earlier
  post-THM-4261/4262 component had 152 edges; exactly six are already owned by
  THM-4266. No other ratio, smaller scale, physical entry, or LRC(14) follows.
source: codex-creative-frontiers-20260827
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4256-uniform-two-three-outsider-ray-endpoint-cocycle-closure
  - THM-4261-semantic-endpoint-band-prefix-union-lift
  - THM-4262-uniform-three-four-outsider-ray-common-deck-closure
  - THM-4266-three-round-learned-carrier-endpoint-descent
  - THM-4267-uniform-four-five-outsider-ray-common-deck-closure
  - THM-4269-uniform-five-six-outsider-ray-common-deck-closure
common_header: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_common.hpp
primary_script: 04-computation/lrc14_four_primitive_rays_common_deck_bridge_thm4270.cpp
primary_output: 05-knowledge/results/lrc14_four_primitive_rays_common_deck_bridge_thm4270.out
primary_script_sha256: 9910c832b9c7a1f3ae1cb089799a9b625bbb72ec00677bafb00368be6c35a48b
primary_output_sha256: d50608147c4214fe9fb5f1b336f165b10ca5e96680ace4c08a758fe2421d87c4
independent_script: 04-computation/lrc14_four_primitive_rays_literal_wall_audit_thm4270.cpp
independent_output: 05-knowledge/results/lrc14_four_primitive_rays_literal_wall_audit_thm4270.out
independent_script_sha256: de94d6dc80ec539fa6ec751e11b3bca5df307d6bbaa4085c3b066a05ae4fcf26
independent_output_sha256: 2514e7d7a0fa8259ce3a1a181cd9455f857fc5a393fcae57a34f6182a992ad36
postprocess_script: 04-computation/lrc14_four_primitive_rays_residual_postprocess_thm4270.py
postprocess_output: 05-knowledge/results/lrc14_four_primitive_rays_residual_postprocess_thm4270.out
postprocess_script_sha256: b8cd34922a769478a29cfacb16fefc56bdd10451e0c5946c64f2d3107b31ba93
postprocess_output_sha256: adb20bedd251a05c2536b48460a6ccdc2266ec7c14fd38c547c6b00085adb49c
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

## 5. Exact current proof-graph contribution

The post-THM-4266 residual is

```text
count=177585,
FNV=6ce05d05eb01daed,
SHA256=009614651bb81e9763b2a9ff4b580497
       bfb6978a6c69d18cf986346e369374d9.                (10)
```

The postprocessor reconstructs THM-4267 semantically from `(10)`, checks its
63-edge `4:5` deletion, and recovers

```text
count=177522,
FNV=33142f955cc93379,
SHA256=d277aebe296153ead14a77207ea1499c
       961c8b06796b7e62f324e34f7a9ef087.                (11)
```

It then reconstructs THM-4269's disjoint 53-edge `5:6` component and recovers

```text
count=177469,
FNV=4d1feae0c1e653d5,
SHA256=289cede32347b364123827e7dea02d728
       b71e8c87d079a9892d3e0492b4a08ae.                (12)
```

Semantic selection of the four proved bridge rays in `(12)` finds

```text
ratio       3:5   7:8   8:9   11:12
new edges    32    44    40      30.                   (13)
```

The four lists are pairwise disjoint. Their 146-edge union has

```text
FNV=3115735824bdb7f5,
SHA256=b2fcda2e0e602a4284d243d486ae426f
       7376de0be429c30d1bb40ee6c455f750.                (14)
```

In the earlier post-THM-4261/4262 order the component had 152 edges. The six
now-owned overlaps are exactly

```text
(602,688), (616,704),
(616,693), (632,711), (640,720),
(638,696),                                               (15)
```

all closed by THM-4266; THM-4267's `4:5` component is ratio-disjoint.
THM-4269's `5:6` component is likewise ratio-disjoint.
Deleting exactly the current 146-edge union gives

```text
count=177323,
FNV=f1dcc8033fa727d9,
SHA256=8c0b1fac01d00bd54784178034f4e5f2
       1c2a29ea95a9cb0ed5a63b06fbc20872,
maximum endpoint=688, unique top edge=(520,688).        (16)
```

The current postprocessor asserts both inherited ledgers, every per-ray
count/FNV/SHA, the pairwise-disjoint union, and the updated ledger. Thus
`(13)--(16)`, rather than the older 152-edge component, are the exact current
proof-graph consequence.

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
