---
id: THM-3313
title: "Projected-k3 z216 third ruler-cost prefix multicover closure"
status: >
  PROVED + VERIFIED-EXACT in the declared projected atlas.  After the two
  independently audited `z1=216` intrinsic-ruler cost prefixes, the next
  deterministic prefix consists of the `gcd72/L72072` singleton row 157 and
  the complete `gcd24/L25872` family at rows 53, 293 and 357.  All four exact
  necessary upper screens are empty: `423=289` crude `+134` status `+0`
  residual.  THM-3308's deterministic compiler presents the status kills as
  133 one-layer modular covers and one support-minimal two-layer circuit.
  Therefore the projected ledger moves `373161 -> 373157`, the `z216` wall
  moves `357 -> 353`, complete live families move `33 -> 31`, and the cap
  remains 216.  Physical entry, arbitrary `k<=1`, the rung and LRC(14) remain
  open.
audit: >
  The exact companion pins the independent two-prefix audit and THM-3308's
  multicover companion, reconstructs the 480-row atlas, removes every older
  closure, reranks all 357 live rows in 33 complete `(gcd(216,L),L)` families,
  and verifies a strict cost boundary before selecting the four rows.  It
  checks every inherited Farkas vector over exact arithmetic but persists only
  deterministic modular covers, reconstructs all three feasible singleton
  tail tables for the unique two-layer circuit, and replays both the row-64
  survivor and the inherited genuine three-layer control.  Normal and
  optimized outputs byte-match the frozen transcript; the source has zero
  assertion nodes and zero floating literals.
source: root/creative-synthesis-recover/2026-08-03
depends_on:
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3281-projected-k3-z216-three-natural-wall-family-screen-descent
  - THM-3308-threshold-chain-modular-multicovers-and-three-layer-status-circuit
script: 04-computation/lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.py
output: 05-knowledge/results/lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.out
script_sha256: 34710b907f3274057d6cb60a0fd7bd72ae2c3879a90b830a3e376d30e5198646
output_sha256: c664acf4e02b40835dfdc76ff135d17332e3d0e57b489e61b19c0f2f7b5998f5
semantic_sha256: 4d8547f53b287d82bb725cd6b3ed78c390ba03a37a3dcb3ab76ac98389e045a3
hash_basis: LF-normalized bytes
---

# THM-3313 -- projected-k3 z216 third ruler-cost prefix multicover closure

**PROVED + VERIFIED-EXACT in the declared projected atlas.**

This theorem applies the support-minimal threshold compiler of THM-3308 to a
new complete-family work prefix and makes direct progress on the live
projected `k=3,z1=216` ledger.

## 1. Deterministic queue

Start from the independently audited state after the first two intrinsic-cost
prefixes: `357` live wall rows in `33` complete ruler families.  Reconstruct
all prior closures and rank each live family by

```text
cost(g,L)=sum_(rows E in the family)L(E)r(E),                (1)
```

where `r(E)` is the exact component count.  Stop after the next nonsingleton
family.  The queue begins

| rank | family | rows | atlas indices | cost |
|---:|---|---:|---|---:|
| 1 | `gcd72/L72072` | 1 | `157` | 1,729,728 |
| 2 | `gcd24/L25872` | 3 | `53,293,357` | 2,121,504 |
| 3 | `gcd24/L76440` | 1 | `141` | 2,293,200 |

The selected prefix contains ranks one and two.  The cost boundary before
rank three is strict.  All three rows of the `L=25872` fibre are retained.

Their labelled bodies and component counts are

```text
157: {1,4,9,11,12,13},       r=24;
 53: {1,2,6,8,11,14},        r=26;
293: {2,4,6,8,11,14},        r=26;
357: {2,6,8,11,12,14},       r=30.                           (2)
```

Intrinsic cost is a deterministic work order, not a safety invariant.

## 2. Exact empty screens

The exact ray/common-status screen gives

| family | states | crude kills | exact status kills | residual |
|---|---:|---:|---:|---:|
| `gcd72/L72072` | 111 | 93 | 18 | 0 |
| `gcd24/L25872` | 312 | 196 | 116 | 0 |
| **total** | **423** | **289** | **134** | **0** |

Rowwise,

```text
157: 111=93+18+0,
 53:  27=17+10+0,
293:  14= 9+ 5+0,
357: 271=170+101+0.                                        (3)
```

Every status dual is checked over exact arithmetic.  No residual passport
remains for a later terminal probe.

The valid logical direction is

```text
empty exact necessary upper screen
  => the selected projected wall row is empty.              (4)
```

Screen feasibility is not used as a converse.

## 3. Modular multicover presentation

THM-3308 recertifies all 134 status kills through pointwise modular majorants
of sums of nested capacity tails.  The new taxonomy is

```text
133 one-layer certificates,
  1 two-layer certificate,
  0 uncompiled status kills.                                (5)
```

The unique multi-layer address is

```text
row 157, divisors (99,616,1001,8008),
q=9009, m=(1365,1287,1287,1287),
histogram=((0,99),(1,1338),(2,2382),(3,1992),(4,2874),(5,324)).             (6)
```

At thresholds `T=(3,4)`, the primitive weight vector `w=(2,1,2,1)` satisfies
the pointwise sixteen-pattern inequality

```text
1[c(P)>=3]+1[c(P)>=4]
 <=2*1[0 in P]+1[1 in P]+2*1[2 in P]+1[3 in P].             (7)
```

The required tails total

```text
H_3+H_4=5190+3198=8388,
```

whereas the weighted marginals total

```text
2*1365+1287+2*1287+1287=7878.                               (8)
```

The contradiction gap is `510`.  Exact feasible common tables exist for
each singleton tail system at thresholds `3`, `4` and `5`; hence support two
is minimum in the full common-table dual cone, not only in the integral
modular template.

The row-64 negative control still retains eight screen residuals, including
an explicit feasible common table.  The inherited row-138 control still
requires three layers.  Thus the compiler neither kills every packet nor
silently caps circuit arity at two.

## 4. Ledger consequence and boundary

The four exact exclusions give

```text
projected ledger:          373161 -> 373157,
z1=216 wall rows:              357 -> 353,
complete live families:          33 -> 31,
projected k=3 cap:              216 unchanged.              (9)
```

The next ranked families are

```text
gcd24/L76440 singleton,           cost 2,293,200;
gcd24/L30576 three-row family,    cost 2,629,536;
gcd72/L7056 nineteen-row family,  cost 3,400,992.             (10)
```

The first two form the next deterministic prefix through a nonsingleton
family.  The nineteen-row fibre is a natural test for reusing capacity
signatures rather than screening each row independently.

Everything here remains inside a necessary projected `k=3,z1=216` atlas.
The theorem restores no physical speed entry, endpoint, owner, phase,
current, arbitrary `k<=1`, rung, or LRC(14) conclusion.
