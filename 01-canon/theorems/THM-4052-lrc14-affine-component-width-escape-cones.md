---
id: THM-4052
title: "LRC(14) affine-component width escape cones"
status: >
  PROVED RELATIVE TO LRCUpTo13 AND THM-4041/4032/4030 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On the physical affine equality boundaries, a
  connected closed pack-safe margin arc cannot fit inside a no-longer open
  spoiled component. For d=2 this gives the exact closure
  1/(42M)>=min(2/(7beta),(alpha+beta-7g)/(7alpha beta)); equivalently every
  hypothetical failure lies in the strict wedge beta<12M and
  alpha beta<6M(alpha+beta-7g). Coarser all-height cones close d=2 when
  E>=12M, d=3 when E>=11M, and d=4 when 3E>=44M. A fully typed rank-eleven
  d=2 row outside THM-4049's residue firewall lies in the new cone. LRC(14)
  remains OPEN.
source: long-precise-frontiers / 2026-08-24
audit: >
  PASS. The primary path checks the THM-4041 component formula and 126,848
  exact scalar-threshold gates, rebuilds three endpoint controls and one
  method hostile, and verifies the finite-box, decoder-tree, crossing,
  residue, width, and direct-clearance gates of a physical rank-eleven row.
  A no-import audit reconstructs strict danger walls and midcell masks on the
  circle, reproduces the complete d=2 component multiset on 780 odd pairs,
  and checks the d/(7E) component bound on all 560 inherited d=3 profiles and
  270 inherited d=4 profiles. The semantic digests agree; normal and
  optimized streams byte-match both frozen outputs.
depends_on:
  - THM-4041-lrc14-d2-affine-defect-edge-boundary
  - THM-4032-lrc14-d3-affine-defect-lattice-boundary
  - THM-4030-lrc14-d4-affine-defect-lattice-boundary
  - LRCUpTo13
related:
  - THM-4049-lrc14-d2-two-phase-residue-firewall
  - THM-4024-lrc14-complete-divisor-incidence-envelope
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
script: 04-computation/lrc14_affine_component_width_escape_thm4052.py
output: 05-knowledge/results/lrc14_affine_component_width_escape_thm4052.out
script_sha256: d0165375b4d3d510549a5ef6d0584c3fd080bd7039bb4ccd0005bd6097628efc
output_sha256: 120b76f5d3fe5a6ae107e4c5579ec07c52d22a097f496e64cc0742bfcc999d91
independent_audit_script: 04-computation/lrc14_affine_component_width_escape_thm4052_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_affine_component_width_escape_thm4052_independent_audit.out
independent_audit_script_sha256: b30309b3885f99786fc9c264b31e458ec19a65b73c8a9af0568c775d65e3801b
independent_audit_output_sha256: 458084ed8923d285d5e3045a7f9655ce38b8eb30c3ff386612763bc26adf4ee4
semantic_sha256: 0d308969f4210b175ee5999bff8fb6d8045da92730b9e7007791ded447442e8b
hash_basis: raw LF bytes
---

# THM-4052 -- affine-component width escape cones

**PROVED RELATIVE TO LRCUpTo13 AND THM-4041/4032/4030 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.** This theorem closes unbounded cones inside all three
remaining affine equality moduli. It does not close the strict complementary
cones and does not prove LRC(14).

Write `||x||` for distance to the nearest integer.

## 1. The exact d=2 theorem

Let `H` be an eleven-element set of distinct positive integers, put
`M=max H`, and let `alpha<beta` be distinct positive odd integers. Write

```text
g=gcd(alpha,beta),       alpha=ga, beta=gb,       gcd(a,b)=1.   (1)
```

For the thirteen-speed row

```text
V=2H union {alpha,beta},                              (2)
```

the following hold.

1. If `a+b<=7`, then `V` is `1/14`-lonely.
2. If `a+b>7`, put

   ```text
   W=min(2/(7 beta),(alpha+beta-7g)/(7 alpha beta)). (3)
   ```

   Then `V` is `1/14`-lonely whenever

   ```text
   1/(42M) >= W.                                     (4)
   ```

Equivalently, every hypothetical failure in this lane must satisfy all three
strict conditions

```text
a+b>7,
beta<12M,
alpha beta<6M(alpha+beta-7g).                        (5)
```

In particular, the simpler cone `beta>=12M` closes. Since `beta` is the
largest exception, this is `E>=12M` in the unified notation below.

## 2. The closed-arc/open-component proof

The inherited eleven-speed result gives a pack phase `y_0` with

```text
||h y_0||>=1/12             for every h in H.         (6)
```

Distance to the nearest integer is one-Lipschitz, so each function
`y -> ||hy||` is `h`-Lipschitz. Therefore the connected closed circle arc

```text
I={y:dist(y,y_0)<=1/(84M)}                            (7)
```

lies in the `1/14`-safe set `G(H)` and has length `1/(42M)`.

For any `y in G(H)`, both lifts

```text
x_j=(y+j)/2,                  j=0,1,                 (8)
```

preserve the pack inequalities, because

```text
||2h x_j||=||h(y+j)||=||hy||.                        (9)
```

Let `Sigma_(alpha,beta)` be THM-4041's open set of pack phases for which the
two odd exceptions spoil both labels. If `y` is outside `Sigma`, one lift in
`(8)` is safe for the whole row. Thus a non-lonely row would force

```text
I subset G(H) subset Sigma_(alpha,beta).              (10)
```

THM-4041 gives every positive-length component of `Sigma` the length

```text
(2/g)L_r,
L_r=min(1/(7b),(a+b-7r)/(14ab)),                     (11)
```

for a positive odd `r` satisfying `7r<a+b`. These lengths decrease with
`r`, so their maximum is exactly `(3)`. Every component is an **open** arc.
A connected closed arc of length at least `W` cannot be contained in an open
arc of length at most `W`; equality is included. Hence `(4)` contradicts
`(10)` and proves the theorem. If `a+b<=7`, THM-4041 says `Sigma` is empty.

Finally, `(4)` says that `1/(42M)` dominates at least one of the two entries
in `(3)`. Clearing positive denominators gives

```text
beta>=12M       or       alpha beta>=6M(alpha+beta-7g), (12)
```

which is exactly the complement of the strict wedge `(5)`.

## 3. Unified d=3 and d=4 corollaries

At the `d=3` and `d=4` boundaries, let `H` be the inherited ten-speed pack,
`M=max H`, and `E` the largest exception. The ten-speed citation supplies
clearance `1/11`, hence the same Lipschitz argument gives a closed pack-safe
arc of length

```text
2(1/11-1/14)/M=3/(77M).                              (13)
```

For `d=3`, the three exceptions are units modulo three. Each can spoil at
most one of the three lifts, so full spoilage has the unique capacity
partition `1+1+1`. For `d=4`, the exceptions have the type

```text
2r,a,b,                  r,a,b odd,                  (14)
```

and full spoilage has the unique partition `2+1+1`. Consequently distinct
label-assignment strata are disjoint open sets. Every connected component of
the fully spoiled set stays inside one assignment and, in particular, inside
one danger tooth of the largest exception. In the pack-phase coordinate that
tooth has length at most

```text
d/(7E).                                               (15)
```

Comparing `(13)` and `(15)`, with the same closed-versus-open equality
argument, proves

```text
d=3:  E>=11M       implies escape,
d=4:  3E>=44M      implies escape.                   (16)
```

The literal independent circle audit attains width ratio one in both bounded
profile universes, so the component bound itself cannot be uniformly
improved without retaining another coordinate.

## 4. Exact endpoint controls and fixed-bank hostiles

The following central cited phase is fully spoiled, yet the margin arc
escapes at the displayed labelled lift:

```text
d   H             exceptions       central y   escape y     lift x       clearance
2   {1,...,11}    (1,133)          1/12        13/154       167/308      1/14
3   {1,...,10}    (1,110,23)       1/11        13/140       51/140       1/14
4   {1,...,10}    (2,185,11)       1/11        137/1540     4757/6160   137/1540
```

For `d=2`, `H={1,...,11}` and exceptions `(1,11)` spoil the central phase
and both guaranteed margin endpoints. This is a hostile to claiming escape
below `(5)`, not an LRC counterexample.

Nor can the two displayed `d=3,4` endpoint phases be frozen into an
all-height bank. The first typed exact hostiles are

```text
d=3: exceptions (1,55,56), phases {1/14,13/140},
d=4: exceptions (26,31,57), phases {15/98,1/14}.      (17)
```

They spoil every lift at both named phases. Their full rows are nevertheless
positive controls, with clearance `1/13` at `x=1/13` and `1/11` at `x=1/11`,
respectively. Thus adaptive component placement, not another fixed two-time
bank, is the live `d=3,4` operation.

## 5. A physical d=2 row outside the residue firewall

This cone is not subsumed by THM-4049. Let

```text
S=(37,43,61,67,73,79,97,103,127),
P=3*5*product_(r in S) r=713721382004055345,
H={P/r:r in S} union {1,4},
alpha=P/5,                    beta=P/3.               (18)
```

Then `M=P/37` and `beta>12M`, so `(12)` closes
`2H union {alpha,beta}`. But

```text
H mod 56=(37,11,5,43,41,23,17,47,39,1,4),           (19)
```

which meets THM-4049's forbidden set in `{11,23}`.

The row is genuinely in the inherited rank-eleven physical lane. Take

```text
s=1, t=2, v=(1,4),
u=(2P/r)_(r in S) direct-sum (P/3,P/5).               (20)
```

It has thirteen distinct speeds, primitive body, `c_2=9`, and

```text
sum n_i=574570283588268864 < (91^6)^2.               (21)
```

An explicit decoder spanning tree uses the edges

```text
37-73, 73-43, 73-97, 43-67, 67-103,
37-79, 37-127, 103-61, 79-(P/3), 37-(P/5),           (22)
```

with admissible primitive sums

```text
110,116,170,110,170,116,164,164,85,47,               (23)
```

while the pair component has sum `5`. Internal relation height is at most
`127`. The exact lower bounds for a crossing relation are

```text
two body + one pair: 54561683510745,
one body + two pair: 1123970680321347,                (24)
```

both greater than `91^6=567869252041`. Thus the decoder graph has components
`11+2`, rank eleven, and no bounded crossing row; the THM-3818 physical typing
is literal rather than heuristic. Independently, `x=1/14` gives full-row
clearance exactly `1/14`.

## 6. Exact scope

THM-4052 closes three unbounded height cones and gives the sharp d=2 scalar
wedge. It does not classify the remaining physical projection, close the
strict complementary cones, or improve the cited lower-dimensional LRC
input. The next exact objects are the intersections of the connected
pack-safe components with the affine spoiled components inside

```text
d=2: beta<12M and alpha beta<6M(alpha+beta-7g),
d=3: E<11M,
d=4: 3E<44M.                                         (25)
```

Those are still infinite physical families. LRC(14) remains open.

**QED.**
