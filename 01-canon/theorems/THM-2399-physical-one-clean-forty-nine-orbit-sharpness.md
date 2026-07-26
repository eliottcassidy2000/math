---
id: THM-2399
title: "Physical one-clean forty-nine-orbit sharpness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. There is
  a primitive, centered, last-lane-typed (1,2,5) scalar packet and a
  strict rational forty-nine-point orbit on which the local scalar
  cover holds at every point while the THM-2396 clean set has exactly
  one root. The common core is literal, q_* is the unique septimal
  depth-one ordinary label, and C_3,c_3 have septimal depth two.
  Hence THM-2396's pointwise multiplicity conclusion N_S(z)>=1 cannot
  be sharpened to N_S(z)>=2 from the same one-orbit incidence and
  valuation data. This does not attain the integrated 66/4459 floor:
  the displayed scalar row has the strict global safe point 3/14 and
  is not an LRC counterexample.
source: codex-2026-07-26-physical-one-clean-orbit
depends_on:
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
related:
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
  - THM-2394-common-core-six-address-transversal-and-labelled-hole-trichotomy
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
script: 04-computation/lrc14_physical_one_clean_orbit_thm2399.py
output: 05-knowledge/results/lrc14_physical_one_clean_orbit_thm2399.out
script_sha256: 2599e39981f7efe9b045f6428d64bef22020d108e7d0729c448b8ccb402666ea
output_sha256: 227add2e7a5eec7f71932a6f13574f0962e0eee055eb44abaee6f78dfb0a94bc
hash_basis: working-tree bytes (LF)
---

# THM-2399 -- physical one-clean forty-nine-orbit sharpness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2396 proves that every high-safe common-core orbit has at least one
clean root. Its search deliberately relaxed the non-core translations,
so it left open whether the physical centered-comb geometry forces two
or more clean roots.

It does not. The obstruction below is physical, rational, strict, and
has all the relevant last-lane valuation labels. It realizes the
one-clean boundary exactly:

```text
one physical 49-orbit
  + the local scalar cover at every orbit point
  + all THM-2396 common-core/high-layer typing
  -> exactly one clean root.                              (1)
```

The word *local* is essential. The same scalar row has a global safe
point, so this is a sharpness theorem for the orbitwise incidence
mechanism, not a counterexample to LRC(14).

## 1. A correctly valued centered packet

As usual, put

```text
D_v={x:||vx||<1/14},

E_H={x:||Hx||<1/7},

C_H=E_H^c.                                             (2)
```

Take the strict rational base

```text
z=8104501/10000019                                     (3)
```

and the scalar data

```text
H=1273,

q_*=77,

(q_1,q_2,q_3,q_4)=(431,643,566,720),

(c_1,c_2,c_3)=(13,169,18193357),

(C_1,C_2,C_3)=(1,13,1399489),       c_i=13C_i.         (4)
```

The exact factorizations are

```text
c_3=49*13^5,

C_3=49*13^4.                                          (5)
```

Consequently the blocker profile is `(1,2,5)`. Moreover:

```text
nu_7(q_*)=1;

nu_7(H)=nu_7(q_i)=0                 for i=1,...,4;

nu_7(C_3)=nu_7(c_3)=2;

(C_1,C_2,c_1,c_2)=(1,13,13,169).                      (6)
```

All six guard/ordinary labels are thirteen-units, the five `q` labels
are distinct, all nine scalar speeds are distinct, and their gcd is
one. Thus (4) has exactly the strict `(1,2,5)` and common-core typing
used in THM-2396.

Writing

```text
V=C_3/49=13^4,
```

direct exact evaluation gives

```text
z notin D_V union D_(13V).                            (7)
```

Hence `C_3` and `c_3` are safe on the whole orbit

```text
O_z={(z+j)/49:j in Z/49Z}.                            (8)
```

## 2. The physical word table

Normalize the orbit index by

```text
n=j+4 mod 49,                    n=r+7s,

r,s in F_7.                                           (9)
```

This sends the physical `D_(C_1)=D_1` word to the horizontal graph

```text
A=(0,0,0,0,0,0,0).
```

Every word below is evaluated from the centered inequalities (2) at
the same base (3); no translate is chosen independently. The other two
common-core words are

```text
B=D_(C_2)=(1,3,5,0,2,3,5),

C=D_(c_2)=(0,4,1,5,1,5,2).                          (10)
```

The guard has the same two addresses in every row,

```text
G=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)),       (11)
```

and the four lower ordinary words are

```text
(0,4,6,3,5,2,4),
(1,2,3,4,4,5,6),
(1,5,2,6,3,6,3),
(1,6,4,2,6,4,2).                                    (12)
```

Finally, `D_(q_*)` is exactly normalized residue row `r=0`.

On each safe row `r=1,...,6`, the two guard addresses and four lower
ordinary addresses are distinct. Their respective holes are

```text
(3,5,5,2,3,5).                                      (13)
```

The corresponding `B` addresses are

```text
(3,5,0,2,3,5).                                      (14)
```

Thus five safe rows have their hole at the actual `c_1=C_2` word.
The sole exception is

```text
(r,s)=(3,5),
```

where the hole is the `C=c_2` address and is outside both quotient
core words `A` and `B`.

## 3. Local cover and the unique clean root

Let

```text
K=1_(E_H)+sum_(i=1)^4 1_(D_(q_i))
```

and use THM-2396's clean set

```text
S={K=0}
   intersection D_(q_*)^c
   intersection D_(c_3)^c
   intersection
     (D_(C_1) union D_(C_2) union D_(C_3))^c.         (15)
```

The word table proves the local scalar cover

```text
O_z intersection C_H
 subset
O_z intersection
 (union_(i=*,1,...,4)D_(q_i)
    union D_(c_1) union D_(c_2) union D_(c_3)).       (16)
```

Indeed, row zero is covered by `q_*`. On every other row, an occupied
address is either guard-dangerous or lower-`q`-dangerous. At the five
holes in (14), `c_1` is dangerous. At the exceptional hole `(3,5)`,
`c_2` is dangerous.

The same inspection gives

```text
#(O_z intersection S)=1.                             (17)
```

In raw indexing the unique point is

```text
j=34,

x=(z+34)/49=348105147/490000931.                     (18)
```

It is the normalized address `(3,5)`. There `K=0`, `q_*` and `c_3`
are safe, and all three quotient blockers are safe, while the actual
`c_2` factor supplies the local scalar cover.

## 4. Strictness and the global hostile control

Across every factor and all forty-nine points, the minimum exact
distance from a danger endpoint is

```text
459267/980001862>0.                                  (19)
```

Therefore the entire word table and the one-clean conclusion persist
on an open interval of base phases. This is not an endpoint
coincidence.

Conversely, the scalar row in (4) is not a global cover. At

```text
x=3/14
```

the guard and every one of the eight ordinary/blocker factors are
strictly safe. Hence the theorem does not construct an LRC
counterexample and does not challenge the conjecture.

## 5. Exact consequence and boundary

For this physical packet,

```text
N_S(z)
 =sum_(j in Z/49Z)1_S((z+j)/49)
 =1.                                                 (20)
```

Therefore THM-2396's pointwise result

```text
N_S(z)>=1
```

is sharp among last-lane-typed centered packets satisfying the scalar
cover on the tested orbit. In particular, no argument using only:

```text
the one-orbit word incidence,
the common-core/high-layer valuations,
and the pointwise scalar-cover clauses on that orbit
```

can replace the right side by two.

This does **not** prove that the integrated bound

```text
delta>=66/4459
```

is globally attained. A stronger global clean-mass theorem could still
use coupling between different base orbits, the fact that the scalar
cover holds on the whole circle, or terminal owner/endpoint phase.
Those are precisely the coordinates absent from this hostile.

No scalar row is excluded, the ledger remains `165`, and LRC(14)
remains open.

## 6. Exact companion

The dependency-free exact companion:

- checks every valuation, distinctness, primitiveness, and common-core
  identity in (4)--(7);
- evaluates all physical masks at the common rational base, reproducing
  (10)--(14) without independent translates;
- verifies the local scalar cover at all forty-nine points;
- finds exactly the one clean point (18);
- computes the strict endpoint margin (19); and
- verifies the strict global safe point `3/14`.

Run

```bash
python3 04-computation/lrc14_physical_one_clean_orbit_thm2399.py
python3 -O 04-computation/lrc14_physical_one_clean_orbit_thm2399.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_physical_one_clean_orbit_thm2399.out
```

after LF normalization. Every assertion remains active under optimized
Python.

## 7. Independent hostile audit

An independent evaluator reconstructed the packet without reusing the
companion's mask tables. It verified `c=13C`, the complete `13`- and
`7`-adic typing, distinctness and primitiveness of the nine scalar
speeds, all seven words, the empty list of local-cover failures, and
the unique clean root. It independently recovered

```text
raw root j=34,
normalized root (3,5),
x=348105147/490000931,
strict margin=459267/980001862,
global safe point=3/14.
```

It also reran the companion under ordinary and optimized Python. Both
LF-normalized transcripts byte-match the stored output, with hashes

```text
source  2599e39981f7efe9b045f6428d64bef22020d108e7d0729c448b8ccb402666ea
output  227add2e7a5eec7f71932a6f13574f0962e0eee055eb44abaee6f78dfb0a94bc.
```

The audit specifically accepted the limiting scope: this is sharpness
of the one-orbit pointwise multiplicity theorem, not attainment of the
integrated clean-mass floor, a global scalar cover, or an LRC
counterexample.
