---
id: THM-2246
title: "Depth-one private joint two-step fibre cap"
status: >
  PROVED + VERIFIED-EXACT + LOCALLY SHARP. In THM-2234's first-depth-one
  private branch, the guard occupancy cap 10 and peeled-complement occupancy
  cap 12 cannot saturate independently: every two-step fibre contains at
  most 112 private roots, not 120. Consequently the second private image has
  measure at least 33709/776160, improving 33709/831600 by the exact factor
  15/14. An exact 169-root phase with five distinct aligned unit masks and a
  labelled depth-two owner attains 112, so no further uniform gain follows
  from the guard, peeled-complement, five unit labels, and one active owner
  alone. The new floor remains below 1/7 and excludes no valuation profile;
  it does not prove LRC(14).
source: codex-2026-07-25-depth-one-joint-private-fibre
depends_on:
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
related:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2239-unrestricted-multicore-signed-dual-profile-exclusion
script: 04-computation/lrc14_depth_one_private_joint_two_step_cap_thm2246.py
output: 05-knowledge/results/lrc14_depth_one_private_joint_two_step_cap_thm2246.out
script_sha256: 938e175990f22feb72b5dbd00ad147279593358b705d0b99af7a157944def791
output_sha256: aac91038011114312c11aeb41be67fa436a0c35bd43a65ae2413230375281a46
hash_basis: working-tree bytes (LF)
---

# THM-2246 -- the private two-step cap is 112, not 120

Retain THM-2234's notation. Thus

```text
T(x)=13x mod 1,
D_b={x:||bx||<1/14},
C_H={x:||Hx||>1/7},
c_1=13a,                         13 does not divide a,

P=C_H minus (union_(i=1)^5 D_(q_i) union D_(c_1)),    (1)
```

and the scalar cover makes `P` a two-owner private remainder:

```text
P subset D_(c_2) union D_(c_3).                       (2)
```

THM-2234 proves

```text
measure(P)>=2593/90090,                               (3)
T(P) subset D_a^c.                                    (4)
```

It obtained its second-image floor by multiplying the separate root caps
`10` and `12`. The present theorem retains the intermediate guard bit and
shows that the exact joint cap is smaller.

## 1. A two-level root invoice

Fix a generic `z` and consider the thirteen intermediate roots

```text
R_z={y:T(y)=z}.                                       (5)
```

Put

```text
d=1_(D_a)(z),             h=1_(C_H)(z),

A=#{y in R_z:y notin D_a},
G=#{y in R_z:y in C_H},
I=#{y in R_z:y notin D_a and y in C_H}.              (6)
```

The exact danger and guard root counts from THM-2234 give

```text
A=11+d,                    G=10-h.                    (7)
```

Since `R_z` has thirteen points, the finite-set inclusion-exclusion bound is

```text
I>=A+G-13.                                              (8)
```

For an allowed intermediate root `y notin D_a`, at most

```text
10-1_(C_H)(y)                                          (9)
```

of its thirteen roots can lie in `P`, because `P subset C_H`. Therefore the
two-step occupancy

```text
n_P^(2)(z)=#{x in P:T^2(x)=z}                         (10)
```

satisfies

```text
n_P^(2)(z)
 <=sum_(y in R_z, y notin D_a)(10-1_(C_H)(y))
 =10A-I
 <=10A-(A+G-13).                                      (11)
```

The four future-bit cases are exact:

```text
 d   h    A    G    I floor    cap in (11)
 0   0   11   10       8          102
 0   1   11    9       7          103
 1   0   12   10       9          111
 1   1   12    9       8          112.               (12)
```

Hence, almost everywhere,

```text
n_P^(2)(z)<=112.                                      (13)
```

This is a co-adapted constraint: the peeled-complement count and guard count
live on the same intermediate root fibre. Treating their separate maxima as
independent loses the intersection debit `I`.

## 2. The improved second-image floor

Haar disintegration over the `169` roots of `T^2` gives

```text
169 measure(P)=integral n_P^(2)(z) dz
              <=112 measure(T^2(P)).                 (14)
```

Combining (3) and (14),

```text
measure(T^2(P))
 >=(169/112)(2593/90090)
 =33709/776160
 =0.0434304782519068....                              (15)
```

THM-2234's sequential value and the exact gain are

```text
33709/831600,

33709/776160-33709/831600
 =33709/11642400.                                     (16)
```

Equivalently, replacing the product cap `10*12=120` by `112` improves the
floor by the factor

```text
120/112=15/14.                                        (17)
```

The deep-owner labels in (2) transport as before. In particular, if both
remaining blockers have depth at least two,

```text
T^2(P)
 subset D_(c_2/13^2) union D_(c_3/13^2)              (18)
```

almost everywhere. The gain in (15) does not discard that two-owner carrier.

## 3. Sharp hostile phase, including all five unit labels

The cap `112` cannot be reduced from the information used above, even if the
five distinct unit masks and one labelled depth-two owner are retained.
Take

```text
H=1,              a=2,
z=325007/700000,

q_i=1+169*700000*i,                 1<=i<=5,

c_1=13a=26,       c_2=169*2=338.                    (19)
```

All five `q_i` are distinct positive thirteen-units. For every two-step root

```text
x=(z+k)/169,                       0<=k<169,          (20)
```

one has

```text
q_i x = x mod 1,                                      (21)
```

because `700000*i*z=325007*i` is integral. Thus all five unit masks align
with `D_1` on this fibre. Since `D_1` is disjoint from `C_1`, they remove no
guard root.

Moreover `z in D_2`, so the labelled blocker `c_2=338` owns every root in
(20). Exact root counting gives

```text
A=12,       G=9,       I=8,

#{x in (20):
    x in C_1,
    T(x) notin D_2,
    x notin union_i D_(q_i)}
 =112.                                                (22)
```

The displayed rational phase lies off every load-bearing guard, peeled, and
unit boundary. Hence the equality persists on a nonempty phase interval;
it is an essential-supremum obstruction, not a null endpoint artifact.

This hostile control need not satisfy the global scalar cover. Its role is
the precise one permitted by the validity gate: it proves that no smaller
two-step fibre cap follows from the local guard, peeled-complement, five
unit-mask, and active-owner data used in (11).

## 4. Exact stopping boundary

The improved floor is still strictly below the mass of one danger comb:

```text
1/7-33709/776160
 =77171/776160>0.                                     (23)
```

Consequently (15) does not force positive mass to remain after either one
of the two owners reaches unit depth and is deleted. It excludes no
first-depth-one valuation profile.

The first failed implication is now exact. A continuation needs at least one
of:

```text
an owner-private floor exceeding the one-comb deletion cost;
a cross-fibre winding restriction excluding the aligned hostile phase;
or a labelled incidence law coupling the two deep owners across fibres.  (24)
```

Another multiplication of the separate local caps cannot supply such a
gain. QED.

## 5. Reproduction

Run

```bash
python3 04-computation/lrc14_depth_one_private_joint_two_step_cap_thm2246.py
python3 -O 04-computation/lrc14_depth_one_private_joint_two_step_cap_thm2246.py
```

Both modes reproduce the stored transcript byte for byte. The companion
checks the four-case cap table, all rational identities, all `169` hostile
roots, exact mask alignment, the active depth-two owner, and boundary
avoidance.
