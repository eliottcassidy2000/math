---
id: THM-2325
title: "Prescribed-target-gain full-lattice 91-unit needle bank"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Under
  THM-2309's primitive scalar and sharp two-coordinate septimal-support
  hypotheses, every nonzero vector in the two-target quotient has at least
  3,134,566,563,840 exact full-lattice relation addresses of controlled
  Bezout height whose nine coordinates are units modulo 91. Consequently
  each of the fourteen projective target directions has at least
  37,614,798,766,080 such addresses. The complete ordered mod-thirteen
  character-pair label of any THM-2321 cubic cell can therefore be copied
  by a genuine full-lattice Fourier needle. No address is proved to
  participate in that current or to belong to the bounded visible carrier;
  no LRC(14) row is excluded.
source: codex-2026-07-25-prescribed-gain-needle-bank
depends_on:
  - THM-2301-essential-affine-arrangement-and-visible-rank-six-address-bank
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
related:
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2318-one-shot-three-prime-mobius-amplifier
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
  - THM-2322-local-hostile-c3-toothpick-ladder
script: 04-computation/lrc14_prescribed_gain_91_unit_needle_bank_thm2325.py
output: 05-knowledge/results/lrc14_prescribed_gain_91_unit_needle_bank_thm2325.out
script_sha256: 7b015639d6f99127fc77202711be79f8887f7c16954fa0eb156cf8e2a04ac7d4
output_sha256: d98d72eae1cb21645c813b60e4e3e95ee735e8f04457246976a3ab960c68648e
hash_basis: working-tree bytes (LF)
---

# THM-2325 -- every target-gain vector has a full-lattice needle bank

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2309 identifies the full mod-thirteen scalar relation kernel modulo
its owner-aligned six-packet as an exact two-target plane. THM-2321
independently identifies the ordered cubic character pair with the same
two-dimensional coordinate plane. The missing statement was not merely
whether all fourteen projective labels exist, but whether a prescribed
**vector** label in that plane contains even one genuine full-lattice
Fourier address.

It contains an enormous bank. The proof has three orthogonal parts:

```text
mod 13:
  count all-unit points in every affine target-vector fibre;

mod 7:
  count all-unit points in the entire scalar kernel exactly;

over Z:
  pair by CRT and correct by one Bezout vector without changing residues.
                                                                  (1)
```

The result closes target-vector **availability**, including the common
scalar on the projective line. It does not couple an address coefficient
to a word current.

## 1. Scalar and owner-packet hypotheses

Use THM-2309's nine scalar labels

```text
E=U disjoint_union {j,a,b},
U={u_0,u_1,...,u_5}.                                (2)
```

Let

```text
w=(w_e)_(e in E) in Z^9
```

be primitive, with every `w_e` nonzero and

```text
13 does not divide w_u,             u in U,
13 divides w_j,w_a,w_b.                              (3)
```

Write

```text
S=sum_e |w_e|,
supp_7(w)={e:7 does not divide w_e},
s=|supp_7(w)|.                                      (4)
```

Assume the sharp septimal condition

```text
s>=2.                                               (5)
```

Choose a Bezout vector `z in Z^9` with

```text
z.w=1,
B(w)=||z||_infinity.                                (6)
```

Fix the selected source `j` and any omitted unit `u_0`. Let `R_13` be
THM-2309's target-grafted owner packet modulo thirteen and put

```text
K_13=(w mod 13)^perp subset F_13^9,
L_13=rowspace(R_13) subset K_13.                    (7)
```

The audited packet properties are:

```text
dim L_13=6,
every one of its nine coordinate functionals is nonzero,

K_13
 =L_13 direct_sum span(e_a,e_b).                    (8)
```

The second line in (8) is complete column brightness. It will be
load-bearing.

## 2. Every target vector has many thirteen-unit addresses

For any nonzero target vector

```text
q=r e_a+s' e_b,             (r,s')!=(0,0),          (9)
```

consider the affine fibre

```text
X_q=q+L_13 subset K_13.                             (10)
```

Parameterize `L_13` by six row coefficients. For coordinate `e in E`,
the condition

```text
x_e=0,                  x in X_q,                   (11)
```

is an affine hyperplane in `F_13^6`. Its normal is the `e`th column of
`R_13`, hence is nonzero by brightness. The nine normals span the full
six-dimensional dual because `R_13` has row rank six.

Apply THM-2301's sharp essential-affine-arrangement bound with

```text
q_field=13,       dimension=6,       hyperplanes=9.
```

It gives, in **every** vector fibre (10),

```text
|X_q intersection (F_13^*)^9|
 >=(13-9+6-1)12^5
 =9*12^5
 =2,239,488.                                         (12)
```

The hyperplanes may coincide; that only increases the complement. Thus no
generic-position assumption is hidden in (12).

## 3. Projectivization retains all twelve scalar fibres

The nonzero quotient vectors split into the fourteen projective
directions

```text
P^1(F_13)
 ={[1:0],[0:1]} union {[1:g]:g in F_13^*}.          (13)
```

For a fixed direction `d=[r:s']`, its twelve vector representatives are

```text
t(r e_a+s' e_b),               t in F_13^*.         (14)
```

Their affine fibres are disjoint: two fibres have the same image in
`K_13/L_13` only when their target vectors agree. Therefore every
projective direction already has at least

```text
12*9*12^5
 =26,873,856                                            (15)
```

all-coordinate-unit mod-thirteen relation vectors.

This scalar factor is structurally important. THM-2321's projective
current sums over exactly the same twelve nonzero character scales. Here
none of those vector scales is absent from the relation lattice.

## 4. The full septimal kernel and the sharp support boundary

At seven, using only a six-packet leaves most of the relation kernel
unused. Put

```text
K_7=(w mod 7)^perp subset F_7^9.                    (16)
```

The number of all-coordinate-unit points in `K_7` depends only on the
support size `s`.

The `9-s` zero-weight coordinates are arbitrary units, contributing
`6^(9-s)`. On the support, rescale `x_e` by the nonzero coefficient
`w_e`. The remaining problem is to count

```text
(y_1,...,y_s) in (F_7^*)^s,
y_1+...+y_s=0.                                      (17)
```

Additive-character orthogonality gives exactly

```text
A_s
 =(1/7)[6^s+6(-1)^s].                               (18)
```

Indeed, the trivial character contributes `6^s`, while for each of the
six nontrivial additive characters the sum over `F_7^*` is `-1`.
Consequently

```text
|K_7 intersection (F_7^*)^9|
 =6^(9-s)[6^s+6(-1)^s]/7.                           (19)
```

For `2<=s<=9`, the minimum occurs at `s=3`:

```text
|K_7 intersection (F_7^*)^9|
 >=5*6^7
 =1,399,680.                                        (20)
```

One quick way to see the minimum is to rewrite (19) as

```text
[6^9+(-1)^s 6^(10-s)]/7.
```

Even `s` lie above the common first term; among odd `s>=3`, the negative
correction decreases in magnitude as `s` grows.

The hypothesis (5) is exact. At `s=1`, (18) is zero; equivalently the
unique supported coordinate of every septimal relation must vanish.
There is then no all-seven-unit address at any height.

## 5. CRT pairing and exact Bezout lifting

Fix one of the mod-thirteen points counted in (12) and one of the
septimal points counted in (20). Coordinatewise CRT gives a unique

```text
y in (Z/91Z)^9                                      (21)
```

with those reductions. Both reductions are scalar relations, hence

```text
y.w=0 mod 91.                                       (22)
```

Every coordinate of `y` is a unit modulo `91`, and its target quotient
modulo thirteen remains the prescribed vector `q`. CRT is a bijection, so
the pair count is injective. Thus every nonzero target-vector fibre has at
least

```text
(9*12^5)(5*6^7)
 =3,134,566,563,840                                 (23)
```

distinct all-`91`-unit modular relation addresses. Every projective
direction has twelve disjoint fibres and therefore at least

```text
12(9*12^5)(5*6^7)
 =37,614,798,766,080.                               (24)
```

Across all `168=13^2-1` nonzero target vectors, the banks contain at
least

```text
168(9*12^5)(5*6^7)
 =526,607,182,725,120                               (25)
```

distinct residues.

Now choose the centered representative of `y`, still denoted `y`, in
`[-45,45]^9`, and define

```text
r_y=y-(y.w)z.                                       (26)
```

Equations (6) and (22) give

```text
r_y.w=0,
r_y=y mod 91,

||r_y||_infinity
 <=45+45 S B(w)
 =45(1+S B(w)).                                     (27)
```

Distinct modular addresses remain distinct exact relations. Hence
(23)--(25) are also counts of controlled-height exact integer relation
addresses.

Because every coordinate is a seven-unit, all nine safe-interval Fourier
factors of length `5/7` or `6/7` are nonzero. These are genuine
full-lattice Fourier needles, not merely formal relations.

## 6. Complete character-pair labels now have lattice witnesses

THM-2321 uses a cubic cell

```text
C_(k,l)=M_k M_l conjugate(M_(k+l)).                 (28)
```

Its complete ordered input-character label maps linearly to the target
quotient:

```text
(k,l)
  -> k e_a+l e_b in K_13/L_13.                     (29)
```

For every allowed nontrivial cell, the vector in (29) is nonzero.
Equation (23) proves that its **exact vector fibre**, not just its
projective gain class, contains a huge bank of full-lattice Fourier
addresses. In particular:

```text
THM-2321 fixed first character:
  prescribe k, select l by row positivity,
  then use the bank in k e_a+l e_b;

THM-2321 fixed projective gain:
  select a positive scalar t on [r:s'],
  then use the bank in t(r e_a+s' e_b);

conditional THM-2318 composition:
  take k to be the collision atom's normalized root colour,
  select l, and retain the same vector label k e_a+l e_b.           (30)
```

The last line also requires the scalar row to satisfy (5).

Thus neither root colour, abstract gain, projective scale, nor existence
of a genuine exact address with the same complete residue label is the
remaining obstruction.

### Why this is still not incidence

The map (29) matches labels. It does not show that one of the addresses
counted in (23) contributes to the coefficient (28). Both cubic input legs
are built from the same root-word function; they are not polarized by the
target coordinates `a,b`. Nor is an address proved to lie in the
degree-`526` visible/Jackson carrier or on THM-2293's canonical `c_3`
edge.

The exact surviving question is now:

```text
given a positive cubic cell C_(k,l),
make one exact relation address in
  k e_a+l e_b+L_13
participate in that very current.                                  (31)
```

This is a coefficient-incidence problem inside an overwhelmingly
nonempty labelled fibre, not a support-existence problem.

## 7. Loss and hostile-boundary ledger

```text
source:
  THM-2309's bright owner packet and target quotient,
  the full septimal scalar kernel,
  THM-2301's essential affine-arrangement theorem;

preserved:
  selected source owner and owner-pivot quotient gauge, both target
  coordinates, all fourteen gains, all twelve nonzero scales on every
  gain, exact relation equations, all nine 91-unit coordinates, and
  controlled exact lift height;

new object:
  a bank in every nonzero vector fibre of K_13/L_13,
  naturally indexed by the complete ordered character pair;

destroyed or unselected:
  bounded visible-carrier membership, an analytic current coefficient,
  current/address incidence, ordinary source frequency, endpoint phase,
  signed amplitude, and canonical c_3-edge membership;

sharp hostile controls:
  |supp_7(w)|=1 gives no all-seven-unit relation anywhere;
  a dark column would turn one affine coordinate condition into an
  identically fixed coordinate and can empty a target fibre;
  THM-2309's bright hostile packet shows that full-lattice abundance does
  not imply bounded selected-owner visibility.                       (32)
```

The theorem makes no sharpness claim for the combined numerical floors.
The septimal support boundary itself is sharp.

No scalar profile is excluded, and LRC(14) remains open.

## 8. Exact verification

The companion reconstructs THM-2309's exact hostile scalar word and
target-grafted packet, checks exact orthogonality, rank six and complete
brightness modulo seven and thirteen, and verifies the two-target quotient.
By exact inclusion-exclusion and modular Gaussian elimination it enumerates
all `168` nonzero mod-thirteen target-vector fibres of that packet; their
actual all-unit counts range from `2,316,060` to `2,526,612`, above (12).

It also evaluates (19) for every support size `1` through `9`, checks all
constants in (23)--(25), and constructs a centered CRT/Bezout exact lift
in the gain-five fibre. Every load-bearing check raises explicitly in
ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_prescribed_gain_91_unit_needle_bank_thm2325.py
python3 -O 04-computation/lrc14_prescribed_gain_91_unit_needle_bank_thm2325.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_prescribed_gain_91_unit_needle_bank_thm2325.out
```

byte-for-byte after LF normalization. QED.
