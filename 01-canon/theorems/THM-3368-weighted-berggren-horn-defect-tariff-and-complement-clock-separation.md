---
id: THM-3368
title: "Weighted Berggren Horn defect tariff and complement-clock separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The THM-3357
  sibling Horn rule has an exact weighted-defect strengthening: the middle
  determinant gate passes whenever the Pythagorean-weighted outer gate debt
  is at most the positive circuit correction K.  This can fire even when
  both outer gates fail.  Complement-clock closure composes lawfully as an
  independent physical-row exit, but cannot be substituted for a Horn seed:
  one fixed THM-3363-closed physical row has saturated Berggren deck
  embeddings whose gate scores tend to minus infinity.  A labelled deck,
  basis gauge, parameter, and support-transfer reconstruction are therefore
  mandatory sidecars.  The result prunes fixed packets; it does not prove
  LRC(14) or propagate a complement certificate through the ternary tree.
source: root/factorial-jacobian-lrc-threebranch-2026-08-14
depends_on:
  - THM-2053-rank-two-parameter-plane-geodesic-terminal
  - THM-2057-scaled-zeta-core-one-tail-closure
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
  - THM-3363-d14-complement-clock-small-lrc-terminal
  - THM-3366-all-sector-complement-clock-completion
related:
  - THM-2055-determinant-gate-normal-fan-and-tangent-sector-reduction
  - THM-2056-kelvin-polar-farey-defect-certificate
script: 04-computation/lrc14_berggren_horn_defect_tariff_thm3368.py
output: 05-knowledge/results/lrc14_berggren_horn_defect_tariff_thm3368.out
script_sha256: 9ec40d5662c4061dd934e17ea99e1d197b19a9bab8a806e084d07b8d9efe721b
output_sha256: a24bfd3f05516e85e416ac9dda5b49fc1676f989301bda2cf5b7c3082cf8457d
hash_basis: working-tree bytes (LF)
---

# THM-3368 -- weighted Berggren Horn debt and clock separation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is a typed composition theorem.  It strengthens the numerical Horn
predicate inside one fixed rank-two gauge, gives an exact clock/Horn packet
that the old pass/pass formulation misses, and proves that a physical-row
clock terminal cannot manufacture the missing gauge data.

## 1. Inherited objects and fixed-gauge convention

Let

```text
u=(m,n)^T,                         0<m<n,                 (1)
x=Lu=(n,2n-m)^T,
y=Mu=(n,2n+m)^T,
z=Ru=(m,2m+n)^T.                                           (2)
```

Write the associated Pythagorean weights

```text
a=n^2-m^2,               b=2mn,               c=n^2+m^2. (3)
```

They are positive and THM-3357 proves the sibling circuit

```text
a x+b z=c y.                                             (4)
```

Fix an ordered finite deck of column vectors

```text
C=(c_i),                         c_i in R^2,              (5)
Delta_C(v)=max_i |det(v,c_i)|.                            (6)
```

For `kappa>=0`, define the determinant-gate score and its negative debt by

```text
G_kappa(v)=||v||^2-kappa Delta_C(v),
A_kappa(v)=max(0,-G_kappa(v)).                            (7)
```

In the LRC application the columns are the labelled integral columns of a
chosen saturated two-plane basis and `kappa=91`.  Both the Euclidean norm and
`Delta_C` are basis-dependent exactly as in THM-2053/2055.  The basis and
the full deck are fixed throughout every use of the tariff below.

## 2. Weighted outer-debt tariff

Put

```text
K(u)=2m(n-m)(3n^2+4mn-m^2)>0.                            (8)
```

Then, for every deck `(5)` and every `kappa>=0`, one has

```text
c G_kappa(y)
 >=a G_kappa(x)+b G_kappa(z)+K(u)                        (9)
 >=K(u)-a A_kappa(x)-b A_kappa(z).                       (10)
```

Consequently,

```text
a A_kappa(x)+b A_kappa(z)<=K(u)
                    ==> G_kappa(y)>=0.                   (11)
```

If the premise in `(11)` is strict, the middle score is strictly positive.
At equality only the weak gate `G_kappa(y)>=0` follows, which is enough for
THM-2053 when `kappa=91`.

### Proof

The maximum of absolute linear functionals in `(6)` is positively
homogeneous and subadditive.  Applying it to `(4)` gives

```text
c Delta_C(y)<=a Delta_C(x)+b Delta_C(z).                 (12)
```

Direct expansion of `(2)--(3)` gives the exact norm correction

```text
c||y||^2-a||x||^2-b||z||^2=K(u).                        (13)
```

Subtract `kappa` times `(12)` from `(13)` to obtain `(9)`.  Since

```text
G_kappa(x)>=-A_kappa(x),
G_kappa(z)>=-A_kappa(z),                                 (14)
```

equation `(10)` follows.  Division by `c>0` proves `(11)`.  QED.

The pass/pass Horn rule of THM-3357 is the special case `A_kappa(x)=
A_kappa(z)=0`.  Formula `(10)`, rather than a Boolean tournament edge, is
the stronger preserved object: the two outer failures may be paid by the
positive circuit correction.

## 3. Exact positive packet beyond pass/pass

Use the saturated AP one-tail deck

```text
c_i=(i,0),       1<=i<=13, i!=12,
c_12=(12,1).                                             (15)
```

Then

```text
Delta_C(r,s)=max(13|s|,|r-12s|).                         (16)
```

Take the primitive opposite-parity Berggren parent

```text
(m,n)=(353,356).                                         (17)
```

Its sibling parameters, weights, and correction are

```text
x=(356,359),       y=(356,1065),       z=(353,1062),
(a,b,c)=(2127,251336,251345),
K=1606017978.                                            (18)
```

At `kappa=91`, exact evaluation of `(16)` gives

```text
(G_91(x),G_91(y),G_91(z))=(-169080,1066,-3893).          (19)
```

Thus both outer gates fail.  Nevertheless their weighted bill is

```text
a A_91(x)+b A_91(z)=1338084208<K,                        (20)
```

and the tariff is sharp on this packet:

```text
c G_91(y)=K-a A_91(x)-b A_91(z)=267933770.               (21)
```

The physical row selected by a positive parameter `(r,s)` in `(15)` is

```text
{r,2r,...,11r,13r,12r+s}.                                (22)
```

THM-2057 proves LRC(14) for every row of this form.  Therefore the two outer
nodes in `(18)` have independent clock/binding exits, while `(11)` certifies
the middle determinant gate.  All three nodes are discharged even though
the original outer-pass/outer-pass rule cannot fire.  The clock exits do not
enter the proof of `(21)`; they close the two physical rows after the
separate numerical tariff has retained their actual gate debts.

## 4. Safety is not a gate seed

The same AP deck at the root `u=(1,2)` gives

```text
(x,y,z)=((2,3),(2,5),(1,4)),
(G_91(x),G_91(y),G_91(z))=(-3536,-5886,-4715).            (23)
```

The three physical tails in `(22)` are respectively

```text
27,                         29,                         16. (24)
```

All three rows are LRC-safe by THM-2057, yet all three determinant gates
fail.  In particular, the tempting replacement

```text
outer rows have clock exits  ==> middle gate passes       (25)
```

is false even inside one fixed saturated plane and one genuine Berggren
sibling packet.  LRC safety is the target predicate; the determinant gate is
only one sufficient certificate for it.

## 5. A fixed THM-3363 row has unbounded negative deck score

There is a sharper obstruction that uses complement-clock completion itself.
Consider the labelled physical row

```text
s=(1,2,3,4,6,12,168,36,60,108,132,156,180).             (26)
```

It is exactly in THM-3363's canonical chart:

```text
F={1,2,3,4,6,12},          L=168,
k=1,                       alpha=1,
D=14,
(d_1,...,d_6)=(14,14,14,14,14,14),
(c_1,...,c_6)=(3,5,9,11,13,15),
(b_1,...,b_6)=(3,5,9,11,13,15).                         (27)
```

Indeed, the aligned tail is `168`, and the six drift speeds are
`168c_i/14`.  THM-3363 therefore proves that `(26)` cannot be an LRC(14)
counterexample.

Now fix the outer Berggren parameter

```text
d=L(1,2)=(2,3).                                         (28)
```

For each component `s_i` of `(26)`, put

```text
epsilon_i=s_i mod 2 in {0,1},
c_i^(0)=((s_i-3epsilon_i)/2,epsilon_i).                  (29)
```

Then `d dot c_i^(0)=s_i`.  The first two columns are `(-1,1)` and `(1,0)`,
whose determinant is `-1`; hence the deck is saturated.  The final column,
corresponding to `s_13=180`, is `(90,0)`.  For an arbitrary integer `q>=0`,
replace it by

```text
c_13^(q)=(90+3q,-2q),                                   (30)
```

leaving all other columns fixed.  Since `d dot (3,-2)=0`, every deck
`C_q` selects the identical physical row `(26)` at the identical Berggren
parameter `(28)`.  The unchanged determinant-`-1` pair proves saturation for
every `q`.

Nevertheless, direct evaluation gives

```text
Delta_(C_q)(d)=270+13q,
G_(91,C_q)(d)=13-91(270+13q)
             =-24557-1183q -> -infinity.                (31)
```

Thus neither the complete physical row nor its THM-3363 complement-clock
certificate determines a gate pass or even a finite lower bound for the gate
score across saturated ambient embeddings.  A fortiori, the coarser
THM-3366 support key cannot supply such a bound.  This is not a defect in the
clock theorem: the physical LRC predicate is basis-free, whereas `(6)--(7)`
are intentionally coordinates of a chosen rank-two proof chart.

## 6. Lawful THM-3366 decision pipeline

For one fixed labelled saturated deck `C` and one Berggren parameter `v`, the
source-to-target map is

```text
v
 -> labelled physical row s(v)=(v dot c_i)_i
 -> witnessed THM-2928 body/tail chart
 -> support key (F,k,D,S_D)
 -> THM-3366 complement-clock exit, when that key is closed. (32)
```

The second arrow in `(32)` is extra data, not an automatic operation on an
arbitrary rank-two row.  It must retain the literal six-body labels, seven
tail labels, aligned multipliers, reduced drift pairs `(c_i,d_i)`, and the
exact reconstruction of the physical speeds.  When the key is in
THM-3366's proved terminal set, the corresponding physical row is LRC-safe
uniformly over the quotient numerators covered by that theorem.

The lawful sibling procedure is therefore:

1. compute `G_91` and the two outer debts in the fixed deck;
2. apply `(10)--(11)` to the middle child;
3. on every child still lacking a determinant certificate, construct an
   actual THM-2928 chart when one exists and query the THM-3366 terminal set;
4. record a clock hit as an independent physical-row exit, never as zero
   gate debt.

This composition can close a sibling packet as in Section 3.  It does not
make clock status hereditary along `L,M,R`, and it supplies no outer gate
seed for another packet.

### Predicate and loss ledger

The combined interface preserves, when its sidecars are kept:

* the fixed-gauge parameter, labelled deck, determinant values, debts, and
  sufficient gate predicate;
* the labelled physical speed row and its body/tail reconstruction;
* THM-3366's pointwise strict-open cover semantics, grid-endpoint owners, and
  exact support-row terminal.

The support quotient by itself forgets:

* the saturated plane, basis gauge, Euclidean parameter norm, determinant
  value, hull owner, and Berggren ancestry word;
* quotient numerators and physical tail locations that THM-3366 deliberately
  eliminates uniformly;
* the lonely time, active runner owner, phase height, and any gate score.

The mandatory join sidecar is therefore

```text
(labelled C, fixed basis, v, sibling role,
 body/tail partition, aligned multipliers,
 reduced drift numerators and denominators,
 exact physical-row reconstruction).                    (33)
```

## 7. Scope boundary

The tariff proves a stronger sufficient determinant certificate inside a
fixed packet.  Complement-clock completion supplies an orthogonal exit for
some physical rows.  Neither predicate is necessary for LRC safety, and
remaining gate- or support-row failures are only unresolved.

No enumeration here shows that every live rank-two parameter has a
THM-2928 chart, that every outer debt satisfies `(11)`, or that every
complement support has a small clock cover.  No owner, phase, or global safe
time propagates through the Berggren tree.  Therefore this theorem does not
prove LRC(14), and it makes no JC, FC(3), tournament-equivalence, or arbitrary
ternary-tree claim.

## Reproduction

```text
python 04-computation/lrc14_berggren_horn_defect_tariff_thm3368.py
python -O 04-computation/lrc14_berggren_horn_defect_tariff_thm3368.py
```

The companion checks the circuit and norm correction symbolically, exercises
the Horn and tariff inequalities on `5,175` exact hostile score cells, and
freezes every integer in Sections 3--5.  Runtime checks remain active under
optimized Python.

**End of proof.**
