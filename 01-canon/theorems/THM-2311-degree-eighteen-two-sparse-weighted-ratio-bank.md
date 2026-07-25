---
id: THM-2311
title: "Degree-eighteen two-sparse weighted-ratio bank"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the genuine nonsplit polynomial
  exact-square-prefix degree-eighteen branch of THM-2262/2297, every
  survivor with exactly two nonzero translation invariants lies in a finite
  bank of 31 weighted-projective ratio points. The C--W plane is already
  empty by THM-2297. On the B--C, B--D, B--W, C--D, and D--W planes, exact
  repeated-branch resultants reduce respectively to 9, 6, 8, 4, and 4
  distinct nonzero algebraic ratios. All displayed ratio polynomials are
  squarefree and pairwise coprime within their plane. This is a finite
  necessary spectral bank, not an emptiness theorem: normalization,
  third-flux, Keller one-form, and whole-Faber sidecars still have to close
  each point. It does not prove JC(2).
source: codex-2026-07-25-degree18-two-sparse-ratios
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
related:
  - THM-2241-monic-transverse-response-depth-and-resultant-nonproper-quotient
  - THM-2285-centered-grid-footprint-and-generic-keller-lines
script: 04-computation/jc2_degree18_two_sparse_ratio_bank_thm2311.py
output: 05-knowledge/results/jc2_degree18_two_sparse_ratio_bank_thm2311.out
script_sha256: 7178b6ba19fe6528412075651d680b7105409e4dd5fcdada11d55de5fa4a553a
output_sha256: 6ee690628aa40697720ac25d932ab5480ecebf850423da6e77e49a5d7cc67dba
hash_basis: working-tree bytes (LF)
---

# THM-2311 -- the two-sparse degree-eighteen locus is a 31-ratio bank

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2297 removes target translation from the degree-eighteen spectral
problem and closes every one-sparse stratum plus the whole `C`--`W` plane.
The remaining weighted cone has coordinates

```text
(B,C,D,W),                  weights (2,3,4,5).       (1)
```

There are five other coordinate planes. On each one the necessary
repeated-branch condition is a weighted homogeneous binary resultant, hence
a one-variable algebraic condition after quotienting by the scaling.
Factoring all five resultants gives a finite bank.

## 1. Inherited necessary condition

After the legal translation of THM-2297, the spectral cubic is

```text
G_0
 =-5878656Wy
  -26040609u^3+49601160Bu^2+1607445u^2y^2
  -20995200B^2u-2857680Buy^2-52907904Du-138915uy^4
  +777600B^2y^2+33592320BD-5598720BCy+78120By^4
  +1959552Dy^2-435456Cy^3+1127y^6.                  (2)
```

Put

```text
Delta_0(y;B,C,D,W)=Disc_u G_0.                      (3)
```

Its `y`-degree is twelve and its leading coefficient is a nonzero constant.
THM-2262 proves that a rational Keller trajectory on a squarefree branch
cover would force a genus-four contradiction. Therefore every survivor in
the retained branch satisfies

```text
Res_y(Delta_0,partial_y Delta_0)=0.                  (4)
```

Equation (4) is necessary, not sufficient. It forgets the component carrying
the trajectory and all differential and polynomial sidecars.

Suppose exactly two of `(B,C,D,W)` are nonzero. The six possible supports are

```text
BC, BD, BW, CD, CW, DW.                              (5)
```

THM-2297 closes `CW`, including both axes and every off-axis critical-value
collision. Its axis arguments also close every one-coordinate boundary of
the other five planes. Hence only the off-axis factors of (4) matter below.

## 2. Weighted ratios and the lifting convention

Use the invariant ratios

```text
BC: t=C^2/B^3,

BD: t=D/B^2,

BW: t=W^2/B^5,

CD: t=D^3/C^4,

DW: t=W^4/D^5.                                      (6)
```

Each is unchanged by

```text
(B,C,D,W)
 ->(rho^2B,rho^3C,rho^4D,rho^5W).                  (7)
```

Over the algebraically closed constant field, the ratio in (6) classifies
the off-axis scaling orbit on the corresponding weighted projective line.

If

```text
t=Y^r/X^s
```

and `p(t)` has degree `d`, write

```text
H_p(X,Y)=X^(sd)p(Y^r/X^s).                          (8)
```

The exact resultant factorization on each plane is a nonzero scalar times
the stated axis monomial and the lifted `H_p` factors with the displayed
multiplicities. Multiplicity is recorded for reproducibility; it is not
being interpreted as a normalization-component count.

## 3. The five exact ratio banks

### 3.1 The B--C plane: nine ratios

Here `t=C^2/B^3`. Apart from the axis factor `B^24`, the resultant is the
product of the lifts of

```text
p_BC,1(t)
 =2000+15309t,                                      multiplicity 1;

p_BC,2(t)
 =125+1134t,                                        multiplicity 2;

p_BC,3(t)
 =3321125-161754894t-2812385772t^2,                 multiplicity 3;

p_BC,4(t)
 =410644531250000
  -18114791748046875t
  -545436228093750000t^2
  -4951165276923468750t^3
  -18946967714644599000t^4
  -26529827304546537363t^5,                         multiplicity 1.
                                                            (9)
```

The degrees sum to

```text
1+1+2+5=9.                                         (10)
```

### 3.2 The B--D plane: six ratios

Here `t=D/B^2`. There is no axis monomial. The factors are

```text
p_BD,1(t)
 =4075-85176t,                                      multiplicity 6;

p_BD,2(t)
 =25-126t,                                          multiplicity 19;

p_BD,3(t)
 =22656250
  -772734375t
  +7600635000t^2
  -30805790400t^3
  +46376717184t^4,                                  multiplicity 2.
                                                            (11)
```

Their degrees give

```text
1+1+4=6.                                           (12)
```

### 3.3 The B--W plane: eight ratios

Here `t=W^2/B^5`. Apart from the axis factor `B^6`, the factors are

```text
p_BW,1(t)
 =5313800000
  +4508659468656t
  -136738899331083t^2,                              multiplicity 3;

p_BW,2(t)
 =5511577600000000000000000000
  +4983290602536960000000000000000t
  -6564822237254419568640000000000t^2
  -3094052863483309848285092659200000t^3
  -81862566455344350924421142159812608t^4
  -744088924275617882256518828471658624t^5
  -2973811237322720333634598763466407943t^6,        multiplicity 1.
                                                            (13)
```

Thus

```text
2+6=8.                                             (14)
```

### 3.4 The C--D plane: four ratios

Here `t=D^3/C^4`. Apart from the axis factor `D^15`, the factors are

```text
p_CD,1(t)
 =22143375+6397664t,                                multiplicity 3;

p_CD,2(t)
 =387420489
  -8964338040t
  +54880100352t^2
  +16544432128t^3,                                  multiplicity 1.
                                                            (15)
```

Therefore

```text
1+3=4.                                             (16)
```

### 3.5 The D--W plane: four ratios

Here `t=W^4/D^5`. Apart from the axis factor `D^3`, the factors are

```text
p_DW,1(t)
 =935886848+430565625t,                             multiplicity 3;

p_DW,2(t)
 =36028797018963968
  +17932072576352256t
  -1448500838400000t^2
  +56162900390625t^3,                               multiplicity 1.
                                                            (17)
```

Again,

```text
1+3=4.                                             (18)
```

## 4. Exact count and weighted-projective consequence

Every polynomial in (9), (11), (13), (15), and (17):

1. has nonzero constant and leading coefficient;
2. is squarefree over `Q`; and
3. is coprime to every other ratio polynomial on its plane.

It therefore has exactly its degree many distinct nonzero roots over `C`,
with no overlap inside the same plane. The complete off-axis candidate
count is

```text
BC  9
BD  6
BW  8
CD  4
DW  4
------
   31.                                              (19)
```

Consequently:

> **Two-sparse weighted-ratio theorem.** Every degree-eighteen survivor in
> the scope of THM-2262/2297 with exactly two nonzero translation invariants
> is represented by one of the `31` roots in (9), (11), (13), (15), and
> (17), on its labelled coordinate plane.

This is an exact finite spectral bank, not merely a degree bound. The plane
label remains part of the state; ratios from different planes are not being
identified.

## 5. What remains at each ratio

For each algebraic ratio, weighted scaling fixes a convenient representative
and leaves an explicit trigonal curve

```text
G_0(u,y)=0.                                         (20)
```

A closure still has to prove at least one of:

```text
the normalization has positive genus;

the rational trajectory is confined to a constant/singular component;

the third flux or Keller one-form has incompatible valuation;

the whole-polynomial Faber sidecar has an uncancellable pole.            (21)
```

The repeated-branch resultant alone supplies none of (21). In particular,
a high factor multiplicity in (11) does not prove a rational component, and
one repeated branch value does not select the component used by a putative
Keller trajectory.

The structural gain is a clean split:

```text
exactly two nonzero invariants -> 31 algebraic cases;

three or four nonzero invariants -> positive-dimensional singular locus.
                                                            (22)
```

The natural next computation is normalization plus the retained differential
sidecars over each ratio number field, beginning with the linear factors in
the `BD`, `CD`, and `DW` banks.

## 6. Connection and loss ledger

```text
source:
  THM-2297's weighted degree-eighteen cone and THM-2262's necessary
  repeated-branch condition;

map:
  restrict to a labelled coordinate plane, quotient the weighted scaling,
  and factor the exact binary repeated-branch resultant;

preserved:
  plane label, invariant ratio, branch-discriminant collision, factor
  multiplicity, and all original sidecar obligations;

destroyed:
  branch point, normalization component, rational trajectory, flux,
  Keller one-form, and whole-polynomial Faber data;

restoring work:
  component-aware normalization and sidecar audit at each algebraic ratio;

hostile control:
  the raw C--W resultant vanishes identically even though THM-2297 closes
  the entire plane by normalization genus.                              (23)
```

The hostile control explains why factorization is a routing theorem rather
than a proof of JC(2).

## 7. Exact reproduction

The companion reconstructs (2), computes `Delta_0`, computes all six plane
resultants, and verifies:

- the five complete lifted factor dictionaries, including axis factors and
  multiplicities;
- the invariant substitutions (6);
- the exact ratio polynomials (9), (11), (13), (15), and (17);
- squarefreeness, pairwise coprimality, and nonzero roots;
- the counts `9+6+8+4+4=31`; and
- the identically vanishing raw `C`--`W` resultant.

Reproduce with

```bash
python3 04-computation/jc2_degree18_two_sparse_ratio_bank_thm2311.py
python3 -O 04-computation/jc2_degree18_two_sparse_ratio_bank_thm2311.py
```

Both runs are byte-identical to the stored output, with all checks active
under optimized Python. The implication from a survivor to the
repeated-branch condition is the mathematical theorem of THM-2262, not a
computer assumption. This theorem remains scoped to the genuine nonsplit,
polynomial exact-square-prefix, reduced degree-eighteen branch and does not
prove JC(2). QED.
