# Strict component-width consumers for the dyadic residuals

**Status: PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED; LRC(14)
OPEN.** The q=2 constants are consequences of THM-4451; the q=4 theorem
below is independently typed and does not pretend that the all-odd q=2 lane
inequality applies to an even tail.

Put

```text
D_n={x in R/Z: ||nx||<1/14},
G_C={y in R/Z: ||cy||>=1/14 for every c in C}.
```

All danger sets are open and all body-safe sets are closed.  A row is called
safe when its maximum clearance is at least `1/14`; therefore its strict
failure set uses `<1/14`.

## 1. Clock two, zero even tails: the factor is two

For three distinct odd tails `T`, let

```text
P_T=(union_(t in T) D_t) intersect
    ((union_(t in T) D_t)-1/2)
```

be the physical `x`-failure set, and let `Q_T` be its quotient under `y=2x`.
Equivalently, `y in Q_T` precisely when both lifts `y/2,(y+1)/2` are killed
by `T`.  The set `P_T` is half-turn invariant.  Its antipodal components pair
under doubling, and every quotient component has twice the physical width.

For a ten-body `C`, strict failure of

```text
2C union T
```

is equivalent to the correctly typed inclusion

```text
G_C subset Q_T.                                           (1)
```

Indeed, over every `y in G_C` neither lift is killed by the doubled body, so
both lifts must be killed by the tails.  Conversely, points outside `G_C`
are already killed by `2C`.

Let `L_C` be the largest length of a positive-length connected component of
`G_C`.  The physical component theorem in THM-4451 and (1) give:

> **q=2 all-odd component entry.**  For every three distinct positive odd
> tails, the row `2C union T` is safe if
>
> ```text
> L_C >= 34/693.                                         (2)
> ```
>
> A strict counterexample necessarily has `L_C<34/693`.

> **q=2 odd-3-unit component entry.**  If no tail is divisible by three, the
> row is safe if
>
> ```text
> L_C >= 38/1001.                                        (3)
> ```
>
> A strict counterexample necessarily has `L_C<38/1001`.

The factors in (2)--(3) are

```text
2*(17/693)=34/693,       2*(19/1001)=38/1001.
```

Using `17/693` or `19/1001` directly against a quotient body component would
lose a factor of two and is invalid.

The inequalities in (2)--(3) are inclusive.  Under failure, a closed
connected body component lies inside one open quotient component.  A closed
interval of length at least the length of an open interval cannot be a
subset of it.  In particular equality cannot be hidden by the fact that the
tail set itself is open.  This is the component analogue of compact-in-open
strictness, not an a.e. measure argument.

More precisely, if `lambda(T)` denotes the actual longest physical component
of `P_T`, then failure forces

```text
L_C < 2 lambda(T).                                       (4)
```

For odd common scale `d`, `P_(dT)=m_d^(-1)(P_T)`, so
`lambda(dT)=lambda(T)/d`.  Thus (4) retains a useful scale coordinate.  On
the two extremal shapes it reads

```text
T=d*(1,9,11):    d L_C < 34/693,
T=d*(1,11,13):   d L_C < 38/1001     (odd 3-unit lane). (5)
```

The universal equalities themselves occur only at scale `d=1`; dilation
shortens components and does not create an equality ray for the universal
bound.

### Body component-count form

Let `kappa_+(C)` be the number of positive-length connected components of
`G_C`; isolated safe points are retained in exact address work but contribute
no Haar mass.  Since one component has length at least
`mu(G_C)/kappa_+(C)`, (2)--(3) imply the inclusive gates

```text
mu(G_C) >= (34/693) kappa_+(C)       for arbitrary odd T,
mu(G_C) >= (38/1001) kappa_+(C)      for odd 3-unit T.    (6)
```

Equivalently a strict failure must satisfy

```text
kappa_+(C) > (693/34) mu(G_C),
kappa_+(C) > (1001/38) mu(G_C),                         (7)
```

in the respective domains.  These are component gates, not assertions of a
universal ten-body mass floor.

## 2. Clock four with one `v_2=1` tail

Let `r,a,b` be positive odd integers with `a!=b`, and define the strict
physical four-sheet tail failure set

```text
P4_(r;a,b)={x: every x+j/4, j=0,1,2,3,
               is killed by one of 2r,a,b}.             (8)
```

An odd tail kills at most one of the four quarter-lifts.  The tail `2r`
kills at most one parity pair, hence at most two lifts.  Covering all four
lifts saturates the exact capacity `2+1+1=4`: `2r` owns one parity pair and
`a,b` own the two remaining lifts, one each.  If

```text
P_ab=(D_a union D_b) intersect ((D_a union D_b)-1/2),
```

then this sheet-coloured saturation gives the exact identity

```text
P4_(r;a,b)
 =[D_(2r) intersect (P_ab-1/4)]
  union [(D_(2r)-1/4) intersect P_ab].                  (9)
```

The two terms in (9) are disjoint.  Also

```text
P_ab=[D_a intersect (D_b-1/2)]
     union[D_b intersect (D_a-1/2)],                    (10)
```

and its two oriented terms are disjoint, because an odd tail cannot be
dangerous on both antipodal sheets.  Consequently every component in (9)
lies in one `D_(2r)` tooth and in one tooth of the larger of `a,b`.  Therefore

```text
lambda_4(r;a,b)
 <= min(1/(14r), 1/(7 max(a,b))).                       (11)
```

This is the q=4 version of lane capacity with the sheet colour retained.
For comparison, the literal scalar analogue uses

```text
R_r=D_(2r) union (D_(2r)-1/4),
Q_n=union_(j=0)^3 (D_n-j/4).
```

These are one-tooth combs of frequency/duty `(4r,2/7)` and `(4n,4/7)`.
For such a comb put, when `kL=m+s` with integral `m` and `0<=s<1`,

```text
Cap_(k,delta)(L)=[m delta+min(s,delta)]/k.
```

The uncoloured lane inequality for a failure interval is only

```text
3L <= Cap_(4r,2/7)(L)
      +Cap_(4a,4/7)(L)+Cap_(4b,4/7)(L).                 (11a)
```

It is not enough: at length `1/110`, (11a) is already an equality for the
genuinely empty carrier `(r,a,b)=(1,1,5)`.  The parity address in (9), not
another Haar sum, is what makes the bound sharp.

### Sharp all-height widths

Equation (11) reduces a violation of `1/98` in the all-odd domain to
`r in {1,3,5}` and `a,b in {1,3,5,7,9,11,13}`.  Including `r=7` classifies
equality.  Exact strict-wall evaluation of these 84 rows gives

> **q=4 all-odd component theorem.**
>
> ```text
> lambda_4(r;a,b) <= 1/98,                              (12)
> ```
>
> with equality, up to swapping `a,b`, exactly at
>
> ```text
> (r;a,b)=(1;7,9), (3;7,13), (5;3,7).
> ```

In the odd-3-unit domain, (11) reduces a violation of `1/110` to

```text
r in {1,5,7},        a,b in {1,5,7,11,13}.
```

Exact strict-wall evaluation of these 30 rows gives

> **q=4 odd-3-unit component theorem.**
>
> ```text
> lambda_4(r;a,b) <= 1/110,                             (13)
> ```
>
> with equality, up to swapping `a,b`, uniquely at
>
> ```text
> (r;a,b)=(5;1,11).
> ```

For visibility, the only nonzero rows in the 30-row residual box are

```text
(1;11,13):17/2002,
(5;1,11):1/110,       (5;1,13):1/910,
(7;1,11):1/539,       (7;1,13):5/637.                  (14)
```

The remaining 25 are empty.  These are exact physical widths, not a finite
height conjecture.  The infinite exterior is closed analytically by (11).

## 3. Clock-four LRC consumers and the factor four

Let `Q4_(r;a,b)` be the quotient of (8) under `y=4x`.  The physical failure
set is quarter-turn invariant and proper; every quotient component has
exactly four times the physical width.  Strict failure of

```text
4C union {2r,a,b}
```

is equivalent to

```text
G_C subset Q4_(r;a,b).                                  (15)
```

Thus (12)--(13), with compact-in-open equality handled exactly as above,
give:

> **q=4 original-body component entries.**  The row is safe whenever
>
> ```text
> L_C >= 2/49                 for arbitrary odd r,a,b,
> L_C >= 2/55                 for odd 3-unit r,a,b.      (16)
> ```

> A strict failure forces the corresponding strict reverse inequality.

The live one-`v_2` residual uses the second line.  It is a direct gate on the
original ten-body and is logically different from applying THM-4153 to the
absorbed body `H=2C union {r}`.  It does not close the residual universally,
because no universal theorem here forces an arbitrary `G_C` component to
have length `2/55`.

If `kappa_+(C)` is as above, the corresponding mass/count gates are

```text
mu(G_C) >= (2/49) kappa_+(C)       (all odd),
mu(G_C) >= (2/55) kappa_+(C)       (odd 3-unit).          (17)
```

For an odd common dilation `d`, quarter sheets are merely permuted and
`lambda_4(dr;da,db)=lambda_4(r;a,b)/d`.  On the odd-3-unit extremal shape,
failure therefore forces

```text
d L_C < 2/55       for (r;a,b)=d*(5;1,11).              (18)
```

Again, equality in the universal physical cap occurs only at `d=1`.

## 4. Exact audit and scope

The audit compares (9) with an independent direct four-sheet wall predicate,
retains deleted strict endpoints, checks every row in both finite boxes, and
uses explicit runtime checks under optimized Python.

```powershell
python -B 04-computation/lrc14_clock_four_one_v2_component_thm4452.py
python -O -B 04-computation/lrc14_clock_four_one_v2_component_thm4452.py
python -B 04-computation/lrc14_clock_four_one_v2_component_thm4452_independent.py
python -O -B 04-computation/lrc14_clock_four_one_v2_component_thm4452_independent.py
```

The outputs are line-identical and end in `PASS`.  The q=2 inputs remain the
separately proved all-height `17/693` and `19/1001` physical theorems in
THM-4451. No q=2 or q=4 component gate is a proof of LRC(14): each is a
sharp, correctly normalized consumer condition on the actual body geometry.

```text
04-computation/lrc14_clock_four_one_v2_component_thm4452_geometry.py
  ca1497c66ec0a7dbfc8c9103409e3cddba283abd9a6149088a9b8f6c914a4d65
04-computation/lrc14_clock_four_one_v2_component_thm4452.py
  03432033c64bdd46f50788f01dbe8545e89e31fccac3b2c2bac23900418548f3
05-knowledge/results/lrc14_clock_four_one_v2_component_thm4452.out
  3f11751c9bff00ecaf39184f811a96d64b700fc14d5b12ad9b6a262bbaa5d8d6
```

The independent audit is
`05-knowledge/results/lrc14_clock_four_one_v2_component_thm4452_independent_audit.md`.
It reconstructs the
four-sheet predicate without importing the primary implementation, verifies
the pointwise identity at walls as well as cells, reproduces the 84/30 boxes,
and directly checks the two-component quotient equality controls.  Its raw
SHA-256 values are

```text
04-computation/lrc14_clock_four_one_v2_component_thm4452_independent.py
  6ead91cfa71c5c31fd212d8fb12cf67fc834e968f74c08c72fe97e163b75197e
05-knowledge/results/lrc14_clock_four_one_v2_component_thm4452_independent.out
  c9e830740e0bf0123244bbf92c2b4bf8fc77bc6213a6195974843a70dd6f3f01
```
