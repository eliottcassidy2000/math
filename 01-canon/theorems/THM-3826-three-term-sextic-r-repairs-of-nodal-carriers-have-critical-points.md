---
id: THM-3826
title: "At-most-three-term sextic R-repairs of nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the c=1
  cubic pseudo-plane, every pure-r
  carrier A=e^2-z/3+r g(e) with deg(g)=6 and at most three nonzero monomials
  has a critical point.  On the first full seven-residue cover b_6=s^3,
  T=s^7, boundary-only support would force one of twenty-two exact support
  cells to satisfy four invariant remainder rows; every saturated cell ideal
  is [1].  A genuine four-term cell survives those four rows, so no denser
  sextic or all-degree conclusion is claimed.
source: jc_sparse_direct_search / first full-seven-residue sextic lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn-boundary, 2026-08-23).
  The audit checked the multiplicity-safe boundary implication, the exact
  Euclidean quotient and remainder degree, surjectivity and support
  preservation on the cubic/seventh-power cover, the complete torus
  saturation universe, the repeated-root and four-term hostile controls,
  and every denominator and finite-root exclusion in the critical-point
  reconstruction.  Normal and optimized runs byte-match the frozen
  transcript and the raw hashes agree.  The deterministic companion
  independently
  expands THM-3820's universal residual, verifies degree thirty-two and its
  leading coefficient, certifies the primitive degree-six logarithmic
  quotient with denominator 15625 b_6^5, constructs the finite cubic cover,
  and back-substitutes all five frozen invariant rows.  It exhausts one
  monomial, six binomial, and fifteen genuine trinomial cells by exact
  saturation, verifies the minimal killing-row table and every previous-row
  hostile survivor, includes a repeated-root control and a genuine four-term
  top-four survivor, and replays the finite-root critical reconstruction.
  A second specialization-first cross-check (lrc_reversible_address plus
  sextic_cell_crosscheck, 2026-08-23) independently recovered the same 22-cell
  universe and the five first-three-row survivors
  `{0},{3},{0,3},{1,3},{3,5}`; in each survivor row 28 forces `T=0`.  Its
  census digest is
  `c88a530caf6b16f552da5886bd86c1f01f3cdd5bb2814d843394442df6a4f4a8`.
  The sharp `{3,5}` control
  `g=e^3-(125/16)e^5-(125/16)e^6` kills rows 31--29 but has nonzero row 28,
  confirming the minimal-row rather than merely final-ideal claim.
depends_on:
  - THM-3820-universal-euler-mod-seven-residual-quadratic-discriminant
related:
  - THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3817-quintic-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_three_term_sextic_r_repair_thm3826.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_three_term_sextic_r_repair_thm3826.out
script_sha256: e798664c9f4bba92f35631e50e533168f2c73803d269cc253a3ea6623853f688
output_sha256: e8e1af9c950fcdebc88f4cdb6099fe499ad3fc2f1d761d476fb119de9d2c5cdb
semantic_sha256: 6a8d50a5a0019c7ac29597f69d4939588cdf77c20a469153000937b7a7367e6e
hash_basis: raw LF bytes
---

# THM-3826 -- every three-term sextic r-repair remains critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  On

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r)
```

use the Poisson law

```text
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.      (1)
```

Let

```text
g(e)=sum_(i=0)^6 b_i e^i,            b_6!=0,
A=e^2-z/3+r g(e).                                      (2)
```

If at most three of `b_0,...,b_6` are nonzero, then `A` has a critical point
on `Y`.  Consequently `A` has no regular Darboux mate.

The scope is sharp for the four-row mechanism used below: a genuine
four-term cell survives all four rows.  Sextics with four or more terms,
degree at least seven, mixed corrections, and arbitrary planar Keller maps
remain **OPEN**.

## 1. Universal residual and multiplicity-safe boundary test

Put

```text
u=re,                         K=1+2u,
P=g u^2-K(2e^3+u e g'),
Q=e^2K^3-729g^3u^2(1+u)^2.                            (3)
```

THM-3820 gives the coefficient-free factorization

```text
Res_u(P,Q)=g(e)^3e^4H(e).                              (4)
```

For the full sextic `(2)`, exact expansion gives

```text
deg H=32,                    LC(H)=53144100 b_6^5.     (5)
```

In particular `H` is nonzero and has a root.  Suppose, for contradiction,
that every root of `H` lies on

```text
V(e g(e)).                                               (6)
```

Factor over `k`, retaining arbitrary multiplicities:

```text
H=mu product_alpha (e-alpha)^(n_alpha).                (7)
```

Every `alpha` in `(7)` is a root of `eg`, so

```text
e g H'/H=sum_alpha n_alpha e g/(e-alpha) in k[e].      (8)
```

Thus boundary-only support necessarily implies

```text
H divides e g H'.                                      (9)
```

This implication uses neither squarefreeness of `g` nor radical support.  At
a boundary point where `eg` has order `m>=1` and `H` has order `n>=1`, the
order of `egH'` is at least `m+n-1>=n`.

## 2. Exact quotient and the first full seven-residue cover

Over `Q(b_0,...,b_6)`, divide `egH'` by `H`.  The quotient has degree six
and can be written uniquely in the normalization

```text
q=N/(15625 b_6^5),                                    (10)
```

where `N` is the primitive integer polynomial frozen in the companion,

```text
deg_e N=6,                 LC_e(N)=500000 b_6^6,
content_Z(N)=1,            gcd(N,b_6)=1.              (11)
```

Define the cleared remainder

```text
R=15625 b_6^5 e g H'-N H=sum_(j=0)^31 r_j e^j.       (12)
```

The degree bound in `(12)` independently certifies that `(10)` is the exact
Euclidean quotient.  If `(9)` held, every `r_j` would vanish.

The leading coefficient has mod-seven weight `3`.  Since `k` is
algebraically closed, choose `s!=0` with

```text
b_6=s^3.                                               (13)
```

Put

```text
C_i=b_i s^(3-i)       (0<=i<=5),
T=s^7.                                                    (14)
```

Equivalently,

```text
b_i=C_i s^(i-3),      b_6=s^3,       T!=0.            (15)
```

This finite cover is surjective on `b_6!=0`.  It preserves exact support:
`b_i!=0` if and only if `C_i!=0`.  The top five pulled-back rows have the
form

```text
r_j|_(15)=lambda_j s^(j-5) F_j(C_0,...,C_5,T),
j=31,30,29,28,27,                                  (16)
```

where every `lambda_j in Q*`.  Thus the Laurent powers are

```text
26,25,24,23,22,                                       (17)
```

and the primitive row term counts are

```text
20,29,41,54,73.                                       (18)
```

The companion back-substitutes `T=s^7` in every row and verifies the
nonzero scalar `lambda_j`; no coefficient equation or `T=0` component is
lost.

## 3. Complete exact-support census

Let `J` be the set of nonzero lower coefficients among `C_0,...,C_5`.  Since
`C_6` is represented by the nonzero cover parameter `s`, the hypothesis
`|supp(g)|<=3` is exactly

```text
J subset {0,...,5},                    |J|<=2.         (19)
```

There are

```text
1+6+binom(6,2)=22                                      (20)
```

such exact-support cells.  For a cell `J` and `1<=m<=4`, set all
`C_i=0` for `i notin J` and form

```text
I_(J,m)=<F_31,...,F_(32-m),
          1-v T product_(j in J) C_j>.                (21)
```

The last equation is the exact affine chart for
`D(T product_(j in J)C_j)`; it simultaneously enforces the cubic-cover torus
and every claimed nonzero support coefficient.  No division is performed.

Exact Groebner reduction over `Q` gives `I_(J,m_J)=[1]` with the following
minimal row count `m_J`:

```text
J       empty  0  1  2  3  4  5
m_J       3    4  2  1  4  1  2

J       01 02 03 04 05 12 13 14 15 23 24 25 34 35 45
m_J      2  1  4  1  3  1  4  3  3  2  3  3  3  4  3. (22)
```

For every entry with `m_J>1`, the companion also checks

```text
I_(J,m_J-1) != [1].                                   (23)
```

Thus `(22)` is not a vacuous saturation or an indexing shift: each previous
row packet has a genuine algebraic survivor on the stated torus, and the
displayed next row is the first exact exit.  All twenty-two cells are empty
once the four necessary rows are imposed.  This contradicts `(9)` and proves
that `H` has a root away from `(6)`.

The census includes repeated roots.  For example,

```text
g=e^4(e-1)^2=e^6-2e^5+e^4                            (24)
```

lies in cell `J={4,5}`.  At `C_4=1,C_5=-2,T=1`, its first four primitive
rows are

```text
30286, -1200098, 4224804, -1419167251,                (25)
```

so no hidden squarefree filter enters the computation.

## 4. A genuine four-term boundary survivor

The support bound is structural, not cosmetic.  In the four-term cell

```text
J={0,1,4},             supp(g)={0,1,4,6},              (26)
```

the saturated ideal of the same top four rows is not the unit ideal.  With
variable order `(v,T,C_0,C_1,C_4)`, its exact Groebner basis is

```text
157464C_4+1953125T^2,
625C_4T+324,
486C_4^2-3125T,
48828125T+209952v,
10C_0-1,
25C_1+4C_4.                                           (27)
```

Indeed `(27)` gives

```text
C_4^3=-10/3,
T=-324/(625C_4),
C_0=1/10,
C_1=-4C_4/25,                                         (28)
```

so every active coordinate and `T` is nonzero over `k`.  This is a genuine
four-term hostile survivor of the complete mechanism used in `(22)`.  The
next row happens to reduce to

```text
F_27 == -911250 C_4 mod I_(J,4),                      (29)
```

and therefore kills this particular cell, but `(29)` does not close the
other four-term cells.  They remain open rather than being inferred from a
bounded table.

## 5. From a residual root to a critical point

Choose a root `eta` of `H` outside `(6)`.  Then

```text
eta g(eta)!=0.                                         (30)
```

By `(4)`, `P,Q` have a common finite root `u_0`: the leading coefficient of
`Q` in `u` is `-729g(eta)^3!=0`, so the resultant has no projective root at
infinity.  Moreover,

```text
Q(0)=eta^2,
Q(-1)=-eta^2,
Q(-1/2)=-729g(eta)^3/16,                              (31)
```

so, on putting `K_0=1+2u_0`, one has

```text
u_0(1+u_0)K_0!=0.                                    (32)
```

Define without choosing a square or cube root

```text
r_0=u_0/eta,
z_0=9g(eta)u_0(1+u_0)/(eta K_0).                     (33)
```

Exact reduction gives

```text
z_0^2-K_0/(9g(eta)) = -Q/(9g(eta)eta^2K_0^2),
r_0^2eta-z_0^3+r_0 = u_0(1+u_0)Q/(eta^3K_0^3),
{A,e}|_0 = -Q/(eta^2K_0^2),
{A,z}|_0 = 3P/eta^2.                                  (34)
```

Thus `(r_0,z_0,eta)` lies on `Y`, and its last two Hamiltonian components
vanish.  Finally the surface-Casimir identity

```text
(1+2re){A,r}-3z^2{A,z}+r^2{A,e}=0                    (35)
```

and `K_0!=0` kill `{A,r}`.  The point is critical, proving the theorem.

## 6. Exact scope and replay

This is the first pure-r result at degree six and therefore the first one in
which all seven Euler residue sectors can occur simultaneously.  Its content
is exactly the sparse support floor `(19)`, not a general sextic closure.
The four-term witness `(27)` records the first lost sidecar explicitly.

Run

```bash
python3 04-computation/jc2_cubic_pseudoplane_three_term_sextic_r_repair_thm3826.py
python3 -O 04-computation/jc2_cubic_pseudoplane_three_term_sextic_r_repair_thm3826.py
```

Both executions must byte-match

```text
05-knowledge/results/jc2_cubic_pseudoplane_three_term_sextic_r_repair_thm3826.out
```

The companion uses no inactive Python `assert`, freezes the five-row packet,
all saturated basis transcripts, the exact support universe, and the semantic
scope, and reports `CHECKS=320` and `RESULT=PASS`.
