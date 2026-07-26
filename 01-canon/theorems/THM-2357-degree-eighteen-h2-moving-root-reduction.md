---
id: THM-2357
title: "Degree-eighteen H2 moving-root and factorization-pivot reduction"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Every
  degree-eighteen H_2 S_5^2 survivor has a unique finite three-cycle
  branch value r at which P(r)=Q(r)=0. Weighted scaling moves it to
  y=1 and gives explicit two-parameter formulas for D and W. After
  removing (y-1)^2, the degree-ten polynomial R_10 satisfies
  R_10=H_2 S_4^2. Its direct coefficient factorization solves four
  variables triangularly and leaves six equations in B,C,s3,s2; the
  first is C-linear with pivot
  -10260B+1771s2-1771. The complementary pivot wall has one explicit
  irreducible degree-five residual. This is a reduction, not an
  emptiness theorem: the pivot wall, K=0 singular-order boundary,
  C=0 boundary, remaining localized ideal, H_4 stratum, and JC(2)
  remain open.
source: codex-2026-07-25-degree-eighteen-h2-moving-root
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - THM-2347-degree-eighteen-double-zero-wall-saturation
related:
  - THM-2335-degree-eighteen-cyclic-square-class-stratum-empty
  - THM-2338-degree-eighteen-deep-common-root-wall-hurwitz-quartet
script: 04-computation/jc2_degree18_h2_factorization_pivot_thm2357.py
output: 05-knowledge/results/jc2_degree18_h2_factorization_pivot_thm2357.out
script_sha256: 110e7c116b9090745cc5434e21afc02747f93fdbfde68e4b02f1bf873ddc82dd
output_sha256: 58ea7044d4d5b0773c2fe4f6f68a77ba0ebf5398b90b7f6c692f4c6038ee3df3
hash_basis: working-tree bytes (LF)
---

# THM-2357 -- move the mixed branch to one and expose its first pivot

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2332 reduces the residual degree-eighteen problem to

```text
F=4P^3+49Q^2=H S^2,             deg(H) in {2,4}.        (1)
```

This theorem treats the mixed case `deg(H)=2`.  Its gain is a canonical
moving-root gauge and a low-degree coefficient system.  It does not prove
that the resulting system is empty.

## 1. The unique three-cycle gives a common root

The `H_2 S_5^2` case has branch signature

```text
(3,2,2).                                             (2)
```

Thus the normalization has one ramification point of index three.  Infinity
is unramified by THM-2332, so its branch value `r` is finite.

At this fibre the depressed cubic

```text
v^3+p(y)v+q(y)=0
```

has a triple root.  A triple root of a depressed cubic is zero, hence

```text
p(r)=q(r)=0,

P(r)=Q(r)=0.                                        (3)
```

If `r=0`, then `P(0)=0`, so the trajectory lies on the wall already
closed by THM-2345.  Every surviving mixed branch therefore has `r!=0`.

The point is smooth on the normalization, but it need not be smooth in the
polynomial-order plane model.  In particular, (3) does **not** by itself
imply `Q'(r)!=0`; that distinction is retained below.

Since `H` records the two transposition branch values, the three-cycle
value occurs in the square factor.  The order/maximal-order discriminant
formula used in THM-2332 adds only an even index order.  Consequently

```text
(y-r) divides S.                                    (4)
```

## 2. Exact common-root incidence

Retain THM-2332's covariants

```text
P
 =245y^4+1890By^2-24300B^2+122472D,

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y.                           (5)
```

For `r!=0`, solving (3) gives exactly

```text
122472D
 =24300B^2-1890B r^2-245r^4,

W
 =-(20/21)BC-(2/27)C r^2
   -r^3(91r^2+1215B)/177147.                        (6)
```

Indeed, after using `P(r)=0`, the equation `Q(r)/r=0` reduces to

```text
7r^3(91r^2+1215B)
 +91854Cr^2+59049(20BC+21W)=0.                     (7)
```

Differentiating on the incidence surface gives the separate plane-model
smoothness factor

```text
Q'(r)
 =14r^2(245r^3+2835Br+26244C).                     (8)
```

The weighted action sends a branch value `r` to `lambda r`.  Since `r!=0`,
choose `lambda=1/r`.  In the resulting gauge,

```text
r=1,

D=(24300B^2-1890B-245)/122472,

W=-(20/21)BC-(2/27)C-(91+1215B)/177147,

K:=245+2835B+26244C,

Q'(1)=14K.                                         (9)
```

The two already-closed parameter walls pull back to

```text
126D-25B^2
 =-(35/972)(54B+7),

20BC+21W
 =-(7/59049)(13122C+1215B+91).                    (10)
```

Thus `54B+7=0` and `13122C+1215B+91=0` carry no Keller
trajectory by THM-2345 and THM-2347.  The factors `C=0` and `K=0` are
separate boundary charts, not closed here.

## 3. The degree-ten mixed residual

Equations (3)--(4) give exact divisions

```text
P=(y-1)p_3,

Q=(y-1)q_5,

F=(y-1)^2 R_10,                                    (11)

R_10=4(y-1)p_3^3+49q_5^2=H_2 S_4^2.               (12)
```

The degree and leading coefficient are

```text
deg(R_10)=10,

lc(R_10)=73060029.                                  (13)
```

On the smooth-order chart `K!=0`,

```text
R_10(1)=49Q'(1)^2=9604K^2 !=0.                    (14)
```

At `K=0`, (12) remains valid, but `R_10` retains an additional root at
one.  That singular-order total-branch boundary remains open.

A direct subresultant encoding of (12) is

```text
deg gcd(R_10,R_10')>=4.                            (15)
```

Hence all ten coefficients of `Sres_0,...,Sres_3` vanish.  The first
companion computes them exactly.  Their total degrees range from `20` to
`37` before localization and from `18` to `34` after stripping the only
detected closed-wall powers.  Four selected coefficients contain
`(54B+7)^3`; none contains a positive power of `C`, `K`, or the
double-zero pullback.  This is an exact encoding, but its monolithic
Groebner basis is not the efficient next representation.

## 4. Direct factorization gives six small equations

Write

```text
R_10=sum_(j=0)^10 a_j y^j,            a_10=L=73060029,

S_4=y^4+s_3y^3+s_2y^2+s_1y+s_0,

H_2=Ly^2+h_1y+h_0.                                  (16)
```

Comparing coefficients from the top solves successively

```text
h_1
 =a_9-2Ls_3,

h_0
 =a_8-L(s_3^2+2s_2)-2h_1s_3,

s_1
 =[a_7-2Ls_3s_2-h_1(s_3^2+2s_2)-2h_0s_3]/(2L),

s_0
 =[a_6-L(2s_3s_1+s_2^2)
    -2h_1(s_1+s_3s_2)-h_0(s_3^2+2s_2)]/(2L).       (17)
```

The remaining coefficients `y^5,...,y^0` give six primitive equations

```text
E_5,...,E_0 in Q[B,C,s_3,s_2].                     (18)
```

Their exact signatures

```text
(total degree, deg_B, deg_C, deg_s3, deg_s2, terms)

(5, 2,1, 5,2, 22),
(6, 3,2, 6,3, 37),
(7, 3,2, 7,3, 43),
(8, 4,2, 8,4, 76),
(9, 4,2, 9,4,105),
(10,5,2,10,5,126)                                  (19)
```

are far smaller than the terminal-subresultant representation.

## 5. The first coefficient exposes a new pivot wall

The first residual is linear in `C`.  Its exact leading coefficient is

```text
coeff_C(E_5)
 =52488(-10260B+1771s_2-1771).                     (20)
```

Put

```text
L_piv=-10260B+1771s_2-1771.                        (21)
```

On `L_piv!=0`, equation `E_5=0` solves `C` rationally.  Substitution into
the other five equations gives primitive numerators in `Q[B,s_3,s_2]`
with signatures

```text
(total degree, deg_B, deg_s3, deg_s2, terms)

(10,5,10,5, 91),
(10,5,10,5,111),
(11,6,10,6,135),
(13,6,13,6,168),
(14,7,14,7,204).                                   (22)
```

Their only denominators are nonzero constants times `L_piv^2`.

On the complementary pivot wall,

```text
B=1771(s_2-1)/10260,                               (23)
```

the constant part of `E_5` is `1127` times the irreducible polynomial

```text
Phi(s_3,s_2)
 =-2057235s_2^2s_3+1925100s_2^2
  -873126s_2s_3^3-1345086s_2s_3^2
  +9986439s_2s_3-6914007s_2
  +672543s_3^5-4035258s_3^4
  +9167823s_3^3-5380344s_3^2
  -6161748s_3+5195399.                             (24)
```

Thus the mixed branch has an exact two-lane continuation:

```text
main lane:
  L_piv!=0, eliminate C and study five equations in B,s_3,s_2;

pivot lane:
  L_piv=0, Phi=0, then restore the remaining five equations.           (25)
```

Neither lane is proved empty here.

## 6. Exact controls and stopping boundary

The factorization companion exhaustively checks all `p^4` tuples
`(B,C,s_3,s_2)` over the good prime fields

```text
p=11,13,17,19,29.                                  (26)
```

There are no raw solutions for `p=11,13,17,29`.  At `p=19` there is one
raw solution, but it lies on the excluded localization divisor

```text
C K (54B+7)(13122C+1215B+91) Disc(H_2)=0.          (27)
```

The prime `23` is intentionally excluded because

```text
23 divides 73060029,
```

so the degree-ten leading coefficient drops and produces irrelevant
modular solutions.

These finite-field counts are exact hostile probes.  They do **not** imply
that the characteristic-zero ideal is empty: closed points can have
nontrivial residue-field degree, and bad reduction must be controlled.

The complete residual obligations are:

```text
K=0 singular-order total-branch boundary;
C=0 boundary;
L_piv=0 with Phi=0;
the localized five-equation main lane;
the separate H_4 S_4^2 stratum.                    (28)
```

Therefore this theorem proves no new degree-eighteen emptiness statement,
and it does not prove `JC(2)` or `DC(2)`.

## 7. Conditional quadratic-ring continuation

Equation (12) also gives the exact norm identity

```text
(7q_5+S_4 sqrt(H_2))(7q_5-S_4 sqrt(H_2))
 =-4(y-1)p_3^3.                                    (29)
```

After normalizing a squarefree quadratic `H_2` to `x^2-1`, its coordinate
ring is

```text
C[x,t]/(t^2-x^2+1)
  =C[s,s^(-1)],

s=x+t,                    conjugation s -> s^(-1). (30)
```

This suggests a much smaller `(1,3,3,3)` Laurent factorization.  It is not
yet a proved consequence: the current off-wall hypotheses do not exclude
additional nonnormal common roots of `p_3,q_5`, overlap above an `H_2`
branch root, repeated roots of `p_3`, or a branch-at-infinity degree drop.
Those gcd and endpoint sidecars must be proved before splitting the norm
factor atomwise.

## 8. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_h2_moving_root_reduction_thm2357.py
python3 -O 04-computation/jc2_degree18_h2_moving_root_reduction_thm2357.py

python3 04-computation/jc2_degree18_h2_factorization_pivot_thm2357.py
python3 -O 04-computation/jc2_degree18_h2_factorization_pivot_thm2357.py

python3 04-computation/jc2_degree18_h2_factorization_pivot_thm2357.py \
  --mode bruteforce --primes 11,13,17,19,29
```

Normal and optimized transcripts agree.  The stored primary transcript is

```text
05-knowledge/results/jc2_degree18_h2_factorization_pivot_thm2357.out
```

The moving-root and finite-field transcripts are

```text
05-knowledge/results/jc2_degree18_h2_moving_root_reduction_thm2357.out

05-knowledge/results/jc2_degree18_h2_factorization_modular_hostiles_thm2357.out
```

The companions reconstruct (5)--(14), every subresultant coefficient and
factor order, the recursion (16)--(19), pivot (20)--(24), all signatures,
generic and `K=0` hostile gcd controls, and every tuple in (26).  No
executable check uses Python `assert`.

Independent audit is pending. QED.
