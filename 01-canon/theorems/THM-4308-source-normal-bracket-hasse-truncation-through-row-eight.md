---
id: THM-4308
title: "Source-normal bracket and depth-filter truncation through row eight"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the fixed normalized
  b=d=0 reduced (2,3) source gauge with residual weight at most twelve, the
  degree-capped Jacobian and nodal-defect equations through G-row eight have
  exactly four successive scalar bracket compatibilities E5--E8. Imposing
  the exact row-eight projections of the depth modules P2 and P3 is
  equivalent to three further parameter equations Delta=896/15,
  Theta=512/75, and zeta3=-3Phi/2. When these seven equations hold, rows four
  through seven are unique and the terminal row-eight solution space is
  affine of dimension seven. The weight-twelve response (U,W,Z) is an
  explicit three-parameter function of (Phi,eta,xi10). An exact finite
  gate-interior control proves that these truncated equations force no
  coefficient wall. This is not a full B2 lift, polynomial-termination,
  seam-entry, exact-M12-existence, JC(2), or DC(2) theorem.
source: root / planar-Jacobian bracket-Hasse continuation, 2026-09-01
audit: >
  PASS. The primary SymPy certificate expands the source directly, derives
  the bracket quotient and its four one-dimensional cokernels, constructs
  all rows, enumerates the two exact depth-projection matrices, evaluates
  all 21 left-null relations, rank-checks the terminal affine fibre, and
  replays a gate-interior finite jet coefficient by coefficient. A separate
  standard-library Fraction/sparse-polynomial verifier reconstructs the
  equations, ranks, response, positive control, and perturbation hostiles
  without importing the primary code. Normal and optimized streams match
  their frozen transcripts exactly.
depends_on:
  - THM-3989-cusp-log-laurent-conductor-and-nondividing-depth-reduction
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4298-weighted-face-source-normal-unimodular-visibility-transform
related:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4005-reduced-two-three-live-seam-invariant-support-atlas
primary_script: 04-computation/jc2_source_normal_bracket_hasse_rows8_thm4308.py
primary_output: 05-knowledge/results/jc2_source_normal_bracket_hasse_rows8_thm4308.out
primary_script_sha256: 3703758f87a628583cf0f2f9e8fb1973f8ee65c875ab679acc20a9867c27e7f1
primary_output_sha256: a32ec035ce23ae496a027b657bbe63596e023f5a0e9dca33f5a220b800485cfc
independent_audit_script: 04-computation/jc2_source_normal_bracket_hasse_rows8_thm4308_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_bracket_hasse_rows8_thm4308_independent_audit.out
independent_audit_script_sha256: 61e634be42bb45476bb25a6e778867ec7c2fba5d90cbd9769f53c5f768ab9a5a
independent_audit_output_sha256: 100f9b7526a00684a565b88ad9b43c4a44186787b0b13b36f0ec53dfa084ae83
hash_basis: raw LF bytes
---

# THM-4308 -- source-normal bracket and depth-filter truncation through row eight

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. THIS IS AN EXACT FINITE
ROW-EIGHT PROJECTION THEOREM. IT DOES NOT ASSERT AN ALL-ROW `B_2` LIFT,
POLYNOMIAL TERMINATION, SEAM ENTRY, OR A CONSEQUENCE FOR `JC(2)` OR `DC(2)`.**

Work over an algebraically closed field `k` of characteristic zero. Use the
fixed `a=1`, `gamma=-1/2`, `I=3/4`, `b=d=0` normalization inherited from
THM-4007 and THM-4230:

```text
p=t(1+x^2t),                 y=xtp,                 u=x^2t,
P(A,C)=C^2-A^3+(3/4)A+1/4,
P(A,C)=G=-u/2+H(p,y).                                      (1)
```

The full residual source of weight at most twelve in this gauge is

```text
H=-3p+(8/3)p^2-(1376/135)p^3
  +K y^2+Phi p^2y+Delta p^4+Theta p y^2
  +eta p^3y+zeta3 y^3+upsilon5 p^5+xi10 p^2y^2
  +alpha11 p^4y+beta11 p y^3
  +U p^6+W p^3y^2+Z y^4,

K=2848/45-(7/6)Delta.                                    (2)
```

No monomial of residual weight greater than twelve is allowed in `(2)`.
That upper bound is part of the finite universe, not a conclusion of this
theorem.

Write

```text
A=sum_(n=0)^8 A_n(x)t^n,             C=sum_(n=0)^8 C_n(x)t^n. (3)
```

The inherited rows are

```text
A_0=1+x^2/4,
C_0=-3x/4-x^3/8,

A_1=4/3+2x^2,
C_1=-4x-(3/2)x^3,

A_2=-32/9-(4/5)x^2,
C_2=(88/15)x-(12/5)x^3,

A_3=2176/135-(Phi/2)x
    +(1088/315-(4/7)K)x^2-(32/15)x^4,

C_3=3Phi/4+(-8128/315+(6/7)K)x+(3Phi/8)x^2
    +(736/105+(3/7)K)x^3+(8/5)x^5.                       (4)
```

The extreme Laurent rows give the sharp necessary caps

```text
deg A_n<=n+1,                   deg C_n<=n+2,       4<=n<=8. (5)
```

As in the earlier source-normal theorems, `(5)` is a consequence of depth
membership. It is not by itself a converse to membership in `B_2`.

## 1. The actual source rows and the visibility firewall

Directly expanding `(2)` under `(1)` gives

```text
G_4=Delta+Phi x+(K-1376/45)x^2+(8/3)x^4,

G_5=upsilon5+eta x+(4Delta+Theta)x^2+3Phi x^3
    +(2K-1376/45)x^4,

G_6=U+alpha11 x+(5upsilon5+xi10)x^2+(4eta+zeta3)x^3
    +(6Delta+3Theta)x^4+3Phi x^5+(K-1376/135)x^6,

G_7=(6U+W)x^2+(5alpha11+beta11)x^3
    +(10upsilon5+4xi10)x^4+(6eta+3zeta3)x^5
    +(4Delta+3Theta)x^6+Phi x^7,

G_8=(15U+5W+Z)x^4+(10alpha11+4beta11)x^5
    +(10upsilon5+6xi10)x^6+(4eta+3zeta3)x^7
    +(Delta+Theta)x^8.                                  (6)
```

Thus THM-4298's weight-twelve flag is indeed

```text
(h_0,h_1,h_2)=(U,6U+W,15U+5W+Z).                        (7)
```

But `(6)` is the firewall between **support visibility** and **row
existence**: rows six through eight contain lower-weight echoes as well as
the three weight-twelve channels. THM-4298 reconstructs the labelled face
from a row flag; it does not supply any of the rows `(3)`.

## 2. The bracket theorem through row eight

For `m>=2`, put `v_n=(A_n,C_n)` and define the lower-row polynomials

```text
B_m=sum_(i=1)^(m-1)
       ((m-i)A_i'C_(m-i)-i A_i C_(m-i)'),

T_m=sum_(i=1)^(m-1) C_i C_(m-i)
    -sum_(i+j+k=m; 0<=i,j,k<m) A_i A_j A_k.              (8)
```

Then the coefficient of `t^(m-1)` in the Keller equation is

```text
m det(v_0',v_m)+B_m=0.                                  (9)
```

Along the boundary row,

```text
grad P(A_0,C_0)=q(-C_0',A_0'),             q=-(x^2+6)/2. (10)
```

Eliminating `v_m` between `(9)` and `[t^m]P=G_m` gives the exact bracket
compatibility

```text
G_m=T_m-(q/m)B_m.                                       (11)
```

Here `A_0'=x/2` and `C_0'=-3(x^2+2)/8`, so the row operator divided by
`m` is

```text
(A_m,C_m) -> (x/2)C_m+(3/8)(x^2+2)A_m.                  (11a)
```

It is onto the polynomials of degree at most `m+3`: first choose the
constant term of `A_m` to match the value at `x=0`, after which the
remainder is divisible by `x` and is absorbed by `C_m`. The caps in `(5)`
are preserved. Coprimality of `x` and `x^2+2` also shows that its
homogeneous kernel is exactly

```text
v_n -> v_n+theta_n v_0',                   deg theta_n<=n. (12)
```

This tangent change fixes `G_n` and changes the next bracket by

```text
delta G_(n+1)=q' theta_n-(q/(n+1))theta_n'.              (13)
```

For `m=n+1`, the map in `(13)` from `k[x]_(<=m-1)` to
`k[x]_(<=m)` is injective with a one-dimensional cokernel. Indeed, a zero
would force `theta_n` to be proportional to `q^m`, which has degree `2m`,
outside the domain unless `theta_n=0`. Integer generators of the four
cokernels, applied to `f=sum f_jx^j`, are

```text
m=5:  21f_0+14f_2+36f_4,
m=6:  77f_0+42f_2+84f_4+360f_6,
m=7:  143f_0+66f_2+108f_4+360f_6,
m=8:  715f_0+286f_2+396f_4+1080f_6+5040f_8.             (14)
```

Applying `(14)` successively to `(6)` gives

```text
E5: 2025upsilon5+9000Delta+1350Theta+184832=0,

E6: 200475U+109350xi10-5593860Delta-529200Theta
    -137763328=0,

E7: 801900W+1782000Delta^2+156163200Delta+868725Phi^2
    +27390480Theta-3434400xi10+12891824128=0,

E8: 21651300Z-225022050Delta^2-59073300Delta Theta
    -9512522400Delta+34749000Phi^2+39092625Phi eta
    +940522560Theta+185376600xi10-1112446017536=0.       (15)
```

> **Bracket conclusion.** Subject to `(1)--(6)`, degree-capped rows
> `A_4,...,A_8,C_4,...,C_8` satisfying
>
> ```text
> [t^r]J_(x,t)(A,C)=delta_(r,0),             0<=r<=7,
> [t^r](P(A,C)-G)=0,                         0<=r<=8       (16)
> ```
>
> exist if and only if the four equations `(15)` hold. When they do,
> `v_4,...,v_7` are unique. The terminal `v_8` is an affine translate of
> the nine-dimensional kernel `(12)`.

This is a triangular statement: `E5` is the obstruction to fixing `v_4`
with `G_5`, then `E6` fixes `v_5`, `E7` fixes `v_6`, and `E8` fixes `v_7`.
The terminal tangent is not tested by a nonexistent `G_9` row.

## 3. Exact row-eight depth projections

Retain THM-3989's filtration

```text
P_d=B_2 intersect tau^(-d)k[s,tau].                     (17)
```

Its symbol theorem says `P_d/P_(d-1)=s^d k[s]`. The representatives
`x^d`, `x^(d-1)u`, and `x^d p^c y^e` for `2c+3e>=2` kill every possible
principal symbol. Induction on `d` therefore gives the exact spanning law

```text
P_d=span_k{x^a u^b p^c y^e: a,b,c,e>=0, a+b<=d}.         (18)
```

For each generator in `(18)`, source normalization gives

```text
x^a u^b p^c y^e
 =x^(a+2b+e)t^(b+c+2e)(1+x^2t)^(c+e).                  (19)
```

Let `pi_8` retain rows `0<=n<=8`. The exact finite matrices obtained from
`(19)` have the following declared universes:

```text
pi_8(P_2): coordinates (n,r), 0<=n<=8, 0<=r<=n+2;
             63 ambient rows, 131 generator columns, rank 51,
             left nullity 12;

pi_8(P_3): coordinates (n,r), 0<=n<=8, 0<=r<=n+3;
             72 ambient rows, 204 generator columns, rank 63,
             left nullity 9.                            (20)
```

The columns are precisely the quadruples in `(18)` with
`b+c+2e<=8`. Thus `(20)` is an exact rational rank calculation in an
explicit finite universe, not a sampled Hasse heuristic.

Write

```text
a_(n,r)=[x^r]A_n,                    c_(n,r)=[x^r]C_n.   (21)
```

Six useful members of the 21-dimensional left-null bank evaluate on the
bracket-compatible family as

```text
-4a_(2,0)+3a_(3,2)-2a_(4,4)+a_(5,6)
 =-(15Delta-896)/45,                                      (22)

-10a_(4,1)+6a_(5,3)-3a_(6,5)+a_(7,7)
 =3(3Phi+2zeta3)/20,                                      (23)

15a_(4,0)-10a_(5,2)+6a_(6,4)-3a_(7,6)+a_(8,8)
 =4(13215Delta+7950Theta-2475c_(8,9)+6075xi10
    -2583808)/7425,                                       (24)

15c_(4,1)-10c_(5,3)+6c_(6,5)-3c_(7,7)+c_(8,9)
 =(6900Delta-13425Theta+4950c_(8,9)-12150xi10
   +3159808)/4950,                                        (25)

-15a_(4,1)+8a_(5,3)-3a_(6,5)+a_(8,9)
 =(27Phi-40c_(8,10)+27eta+45zeta3)/30,                   (26)

-6c_(3,0)+5c_(4,2)-4c_(5,4)+3c_(6,6)-2c_(7,8)+c_(8,10)
 =(40c_(8,10)-27eta-27zeta3)/40.                         (27)
```

Equations `(22)--(27)` force

```text
Delta=896/15,                 Theta=512/75,
zeta3=-3Phi/2,                                           (28)

c_(8,9)=(1215xi10-348032)/495,
c_(8,10)=27(eta+zeta3)/40.                              (29)
```

For example, eliminating `c_(8,9)` from `(24)--(25)` gives

```text
33330Delta+2475Theta-2007808=0,                         (30)
```

which yields the stated `Theta` after `(22)`.

The exact audit evaluates all twelve `P_2` and all nine `P_3` left-null
relations, not only the six displayed witnesses. After `(28)`, their
coefficient matrix on the nine terminal tangent coordinates has rank two;
`(29)` solves those two independent rows, and every remaining null relation
then vanishes identically. Therefore:

> **Depth-filtered finite iff.** Rows satisfying `(16)` and also
>
> ```text
> pi_8(A) in pi_8(P_2),                 pi_8(C) in pi_8(P_3) (31)
> ```
>
> exist if and only if `(15)` and `(28)` hold. Their terminal solution
> space is affine of dimension `9-2=7`. Equations `(29)` fix two terminal
> coordinates; the other seven tangent coordinates remain free.

The assertion in `(31)` is membership in an exact **projection** of the
depth modules. It is not membership of an infinite formal series, much less
the existence of polynomials in `B_2` satisfying every later row.

## 4. The three-parameter weight-twelve response

Substituting `(28)` into `(2)` and `(15)` gives

```text
K=-32/5,                 upsilon5=-731648/2025,
zeta3=-3Phi/2,                                             (32)

U=(475515904-109350xi10)/200475,

W=-(4343625Phi^2-17172000xi10+143826305024)/4009500,

Z=(12506118074368-173745000Phi^2-195463125Phi eta
   -926883000xi10)/108256500.                              (33)
```

Thus the top face factors through the three parameters
`(Phi,eta,xi10)`. The coefficients `alpha11,beta11` remain live in the
lower-weight rows but do not enter `(33)`. This is a response calculation,
not a proof that any point of the response image extends to an exact-M12
Keller pair.

## 5. Exact positive and hostile controls

For a full finite positive control, take

```text
Phi=eta=xi10=zeta3=alpha11=beta11=0,
Delta=896/15,                    Theta=512/75.             (34)
```

Equations `(32)--(33)` then give

```text
U=475515904/200475,
W=-35956576256/1002375,
Z=3126529518592/27064125,                                  (35)

Lambda=U+W+Z=443979321344/5412825,
D=W^2-4UZ
 =5173344945126466650112/27128402296875.                   (36)
```

Choose the terminal tangent so that

```text
c_(8,9)=-348032/495,                  c_(8,10)=0,          (37)
```

and set the other seven terminal tangent coordinates to zero in the
certificate's fixed basis. Direct expansion then verifies every equation in
`(16)` and both column-space memberships in `(31)`. All four quantities
`U,Z,D,Lambda` are nonzero. Hence this exact finite bracket/Hasse jet lies in
the interior of THM-4298's four coefficient walls. The row-eight theorem
cannot itself force any one of those walls.

The hostile bank is coordinatewise and exact:

1. after each compatible predecessor, increasing respectively
   `upsilon5,U,W,Z` by one changes the left sides of `E5,E6,E7,E8` by
   `2025,200475,801900,21651300`; the corresponding tangent system becomes
   inconsistent;
2. increasing `Delta` from `896/15` by one makes `(22)` equal `-1/3`;
3. increasing `zeta3` from `-3Phi/2` by one makes `(23)` equal `3/10`;
4. after fixing `Delta`, increasing `Theta` from `512/75` by one changes
   the eliminant `(30)` by `2475`;
5. the THM-4298 observer packets `(1,6,15)` for `p^6` and `(0,1,5)` for
   `p^3y^2` remain distinct, but they test visibility only.

For items 2--4, the bracket equations `E5--E8` are re-solved after the
perturbation. The direct Jacobian and `P-G` rows still pass, so these are
hostiles to depth-projection membership rather than accidental bracket
failures.

## 6. Scope and reproduction

The theorem proves exactly the row-eight projected iff `(31)`. It does not
prove any of the following:

1. that a finite solution extends to row nine or to all rows;
2. that its coordinates come from actual elements of `P_2,P_3` rather than
   their row-eight projections;
3. polynomial termination or the existence of a Keller pair in `B_2`;
4. entry into the reduced seam by an arbitrary hypothetical counterexample;
5. that residual weight is at most twelve outside the declared universe;
6. any coefficient wall, exact-M12 boundary case, or emptiness result;
7. `JC(2)` or `DC(2)`.

Reproduce from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_bracket_hasse_rows8_thm4308.py
python3 -B -O 04-computation/jc2_source_normal_bracket_hasse_rows8_thm4308.py
python3 -B 04-computation/jc2_source_normal_bracket_hasse_rows8_thm4308_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_bracket_hasse_rows8_thm4308_independent_audit.py
```

The normal and optimized streams must byte-match their corresponding frozen
outputs. The two implementations share only the theorem's displayed gauge
and equations; the independent path imports no code from the primary path.
