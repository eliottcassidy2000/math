---
id: THM-2918
title: "Fixed-low-pair high-gap cube-root torsor and sign-blind quartic exit"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  Every fixed independent real mean-zero factorial low
  pair whose binary cubic has nonzero remainder modulo its positive
  quadratic Gram polynomial, together with a fixed mean-one subtraction,
  acquires exactly three locally unique real cubic-null branches after one
  sufficiently remote shared atom is added at the natural cube-root scale.
  The limiting branches form a free C3 torsor.
  All three have nonzero quartic remainder, certified by a positive
  quadratic field norm, while every real linear projection of their
  limiting quartic remainders sums to zero.  Arbitrary factorial
  three-slot low supports satisfy the hypothesis by THM-2824 and have
  exactly one positive local sheet; its scalar endpoint sign is not
  support-uniform.
source: root/fixed-low-triple-high-gap-trichotomy-2026-07-29
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2853-gamma-adjacent-tensor-cycle-weighted-positivity
related:
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
  - THM-2910-nonconsecutive-cubic-null-endpoint-holonomy-sign-reversal
  - THM-2914-eventual-high-gap-cubic-null-positive-holonomy-branch
script: 04-computation/gmc_fixed_low_pair_cube_root_torsor_thm2918.py
output: 05-knowledge/results/gmc_fixed_low_pair_cube_root_torsor_thm2918.out
script_sha256: 16d054113c0410d05eb86aa25f2e553a6017a59221c2ed259e877eebd6d79908
output_sha256: 4b4369ce5cc507261c99af7995afc03b75e524fa5ceb7afa207c422d75052df5
hash_basis: LF-normalized bytes
---

# THM-2918 -- fixed-low-pair high-gap cube-root torsor

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

The three limiting branches in THM-2914 are not an accident of the
consecutive low support.  They are the three cube roots of one nonzero
element in a quadratic quotient field.  This observation gives both a
strictly broader positive theorem and the exact reason that one scalar
endpoint orientation cannot control all branches.

## 1. Fixed low data and the remote shared atom

Let

```text
L(s^n)=n!,                         f_n=s^n/n!.          (1)
```

Fix real polynomials

```text
E_0,E_1,D in R[s],                 E_0,E_1 independent,

L(E_0)=L(E_1)=0,                   L(D)=1.              (2)
```

Put

```text
E(z)=E_0+zE_1,
q(z)=L(E(z)^2)=g_0+2g_1z+g_2z^2,
c_0(z)=L(E(z)^3).                                      (3)
```

The factorial functional is a positive inner product on real
polynomials:

```text
L(P^2)=integral_0^infinity P(s)^2 exp(-s) ds.           (4)
```

Consequently

```text
g_0>0,                 g_2>0,
Delta=g_0g_2-g_1^2>0.                                  (5)
```

In particular `q` is irreducible over `R`.  Let

```text
A=R[z]/(q),                    rbar=[c_0] in A,          (6)
```

and assume throughout the main theorem that

```text
rbar!=0.                                                (7)
```

For integers `R` larger than the fixed degrees, define

```text
B_R=f_R-D,
T_k(R)=(kR)!/(R!)^k,
epsilon_R=T_3(R)^(-1/3),
H_R=T_4(R)/T_3(R)^(4/3)>0.                              (8)
```

For real parameters `(xi,eta)`, write

```text
ell(z)=xi+eta z,
W_R(z)=E(z)+epsilon_R ell(z)B_R,                        (9)

Q_R(z)=L(W_R(z)^2),
C_R(z)=L(W_R(z)^3),
A_R(z)=L(W_R(z)^4).                                   (10)
```

Here `C_R` is a binary cubic moment, not the shared atom `B_R`.
The normalizations in `(2)` give

```text
L(B_R)=0,                         L(W_R(z))=0.           (11)
```

Thus a common projective zero of `Q_R` and `C_R` is genuinely a
first-three-moment null direction, not merely a quadratic/cubic
multipole coincidence.

## 2. Exact factorial tail

If `P` is any fixed polynomial and

```text
D^(k-j)=sum_a d_(k-j,a)s^a,
```

then direct expansion gives

```text
L(P B_R^k)
 =sum_(j=0)^k (-1)^(k-j) binom(k,j)
    sum_(a,d) d_(k-j,a)p_d
      (jR+a+d)!/(R!)^j.                               (12)
```

For

```text
M_(j,b)(R)=(jR+b)!/(R!)^j,                             (13)
```

every term in `(12)` belongs to a finite list, because `E_0,E_1,D`
are fixed.  After the scaling in `(8)`, the successive-ratio limit of
an order-three term

```text
T_3^(-k/3)M_(j,b),             0<=j<=k<=3,             (14)
```

where the error list has `1<=k<=3`, is `j^j/3^k`, with base `1`
for `j=0,1`.  The sole nondecaying order-three term is
`(k,j,b)=(3,3,0)`.  Every other limit is at most `4/9`.

For the fourth moment divided by `H_R`, lower terms have the forms

```text
T_3^((4-k)/3)M_(j,b)/T_4,       j<=k<=3,
M_(j,b)/T_4,                    j<4, k=4.              (15)
```

Their largest successive-ratio limit is `81/256`.
The finite ratio test therefore gives, uniformly with one parameter
derivative on every compact `(xi,eta)` set,

```text
Q_R -> q,
C_R -> c_0+ell^3,                                      (16)

A_R/H_R -> ell^4.                                      (17)
```

All coefficient errors in `(16)--(17)` are `O(2^(-R))` after
increasing the fixed threshold.  This is an exact factorial-ratio
argument, not an interchange of infinite series.

For large `R`, `Q_R` remains a positive quadratic Gram polynomial.
Let

```text
mathcal R_R(xi,eta)
 = coefficient vector of rem_(Q_R)(C_R) in the basis (1,z). (18)
```

Polynomial division is smooth while the leading quadratic coefficient
stays nonzero.  Equations `(16)` give

```text
mathcal R_R -> mathcal R_infinity
             = coefficient vector of rbar+[ell]^3 in A              (19)
```

in `C^1` on compact sets.

## 3. The exact `C_3` torsor

The algebra `A` is a copy of `C`.  More explicitly,

```text
i=(g_2z+g_1)/sqrt(Delta),                 i^2=-1,       (20)

omega=(-1+sqrt(3)i)/2,                    omega^3=1,
1+omega+omega^2=0.                                      (21)
```

Because `rbar!=0`, the equation

```text
[ell]^3=-rbar                                           (22)
```

has exactly three solutions in `A`.  If one is `ell_0`, all of them
are

```text
ell_j=omega^j ell_0,                    j=0,1,2.         (23)
```

Every element of `A` has unique real coordinates `xi+eta z`.
Thus `(22)` gives three distinct real points

```text
p_j=(xi_j,eta_j) in R^2.                               (24)
```

They form a free torsor under the group `mu_3={1,omega,omega^2}`.

The derivative of `(22)` at `ell_j` is multiplication by
`3ell_j^2`.  For `ell=a+bz`, define

```text
N_q(a,b)=g_2a^2-2g_1ab+g_0b^2.                        (25)
```

This is positive definite by `(5)`, and the real determinant of
multiplication by `3ell^2` is

```text
9(N_q(a,b)/g_2)^2>0.                                  (26)
```

Hence all three limiting zeros are simple.

Choose disjoint closed neighborhoods of the three points in `(24)`.
The ordinary implicit-function theorem, or the quantitative contraction
argument used in THM-2914, applied to the `C^1` convergence `(19)`
gives an integer `R_0` such that for every `R>=R_0`:

1. `mathcal R_R` has one and only one zero `p_(R,j)` in the `j`th
   neighborhood;
2. `p_(R,j)=p_j+O(2^(-R))`; and
3. these are three distinct real locally unique branches.

At a zero, `Q_R` divides `C_R` exactly.  There is no degeneracy hidden
in this statement: the two coefficient polynomials in `(9)` are
independent.  Indeed, cancellation of their `s^R` coefficients first
gives the same linear relation among the high coefficients, after
which the fixed independence of `E_0,E_1` forces both scalars to vanish.
Their quadratic Gram determinant is therefore positive.

## 4. Sign-blind quartic exit

On the `j`th finite branch put

```text
S_(R,j)=rem_(Q_R)(A_R/H_R).                            (27)
```

Equations `(17)` and the branch convergence give

```text
S_(R,j) -> [ell_j^4] in A.                             (28)
```

The limit is nonzero because `A` is a field and `ell_j!=0`.  Equivalently,
if `[ell_j^4]=a_j+b_jz`, then

```text
N_q(a_j,b_j)
 =N_q(xi_j,eta_j)^4/g_2^3>0.                         (29)
```

The corresponding finite quadratic quotient norm is positive definite
and converges to `(29)`.  Thus `S_(R,j)!=0` for every sufficiently
large `R`.  No one of the three cubic-null branches is quartic-null.

This is stronger than a sign choice for one endpoint coordinate: the
certificate retains both coefficients of the quartic remainder and
then takes their positive field norm.

## 5. Why every linear endpoint orientation loses a branch

The same torsor gives an exact stopping rule for scalar projections.
Since `omega^4=omega`,

```text
[ell_j^4]=omega^j[ell_0^4],

sum_(j=0)^2 [ell_j^4]=0 in A.                         (30)
```

Therefore every real linear functional `lambda:A->R` satisfies

```text
sum_(j=0)^2 lambda([ell_j^4])=0.                      (31)
```

THM-2872's limiting endpoint-holonomy determinant is such a linear
functional: it vanishes on every quartic multiple of `q`, and hence
factors through the quotient `A`.  Consequently, unless one limiting
endpoint value is zero, its three signs must be mixed.  THM-2910's
support-dependent sign reversal is a finite shadow of `(31)`, not a
failure of the full quartic remainder.

Equation `(31)` is a limiting identity.  No exact zero-sum assertion is
made for the three finite-`R` continuations.

There is also a precise `C_2/S_3` boundary.  Conjugation in `A` sends
`omega` to `omega^(-1)`, but it sends the fibre

```text
ell^3=-rbar
```

to the fibre `ell^3=-conj(rbar)`.  Unless `rbar` is conjugation-fixed,
or an extra phase identification is chosen, conjugation is not a
canonical involution of the same fixed torsor.  Thus `C_3` is the
canonical fixed-fibre symmetry.  The semidirect product
`C_3 semidirect C_2 isomorphic to S_3` belongs to the conjugate-family
groupoid, or to a noncanonically identified fibre, not automatically
to each low pair.

## 6. Every fixed factorial three-slot support

Fix arbitrary integers

```text
0<=a<b<c
```

and specialize

```text
E_0=f_b-f_a,
E_1=f_c-f_b,
D=f_c.                                                (32)
```

Then the finite branches live on the four slots

```text
{a,b,c,R}.                                            (33)
```

The low pair in `(32)` spans the mean-zero plane on the three fixed
slots.  THM-2853 expands all of its quadratic and cubic tensor entries
as strict positive sums of adjacent-difference tensors.  More
decisively, THM-2824 proves that no nonzero member of this plane can
have its first three factorial moments all zero.  Since a common
nonreal root of `q` and `c_0` would provide exactly such a member,

```text
[c_0]!=0 in R[z]/(q).                                 (34)
```

Thus every fixed triple `(a,b,c)` satisfies the hypothesis of the
theorem.  Adjoining the remote slot `R` produces three locally unique
real cubic-null branches at scale

```text
((3R)!/(R!)^3)^(-1/3),                                (35)
```

and the fourth factorial moment exits on every branch.

### 6.1. Exactly one positive sheet

The factorial specialization has more orientation than the abstract
theorem requires.  THM-2853 gives

```text
g_0>0,                  g_1>0,                 g_2>0,
rho=g_1/sqrt(g_0g_2) in (0,1).                         (36)
```

THM-2824 proves that its second divisibility invariant satisfies

```text
I_2(c_0)<0.                                             (37)
```

For a pure cube `ell^3`, with `xi,eta>0`, put

```text
sigma=xi sqrt(g_2)/(eta sqrt(g_0))>0.
```

Direct substitution into `I_2` gives

```text
I_2(ell^3)
 =g_0^(3/2)g_2^(1/2)eta^3
    (-sigma^3+3sigma-2rho).                            (38)
```

The cubic

```text
p_rho(sigma)=-sigma^3+3sigma-2rho
```

has one root in `(0,1)`, one in `(1,2)`, and is positive
exactly between those two positive roots.  At a limiting null,
`I_2(ell^3)=-I_2(c_0)>0`.

Geometrically, the inverse image of this strict half-plane under the
cube map on `A` has three connected angular sectors, permuted freely by
`mu_3`.  Equation `(38)` says that its intersection with
`xi>0,eta>0` is exactly one of those sectors.  Hence exactly one member
of the limiting torsor has both coordinates positive.  Its finite
continuation remains positive, and the other two local continuations
remain outside the positive quadrant, for every sufficiently large
`R`.  Indeed, on `eta=0` the required positive value of
`I_2(ell^3)` forces `xi<0`, while on `xi=0` it forces `eta<0`; no other
limiting sheet sits on the boundary of the positive quadrant.

For `(a,b,c)=(0,1,2)`, this positive member is the branch in THM-2914.

### 6.2. A positive-sheet scalar-orientation reversal

Even after selecting the unique positive sheet, its scalar endpoint
orientation depends on the fixed low triple.

For `(a,b,c)=(0,1,2)`, THM-2914 proves that the limiting endpoint
determinant is positive.  For `(a,b,c)=(0,8,10)`, exact arithmetic gives

```text
q=12869+61776z+110110z^2,

rbar=2(157110254569762+297335210491143z)/2695.         (39)
```

Write the positive cube root as `ell=eta(T+z)`.  Its projective
coordinate is the unique root in

```text
341/1000<T<683/2000                                   (40)
```

of

```text
425189351002334490T^3
-674002992104278980T^2
+229061153084068755T
-16579594004575622=0.                                 (41)
```

On `(40)`, the constant coefficient of
`rem_q((T+z)^3)` is negative.  Therefore
`eta^3=-rbar_0/a(T)>0`, so `(40)` is the positive sheet selected in
Section 6.1.  Its limiting endpoint determinant factors as

```text
J_(0,8,10)(T)
 =-572(110110T^2-12869)
      (5945940T^2-4954565T+694926)<0.                 (42)
```

Sturm isolation proves `(40)--(42)` exactly.  By continuation, the
finite positive sheet has negative endpoint determinant for every
sufficiently large remote slot.  Thus

```text
(0,1,2): positive J,             (0,8,10): negative J. (43)
```

The full quotient norm `(29)` stays positive in both cases.  It is
therefore genuinely load-bearing; it is not ornamental strengthening
of a scalar sign proof.

## 7. Boundary and scope

The hypotheses and losses are sharp:

1. if `E_0,E_1` are dependent, `q` does not define the quadratic field;
2. without `L(E_0)=L(E_1)=0` and `L(D)=1`, the same quotient algebra
   still gives quadratic/cubic multipole branches, but not the claimed
   first-three-moment null directions;
3. if `[c_0]=0`, equation `(22)` has only the triple root `ell=0`, the
   derivative vanishes, and another scale is required;
4. the theorem constructs only three local branches at the natural
   cube-root scale and does not exclude escaping or boundary branches;
5. it gives no explicit numerical `R_0`;
6. it does not classify complex coefficient branches away from these
   real continuations; and
7. it proves no arbitrary-support GMC or SFC theorem.

The preserved object is the full two-coordinate quotient remainder.
The destroyed object in a scalar endpoint proof is one transverse
coordinate of that remainder.  The positive norm `(25)` is the exact
sidecar that repairs this loss.

## 8. Exact companion

The exact companion verifies:

1. the generic quadratic-field relation and norm/determinant formulas;
2. the `mu_3` identities and the free three-sheet orbit;
3. multiplicative norm invariance and the zero-sum fourth-power orbit;
4. annihilation of arbitrary quartic multiples of `q` by the
   THM-2872 endpoint functional;
5. the `(0,1,2)` quotient remainder, ninth-degree cube-root eliminant
   and scalar zero-sum control;
6. all factorial rate bases in `(14)--(15)`;
7. all `84` triples `0<=a<b<c<=8` as a finite hostile control of
   `(34)`, including strict quadratic/cubic tensor positivity; and
8. the exact positive-sheet bracket and negative endpoint factor for
   `(0,8,10)`, paired with the positive `(0,1,2)` endpoint.

It uses explicit `require` gates in ordinary and optimized mode.  Run

```text
python3 04-computation/gmc_fixed_low_pair_cube_root_torsor_thm2918.py
python3 -O 04-computation/gmc_fixed_low_pair_cube_root_torsor_thm2918.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_fixed_low_pair_cube_root_torsor_thm2918.out.
```

The companion checks the exact algebraic core.  The uniform tail and
finite-`R` continuation are proved in Sections 2--4; no numerical cutoff
is inferred from the finite census.
