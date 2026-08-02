---
id: THM-3110
title: "Arbitrary-anchored product-Gamma dominant tail and low-histogram reduction"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  On every
  anchored support {0,a,b}, 1<=a<b, the two cleared width-three invariants
  have universal 24- and 25-row response banks.  After an exact
  chamber-dependent common-root deletion, one positive Ferrers row
  majorizes every negative row and gives an explicit all-layer tail
  threshold.  Dual Cauchy and exact chamber-Newton certificates prove every
  Schur response coefficient through total degree twelve, and a separate
  two-active-layer certificate reaches total degree fourteen.  The residual
  finite histogram bank remains OPEN.  Its labelled Ewens zeta transform is
  a signed current supported exactly in atomic rank four, reducing the
  missing argument to a relative rank-four cycle/edge-slide problem.
source: root/multiscale-newton-flag-2026-08-02
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3107-finite-layer-product-gamma-consecutive-width-three-orientation
related:
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3100-product-gamma-response-monge-compactification-and-bad-prefix-spectrum
script: 04-computation/gmc_product_gamma_arbitrary_anchored_tail_thm3110.py
output: 05-knowledge/results/gmc_product_gamma_arbitrary_anchored_tail_thm3110.out
script_sha256: 4bc0a6e39ac4034f5f702e4eb4ffca20349ebed05e57a160dba58524318f9244
output_sha256: b1c2927ce73d866f1b52faf8d916dda4c80f2bf4f1da0983b71f1ccf12b5567b
secondary_script: 04-computation/gmc_product_gamma_arbitrary_anchored_schur_thm3110.py
secondary_output: 05-knowledge/results/gmc_product_gamma_arbitrary_anchored_schur_thm3110.out
secondary_script_sha256: 15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f
secondary_output_sha256: 5f85b5fc3ace0d2be50e9278a8c64a3543869b37c5ea16b6b5a166c219956882
hash_basis: LF-normalized bytes
---

# THM-3110 -- arbitrary-anchored product-Gamma dominant tail and low-histogram reduction

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3107 proves every anchored arithmetic progression `{0,d,2d}` for an
arbitrary finite product of positive Gamma shapes.  The same response
method has a uniform, but not yet complete, extension to every anchored
triple.  The large-layer tail is controlled by one Ferrers-dominant row.
The remaining obstruction is finite and begins only beyond a substantial
Newton-positive front face.

This theorem records that reduction exactly.  It does **not** claim that
every anchored triple is good.

## 1. The two universal response banks

Fix integers

```text
1<=a<b                                                        (1)
```

and a positive moment sequence `(w_n)`, with `w_0=1`.  Put

```text
f_n=s^n/w_n,             U=f_a-f_0,             V=f_b-f_a.    (2)
```

Define `g_ij,t_ijk,I_1,I_2` exactly as in THM-3107:

```text
g11=L(U^2),  g12=L(UV),  g22=L(V^2),
t111=L(U^3), t112=L(U^2V), t122=L(UV^2), t222=L(V^3),

I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22,
I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2.          (3)
```

Use the shorthand

```text
A=w_a, B=w_b, A2=w_(2a), AB=w_(a+b), B2=w_(2b),
A3=w_(3a), A2B=w_(2a+b), AB2=w_(a+2b), B3=w_(3b).   (4)
```

After clearing the positive denominators, set

```text
B1=-A^5 B^3 I1,                  B2=-A^5 B^4 I2.     (5)
```

Direct multilinear expansion gives the following complete banks:

```text
B1=
 +A^5 B3 -3A^4 B AB2 -2A^4 AB B2
 +3A^3 B^2 A2B +2A^3 B A2 B2 +4A^3 B AB^2
 -2A^3 A2 B3 +3A^3 B2 A2B -A^2 B^3 A3
 -6A^2 B^2 A2 AB +6A^2 B A2 AB2 -6A^2 B AB A2B
 -3A^2 B B2 A3 +2A B^3 A2^2 -3A B^2 A2 A2B
 +6A B^2 AB A3 +A A2^2 B3 -3A A2 B2 A2B
 +2A AB B2 A3 -B^3 A2 A3 -3B A2^2 AB2
 +6B A2 AB A2B +B A2 B2 A3 -4B AB^2 A3;            (6)

B2=
 -A^5 B2^2 +4A^4 B AB B2 -2A^4 AB B3 +3A^4 B2 AB2
 -2A^3 B^2 A2 B2 -4A^3 B^2 AB^2 +2A^3 B A2 B3
 -6A^3 B B2 A2B +4A^2 B^3 A2 AB -3A^2 B^2 A2 AB2
 +6A^2 B^2 AB A2B +3A^2 B^2 B2 A3 +2A^2 A2 AB B3
 -3A^2 A2 B2 AB2 +A^2 B2^2 A3 -A B^4 A2^2
 -4A B^3 AB A3 -2A B A2^2 B3 +6A B A2 B2 A2B
 -4A B AB B2 A3 +B^4 A2 A3 +3B^2 A2^2 AB2
 -6B^2 A2 AB A2B -B^2 A2 B2 A3 +4B^2 AB^2 A3.    (7)
```

Thus the generic banks contain `24` and `25` collected rows.  When two
linear subscripts coincide, as at `b=2a`, substitution in `(6)--(7)`
combines the rows; no genericity assumption is being made.

## 2. Product-Gamma roots and the chamber divisor

Now specialize to a finite nonempty product of positive Gamma shapes:

```text
w_n=product_(nu=1)^J (theta_nu)_n,       theta_nu>0.  (8)
```

Put `x_nu=1/theta_nu` and remove the common geometric factor.  For one
Gamma layer define

```text
P_n(x)=product_(r=1)^(n-1)(1+r x),
P_lambda(x)=product_i P_(lambda_i)(x).                 (9)
```

Every row of `(6)` has weighted sum `5a+3b`; every row of `(7)` has
weighted sum `5a+4b`.  Hence the geometric factors cancel.  Put

```text
M=max(2a,b),
C1=(a,a,a,M),                 C2=(a,a,a,M,M).         (10)
```

The Ferrers-column minimum over every row is exactly

```text
gcd_lambda P_lambda=P_(C1)       in B1,
gcd_lambda P_lambda=P_(C2)       in B2.               (11)
```

Here is a symbolic proof of exactness.  Below root level `a`, every row has
at least four, respectively five, parts, and `Q_1,Q_2` attain those counts.
For `a<=r<M`, every `B1` row has at least one part of size at least `M`
and every `B2` row has at least two; the rows `A^5B3` and `A^5B2^2`
attain one and two.  At and above `M`, the rows `AB^3A2^2` and
`AB^4A2^2` have no part larger than `M`.  These three root ranges give
exactly the column multiplicities of `(10)`.  This proof covers the wall
`b=2a` as well as both open chambers.

Thus `b>=2a` deletes `P_a^3P_b` and `P_a^3P_b^2`, while `b<2a`
deletes `P_a^3P_(2a)` and `P_a^3P_(2a)^2`.  For a bank row `lambda`, write

```text
R_(j,lambda)(x)=P_lambda(x)/P_(Cj)(x).                 (12)
```

After the positive deletion, the cleared invariant is the signed product
kernel

```text
K_j(x_1,...,x_J)
 =sum_(lambda in Bj) c_lambda product_nu R_(j,lambda)(x_nu).  (13)
```

Consequently `K_j>0` implies `I_j<0`.

## 3. Ferrers dominance and all-layer tails

For a partition `lambda`, let its root multiset be

```text
E(lambda)=disjoint_union_i {1,...,lambda_i-1}.         (14)
```

The unique positive rows which majorize every negative row are

```text
Q1=(3b,2a,2a,a),               c_(Q1)=1,
Q2=(3b,a+b,2a,a,a),            c_(Q2)=2.              (15)
```

The assertion follows by sorting the nine possible linear subscripts in
the four exact ratio intervals

```text
1<b/a<3/2,  3/2<=b/a<2,  2<=b/a<3,  3<=b/a,          (16)
```

and comparing every prefix sum.  All differences are linear in `a,b` and
are nonnegative at the two endpoints of the relevant interval.  Equality
does not identify a second dominant response row.

Majorization implies that, after subtracting the common multiset
`E(C_j)`, the decreasing list of roots of every negative row is
coordinatewise at most the corresponding list for `Q_j`, after padding by
zeros.  If positive vectors `q,s` satisfy `s_i<=q_i`, then successive
one-coordinate decreases give

```text
e_k(s)/e_k(q)<=e_1(s)/e_1(q),              k>=1.      (17)
```

Indeed, when only coordinate `i` decreases by `delta`, the relative loss
in `e_k` is

```text
delta e_(k-1)(q without q_i)/e_k(q),                   (18)
```

which is at least `delta/e_1(q)` because
`e_1(q without q_i)e_(k-1)(q without q_i)>=e_k(q without q_i)`.
Multiplying the successive relative inequalities proves `(17)`.

Let

```text
T(n)=binom(n,2),
D1=T(3b)+2T(2a)-2T(a)-T(M),
D2=T(3b)+T(a+b)+T(2a)-T(a)-2T(M).                    (19)
```

These are the degree-one coefficients of the two dominant residual rows.
For a negative row `lambda`, put

```text
Delta_(j,lambda)=T(Q_j)-T(lambda),
rho_(j,lambda)=1-Delta_(j,lambda)/D_j,                 (20)
```

where `T(lambda)=sum_i T(lambda_i)`.  Then `(17)` proves that the exact
tail threshold is

```text
L_j(a,b)=min { ell>=0:
  sum_(c_lambda<0) |c_lambda| rho_(j,lambda)^ell<c_(Qj) }.  (21)
```

Every present response coefficient whose histogram has at least `L_j(a,b)`
active Gamma layers is strictly positive, with every per-layer degree lying
within the dominant residual degree.  Degrees above that bound vanish in
every row.  This is an exact finite formula, not an asymptotic statement.

For reference,

```text
D1=3a^2-a+4b^2-b                    if b>=2a,
D1=(2a^2+9b^2-3b)/2                 if b<2a,
D2=2a^2+ab-a+4b^2-b                 if b>=2a,
D2=-2a^2+ab+a+5b^2-2b               if b<2a.          (22)
```

The negative coefficient masses are `37` and `39`, and the closest
negative rows have gaps

```text
delta1=min(a^2,2b(b-a)),
delta2=a min(a,b-a).                                     (23)
```

Thus a convenient closed upper bound follows by replacing every ratio in
`(21)` by `1-delta_j/D_j`:

```text
ell>log(37)/[-log(1-delta1/D1)]       for B1,
ell>log(39/2)/[-log(1-delta2/D2)]     for B2.          (24)
```

At `(a,b)=(1,2)`, the exact sums in `(21)` give the sharp dominant-mode
certificate thresholds `14` and `16` from THM-3107.

## 4. Collision divisor and the Schur-positive initial face

Write

```text
Phi_j(e_(k1)...e_(kell))
 =sum_lambda c_lambda product_(r=1)^ell
    e_(kr)(E(lambda) minus E(C_j)).                    (25)
```

This is exactly the coefficient of a monomial using `ell` active Gamma
layers with per-layer degrees `k_1,...,k_ell`.  Every coefficient of total
degree below five vanishes.

The polynomial dependence on the support indices carries a universal
degree-nine divisor.  Work in the formal response ring `Q[a,b][[X,Y]]` and
introduce the continuous-index moment response

```text
F_t(X)=exp(sum_(r>=1)(-1)^(r-1) p_r(t)p_r(X)/r),
p_r(t)=sum_(j=1)^(t-1)j^r.                            (26)
```

Faulhaber's formula makes every coefficient of `F_t` a polynomial in `t`
and gives `F_0=F_1=1`.  Adjoin the group-like formal exponentials
`y_t=exp(tY)`, define the moment functional by

```text
L_X(y_t)=F_t(X),                 f_t=y_t/F_t(X).       (27)
```

At positive integer `t`, `(26)` is exactly
`product_(x in X)P_t(x)`, so `(27)` is the normalized monomial of the
actual response model.  Since `F_t` has constant coefficient one, it is a
unit.  Formal Taylor divisibility, applied to both `y_t` and `F_t^(-1)`,
therefore gives

```text
f_(t+h)-f_t in (h).                                   (28)
```

Thus the collision differences below are coefficientwise divisibilities
in the normalized response algebra, not a heuristic interpolation.

Each term of `I1` is multihomogeneous of degree four in `U` and degree
three in `V`; each term of `I2` has degrees three and four.  Therefore

```text
U=f_a-f_0       supplies a^4 in I1 and a^3 in I2,
V=f_b-f_a       supplies (b-a)^3 in I1 and (b-a)^4 in I2.  (29)
```

There are also two factors at the third collision.  Put `W=U+V=f_b-f_0`
and replace `V` by `epsilon W-U`.  With the evident polarized notation,
direct multilinearity gives

```text
I1(U,epsilon W-U)=epsilon^2 A+epsilon^3 B,
I2(U,epsilon W-U)=-epsilon^2 A+epsilon^3 C+epsilon^4 D, (30)

A=3gUU^2 tUWW-6gUU gUW tUUW-gUU gWW tUUU+4gUW^2 tUUU,
B=-gUU^2 tWWW+3gUU gWW tUUW-2gUW gWW tUUU,
C=2gUU^2 tWWW-6gUU gWW tUUW+4gUW gWW tUUU,
D=-2gUU gUW tWWW+3gUU gWW tUWW-gWW^2 tUUU.
```

Here `W=f_b-f_0` lies in `(b)`, so `(30)` supplies `b^2` in both
invariants.  Multiplying `I1,I2` by their clearing factors
`F_a^5F_b^3,F_a^5F_b^4` only multiplies by units and preserves all three
ideals.  The polynomials `a,b,b-a` are pairwise coprime.  Dividing by the
common response polynomial in `(11)` likewise preserves these coefficient
ideals because that polynomial has constant coefficient one: formal series
division is triangular.  Consequently

```text
base1=a^4b^2(b-a)^3,             base2=a^3b^2(b-a)^4  (31)
```

divide every response multicoefficient, and a coefficient of total response
degree `N`, divided by `base_j`, is a chamber polynomial of total degree at
most `2N-9`.  This is the load-bearing fact which makes the finite Newton
tables global identities rather than fitted grids.

There is a stronger symmetric-function organization.  Let `X` be the
alphabet of Gamma-layer variables and let `S_R=E(R)\E(C_j)` be a residual
root alphabet.  Dual Cauchy gives

```text
product_(x in X) product_(r in S_R)(1+rx)
 =sum_mu s_mu(X)s_(mu')(S_R).                          (32)
```

Therefore

```text
K_j(X)=sum_mu A_(j,mu)s_mu(X),
A_(j,mu)=Phi_j(s_(mu')).                               (33)
```

Schur functions have nonnegative monomial coefficients.  It is therefore
enough, degree by degree, to orient the finite family `A_(j,mu)` rather than
every layer histogram separately.

For one active layer, the degree-five coefficients factor as

```text
Phi_1(e5)=a^4 b^2 (b-a)^3(3a+5b-5)/2,
Phi_2(e5)=a^3 b^2 (b-a)^4(3a+4b-5)/2.                 (34)
```

More generally, for every partition `mu` of five, let `f^mu` be the number
of standard tableaux and put

```text
kappa_mu=sum_i mu_i(mu_i-2i+1).                        (35)
```

Exact character algebra gives

```text
A_(1,mu)=base1 f^mu/2 (3a+5b-kappa_mu/4),
A_(2,mu)=base2 f^mu/2 (3a+4b-kappa_mu/4).              (36)
```

Since `kappa_mu<=5*4=20`, every degree-five Schur coefficient is strictly
positive.

For degrees six through twelve use the nonnegative integer chamber
coordinates

```text
wide:   a=A+1,       b=2a+D,
tight:  a=A+D+2,     b=A+2D+3.                         (37)
```

For every partition `mu` of size `N`, `5<=N<=12`, the exact binomial
Newton expansion

```text
A_(j,mu)/base_j
 =sum_(r+s<=2N-9) c_(r,s) binom(A,r)binom(D,s)         (38)
```

has `c_(r,s)>=0` in both chambers.  The complete exact census is

```text
N       5     6      7      8      9      10      11      12
slots  84   440   1260   3168   6600   13104   23520   41888. (39)
```

There are `90,064` slots in all: `89,916` are positive, `148` are zero, and
none is negative.  The zeros are structural boundary/degree deficiencies.
By `(32)--(33)`, every layer-histogram coefficient of total degree between
five and twelve is therefore nonnegative, for an arbitrary number of
Gamma layers and every anchored support.

As independent shadows of this Schur certificate, the one-active-layer
Newton bank through degree eight has `276` positive and four structural-zero
slots, while the `45` two-active pairs `1<=k<=l`, `5<=k+l<=14`, have
`18,760` nonnegative Newton slots and no negative one.  The latter is a
separate exact arbitrary-support result beyond the Schur range proved here.
It does not cover histograms with three or more active layers in degrees
12--14.

A finite exact hostile scan on the six supports

```text
(1,2),(1,3),(2,3),(1,4),(2,5),(4,5)                  (40)
```

finds no negative Schur coefficient through degree fourteen (`5,952`
values, of which `2,832` lie in the still-open degrees thirteen and
fourteen).  This is evidence only outside the uniform degree-twelve
theorem.

## 5. The labelled Ewens current and the exact remaining topology

There is a useful reformulation of the open finite bank.  Start with five
disjoint labelled packets of `a` microletters and three disjoint labelled
packets of `b` microletters for `B1`, or four `b`-packets for `B2`.  A set
partition `pi` of the macro packets induces the corresponding partition of
their `5a+3b` or `5a+4b` microletters into block unions.  Distribute each
collected row coefficient uniformly over the labelled macro partitions of
that colour type.  If `omega(pi)` is the resulting signed weight, define
its upper zeta transform

```text
W(kappa)=sum_(pi>=kappa) omega(pi).                    (41)
```

For a permutation `sigma` of the microletters, let `orb_macro(sigma)` be
the least macro coarsening whose packet-block unions are `sigma`-invariant;
equivalently, join macro packets whenever a cycle meets both.  The cycle
enumerator identity

```text
P_pi(x)=sum_(sigma: orb_macro(sigma)<=pi) x^(defect(sigma))  (42)
```

shows that for any number of Gamma layers each coefficient is a positive
weighted sum of

```text
W(join_nu orb_macro(sigma_nu)).                        (43)
```

The tempting claim `W>=0` is false, but the failure is extremely rigid.
Exact enumeration of all `Bell(8)=4140` and `Bell(9)=21147` labelled
macro partitions gives

```text
                         positive   negative   zero
B1 rank-four current        285        195     3660
B2 rank-four current        720        900    19527.  (44)
```

Every nonzero value of `W` has atomic rank exactly four.  For `B1` its
positive minimum is `1/30` and its negative minimum is `-1/15`; for `B2`
they are `1/60` and `-1/30`.  The zero column in `(44)` counts all Bell
partitions; among rank-four partitions alone the zero counts are `1221`
and `5331`.

The positive and negative collected masses are balanced, but a
first-order transport using partition/root dominance has worst-chamber
capacity only `28`.  Thus any uniform argument of that kind leaves residual
masses `9` for `B1` and `11` for `B2`; special supports can add edges but do
not give a uniform full transport.  Bare partition-lattice positivity, and
every symmetric relabelling of it, fails.  At total defect four the signed
rank-four tree current cancels.  Degree five is the first place where a
rank-four join carries a redundant insertion edge.  The remaining problem
is therefore a relative rank-four cycle/edge-slide positivity theorem, not
an unstructured search through histograms.

The zeta expansion `(41)--(44)` applies to the undeleted bank.  Division by
the common response factor in `(11)` is triangular with constant term one,
but it is not a positive zeta operation; the rank-four current is therefore
a topological reduction of the residual `Phi`, not already a positive
decomposition of it.  This also explains why the proved Schur and
two-layer faces do not yet prove the full result: they test large explicit
faces of `(43)`, but do not provide the required cycle-space injection for
an arbitrary tuple of layers.

## 6. Exact scope

The theorem proves:

1. the universal `24/25` banks and their exact common-root quotients;
2. a unique Ferrers-dominant response and an explicit all-layer tail for
   every anchored triple;
3. every Schur/histogram coefficient through total degree twelve;
4. the complete two-active-layer face through total degree fourteen; and
5. the rank-four Ewens-current reduction, the irreducible `9/11`
   first-order transport deficit, and its sharp negative hostiles.

It does **not** prove positivity of every remaining finite histogram.  In
particular it does not yet prove arbitrary anchored product-Gamma
width-three goodness, the arbitrary three-slot case, SFC in width three,
GMC(2), NC2, LRC(14), JC(2), or DC(2).

## 7. Exact companions

The primary companion is reproduced by

```text
python 04-computation/gmc_product_gamma_arbitrary_anchored_tail_thm3110.py
python -O 04-computation/gmc_product_gamma_arbitrary_anchored_tail_thm3110.py
```

It derives the two multilinear banks, the polarized collision orders, the
four-chamber Ferrers comparisons, the exact common-root counters, the tail
gaps, the one- and two-layer Newton faces, and the complete labelled Ewens
zeta census.  Both executions byte-match the stored output.

The independent Schur/Newton route is reproduced by

```text
python 04-computation/gmc_product_gamma_arbitrary_anchored_schur_thm3110.py
python -O 04-computation/gmc_product_gamma_arbitrary_anchored_schur_thm3110.py
```

It rederives the banks and common divisors, evaluates the Schur coefficients
by exact power-sum/Newton and Jacobi--Trudi arithmetic, proves all `90,064`
generic chamber slots through degree twelve, replays the six-support
degree-fourteen hostile scan, and independently reconstructs the sharp
first-order transport capacity `28`.  Both executions byte-match its stored
output.  Every count and sign in this theorem is therefore an exact integer
or rational identity; the remaining all-degree step is mathematical, not a
floating-point extrapolation.

**QED candidate, pending independent hostile audit.**
