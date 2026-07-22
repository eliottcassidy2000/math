---
id: THM-2112
title: "Effective guard/speed ratios and all-height finiteness at rank seven"
status: >
  PROVED by a connected two-torus strip-tiling obstruction, made effective
  through the rank-six flat-cut parameter lattice.
  There is an absolute R_7 such that if the closed safe set of seven speeds
  lies in one open 1/7 guard comb, then some speed q has both entries of the
  reduced ratio (h/gcd(h,q),q/gcd(h,q)) at most R_7. One explicit choice is
  R_7=5*28^8*(7*57^42)^17. It strengthens the support-three relation
  alternative THM-2083 in arity, while THM-2085 remains stronger in giving
  the explicit coefficient height 57 for its three-term relation. The first-
  ratio statement alone does not close the resulting templates or LRC(14).
  Combined with
  THM-2087's complete short-relation cut, it also proves that only finitely
  many primitive contained packets can lie in the branch where every
  guard/speed relation has height greater than 57. Thus every infinite
  primitive obstruction family must stay on one of the explicit reduced
  guard ratios of height at most 57. For odd guards, a recursive residual-
  capacity argument then bounds all seven reduced ratios and proves that only
  finitely many primitive seven-speed guard containments exist. An explicit
  BV/Fourier recursion bounds all seven ratios and then the primitive terminal
  height. THM-2077/2092 then give the full depth-four original-row box
  max(S)<=floor((128/3)L_7^8). The resulting finite audit is numerical but far
  too large for a direct box run; it remains undone and LRC(14) is not closed.
source: codex-2026-07-22-LRC14-two-torus-strip-tiling
depends_on:
  - THM-2080
  - THM-2081
  - THM-2087
  - THM-2088
  - THM-2089
related:
  - THM-2082
  - THM-2083
  - THM-2085
  - THM-2086
  - THM-2090
  - THM-2091
  - THM-2092
  - THM-2094
  - THM-2095
  - THM-2096
  - THM-2097
  - THM-2098
  - THM-2099
  - THM-2114
---

# THM-2112 -- an effective guard/speed-ratio alternative

Put

```text
D_q={t in R/Z:||qt||<1/14},
E_h={t in R/Z:||ht||<1/7},
C_h=(R/Z) minus E_h,
G_Q=intersection_(q in Q)((R/Z) minus D_q).             (1)
```

For positive integers `h,q`, define the two-term relation height

```text
rho(h,q)
 =min{max(|a|,|b|):(a,b) in Z^2 minus {(0,0)}, ah+bq=0}
 =max(h/g,q/g),                         g=gcd(h,q).       (2)
```

Thus `rho(h,q)` is exactly the height of the reduced rational ratio `q/h`.

> **Theorem.** Put
>
> ```text
> R_7=5*28^8*(7*57^42)^17.                             (3a)
> ```
>
> Then, for
> every seven-element set `Q` of positive integers,
>
> ```text
> G_Q subset E_h     implies     min_(q in Q)rho(h,q)<=R_7.  (3)
> ```

No parity, divisor-completeness, hereditary-primitivity, or terminal-height
assumption is needed. In the THM-2073 depth-four application those additional
coordinates remain available after (3).

Thus a rank-seven obstruction cannot send **all seven** reduced guard/speed
ratios to increasing height. At least one speed belongs to the explicit finite
list of rational multiples of the guard of reduced height at most `(3a)`.

## 1. The compact connected-group limit

Suppose, for contradiction, that (3) has no uniform constant. Choose contained
packets

```text
G_(Q_n) subset E_(h_n),
Q_n={q_(n,1),...,q_(n,7)},
min_i rho(h_n,q_(n,i))->infinity.                       (4)
```

Label each set in increasing order and put

```text
v_n=(h_n,q_(n,1),...,q_(n,7)) in Z^8.                  (5)
```

Let `mu_n` be the pushforward of Haar measure on `R/Z` under

```text
t -> v_n t in (R/Z)^8.                                 (6)
```

For a character `m in Z^8`,

```text
mu_n_hat(m)=1 if m.v_n=0, and 0 otherwise.              (7)
```

Enumerate `Z^8` and take a diagonal subsequence on which every value in (7)
is eventually constant. Let

```text
L={m in Z^8:m.v_n=0 eventually},
K=L^perp subset (R/Z)^8.                                (8)
```

The limiting Fourier transform is `1_L`, so

```text
mu_n -> Haar measure m_K weakly.                        (9)
```

The subgroup `L` is saturated: if `k m in L` for a nonzero integer `k`, then
`k(m.v_n)=0` eventually, hence `m.v_n=0` eventually and `m in L`. Therefore
`Z^8/L` is torsion-free and `K` is connected.

Write `chi_0,...,chi_7` for the coordinate characters restricted to `K`.
Condition (4) says that `L` has no nonzero vector supported on `{0,i}`.
By duality, every map

```text
(chi_0,chi_i):K -> (R/Z)^2                              (10)
```

is surjective. In particular its pushforward of `m_K` is two-dimensional
Haar measure.

## 2. Containment becomes an exact seven-strip tiling

In `(R/Z)^8`, let

```text
C_0={x:||x_0||>=1/7},
B_i={x:||x_i||<1/14},
A_i=C_0 intersect B_i.                                  (11)
```

The set

```text
B=C_0 intersect intersection_i((R/Z)^8 minus B_i)       (12)
```

has `mu_n(B)=0` by (4). Its boundary lies in the finite union of coordinate
level sets `||x_0||=1/7` and `||x_i||=1/14`. Every coordinate marginal of
`m_K` is Haar, so this boundary has zero `m_K`-measure. Thus `B` is a
continuity set for (9), and

```text
m_K(B)=0.                                               (13)
```

Consequently the seven `A_i` cover `C_0` almost everywhere. Surjectivity in
(10) gives

```text
m_K(A_i)=(5/7)(1/7)=5/49,
sum_i m_K(A_i)=5/7=m_K(C_0).                            (14)
```

Equality of the sum and union measures forces multiplicity exactly one:

```text
1_(C_0)=sum_(i=1)^7 1_(A_i)       m_K-almost everywhere. (15)
```

In particular,

```text
m_K(C_0 intersect B_i intersect B_j)=0                 (16)
```

for every `i!=j`.

If the three characters `chi_0,chi_i,chi_j` were rationally independent,
their map from `K` onto the three-torus would be surjective. The measure in
(16) would then equal

```text
(5/7)(1/7)(1/7)=5/343>0,                               (17)
```

a contradiction. Hence every such triple is dependent. Fixing `i=1`, the
pair `chi_0,chi_1` is independent and every `chi_j` belongs to their rational
span. The eight coordinate restrictions generate the character group of
`K`, so

```text
rank(K)=2.                                              (18)
```

Thus any counterexample sequence to (3) would converge to an exact tiling of
the guard complement on one connected two-torus by seven centered character
strips.

## 3. Centered character strips cannot make that tiling

Let `Gamma=K_hat`, a free abelian group of rank two. Partition
`chi_1,...,chi_7` into rational rays in `Gamma tensor Q`, identifying opposite
rays because `||-x||=||x||`. For a ray `R`, choose its primitive character
`psi_R` and write

```text
chi_i=m_i psi_R,                    m_i in Z minus {0},
S_R(z)=sum_(i in R)1_{||m_i z||<1/14}.                  (19)
```

The integer-valued circle step function `S_R` is not constant. It equals
`|R|` near `z=0`, whereas multiplication by each nonzero `m_i` preserves Haar
measure, so

```text
integral S_R(z) dz=|R|/7.                              (20)
```

Therefore `S_R` has a boundary value `z_R` across which its net jump is
nonzero.

The characters `chi_0` and `psi_R` are rationally independent by (10), so

```text
(chi_0,psi_R):K -> (R/Z)^2                             (21)
```

is surjective. The boundary circle `psi_R=z_R` consequently meets the open
guard complement `||chi_0||>1/7` in open arcs. Boundary circles belonging to
other rational rays meet it in only finitely many points; the two guard
boundary circles do likewise. Choose a point on one of those open arcs away
from every such intersection.

In a small neighbourhood of that point, crossing `psi_R=z_R` changes `S_R`
by its nonzero jump, while every other ray contribution is locally constant
and both sides remain in the guard complement. The total multiplicity

```text
sum_R S_R(psi_R(x))=sum_i 1_(B_i)(x)                   (22)
```

therefore takes different constant integer values on the two adjacent open
cells. Equation (15) says it equals one almost everywhere on both cells, a
contradiction.

This rules out the sequence (4). The usual bad-sequence contradiction now
supplies one absolute `R_7`, proving (3). QED.

## 4. Quantitative chambers and the explicit constant

The compactness proof can be effectivized because THM-2087--2089 bound the
characters in the only branch where the cut matrix has a two-dimensional
kernel.

### A quantitative strip-arrangement lemma

Let

```text
a_0,a_1,...,a_7 in Z^2,
det(a_0,a_i)!=0,                  ||a_j||_infinity<=M,  (Q1)
```

and define the strict safe region

```text
U={x in (R/Z)^2:
     ||a_0.x||>1/7 and ||a_i.x||>1/14 for every i}.    (Q2)
```

Put

```text
D=28M^2.                                               (Q3)
```

Then `U` contains a Euclidean ball of radius strictly greater than

```text
1/(5D^8).                                              (Q4)
```

First, `U` is nonempty. If it were empty, the seven danger strips would cover
the guard complement almost everywhere. Every map

```text
x -> (a_0.x,a_i.x)
```

is onto the two-torus by `(Q1)`, so each restricted danger strip has measure
`5/49`. Their total measure is `5/7`, exactly the guard-complement measure.
The cover would therefore have multiplicity one almost everywhere. Partition
the `a_i` into rational rays. On a ray with primitive character `psi`, write
`a_i=m_i psi`. The circle step function

```text
S(z)=sum_i 1_(||m_i z||<1/14)                          (Q5)
```

is not constant: it equals the number of indices on the ray near zero, while
its mean is one seventh of that number. Hence some boundary value has a
nonzero net jump. Since `det(a_0,psi)!=0`, its boundary fiber meets the strict
guard complement in open arcs. Choose a point there away from the finitely
many other ray and guard boundaries. Crossing the fiber changes the total
multiplicity, contradicting multiplicity one on both adjacent open cells.

It remains to quantify one cell. Lift a nonempty chamber of `U` to `R^2` and
intersect its closure with a translated unit square. This also handles an
unbounded strip chamber: the square sides supply vertices. The lifted
arrangement lines have equations

```text
a_j.x=k plus-or-minus 1/14,
a_0.x=k plus-or-minus 2/14,                            (Q6)
```

and the four square sides have integral equations. At a vertex, Cramer's rule
uses two nonparallel integer normals of sup-height at most `M`; each vertex
coordinate consequently has reduced denominator at most

```text
14|det|<=28M^2=D.                                      (Q7)
```

The resulting positive-area convex polygon has a nondegenerate triangle whose
vertices are polygon vertices. A difference of two vertex coordinates has
denominator at most `D^2`, so the doubled area of that triangle, being a
nonzero rational number, has absolute value at least `D^(-8)`. Its area is at
least `1/(2D^8)`, while its perimeter is at most `3 sqrt(2)`. Its incircle
therefore has radius at least

```text
1/(3 sqrt(2) D^8)>1/(5D^8).                            (Q8)
```

The open incircle lies in the lifted chamber and proves `(Q4)`.

### Applying the lemma to a flat cut

Assume first that THM-2087's two-term branch does not occur. Choose its
complete-cut relations and form the THM-2088 matrix. If its rank is seven,
divide the whole packet by its common gcd. The primitive kernel vector is a
maximal-minor vector; common dilation preserves containment because circle
multiplication is onto, and it preserves every reduced pair height. THM-2088
then gives

```text
max(h,max Q)<=91421508108581.                          (Q9)
```

Suppose instead that the cut matrix has rank six. THM-2089 supplies

```text
q_i=(N_i z+R_i h)/D_i,
0<|N_i|,|D_i|<=A=57^6,        |R_i|<=6A,              (Q10)
```

on the integral parameter lattice `Lambda`, whose index satisfies

```text
J=[Z^2:Lambda]<=57^36.                                 (Q11)
```

Choose a Hermite-normal-form basis of `Lambda`. Every basis entry has
absolute value at most `J`. In its coordinates `n in Z^2`, the guard and
speeds are integral linear forms

```text
h=ell_0(n),                    q_i=ell_i(n),            (Q12)
```

and `(Q10)--(Q11)` give

```text
||ell_j||_infinity<=7AJ<=M:=7*57^42.                  (Q13)
```

Indeed, on either basis vector the numerator in `(Q10)` has magnitude at most
`(|N_i|+|R_i|)J<=7AJ`; division by `D_i` only decreases it, and it is integral
because the basis vector belongs to `Lambda`.

The guard form is independent of every speed form. Otherwise a nonzero row
supported on `(h,q_i)` would lie in the rank-six cut row space. But the
speed-only kernel vector `(0,r)` from THM-2088 has every `r_i!=0`; vanishing
on it kills the speed coefficient, and vanishing on the positive packet then
kills the guard coefficient. Thus the characters in `(Q12)` satisfy `(Q1)`.

Replace `n` by its primitive direction. This scales every frequency by the
same positive integer, so it changes neither containment nor any reduced pair
height. The geodesic

```text
L_n={tn:t in R/Z} subset (R/Z)^2                       (Q14)
```

has covering radius at most `1/(2||n||_2)`: it is the kernel of the primitive
character `(-n_2,n_1)`, and the distance to that kernel is at most the circle
distance `1/2` divided by `||n||_2`. Under containment, `L_n` misses the ball
in `(Q4)`. Hence

```text
||n||_2<=(5/2)D^8.                                    (Q15)
```

Every coordinate form in `(Q12)` is then at most `2M||n||_2`, so for every
speed in this flat-cut branch,

```text
rho(h,q_i)<=5MD^8
             =5*28^8*(7*57^42)^17.                   (Q16)
```

This dominates `57` and `(Q9)`. Together with THM-2087's pair branch and
THM-2088's rank-seven branch, `(Q16)` proves the explicit value `(3a)`. QED.

## 5. The no-short-pair branch is finite

There is a stronger consumer once THM-2087's complete cut is retained.
Normalize a packet by

```text
gcd(h,q_1,...,q_7)=1                                   (23)
```

and assume the no-short-pair branch

```text
rho(h,q_i)>57                    for every i.           (24)
```

THM-2087 supplies a nontrivial cut `A|B` of the seven labels and, for every
cross pair `i in A`, `j in B`, a relation

```text
a_ij h+b_ij q_i+c_ij q_j=0,
max(|a_ij|,|b_ij|,|c_ij|)<=57.                         (25)
```

Both `b_ij` and `c_ij` are nonzero by (24). Choose one such relation for each
cross pair and let `M` be their coefficient matrix on the eight columns
`(h,q_1,...,q_7)`. There are only finitely many possible matrices `M`.

Delete the guard column and call the resulting matrix `N`. It is a weighted
edge matrix of the connected complete bipartite graph on `A|B`, with both
endpoint weights on every edge nonzero. If `x in ker(N)`, the value at any
one vertex determines the value at every neighbour, and connectivity
propagates this determination through the graph. If the initial value is
zero, every value is zero. Hence

```text
dim ker(N)<=1,                  rank(N)>=6.             (26)
```

The actual positive packet lies in `ker(M)`, so `rank(M)<=7`, while (26)
gives `rank(M)>=6`.

If `rank(M)=7`, its kernel is a rational line and contains at most one
primitive positive integer vector. Thus a fixed `M` contributes at most one
packet.

Suppose `rank(M)=6`. Then (26) forces `rank(N)=6`. Let `0!=r in ker(N)`.
Every coordinate of `r` is nonzero: one zero again propagates across the
connected cut. Therefore

```text
(0,r) in ker(M),                                      (27)
```

and together with the actual vector `(h,Q)` it spans the two-dimensional
kernel of `M`.

No nonzero two-term guard/speed relation can belong to the row space of this
rank-six `M`. Indeed, a row supported on `(h,q_i)` that vanishes on `ker(M)`
must vanish on `(0,r)`. Its speed coefficient times `r_i` is then zero; since
`r_i!=0`, the speed coefficient is zero, and vanishing on the positive packet
then kills the guard coefficient as well.

THM-2112 supplies some label `i` and a nonzero pair row `ell` of height at
most `R_7` that annihilates the actual packet. There are only finitely many
choices of `(i,ell)`. The preceding paragraph makes `ell` independent of the
row space of `M`, so

```text
rank([M;ell])=7.                                       (28)
```

Again the common kernel is one rational line and contains at most one
primitive positive integer vector. Since both the cut matrices and the pair
rows come from finite sets, only finitely many primitive packets satisfy
(24) and guard containment.

Consequently there is an absolute effective `B_7` such that

```text
G_Q subset E_h and gcd(h,Q)=1
  implies
  [min_i rho(h,q_i)<=57] or [max(h,Q)<=B_7].            (29)
```

This proves the rank-six/rank-seven matrix input targeted by THM-2088 and
then consumes the rank-six plane using the two-torus ratio theorem. For
example, choose six independent cut-tree rows and append the pair row of
height at most `R_7`. Hadamard's inequality makes the primitive common kernel
vector at most

```text
B_7=2*3^3*57^6*R_7.                                   (29a)
```

Thus this no-short-pair finite box is numerical, although enormous.

## 6. Odd guards have only finitely many primitive contained packets

For the terminal application `h` is odd, and THM-2080 makes the preceding
one-ratio result recursive.

Order the seven values `rho(h,q_i)` increasingly. We prove by induction on
`k=1,...,7` that the `k`th value is bounded by an absolute constant over all
odd-guard contained packets. The case `k=1` is (3).

Suppose the result is known through `k`, where `1<=k<=6`, but no uniform
bound exists for the `(k+1)`st value. Choose contained packets for which the
first `k` heights stay bounded and every remaining height tends to infinity.
There are only finitely many labels and reduced ratios in the bounded bank,
so after passing to a subsequence we may assume that the same `k` labels have
fixed ratios

```text
q_j/h=r_j/s_j,                    j=1,...,k.             (30)
```

In the compact-group limit of Section 1, the characters
`chi_0,chi_1,...,chi_k` therefore lie on one rational character ray. Choose
its primitive character `psi`; all these characters are fixed nonzero integer
multiples of `psi`.

Let `R_k` be the part of the guard complement that is safe for these `k`
fixed danger characters:

```text
R_k={||chi_0||>=1/7 and ||chi_j||>=1/14 for j<=k}.      (31)
```

This is the pullback under `psi` of one fixed circle step set. Its Haar
measure `c_k` is positive. Indeed, after clearing the fixed denominators in
(30), the guard and the `k<=6` speeds are a common dilation of one integer
packet with odd guard. THM-2080 says that the safe set of at most six distinct
speeds is not contained in that guard, and in fact leaves positive measure
outside it. Common dilation preserves Haar measure, so

```text
c_k=measure(R_k)>0.                                    (32)
```

For every remaining label `i`, its pair relation height with the guard tends
to infinity. Hence `chi_i` is rationally independent of `chi_0`, and therefore
of `psi`; the map `(psi,chi_i)` is onto the two-torus. It follows that

```text
measure(R_k intersect {||chi_i||<1/14})=c_k/7.          (33)
```

As in Section 2, all threshold boundaries have zero limiting Haar measure,
so exact containment passes to the limit and says that the `7-k` remaining
danger strips cover `R_k`. The union bound and (33) would give

```text
c_k<=((7-k)/7)c_k<c_k,                                 (34)
```

which is impossible. This proves the induction step. Thus all seven reduced
guard/speed ratios belong to finite absolute banks.

For any fixed labelled ratio pattern, clear its denominators:

```text
h=d H,                     q_i=d Q_i.                  (35)
```

The integer tuple `(H,Q_1,...,Q_7)` is fixed by the pattern. If `Q` is
primitive, then

```text
1=gcd(q_1,...,q_7)=d gcd(Q_1,...,Q_7),                 (36)
```

so `d=1`. There are only finitely many patterns, hence only finitely many
primitive odd-guard packets satisfying `G_Q subset E_h`. Equivalently, an
absolute finite terminal-height bound exists.

This induction also has a numerical form. Define an increasing sequence
`L_1,...,L_7` by

```text
L_1=R_7,                                                (R1)

L_(k+1)=ceil( [98(7-k)/(5k(6-k))]
                  L_k^(2k)(1+kL_k) ),       1<=k<=5,    (R2)

L_7=ceil( (686/3)L_6^12(1+6L_6) ).                    (R3)
```

Then the `k`th smallest reduced guard/speed height is at most `L_k`.
To prove the induction step, suppose the first `k` ratios are at most `L_k`
and write them in lowest terms as in `(30)`. Put

```text
H=lcm(s_1,...,s_k),        h=dH,
Q_j=H r_j/s_j,             q_j=dQ_j.                  (R4)
```

Then

```text
H<=L_k^k,              Q_j<=H L_k.                    (R5)
```

On the base circle define

```text
F_k(x)=1_(||Hx||>=1/7)
       product_(j<=k)1_(||Q_jx||>=1/14),
c_k=integral F_k,
V_k=Var(F_k)<=2H+2 sum_j Q_j<=2H(1+kL_k).             (R6)
```

For a remaining speed `q`, reduce

```text
d=ga,                 q=gb,                 gcd(a,b)=1. (R7)
```

The BV Fourier estimate `|F_k_hat(m)|<=V_k/(2 pi |m|)`
and the danger coefficient bound `|s_n|<=1/(pi|n|)` give the exact
correlation estimate

```text
| measure(F_k(dt)1_(D_q)(t))-c_k/7 |
 <=sum_(ell!=0) V_k/(2 pi |b ell|) * 1/(pi |a ell|)
 =V_k/(6ab).                                           (R8)
```

Here multiplication by `g` preserves Haar measure, and the only Fourier
relations are `ma+nb=0`, namely `(m,n)=(b ell,-a ell)`. If
`s=gcd(H,b)`, then

```text
rho(h,q)=max(aH/s,b/s)<=Hab,                           (R9)
```

so a lower bound `rho(h,q)>L` makes the error in `(R8)` at most
`V_k H/(6L)`.

Containment says the `7-k` remaining dangers cover the residual set. If all
their reduced heights exceeded `L`, the union bound and `(R8)--(R9)` would
give

```text
c_k<=(7-k)(c_k/7+V_kH/(6L)).                          (R10)
```

Thus some remaining height is at most

```text
7(7-k)V_kH/(6k c_k).                                  (R11)
```

For `k<=5`, THM-2080's guard-centred Hunter star and its overlap floor give

```text
c_k>=5(6-k)/42.                                       (R12)
```

For `k=6`, the same theorem gives a uniform strict floor. Its exact proof
shows that every overlap other than the unique equality ratio `q=6h` is at
least `1/42+1/2940`. Six distinct speeds therefore satisfy

```text
c_6>=sum_(j<=6)measure(E_h intersect D_(q_j))-1/7
    >=5/2940=1/588.                                   (R13)
```

Substituting `(R5)--(R6)` and `(R12)` into `(R11)` gives `(R2)`; substituting
`(R13)` at `k=6` gives the coefficient

```text
7*588/(3*6)=686/3,                                    (R14)
```

which is `(R3)`. This proves the effective induction.

Finally, if `Q` is primitive, let `S=lcm(s_1,...,s_7)`. Clearing every ratio
simultaneously gives `h=eS` and `q_i=eS r_i/s_i`. Since `gcd(Q)=1`, `e=1`.
Consequently

```text
h<=L_7^7,                         max Q<=L_7^8.         (R15)
```

This is an explicit all-height finite box for primitive odd-guard rank-seven
containments. In the depth-four tower, THM-2077's terminal-to-original height
transfer, recorded explicitly in THM-2092, now gives

```text
max(S)<=floor((128/3)max(Q))
      <=floor((128/3)L_7^8).                           (R16)
```

Thus the entire original thirteen-speed row is in an explicit finite box on
this branch. QED.

## 7. Frontier effect and scope

THM-2083 proves that some guard/two-speed triple has a uniformly bounded
support-three relation, and THM-2085 gives the explicit height `57` for that
statement. The present theorem reaches lower arity effectively, at the price
of a vastly larger constant: some single reduced ratio

```text
q/h=b/a,                    1<=a,b<=R_7                 (37)
```

comes from a finite list. This is exactly the ratio coordinate exposed by
THM-2082's translated-prime-grid branch, so the two results can now be joined
without pretending that rank-one code support remembers projective residues.

For primitive odd-guard packets the theorem now bounds every ratio and the
common scale explicitly by `(R1)--(R15)`. It therefore turns the depth-four
rank-seven terminal, and by `(R16)` its full original row, into a numerical
finite audit. No THM-2090 splice is needed merely to prove this finiteness;
its global templates can still organize the exact audit. The bound is
deliberately coarse and far beyond a direct height-box run, so deciding every
row remains necessary. Thus this theorem does not close the depth-four
terminal or LRC(14).

The challenged assumption is that the relevant vertices are runners. In the
limit proof the faithful vertices are **rational character rays**; their
boundary circles preserve the exact outside-guard tiling predicate. Pairwise
outside-guard intersection is the observable, while rational proportionality
is the gauge. Orienting these rays as a tournament would discard the jump
sign that proves the contradiction. The faithful carrier is instead the
rank-two relation matroid together with its signed boundary arrangement.
