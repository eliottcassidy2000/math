# Finite scale certificates form a convolution calculus

Status: **PROOF CANDIDATE under root audit; finite-exact controls retained.**
This extends the exact global two-step pairing tree, not the original
Collatz map or longer-horizon pairing constraints.

## 1. Inheritance, map, and concept board

[THM-4491](../../01-canon/theorems/THM-4491-pairing-two-step-tree-extra-density.md)
gives the exact rooted constraint tree and its minimum lower density alpha.
[THM-4492](../../01-canon/theorems/THM-4492-pairing-two-cutoff-density-separation.md)
proves that finite weighted cutoff objectives have exact fringe limits,
refutes natural-density attainment of alpha, and gives upper density at
least `821510388809/2677850419968=0.30677978974599224...`.
The closest hostile example is `C(2)+C(5)=1` versus joint minimum two.
The corrected near miss is independent prefix optimization followed by
unjustified gluing at comparable scales. The least-used sidecar is the
exact common-assignment functional for an arbitrary finite positive kernel.

The live board is: common assignments; probability measures on scales;
convolution; finite-depth boundary loss; logarithmic averaging; and global
realization. Anchor: the exact dual calculus of cutoff certificates.
Niche: lower logarithmic density. Wildcard: the user's ten-vertex prompt
as a declared small legal-constraint carrier, with a rounding hostile test.

For a global two-step pairing let `F` be its flipped pair indices and
`A_F(X)=|F intersect[1,X]|`. Its exact clauses on indices at least two are

```
epsilon_(2k)+epsilon_(3k)=1,
epsilon_(3k+1)<=epsilon_(2k+1)<=epsilon_(3k+2).
```

Every legal finite prefix extends globally. Put `r=2/3`. A nonzero finite
nonnegative kernel `a=(a_0,...,a_K)` has

```
J_a(X)=min_F sum_k a_k A_F(floor(r^k X)),
Lambda_a=sum_k a_k r^k,
gamma_a=lim_(X->infinity) J_a(X)/X.
```

The limit and its exact nonnegative fringe series are THM-4492. Normalize
the kernel to a probability distribution `p_k=a_k r^k/Lambda_a` and define

```
B(p)=gamma_a/Lambda_a.
```

Conversely a finite rational probability vector gives `a_k=p_k r^(-k)`,
so the map loses only multiplication of all original prices by a positive
constant. It retains all scale locations and the common assignment.
Passing from the functional to its scalar value loses the minimizing
assignments; no subsequent global realization claim is made.

## 2. Concavity, translation, and convolution

**Proposition 1.** On finite rational probability distributions, `B` is
concave, invariant under translation by a nonnegative integer, and
Lipschitz in total variation:

```
|B(p)-B(q)| <= TV(p,q) = (1/2) sum_k |p_k-q_k|.
```

**Proof.** With `d_F(X)=A_F(floor X)/X`, one has exactly

```
J_a(X)/(Lambda_a X)=min_F sum_k p_k d_F(r^k X).
```

Each profile lies in `[0,1]`. Minimum of linear functionals is concave;
two probability vectors have equal total mass, so their difference against
any such profile is at most total variation. Both properties survive the
existing limits. For translation by `s`, compare the shifted objective at
`X` with the original objective at `Y=floor(r^s X)`. Nested and fused
floors differ by at most one entry for each fixed coefficient. Their total
cost difference is bounded independently of `X`; division by `X`, followed
by the normalization factor `r^s`, proves translation invariance.

The same argument extends `B` uniquely and continuously to finite real
probability vectors: approximate by rational vectors, using the uniform
finite-`X` total-variation estimate before taking the limit. No external
real-weight theorem is required.

**Proposition 2 (convolution domination).** For finite probability vectors
`p,q`,

```
B(p*q) >= max(B(p),B(q)).
```

Indeed `p*q` is the convex combination with weights `q_j` of the translates
of `p`. Concavity and translation invariance prove the first bound;
commutativity proves the other. In unnormalized prices this is

```
gamma_(a*b) >= Lambda_b gamma_a,
Lambda_(a*b)=Lambda_a Lambda_b.
```

There is also a useful exact finite-cutoff inequality:

```
J_(a*b)(X) >= sum_j b_j J_a(floor(r^j X)).
```

Nested floors never exceed the fused floor. Therefore the left objective
for each common assignment dominates the corresponding sum of shifted
objectives; minimum of the sum is at least sum of the separate minima.
The inequality can be strict from floor effects, common-assignment
incompatibility, or both. The controls below separate these mechanisms.

## 3. Uniform adjacent scales exhaust every finite positive kernel

Let `u_m` be uniform probability on `0,...,m-1`, and put
`B_*=sup_p B(p)` over all finite positive cutoff kernels.

**Proposition 3 (uniform completeness).**

```
lim_(m->infinity) B(u_m)=B_*,
B(u_m) >= B(p)-sum_k p_k min(k/m,1).
```

For a fixed `p`, convolution domination gives `B(p*u_m)>=B(p)`. A
translated uniform window differs in total variation from the original by
`min(k/m,1)`, hence

```
TV(p*u_m,u_m)<=sum_k p_k min(k/m,1).
```

Proposition 1 proves the displayed quantitative bound. For every fixed
finite `p` its error tends to zero. Thus `liminf_m B(u_m)>=B(p)` for every
`p`, while `B(u_m)<=B_*`; take the supremum only after the fixed-kernel
limit. This proves completeness without optimizing an unbounded list of
unrelated coefficients.

Moreover `m B(u_m)` is superadditive: split the uniform window of length
`m+n` into a window of length `m` and a translate of one of length `n`,
then use concavity. In particular `B(u_(jm))>=B(u_m)` for integer `j>=1`.
Neither this argument nor the experiments claim that `B(u_m)` is monotone
at every consecutive value of `m`.

The integer representative of `u_m` is precisely
`a_k=3^k 2^(m-1-k)` for `0<=k<m`. Thus the previously empirical adjacent-
scale family is a complete family for the entire finite positive-kernel
method. Completeness concerns certificates, not globally realizable
assignments or attainment of their supremum.

## 4. Lower logarithmic density also obeys every certificate

For each global pairing define

```
lower_log_density(F)=liminf_(X->infinity)
    (1/log X) sum_(i<=X, i in F) 1/i.
```

**Proposition 4.** Every global two-step pairing satisfies

```
lower_log_density(F) >= B_*
    >= 821510388809/2677850419968
     = 0.30677978974599224... .
```

Fix `p` and `epsilon>0` first. The exact minimum limit implies that for
all sufficiently large real `X`, every global `F` obeys
`sum_k p_k d_F(r^k X)>=B(p)-epsilon`. Put `X=exp(t)` and integrate over
a long interval of `t`. Each finite shift `k log r` changes the integral
of the bounded profile by only a bounded endpoint error. After dividing
by the interval length, the weighted sum has the same lower time average
as `d_F(exp(t))`. Let the interval length tend to infinity and then let
`epsilon` tend to zero. Abel summation gives

```
sum_(i<=X, i in F) 1/i
    = A_F(X)/X + integral_1^X A_F(x)/x^2 dx.
```

The first term is bounded by one; the integral is the logarithmic time
average just obtained. Thus lower logarithmic density is at least `B(p)`.
Only now take the supremum over the fixed kernels. The exact rational
bound is the inherited eight-scale, depth-twenty certificate, independently
recomputed here using full fixed-root cost arrays.

This strengthens the density consequence of the finite-scale method.
It does not prove `B_*` is the minimum logarithmic, natural, or upper
density: each equality would require a matching globally legal construction.

## 5. Finite depth has a precise delay, not shift invariance

Let `B_H(p)` be the fringe series through depth `H`, with `B_H=0` for
negative `H`. If `S^s p` translates the probability kernel by `s`, then

```
B_H(S^s p)=B_(H-s)(p).
```

Before depth `s` all root weights vanish. Afterwards the normalized root
weights equal `r^(-s)` times the old pattern with depth reduced by `s`.
The same is true of differences and tolls; the residue period is repeated
`2^s` times. Hence the shifted toll sum is `3^s` times the old one, exactly
canceling the denominator's depth shift in the series.

The finite partial functional is itself concave. Indeed, writing
`T_H=sum_(i mod 2^H) min_b F_H^b(i)`, the exact period count gives
`T_H=B_H_toll+3T_(H-1)`: every child residue occurs three times.
Thus the partial series is `T_H/3^(H+1)` before normalization, a sum of
minima of linear cost functionals. Consequently the correct finite-depth
version of convolution domination is

```
B_H(p*q) >= sum_j q_j B_(H-j)(p).
```

A kernel translated beyond depth `H` has zero retained certificate even
though its exact infinite-depth value is unchanged. This is the cheapest
hostile example against reading fixed-depth comparisons as exact ordering
of kernel strengths. The lost coordinate is the depth still needed to
reach the translated cost, not failure of exact convolution domination.

## 6. Exact ten-vertex and compatibility controls

The legal graph on pair indices `2,...,11` has ten vertices and nine edges:

```
complement: (2,3), (4,6), (6,9);
order:      4<=3, 3<=5, 7<=5, 5<=8, 10<=7, 7<=11,
```

where the order is between bits at those indices. Exhausting all 1,024
assignments leaves exactly sixteen legal assignments, each globally
extendible. With `a=(1,0,1)`, at `X=11` the separate objectives have costs

```
min(A(11)+A(4))=3,       min(A(4)+A(1))=1.
```

Their common nested objective still has optimum four. The fused
convolution objective has minimum five because its last cutoff is two:

```
floor(r^2 floor(r^2*11))=1,       floor(r^4*11)=2.
```

The apparent strict gain in this ten-vertex example is entirely a floor
boundary effect. At the exactly aligned cutoffs `486,216,96`, the same
kernel instead gives

```
min(A(486)+A(216))=210,
min(A(216)+A(96))=94,
min(A(486)+2A(216)+A(96))=306 > 210+94.
```

Here there is no floor discrepancy: the extra two are a genuine cost of
one common assignment. No global minimality of this witness is claimed.
The ten-vertex object is a constraint tree, not a tournament: the native
relations are complement and partial order clauses, and most vertex pairs
have no relation. Orienting all absent pairs would add unsupported data.

Reproduce with `python3 04-computation/experiments/crossroads10_20260926_flow.py`
or with `-O`. The [output](crossroads10_20260926_flow.out) also records exact
finite convolution checks for cutoffs 1 through 1,000, finite-depth shift
and convolution controls, uniform-window total-variation checks, and the
independent inherited rational certificate. This is structural progress
in the finite positive-certificate method; its global realization gap and
the longer-horizon legal-constraint carrier remain open.

## 7. The minimum lower logarithmic density is attained by a global pairing

This is an attainment theorem for a precisely defined constant `h`, not
an identification of `h` with `B_*`. Work at integer cutoffs and set the
unconstrained bit `epsilon_1=0`; changing that single bit affects harmonic
cost by at most one and no asymptotic density. Define

```
H(X)=min_(global P_2 F) sum_(2<=i<=X, i in F) 1/i,
h=liminf_(integer X->infinity) H(X)/log X.
```

**Proposition 5 (harmonic prefix extension and attainment).** A fixed legal
prefix through an integer `Y>=2` extends through every integer `X>Y` with
harmonic cost at most

```
H(X) + sum_(2<=i<=Y) epsilon_i/i + 3
    <= H(X)+log Y+3.
```

Consequently a global two-step pairing attains lower logarithmic density
`h`, and every global two-step pairing has lower logarithmic density at
least `h`. In particular `h>=B_*`; the reverse inequality remains open.

**Proof of the extension bound.** In the finite tree cut at `X`, assign
weight `1/i` to vertex `i`. Let `Delta_i` be the difference of the optimum
subtree costs with root bit one and zero. Then

```
|Delta_i| <= 3/(i-1)                 (2<=i<=X).
```

The leaf case is immediate. Each child `j` has
`j-1 >= (3/2)(i-1)`, so its gap is at most `2/(i-1)` by induction.
For an even root the recurrence is `Delta_i=1/i-Delta_j`. For an odd
root it is `Delta_i=1/i+min(Delta_left,0)+max(Delta_right,0)`.
The negative and positive terms have opposite signs, so the absolute
value of their sum is bounded by the *maximum* child bound, not the sum
of two child bounds. Adding `1/i<=1/(i-1)` proves the claim.

Remove the fixed prefix. The remaining forest roots are exactly children
`j>Y` of parents at most `Y`. Each satisfies

```
Y < j <= floor((3Y+1)/2).
```

The fixed parent allows at least one child bit, independently for each
child. Requiring a particular allowed bit costs at most `|Delta_j|`
above the free optimum of that component. The unrestricted optimum
`H(X)` is at least the sum of these free component optima, since all
omitted prefix costs are nonnegative. Hence the entire prescribed prefix
costs at most its own harmonic weight plus

```
sum_(frontier j) 3/(j-1)
    <= 3(Y+1)/(2Y) <= 9/4 < 3
```

above `H(X)`. This proves the displayed extension inequality.

**Proof of attainment.** Every fixed global assignment has cost at least
`H(X)` for each integer `X`, hence lower logarithmic density at least `h`.
Starting from a finite legal prefix through `Y`, choose an arbitrarily
large integer `X` on the defining liminf sequence for `H`, so that
`H(X)/log X<=h+epsilon` and `log Y/log X<=epsilon`. Extend the old prefix
using the proved conditional bound, then repeat with `epsilon` tending
to zero and increasing cutoffs. The nested legal prefixes define one
global legal assignment. At the chosen cutoffs its harmonic ratio tends
to `h`; the universal lower bound gives equality of the liminf. Integer
and real cutoffs have the same asymptotics by rounding.

Thus global lower-logarithmic attainment needs no additional minimax
assumption. The remaining sharpness question has become the concrete
finite-optimization obligation

```
liminf_X H(X)/log X <= B_*.
```

The script checks the harmonic root-gap bound with exact fractions at
every node through cutoff 233, checks the cutoff-eleven harmonic optimum
against all sixteen legal assignments, and verifies all sixteen prefixes
extend through 233 within the exact proved bound.

## 8. A shared phase blocks a generic minimax shortcut

Compactness, finite prefix extendibility, and the kernel identities above
do not by themselves imply that `B_*` equals the minimum logarithmic
density. A family of actual subsets of the positive integers provides a
hostile control. Put `R=3/2`, `lambda=sqrt(R)`, and for `theta` modulo one
let

```
F_theta = {n>=1: fractional_part(log_R n-theta) lies in [0,1/2)}.
```

Equivalently its active real bands are
`[R^(j+theta),lambda R^(j+theta))`. One may take the closure of the
indicator sequences in the product topology. The closure only adds
boundary choices: for each finite range, a convergent phase subsequence
stabilizes membership away from band boundaries. At a fixed phase the
set of integer boundary points has size `O(log X)` through `X`, and
their reciprocals have bounded total sum. Thus every added sequence has
the same limiting logarithmic density as its phase representative. Every
allowed finite prefix extends by definition.

For `z=X/R^(j+theta)` in `[1,R)`, direct geometric summation gives,
uniformly in phase up to `O(log X/X)`, the prefix density profile

```
a=(lambda-1)/(R-1)=sqrt(6)-2,
g(z)=1-(1-a)/z          for 1<=z<lambda,
g(z)=R a/z             for lambda<=z<R.
```

Its minimum is `a`. Scaling `X` by any integer power of `r=1/R` leaves
this phase unchanged. Therefore for every finite probability kernel,

```
lim_(X->infinity) min_theta sum_k p_k d_(F_theta)(r^k X)
    = min_z g(z) = sqrt(6)-2.
```

All finite kernels have the same exact certificate. Yet every member of
the family has logarithmic density one half: each active band occupies
exactly half of its logarithmic period, and integer endpoint errors have
summable reciprocal size. Equivalently,
`integral_1^R g(z) dz/z=(1/2) log R`. Thus

```
sup_(finite kernels) B(p)=sqrt(6)-2 < 1/2
    =minimum lower logarithmic density in this family.
```

The lost coordinate is a common phase within a multiplicative period.
Translations by whole scale steps cannot observe it. This is a hostile
family for a proposed generic inference, not a counterexample within the
two-step pairing tree. To prove `h=B_*` for that tree one must exploit its
specific freedom to assemble assignments across such phases, or otherwise
prove the remaining harmonic-optimum bound. Replacing that obligation by
compactness or by a finite-kernel minimax slogan would leave the same gap.
