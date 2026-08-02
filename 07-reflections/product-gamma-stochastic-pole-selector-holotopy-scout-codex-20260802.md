# Product-Gamma stochastic pole-selector holotopy: exact repair polytope and portability wall

**Status:** `FINITE-EXACT REFLECTION / THEOREM-GRADE CANDIDATE / NO GMC OR NC2 CLOSURE`

**Date:** 2026-08-02

**Exact companion:**
[`04-computation/gmc_stochastic_pole_selector_polytope_scout.py`](../04-computation/gmc_stochastic_pole_selector_polytope_scout.py)

**Stored transcript:**
[`05-knowledge/results/gmc_stochastic_pole_selector_polytope_scout.out`](../05-knowledge/results/gmc_stochastic_pole_selector_polytope_scout.out)

## 1. Inheritance and the failed deterministic implication

The starting point is the fixed $Q$, one-letter virtual-prefix current from
the pole/partition bicomplex scout.  At support $(a,b)=(1,3)$, bank $I_2$,
every deterministic choice of the first physical pole
$r\in\{1,\ldots,8\}$ eventually fails the partition-upset inequalities by
degree $8$ or $9$.
THM-3128 explains why the first failure cannot be repaired inside the full
labelled-deletion fibre: the load-bearing top upset is invariant on that
fibre.

The first failed implication is therefore precise:

> positivity of the scalar pole flag does not imply that a deterministic
> first virtual letter has a Hasse-positive Young-type current.

The present scout changes the selector object rather than changing the
deletion fibre.  It replaces one pole $r$ by a probability law
$\lambda=(\lambda_1,\ldots,\lambda_8)$ and averages the eight already-defined
one-letter currents.  This is a finite convex relaxation, not a pole ordering.

## 2. The exact fibre polytope

For degrees $N=5,\ldots,9$, the partition lattices have respectively

\[
10, 27, 47, 168, 573
\]

upsets, including the empty upset.  Of the resulting 825 exact row
inequalities, 815 are nonzero and pairwise primitively distinct.  Exactly 812
are coordinatewise nonnegative on the pole simplex.  The three remaining
primitive rows are

\[
\begin{aligned}
A={}&(83547971350,15688221032,-569502791,-2926894294,\\
&-14620631052,-44124052840,-85206825873,-116934942806),\\
B={}&(77268972364224,14817264695349,-301139382544,\\
&-1779540733230,-11604925470240,-38159896577461,\\
&-74842829308296,-98601655955064),\\
C={}&(-133125723152,-49818680888,2688744487,18251699828,\\
&112450828200,398996381528,867825875007,1262919484552).
\end{aligned}
\]

The apparent middle wall is redundant.  With

\[
\alpha=\frac{8315869923144}{9467425097},
\]

every coordinate of $B-\alpha A$ is nonnegative.  Consequently the full
finite selector polytope is exactly

\[
\mathcal P_{1,3;I_2}^{5:9}
=\Delta_7\cap\{A\lambda\ge0\}\cap\{C\lambda\ge0\}.
\]

Exact active-set enumeration gives 24 vertices.  Every vertex has support
exactly two.  More sharply, there are exactly two boundary points on every
cross edge

\[
\{r_i,r_j\},\qquad i\in\{1,2\},\quad j\in\{3,4,5,6,7,8\},
\]

and no other vertices.  No deterministic pole lies in the polytope, so support
two is minimal.  The stored transcript gives all 24 rational endpoints.

On the $r=1,r=8$ edge the exact feasible interval for the $r=1$ mass is

\[
\frac{3077235337}{5275866162}
\le p\le
\frac{157864935569}{174505650963}.
\]

In particular

\[
\lambda=\frac35\delta_1+\frac25\delta_8
\]

has exact saturated Hasse transport in every degree $5,\ldots,9$.

## 3. Why the convex repair is not a synthetic pole

On the scalar row polynomial, averaging the one-pole cofactors only sees the
mean:

\[
\sum_r\lambda_r(1-rt)=1-t\,\mathbb E_\lambda r.
\]

That makes it tempting to replace the law by one synthetic barycentric pole
$\bar r=\mathbb E r$.  Young types retain the whole moment sequence instead:

\[
\mathbb E\,p_k[X-r]=p_k[X]-\mathbb E(r^k),
\qquad
p_k[X-\bar r]=p_k[X]-\bar r^k.
\]

For the $(3/5,2/5)$ law,

\[
\bar r=\frac{19}{5},\qquad
\mathbb E(r^2)=\frac{131}{5},\qquad
\operatorname{Var}(r)=\frac{294}{25}.
\]

The averaged current passes all five degrees.  The deterministic virtual pole
$r=19/5$, which has the same scalar cofactor, passes $N=5,6,7$ but fails
$N=8,9$.  Thus variance is not cosmetic: the convex selector repairs a
scalar-equivalent but type-distinct current.

This is the cleanest conceptual content of the computation.  The Young
refinement is a moment detector for pole laws.

## 4. Portability: a positive chamber result and a support wall

The repair is not merely an $I_2$-only accident.  At the same support
$(1,3)$, the simple law

\[
\lambda_1=\frac9{500},\qquad
\lambda_3=\frac{487}{500},\qquad
\lambda_8=\frac4{500}
\]

passes every exact upset and max-flow test for both banks $I_1,I_2$ and all
degrees $5,\ldots,9$.

It is nevertheless not a universal law across supports.  At support $(1,2)$,
bank $I_2$, let $\lambda$ be any probability law on the five physical pole
values $1,\ldots,5$.  Two necessary upset rows, one from $N=6$ and one from
$N=9$, are

\[
\begin{aligned}
R_6&=(2469992,-986920,-2955376,-3435376,-2426920),\\
R_9&=(-60532076544,1282401120,16315005312,
      -34474521120,-48239160000).
\end{aligned}
\]

But

\[
10000R_6+R_9=
(-35832156544,-8586798880,-13238754688,-68828281120,-72508360000)
\]

is strictly negative in every coordinate.  If both required inequalities were
nonnegative on $\lambda$, their displayed positive combination would be
nonnegative too, a contradiction.  Hence **no probability law on the physical
one-pole deletions works even in this single support/bank fibre** through
degrees $5,\ldots,9$.

A seductive portable rule illustrates how sharp the wall is.  Choosing the
smallest physical pole with weight $3/5$ and the largest with weight $2/5$
passes 1,141 of the 1,150 exact support/bank/degree cases in the THM-3120
universe.  The remaining nine failures are recorded exactly in the transcript.
Near-universality is not universality.

## 5. Degree ten empties the selector space

The cutoff at degree nine is load-bearing.  At the original support
$(1,3)$, bank $I_2$, degree ten has 42 partitions and 3,588 upsets.  Only six
upset rows have a negative coordinate, but two nested rank-tail upsets already
exclude every probability law on the eight physical poles.  Let

\[
U_2=\{\mu\vdash10:\ell(\mu)\le8\},\qquad
U_1=\{\mu\vdash10:\ell(\mu)\le9\}.
\]

Their raw pole rows are

\[
\begin{aligned}
R={}&(-105427024770557952,-41110096821098400,1323165981189312,\\
&10490179768896384,84003185593589760,328063977223727328,\\
&734662852146861120,1034128690334272512),\\
S={}&(2135515427122176,1594234834827264,-128048320760256,\\
&-1149608741796864,-8750331832665600,-36573564507761664,\\
&-89737341671630016,-137769965901950976).
\end{aligned}
\]

The positive Farkas combination

\[
R+11S=(-81936355072214016,-23573513637998496,-85365547173504,
-2155516390869120,-12250464565731840,-74245232361650976,
-252447906241069056,-481340934587188224)
\]

is strictly negative in all eight coordinates.  Hence the degree-ten
selector polytope is empty: no probability mixture of the physical one-pole
currents survives.  In particular the $3/5,2/5$ law first fails at degree ten,
with exact deficit $269133385522535424/5$ and first unresolved type
$(7,1,1,1)$.

Two cheap augmentations also fail sharply.  The empty/no-deletion state has
zero coordinate on both $U_2$ and $U_1$, so the same Farkas combination forces
all positive one-pole mass to vanish; the lazy selector polytope is only the
trivial all-empty point.  And for the 33 physical two-pole submultisets, the
single upset

\[
U_3=P_{10}\setminus\{(1^{10}),(2,1^8),(2,2,1^6)\}
\]

has strictly negative response on every state, ranging from
$-2197919641631883264$ to $-64486707449419008$.  Thus jumping directly to a
random physical two-pole prefix does not repair degree ten either.  Any viable
augmentation must mix prefix depths with a genuinely new reference state, or
leave the positive physical-prefix simplex.

## 6. Holotopy interpretation and exact typing boundary

The probability simplex is a useful “filled” selector space: deterministic
pole deletions are its vertices, while a law is a measure-valued barycentre.
The two exceptional upset functionals cut this selector simplex by two walls;
the 24 two-pole vertices are its sharp edge contacts.  In that limited sense
the repair is a convex or measure-valued holotopy between failed deterministic
selectors.

What is preserved:

- the support, bank, dominant quotient alphabet $Q$, and degree;
- the already-defined one-letter virtual-subtraction currents;
- zero total mass and every exact partition-upset inequality;
- hence the finite Hasse-transport certificate after averaging.

What is destroyed or absent:

- a deterministic pole flag or ordered factor removal;
- a proof that the weights arise from a positive decomposition of the original
  scalar response or of a PSD operator;
- consistency between different supports or between successive prefix depths;
- a canonical coupling of the random first pole to later pole choices;
- any implication for arbitrary radial polynomial coefficients.

The comparison across $I_1,I_2$ is an exact arithmetic comparison of currents
at a fixed support.  It is not, by itself, a common physical selector until a
lawful common-root/operator typing is reconstructed.  Likewise, a law must not
be silently replaced by its mean pole: the variance hostile proves that this
changes the Young current.

## 7. Consequence for the GMC/NC2 frontier

This result opens one route and closes another.

The opening is that deterministic virtual-prefix positivity is stronger than
necessary at a fixed fibre.  A positive measure over first-pole removals can
restore the complete finite Hasse cone even when every atom fails.  If a future
product-Gamma representation supplies these weights as an honest positive
mixture, the partition-current obstruction is not terminal.

The closure is that no support-independent stochastic first-pole recipe can be
assumed.  The exact $(1,2),I_2$ Farkas pair excludes every law supported on its
physical poles.  Therefore this scout does not prove NC2, does not handle
arbitrary radial coefficients, and does not weaken the scalar product-Gamma
positivity theorem.  It identifies a new missing sidecar: a representation
whose selector measure may use more structure than a probability law on one
physical pole.

## 8. Cheapest next decisive tests

1. **Two-step selector coupling.**  Ask whether a feasible first-pole law at
   $(1,3),I_2$ extends to a probability law on ordered pairs of distinct pole
   labels whose first and second averaged currents are both Hasse-positive.
   This is a finite martingale/coupling problem and is the natural test of
   whether the convex repair can become a genuine stochastic pole flag.

2. **Signed or augmented support at the Farkas wall.**  The $(1,2),I_2$
   certificate rules out positive laws on the five physical one-letter maps.
   Test the cheapest enlargement: the empty deletion, a two-letter deletion,
   or a Schur-positive virtual letter.  The strict Farkas vector tells exactly
   how much new direction is required.

3. **Operator realization.**  Search the product-Gamma integral for a latent
   random variable whose conditional first-pole law produces the exact
   $A,C$ halfspaces.  A scalar identity depending only on the mean is
   insufficient; at least a second-moment sidecar must survive.

4. **Minimal degree-ten augmentation.**  Since $R+11S$ is strictly negative
   on every physical one-pole direction, add the empty deletion and then the
   physical two-letter deletions.  Compute the smallest new support and exact
   mass needed to cross this separating wall.
