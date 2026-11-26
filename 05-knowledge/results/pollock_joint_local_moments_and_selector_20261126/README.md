# Pollock's joint moment filter and a positive-window selector

Extend the [cubic-defect packet](../pollock_five_clebsch_and_defect_20261126/) and the [marked Richelot golden pair](../pollock_ogg_richelot_golden_reduction_20261126/). For five indices $x_i=c+y_i$ at an integer center, with $\sum_i y_i=0$, put
\[
Q=\frac12\sum_i y_i^2,\qquad B=\frac16\sum_i y_i^3,\qquad
N=5F(c)+cQ+B,\quad F(x)=(x^3-x)/6.
\]
Both moments are integral for integer offsets: use parity for $Q$ and $y^3\equiv y\bmod6$ for $B$. Equivalently $Q=-e_2(y)$ and $B=e_3(y)/2$ on sum zero. The previous uncentered defect is $D=150B$ here.

**Results.** There is a complete joint test over $\mathbb Z_p$ for $p\ge5$, and the scalar target at a fixed integer center is locally unrestricted at every prime including 2 and 3. For the actual bounded integer problem, an elementary two-square inversion finds or excludes a prescribed moment pair. In sharp contrast with local universality, almost all targets in an $o(c^{5/3})$ growing window above the central value cannot be represented **at that center**. This is not a statement about unrestricted Pollock representations.

## 1. Complete $p$-adic joint test for $p\ge5$

First reduce $(Q,B)$ to $(q,b)$ mod $p$. The missing pairs are:

| $p$ | Missing $(q,b)$ |
|---|---|
| 5 | $(0,b)$ with $b\ne0$ |
| 7 | $q\in\{0,3,5,6\}$ and $b=\pm2$ |
| 11 | $(1,\pm5),(3,\pm2),(4,\pm4),(5,\pm1),(9,\pm3)$ |
| 13 | $(0,\pm2),(0,\pm3)$ |
| $p\ge17$ | none |

Every present pair lifts to all its higher digits except $(0,0)\bmod p$ at $p=7,19,23,47$. For those four exceptions every solution is zero in all coordinates mod $p$. Require $p^2\mid Q$ and $p^3\mid B$, divide the pair by those powers respectively, and recurse at the same prime. The identically zero pair always has the all-zero solution. For a nonzero pair the descent terminates. This is necessary and sufficient in $\mathbb Z_p$, without positivity. Neither absence list nor exceptional recursion is a classification of Ogg primes.

A residue tuple is regular exactly when it has at least three distinct entries: the Jacobian minor on three coordinates for $(\sum y_i,Q,B)$ is half a Vandermonde determinant. Each regular tuple therefore lifts to an arbitrary higher-digit pair, with $p^2$ extensions per level after retaining the sum condition. The replay exhausts all sorted centered residue multisets for primes $5\le p\le61$: enumerate $a\le b\le c\le d$, force the last residue by the sum and discard if the forced residue is smaller than $d$. This finds the table and verifies that the four singular-only pairs at the stated primes have just the zero tuple. The following argument makes the remaining range unconditional.

### K3 over a nonzero pair, Bring curve over the zero pair

For $p>5$ use the compactification cut out by $\sum y_i=0$, $Q(y)=qt^2$, $B(y)=bt^3$ in $\mathbb P^5$ (equivalently a quadric and cubic in $\mathbb P^4$). At infinity is the smooth genus-four complete-intersection [Bring curve](https://link.springer.com/article/10.1007/s40879-023-00706-0). Indeed a gradient dependence forces at most two coordinate values. A centered two-value nonzero vector has values $(5-m)d,-md$ with multiplicities $m,5-m$ and $Q=\frac52m(5-m)d^2\ne0$ for $p>5$; a constant centered vector is zero. Thus the infinity curve has no such point and is smooth; connectedness and genus follow from the complete-intersection formulas. Smoothness there also implies smoothness of the projective surface at infinity.

For $(q,b)\ne(0,0)$ any singular point in the affine chart must be a two-value vector. Up to swapping the groups the possibilities are:
* $(4d,-d,-d,-d,-d)$ gives $(q,b)=(10d^2,10d^3)$ and five points on $q^3=10b^2$.
* $(3d,3d,-2d,-2d,-2d)$ gives $(q,b)=(15d^2,5d^3)$ and ten points on $q^3=135b^2$.

The two cusps are disjoint for nonzero pairs at $p>5$ because $135-10=125$, and their parameter is recovered from $b/q$. At a singular point with values $a,a'$, the quadric tangent space consists of vectors whose sums within each of the two blocks vanish. Subtract $(a+a')/2$ times the quadratic equation from the cubic to cancel the linear part on the sum-zero hyperplane. The Hessian restricted to the quadric tangent is proportional to the difference of the block squared norms. Its nonzero blocks are nondegenerate since both block sizes are units. Thus every singular point is an ordinary double point. The quadric is smooth on the intersection at these affine points; away from them regularity or the infinity check applies.

The projective intersection has dimension two: the smooth curve at infinity has dimension one, so no component can have larger dimension, and the equation count gives the lower bound. The complete-intersection Koszul resolution gives geometric connectedness, trivial canonical and $H^1(\mathcal O)=0$; isolated singularities and the complete-intersection depth condition give normality. Resolving the ordinary double points preserves canonical and $H^1(\mathcal O)$, so the resolution is a K3 (cf. [Auel–Scheible's complete-intersection model](https://math.dartmouth.edu/~auel/papers/docs/AS_ciK3.pdf)). Its Weil lower bound is $p^2+1-22p$. Remove at most $p+1+8\sqrt p$ rational points at infinity and $10(p+1)$ over affine nodes, since an exceptional conic has at most $p+1$ rational points. A positive lower bound $p^2-33p-10-8\sqrt p>0$ for $p\ge37$ produces a regular affine tuple for every nonzero pair. For the zero pair one instead needs a Bring point. The genus-four Hasse–Weil lower bound is positive for $p\ge67$. The finite replay checks the remaining primes and also compares direct counts to $p+1-4a_p$ for the conductor-50 class, using $Y^2+XY+Y=X^3-126X-552$ and the [Bring isogeny description, Proposition 4.5](https://link.springer.com/article/10.1007/s40879-023-00706-0). That comparison is a control, not an assumption in the exhaustion.

## 2. Golden line vs node in characteristic five

The quadratic form mod 5 has the constant vector as radical. On the slice $y_0=0$ its matrix is $(I+J)/2$ in three independent coordinates, with determinant 3. A slice enumeration gives 25 vectors of norm zero; adding the constant translates gives 125 tuples, exactly the 120 permutations of all residues and the five constant tuples (all have both moments zero mod 5). For a nonzero $q$ the slice has 20 vectors per square $q$ and 30 per nonsquare $q$. Translation by $r$ preserves $Q$ and changes $B$ by $rq$ mod 5, so the same counts hold at each $b$. There cannot be exactly two residue values in a centered tuple, since the first block size is invertible mod 5 and the other is its negative. Thus all nonconstant tuples are regular.

The Richelot-selected nodal line identifies the five coordinates with $i\in\mathbb F_5$ in the selected cyclic order. Its 20 tuples $y_i\equiv ai+r$ with $a\ne0$ all lift with independently prescribed $Q,B\in5\mathbb Z_5$. The Vandermonde step is the joint lift, not just a phase match of $cQ+B$. The replay exhausts all $25^4=390625$ centered tuples mod 25 and checks both the full regular counts and the designated line seed.

In contrast a nonzero constant residue tuple cannot lift to the Clebsch cubic mod 25. Write $y_i=r+5z_i$, so $\sum z_i=-r$. Direct expansion gives $Q\equiv-5r^2/2$, $B\equiv-5r^3/3\bmod25$; in particular $v_5(B)=1$ and $v_5(D)=3$. Thus the scheme point at the node and the regular points of the specialized golden line have different lifting behavior.

## 3. Fixed-center scalar universality, but almost no nearby exact targets

For every $p\ge5$, the affine line $b=N-5F(c)-cq$ in the residue pair plane meets the regular locus. This follows from the table: at $7$ a horizontal line only loses at most four entries (plus at most one singular-only entry on the zero row), and a nonhorizontal line loses at most two from the missing $b$ values and at most the zero pair. At every other prime fewer than $p$ pairs are excluded from the full plane. Choose one regular pair, prescribe its higher digits with $cQ+B$ equal to the target difference, and lift.

At 3 use the integral polynomial $cQ+e_3/2$ on sum zero. The derivative in a transfer between two coordinates is the difference of $cy+y^2/2$ evaluated at them. For each $c\bmod3$ and each target residue, the short regular seeds are given and verified in the replay; Hensel solves the scalar equation. At 2 take offsets $(2t,1-2t,-1,0,0)$. The change in $F(x)$ for a step $2t$ is the $\mathbb Z_2$-polynomial $t(3x^2-1+6xt+4t^2)/3$; for the two adjacent base indices the combined change has slope one mod 2. Hence this one variable suffices for a scalar 2-adic lift. No complete joint classification at 2 or 3 is claimed.

For a concrete bounded residue consequence, for any modulus $M$ combine the centered local witnesses by CRT on the offsets with period $L=M$ if $M$ is odd and $L=2M$ otherwise. In the sum-zero hyperplane the polynomials $Q=-e_2$, $B=e_3/2$ show these periods suffice. Choose representatives in $[0,L-1]$, whose sum is $kL$ for $0\le k\le4$. Subtract $L$ from the $k$ largest. The resulting offsets sum to zero and have absolute value at most $4L/5$: the largest unshifted value is at most $kL/(k+1)$ by summing the top $k+1$ residues; applying the same bound to the deficits controls the negative end. So for $c>4L/5$ this yields five positive indices of mean $c$ hitting any residue mod $M$. The unknown integer quotient is not controlled.

Here is a global obstruction at a fixed center. The sharp unconstrained real bound is
\[
10B^2\le Q^3 .
\]
For $Q>0$ maximize or minimize $B$ on the mean-zero sphere. Lagrange multipliers force at most two values, and the two-value calculations above give $B^2/Q^3=1/10$ or $1/135$. So the five-node cusp is the outer real envelope, while the ten-node cusp is interior. If $c\ge2$ and $a=c-1$, positivity requires $y_i\ge-a$ and gives $Q\le10a^2$ by convexity on the translated simplex. In addition
\[
 B\ge -aQ/3 + 2Q^2/(15a),
\]
by summing $(y_i+a)(y_i-t)^2\ge0$ and taking $t=2Q/(5a)$. In particular $\Delta=N-5F(c)=cQ+B\ge(2c+1)Q/3$. For real $H\ge1$ this leaves only $Q\le3H/(2c+1)$ for $\Delta\in[0,H]$, and hence at most
\[
 O\!\left(1+(H/c)^{5/2}\right)
\]
distinct integer $\Delta$ at center $c$, by summing $O(1+Q^{3/2})$ cubic defects per $Q$. For $H\longrightarrow\infty$ with $H=o(c^{5/3})$ this is $o(H)$. For example $5F(c)+1$ is locally solvable everywhere at that mean but is not exactly representable there for any $c\ge2$. This does not deny representations with other means. The local filter should be paired with an explicit real window rather than treated as a fixed-center Hasse principle for positive integers.

## 4. Exact selector for $(c,Q,B)$ using two squares

Choose one singleton offset $u$, two pair sums $s,t=-u-s$ and the corresponding differences $a,b$. After exchanging the pairs assume $s>t$. Put
\[
 A=4Q-2u^2-s^2-t^2,\qquad C=(24B-4u^3-s^3-t^3)/3 .
\]
Then
\[
 a^2=(C-tA)/(s-t),\qquad b^2=(sA-C)/(s-t).
\]
This comes from the fact that on each pair the squared difference enters the quadratic moment with a constant weight and the cubic moment with weight proportional to the pair sum. Require $C$ integral and both fractions to be nonnegative integral squares. Then require $a\equiv s,\ b\equiv t\bmod2$ and check each recovered offset $u,(s\pm a)/2,(t\pm b)/2$ is at least $1-c$.

This is exhaustive: unless all five entries agree, one can leave four entries not all equal after removing a singleton. Of the three pairings of four entries some have unequal sums, since equality for all three would force all four entries equal. The all-equal centered tuple solves only $(Q,B)=(0,0)$. Conversely every acceptance is an actual positive witness, not a modular test. The choice variables are bounded by $|u|\le\sqrt{2Q}$ and $|s|\le2\sqrt{2Q}$ and the positive simplex of side $5(c-1)$, giving $O(\min(c^2,Q+1))$ square tests per specified pair. For a scalar target iterate candidate $Q$ with $B=N-5F(c)-cQ$, applying the real and local guards before the scan. This is an exact selector at a chosen center, not a coverage proof. The replay checks all 22,498 realized pairs at centers 1 through 13 against sorted exhaustive enumeration, with additional perturbed-pair controls.

Replay: `python replay.py`.
