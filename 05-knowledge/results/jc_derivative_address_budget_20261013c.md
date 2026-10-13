# Finite derivative frontiers and source-degree budgets

2026-10-13. **PROVED:** the scoped frontier, strand-growth and weighted-budget
criteria below. **FINITE-EXACT:** the declared source-polynomial controls.
**OPEN:** finding a successful criterion for every planar Keller map. This
does not prove JC2 or use finite experiments as a universal closure bound.

## 1. Inheritance and the change of object

The closest mechanism is the
[finite dual-field inverse compiler](jc_bilinear_reconstruction_20261013a.md):
an exactly closed finite-dimensional space of actual source polynomials,
containing the coordinates, yields a polynomial inverse. The current
[THM-4626, generic-line and local-nilpotence equivalence](../../01-canon/theorems/THM-4626-surjective-dual-fields-are-locally-nilpotent-and-hyp-9255-is-the-planar-jacobian-conjecture.md)
already removes the supposed nonzero torsion-cokernel escape. We do not
reintroduce that route or repeat its proof.

The incoming [transverse jet identity](keller_transverse_jet_length_20261013.md)
shows why thickening the same ordinary global cokernel cannot reduce its
obstruction: its length is exactly $b\mu$, and a fixed-component mate shear
preserves the connection itself. That carrier is different from the
polynomial derivative addresses used here. Its punctured-line and Pruefer
controls are the corrected near misses: local or quotient success need not
preserve the global source. The canonical hostile used below is the actual
three-variable Keller map of
[THM-4621, the THM-1300 coordinate-pair fibrations](../../01-canon/theorems/THM-4621-thm1300-coordinate-pair-fibrations-are-rational-functions-minus-critical-points.md),
whose duals are not locally nilpotent. It is not a planar counterexample.

Our least-used sidecar is the **minimal excluded derivative address**, together
with a literal polynomial zero at that address. The board is: well-founded
addresses, source-polynomial identities, degree capacity, constant coefficient
recurrences, and evaluated jets. The anchor is an inverse certificate; the
niche is a necessary degree-growth test for nonclosure; the wildcard is how
a seemingly slow recursion can hide its first nonzero point-jet coefficient.
No identity with the Conway or Hofstadter recurrence is asserted.

## 2. A finite address boundary suffices

Let $k$ have characteristic zero, let $F:k^n\to k^n$ be a polynomial map
with nonzero constant Jacobian, and let $D_1,\ldots,D_n$ be its polynomial
dual derivations. Thus $D_jF_i=\delta_{ij}$ and all the $D_j$ commute.
These are the duals of the **complete supplied source map**.

Let $S\subset\mathbb N^n$ be a nonempty finite downward-closed set. Write
$\partial S$ for the minimal elements of its complement under coordinatewise
order. They are computed by taking one-step neighbours outside $S$ and
keeping those whose immediate predecessors all lie in $S$.

**Frontier criterion.** If

$$D^\alpha x_i=0\quad\text{as a source polynomial for every }
  \alpha\in\partial S\text{ and every }i,                 \tag{1}$$

then $F$ is a polynomial automorphism. At any retained source point $a$ its
inverse is the finite expression

$$G_i(T)=\sum_{\alpha\in S}
 \frac{(D^\alpha x_i)(a)}{\alpha!}(T-F(a))^\alpha.          \tag{2}$$

In particular its total degree is at most $\max_{\alpha\in S}|\alpha|$;
its degree in $T_j$ is at most $\max_{\alpha\in S}\alpha_j$.

**Proof.** Every excluded address dominates a minimal excluded one. Since
the dual derivations commute, (1) kills all addresses outside $S$. Hence each
$D_j$ is locally nilpotent on every coordinate, and Leibniz gives local
nilpotence on the whole polynomial ring. The commuting exponentials define
a polynomial additive action taking $F$ to $F+t$. Evaluate that action at
$a$ to obtain (2) and $F\circ G=\mathrm{id}$. The coordinates of $G$ are
algebraically independent because their algebra contains the target
coordinates through $F\circ G$. Its image is therefore dense, and
$G\circ F\circ G=G$ implies $G\circ F=\mathrm{id}$. This is the inherited
inverse argument with a smaller, explicitly marked stopping interface.

Conversely, if $G=F^{-1}$ is polynomial, then

$$D^\alpha x_i=(\partial^\alpha G_i)(F).$$

Thus the set of addresses where at least one coordinate derivative is a
nonzero polynomial is a finite downward-closed set. It is exactly the
downward closure of the monomial supports of the inverse coordinates.
Its largest total address is exactly $\deg G$. The frontier criterion is
therefore complete as a semidecision system, without any asserted uniform
search cap. A supplied larger set is also sufficient; minimality is optional.

A rectangle gives the convenient sufficient test
$D_j^{m_j}x_i=0$ for every $i,j$. Its inverse degree bound is
$\sum_j(m_j-1)$, **not** $\max_j(m_j-1)$. A simplex retains a total-order
bound; arbitrary staircases retain anisotropic support. The implementation
checks downward closure rather than assuming independent coordinate caps
are a total-degree cap.

For $F=(x-y^2,y)$ the minimal set is
$(0,0),(0,1),(0,2),(1,0)$ and its excluded boundary is
$(0,3),(1,1),(2,0)$. The inverse is $(X+Y^2,Y)$. A nonminimal $2\times3$
rectangle still works, but does not give the sharp support record.

## 3. Nonterminating dual strands must occupy new polynomial dimensions

The no-nonzero-eigenvalue lemma for characteristic-zero Keller duals is
inherited from THM-4620 and the previous inverse compiler. Briefly, if
$Dh=\lambda h$ for nonzero $\lambda$ over an algebraic closure, differentiation
of the minimal polynomial of $h$ over $k(F)$ makes each lower coefficient
satisfy a rational equation $a'=c a$ with nonzero $c$. Such a rational
function is zero, by its poles and degree. The resulting polynomial $T^m$
cannot be the minimal polynomial of nonzero $h$.

**Strand theorem.** For one dual $D$ and a polynomial $h$, the nonzero
initial terms

$$h,Dh,\ldots,D^r h$$

are linearly independent over $k$. More precisely, they are independent
whenever $D^r h\ne0$. If a constant-coefficient recurrence first occurs,
the resulting finite cyclic space is $D$-stable. Its only possible eigenvalue
is zero, so its minimal polynomial is $T^m$. The first recurrence is therefore
termination, not a persistent oscillation or a nonzero geometric mode.

Let $d_r=\max_{0\le j\le r}\deg(D^j h)$. Nontermination through $r$ forces

$$\binom{n+d_r}{n}\ge r+1.                               \tag{3}$$

This is simply the dimension of the full degree-$d_r$ polynomial space.
In the plane, $d_r$ must be at least the smallest integer $d$ with
$(d+1)(d+2)/2\ge r+1$. A uniform bound $d_r=o(r^{1/n})$ would therefore
force termination. Such a bound for each coordinate and each dual would
produce a finite rectangle and hence an inverse by section 2. **No such
bound for arbitrary planar Keller inputs is proved here.**

One finite formulation is useful: if all iterates through
$r=\binom{n+B}{n}$ have degree at most $B$, then the last one must be zero.
It is an exact capacity test, not a claim that computing all those degrees
is faster than testing the final derivative directly.

The coefficient field is essential. Every nonzero strand is trivially
one-dimensional over $k(x)$, where $Dh=(Dh/h)h$. For the inherited THM-1300
map the first five $D_3$-iterates of $x$ have constant-field rank five,
although this rational-function recurrence holds. Its coefficient is not
constant. Likewise the Euler derivation $x\partial_x$, **outside the Keller
dual category**, has bounded degree and a nonzero eigenvalue; it violates
the conclusion if that hypothesis is dropped.

## 4. A directly checkable weighted budget

Choose positive integer source weights $w_1,\ldots,w_n$. Suppose every
monomial $x^\beta$ with nonzero coefficient in every $D_jx_i$ satisfies

$$\beta\cdot w\le w_i.                                  \tag{4}$$

Then each dual preserves weighted degree. Put $L=\max_i w_i$ and

$$M=\#\{\beta\in\mathbb N^n:\beta\cdot w\le L\}.$$

The complete weighted polynomial space of degree at most $L$ has dimension
$M$, contains the unit and all source coordinates, and is exactly stable
under every dual. On this space the commuting dual matrices are nilpotent.
They have a common strict triangular flag, so every mixed product of length
$M$ vanishes. Consequently the inverse has total degree at most $M-1$ and
the frontier search is guaranteed to close by that order.

This gate permits **nonincreasing**, not necessarily strictly decreasing,
weights. Nilpotence on each equal-weight layer comes from the actual Keller
no-eigenvalue lemma. Inequalities (4) alone are linear rational inequalities
in the weights; the supplied-weight consumer checks them exactly. It is not
an implementation of a complete weight-feasibility solver.

For every integer $m\ge2$, the triangular family $F=(x-y^m,y)$ admits
$w=(m-1,1)$. Here $M=m+1$, the inverse degree is exactly $m$, and
$D_2^m x=m!$, $D_2^{m+1}x=0$. Thus the general bound is attained. The code
also verifies the symbolic invoice $m=1024$ without generating a long
derivative sequence. These familiar triangular automorphisms are controls,
not a newly solved Jacobian family.

The gate is not necessary. The automorphism
$F=(x+y^2,y+(x+y^2)^2)$ has an $xy$ monomial in $D_1x$; (4) would require
$w_x+w_y\le w_x$, impossible for positive weights. Its exact frontier still
closes. Failure of a supplied weight, of all positive weights in this
coordinate presentation, or of a capped frontier search is not a
noninvertibility certificate.

## 5. Evaluation, well-founded addresses and the next obligation

For $F=(x-y^m,y)$ at source $(0,0)$, the evaluated $D_2^r x$ is zero for
$1\le r<m$, while its actual polynomial is nonzero. The first derivative
$my^{m-1}$ already exposes the obstruction to claiming premature closure.
At $m=1024$, every evaluated value on the proposed degree-one simplex
boundary is zero, so even checking the whole wrong boundary at that point
would accept it. The literal polynomial checker rejects the same packet.
Arbitrarily many vanishing **point values** do not replace one boundary
polynomial identity. This is the precise address-recursion lesson: both
the predecessors and their values in the original carrier matter.

| Source | Preserved predicate | Missing coordinate / hostile |
|---|---|---|
| Finite address ideal | All descendants of a zero boundary remain zero | Commuting actual duals and complete boundary |
| Full polynomial derivative | Literal equality on the source | Evaluation at one point can hide arbitrarily late support |
| Constant-field strand | Independent new direction until termination | Rank over the source function field collapses to one |
| Positive source weights | Stable finite polynomial capacity | A failed weight gate can still have a successful exact frontier |
| Ordinary transverse jets | Their exact linear length | They are a different global module, not this derivative address set |

The concrete next obligation is to prove (4), a comparably effective
filtration after a **certified source coordinate change**, or the subcritical
strand-growth bound for a specified unresolved Keller class. The verifier
and inverse compiler now accept those statements as exact finite inputs.
No recursion analogy, scalar slow-increment pattern, or quotient length
automatically supplies them.

## 6. Reproduction

[Script](../../04-computation/experiments/jc_derivative_address_budget_20261013c.py)
and [saved output](jc_derivative_address_budget_20261013c.out):

```text
python -B -X utf8 04-computation/experiments/jc_derivative_address_budget_20261013c.py
python -B -O -X utf8 04-computation/experiments/jc_derivative_address_budget_20261013c.py
```

Declared controls: four ordered three-shear maps at two anchors; exact
boundary and two-sided inverse identities; nine integer points per receipt;
five sharp triangular weight/strand budgets; the degree-1024 delayed-jet
control; the actual inherited three-variable nonnilpotent strand; and exact
type, missing-boundary, wrong-source and cap hostiles. All checks are active
under `-O`. Normal and optimized runs both match the saved output at
**189 exact checks** (about ten seconds each). LF SHA256:
`ce52fafb43ed3259c15f71b0566e34f0b9e800d2c05702000e5b4618e955bf16`.
