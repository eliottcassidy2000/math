# Quadratic cycles and exact polynomial-observable return budgets

2026-10-12. **PROVED, scoped:** the critical-value filtration, repeated-value observable ring, and its sharp polynomial-return compiler below. **INHERITED:** the local observable argument and separated rational Keller construction. **FINITE-EXACT:** the declared controls in the companion script. **OPEN:** the remaining planar polynomial-source problem and JC(2). A rational symplectic map of degree greater than one is not a polynomial Keller map or a birational source change.

## 1. Inheritance and the actual next interface

The current [time-form divisor note](timeform_ramification_20261011.md) already recovers the all-degree obstruction for its genus-two carrier; [source frames](source_frames_eit_20261011.md) and [Möbius source returns](mobius_affine_source_return_20261011.md) close their specified one-flow and changing-generator cases. The latter retains an affine offset and does construct an elliptic involution. Its fixed-invariant classification is not an invitation to repeat a polynomial ansatz in the same closed carrier.

The closest exact mechanism here is [THM-3567, separated rational Keller maximal observables](../../01-canon/theorems/THM-3567-separated-rational-keller-maximal-observable-nodal-completion-no-go.md). Its stated hypotheses include distinct critical values at different critical points. Iteration violates that hypothesis dramatically. We extend its local proof to repeated critical values and apply it to the two specified quadratic maps. The new outcome is an all-depth filtration with an exact cancellation budget, rather than a new claim about arbitrary planar sources.

The canonical hostile is a unit-Jacobian rational map with two source sheets and a removed ramification divisor. The corrected near miss is interpreting a periodic orbit, or even a whole returning cotangent fibre, as polynomial return on the complete affine plane. The least-used sidecar is the *critical-value divisor with its local ramification orders*. The concept board is: **marked source / cotangent lift / postcritical orbit / local ramification / polynomial observables / finite-depth return**.

The two dynamics are already distinguished in [difference visibility](difference_visibility_norms_20261004.md), [THM-4139, rational three-cycle lift](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md), and [THM-4146, horizontal divisor/fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md). No new arithmetic-dynamics literature theorem is imported here.

## 2. Cotangent return keeps a multiplier and a source sheet

Let $f_c(x)=x^2+c$. Its separated cotangent lift is

$$\widehat f_c(x,p)=\left(f_c(x),\frac{p}{2x}\right),\qquad x\ne0.$$

It preserves $dx\wedge dp$. With $F_n=f_c^{\circ n}$, the actual iterate is

$$\widehat f_c^{\circ n}(x,p)=\left(F_n(x),\frac{p}{F_n'(x)}\right),
\qquad F_n'(x)\ne0. \tag{1}$$

The guard is precisely avoidance of zero at every preceding base step. It cannot be removed by noting that some chosen observable has a polynomial extension.

At a base cycle of period $r$ with multiplier $\mu=\prod_{i=0}^{r-1}2x_i\ne0$, the returning momentum is $p/\mu$. Therefore a specified cotangent point returns iff $p=0$ or $\mu=1$. This proves an actual connection between the scalar multiplier and symplectic return; the secant product $\prod_i(x_i+x_{i+1})=1$ is a different invariant and does not supply this test.

For $c=-7/4$, let

$$P(x)=x^3+\tfrac12x^2-\tfrac94x-\tfrac18.$$

The inherited factorization, checked exactly again, is

$$F_3(x)-x=(f_c(x)-x)P(x)^2.$$

Here $\operatorname{disc}(P)=49$, $\gcd(P,f_c-x)=1$, and $F_3'=1$ modulo $P$. Thus its three distinct real roots form an exact 3-cycle and every momentum returns over each of them. For $c=-29/16$, the rational cycle

$$-7/4\longmapsto5/4\longmapsto-1/4\longmapsto-7/4$$

has $\mu=35/8$, so only momentum zero returns, while every momentum is rescaled by $8/35$.

Neither example is a birational map of the plane. The function-field degree of (1) is $2^n$: adjoining $x$ to $k(F_n,p/F_n')$ has the same degree as adjoining $x$ to $k(F_n)$, and then $p=(p/F_n')F_n'(x)$. Already $(1,2)$ and $(-1,-2)$ are distinct legal sources with the same image $(-3/4,1)$ for $c=-7/4$. A marked inverse branch is missing from an unmarked endpoint.

The number 63 has a separate exact critical-orbit role:

$$0\mapsto-7/4\mapsto21/16\mapsto-7/256,
\qquad (-7)(-7+4)^2+4^3=64-63=1.$$

This is the inherited critical numerator/norm computation, not the algebraic 3-cycle starting at a root of $P$.

## 3. Exponentially many critical points, linearly many retained values

Work over $\mathbf Q$, or its algebraic closure for local orders, and assume $v_2(c)<0$. Set $v_i=f_c^{\circ i}(0)$ and

$$\Delta_n(T)=\prod_{i=1}^n(T-v_i). \tag{2}$$

**Proposition.** The values $v_i$ are nonzero and pairwise distinct. Every critical point of $F_n$ is simple; the critical points mapping to $v_i$ number $2^{n-i}$. Consequently $F_n$ has $2^n-1$ critical points but exactly $n$ critical values.

**Proof.** Put $b=v_2(c)<0$. The two terms in $v_i=v_{i-1}^2+c$ have distinct valuations after the first step, hence

$$v_2(v_i)=2^{i-1}b.$$

This proves nonzero/distinct values and, in particular, proves that the critical orbit is not preperiodic; no numerical orbit evidence is needed. The chain rule gives

$$F_n'(x)=2^n\prod_{j=0}^{n-1}F_j(x),\qquad F_0(x)=x. \tag{3}$$

Each $F_j$ has simple roots: a repeated root would encounter zero at an earlier step and return to zero, contrary to the preceding valuation argument. The same argument makes these root sets pairwise disjoint. A root of $F_{n-i}$ maps under $F_n$ to $v_i$, giving the count and simple ramification index two. QED.

The fibre $F_n=v_i$ contains $2^{n-i}$ double roots and $2^n-2^{n-i+1}$ simple roots. In particular, remembering only the value loses the number and location of the source sheets. Repetition of a critical *value* does not make a critical *point* nonsimple.

## 4. Exact all-degree polynomial-return compiler

Keep the common **target** polynomial ring $k[X,Y]$ fixed and define

$$\mathcal A_n=\{H(X,Y)\in k[X,Y]:H(F_n(x),p/F_n'(x))\in k[x,p]\}.$$

This is polynomial regularity of an observable after rational pullback. It is neither equality $H\circ\widehat f^n=H$ nor regularity of the full map (1).

Write $H=\sum_{j\ge0}a_j(X)Y^j$. Then

$$H\in\mathcal A_n\quad\Longleftrightarrow\quad
\Delta_n(X)^{\lceil j/2\rceil}\mid a_j(X)\quad\text{for every }j. \tag{4}$$

**Proof.** Distinct powers of $p$ cannot cancel. At a simple critical point $\alpha$ with $F_n(\alpha)=v_i$, the local orders of $F_n-v_i$ and $F_n'$ are two and one. Thus the coefficient $a_j(F_n)/(F_n')^j$ is regular there iff $2\operatorname{ord}_{v_i}a_j\ge j$. All critical points over the same value impose the same condition. There are no other finite poles. This proves both directions even though different critical points share a value. QED.

In particular, with target generators $A=\Delta_n(X)Y$ and $B=\Delta_n(X)Y^2$,

$$\mathcal A_n=k[X,A,B]\cong k[X,A,B]/(A^2-\Delta_n(X)B). \tag{5}$$

For even $j=2m$ the minimum monomial is $B^m$; for odd $j=2m+1$ it is $AB^m$. Algebraic independence of $X,Y$ shows that the displayed relation is the full kernel. The surface has exactly $n$ ordinary double points $(v_i,0,0)$: the gradient vanishes precisely there, and its local quadratic part is $A^2-\Delta_n'(v_i)(X-v_i)B$, nondegenerate.

The **full-field** extension of inherited THM-3567 also holds. If a rational target $H\in k(X,Y)$ pulls back polynomially, first regard it as a reduced rational function in $Y$ over $k(X)$. Scalar extension to $k(x)$ and the invertible linear rescaling of $p$ force it to be polynomial in $Y$. A finite pole of any coefficient at $X=z$ remains a pole at any root of $F_n(x)=z$; division by $F_n'$ cannot cancel it. Hence those coefficients are polynomials, and (4) applies. Therefore the image of (5) is exactly

$$k(F_n,p/F_n')\cap k[x,p].$$

**Variable-identification boundary.** The filtration below is in the fixed target ring $k[X,Y]$. Its images live in different embedded source subfields as $n$ changes. We do not assert those source intersections are nested.

**Sharp recursion and budget.** Since $\Delta_{n+1}=(X-v_{n+1})\Delta_n$,

$$\mathcal A_{n+1}\subsetneq\mathcal A_n,\qquad
A_{n+1}=(X-v_{n+1})A_n,\quad B_{n+1}=(X-v_{n+1})B_n.$$

A nonzero momentum-$j$ coefficient owes at least $n\lceil j/2\rceil$ base degree. The observable $Y\Delta_N(X)$ returns polynomially at every depth $1\le n\le N$ and fails at depth $N+1$. Thus arbitrarily deep finite regularity exists but never implies all-depth regularity. Indeed,

$$\bigcap_{n\ge1}\mathcal A_n=k[X]. \tag{6}$$

Each nonzero coefficient with $j>0$ would otherwise have infinitely many distinct roots. This is a stronger, exact statement than failure of a bounded search, but its parameter and source scope remain essential.

The new positive operation is `compile_return`: it authenticates (4) and returns the unique polynomial extension on the original source. Parameter, depth, target polynomial and extension are kept in a checked receipt. The program separately divides the substituted numerator by $(F_n')^j$, rather than treating the critical-value test alone as the computed output. Point evaluation of the original rational lift still rejects $F_n'=0$.

## 5. What cancellation buys, and what it does not

The minimal linear-momentum completion is

$$P=F_n(x),\qquad Q=\frac{\Delta_n(F_n(x))}{F_n'(x)}p.$$

It is polynomial, but its Jacobian is $\Delta_n(F_n(x))$, not one. The same critical-value divisor paying the denominator now appears in the area multiplier. A polynomial observer is useful information; changing the observer into a coordinate has a retained ramification bill.

**Stronger all-observable obstruction.** No two members of the full-field polynomial-observable ring can have nonzero constant source Jacobian. This extends the inherited nodal-completion mechanism to the repeated-critical-value setting, rather than testing only the displayed minimal pair.

Indeed write $h(x)=\Delta_n(F_n(x))/(F_n'(x))^2$, a polynomial coprime to $F_n'$. The three source generators are

$$P=F_n(x),\qquad A=F_n'(x)h(x)p,\qquad B=h(x)p^2.$$

At each actual source point $(\alpha,0)$ with $F_n'(\alpha)=0$, all three differentials vanish. By the chain rule the differential of *every* polynomial in these generators vanishes there. Section4 proved that these are every rational target observable whose pullback is polynomial. Therefore every pair of such observables has zero source Jacobian at $(\alpha,0)$; a nonzero constant is impossible. Although the rational map itself is undefined at these points, its proposed polynomial extensions are defined there and must pass the standard-area test there. Deleting these points changes the source problem.

**General cotangent corollary, without squarefree ramification.** Let $F\in k[x]$ have degree at least two in characteristic zero. Then any two elements of $k(F,p/F')\cap k[x,p]$ have zero Jacobian at $(\alpha,0)$ for every critical point $\alpha$ of $F$. Consequently no rational target reencoding of this cotangent lift into two source-polynomial observables is a polynomial Keller pair. This is an elementary extension of the recovered separated-carrier obstruction, not a statement about arbitrary rational symplectic maps or arbitrary planar sources.

To prove it, the full-field pole argument from Section4 still writes each observable uniquely as $\sum_j a_j(F)(p/F')^j$ with $a_j\in k[X]$. Put $e=\operatorname{ord}_\alpha(F-F(\alpha))\ge2$. For $j=0$ the source $x$ derivative contains $F'$. The $p$ derivative at $p=0$ can come only from $j=1$. Polynomiality forces $e\operatorname{ord}_{F(\alpha)}a_1\ge e-1$, hence $\operatorname{ord}_{F(\alpha)}a_1\ge1$. The resulting coefficient $a_1(F)/F'$ therefore vanishes at $\alpha$, with order at least $e-(e-1)=1$. Terms $j\ge2$ have both derivatives zero at $p=0$. Thus every observable differential vanishes there. Critical points exist over the algebraic closure; that suffices to exclude a nonzero constant Jacobian over $k$. The degree-one boundary is real: $(2x,p/2)$ is a polynomial symplectic map.

Three sharp controls distinguish the hypotheses.

1. At $c=0$, $F_2=x^4$ has ramification index four. The squarefree-critical-point rule would incorrectly accept $X^2Y^3$, whose pullback is $p^3/(64x)$. In general a point of ramification index $e$ requires $\operatorname{ord}_{F(\alpha)}a_j\ge\lceil j(e-1)/e\rceil$, and points sharing a value impose the maximum of these bounds.
2. At $c=-2$, the postcritical set is the finite set $\{-2,2\}$. The nonzero momentum observable $(X^2-4)Y^2$ returns polynomially at every depth. Equation (6) requires the proved infinite critical orbit, not merely a quadratic map.
3. Polynomial return is not a polynomial invariant. For any polynomial $f$ of degree $d\ge2$, a polynomial $H=\sum a_j(x)p^j$ invariant under $(f,p/f')$ must be constant. For $a_j\ne0$, coefficient equality gives $a_j(f)=(f')^ja_j$. Comparing degrees forces $\deg a_j=j$, and comparing leading coefficients then forces $1=d^j$, impossible when $j>0$. The $j=0$ coefficient is constant by degree. This does not contradict the returning fibres in Section2 or finite-depth compiler in Section4.

| Source | Target/map and preserved predicate | Lost coordinate / mandatory sidecar | Decisive control |
|---|---|---|---|
| Marked quadratic cycle | Cotangent cycle via (1) | Derivative multiplier, not only secant product | $1$ versus $35/8$ |
| Rational cotangent map | Polynomial observable via (4) | Full-map poles, area multiplier and source sheet | $(1,2)$ and $(-1,-2)$ collide |
| Critical points | Distinct-value divisor $\Delta_n$ | Sheet multiplicities and local index | $2^n-1$ points but $n$ values |
| Depth-$N$ regular observable | Deeper polynomial-return test | New critical-value factor | $Y\Delta_N$ fails exactly at $N+1$ |
| Completed polynomial coordinates | Standard-area Jacobian test | Multiplier $\Delta_n(F_n)$ | Polynomiality does not make it constant |

**Next falsifiable obligation, rather than another closed-carrier ansatz.** A proposed changing-source construction should report its exact pullback of the standard area form, every source sheet and deleted divisor, and the regularity of its full coordinates—not just an observable or a selected cycle. Within this quadratic cotangent family the generic covering degree and completion multiplier already settle those tests negatively. The remaining planar return problem must change that mechanism; this packet supplies an exact filter and positive observable compiler, not a candidate JC(2) solution.

## 6. Reproduction

[Exact script](../../04-computation/experiments/quadratic_clock_return_20261012a.py), [saved output](quadratic_clock_return_20261012a.out):

```text
python -B -X utf8 04-computation/experiments/quadratic_clock_return_20261012a.py
python -B -O -X utf8 04-computation/experiments/quadratic_clock_return_20261012a.py
```

The declared universe is $c\in\{-7/4,-29/16,3/2\}$ at depths1–4 for structural identities and 48 positive receipts (momentum degree0–4 at depths1–2, degree0–2 at depths3–4); two parameters at depths1–3 for 288 independent bounded coefficient controls; eight direct simultaneous-substitution controls; separate higher-ramification, finite-postcritical, cycle, collision and thirteen exact-type/forged-receipt hostiles. These finite controls support the formulas; the arbitrary-depth and arbitrary-degree conclusions are proved above.

The script also verifies the six generator derivatives and a nonlinear target reencoding at every declared parameter/depth, and the general local-index formula for $F=x^e$, $2\le e\le6$, momentum degrees1–6. Normal and optimized runs agree with saved output at **1,017 always-active exact checks**, LF SHA256 `fb42d254cc679839239dd41078eb39693b3926f9ba6ab181243c86f5715af08f`.
