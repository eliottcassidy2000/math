# Six supplied triples: recurrence slots, golden signatures, and singular carriers

**Status:** PROVED elementary identities and a rational-scalar rigidity lemma;
CITED classical recurrence families and transforms; FINITE-EXACT replay of the
declared universes below. The two supplied noncanonical triples remain literal
objects throughout. Neither a classification of all integral recurrences nor
a new modular parametrization is claimed.

## 1. Inheritance and the exact parameter convention

The live objects are a recurrence with its initial jet, its differential
operator, the four singular locations, a golden power signature, and integer
roles in earlier packets. Their maps preserve different data.

The current proved starting point is
[THM-4711, Lucas-carry Wronskians](../../01-canon/theorems/THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md),
especially §7 and its six correctly named families. The nearby hostile is its
excluded D prime 5, where the singular polynomial degenerates. The corrected
near miss is [MISTAKE-638](../../01-canon/MISTAKES.md): factoring a truncated
polynomial does not identify its actual divided carry. The least-used sidecar
here is the accessory parameter $b$: singular-location data forget it.

Use the convention

$$
(n+1)^2u_{n+1}=(an(n+1)+b)u_n-cn^2u_{n-1},
\qquad u_0=1,\quad u_1=b.                         \tag{1}
$$

In [Zagier's original paper, Table 1 and §4](https://people.mpim-bonn.mpg.de/zagier/files/tex/AperylikeRecEqs/fulltext.pdf),
the columns are $(A,B,\lambda)$; these become $(a,c,b)$ here. His standard
six nondegenerate sporadic rows, also used in THM-4711, are:

| Name | $(a,b,c)$ | $a^2-4c$ | $u_0,u_1,u_2$ |
|---|---|---:|---|
| A, Franel | $(7,2,-8)$ | 81 | $1,2,10$ |
| B | $(9,3,27)$ | $-27$ | $1,3,9$ |
| C | $(10,3,9)$ | 64 | $1,3,15$ |
| D, Apéry $\zeta(2)$ | $(11,3,-1)$ | 125 | $1,3,19$ |
| E | $(12,4,32)$ | 16 | $1,4,20$ |
| F | $(17,6,72)$ | 1 | $1,6,42$ |

The supplied list contains B, D, C, E and two different objects:

$$
u_2(11,5,125)=\frac52,\qquad
u_2(17,5,72)=\frac{123}{4}.                       \tag{2}
$$

This follows directly from $4u_2=(2a+b)b-c$. In particular, putting 5 in F's
accessory slot destroys integrality at the first possible test. The occurrence
of $(11,5,125)$ has a genuine golden interpretation in §4, but with different
slot meanings. No input has been silently replaced.

The earlier [sporadic CM detector](lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md)
retains rational cusp slopes, modular frames, and CM restrictions. Those
objects are not inferred from the elementary carrier constructed next. Its
correction of a level-9 twist and the older label repair in
[the factorization bridge](factorization_vertex_supersymmetric_bridge_20261018.md)
are additional reasons to retain names and parameters together.

## 2. A useful carrier, and exactly what it loses

Writing $f(x)=\sum u_nx^n$ and $h(x)=1-ax+cx^2$, equation (1) becomes

$$L_bf:=(xh(x)f')'+(cx-b)f=0.                     \tag{3}$$

Assume characteristic zero and $c(a^2-4c)\ne0$. The four singular locations
are $0,\infty$ and the two roots of $h$. In the reciprocal coordinate they
are $0,\infty$ and the roots of $X^2-aX+c$. Their double cover is the smooth
auxiliary elliptic curve

$$E_{a,c}:Y^2=X(X^2-aX+c).                         \tag{4}$$

**This is the double cover of the singular-location set. It is not asserted
to be the elliptic fibre whose periods satisfy (3).** In particular, its
$j$-invariant is not a modularity certificate for (1).

The usual Weierstrass invariants, or direct substitution into the cubic
discriminant formula, give

$$
\Delta_E=16c^2(a^2-4c),\qquad c_4=16(a^2-3c),\qquad
J(a,c)=\frac{256(a^2-3c)^3}{c^2(a^2-4c)}.          \tag{5}
$$

Thus A, C, F give $1556068/81$, B gives 0, E gives 1728, and D gives
$488095744/125$. Distinct values prohibit a Möbius identification of the full
four-point singular sets over an algebraic closure. Equal values do not
restore the accessory, a rational twist, or a marked differential.

The harmonic and equianharmonic cases have short exact tests:

$$J=0\iff a^2=3c,\qquad
J=1728\iff a=0\ \text{or}\ 2a^2=9c.               \tag{6}$$

For the second test use
$4(a^2-3c)^3-27c^2(a^2-4c)=a^2(2a^2-9c)^2$.
The nondegeneracy assumption is essential in both equivalences.

The pair $(17,6,72)$ and $(17,5,72)$ is a minimal useful loss witness: the
same entire four-point set and same $J$ accompany different integrality.
The script's carrier receipt therefore contains the original $(a,b,c)$,
and its consumer rejects a valid F carrier presented for the altered source.
Merely comparing invariant fields would wrongly accept it.

## 3. Exact transformations and a rigidity boundary

For a nonzero scalar $s$, setting $v_n=s^nu_n$ transforms the parameters by

$$ (a,b,c)\longmapsto(sa,sb,s^2c).                 \tag{7}$$

The root discriminant is multiplied by $s^2$ and $J$ is unchanged. On the
carrier, $X=sx,\ Y=s\sqrt{s}\,y$ identifies the curves over a field containing
$\sqrt{s}$. For $s=r^2$ the map is rational: $(X,Y)=(r^2x,r^3y)$.
For $s=-1$, alternating coefficients retain a quadratic twist. Integer $s$
preserves coefficient integrality in this direction; division by $s^n$ needs
its own divisibility condition. The implemented scalar API uses nonzero
integers and retains the original source.

**Rational-scalar rigidity.** Fix $a,c$ and the independent variable. If a
nonzero rational function $R(x)$ multiplies the entire solution space of
$L_b$ into that of $L_{b'}$, then $R$ is constant and $b=b'$.

Indeed, with $p=xh$, direct differentiation gives

$$
L_{b'}(Rf)-R L_bf
=2pR'f'+\bigl(pR''+p'R'+(b-b')R\bigr)f.           \tag{8}
$$

At an ordinary point, solutions of $L_b$ have arbitrary independent values
of $(f,f')$. Hence both rational coefficients vanish. As $p\ne0$ in the
function field, $R'=0$, so characteristic zero gives constant $R$, followed
by $b=b'$. This concerns the whole two-dimensional solution space, not the
ratio of two chosen special solutions. It also excludes neither a change
of independent variable nor a higher-order differential intertwiner.

Consequently the wrong F accessory cannot be repaired by a rational scalar
gauge at fixed coordinates. This is stronger than merely observing (2),
but does not pretend to classify all transformations between equations.

For a concrete, source-retaining inverse of the initial coefficient map,
if $6u_2-8b^2\ne0$ then

$$
a=\frac{9u_3-17bu_2+4b^3}{6u_2-8b^2},\qquad
c=(2a+b)b-4u_2.                                  \tag{9}
$$

This follows by eliminating $c$ from the first two nontrivial recurrences.
It formalizes the initial-jet search already used by Zagier. The denominator
guard is real: with $b=c=0$, every $a$ gives the constant series 1. A
degenerate decoder must report missing information, not nonexistent data.

There is a stronger classical coordinate transformation: binomial
reanchoring. [Ji-Cai Liu, Theorem 2.5](https://arxiv.org/html/2406.18059v1)
gives its generic differential operator, with a cubic tail that disappears
at a singular anchor. The separately proved
[reanchoring kernel and carry packet](reanchoring_kernel_20261018c.md)
specializes it to A, C, F and keeps the divided-carry correction. In
particular $C_n=\sum_k\binom nk A_k$, and $C_4=639$. This is an actual
sequence map; the equality of this coefficient to a word in the prior
digit-tournament packet does not supply a map of those marked objects.

## 4. The golden signature explains one supplied triple

Let $F_m,L_m$ be the Fibonacci and Lucas integers, and put
$\phi=(1+\sqrt5)/2$, $\psi=(1-\sqrt5)/2$. The retained signed identity is

$$L_m^2-5F_m^2=4(-1)^m,\qquad
\phi^m,\psi^m=\frac{L_m\pm F_m\sqrt5}{2}.         \tag{10}$$

See the corrected source and sign history in
[difference families, §4](abc_difference_families_20261011c.md). For every
$m\ge1$, the polynomial

$$X^2-L_mX+(-1)^m$$

has these two slopes and discriminant $5F_m^2$. In particular,

$$ (L_5,F_5,5F_5^2)=(11,5,125),\qquad c=-1.       \tag{11}$$

So the supplied triple is exactly a meaningful signature
$(\text{slope sum},\text{square-root coefficient},\text{discriminant})$.
The inverse $c=(a^2-\Delta)/4$ recovers $-1$; it never supplies the accessory
$b=3$. Reading the same three integers as recurrence slots instead gives
the nonintegral row in (2), with discriminant $-379$ and
$J=4195088384/5921875$.

This all-power construction preserves the quadratic splitting field and
the signed norm. It does not preserve $J$: at $m=1$ the carrier has
$J=16384/5$, while at $m=5$ it has D's value in §2. Nor is taking powers of
the slopes an invertible reparametrization of (3). There is no claim that
arbitrary $m$ supplies an integral sporadic recurrence.

## 5. Bounded contextual atlas of every supplied integer

This is a route-level atlas of the thirteen distinct supplied integers,
not an enumeration of every occurrence of a small numeral in the repository.
The search used exact triples, named recurrence mechanisms, and the earlier
golden, digit, divisor, and elliptic packets. Each row distinguishes a
mathematical role from an unproved proposed identification.

| Integer | Exact role here | Recovered context and preserved boundary |
|---:|---|---|
| $-1$ | D's $c$; odd golden norm | [Golden difference families](abc_difference_families_20261011c.md): retain the sign. |
| 3 | B,C,D accessory | [CM detector](lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md): accessory, characteristic and level are different coordinates. |
| 4 | E accessory; E slope | [Divisor recursion](divisor_equation_recursions_20261011f.md): also a strict record input, not an operator transformation. |
| 5 | D field and $F_5$ | Equation (11) explains the supplied signature, not a replacement accessory. |
| 9 | B's $a$, C's $c$ and slope | [Digit C9](digit_idoneal_difference_20261018b.md): no graph map follows from a shared integer. |
| 10 | C's $a$ | [Ten-edge Paley cores](paley_next_cores_20261014i.md): an edge budget is not a recurrence parameter. |
| 11 | D's $a=L_5$ | Its splitting field is $\mathbb Q(\sqrt5)$, not the $\sqrt{-11}$ field in [the HN bridge](hn_elliptic_bridge_20261016c.md). |
| 12 | E's $a$ | Strict divisor record; the degree-12 marking in [the elliptic bridge](jacobian_elliptic_bridge_20261016c.md) retains a different source. |
| 17 | F's $a$ | [Fermat/order packet](fermat_wieferich_order_20261014e.md): primality does not determine $b$; canonical F has $b=6$. |
| 27 | B's $c=3^3$ | A coefficient and a prime-power field in [the wheel packet](wheel_apery_extension_20261014h.md) are separate typed inputs. |
| 32 | E's $c=2^5$ | Being a prime power does not make it a highly composite record. |
| 72 | F's $c=8\cdot9$ | The singular slopes 8,9 survive the failed $b=5$ perturbation. |
| 125 | D's discriminant $5^3$ | Putting this discriminant in the $c$ slot changes the curve and recurrence. |

For an exact comparison with the recovered highly-composite thread, the
strict divisor-count records through 125 are
$1,2,4,6,12,24,36,48,60,120$. Among the supplied positive integers only
4 and 12 are record **inputs**. The record **values** intersect the supplied
set in $3,4,9,10,12$. These are different predicates; for instance
$\tau(72)=\tau(60)=12$, so 72 is not a new record. This classification is
finite and complete to the stated cutoff.

## 6. Reproduction and the remaining obligation

The [standard-library checker](../../04-computation/experiments/zagier_six_triples_20261018c.py)
and [saved JSON](zagier_six_triples_20261018c.json) retain the exact source and
rational coefficients. The replay covers all 729 triples in $[-3,5]^3$
through index 8, all six canonical prefixes through 24, four independent
classical binomial formulas through 16, golden powers 1 through 24, and
divisor records through 125. The operator is independently expanded; the
gauge commutator is checked after full denominator clearing. Tests include
invalid types, incorrect source labels, smoothness failures, a degenerate
jet, and both literal noncanonical triples. Assertions remain active under
optimization.

Normal and independently invoked optimized replays agree exactly with the
saved JSON: **9,938 checks**; LF SHA256
`7c6e484389961f8705e09bbf941018a7ec2fab1e0a6a8018bb5ba3e170c1dae2`.
An independent peer read the full proof and production interfaces and replayed
the optimized checker with the same 9,938 checks: **PASS**, no correction.

```powershell
python -B -X utf8 04-computation/experiments/zagier_six_triples_20261018c.py --out 05-knowledge/results/zagier_six_triples_20261018c.json
python -B -O -X utf8 04-computation/experiments/zagier_six_triples_20261018c.py
```

The useful next question is marked: which changes of variable together with
an explicit differential or carry transformation preserve a selected
arithmetic recurrence? Matching a discriminant, splitting field, integer
alphabet, or auxiliary $J$ is insufficient. The adjacent reanchoring packet
answers one such transformation problem for A/C/F; general accessory or
modular equivalence remains a separate obligation.
