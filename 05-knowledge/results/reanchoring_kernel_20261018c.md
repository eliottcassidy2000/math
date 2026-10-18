# Reanchoring a recurrence without losing its carry or group action

Status: **PROVED**, the scoped identities and action obstruction below;
**CITED**, the classical binomial-transform mechanism; **FINITE-EXACT**, the
declared tests. This is a transfer tool, not a Collatz or plane-colouring proof.

## 1. Inheritance and the live comparison

The closest proved mechanism is [THM-4711, Legendre Wronskian and binomial-cube
Lucas carry divisibility](../../01-canon/theorems/THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md),
especially its six-family first-block carry identity. Its hostile is loss of
the carry polynomial when only the Lucas polynomial is retained; its corrected
near miss is the characteristic-p mistake that a polynomial with zero
derivative must be constant. The least-used sidecar here is the actual
coordinate-change law for that carry.

The anchor is the literal six-colour boundary-support problem. The niche is
Zagier recurrence reanchoring; the wildcard is the difference between the
two six-element S4 sets. The concept board is: pinned source; actual incidence;
relative six-core; residual branch phases; recurrence accessory parameter;
first-block carry. A connection is useful only if its map preserves its
consumer's predicate.

For geometric status use [the coupled fractional packet](coupled_fractional_geometry_20261018c.md)
and [the joint branch compiler](coupled_support_branches_20261018c.md).
For the supplied number/triple atlas use [the Zagier packet](zagier_six_triples_20261018c.md).
Shared small integers or splitting fields are prompts for these maps, not
substitutes for them.

## 2. Classical transform, with the root condition made explicit

Use the convention

\[
 (n+1)^2u_{n+1}-(an^2+an+b)u_n+cn^2u_{n-1}=0,
 \qquad u_0=1,\quad u_1=b.
\]

For \(f(x)=\sum u_nx^n\), put \(h(x)=1-ax+cx^2\). Its equation is

\[
 (xh f')'+(cx-b)f=0.
\]

In the reciprocal coordinate \(t=1/x\), with
\(Y(t)=f(1/t)/t\), this becomes the self-adjoint equation

\[
 (P(t)Y')'+(t-b)Y=0,\qquad P(t)=t(t^2-at+c).
\tag{1}
\]

Let \(d\) be a root of \(P\), and let \(\epsilon\ne0\). Substitute
\(t=d+\epsilon z\), \(Z(z)=\epsilon Y(d+\epsilon z)\). Then (1) has exactly
the same form, with

\[
 \boxed{a'=(a-3d)/\epsilon,\quad b'=(b-d)/\epsilon,
 \quad c'=P'(d)/\epsilon^2.}
\tag{2}
\]

Indeed \(P(d+\epsilon z)/\epsilon^3
=z^3-a'z^2+c'z\). If \(P(d)\ne0\), an additional constant
\(P(d)/\epsilon^3\) remains in the leading polynomial. Thus the root guard
is necessary and sufficient for this change of the **whole operator** to
remain in this family. This is not a claim that an exceptional single solution
cannot satisfy some additional lower-order relation.

The corresponding normalized sequence is

\[
 u'_n=\epsilon^{-n}\sum_{k=0}^n\binom nk(-d)^{n-k}u_k,
\quad
 g(x)=\frac1{1+(d/\epsilon)x}
 f\!\left(\frac{x/\epsilon}{1+(d/\epsilon)x}\right).
\tag{3}
\]

This binomial transform is classical, not a new theorem of this packet.
Zagier gives the sporadic table and related transforms in
[his original paper](https://people.mpim-bonn.mpg.de/zagier/files/tex/AperylikeRecEqs/fulltext.pdf).
Ji-Cai Liu's [arXiv:2406.18059v1](https://arxiv.org/html/2406.18059v1), Lemma 2.4
and Theorem 2.5, gives the general generating-function/operator transform;
its cubic tail vanishes precisely at the same root condition. Equations
(1)--(3) give a direct derivation in the present convention.

For Franel A, \((a,b,c)=(7,2,-8)\), the three roots are \(0,-1,8\).
Choosing \((d,\epsilon)=(-1,1)\) gives C, and choosing \((8,-1)\) gives F:

\[
 C_n=\sum_k\binom nk A_k,\qquad
 F_n=\sum_k\binom nk8^{n-k}(-1)^kA_k.
\]

Their parameters are \((10,3,9)\) and \((17,6,72)\), respectively.
In particular \(C_4=639\), an exact occurrence of one of the earlier
nine digit addresses. That single coincidence adds no further invariant.
The six choices of root and orientation give three signed pairs of
presentations. They form a coordinate-change groupoid; do not call them
six isometries of the original unevenly spaced root set.

The literal supplied variants stay separate:
\((17,5,72)\) gives \(u_2=123/4\), and \((11,5,125)\) gives \(u_2=5/2\).
The latter also has a useful, different type as the golden discriminant
signature \((L_5,F_5,5F_5^2)\). Neither replacement silently preserves the
recurrence's accessory parameter.

## 3. The first carry has an indispensable correction

Here is a general integer identity; it does not require a Zagier recurrence.
Let \(u_n\) be integers, let \(b,s\in\mathbb Z\), let \(p\) be an odd prime,
and assume just the first-block Lucas relation

\[
 u_{p+j}\equiv b u_j\pmod p\qquad(0\le j<p).
\]

Define

\[
 v_n=\sum_{k=0}^n\binom nk s^{n-k}u_k,
 \quad V_u(j)=\frac{u_{p+j}-bu_j}{p},
 \quad V_v(r)=\frac{v_{p+r}-(b+s)v_r}{p}.
\]

The last quotient is integral, and for \(0\le r<p\),

\[
 \boxed{V_v(r)\equiv
 \sum_{j=0}^r\binom rj s^{r-j}V_u(j)
 +\mathcal C_{p,r}(u,s)\pmod p,}
\tag{4}
\]

where the following numerator is divisible by \(p\):

\[
\begin{split}
 p\mathcal C_{p,r}(u,s)={}&
 \sum_{k=0}^{p-1}\binom{p+r}{k}s^{p+r-k}u_k\\
 &+b\sum_{j=0}^r\binom{p+r}{p+j}s^{r-j}u_j
 -(b+s)\sum_{j=0}^r\binom rj s^{r-j}u_j.
\end{split}
\tag{5}
\]

Proof: split the defining sum for \(v_{p+r}\) at \(k=p\), substitute
\(u_{p+j}=bu_j+pV_u(j)\), and keep the division by \(p\) exact. Lucas'
binomial congruence gives zero for the lower-block terms with \(k>r\),
and gives \(\binom rj\) for \(k=j\) or \(p+j\). Fermat gives
\(s^p\equiv s\), including \(p\mid s\). These prove the divisibility in
(5). Only after this division may the coefficients of \(V_u(j)\) be reduced
modulo \(p\), giving (4).

For an orientation \(\epsilon=\pm1\), put
\(u'_n=\epsilon^{-n}v_n\), \(b'=(b+s)/\epsilon\). Since \(p\) is odd,

\[
 V_{u'}(r)\equiv\epsilon^{-r-1}
 \left[\sum_{j=0}^r\binom rj s^{r-j}V_u(j)
 +\mathcal C_{p,r}(u,s)\right]\pmod p.
\tag{6}
\]

Arbitrary rescalings are intentionally not covered by this integer carry
formula: additional Fermat-quotient data can enter. Equation (2) itself is
valid for any nonzero scale in characteristic zero.

At zeroth order the first-block polynomials satisfy

\[
 F_v(x)=(1-sx)^{p-1}F_u\!\left(\frac{x}{1-sx}\right)\pmod p.
\tag{7}
\]

Interpret this as a polynomial identity after cancelling each term's
denominator. It follows from
\(\binom{p-1-k}{n-k}(-1)^{n-k}\equiv\binom nk\).
This elegant zeroth-order transport does **not** justify discarding (5).
That correction contains the middle of the first block, binomial carries,
and power carries; it is the precise missing information.

### A cheap hostile and its survivor

For A to C, at \(p=7,r=1\), the actual carry residue is 5; the naive
triangular transform gives 6; the correction is 6, restoring \(6+6=5\).
For A to F the same first failing index gives actual 2 versus naive 6,
with oriented correction 3. The first failed implication is that a transform
of Lucas data also transforms divided differences without a correction.
The strongest survivor is the exact, source-specific affine law (4)--(6).

This is the proposed reusable move: carry the first defect of a compressed
certificate through the change of coordinates, rather than recomputing or
discarding it. Its colouring analogue is only a **research question**: which
residual phase constraints remain after branch contractions? The answer must
come from the literal edge compiler, not from the recurrence analogy.

## 4. Two six-element objects are not the same S4 object

The earlier [K4 duality packet](digit_idoneal_difference_20261018b.md) models
six colours by the six edges of K4. The vertex permutation action of S4 on
these edges is faithful. In contrast, permutation of four marked branch
points acts on a generic six-value cross-ratio orbit through S3; its kernel
is the normal Klein four group.

A decisive witness is the double transposition \((01)(23)\). It fixes all
three partitions into opposite pairs and hence every cross-ratio expression,
but swaps K4 edges \(02\leftrightarrow13\) and
\(03\leftrightarrow12\). Thus no S4-equivariant bijection can identify the
two six-element sets. Genericity matters: special cross-ratios can have
smaller orbits, which cannot repair this six-element identification.

The positive connection is a common space of 24 ordered frames with
**different** four-element stabilizers. Forgetting a frame in one construction
loses different information from forgetting it in the other. Retaining that
frame permits a correspondence, but does not supply an invariant bijection
or a new geometric edge. The opposite-edge colour involution from the older
packet remains useful and is not contradicted by this action obstruction.

## 5. Reproduction and what this changes

Run:

```text
python -B -X utf8 04-computation/experiments/reanchoring_kernel_20261018c.py --out 05-knowledge/results/reanchoring_kernel_20261018c.json
```

The recorded universe is A/C/F indices 0 through 193, every prime 7 through
97 and every carry index \(0\le r<p\), plus all 24 S4 permutations.
There are 10,704 explicit checks, also active under `python -O`. Franel's
independent binomial-cube formula supplies the source oracle; both transformed
sequences are checked against their own recurrences and inverse transforms.
The nonroot \(d=1\), the two literal nonintegral variants, naive carry deletion,
and the double transposition are hostile controls. The program does not prove
Lucas congruences for all primes; that inherited dependency is THM-4711 and
its cited sources. Equations (4)--(6) are proved above for every sequence
satisfying the stated first-block hypothesis.

An independent proof review and optimized replay agreed with the saved JSON;
its LF SHA256 is
`6b352ffdc8d4d6ae272aea4b1af682e7ccbf7f48171bb4c6a38b4636793cfd24`.

The useful next target is a source-retaining transport law, not a numerical
identification of triples. In the geometric anchor, retain a nonempty relative
six-core, all actual edges, all boundary patterns, and all residual component
phases. A successful reanchoring is an exact change of presentation; a
successful obstruction additionally has to eliminate every separating
colouring while retaining an equal-endpoint witness.

## 6. Incoming companion descent: a precise adjacent mechanism

The incoming [carry-companion descent packet](carry_descent_transport_branches_20261018l.md),
commits `e8a1fac0ed` and `f70032d1a6`, extends THM-4711 to the all-prime
Apéry zeta(3) carry divisibility statement. Its descent retains both sheets
of a quadratic map, the pole divisor, and the normalized companion before
identifying the actual carry. This is the relevant connection to the present
work: transformations act on a period **and its first response**, not on a
period alone. Our binomial correction (5) is the response of a different
operation and is proved directly; it is not inferred from that descent.

The incoming note also supplies witnesses for every one of its geometry's
375 local branches. Thus the new CRT compiler is a way to impose and certify
additional joint constraints, not evidence that the old branches have been
eliminated. Its artificial five-rich contact control must retain that scope.

## 7. The next three-state object, and a four-vertex warning

The complete geometric test now retains all 570 boundaries. Its fixed-anchor
triangle already defeats anchor-only pair-label charts; increasing the number
of such chart failures is therefore not the next target. The correct larger
object is a **variable three-state labelling with pinned bipartite fibres**.

Explicitly, write an actual colour uniquely as
\(c_v=3s_v+4\ell_v\pmod6\), where \(s_v\in\mathbb F_2\) and
\(\ell_v\in\mathbb F_3\). For a proposed labelling \(\ell\) agreeing with
every pinned colour modulo 3, retain exactly
those edges whose endpoints have the same \(\ell\). A six-colour lift exists
if and only if this retained graph is bipartite and every pair of pinned
vertices in one component has the prescribed shade difference equal to its
path parity. This is the elementary pair/shade equivalence used by the
inherited compiler: edges with different labels are already satisfied;
edges with equal labels demand opposite bits. It concerns arithmetic CRT
labels, not an automatic additive structure on K4's three directions.

The initial label agreement is essential: one isolated vertex pinned to
colour 0 and proposed label 1 passes every cycle/path test vacuously, but
neither shade lifts it to colour 0. The independent review caught this omitted
prose guard before promotion; the compiler already retained both coordinates
of every pin. This smallest hostile is included in the checker.

Without pins, one can express the condition as follows: **no odd cycle may
be monochromatic in the three-state labelling**. With pins, the path-parity
conditions must also be retained. Thus the next obstruction is a constraint
system on odd cycles and marked paths; it is not a proper three-colouring
problem for the original graph.

The earlier four-vertex warning reappears exactly. In the abstract K4, take
the cycle basis \(012,013,023\), and label vertex 0 by 1 and vertices 1,2,3
by 0. Every basis cycle has mixed labels, but the dependent triangle 123 is
monochromatic. It is the mod-2 sum of those three basis cycles. Monochromaticity
is not a linear predicate on the cycle space. Therefore a fixed cycle basis
cannot replace the full odd-cycle condition. This control is abstract and
does not assert a unit-distance embedding of K4 in the plane.

The script checks the witness and all 81 three-state labellings of K4 against
all 16 binary shade assignments. This supplies a cheap hostile for any future
compression of the variable-label search. A useful compact certificate must
preserve the monochromatic-cycle and pinned-path predicates, rather than
only preserve the linear span of cycles.
