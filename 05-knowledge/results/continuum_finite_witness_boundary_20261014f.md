# Continuum changes and finite exact witnesses

**Status: PROVED elementary finite-witness and local-sieve statements; CITED classical transfer/independence inputs; FINITE-EXACT checker controls. No solution of CH, the six/seven plane-colouring question, the Jacobian conjecture, or Guy D19 is supplied.**

The useful common structure is an exact finite witness, with its arithmetic domain and source data retained. A forcing extension can add real numbers and global functions while leaving every such finite witness unchanged. This gives a precise restriction on a continuum-based approach; it does not give a practical search bound or a general theorem that stronger axioms cannot help prove arithmetic.

## 1. Inheritance and the maps being compared

The closest existing construction is [the finite algebraic plane-witness packet](plane_algebraic_witness_20261011.md): compactness followed by real-closed-field transfer, with exact coordinates, distinctness and all edge/nonedge equations retained. Its positive control is the seven-vertex four-chromatic Moser spindle, not a six-chromatic graph. The [finite JC derivative-closure packet](jc_bilinear_reconstruction_20261013a.md) instead consumes polynomial identities and exports an inverse. The [marked Bott packet](collatz_bott_marked_periodicity_20261007d.md) separates an injective marked matrix carrier from a stable or homotopy quotient. The [parity continued-fraction packet](collatz_parity_cf_boundary_20261014a.md) shows why a common irrationality exponent can forget an otherwise lossless source code.

The current canon already states the HN arithmetic/choice boundary in [THM-4622, logical-status dictionary](../../01-canon/theorems/THM-4622-logical-status-dictionary-absoluteness-choice-and-automorphisms-for-jc2-collatz-hadwiger-nelson.md). Its MISTAKE-594 correction restricts the unparameterized field claim to computably presented algebraic fields and separates it from a fixed countable field shared by two models. We retain that correction: no theorem about arbitrary unspecified noncomputable field presentations is used here. [HYP-2232, cardinal-shadow carrier](../hypotheses/HYP-2232-continuum-hypothesis-cardinal-shadow-carrier.md) is earlier method provenance, not an additional set-theoretic theorem. Our contribution is the explicit finite-verifier interface, the bounded JC dual search, and the source-typed D19 sieve/hostiles, not a new HN absoluteness theorem.

The live concepts here are finite codes, field restrictions, global models, local congruence filters, and source markers. Their interfaces are:

| Source | Target/map | Preserved predicate | Data that must remain |
|---|---|---|---|
| Finite real unit graph | Algebraic realization of the same finite formula | All specified distances, distinctness, and colour obstruction | Labelled vertices, metric equations, palette and refutation |
| A complex Keller counterexample of fixed finite degree | An algebraic solution of its coefficient/collision equations | Normalized determinant and a marked collision | Dimension, degree bound, characteristic zero, differing-coordinate witness |
| Rational distances to four square corners | Primitive integer scaling | Four exact square equations | The *same* side length in both directions; all four distances |
| An exact finite witness in a ground model | The same finite code in a forcing extension | Acceptance by the fixed verifier | Fixed finite program and finite parameters |

The decisive hostile examples are an algebraic distance that is not rational, a characteristic-$p$ Keller collision that is not characteristic-zero Keller, and a rectangle witness that is not a square witness.

## 2. What the real-number classifications do and do not say

There are inclusions
\[
\mathbb Q\subsetneq\mathbb R_{\rm alg}\subsetneq\mathbb R.
\]
The irrational reals include both algebraic irrationals and transcendental reals. Thus “real versus transcendental” is not a complementary classification. A Dedekind cut represents a real number; it does not by itself give a finite certificate of its arithmetic type.

Every nonempty rational interval contains a rational, an algebraic irrational and a transcendental number. For the middle case, take a rational point safely inside the interval and add a sufficiently small nonzero rational multiple of $\sqrt2$. For the last case, a real interval is uncountable and the algebraic numbers are countable. Consequently no observation consisting only of a nonempty rational interval can certify which of these three classes contains the unknown point.

This does not make exact classification impossible for a suitably specified source: a defining polynomial with an isolating interval is additional information. Nor does an irrationality exponent supply that information. The inherited continued-fraction constructions give different exactly marked sources with the same exponent $2$; the exponent is an invariant with large fibres, not a reconstruction map.

Likewise, an equidistribution statement needs a sequence and a measure; an equidecomposition needs a permitted group of motions, pieces and regularity conditions. Equal cardinality supplies neither. For example intervals of lengths one and two are bijective but cannot be finitely partitioned into measurable pieces that match by translations: finite additivity and translation invariance preserve length. This elementary obstruction has no dependence on CH.

## 3. Finite-code absoluteness

**Theorem A.** Let $M\subseteq N$ be transitive models of ZFC with the same natural numbers. Let $V$ be a fixed terminating exact verifier specified by a finite program and finite integer parameters. Then the set of finite strings accepted by $V$ is the same in $M$ and $N$. In particular
\[
 M\models\exists w\ V(w)=1
 \quad\Longleftrightarrow\quad
 N\models\exists w\ V(w)=1.                         \tag{1}
\]

**Proof.** Both models have exactly the same finite strings over a finite alphabet. On a fixed input, every step of a finite integer computation is the same, so its final acceptance decision is the same. The quantifier in (1) therefore ranges over the same codes and the same predicate. Equivalently one can include a finite accepting computation in the witness; then the verification of that computation is bounded. No absoluteness theorem about arbitrary real quantifiers is needed. $\square$

The same argument extends to every first-order arithmetic statement with fixed finite integer parameters: induction on its formula preserves atomic arithmetic and Boolean operations, and each number quantifier ranges over the same natural numbers. Thus the usual universal Collatz statement, written $\forall n>0\;\exists k\;T^k(n)=1$ with its fixed integer iteration rule, is also unchanged between these models. This proves no instance or universal convergence theorem; it says that changing CH through forcing is not a way to change its truth value.

A forcing extension is a principal instance: it adds no natural numbers or finite strings. **CITED:** Gödel's constructible-universe result and Cohen's forcing result give the relative independence of CH from ZFC, assuming ZFC is consistent; see [Gödel, PNAS 1938](https://doi.org/10.1073/pnas.24.12.556) and [Cohen, PNAS 1963, pp.1143–1144](https://www.pnas.org/doi/pdf/10.1073/pnas.50.6.1143). These are not unconditional proofs that a transitive model of ZFC exists. Theorem A is the elementary finite-code part of this discussion, proved above.

In particular, forcing CH to hold or fail cannot change the truth of a problem already reduced to (1). This is a semantic statement about these models and this fixed problem. It is not a claim about arbitrary nonstandard models, arbitrary stronger axioms, or new real parameters that were absent from the original problem. It also gives no bound on the first accepted code: a complete enumeration is a semidecision procedure, and a finite unsuccessful search proves no global nonexistence.

### 3.1 Plane-colouring lower bounds

For a fixed integer $k$, the statement that the Euclidean plane has no proper $k$-colouring is equivalent, in ZFC, to existence of a finite non-$k$-colourable unit-distance graph. This is the compactness step inherited from [the prior plane packet](plane_algebraic_witness_20261011.md). Its finite geometric formula includes
\[
 (x_i-x_j)^2+(y_i-y_j)^2=1
\]
for edges, inequality from $1$ for specified nonedges, and positive squared distance for distinct vertices. Real-closed-field transfer gives a realization in the real algebraic numbers. A finite algebraic-number encoding plus a finite exhaustive colour refutation is an exact finite witness.

Thus each finite-palette lower bound, and hence the usual finite plane chromatic number, is unchanged by forcing. Even though the larger model has new real points, any newly alleged finite obstruction would have an algebraic finite code already available in the smaller model. Conversely the old obstruction remains valid. Global colourings may be different functions; compactness preserves existence at the same palette size. This assertion concerns unrestricted proper colourings under ZFC, not measurable, definable, or choice-free variants.

Conditional on a non-five-colourability premise, fair algebraic-witness enumeration eventually finds a finite obstruction and vertex-critical deletion gives a six-chromatic subgraph. It supplies neither its size nor a deadline. A finite six-chromatic witness still does not decide whether the full-plane answer is six or seven. The existing finite-field stabilization theorem remains useful: some finite real-number field attains the full finite chromatic number, but no field or degree bound is selected by this argument.

### 3.2 A finite algebraic JC failure predicate

Fix dimension $n$ and a finite degree bound $d$, in characteristic zero. A constant nonzero Jacobian determinant can be normalized to $1$ by scaling one target coordinate. **CITED:** [Ax, 1969, Corollary 1, p.1](https://msp.org/pjm/1969/31-1/pjm-v31-n1-p01-p.pdf) gives surjectivity of injective polynomial self-maps over algebraically closed fields; the characteristic-zero automorphism conclusion is also stated and proved in [Kaliman, Lemma 1, pp.1–2](https://www.math.miami.edu/~kaliman/library/Ka.PAMS.2005.ax.pdf). Therefore a nonautomorphic complex Keller self-map has a collision.

Introduce variables for the coefficients of $F=(F_1,\ldots,F_n)$, two points $a,b$, a coordinate $j$, and a scalar $u$. The finite equations are
\[
 \det JF\equiv1,\qquad F_i(a)=F_i(b)\ (1\le i\le n),
 \qquad u(a_j-b_j)=1.                              \tag{2}
\]
The determinant equation means every coefficient identity, not evaluation at one point. The last equation records distinctness without an informal inequality.

If (2) has a complex solution, its ideal over $\mathbb Q$ is proper, hence it has a solution over $\overline{\mathbb Q}$ by the weak Nullstellensatz. All finitely many entries belong to a number field and admit a finite exact encoding. Conversely any verified algebraic solution is a characteristic-zero Keller collision. Taking the enumerable union over $d$ and $j$ proves that failure of JC in a fixed dimension has a finite-witness formulation of the form (1). The argument is a logical interface and makes no claim to settle any dimension. It may replace an initially transcendental-coefficient counterexample by a different algebraic-coefficient one; it does not reconstruct or preserve arbitrary extra properties of the original map.

There is a positive bounded decision target as well. For each fixed $(n,d,j)$, write the finite equations (2) as $E_1=\cdots=E_s=0$. Dovetail two exact searches: algebraic solutions, and rational polynomial identities
\[
 1=\sum_{i=1}^s H_iE_i.                              \tag{2a}
\]
If there is no algebraic solution, the weak Nullstellensatz supplies the identity; its coefficients can be taken rational because any finite coefficient comparison is a consistent linear system over $\mathbb Q$. If there is a solution, an exact algebraic encoding eventually appears. The two outcomes exclude each other and one search halts. Thus the fixed degree collision problem is decidable, with a finite positive or negative certificate. Running through unbounded degrees does not produce a universal stopping rule or a proof of JC. This is a concrete surviving target for computation, not a new degree bound.

The characteristic marker is essential. In characteristic $p$, $f(t)=t-t^p$ has derivative $1$ and $f(0)=f(1)=0$. The identical integer polynomial in characteristic zero has derivative $1-pt^{p-1}$, so it is not Keller. The checker tests this distinction for $p=2,3,5,7,11$ and refuses positive-characteristic packets at its characteristic-zero interface. Its positive controls are supplied rational polynomial inverse pairs, checked in both composition directions. It supplies no characteristic-zero collision.

## 4. Guy D19: a rational witness cannot be replaced by an algebraic one

**CITED source identification.** Berry's [*Points at Rational Distance from the Corners of a Unit Square*, 1990, pp.505–506](https://www.numdam.org/item/ASNSP_1990_4_17_4_505_0.pdf), explicitly identifies Guy §D19 and distinguishes three-distance families from the four-distance question. His question permits a point anywhere in the plane of the square. An inside-square requirement is a separate additional condition. The relevant two pages and Fig.1 were inspected visually. This note does not infer a current solution from that historical paper.

For a unit square let the four rational distances be $a,b,c,d$, with $a,b$ corresponding to $(0,0),(1,0)$ and $c$ to $(0,1)$. Subtracting squared distances gives
\[
 x=(a^2-b^2+1)/2,\qquad y=(a^2-c^2+1)/2,
\]
so the point coordinates are automatically rational. Clearing denominators and dividing their common content gives the equivalent finite integer system
\[
\begin{aligned}
 A^2&=X^2+Y^2,& B^2&=(X-L)^2+Y^2,\\
 C^2&=X^2+(Y-L)^2,& D^2&=(X-L)^2+(Y-L)^2,             \tag{3}
\end{aligned}
\]
with $L>0$, $X,Y\in\mathbb Z$, $A,B,C,D\ge0$, and $\gcd(X,Y,L)=1$. Dividing by the coordinate content also divides every distance: $g^2\mid A^2$ implies $g\mid A$, and similarly for the others. Conversely (3), divided by $L$, gives the desired rational point and distances. For the strict interior variant add $0<X,Y<L$; no positivity restriction on $X,Y$ belongs to the unrestricted question.

This is another finite-witness existence problem governed by Theorem A. Rationality remains indispensable: the centre of the unit square has four distances $\sqrt2/2$. It is a real-algebraic solution of the distance equations, but not a rational-distance solution. Real-closed-field transfer cannot remove this distinction.

### 4.1 A complete cheap local guard

**Theorem B.** Every primitive integer solution of (3) has
\[
 12\mid L,\quad
 \text{one of }X,Y\text{ is divisible by }4\text{ and the other is odd},
 \quad \text{exactly one of }X,Y\text{ is divisible by }3.           \tag{4}
\]
In particular all four distances are coprime to $6$.

**Proof.** If $L$ is odd, each of the pairs $X,X-L$ and $Y,Y-L$ contains an odd number. One corner would then have squared distance $2\bmod4$. Thus $L$ is even. Both coordinates cannot be odd (the same obstruction), nor both even (primitivity), so exactly one is even. In a sum of an odd square and an even square that is a square, the even coordinate must be $0\bmod4$, since $1+4=5\bmod8$ is not a square. Applying this to that coordinate and its translate by $L$ gives $4\mid L$.

If $3\nmid L$, each coordinate pair contains a nonzero residue modulo $3$. A corner then has squared distance $2\bmod3$, impossible. Hence $3\mid L$. The coordinates cannot both be nonzero modulo $3$, nor both zero by primitivity. This proves (4). Every corner then has odd squared distance $1\bmod3$, proving the last assertion. $\square$

There is an exact finite converse *only at this local level*. For triples $(X,Y,L)\bmod24$ satisfying $\gcd(X,Y,L,24)=1$, all four corner norms are square residues modulo $24$ if and only if (4) holds modulo $24$. In fact all four norms are then $1\bmod24$. Necessity is the preceding proof using only the $2$- and $3$-parts of primitivity. Sufficiency follows directly modulo $8$ and $3$. There are exactly $128$ such triples: two choices for $L\bmod8$, sixteen coordinate pairs with the stated parity condition, and four coordinate pairs modulo $3$, with $L=0\bmod3$. Equivalently $L=0$ or $12\bmod24$, with $64$ coordinate pairs for each.

This does not make the local guard sufficient for (3). The primitive interior point
\[
 (X,Y,L)=(3,4,12)
\]
has corner norms $(25,97,73,145)$, all $1\bmod24$, while only the first is an integer square. The checker retains all four exact equations. A separate positive rectangle control, sides $6,8$ with point $(3,4)$ and distances $(5,5,5,5)$, is rejected by the square interface because its two side lengths differ.

## 5. What remains a productive target

The invariant statement is about existence of a correctly encoded finite witness. It does not turn a local admissibility table into a witness, a dense class into an effective classification, or a stable homotopy class into a marked coordinate object. The surviving tasks are concrete:

* HN: build an exact algebraic unit graph with a complete five-colour refutation, retaining all cross-piece metric constraints.
* JC: prove exact polynomial closure/inverse identities, or provide a characteristic-zero coefficient-and-collision packet satisfying (2).
* D19: solve all four equations (3) simultaneously on a primitive source surviving (4); an ever larger local sieve remains only a filter unless a valid global lifting theorem is supplied.

None of these obligations is discharged by choosing CH or its negation. The present theorem says forcing cannot create their finite witnesses. A proof using additional principles must still establish the exact source predicate; this note does not make a blanket claim about the proof-theoretic strength of arbitrary new axioms.

## 6. Reproduction and scope of the checker

Files: [exact checker](../../04-computation/experiments/continuum_finite_witness_boundary_20261014f.py) and [saved JSON](continuum_finite_witness_boundary_20261014f.json).

```powershell
python -B -X utf8 04-computation/experiments/continuum_finite_witness_boundary_20261014f.py
python -B -O -X utf8 04-computation/experiments/continuum_finite_witness_boundary_20261014f.py
```

The finite universe is all $24^3=13,824$ local residue triples, 100 explicitly specified two-shear polynomial automorphisms, five characteristic-change hostiles, and the inherited exact algebraic Moser certificate. The production rational polynomial checker rejects a pointwise Jacobian test, an unmarked coincident pair, wrong characteristic and malformed exact-number inputs. The D19 checker allows exterior points by default and retains an explicit optional interior guard. No actual D19 square receipt is supplied. The code is not a general algebraic-number JC enumerator or a computation of set-theoretic absoluteness; those claims are the proved reductions above.

Normal and optimized runs agree: **14,096 checks**. Saved JSON LF SHA256: `71f9ab41448332bc9d0a813b34fa43fa4a970a1fba18f1bee299896df9db0a72`.
