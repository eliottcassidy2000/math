# Actual-source conductor, collision multiplicity, and the minimal-defect cylinder

20 October 2026. **PROVED scoped algebraic statements; independent review pending.**

## Outcome and boundary

The preceding collision packet identified the support of the singular graph. This note restores its multiplicities. In the actual primitive monic graph gauge, one polynomial is simultaneously

- the classical divided-difference D-resultant;
- a generator of the normalization conductor on the source;
- the zeroth Fitting determinant of the finite ordered-collision module;
- the Jacobian-normalized graph derivative `R_Z(P,Q,x)/J`.

Its norm is the square of an explicit order-index determinant. Every source-line specialization has degree twice its total affine delta invariant. Thus an odd specialized degree, or a nonempty collision-free specialization, excludes a proposed noninvertible Keller realization. The smallest possible nonzero generic degree, two, forces the entire graph image ring to be a polynomial cylinder over the ordinary nodal curve.

This is an exact refinement and synthesis of classical resultant/conductor mechanisms, **not a claim of historical novelty**. The D-resultant and its derivative identities are already in van den Essen--Yu and Campbell--van den Essen. No universal contradiction, new global function-field-degree exclusion, or proof of the complex Jacobian conjecture is asserted. In particular, no quartic target is reopened.

## 1. Inheritance, setup, and concept board

The relevant preceding packets are:

- [Finite selector/Euler gate](../inverse_trace_lattice_20261020/FINITE_SELECTOR_EULER_GATE.md);
- [Source-line specialization](../inverse_trace_lattice_20261020/SOURCE_LINE_AND_IRREDUCIBLE_BOUNDARY.md);
- [Collision gap and singular locus](../inverse_trace_lattice_20261020/COLLISION_GAP_AND_SINGULAR_LOCUS.md), and its [independent audit](../inverse_trace_lattice_20261020/collision-review/AUDIT.md).

Closest mechanism: the finite graph normalization and the collision-gap unit. Canonical hostile: the finite immersed nodal-cylinder family that has a collision on every source line but a nonconstant Jacobian. Corrected near miss: singular support, a unit gap, or surjective scalar projection does not force properness. Least-used sidecar: the multiplicity of the finite collision algebra, rather than its reduced points.

Concept board: finite normalization; monogenic trace-dual lattices; divided-difference resultants; affine delta and contact orders; polynomial equalizer rings. The key move is to use the *actual algebra order and its multiplication*, rather than a repaired trace bilinear form alone.

Work over `C`. Suppose

    F=(P,Q): A²_(x,y) -> A²_(u,v),
    J=P_x Q_y-P_y Q_x in C*,
    P monic in y of degree n>=1,
    C(x,P,Q)=C(x,y).

The last condition says that the selected source coordinate x is primitive over `C(P,Q)`. Let

    B=C[x,y],       A=C[x,u] with u=P(x,y),
    S=C[x,P,Q] subset B.

The actual graph `G=(P,Q,x)` is finite birational onto the hypersurface `W=V(R)` because y is integral over S and the fraction fields agree. We use the unique normalization of the irreducible equation R for which it is **monic in v**. This is possible: Q is integral over A, its minimal polynomial over Frac(A) has coefficients in the normal ring A, and its kernel generator differs from R by a nonzero complex scalar. Its v-degree is n. If an earlier graph equation is R'=lambda R for lambda in C*, its derivative quotient is c'=lambda c, whereas the divided-difference resultant below remains c. The norm identity then reads Norm(c')=lambda^n I². The conductor ideal and degree tests are unchanged; the exact scalar identities below always use monic-v R.

No low-trace hypothesis is needed in Sections 2--7. When it is imposed, the previous packet additionally places all conductor support off the nonproperness boundary H=0.

For `a in C`, write

    P_a(y)=P(a,y), Q_a(y)=Q(a,y),
    D_a=F({x=a}), S_a=C[P_a,Q_a] subset B_a=C[y].

The preceding source-line proof supplies that `A¹ -> D_a` is its finite birational normalization with nowhere-zero derivative. Briefly: monicity makes the parametrization finite; normalization of its image is A¹; the induced polynomial cover A¹ -> A¹ has degree one since otherwise its derivative would vanish, contrary to the Keller condition. Thus Q_a is primitive over C(P_a), and its minimal polynomial has degree n.

The specialized graph equation can be identified directly by linear algebra. Let T_Q be the multiplication-by-Q matrix on the free A-module `B=A[y]/(P-u)`. Global primitivity makes its degree-n characteristic polynomial equal to the monic-v equation R. Specializing x=a commutes with taking this determinant, so `R(u,v,a)=det(vI-T_Q(a,u))` is the characteristic polynomial of Q_a on the degree-n extension C(y)/C(P_a). Slice primitivity makes its minimal polynomial degree n as well; hence it equals this characteristic polynomial. Therefore R(u,v,a) is itself the monic minimal polynomial, reduced and irreducible. No analytic reduced-section argument or blanket base-change property for conductors is needed.

## 2. The exact conductor generator

Define, by polynomial substitution,

    c(x,y)=R_Z(P(x,y),Q(x,y),x)/J.                  (1)

It is a polynomial. Differentiating the graph identity gives

    R_u P_x+R_v Q_x+R_Z=0,
    R_u P_y+R_v Q_y=0.

Since J is nonzero constant, solving this pair gives

    R_v(P,Q,x)=c P_y,
    R_u(P,Q,x)=-c Q_y.                             (2)

The polynomial c is not zero: R is separable as a monic v-polynomial and P_y is nonzero.

**Theorem A.** The conductor of S in its normalization B is

    (S:B)={b in B : bB subset S}=cB.               (3)

For every a, separately,

    (S_a:B_a)=c(a,y)B_a.                           (4)

In particular c(a,y) is never identically zero. Statement (4) is proved independently for the specialized rings; it does not assume that arbitrary conductor ideals commute with base change.

### Proof by monogenic trace-duals

Both B and S are finite free A-modules of rank n, with bases

    1,y,...,y^(n-1)     and     1,Q,...,Q^(n-1).

For a lattice L in their common separable fraction field, set

    L^dagger={z : Tr(zL) subset A}.

The elementary monogenic trace-dual formula gives

    B^dagger=(1/P_y)B,
    S^dagger=(1/R_v(P,Q,x))S.                       (5)

For completeness, if f is monic of degree n with root theta, Lagrange interpolation says `Tr(h(theta)/f'(theta))` is the coefficient of `T^(n-1)` in the remainder of h modulo f. The pairing on the power basis has zeros below its anti-diagonal and ones on that anti-diagonal. Its determinant is a unit. This proves (5) over A, not merely over its fraction field, and also proves `L^(dagger dagger)=L` for these free lattices.

For b in B, dualizing reverses lattice inclusion, so

    bB subset S
      iff b S^dagger subset B^dagger
      iff (b P_y/R_v)S subset B
      iff b/c in B.

The last equivalence uses `1 in S subset B`. Since c is already polynomial, this proves (3). Apply exactly the same calculation over `C[u]` to `P_a,Q_a,R(u,v,a)` to get (4). Their field degrees are all n, as justified in Section 1. This also proves c(a,y) is nonzero. QED.

## 3. Divided differences and the finite collision module

Use a second source parameter t and put

    U(x,y,t)=(P(x,t)-P(x,y))/(t-y),
    V(x,y,t)=(Q(x,t)-Q(x,y))/(t-y).

These are polynomial divided differences; U is monic in t of degree n-1. When n=1, U=1 and the collision module below is zero, with determinant convention one.

**Theorem B.** With R normalized monic in v,

    Res_t(U,V)=c(x,y).                              (6)

Consequently, for

    M_col=C[x,y,t]/(U,V),

viewed as a finite B-module,

    Fitt_0^B(M_col)=(c).                            (7)

### Proof

Over the generic splitting field, let t_1,...,t_(n-1) be the other roots of `P(x,t)=P(x,y)`. Because Q is primitive, the monic minimal polynomial of Q has the n corresponding distinct generic Q-values as its roots. Thus

    R_v(P,Q,x)=product_i (Q(x,y)-Q(x,t_i)),
    P_y(x,y)=product_i (y-t_i).

Their quotient is the product of V over the roots of the monic U, namely the resultant in (6). Equality over the fraction field is a polynomial identity. The proof also covers n=1 through the empty-product convention.

The free B-module `E=B[t]/(U)` has rank n-1. Multiplication by V has cokernel M_col and determinant `Res_t(U,V)`. This is a square presentation, so its determinant generates the zeroth Fitting ideal. QED.

The earlier exact gap identity makes t-y a unit in M_col. Therefore this module has no diagonal or ramification contribution for a Keller map: it records *ordered distinct source-line collisions*, including their nonreduced contact multiplicity. The resultant itself is classical. The present connection identifies its precise role in the actual-source graph/conductor packet.

## 4. The norm-square and exact delta identity

Let M be the n by n inclusion matrix expressing `1,Q,...,Q^(n-1)` in the power basis of B over A. It is computed by reducing powers of Q modulo `P(x,y)-u`. Put

    I(x,u)=det M !=0.

**Theorem C.** In A,

    Norm_(B/A)(c)=I(x,u)^2.                         (8)

For every a, let

    delta_a=dim_C(B_a/S_a),

the sum of the affine local delta invariants of D_a. Then

    deg_y c(a,y)
       =2 delta_a
       =dim_C C[y,t]/(U(a,y,t),V(a,y,t)).            (9)

The last dimension is scheme length, not the number of distinct ordered pairs.

### Proof

The two trace Gram matrices satisfy

    Gram_S=M^t Gram_B M.

For a monic degree-n polynomial, the power-basis discriminant is `(-1)^(n(n-1)/2)` times the norm of its derivative. The two signs cancel. Taking determinants and using (2) proves (8) exactly in the monic-v normalization.

Specialize independently to x=a. The matrix M(a,u) still has nonzero determinant. Over the PID C[u], Smith normal form gives

    dim_C(B_a/S_a)=deg_u I(a,u).

For every nonzero polynomial h(y),

    deg_u Norm_(C[y]/C[P_a])(h)=deg_y h.            (10)

Indeed, factor h into linear factors: the norm of each `y-r` is a nonzero constant times `P_a(r)-u`, so has u-degree one. Apply (10) to c(a,y) and specialize (8) to obtain the first equality in (9).

For the second, specialize the multiplication-by-V presentation of Section 3 over C[y]. Its determinant is the nonzero polynomial c(a,y). Its cokernel has finite length equal to the degree of that determinant, again by Smith normal form. QED.

The local geometric refinement is

    ord_(y=s_i)c(a,y)=sum_(j != i) I_p(C_i,C_j),    (11)

where the smooth normalization branches C_i at p are indexed by the source parameters mapping to p, and I_p is plane-curve intersection multiplicity. To see this directly, complete at p and write its reduced local equation as the product of its distinct smooth branch equations. An element supported only on branch i must be divisible by the product of all other branch equations. Its restriction to branch i has valuation exactly the right side of (11). These supported elements characterize the conductor in the product of the branch normalization rings. This proves (11) without a transversality assumption.

## 5. Necessary specialization and coefficient gates

If F is noninvertible, the one-line injectivity theorem of Gwozdziewicz forces a collision on every x=a. Consequently

    deg_y c(a,y) is even and at least 2, for every a. (12)

For an automorphism, S=B and c is a nonzero constant instead. Thus a zero-degree specialization of c already forces F to be an automorphism. An odd-degree specialization cannot occur for any Keller map in this gauge.

Write `c=sum_(j=0)^N c_j(x)y^j`, with c_N nonzero. In the noninvertible case:

- N is even and N>=2.
- For every odd j, `V(c_N,...,c_(j+1)) subset V(c_j)` on the x-line.
- `V(c_N,...,c_1)` is empty.
- The degree at a special parameter may drop only by an even amount and never below two.

These are coefficient tests on an actual proposed source pair or graph. They are necessary, not realization certificates. In particular, an even-degree abstract polynomial need not be a conductor, and a norm-square need not arise from a Keller map.

The norm-square (8) is stronger than pointwise parity: it retains the whole inclusion determinant over C[x,u]. It is the multiplication-aware lattice sidecar missing from a trace-unimodular lattice repair.

## 6. The degree-two sector is a nodal cylinder

**Theorem D (minimal-defect classification).** Suppose F is noninvertible and the generic y-degree of c is two. Then, after a polynomial shear and constant scaling of y that keep x fixed,

    (c)=(y²-1),
    S=C[x]+(y²-1)C[x,y].                            (13)

Equivalently,

    S is isomorphic over C[x] to
    C[x,w,z]/(z²-w²(w+1)),
    w=y²-1, z=y(y²-1).                             (14)

Every D_a has exactly one ordinary node and no other affine singularity. Its two normalization preimages are the two polynomial conductor sections.

### Proof

By (12), every specialization has degree exactly two. Therefore the leading coefficient c_2(x) has no complex zero and is a nonzero constant. By (9), delta_a=1 for every a.

All normalization branches are smooth. Equation (11) shows that each branch over a point with r>=2 branches contributes at least r-1 to the conductor degree. Total degree two therefore allows only one point with two branches, whose pairwise intersection multiplicity is one. Thus this point is an ordinary transverse node; c(a,y) has two distinct simple roots for every a.

The quadratic discriminant `c_1(x)^2-4c_2c_0(x)` is therefore nowhere zero, hence a nonzero constant. Complete the square and rescale to obtain the ideal `(y²-1)`. The allowed change of coordinates is `y -> alpha*y+q(x)` with alpha nonzero constant, so it is a polynomial source automorphism with constant Jacobian and preserves the x-slices.

The two roots y=+1,-1 are the preimages of the unique node on every slice. Consequently

    P(x,1)=P(x,-1), Q(x,1)=Q(x,-1)

as polynomial identities. Hence S lies in the equalizer

    E={f in C[x,y]: f(x,1)=f(x,-1)}.

Dividing f by y²-1 leaves a remainder `a(x)y+b(x)`. The equalizer condition is exactly a(x)=0, so `E=C[x]+(y²-1)B`. Conversely Theorem A gives `(y²-1)B subset S`, and C[x] is already contained in S. Thus S=E, proving (13).

Finally `(y²-1)B` is generated as an algebra together with C[x] by w and z: use `y²=w+1` to reduce even and odd powers. Their displayed irreducible relation proves (14). QED.

This theorem is a classification/reduction, **not yet an exclusion**. It reduces this minimum-conductor sector to the precise question whether the nodal equalizer admits P,Q with

    {P,Q}=J in C*,       C[x,P,Q]=C[x]+(y²-1)B,

together with the original source/primitive assumptions. Merely taking P,Q in the equalizer is weaker than requiring the generation identity. An area integral between y=-1 and y=1 gives only a derivative identity and no contradiction. No claimed bound `delta_a >= generic_degree(F)-1` is used; no such bound was established in the checked sources.

## 7. Controls and what they distinguish

`verify_conductor.py` runs **293 exact checks** over 30 parameterized source pairs, plus coefficient and equalizer controls. It uses explicit exceptions rather than Python assert, so optimized Python still executes every check.

Positive Keller controls are polynomial automorphisms

    P=x+y^r, Q=k y+P^m,
    r in {1,2,3,4}, k in {1,2,-1}, m in {2,3}.

With monic-v graph equation `R=(v-u^m)^r+k^r(Z-u)`, they have J=k, conductor generator `c=k^(r-1)`, and index determinant `k^(r(r-1)/2)`. This checks the scalar normalization, norm-square, and direct divided-difference resultant rather than silently setting every Jacobian to one. The n=1 empty-collision-module convention is checked as well.

The collision-bearing controls are **not Keller**:

    P=y²,
    Q=x+y(y²-1)^k,
    R=(v-Z)²-u(u-1)^(2k),          1<=k<=6.

Every slice is a finite immersed rational curve with the same pair y=+1,-1 colliding, with contact order k. Its conductor is `(y²-1)^k`, its total delta is k, and its index determinant is `(u-1)^k`. However J=-2y and

    R_Z(P,Q,x)=-2y(y²-1)^k

has odd y-degree 2k+1. Thus the unnormalized graph derivative cannot be used as the conductor outside the Keller hypothesis. The k=1 case is an ordinary nodal cylinder; k>=2 shows why two reduced collision points do not mean conductor degree two.

The cusp control `P=y²,Q=x+y³` has conductor y² and passes degree parity and norm-square, but violates immersion. The node-to-cusp family `P=y²,Q=y(y²-x)` has conductor `y²-x`, of degree two everywhere, but its discriminant 4x vanishes. At x=0 the slice derivative vanishes. These show why the immersion hypothesis is indispensable for Theorem D.

Finally, `x y²+y+1` has a forbidden odd degree drop at x=0, whereas `x² y⁴+y²-1` passes the numerical parity/drop test without being claimed realizable.

Reproduce:

    python3 verify_conductor.py
    python3 -O verify_conductor.py

## 8. Attribution, exact scope, and next target

Primary sources checked 20 October 2026:

1. A. van den Essen and J.-T. Yu, *The D-resultant, singularities and the degree of unfaithfulness*, Proc. AMS 125 (1997), 689--695. [Author-institution PDF](https://hub.hku.hk/bitstream/10722/156137/1/content.pdf). Section 1 defines the divided-difference resultant; Proposition 1.1 gives the derivative identities. The monic-v normalization used here absorbs the source's resultant sign. This paper also credits Abhyankar's characteristic-zero Taylor resultant.
2. L. A. Campbell and A. van den Essen, *Jacobian pairs, D-resultants, and automorphisms of the plane*, J. Pure Appl. Algebra 104 (1995), 9--18. [Publisher DOI](https://doi.org/10.1016/0022-4049(95)00116-E). The publisher abstract was checked; its principal constant-Jacobian/constant-D-resultant automorphism criterion is classical. The full-paper mirror was blocked, so no uninspected theorem is used.
3. J. Gwozdziewicz, *Injectivity on one line*, 1993, Theorem 1.1. [Primary arXiv PDF](https://arxiv.org/pdf/alg-geom/9305008). This is the sole non-elementary external input for the strict lower bound in (12); the conductor, norm-square, and delta identities do not depend on it.

The work advances the earlier requested singular-locus target by giving exact contact multiplicity and a polynomial elimination obstruction. It does **not** show the required conductor has an odd or zero-degree specialization in every hypothetical Keller map. That remaining existence step is genuinely open here.

Two precise continuations are now available:

- For a named low-trace/source family, compute the D-resultant or `R_Z(P,Q,x)`, then seek an odd/zero specialization, a failed norm-square, or an incompatible branch-contact divisor.
- In the minimum-defect sector, analyze constant Poisson brackets and the full generator condition in the explicit nodal-cylinder equalizer. A proof that no such pair exists would close that sector. Neither arbitrary pairs nor all higher conductor degrees have been excluded.

The transfer contract is: actual source order S subset B -> collision module -> conductor divisor -> slice delta. The preserved information is normalization defect including contact multiplicity. Reducing to the singular support destroys the weights; reducing to even degree destroys locations and contact pairing; reducing to the norm square destroys the full source algebra. The actual inclusion matrix and equalizer are the required sidecars. There is no transfer of a conjecture solution from a shared trace or determinant alone.
