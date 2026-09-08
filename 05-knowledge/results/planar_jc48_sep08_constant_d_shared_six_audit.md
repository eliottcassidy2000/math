# Independent audit: the constant-D shared finite-six class

**Status: PROVED / INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.**
Root accepts [constant_d_shared_six](planar_jc48_sep08_constant_d_shared_six.md), with its full every-finite-position coefficient class: no polynomial Jacobian mate exists, and after the necessary local entry its five-case table completely classifies rational-mate existence. The fixed DG surface, square-prefix representation, constant-D hypothesis and finite-six/infinity-two boundary divisor are essential. JC(2) remains open.

Root supplied the initial coefficient target and one explicit rational submersion. The geometry producer independently found and wrote the complete case classification, including the conic residue argument, every degenerate coefficient case and the original-source pole obstruction. Root then independently read the whole resulting proof and standalone implementation, rederived the residue and regularity arguments, and executed normal and optimized replays. Both reproduce the frozen 68-gate output. A further direct original-source calculation checks all four rationally accepted cases and the source divisors used for polynomial exclusion.

## Full globality and the necessary active entry

The numerator-box calculation is complete. With N=(x-p)^6, deg A,B,C<=4 and N=Ax^4+Bx²+C, one must have deg A<=2. The top two coefficients of B are forced, while its three remaining coefficients and the three coefficients of A are free. The displayed shifted rows follow, retaining p and all six parameters; their converse lies in the full global numerator box. The separate L1 box is treated in the same way.

On the actual infinity divisor D, H and L restrict to affine functions with slopes 2-A4 and -B4. Constancy of H²+L forces both slopes to vanish, because an affine L cannot cancel a nonzero quadratic term of H². This gives A4=2, B4=0 exactly. The actual constants H|D=4p²+2pc+k and L|D=2pn+e retain the dependence on the original finite location. Their disappearance from the ordinary source formulas does not assert that translation is a surface automorphism.

For the full constant-D rows, the complete infinity numerator has leading boundary term r² times a unit, and its normal constant and the lower numerator's constant both vanish. The rescaling s=rZ has the quartic displayed in the primary. It has nonzero constant one; a repeated nonzero root must satisfy the fixed eliminant ZP0'-4P0, whose constant is -4. Thus, outside finitely many fibre values, all four roots are nonzero and simple and exhaust the local Weierstrass degree. The weighted relative form has order 2+2-3=1. There are no primitive poles at infinity. This check includes the canonical r² and does not import the finite logarithmic statement there.

At the shared finite point, a normal unit makes all branches regular. Since infinity has no primitive poles, a proposed rational primitive would then be constant on every compact generic component, a contradiction. Otherwise the finite first-jet theorem produces nonzero logarithmic residues unless beta1=m1=0. Thus beta0=beta1=m1=0 is a necessary rational entry. The inherited first-jet statement is a rational obstruction, not a polynomial-only source critical-point test. Outside this entry there is no rational mate; inside it the full remaining coefficient family is exactly the one classified in the table.

Generic components meet the original source, where the relative form is nonzero. No generic connectedness is assumed for this entry argument. In particular the later composition cases are not incorrectly removed at this step.

## Actual fields and the five exhaustive coefficient cases

Put q=1+u²t, v=uq, h=H=v²+cv+bq+k, L=dq+nv+e. The inverse u=v/q, t=(q-1)q²/v² proves C(u,t)=C(v,q). The Jacobian J(v,q)=v²/q gives the exact source volume q dv wedge dq/v². When b is nonzero, h,v are rational field coordinates and the volume is q dv wedge dh/(b v²). These identities are used for rational exactness only; original polynomiality is checked separately.

For b*d nonzero, eliminating q gives a conic F=h²+(d/b)h-(d/b)v²+(n-dc/b)v+e-dk/b. At its two generic infinities, write w=1/v and xi=h/v. The local equation has xi²=d/b at w=0, with two nonzero simple roots. The complete transformed differential has residues -epsilon/(2b²sqrt(d/b)). They are nonzero for every free remaining parameter. Thus this whole case has no rational mate. There is no unproved outer-residue condition or sampled coefficient restriction in this implication.

For b=0,d nonzero, C(F,v) is the full source field. Its fibre differential is -(f-P(v)-e)dv/(d²v²), where P=(v²+cv+k)²+nv. The only finite residue is (2kc+n)/d². Its vanishing is necessary and sufficient for rational integration, and the displayed primitive is correct. On the original line u=0, F_t=0 and F_u=2kc+n. Consequently every rationally accepted coefficient choice has an actual affine critical line and cannot have a polynomial mate. This retains the source line which the rational inverse chart deletes.

For b=d=0, F=P(v)+e with outer degree four. Its derivative factor is nonconstant, excluding polynomial mates. The identity J(v,1/(2u²))=1 supplies a rational mate after division by P'(v), for every parameter. The proof does not incorrectly claim that the constant field is C(F) in this composite case, and it never divides by n here.

For b nonzero,d=0,n nonzero, C(F,h) is the full source field. Substituting v=(f-h²-e)/n in the actual form gives two generic residues epsilon(nk+2c rho²)/(4b²rho³), rho²=f-e. Their vanishing for the independent fibre parameter forces c=k=0, and conversely the explicit G0=1/(2b²v)-H/(n b²) has source Jacobian one. In this noncomposite field C(F)(h), the constants of the nonzero Hamiltonian derivation are exactly C(F), so every rational mate is G0+R(F). This kernel is paid by the explicit field coordinate, not inferred merely from generic connectedness.

For b nonzero,d=n=0, F=H²+e. Multiplying a mate by 2H or dividing by it shows rational-mate equivalence with H. At fixed H, the relative form has sole residue c/b². Thus c=0 is the exact rational criterion, with the displayed primitive. The polynomial factor 2H excludes all polynomial mates, irrespective of whether that residue vanishes.

These five disjoint cases exhaust the parameters, including every zero-coefficient degeneration. All introduced denominators are used only in their declared nonzero case. The constants k,e and the arbitrary original position p are retained.

## The noncritical rational family and the full source pole obstruction

The case b*n nonzero,d=c=k=0 is an actual polynomial submersion, not a critical-point near miss. Its displayed rational mate is regular wherever v is nonzero, which proves the source gradient nonzero there. The complement is v=uq=0: directly F_u=n on u=0, and F_t=nu³ on q=0, where u is nonzero. Both are nonzero. Thus every affine point is covered.

Nevertheless every rational mate has an unavoidable source pole. The primary's exact factorization

    F-(b²+e)=u K,
    K=u(q²+bt)(H+b)+nq

is correct. The cofactor satisfies K(0,t)=n, K(u,0)=u³+2bu+n and K(u,-1/u²)=-b²/u. It is nonconstant and has a nonempty curve divisor over C; that divisor avoids both u=0 and q=0. Hence G0 is regular on every one of its components.

On u=0, G0 has a simple pole and F-(b²+e) has order one. Any invariant correction R(F) canceling it must have a pole at b²+e. The same correction then has a pole along every component of K=0, where G0 is regular. Cancellation is impossible. A higher-order pole of R would already produce an uncancelled higher-order pole on u=0, so it cannot evade this argument. This proves absence of a polynomial mate without making any claim that rational exactness or source submersion alone is impossible.

I independently read the canon antecedent [THM-3978, linear-seam-submersion-rational-mate-pole-obstruction](../../01-canon/theorems/THM-3978-linear-seam-submersion-rational-mate-pole-obstruction.md) and the proved [octuple transport, Section 5](planar_jc48_sep08_octuple_transport.md). They supply the recovered mechanism; the new proof correctly rederives the actual rational field, special fibre and pole divisors before using it. No external novelty claim is part of this audit.

## Relation to the earlier M-unit class

The older fixed-point M-unit theorem remains true within its own hypothesis. Its all-finite-position extension is [const_d_translation](planar_jc48_sep08_const_d_translation.md), which pays the changed global section intersection rather than an unsupported surface translation. The new shared theorem is its complementary M(p)=0 case under the same constant-D and boundary-divisor hypotheses. Once this shared theorem is promoted, the two proved results therefore exclude polynomial mates for the complete constant-D finite-six/infinity-two class at every finite position. Rational mates persist in the shared branch, including the nonconstant-L submersion; the union must not be described as a rational exclusion.

This combined consequence still requires constant D. It says nothing about the remaining nonconstant-D finite-six/infinity-two class, other root partitions, arbitrary quartics or JC(2).

## Independent literal source controls

The following separate calculation constructs the original polynomials and differentiates every rationally accepted split, along with the critical-line and same-fibre controls. It imports no producer code:

```python
import sympy as S
u,t,b,c,k,d,n,e=S.symbols('u t b c k d n e')
q=1+u*u*t;v=u*q
H=v*v+c*v+b*q+k
F=H*H+d*q+n*v+e

def jac(A,B):return S.cancel(S.diff(A,u)*S.diff(B,t)-S.diff(A,t)*S.diff(B,u))
def check(name,v):
 if S.cancel(v)!=0:raise RuntimeError(name)
 print(name,'PASS')
F1=F.subs({b:0,n:-2*k*c})
G1=(v**3/S.Integer(3)+c*v*v+(c*c+2*k)*v+(F1-e-k*k)/v)/(d*d)
check('b0 accepted rational family',jac(F1,G1)-1)
check('b0 original critical line u derivative',S.diff(F1,u).subs(u,0))
check('b0 original critical line t derivative',S.diff(F1,t).subs(u,0))
F2=F.subs({c:0,k:0,d:0});H2=H.subs({c:0,k:0})
G2=1/(2*b*b*v)-H2/(n*b*b)
check('nonconstant L rational submersion',jac(F2,G2)-1)
K=S.cancel((F2-b*b-e)/u)
check('same fibre cofactor at u0',K.subs(u,0)-n)
check('same fibre cofactor at q0',K.subs(t,-1/u**2)+b*b/u)
check('submersion at u0',S.diff(F2,u).subs(u,0)-n)
check('submersion at q0',S.diff(F2,t).subs(t,-1/u**2)-n*u**3)
H3=H.subs(c,0)
G3=(v+(H3-k)/v)/(2*b*b*H3)
check('constant L accepted rational family',jac(H3*H3+e,G3)-1)
w=S.symbols('w');P=(w*w+c*w+k)**2+n*w+e
G4=1/(2*u*u*S.diff(P,w).subs(w,v))
check('composite accepted rational family',jac(P.subs(w,v),G4)-1)
print('Independent original-source shared-six identities PASS')
```

All ten independently constructed identities pass. These supplement the full analytic case proof and the producer's all-parameter replays; they are not a finite parameter census.

## Frozen source and final acceptance

I read every producer gate. The global numerator box, full shifted chart, complete infinity tangent, actual field inverses and Poisson coefficients, three residue mechanisms, all rational positives, polynomial derivative factors and the original-source submersion/pole divisors are checked symbolically. Branch exhaustion, field kernels and the all-parameter case distinction are the written analytic proof; simple order-bookkeeping gates do not replace them.

Both commands run independently from the worktree root:

    python3 04-computation/planar_jc48_sep08_constant_d_shared_six.py
    python3 -O 04-computation/planar_jc48_sep08_constant_d_shared_six.py

Each passes 68 always-active gates and reproduces all 392 frozen output bytes.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source | 9140 | `798b7824988fa6cbdcaaccd2689997d29ad57473e64206ca8b7f42a1e31d4aa6` |
| Output and both replays | 392 | `ecd0a552a4bde15e90bebc8dd5764f6105dc24d05463107e7d30404b048c410c` |
| Primary before promotion | 15054 | `352e67471e66848acb37c833c6d9b05c16ff4eb29f613880fc403290b0c57064` |

No mathematical or source correction was required. The source/output remain frozen. This audit accepts promotion of the complete shared theorem and rational table, and the precisely scoped constant-D polynomial corollary using the separately proved M-unit supplier.
