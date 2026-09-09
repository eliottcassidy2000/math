# The two-constant boundary closes geometric reconstruction for quintic mutations

**Status: PROVED ANALYTICALLY + FINITE-EXACT; independently accepted
analytically by the root referee.** This is a separate complete q=2
argument. It closes the precise boundary left outside the three-constant
proof in [Geometric fibre reconstruction](continuing14_20260908_geometric_fibre_reconstruction.md),
without changing that frozen packet.

## 1. Whole quadratic-derivative normalization

For every q=2 member of the continuing13 mutation supplier,

    f(z)=a z(z-beta),  a*beta!=0, Q'=f^2, Q(0)=0.

In coordinate w=z/beta one has

    Q(beta w)=C S(w),
    C=a^2 beta^5/30 !=0,
    S(w)=6w^5-15w^4+10w^3.                             (1)

Thus this is the whole q=2 supplier family, not a sample of primitive
polynomials. The source parameters h and lambda already satisfying the
mutation hypotheses do not occur in its generic affine fibre.

Let L be an algebraic closure of C(T), with the SAME target T throughout.
The generic normalized source curve is P1_L minus the two constant
points 0,1, the five roots of C S(w)=T, and infinity.

The independently proved cross-ratio count lemma applies also to q=2:
each constant finite point belongs to fifteen degenerating quadruples,
and each of the six outer points belongs to five. Consequently any
L-isomorphism of two such curves, extended to their P1 completions,
preserves the unordered pair {0,1}. This deduction does not claim that
two constant points alone determine the Mobius map.

## 2. Exhausting the remaining Mobius maps

Pass to v=w/(w-1). The distinguished constant pair becomes {0,infinity}.
Every Mobius map preserving it has exactly one of the forms

    v -> A v       or       v -> A/v,     A in L^*.

The remaining six punctures are exactly the simple roots of

    B_C(v)=(v-1)[C v^3(v^2-5v+10)-T(v-1)^5]
          =(C-T)v^6-6(C-T)v^5+15(C-T)v^4
            +(20T-10C)v^3-15T v^2+6T v-T.             (2)

The factor v-1 is the old infinity point. It is simple, because
B_C'(1)=6C!=0. The five remaining roots are distinct at the generic
target T and avoid v=1,0,infinity. In particular (2) is a degree-six
squarefree polynomial over L, with nonzero leading and constant terms.

For an isomorphism from C_1 to C_2 in the scaling case, equality of the
six-point divisors gives

    B_(C2)(A v)=k B_(C1)(v),   k in L^*.

Constant coefficients give k=1. Linear coefficients then give A=1,
and leading coefficients give C_2=C_1.

In the inversion case the corresponding equality is

    v^6 B_(C2)(A/v)=k B_(C1)(v).

The leading two coefficients give

    -T=k(C_1-T),   6T A=-6k(C_1-T),

so A=1. The constant coefficient then requires

    (C_1-T)(C_2-T)=T^2,
    C_1 C_2-(C_1+C_2)T=0.                            (3)

Because T is transcendental over C and C_1,C_2 are nonzero constants,
(3) is impossible: its constant coefficient would force C_1 C_2=0.
Thus inversion gives no isomorphism. This proof allows A and k to have
arbitrary algebraic degree over C(T); they were not assumed constant.

We have proved the complete normalized iff: these q=2 geometric generic
curves are isomorphic if and only if C_1=C_2. Undoing normalization in
(1) gives exactly an affine identity Q_2(az+b)=Q_1(z). Conversely any
such identity maps all punctures, so supplies the required isomorphism.

## 3. Combined theorem and scope

Combining this complete q=2 pencil calculation with the independently
proved q>=3 cross-ratio argument gives the geometric reconstruction
theorem for the ENTIRE mutation supplier, q>=2:

    generic original affine source curves are isomorphic
    over an algebraic closure of C(T)
    iff Q_2(az+b)=Q_1(z) for a!=0,b in C.

The previous C(T) theorem is strengthened. The q>=3 supplement correctly
identified the limitation of its own three-point argument; the present
coefficient exhaustion resolves that remaining case by a different
mechanism. No old proved statement is retracted and no previous frozen
artifact is changed.

This theorem concerns the specified generic original source curves and
the same target coordinate. It does not classify arbitrary quintic
first functions on W2, arbitrary rational affine fibrations, or total
fibrations with source parameters included. The full pointed torsion
comparison in the explicit degree-fifteen pair remains unchanged.

## 4. Exact controls

The independent engine checks the all-parameter normalization, the exact
six-point polynomial and discriminant, the full scaling and inversion
coefficient implications, the nonzero-constant obstruction, and lawful
changes of the original beta coordinate. The reflection S(1-w)=1-S(w)
is retained as a hostile to silently changing the target. No producer
engine is imported and no finite enumeration supplies the all-parameter
proof above.

The frozen engine passes 50 always-active exact gates in both modes.

Reproduce with `python continuing14_20260908_quintic_geometric_boundary.py`
and the same command with `python -O`. Output and certificate are frozen
byte-identically with LF newlines; the source writes its certificate
beside itself or into repository 05-knowledge/results after relocation.
