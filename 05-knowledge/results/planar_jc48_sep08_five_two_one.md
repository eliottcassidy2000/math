# The binary (5,2,1) class from exact inverse coefficients

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a complete coefficient-family argument in the fixed DG surface.
It is not a finite mate-degree calculation or a general quartic closure.

## 1. Statement, inheritance and full global coefficient space

Use `W=(P1_x x P1_z)\{z=x^2}`, `t=1/(z-x^2)`, and the original
volume `dx wedge dt`. Let `H in L2`, `L in L1`, `deg_t H=2`, and
`F=H^2+L`. Assume its leading binary octic has three distinct zeros of
multiplicities five, two and one.

**Theorem.** No polynomial mate of any degree exists. When the
double point is infinity, a rational mate forces `L` constant. When the
double point is finite, no rational mate exists, irrespective of `L`.
Section 7 gives exact rational mates with constant `L` in the permitted
infinity placement, for every position of the finite fivefold point.

The finite-double alternative follows directly from the
[proved leading differential gate](planar_jc48_sep08_leading_exactness.md):
`du/sqrt(N)` has a nonzero simple pole at a finite root of multiplicity two.
Thus it suffices to put the double point at original infinity. Actual
source/surface scaling normalizes separation and leading scalar; write
`u=x-p` and

    N=u^5(u-1).

The translation is only a polynomial source coordinate, not an assumed
automorphism of `W`. The complete global rules in the original `x`, from
[the full DG filtration](planar_jc48_sep08_dg_quadratic.md), give

    P=a u^4+b u^3+c u^2+P1 u+P0,
    Q=(a-1)u^2+(b-2pa+4p+1)u+d,
    M=M4 u^4+M3 u^3+M2 u^2+M1 u+M0,
    R=M4 u^2+(M3-2pM4)u+e,
    H=N t^2+P t+Q,  L=M t+R.                         (1)

Every lower row in (1) is retained. In particular, the `4p+1` term in
`Q` comes from the original leading polynomial's fifth coefficient.
We absorb the additive constant `e` into the generic value of `F`.

The inherited mechanism is the original-fibre formal inverse and the
[finite first-jet gate](planar_jc48_sep08_shared_roots.md). The closest
corrected near miss is to search for a primitive of low pole degree
before exhausting exact coefficient consequences. The live concepts are
full global sections, moving-root residues, the algebraic leading field,
Ferrari pair sums, and the distinction between formal necessities and
rational algebraization. The source-to-target operation below retains
the complete original `F` and all induced rows; coefficient extraction
loses algebraization of a compatible infinite collection. No sufficiency
of the hierarchy is asserted.

## 2. Active jets and the two exhaustive strata

At the finite simple point all normalized generic relative forms are
regular. At the original infinity point of multiplicity two they are also
regular for every `M`: the actual chart has

    x=1/r, t=-r^2-r^4b_D, omega=r^2 dr wedge db_D.

An M-unit or normal unit is already regular. In the remaining shared
case, the complete generic tangent quartic has four simple nonzero roots;
its unweighted forms are logarithmic, and multiplication by `r^2` makes
each regular. This is the complete low-multiplicity classification in the
[shared-root theorem](planar_jc48_sep08_shared_roots.md), including its
explicit infinity qualification.

At the finite fivefold point an M-unit is unbalanced (`5!=3j`) and regular
by [the unit-M analysis](planar_jc48_sep08_quartic_common_root.md); a normal
unit is regular. If a rational mate exists, the fivefold point must
therefore be active, since otherwise a rational primitive would have no
poles on each compact generic component although the relative form is
nonzero there. The finite first-jet theorem requires

    P0=P1=M0=M1=0.                                  (2)

No connectedness assumption is needed for this preliminary argument.

The complete inverse coefficient `T2=-M/(4N)` implies that `M du/N`
is rational-exact, either directly or by normalized trace from its
quadratic field. Its residue at the finite simple point requires
`M(1)=0`. Hence (1)--(2) become

    P=u^2(a u^2+b u+c),
    M=u^2(u-1)(lambda u+mu),
    R=lambda u^2+(mu-lambda-2p lambda)u.              (3)

Put `w=u^2t`, with actual volume `u^-2 du wedge dw`. The first two
coefficients of the complete `F` are

    f(w)=(cw+d)^2-mu w,
    g(w)=2(cw+d)(-w^2+bw+b-2pa+4p+1)
          +(mu-lambda)(w+1)-2p lambda.               (4)

At every simple root of the original `f(w)=zeta`, rational exactness
requires the moving-root residue to vanish. Thus `(g/f')'=0` as an
identity, not merely at one fixed level. If `c!=0`, the cubic coefficient
`-2c` of `g` cannot be a constant multiple of the linear `f'`. Therefore
`c=0`.

If `mu!=0`, `f` is now linear. Exactness forces `g` constant; its
quadratic and linear coefficients give `d=0` and `mu=lambda`. This is
stratum **n2**:

    M=lambda u^2(u^2-1), R=lambda(u^2-2pu),
    c=d=0, lambda!=0.                               (5)

If `mu=0`, nonconstant `L` requires `lambda!=0`. This is stratum **n3**:

    M=lambda u^3(u-1), R=lambda u(u-1-2p),
    c=0, lambda!=0.                                 (6)

These alternatives exhaust every nonconstant `L`; they are not selected
parameter families. If `L` is constant, the polynomial conclusion is
already immediate from the factor `2H` in its Jacobian.

## 3. Every nonzero-index inverse coefficient is an exact differential

This section derives the additional coefficient tests directly, keeping
the complete original fibre. Let

    K=C(u)(sqrt(N)), v=F^(-1/4),
    t=T(u,v)=sum_(k>=-1) T_k(u) v^k,
    T_-1=1/sqrt(N).

The formal implicit-function theorem gives this unique branch in `K((v))`.
A rational mate composes to a Laurent series `Gtilde(u,v)` over the same
field. Differentiating `F(u,T)=v^-4` and using `J(F,G)=1` gives

    partial_u Gtilde |_v = (1/4)v^5 partial_v T,
    partial_u [v^(k+4)]Gtilde = (k/4) T_k.           (7)

Consequently `T_k du` is exact in `K` for every `k!=0`. For a rational
coefficient, normalized field trace gives exactness already in `C(u)`.
For an algebraic coefficient, residues on the normalized quadratic curve
must vanish. No implication from all these necessities to existence of a
rational mate is used.

Center the quadratic by `t=y-P/(2N)` and write

    D0=Q-P^2/(4N), E0=R-MP/(2N),
    F=(N y^2+D0)^2+M y+E0.

Set `W=sqrt(N) v y`. Its exact equation is

    W^4+2D0 v^2 W^2+(M/sqrt(N))v^3 W
        +(D0^2+E0)v^4=1.                            (8)

Coefficient comparison through order six gives

    T1=-D0/(2sqrt(N)),
    T2=-M/(4N),
    T5=-(2D0^3+4D0E0+M^2/N)/(32sqrt(N)).             (9)

The translation affects only `T0`, which is not constrained by (7).

The opposite-root pair yields the rational coefficients particularly
efficiently. If `z=v^4`, write

    T_(4j+2)=-M/(4N) * [z^j]V(z),
    (1-E0z)V^2 + D0 M^2 z^2 V^4/(4N)
                  +M^4 z^3 V^6/(64N^2)=1,
    V(0)=1.                                        (10)

One may obtain (10) by factoring the depressed quartic into the two
quadratics containing opposite inverse roots; it is the same exact
pair-sum identity proved in Section 5 of
[the all-finite (5,3) theorem](planar_jc48_sep08_five_three.md).
For completeness the needed coefficients are

    V1=E0/2,
    V2=(3E0^2-D0 M^2/N)/8,
    V3=(40E0^3N^2-40D0E0M^2N-M^4)/(128N^2),
    V4=7(2D0^2M^4-20D0E0^2M^2N+10E0^4N^2-E0M^4)
          /(256N^2).                               (11)

The exact source independently substitutes (11) in (10), and derives
(9) directly from (8). This controls normalization and pair-sum factors
through `T18` without inferring them from a finite sample of polynomials.

## 4. Stratum n2: four rational residue rows give a contradiction

In (5), the `T6` residue at zero gives

    b=2-4p-2a.

Put `A=a-2`, so the complete retained coefficients are

    P=u^3[(A+2)u-2-4p-2A],
    Q=(A+1)u^2+(-1-2A-2pA-4p)u.

Define `qA=A^2+6pA+12p^2`. Exact substitution in (11) gives

    Res_0 T10 = -lambda^3 qA/32.                     (12)

Modulo `qA=0`, the next residue is

    Res_0 T14 = -5lambda^4[lambda+32p^2(A+3p)]/512.  (13)

Thus `lambda=-32p^2(A+3p)`. In particular `p!=0`, since `lambda!=0`.
After these two relations the next residue reduces to

    Res_0 T18 = -77lambda^5 p^3(A+4p)/64.            (14)

It forces `A=-4p`. But substitution in `qA` gives `4p^2`, which is
nonzero. This contradicts (12). All reductions are polynomial identities
with denominators only nonzero rational constants and the explicitly
nonzero `lambda`; no exceptional coefficient stratum was divided away.
Hence the whole n2 stratum has no rational mate.

The earlier rows alone do not suffice. For example `p=A=0` satisfies the
normalized `T6,T10` equations, whereas `T14` is nonzero for every nonzero
`lambda`. The exact source retains this staged hostile rather than
mistaking a partial residue test for closure.

## 5. Stratum n3: two algebraic infinity residues finish the reduction

In (6), the rational residues give

    Res_0 T6=lambda^2(b+4p+2)/16,
    Res_0 T10=-d lambda^3/32.

Therefore `b=-4p-2`, `d=0`. Again put `A=a-2`; then

    D0=-u(A^2u^2-8Ap+16p^2)/[4(u-1)],
    E0=-A lambda u^2/2.                             (15)

The leading quadratic field is the genus-zero curve

    y^2=u(u-1),  sqrt(N)=u^2 y.

Using (9),

    T1=(A^2u^2-8Ap+16p^2)/(8y^3).                   (16)

There are two infinity points, with `y/u` tending respectively to
`epsilon=+1,-1`. In their actual local parameter `z=1/u`, the residue
of `T1 du` is `-epsilon A^2/8`. Exactness forces `A=0`.

Now `D0=-4p^2u/(u-1)` is bounded at those infinities and `E0=0`.
Moreover

    M^2/N=lambda^2u(u-1).

In `T5` the terms involving `D0` are `O(u^-3)`, while the last term is

    -lambda^2/(32 epsilon u)+O(u^-2).

Thus the actual residues of `T5 du` are
`epsilon lambda^2/32`, nonzero. This contradicts (7) for `k=5`.
There is no unproved genus or primitive-space premise here: both residues
are computed on the explicit normalized leading field, using its two
infinity points. The entire n3 stratum has no rational mate.

This also shows why continuing only the rational even-index rows can
miss a short obstruction. At the tuned coefficients (15), `T14,T18`
have zero rational residue at zero; the algebraic `T1,T5` rows supply the
missing information.

## 6. Scope of the complete conclusion

Sections 2--5 exhaust all nonconstant `L` when the double root is infinity.
A rational mate forces `L` constant, and then
`J(H^2+constant,G)=2H J(H,G)` cannot equal one for a polynomial `G`.
At every other placement the finite-double leading residue already
excludes rational mates. The all-location polynomial statement follows.

The proof retains every induced global coefficient and every location.
It imposes no bound on mate degree and makes no numerical or sampled
parameter inference. Conversely the inverse hierarchy is used only in
its necessary direction. Formal coefficient primitives need not
algebraize to a rational function, and no general exactness equivalence
is claimed.

## 7. Sharp rational boundary and exact verification universe

For every `p`, let `u=x-p` and

    Z=u^2t+1,
    H=u(u-1)Z^2+kappa,
    G_H=-(1+2u)/(3u^2 Z).

These satisfy `J(H,G_H)=1` exactly. Their complete global coefficients
are

    N=u^5(u-1),  P=2u^3(u-1), Q=u(u-1)+kappa.

They obey the original filtration for every `p`; no surface translation
has been assumed. Hence

    F=H^2+e,  G=G_H/(2H),  J(F,G)=1

is a rational-mate family in this exact partition. It is the correct
boundary to the rational conclusion and cannot be polynomial because
of the factor `2H`.

The reproducible source is
[planar_jc48_sep08_five_two_one.py](../../04-computation/planar_jc48_sep08_five_two_one.py).
Its named universe comprises the complete all-p section equations,
both first-row alternatives, direct formal inverse coefficients through
order six, the exact Ferrari recurrence through order four, every
symbolic residue and ideal reduction in Sections 4--5, the two infinity
signs, and the sharp constant-L family. It separately checks the full
original chart and field-normalization factors. Its gates are
always-active under `-O`.

The local/global necessity and the all-parameter implication are the
analytic proof, not an extrapolation from the gate count. The source and output are frozen; independent review remains pending.


Reproduce from the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_five_two_one.py
python3 -B -O 04-computation/planar_jc48_sep08_five_two_one.py
```

Both executions pass **59 always-active exact gates** and reproduce all
373 bytes of [the frozen output](planar_jc48_sep08_five_two_one.out).

| Artifact | SHA-256 |
| --- | --- |
| Source, 7,045 bytes | `4735e3d4c91c86218e192c4aa88a4b7a9f082461622b3141f0870fd09d42ac2e` |
| Frozen output and both replays, 373 bytes | `bdccd417a403061e0362efe8e1a5b8ef711f2cd81f480bbff1bac60a45f38dca` |
| Semantic record | `9b6d82638dc683f69afbf67aa90336489cff490850e6cc69edd401731513ac10` |


## Accepted independent audit

The [complete independent audit](planar_jc48_sep08_five_two_one_audit.md)
accepts the analytic proof, full coefficient scope and actual source
coordinates. Both normal and optimized replays reproduce all59
gates and the frozen output. The candidate is now in the proved graph.
