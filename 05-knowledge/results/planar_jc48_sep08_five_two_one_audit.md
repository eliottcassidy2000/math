# Independent audit of the complete binary (5,2,1) obstruction

**Status: PROVED / INDEPENDENT AUDIT PASS.**
The entire analytic proof, frozen source, and independent normal and
optimized replays of [the primary](planar_jc48_sep08_five_two_one.md)
are accepted. Separate coefficient recurrences, rational-series
residue extraction, and explicit monic quadratic reductions reconstruct
the decisive higher-residue identities. No correction was needed.
Producer artifacts were not edited; root owns status promotion and
integration.

## 1. Exact statement, locations, and coefficient universe

The theorem concerns actual global `H in L2`, `L in L1` on the
fixed DG surface, `F=H^2+L`, with leading binary multiplicities
five, two, one at three distinct points. It excludes every
polynomial mate, without a bound on its degree. If the double
point is finite, leading exactness excludes every rational mate.
If it is the original infinity, a rational mate forces `L`
constant. The latter boundary is sharp, as shown by an actual
all-position constant-`L` rational family.

The finite double-root argument uses a nonzero simple pole of
`du/sqrt(N)` on the actual leading field. It remains applicable
if the fivefold or simple root is at infinity. With the double
root at infinity, actual scaling normalizes the nonzero separation
and leading scalar; the remaining arbitrary finite position is
retained in `u=x-p`. This coordinate is not asserted to act on
the compactification.

For `N=u^5(u-1)`, the complete original global rows are

    P=a u^4+b u^3+c u^2+P1u+P0,
    Q=(a-1)u^2+(b-2pa+4p+1)u+d,
    M=sum_(i=0)^4 M_i u^i,
    R=M4u^2+(M3-2pM4)u+e.

I checked the original-coordinate high coefficients and the
full polynomial expressions in the second chart. In particular,
the `4p+1` term is induced by the actual leading polynomial,
and the arbitrary `P0,P1,M0,M1` remain until the local theorem
forces their vanishing. These are complete section spaces, not
a selected coefficient family. The additive lower constant is
absorbed only into the original generic value.

## 2. Local entry and exhaustive moving-root alternatives

The finite simple root contributes regular relative forms. At
original infinity the multiplicity is two and the actual
canonical multiplier is `r^2`. Normal-unit and M-unit cases
are regular; the remaining shared tangent quartic has four
generic nonzero simple determinations whose unweighted
logarithmic order becomes regular after multiplication by
`r^2`. These are the already-proved complete low-multiplicity
local alternatives, including the original infinity chart.

At the finite fivefold root, the M-unit case is unbalanced,
and the normal-unit case is also regular. Unless that root
is active, a rational primitive would have no pole on any
compact normalized generic component and would be constant,
contrary to its nonzero differential. The shared first-jet
gate therefore forces `P0=P1=M0=M1=0`. This preliminary
argument is componentwise and does not assume connectedness.

The rational coefficient `T2=-M/(4N)` is exact after normalized
trace, so its simple-root residue requires `M(1)=0`. The entire
remaining lower coefficient is therefore

    M=u^2(u-1)(lambda u+mu),
    R=lambda u^2+(mu-lambda-2p lambda)u.

There is no unexamined fourth-degree lower polynomial outside
this factorization. In the actual chart `w=u^2t`, the complete
constant and first rows of `F` are the primary's `f(w),g(w)`.
The moving-root residue is proportional to `g f''-g'f'` at
a simple root of the original `f(w)=zeta`. Since that root
is generic, its vanishing forces the rational identity
`(g/f')'=0`, not only a relation on a sampled level.

If `c!=0`, the cubic term `-2c w^3` in `g` cannot be a
constant multiple of the linear `f'`, so `c=0`. If `mu!=0`,
the now-linear `f` forces `g` constant; its quadratic and
linear coefficients give `d=0,mu=lambda!=0`. This is the
complete n2 stratum. If `mu=0`, nonconstant `L` requires
`lambda!=0`, giving exactly n3. If both lower parameters
vanish, `L` is constant and has already left the nonconstant
branch. No nonconstant parameter case is discarded.

## 3. Universal inverse equations and exactness directions

For the actual inverse `F(u,T)=v^-4`, the formal chain rule
with `J(F,G)=1` gives

    Gtilde_u=(1/4)v^5 T_v.

For every existing coefficient index other than zero, the
coefficient `T_k du` must therefore be exact in
`C(u)(sqrt(N))`. Index zero has zero multiplier and is not
constrained. A rational coefficient becomes exact already
in `C(u)` by normalized field trace; an algebraic coefficient
has zero residues on the normalized leading field. No
converse algebraization assertion is made.

Centering by `t=y-P/(2N)` gives

    D0=Q-P^2/(4N), E0=R-MP/(2N),
    W^4+2D0v^2W^2+(M/sqrt(N))v^3W
                  +(D0^2+E0)v^4=1,
    W(0)=1.

The coefficient of each new `W` coefficient is four. I
independently reconstructed this recurrence through the
needed order and obtained exactly

    T1=-D0/(2sqrt(N)), T2=-M/(4N),
    T5=-(2D0^3+4D0E0+M^2/N)/(32sqrt(N)).

For the rational even coefficients, factor the centered
monic quartic into the two opposite-pair quadratics. The
pair-sum square `X` satisfies the universal resolvent
`X(B+X)^2-C^2-4DX=0`. Substitution of the opposite sum
`-M v^2V(v^4)/(2N)` gives precisely

    (1-E0z)V^2+D0M^2z^2V^4/(4N)
                         +M^4z^3V^6/(64N^2)=1.

Each new `V` coefficient has coefficient two. The primary's
four displayed `V` coefficients were checked by direct
substitution before any family specialization. Their
normalizations yield `T6,T10,T14,T18` after multiplication
by `-M/(4N)`.

No division by a potentially zero coefficient restricts the
final formula. The universal identities can be derived on
the nonzero generic coefficient field and then cleared as
polynomial coefficient identities; the zero lower row is
included by specialization. In the two nonconstant strata
actually used, `lambda` is explicitly nonzero. The entire
original `F`, including its induced `Q,R`, is retained.

## 4. Independent n2 residue reconstruction and contradiction

The first residue in the unrestricted active family is

    Res T6=lambda(2a mu+b lambda+4p lambda
                              +2lambda-4mu)/16.

At `mu=lambda` it forces `b=2-4p-2a`. Put `A=a-2`.
I independently extracted the subsequent residues and obtained

    Res T10=-lambda^3(A^2+6pA+12p^2)/32.

Writing `qA=A^2+6pA+12p^2`, reduction modulo this monic
quadratic gives

    Res T14=-5lambda^4[lambda+32p^2(A+3p)]/512.

Thus `lambda=-32p^2(A+3p)`. Since `lambda!=0`, this
already excludes `p=0`. After this substitution and the
same monic reduction, the next residue is

    Res T18=-77lambda^5p^3(A+4p)/64.

It forces `A=-4p`, whereas `qA(-4p)=4p^2!=0`.
This proves the contradiction for every complex parameter,
including every boundary value that might have disappeared
under an unjustified division.

For independence, my residue implementation did not call
the producer or SymPy's residue function. It canceled the
literal rational expression, removed the denominator's
power of `u`, and generated the needed reciprocal Taylor
coefficients by their elementary linear recurrence. All
remaining denominator constant terms were checked to be
nonzero rational constants. I also implemented explicit
descending substitution of `A^2=-6pA-12p^2`, without using
the producer's polynomial remainder function. This recovers
the four identities above. All resulting finite residues
have constant-only denominators; neither `p` nor a hidden
linear factor was inverted.

The staged hostile `p=A=0` passes the earlier tuned `T6,T10`
conditions but fails `T14` for every nonzero `lambda`. It is
properly retained as a warning about partial tests, not a
purported surviving rational mate.

## 5. n3 and the two actual algebraic infinity residues

For `mu=0`, the independently reconstructed rational
residues are

    Res T6=lambda^2(b+4p+2)/16,
    Res T10=-d lambda^3/32.

Thus `b=-4p-2,d=0`. Direct substitution gives

    D0=-u(A^2u^2-8Ap+16p^2)/[4(u-1)],
    E0=-A lambda u^2/2.

The actual leading field is `y^2=u(u-1)` with
`sqrt(N)=u^2y`. There are two unramified infinity points,
and `y/u` tends to `epsilon=+1` or `-1` respectively.
The coefficient

    T1=(A^2u^2-8Ap+16p^2)/(8y^3)

has residue `-epsilon A^2/8` in the actual local
coordinate `z=1/u`. Both signs agree with direct Laurent
expansion, including the minus sign from `du=-z^-2dz`.
Exactness forces `A=0` over `C`.

Then `D0=-4p^2u/(u-1)` is bounded and `E0=0`.
The terms involving `D0` in `T5` are `O(u^-3)`, while

    M^2/N=lambda^2u(u-1),
    T5=-lambda^2/(32epsilon u)+O(u^-2).

The two residues are consequently `epsilon lambda^2/32`,
nonzero because `lambda!=0`. This contradiction uses
the original inverse coefficient on the normalized
leading field; it does not require a compact generic-
fibre genus formula or a chosen primitive space.

The later rational even-row residues can vanish in this
tuning. The proof correctly retains the algebraic odd
rows, rather than equating rational-coefficient exactness
with the full inverse hierarchy.

## 6. Actual conclusion and rational sharpness

The n2/n3 alternatives exhaust every nonconstant `L`
with the double point at infinity, so a rational mate
forces `L` constant there. For constant `L`, the
polynomial bracket `J(H^2+constant,G)` has nonconstant
factor `2H` and cannot equal one. The finite-double
leading residue covers every remaining root placement.
This proves the full all-location polynomial statement
with no mate-degree bound.

The sharp family is literal and correctly normalized:

    Z=u^2t+1, H=u(u-1)Z^2+kappa,
    G_H=-(1+2u)/(3u^2Z).

I checked its exact Jacobian and the full original
globality. Its leading coefficient is `u^5(u-1)`;
`P=2u^3(u-1)` and `Q=u(u-1)+kappa` satisfy the
induced rows for every original `p`. Therefore
`H^2+e` has rational mate `G_H/(2H)`. The denominator
is essential and preserves the correct limit of the
rational conclusion. No surface-translation or
arbitrary-Keller entry claim is introduced.

## 7. Independent replays and frozen pins

The entire source was read, including direct inverse
recursion, literal Ferrari substitution, the full global
chart, scalar denominator checks, all residue reductions,
both infinity signs, and the staged hostiles. Its
59 gates are explicit exception checks and remain
active under Python optimization. The coefficient
universe is symbolic; no numerical parameter bank
is substituted for the all-parameter proof.

The separate referee source reconstructs the universal
`V` coefficients rather than copying their formulas,
extracts finite residues by a reciprocal series, and
uses explicit monic long reduction. Its 16 primary
checks recover all n2 residues and the odd inverse
coefficients. Additional direct checks recover the
untuned `T6`, n3 `T10`, and constant-only residue
denominators.

Independent full commands were

```sh
python3 -B 04-computation/planar_jc48_sep08_five_two_one.py > /tmp/five-two-one-audit-normal.out
python3 -B -O 04-computation/planar_jc48_sep08_five_two_one.py > /tmp/five-two-one-audit-optimized.out
python3 -B /tmp/five_two_one_independent.py > /tmp/five-two-one-independent.out
```

Both producer replays completed successfully with
**59 always-active gates** and reproduce the frozen
373-byte output byte-for-byte. The full primary,
source and final replay pins were checked at acceptance.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Primary before promotion | 12,980 | `c37a4827199bc4b3dfa649c50ac90b369b8911ecdc4e3acc92a4f3becdb866f4` |
| Frozen source | 7,045 | `4735e3d4c91c86218e192c4aa88a4b7a9f082461622b3141f0870fd09d42ac2e` |
| Frozen output and each independent replay | 373 | `bdccd417a403061e0362efe8e1a5b8ef711f2cd81f480bbff1bac60a45f38dca` |
| Independent recurrence/reciprocal source | 2,604 | `19c52274d101e939647814a2d9a8d72794352ba7a72b9a9d7ec267666879c335` |
| Independent output | 198 | `372e32f1a05745816754e127751f9ed08fe31d2952da09ed199fd24ffd94495b` |

The producer semantic digest is
`9b6d82638dc683f69afbf67aa90336489cff490850e6cc69edd401731513ac10`.
The independent control digest is
`01115e68fc092efdf74268d99f19e3cdd71b28c666f030505df7ddc5c85146e0`.
No correction or unresolved mathematical qualification
remains. The accepted result has precisely the stated
all-location polynomial scope and constant-`L`
rational boundary.
