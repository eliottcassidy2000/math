# Independent audit: the universal degree-six constant-boundary gate

**Status: PASS — analytic proof, exact source, normal/optimized replays.**

Auditor: root, independently of the producer `three_ray_geometry`.
The [primary theorem](planar_jc48_sep08_degree_six_constant_d.md) was
audited before status promotion. This audit does not enlarge its scope.

## 1. Accepted statement and exact hypotheses

For the specified DG surface with original source volume `dx wedge dt`,
let `F=H^2+L`, `H in L2`, `L in L1`, and let the leading polynomial
`N` of `H` have degree exactly six. A rational constant-Jacobian mate
forces both `H|D` and `L|D` constant. No assumption about finite roots,
their jets, genus, or the degree or pole divisor of the mate is used.
The short complete `(5,2,1)` consequence is also accepted.

The theorem is necessary, not sufficient for a rational mate. It is not
a theorem about arbitrary quartics with degree-six leading square root,
and supplies no reduction of arbitrary planar Keller maps to this surface.
Constant scaling makes `N` monic and preserves all relevant predicates.

## 2. Independent analytic reconstruction

I reconstructed the original, unshifted global section rows directly from
the proved numerator boxes. They are

    N=x^6+n5*x^5+...+n0,
    P=(A+2)*x^4+b*x^3+c*x^2+P1*x+P0,
    Q=(A+1)*x^2+(b-n5)*x+d,
    M=lambda*x^4+mu*x^3+nu*x^2+M1*x+M0,
    R=lambda*x^2+mu*x+e.

These have all eighteen independent coefficient parameters. The change
`u=x-p` gives the primary's full shifted family, including `-2pA` in
the linear `Q` row and `-2p lambda` in `R`. Translation is used as a
source coordinate only; it is not spent as a projective automorphism.
The actual boundary slopes are `-A` and `-lambda`.

Polynomial division in the original `x`, independently of the producer's
normalized local numerators, gives

    D0=Q-P^2/(4N)=-A^2*x^2/4+O(x),
    A=0 => D0=O(1), E0=R-MP/(2N)=O(x),
    M^2/N=lambda^2*x^2+O(x).

The vanishing of the entire linear term of `D0` after `A=0` is essential
and was checked symbolically, with all low coefficients retained.

For an independent inverse-coefficient path, write the centered equation
in the scaled fibre variable as

    phi(z)=1+2D*z^2+B*z^3+(D^2+E)*z^4.

Lagrange inversion gives its inverse coefficient at index `j>0` as

    -(1/j) [z^(j+1)] phi(z)^(j/4).

Direct binomial extraction, without the producer's implicit recurrence,
gives `-D/2`, `-B/4`, and `-(2D^3+4DE+B^2)/32` at indices 1,2,5.
Substituting `B=M/kappa` and dividing by `kappa` recovers all three
displayed `T_j`, including the coefficient and sign of `M^2/N` in `T5`.

The formal inverse exists over the actual field `K=C(u)(sqrt(N))`:
its scaled derivative at the selected constant root is `4 sqrt(N)`,
nonzero in the field. Rational substitution is injective by the unique
lowest valuation of a polynomial's highest fibre-degree term. The chain
rule with the original Jacobian gives

    partial_u Gtilde=(1/4)v^5*T_v,
    partial_u [v^(j+4)]Gtilde=(j/4)T_j.

Thus the required coefficient differentials are exact. No convergence
or algebraization of an arbitrary formal solution is claimed. The
index-zero exception is correctly excluded.

With `z=1/u`, the local square root has form
`epsilon*z^-3*(1+O(z))`. For nonsquare `N` there are two unramified
infinity points; for square `N` the chosen field is `C(u)` with one.
The respective residues are

    Res(T1 du)=-epsilon*A^2/8,
    A=0 => Res(T5 du)=epsilon*lambda^2/32.

Exactness forces `A=lambda=0`. Repeated finite roots change neither
the local unit calculation nor this field distinction. The proof uses
leading-field coefficient differentials, so no actual quartic-fibre
genus or irreducibility premise is needed. The source volume was used
before coefficient extraction; this is not an unweighted replacement
of the second-chart volume.

## 3. The shorter (5,2,1) consequence

I checked its logical order separately. A finite double root has the
nonzero leading logarithm. For a double root at original infinity,
the universal theorem is applied before the active finite-root suppliers.
The M-unit and normal-unit cases at the finite fivefold point are regular;
if it were inactive, every generic compact component would support a
nonzero exact holomorphic differential, a contradiction. The first-jet
and `T2` simple-root conditions give the complete displayed active rows.

In the actual chart `w=u^2 t`, the original form is
`u^-2 du wedge dw`. Expanding the implicit moving root of
`F=f(w)+u g(w)+...` gives residue `-(g/f')'/f'` whenever `f` is
nonconstant. All generic values are required to turn vanishing at the
roots into an identity. A cubic `g` cannot be a multiple of linear `f'`,
forcing `c=0`; with nonconstant `L`, the ensuing linear `f` first forces
`d=0`, then forces the nonzero coefficient `mu` to vanish. The omitted
case `c=mu=0` already has constant `L` and needs no moving-root argument.
The final polynomial obstruction is the nonunit factor `2H`.

## 4. Controls and source audit

I read the full standalone source. Its runtime checks remain active
under `-O`, the complete symbolic family precedes finite controls, and
there is no hidden root filter, parameter census, mate-degree cap, or
producer-file import used to establish the result.

The named positive control `H=u^2(1+u^2t)^2+(1+u^2t)`,
`L=u(1+u^2t)`, `G=1/(2u(1+u^2t))-H` has nonconstant `L` and an actual
rational mate, preserving the sharp constant-boundary limit. The `T1`
only hostile has a nonzero `T5` residue. The constant-boundary finite
double hostile has a nonzero residue at zero. The nonglobal pure
monomial has the exact mate `1/(8u^11t^3)` and a genuine boundary pole,
so global lower-row constraints cannot be removed. The `(5,2,1)`
constant-`L` family is global for every original finite shift and has
the stated rational mate. These controls support, rather than replace,
the universal coefficient proof.

Independent complete replays:

    python3 -B 04-computation/planar_jc48_sep08_degree_six_constant_d.py
    python3 -B -O 04-computation/planar_jc48_sep08_degree_six_constant_d.py

Both exit successfully, reproduce **90 gates**, and match the frozen
455-byte output byte for byte. The additional Lagrange extraction and
original-coordinate polynomial divisions described above also pass.

Frozen producer pins before status promotion:

    source bytes 9164
    source SHA256 87d48d9df51f2a4601132e2094d46cccc13d783e16e977a2bbf5ceb350365727
    output bytes 455
    output SHA256 061bccb24949d5efe2381ba56b9dc8c7c85500d02eb35f5e63975dad3be28863
    primary bytes 16051
    primary SHA256 c1953d97ffd0e0b5e4ecdb057e33c9579d6c5ea2fbf9249b9b7d02dc6a078969

No correction remains. Parent-owned promotion may change the primary
status and add this audit link; the source and output stay frozen.
