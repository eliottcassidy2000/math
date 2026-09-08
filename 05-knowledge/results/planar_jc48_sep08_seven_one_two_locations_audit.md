# Independent audit of the two sevenfold-finite (7,1) placements

**Status: INDEPENDENT FULL ANALYTIC / SOURCE / REPLAY AUDIT PASS.**
The [primary proof](planar_jc48_sep08_seven_one_two_locations.md), its
standalone producer, and its frozen output pass. This sidecar records
acceptance of the complete polynomial-mate exclusion when the sevenfold
point is finite. It does not promote the primary's status, alter any
producer artifact, exclude every rational mate, or settle JC(2).

## 1. Accepted statement and proof boundary

The actual surface and volume are

    W=(P1_x x P1_z) minus {z=x^2},
    t=1/(z-x^2), x=1/r, t=-r^2-r^4 b_D,
    omega=dx wedge dt=r^2 dr wedge db_D.

For `H in L2`, `L in L1`, `deg_t H=2`, put `F=H^2+L`.
If the complete binary octic leading section has two distinct zeros of
multiplicities seven and one, with the sevenfold zero finite, the
primary proves that no original polynomial `G(x,t)` has nonzero constant
Jacobian with `F`. The degree of `G` is unrestricted. Multiplying `G` by
a nonzero constant reduces its Jacobian to one.

The proof covers both an arbitrary finite simple point and the original
point at infinity. It retains the arbitrary original finite position
`p` throughout. The placement with the sevenfold point at infinity is
outside this theorem. In particular, “OPEN” in the primary's closing
scope means that this argument has not closed that placement.

The conclusion is polynomial. Some necessary conditions below hold for
rational mates, and the final residual families have rational
obstructions, but their entry uses original polynomial support. The two
explicit rational positive families satisfy the exact same leading
partitions and are global `L2` sections for every `p`. They rule out
silently strengthening the final theorem to rational mates.

The main proof move is to impose the polynomial first-jet constraint
before changing rational coordinates, then retain the actual field and
the actual source volume while computing generic residues. No bounded
search in the degree of a mate is used. No generic irreducibility
assumption substitutes for the field or component arguments below.

## 2. Inherited entry and completeness of the original rows

I reread the relevant statements in the following current proved files:

* [DG filtration](planar_jc48_sep08_dg_quadratic.md), especially its full
  quadratic coefficient space and actual second chart;
* [M-unit local analysis](planar_jc48_sep08_quartic_common_root.md),
  including all unbalanced normal orders and the infinity multiplier;
* [shared-root first jets](planar_jc48_sep08_shared_roots.md), including
  simple-root and normal-unit regularity and the finite logarithmic test;
* [leading inverse exactness](planar_jc48_sep08_leading_exactness.md),
  including the radical coefficient field, faithful substitution, and
  coefficientwise primitive equation.

Each is promoted and independently audited on the inspected worktree.
The primary reproduces the needed polynomial weight argument directly;
its separate all-m extension is not a provisional dependency.

At a finite sevenfold point, M-unit regularity is indeed unbalanced:
`7=3j` has no integral solution, and the cases `j=0`, finite positive
`j`, and infinite normal order are all included by the local supplier.
At a shared root with nonzero normal derivative the form is regular on
every normalized branch. At the other, simple root the M-unit and shared
analyses are regular as well. If that root is the original infinity,
the additional canonical factor only increases differential orders.

Thus an inactive sevenfold point would leave the relative form
holomorphic on every compact normalized generic component. Exactness
would make its primitive constant, whereas the actual form is nonzero
on the component's source part. Generic components cannot be contained
in the fixed boundary divisors; a generic level excludes their fixed
values. This gives the stated componentwise contradiction without
assuming geometric generic integrality. At an active shared sevenfold
point, the finite first-jet supplier gives a nonzero residue unless

    P(0)=P'(0)=M(0)=M'(0)=0.

This use is finite; no logarithmic obstruction is transported through
the canonical zero at infinity.

The full fixed-leading `L2` space has six free coefficients. One can
recover this directly from the fifteen-dimensional filtration: its map
to the nine coefficients of the leading binary octic is surjective,
with a six-dimensional kernel. Alternatively, the three numerator rows

    Q, P-2Qx^2, N-Px^2+Qx^4

must separately have degree at most four. Solving their coefficients
with `x=u+p` gives exactly the primary's two displayed six-parameter
spaces, with free coordinates `a,b,c2,c1,c0,d`. Their independence is
visible in the five free coefficients of `P` through degree four and
the final constant of `Q`; there is no special-p rank loss. The
analogous two linear rows give

    M=sum_(i=0)^4 m_i u^i,
    R=m4 u^2+(m3-2p m4)u+e.

The frozen producer checks the complete symbolic post-entry carriers,
not just the `p=0` specialization. Completeness of the pre-entry spaces
is provided by the filtration and the displayed coefficient solution,
rather than by a numerical rank claim in the producer.

For two finite points, scaling `x` and `z` by compatible nonzero powers
preserves the graph and normalizes their separation; the finite point
becomes a still-arbitrary `p`. A constant rescaling of `H` and the
corresponding square rescaling of `F,L` normalizes its leading scalar.
These have nonzero constant source/target Jacobians. Writing `u=x-p`
subsequently is an ordinary original coordinate change, not a claimed
automorphism of the compactification. The source does not discard the
resulting p-dependent coefficients.

## 3. Polynomial entry and the inverse coefficient at the simple root

Under the active rows, put `w=u^2t`. The image of an original polynomial
mate is a finite Laurent sum `sum u^ell g_ell(w)`; a nonzero negative
row is divisible by `w`. The exact bracket is `u^2 J_(u,w)`.

For nonconstant `f0(w)`, the least negative mate row gives the nonzero
term `-ell f0' g_ell u^(ell+1)`. Every higher row and every positive-u
term of `F` has strictly higher weight, so cancellation at that first
weight is impossible. Only `ell=-1` can yield weight zero; the factor
`w` in `g_ell` prevents a constant coefficient. If the mate has no
negative rows, all nonzero bracket terms have positive weight. This
includes `ell=0`, where the nominal leading coefficient vanishes.

Consequently `f0` is constant. In the square-prefix expression this
forces `P2=M2=0`. If `f1` is nonconstant of degree `s>=1`, a least
negative mate coefficient of degree `j>=1` gives

    f1 g_ell' - ell f1' g_ell,

whose leading coefficient has nonzero factor `j-ell s` and degree
`s+j-1>=1`. Again only `ell=-2` could have weight zero, and its
coefficient cannot be constant. Nonnegative mate rows give positive
weight. Therefore `f1=beta` is constant. If `beta=0`, substituting
back to the original ring makes `F-f0` divisible by `u^2`, so both
partials vanish on the source line `u=0`. This establishes the required
nonzero-beta condition, not merely the constancy of a first row.

These arguments yield precisely

    P2=M2=0, M3=-2P3 Q0, 2Q0 Q1+R1 != 0.

They use the complete polynomial image constraint, without a degree
bound. They do not apply to a general rational mate.

I also reconstructed the additional all-finite simple-root constraint.
Center the original quadratic `H` by writing

    H=N s^2+D, F=(N s^2+D)^2+M s+E.

With `alpha^2=N` and `v=F^(-1/4)`, formal inversion begins

    s=alpha^-1 v^-1-D/(2alpha) v-M/(4N) v^2+O(v^3).

Cancellation of the four Laurent coefficients through order `v^-1`
recovers the displayed `v^2` coefficient directly. The original
centering is independent of v and does not change it. The exact
formal equation `Gtilde_u=(1/4)v^5 T_v` then makes `T2 du` exact in
the chosen radical coefficient field. Since `T2=-M/(4N)` already
lies in `C(u)`, normalized field trace descends the primitive to
`C(u)`. Equivalently, all local ramified residues must vanish.
The simple pole at `u=1` therefore forces `M(1)=0`. This remains
valid for a proper-degree radical extension.

Combining the full rows, first jets and this simple-root condition
gives exactly (7) in the primary, with `lambda=2bd`. If lambda is
nonzero, both b and d are nonzero. The other placement gives exactly
(8), with `mu=-2bd`. No division by p or a possibly vanishing beta
expression is used later. A zero beta would already exclude the
polynomial mate.

## 4. All-finite actual field and its generic residue obstruction

I independently expanded the original carriers after

    v=u+u^3t-2p, z=uv, A=a-4p, k=d+2bp.

They give the complete `H=z^2+A z+k+v(b-z)` and
`L=lambda(z-v)+e-2p lambda`. The constant shift of the target is
legitimate. The inverse `u=z/v`, `t=(v+2p-u)/u^3` is rational, and
the volume follows directly from

    omega=u^-3 du wedge dv=(v^2/z^3) dz wedge dv.

In this residual case lambda and b are nonzero. For `h=H,f=F`, the
displayed `z=Q/E`, `v=W/(lambda E)` invert both equations literally.
The denominator E has nonzero f coefficient, Q has nonzero f
coefficient `-b`, and W has leading f coefficient one in degree two.
The dominance derivative given in the primary is nonzero as a rational
function. Hence the field is actually `C(f,h)`, with f transcendental,
not merely a selected component of a potentially disconnected model.
The transported volume is exactly

    eta=-W^2 dh/(lambda^2 Q^3).

This sign and every factor come from the actual source volume. The
inverse maps make no claim to preserve original polynomiality.

The final primary includes the corrected term
`-A(h-k)/b` in B0. It is necessary for the complete identity

    W=w^2+lambda B1 w+lambda^2 B0, Q=bw,
    w=h^2-(lambda/b)h+lambda k/b-f.

The earlier scout omission was repaired before the source freeze. No
correction remains. The polynomial part of `W^2/w^3`, including
its f-dependent coefficient, has a rational primitive in the h-line
over `C(f)`.

For the remaining simple, double and triple w-poles, expansion in
the uniformizer `w=P(h)-f` gives the three terms of the primary's
residue formula (16), with derivative `D_P=(1/P')d/dh` and its
factor `1/2` for the triple pole. At a generic moving root, h is
transcendental over the original parameter field. Thus vanishing
of its residue requires a rational identity in h; it is not a
condition checked at a finite sample of levels. The critical points
of the fixed degree-two P give finitely many exceptional levels and
cannot change that identity.

Independently taking numerator and denominator degrees, rather than
the producer's symbolic limit, gives degree difference one and
leading coefficient `3/b^2`. The first term contributes this slope;
the derivative of the double-pole term is bounded and the final
term tends to zero. No parameter specialization can kill it while
b is nonzero. This proves the claimed rational nonexactness for
every polynomial candidate that enters this residual family.

## 5. The other placement: actual quadratic trace and all degeneracies

Direct original-carrier expansion with `A=a+4p`, `k=d+2bp`
gives

    H=u v(v+A)+bv+k,
    L=lambda u v+mu v+e,

where e absorbs the actual constant `2p mu`. Set `h0=bv+k`.
The inverse `u=(h-h0)/(v(v+A))` shows that the actual source
field is `C(v,h)`, and direct differentiation gives

    omega=v^2(v+A)^2 dh wedge dv/(h-h0)^3.

The polynomial equation in h over `C(f,v)` has discriminant linear
in f with coefficient four. It is therefore nonsquare in
`C(v)(f)`, and the extension is genuinely separable quadratic.
This pays the trace even if lambda is zero. The denominator
`v(v+A)` is a nonzero rational function also when A is zero.
No quotient map is assumed birational unless the displayed inverse
establishes it.

I checked the full trace independently by replacing the second root
of `y^2+B y+w` with `-B-y`, adding the two rational expressions,
and reducing the numerator modulo the quadratic. The result is

    Tr(1/[y^3(2y+B)])=1/w^2-B^2/w^3.

It agrees with the producer's separate companion-matrix calculation.
The derivation holding f fixed extends uniquely across a separable
algebraic extension, so trace commutes with differentiation. A
primitive in the actual source field would therefore give a rational
primitive of the full trace in `C(f,v)`. Normalizing the trace
would multiply it by `1/2`, without changing the conclusion.

When b is nonzero, the moving polynomial `P=(bv+k)^2+mu v+e`
has degree two. Residue vanishing is precisely

    (1/2) D_P^2 L1-D_P L2=0,
    L1=v^2(v+A)^2 B^2/P', L2=v^2(v+A)^2/P'.

It is a rational identity in v, by the same generic moving-root
argument. Its first integration gives

    L1'-2v^2(v+A)^2=c_* P'.

Here the expression being differentiated already depends only on v
and the original parameters; its derivative-zero constant belongs to
that constant field. There is no unaccounted f-dependent function.
The left side has degree four and leading coefficient eight, whereas
the right side has degree at most one. I recovered these degree and
leading-coefficient statements by rational numerator/denominator
comparison independently of the producer's limit. They include every
mu, lambda zero or nonzero, and A zero or nonzero. Generic moving
roots avoid the fixed poles at v=0 and v=-A; their possible coincidence
with each other when A=0 creates no missing parameter branch.

## 6. The b=0 case, all-p rational positives, and conclusion

When b is zero the polynomial entry makes mu zero. Constant L is
handled at the end. Otherwise lambda is nonzero, and

    z=(f-h^2)/lambda, v=(h-k)/z-A

is an actual rational inverse. The relative differential calculated
from the same volume is

    eta=(h-k-Az)^2 dh/(lambda z^6).

Passing to `f=rho^2` is an algebraic extension of the generic
constant field, not a special numerical fibre. Exactness persists
under this extension. Removing the nonzero common factor lambda
cubed leaves

    [lambda(h-k)-A(rho^2-h^2)]^2 dh/(rho^2-h^2)^6.

I independently computed the residue by setting `h=rho+epsilon`.
The numerator is a polynomial of degree four in epsilon, and the
coefficient of epsilon to the fifth power in its product with
`(2rho+epsilon)^-6` gives

    -[80A^2 rho^4+140Ak lambda rho^2
                  +63k^2 lambda^2-7lambda^2 rho^2]
                         /(512rho^11).

This calculation uses a finite binomial convolution, independently
of the producer's residue routine. Vanishing as an identity first
forces A=0 from degree four and then lambda=0 from degree two.
There is no division by k and no exceptional k branch. This
contradicts the nonconstant-L branch under consideration.

For constant L, `J(H^2+L,G)=2H J(H,G)` cannot be a nonzero
constant in the original polynomial ring. This final factor argument
does not apply to a rational G, which is precisely the sharp boundary.

For both positives, write `H=A0(u)v^2+d`, `G_H=B0(u)/v`.
The actual source bracket reduces to the scalar identity

    A0' B0+2 A0 B0'=-u^-3.

The displayed pairs

    A0=u(u-1), B0=-(8u^2+4u+3)/(15u^3),
    A0=u,      B0=1/(5u^3)

satisfy it exactly. This independently gives both Jacobians equal
to one without copying the producer's original-variable expansions.
The full numerator-row checks show the associated H are actual
global sections for every p. Dividing `G_H` by `2H` then gives
the advertised rational mates of `H^2+e`. The poles are retained;
neither example asserts a polynomial or globally regular mate.

## 7. Exact replay, independent path, and frozen pins

I read the complete standalone
[producer](../../04-computation/planar_jc48_sep08_seven_one_two_locations.py).
It imports SymPy and hashing, but no inherited mathematical program.
Its gates use explicit exceptions rather than removable assertions.
Its symbolic universe consists of both complete post-entry carrier
families with arbitrary original p, their exact field and volume
maps, the entire residue expansion, the full quadratic trace, both
generic nonzero leading coefficients, the exceptional residue, and
the same-class rational positives. These are literal identity
controls supporting the analytic proof. They are not a finite
enumeration of all polynomial mates or a substitute for the entry
and generic-fibre arguments.

Fresh independent normal and optimized runs both exited zero and
reproduced **46 gates**. Each output is byte-for-byte equal to the
**322-byte** frozen producer output. The inspected pins are:

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Primary, before status promotion | 15143 | `2c0657af853fee8fe3f90372ffac93d81cca9e79bc0625a669fef8f0abd7b05c` |
| Source | 5925 | `62aac3b8ce69b6781b0434c32bf40d92a21fd36d30249866285db25d746c8dd1` |
| Frozen output; independent normal; independent optimized | 322 each | `ee0d52487c625eb31ea960a20b5fed7f515f9696da162e72464b65c9a68b0431` |

Semantic gate SHA256:
`d6888a5c8c4ad87d99944d0224fe9b9098458fb83c8dbd23af7bc1d54a9b6d66`.

Reproduction commands, from the worktree root, are

```sh
python3 -B 04-computation/planar_jc48_sep08_seven_one_two_locations.py
python3 -B -O 04-computation/planar_jc48_sep08_seven_one_two_locations.py
```

The independent outputs were retained temporarily as
`/tmp/planar_jc48_sep08_seven_one_two_locations_audit_normal.out`
and the matching `_optimized.out`. An additional independent
temporary calculation passed **15 alternative controls**: direct
centered inverse coefficient cancellation, source and chart volume
derivatives, paired-root trace reduction, rational degrees and leading
coefficients without limits, the exceptional binomial residue without
a residue engine, and the two scalar ODEs for the rational positives.
The complete formulas for these alternative paths are recorded above;
the temporary script is not a new required repository dependency.

**Final acceptance:** complete analytic scope, inherited entry,
all-position section completeness, original polynomial gate, actual
field maps and volumes, separable trace, generic residue identities,
all displayed parameter branches, source and both replay modes PASS.
The pre-freeze B0 correction is present. No correction or further
producer edit is required. Parent owns promotion and checkpointing.
