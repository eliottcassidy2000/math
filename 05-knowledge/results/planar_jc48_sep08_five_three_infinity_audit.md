# Independent audit: both infinity placements of binary (5,3)

**Status: PROVED / INDEPENDENT AUDIT PASS.**
The full analytic argument, frozen source, and independent normal and
optimized replays of [the primary](planar_jc48_sep08_five_three_infinity.md)
are accepted. No mathematical or source correction was needed. This
referee edited neither the producer artifacts nor an inherited bundle;
root owns promotion and integration.

## 1. Precise accepted scope

For actual global sections `H in L2`, `L in L1` on the fixed DG surface

    W=(P1_x x P1_z)\{z=x^2}, t=1/(z-x^2),
    x=1/r, t=-r^2-r^4 b_D, omega=r^2 dr wedge db_D,

assume that the leading binary octic has two distinct roots, with
multiplicities five and three, and one of these is the original
infinity. The primary excludes polynomial mates of `F=H^2+L` of
every degree, for both placements and every original finite position.
It does not exclude all rational mates. Both stated nonconstant-`L`
rational families are actual global controls; their pole-repair
arguments exclude polynomial mates without denying the rational ones.

The [all-finite proof](planar_jc48_sep08_five_three.md) and its
[independent audit](planar_jc48_sep08_five_three_audit.md) are separate
proved inputs. Their combination with this accepted theorem covers
the entire binary `(5,3)` polynomial-mate stratum. It is not an entry
theorem for arbitrary Keller maps, a theorem for all quartic
polynomials, or a resolution of JC(2).

I read the entire primary and source, checked the complete weighted
local classification, and independently reconstructed the decisive
field maps, critical-point alternatives, Jacobians, and principal
parts. The proof keeps the actual generic fibre constant throughout.
Its inherited inputs are the full global filtration, the shared-root
first-jet obstruction, the M-unit Newton/Morse classification, the
complete inverse-coefficient exactness identities, and the precise
noncomposition/geometric-integrality criterion. Their use is consistent
with the scope stated in the primary.

## 2. Complete global rows and retained original position

After normalizing the leading scalar, let `u=x-p`, with arbitrary
original finite location `p`. This is a polynomial source coordinate;
the original compactification is still used to test globality. For
`N=u^d`, `d=3,5`, the complete rows are

    P=sum A_i u^i,
    Q=A4 u^2+(A3-2p A4-e_d)u+C0,
    M=sum B_i u^i,
    R=B4 u^2+(B3-2p B4)u+E0,
    0<=i<=4, e_3=0, e_5=1.

I checked these against the original numerator row bounds. In
original `x`, the relations are `Q2=P4`, `Q1=P3-N5`, and
`R2=M4`, `R1=M3`. Expanding `x=u+p` gives exactly the displayed
coefficients; there is no omitted `p` term. Conversely these complete
relations put both sections in their required global spaces.
The source verifies the original-coordinate coefficients and direct
polynomiality after substituting `u=1/r-p,t=-r^2-r^4 b_D`.

Thus later disappearance of `p` in some residuals is a consequence
of intersecting the full coefficient space, not an assumed surface
translation. The polynomial scalings used to normalize the conic
exception preserve polynomial-mate existence, and are applied after
the actual global rows have been paid.

## 3. Complete weighted infinity lemma

The load-bearing local form is

    calN=r^m a(r)+s B(r)+s^2 C(r,s), a(0)!=0,
    calM=D(r)+s E(r,s),
    Phi=calN^2+s^3 calM-zeta s^4,
    eta=(unit) r^2 s^2 dr/Phi_s,

with `m=3` or `5`. This includes the actual canonical multiplier
`r^2`. I checked the regularity claim for every order pair
`j=ord B`, `n=ord D`, including vanishing functions and higher
contact. A normal unit is covered by the proved analytic-centre
lemma. An M-unit in multiplicity five is unbalanced and regular;
in multiplicity three, the balanced `j=1` case has at worst the
already-proved unweighted logarithmic alternatives, which become
regular after multiplication by `r^2`.

For `j,n>=1`, the local Weierstrass degree is four. The remaining
cases in the primary exhaust those four determinations as follows.

| Case | Normalized branches or determinations | Weighted differential orders |
| --- | --- | --- |
| `j=1`, either `m` | Two simple `s~r` branches; two determinations centered at order `m-1` | The simple branches have order one; the cancelling expression is `r^((m+1-ell)/2) dr`, `1<=ell<=m-1`, hence regular |
| `m=3,j>=2,n=1` | One `s~r` branch and one ramified branch with `(ord r,ord s)=(3,5)` | One and five |
| `m=3,j>=2,n>=2` | Four determinations `(r,s)=(tau^2,tau^3 Z)`, paired into two branches | Two on each |
| `m=5,j=2,n>=2` | Two `s~r^2` branches and two centered at order three | Zero on the simple branches; cancelling expression `r^((3-ell)/2) dr`, `1<=ell<=3` |
| `m=5,j>=3,n=2` | One `s~r^2` branch and one branch of orders `(3,8)` | Zero and two |
| `m=5,j>=3,n>=3` | Four determinations `(r,s)=(tau^2,tau^5 Z)`, paired into two branches | Zero on each |
| `m=5,n=1,j>=3` | One `s~r` branch and three simple `s~r^3` determinations | One on each |
| `m=5,n=1,j=2` | Same low branch, with the three high determinations governed by the full cubic below | All regular, including its double-root contact |

For the cancelling centres, the generic term `-zeta s` bounds the
contact by the stated centre order. A nonnegative rational exponent
of `r` multiplying `dr` remains regular after normalization; the
ramification contribution to `dr` must be added and is nonnegative.
For the fractional faces, direct evaluation of
`2 ord r+2 ord s+ord(dr)-ord(Phi_s)` gives the listed orders.

The special face for `m=5,j=2,n=1` is

    (a0+B2 Z)^2+D1 Z^3.

Its constant is nonzero. It has no triple root; a double root is
necessarily `Z=-3a0/B2`, `D1=4B2^3/(27a0)`, and its second
derivative is nonzero. I checked those identities independently.
The generic fibre enters two weights above this initial face, so
the derivative of the Morse critical value with respect to `zeta`
has exact order two and nonzero leading coefficient. Contact one
gives a ramified branch with `eta=(unit) r dr/v` of order two;
contact two gives two smooth branches of order zero. Neither case
may be deleted or mislabeled as an unweighted finite obstruction.

The other nonzero faces are generically simple: their repeated-root
eliminants have a nonzero constant and are independent of the generic
fibre parameter. The monomial fractional faces have nonzero simple
roots directly. Higher terms either have strictly higher weight or
enter the explicitly bounded centre/Morse calculation. No branch
count is inferred from a numerical Puiseux sample.

All relative forms at original infinity are therefore regular.
At points of retained `D` they are regular as well, and ordinary
smooth points on the generic affine fibre have unit relative form.
This pays the full compactification needed for the finite pole
budgets used below.

## 4. Finite triple, quintuple infinity: entry and genus obstruction

For `N=u^3`, a rational mate with nonconstant `L` must have an
active finite triple. Otherwise the preceding lemma and the finite
unit alternatives leave no primitive poles on any compact generic
component, contradicting its nonzero differential. The first jets
therefore give `P0=P1=M0=M1=0`. Four nonzero simple determinations
at `u=tau^2,s=tau^3 Z` give two actual double differential poles,
with total primitive pole capacity two.

A nonconstant restriction to original `D` would give a generic
differential zero of order two and hence primitive local degree
three. The pole-capacity contradiction applies on its own component,
without using irreducibility. Thus `F|D` is constant. The original
slopes are `-A4,-B4`, so constant `F|D` forces `A4=B4=0`.
The traced inverse coefficient `T2` then forces `B2=0` by residue.
The complete nonconstant-`L` residual is

    H=u^3t^2+u^2(a u+b)t+a u+c,
    Z=u^3t+u, L=lambda Z, lambda!=0.

For `a!=0`, direct infinity coefficients are

    calN(r,0)=r^5(1-pr)^3,
    [s]calN=-a r+O(r^2), calM(r,0)=-lambda r+O(r^2).

The two low branches `s~r` have differential zeros of order one.
The cancelling centre has order four and contact exactly one,
because the leading M term is nonzero. Its two determinations form
one branch with `r=tau^2`, differential order six. Together with
the two finite double poles, and no generic point on the constant
divisor `D`, the canonical degree is `1+1+6-2-2=4`.

The noncomposition proof is valid: outer degree four would make
the leading `N=u^3` a square up to scalar; outer degree two, after
completing the square, forces the entire square prefix and then
`L` constant. The precise characteristic-zero closed-polynomial
criterion cited in the primary supplies geometric generic
integrality. The normalized generic curve consequently has genus
three.

The second ingredient is an actual degree-three map, not a trial
function count. Put `U=1/u,z=Z`; the inverse is
`u=1/U,t=(z-1/U)U^3`, so this is a birational field change and

    h=H=z^2 U^3-2z U^2+(1+bz)U+a z+c-b,
    omega=-U dU wedge dz,
    z=(zeta-h^2)/lambda.

I reconstructed this identity and its leading cubic coefficient
independently. Over `C(z)`, the polynomial in `U` has degree three,
so `C(U,z)/C(h,z)` has degree exactly three. The geometrically
integral generic curve retains this degree on extending the generic
constant field.

A rational primitive is a nonconstant map of degree at most two.
Degree one would make the curve rational. At degree two, the fields
of the primitive and of `h` generate the full curve field: the
remaining extension degree divides both two and three. The map
to `P1 x P1` is therefore birational onto a curve of bidegree
`(2,3)`, whose arithmetic genus is two. Normalization cannot raise
it to three. This proves the claimed rational exclusion in this
subcase. No unsupported cubic trace or guessed primitive basis is
used.

## 5. The finite-triple a-zero branch and its conic control

For `a=0`, the genuine chart `v=ut`, `u!=0`, gives

    H=u v(v+b)+c,
    F=H^2+lambda(u^2v+u).

I independently differentiated this polynomial and recovered the
critical formulas

    u=-(2v+b)/[v(3v+b)], H=lambda/[2v(3v+b)],
    C(v)=-4v^3+6(c-b)v^2+2b(c-b)v-lambda.

Every zero of `C` outside `0,-b/3,-b/2` gives a finite source
point with `u!=0` and both coordinate derivatives zero. The chart
Jacobian is nonzero there, so it is an actual source critical point.
Zero is not a root because `lambda!=0`. If `b=0`, all three
forbidden addresses coincide with zero, so a usable root exists.
For `b!=0`, concentration of all three roots at the other two
addresses has exactly the four patterns listed in the source.
Direct coefficient comparison leaves only

    c=b/3, lambda=4b^3/27,
    C(v)=-4(v+b/3)^3.

This is a complete cubic root-concentration argument, not a bounded
parameter probe. The explicit scaling `u=3U/b,t=b^2 T/9` and
constant output factor reduce this exception to the stated literal
without changing existence of a polynomial mate.

For that literal, direct original-coordinate calculation gives

    H=u^3t^2+3u^2t+1,
    A=1+u(ut+1)(ut+4), B=1+u^2t(ut+1),
    F=H^2+4(u^3t+u)=AB,
    G0=1/[2A u(ut+1)], J(F,G0)=1.

The displayed inverse from `A,B` to `u,v` is rational, so the
function field is `C(F)(A)` and its fibre derivation has constants
exactly `C(F)`. Every rational mate is `G0+K(F)`. This constants
claim is paid by a birational field map, not by a generic-fibre
picture or an assumption about connectedness.

On the two actual disjoint irreducible divisors `E={u=0}` and
`Gamma={1+ut=0}` in `F=1`, I independently checked

    [(F-1)G0]|E=2, [(F-1)G0]|Gamma=1.

The fibre has simple vanishing at the generic point of each divisor.
A single rational function of `F` cannot cancel both distinct simple
principal parts, and any higher pole of that function survives.
This establishes polynomial exclusion while preserving the exact
nonconstant-`L` rational mate.

## 6. Finite quintuple, triple infinity: residues and source criticality

For `N=u^5`, regularity at original infinity and the finite unit
alternatives again force the active finite jets. The residue of
`M/u^5` forces `B4=0`. The resulting full coefficients are

    H=u^5t^2+u^2(a u^2+b u+k)t
        +a u^2+(b-2pa-1)u+c,
    L=n u^2t+m(u^3t+u).

In the actual chart `w=u^2t`, the constant and linear coefficients
in `u` are the displayed `f(w),g(w)`. At a simple root of
`f(w)=zeta`, the relative-form residue is proportional to

    g f''-g' f'.

Its vanishing on the original generic fibre forces the polynomial
identity `(g/f')'=0`, hence `g=C f'`. It is stronger than a
coincidence on one special fibre. If `k!=0`, the nonzero cubic
coefficient `2k` of `g` is incompatible with the linear `f'`.
Thus `k=0`. If `n!=0`, comparing the quadratic and linear terms
forces `c=m=0`; the original polynomial `F` is then divisible by
`u^2`, and the entire source line `u=0` is critical. This already
excludes polynomial mates, independently of any additional rational
residue conditions.

For `n=0,m!=0`, the higher inverse identities are used at the
correct level: they belong to the full original square-prefix
inverse, and normalized trace makes their derivatives rational.
Their residues are

    Res T6=m(2an+bm-2m)/16,
    Res T10 |(n=0,b=2)=c m^3/32.

They force `b=2,c=0`, retaining every `p`. With
`v=u(1+u^2t)` the actual residual becomes

    H=v^2/u+a u(v-2p), L=m v.

For `a!=0`, its fold is `a(v-2p)u^2=v^2`. On this fold,
`H=2v^2/u`, `F_u=0`, and `F_v=12a v^2-16ap v+m`.
I checked the displayed derivatives and their fold reduction.
The quadratic has a root outside `0,2p`: zero is forbidden by
`m!=0`, and concentration at `2p` would force first `p=0` and
then `m=0`. At an allowed root, the square-root equation for
`u` has nonzero solutions. The chart Jacobian from `(u,t)` is
`u^3`, so these are actual source critical points. This argument
covers all nonzero `a,m` and every complex `p`.

## 7. The last rational family and the second principal coefficient

In the remaining case `a=0`, direct original-coordinate calculation
gives

    q=1+u^2t, v=u q, H=u q^2,
    F=H^2+m v, G0=1/(6v^3), J(F,G0)=1.

Its complete global rows hold for all original finite positions.
The field map `(u,t)->(h,v)` is birational: `u=v^2/h`, followed
by the rational formula for `t`. Since `F=h^2+m v`, the field
is `C(F)(h)`. Hence every rational mate is `G0+K(F)`.

The two actual source divisors in `F=0` are `E={u=0}` and
`Gamma={1+u^2t=0}`. They are disjoint and irreducible. At the
generic point of `E`, `v=u+u^3t` and

    H^2=v^2+O(v^4), F=m v+v^2+O(v^4).

At the generic point of `Gamma`, `u` is a unit and
`H^2=v^4/u^2`, so `F=m v+O(v^4)`. I checked the first expansion
also by inverting the original polynomial with its generic source
coordinate `t` retained, rather than substituting a sampled point.
These give respectively

    G0=m^3/(6F^3)+m/(2F^2)+O(F^-1),
    G0=m^3/(6F^3)          +O(F^-1).

The nonzero difference in the second principal coefficient survives
any common rational correction `K(F)`. A higher pole of `K` could
not disappear either. Thus this exact rational family has no
polynomial mate. The proof uses actual source divisors, not just
places of an auxiliary rational model.

When `L` is constant in either placement, the bracket of a
polynomial with `F=H^2+constant` has the nonconstant factor `2H`.
This completes the polynomial exclusion and is consistent with
both rational sharpness families.

## 8. Independent exact controls, replays, and final pins

The source was read in full. Its exact universe is symbolic:
complete shifted section rows; every declared local face regime;
all four possible cubic root concentrations; the traced higher
residues; actual critical-point eliminations; birational field
inverses; two exact rational Jacobians; and the distinct same-fibre
principal parts. The order arithmetic controls do not replace the
analytic normalization argument, and no finite coefficient scan
is treated as a universal theorem.

Separate direct calculations, without importing the producer, recovered
the cubic field map, the critical compatibility and both critical
derivatives, all four root-concentration patterns, the conic factor
and Jacobian, both simple principal coefficients, the second exact
rational Jacobian, and the original-source inversion yielding the
third and second principal coefficients along `E`. A local check
script initially requested a negative SymPy series truncation order;
using the supported order zero produced the complete Laurent
principal part and passed. This was a referee command issue, not a
producer identity or theorem correction.

Independent full replays were run from the declared worktree:

```sh
python3 -B 04-computation/planar_jc48_sep08_five_three_infinity.py > /tmp/five-three-infinity-audit-normal.out
python3 -B -O 04-computation/planar_jc48_sep08_five_three_infinity.py > /tmp/five-three-infinity-audit-optimized.out
```

Both completed successfully with **107 always-active gates** and
reproduce the frozen output byte-for-byte. No assertion is removed
by optimization. The entire final primary, source and output were
reread and hashed at acceptance.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Primary before status promotion | 19,214 | `ee9f82c3ec006994e4e3454999c2e2596c5b4bf0aa77af848633b9004c7e4c3e` |
| Frozen producer source | 8,514 | `a693f7ad98a082d6897fa247316904f89056a146a771fe04862eccd3067c15a1` |
| Frozen output and each independent replay | 349 | `05ee85e0d98edf3f36d14b8cfe7116aff7824b814c5e7b0abfe0262d3088cceb` |

The semantic record digest is
`d2c7bd4f6fead41fa623987705af5b830f098bf14243d9a1f53d9b44d70fcf73`.
All stated analytic and source obligations pass. No remaining
repair is requested, and the accepted rational/polynomial scope
is exactly the one stated in Section 1.
