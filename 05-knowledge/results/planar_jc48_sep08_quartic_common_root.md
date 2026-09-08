# Every rational mate of a global quartic square prefix requires a common boundary root

**Status: PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT controls PASS.**
The statement concerns the fixed DG surface and the full rational source
field. It is not a proof of JC(2), nor a construction from the remaining
necessary boundary condition.

## 1. Actual claim, retained data, and dependencies

Use `W=X\S`, `X=P1_x x P1_z`, `S={z=x^2}`, with

    t=1/(z-x^2),  r=1/x,  t=-r^2-r^4b,
    D={r=0} in W,  omega=dx wedge dt=r^2 dr wedge db.

The [full global filtration](planar_jc48_sep08_dg_quadratic.md) identifies
global functions `H in L_2`, `L in L_1` with sections

    H=N/s^2,   L=M/s,   s=z-x^2,
    N in H^0(X,O(4,2)),   M in H^0(X,O(2,1)).

Assume `deg_t H=2`, so the binary octic `N|S` is nonzero. Its zero
multiplicities, including infinity, sum to eight. Put `F=H^2+L`.

**Theorem.** If `N|S` and `M|S` have no common zero, then
there is no `G in C(x,t)` with `J_(x,t)(F,G)=kappa!=0`.
In particular, any hypothetical global quartic Keller pair with this
actual square-prefix entry must have a common boundary root.

The source is this fixed `F`, its compact generic fibre, and its actual
relative form. The map restricts `omega/dF` to each normalized fibre
component. It preserves exactness of a proposed rational mate and all
labelled boundary valuations. A residue test alone loses the higher
principal parts. The needed sidecar is the degree of the primitive on
the same compact component, not just the genus or a local primitive.

Closest mechanisms are the safe local proof and the global higher-pole
hostiles in [the quartic differential theorem](planar_jc48_sep08_quartic_differential.md),
and the componentwise pole-degree argument in
[the pole-degree note](planar_jc48_sep08_pole_degree.md).
The argument below also spells out the degree consumer directly; that
companion has now passed its independent audit. The global
square-prefix supplier is [the quartic boundary theorem](planar_jc48_sep08_quartic_boundary.md);
it is a consumer route, not an assumption that arbitrary source
coordinates already extend to this surface.

Live concepts: a binary octic's total multiplicity; the three-term
Newton polygon; critical values of a double-root family; primitive
pole degree; and the actual zero divisor on `D`. The corrected near
miss is to claim that `M`-unit makes the relative form regular or
logarithmic. The multiplicity-six second-kind example refutes both
shortcuts. The present proof retains its pole budget instead.

## 2. Local classification when M is a unit

Near a boundary point use a local coordinate `w` on `S` and equation
`s=0` for `S`. The generic fibre is

    E=N(s,w)^2+s^3 M(s,w)-c s^4=0,
    N=A(w)+s B(w)+s^2 C(s,w),                            (1)

where `M(0,0)!=0`. Set `m=ord_w A>=1` and `j=ord_w B`, allowing
`j=infinity` when `B=0`. Write `a,b` for the nonzero leading
coefficients of `A,B` when relevant. Up to a holomorphic multiplier,
the relative form is

    eta=s^2 dw/E_s.                                     (2)

Away from `Dbar` that multiplier is a unit. At `S intersect Dbar`
it has an extra factor `w^2`, so the estimates below only improve.
The full infinity-chart computation and transitions are given in
Section 2 of the audited quartic differential note.

The following cases exhaust the local branches. One can read them
from the lower Newton polygon of (1), or obtain them by the indicated
rescalings. At `w=0`, (1) has Weierstrass degree two in `s` when
`j=0` and degree three when `j>0`, including `j=infinity`.

### 2.1 The case m<3j

Here every branch has `ord_w s=2m/3`, with leading balance
`a^2 w^(2m)+M(0,0)s^3=0`. The three nonzero roots of this leading
cubic are simple, before allowing the necessary Puiseux ramification.
The order of `N` is `m`; the order of `N_s` exceeds `m/3` because
`N_s=B+2sC+s^2C_s`. Consequently

    E_s=3M(0,0)s^2(1+terms of positive order).

Formula (2) is a unit times `dw`, hence regular on every normalized
branch. This includes `B=0` and every sufficiently large finite `j`.

### 2.2 The case m>3j

First consider `j>0`. There is one simple branch with
`s~-(b^2/M(0,0))w^(2j)`. On it the terms `sB` and `s^3M`
give the leading equation and `E_s` has the same order `4j` as
`s^2`; its differential is regular.

The other two Puiseux determinations, also the entire boundary
part when `j=0`, have the cancellation balance

    s~-(a/b)w^n,       n=m-j>2j.

On either branch, `N^2=-s^3(M-cs)` forces `ord_w N=3n/2`.
Moreover `N_s~b w^j`, and the term `2N N_s` strictly dominates
the other terms of `E_s`, since `j<n/2`. Therefore

    eta=(unit) w^((m-3j)/2) dw.                         (3)

These determinations may belong to one ramified normalized branch.
The form is regular even after that possible quadratic normalization:
the exponent is positive, and ramification adds the nonnegative
order of `dw`. There are no omitted higher-contact branches:
the simple branch plus this pair exhaust the stated Weierstrass
degree. For `j=0` the third putative Newton root has `s`-order zero
and does not pass through the boundary point.

### 2.3 The balanced case m=3j=3k

Here `k>=1` and every boundary branch has `s=w^(2k)Z` with a
nonzero leading value of `Z`. Divide (1) by `w^(6k)`:

    Q(w,Z,c)=Q0(w,Z)-c w^(2k)Z^4,
    Q(0,Z,c)=P(Z)=(a+bZ)^2+M0 Z^3.                     (4)

The division has no negative `w` powers; in particular the higher
`s^2 C` term contributes only positive powers after rescaling.
The cubic has no zero root. It cannot have a triple root: comparing
with `M0(Z-Z0)^3` would give

    b^2=-3M0 Z0,   2ab=3M0 Z0^2,   a^2=-M0 Z0^3.

The identity `(2ab)^2=4a^2b^2` would then require
`9M0^2 Z0^4=12M0^2 Z0^4`, impossible. A simple root of `P`
has `Q_Z` a unit and hence gives a regular differential.

It remains to quantify a double root `Z0`. Because `Q_ZZ` is a
unit there, the implicit function theorem gives its critical centre
`Z=z(w,c)` satisfying `Q_Z=0`, `z(0,c)=Z0`. Let

    R(w,c)=Q(w,z(w,c),c).

The chain rule gives the exact identity

    partial_c R = -w^(2k) z(w,c)^4.                    (5)

Thus coefficients of orders below `2k` are independent of `c`, and
the coefficient of `w^(2k)` has nonzero linear slope `-Z0^4`.
Over `C(c)`, or for `c` outside a finite exceptional set,

    1<=lambda=ord_w R<=2k.                             (6)

This is the required genericity proof; it is not inferred from
one perturbation or an endpoint sign. Taylor expansion around the
critical centre and an analytic square root of the quadratic unit
put the equation into the form `v^2+R=0`. The coordinate change
has unit `Z` derivative, so `Q_Z` equals a unit times `v` on it.
Since `E_s=w^(4k)Q_Z` and `s^2=w^(4k)Z^2`, formula (2) becomes

    eta=(unit) dw/v.                                  (7)

If `lambda=2l`, there are two smooth branches with differential
pole order `l` each. If `lambda=2l+1`, one normalized branch has
`w=tau^2`, `v=(unit)tau^(2l+1)`, and differential pole order
`2l=lambda-1`. The latter includes the regular case `lambda=1`.
In particular higher than simple poles are allowed and retained.

## 3. The complete primitive budget for the binary octic

For a differential of order `-p<0`, a rational primitive can have
only a pole of order `p-1`; a nonzero simple pole instead forbids
an exact primitive outright. Define the possible primitive pole
budget as `max(p-1,0)`, summed over normalized branches. The even
and odd alternatives in (7) both give exactly

    local budget = max(lambda-2,0) <= 2k-2.             (8)

Any additional zero of the canonical numerator can only decrease
this budget. Unbalanced roots and simple roots of the balanced
cubic contribute zero.

Because all boundary multiplicities sum to eight, the only
balanced possibilities are `m=3` and `m=6`. At `m=3`, (8) is zero.
At `m=6`, it is at most two, and there can be at most one such
boundary point. All other multiplicities `1,2,4,5,7,8` contribute
zero by Section 2. Hence the total possible pole degree of a
rational primitive, across every component of the compact generic
fibre, is at most two. If any simple pole occurs, no primitive
exists and the desired conclusion already follows.

## 4. An actual point on D forces degree at least three

The restrictions `H|D` and `L|D` are polynomials in `b` of degrees
at most two and one. If `F|D=H|D^2+L|D` were constant, both
restrictions would be constant: the leading term of a nonconstant
square cannot be cancelled by a polynomial of degree at most one.

In the complete chart at infinity, with `q=1/z` and
`sigma=q-r^2`, let

    Ninf=r^4q^2N(1/r,1/q),  Minf=-r^2qM(1/r,1/q).

On `Dbar`, `H=Ninf(0,q)/q^2` and `L=Minf(0,q)/q`. Constancy
would force `Ninf(0,0)=Minf(0,0)=0`. This is a common zero at
`S intersect Dbar`, contrary to the hypothesis. Therefore `F|D`
is nonconstant.

Choose `c` generic also for this restriction. There is a point
`(r,b)=(0,b0)` on the fibre with `F_b(0,b0)!=0`. The fibre is
smooth there, `r` is its parameter, and

    eta=(unit)r^2 dr.                                  (9)

All generic components are treated separately; no irreducibility
or rational constants-field hypothesis is imposed. The generic
fibre is smooth in `W` after excluding its finitely many critical
values. The form is regular there, and nonzero on each component
meeting `W` away from `D`. Neither `S` nor `Dbar` is a generic
component.

Suppose a rational mate existed. Its restriction to each compact
normalized component is meromorphic for a generic `c`, and
`dG=kappa eta`. This equality prevents any pole of that restricted
function in `W`, regardless of whether `G` was globally regular
on the surface. At the boundary, Section 3 bounds the total pole
degree by two. On the component containing the point (9), the
function is nonconstant and has local ramification index three:
`G-G(b0)` starts with a nonzero multiple of `r^3`. A meromorphic
function's local degree cannot exceed its degree as a map to `P1`,
which equals its total pole degree on that same component. This
gives `3<=degree G<=2`, a contradiction. This proves the claim.

## 5. Hostiles, exact controls, and the remaining entry

The globally admissible multiplicity-six example in the quartic
differential note has residue-free double poles. It survives local
residue tests but fits the sharp budget of two. Its affine critical
points already exclude polynomial mates, so it is not advertised
as a new noncritical Keller class. Its role here is to refute a
false local regularity claim and verify why (8), rather than a
logarithmic-only argument, is necessary.

The exact producer checks all `m=1..8` with representative finite
`j` and the `j=infinity` case, both normalizations in (8), the cubic
no-triple-root identity, and the exact critical-value derivative.
It independently reconstructs the actual multiplicity-six equation
and polar coefficient from its global section. The all-parameter
proof is Sections 2--4, not an extrapolation from this finite bank.

The strongest remaining necessary sidecar is now an actual common
zero of the binary octic and quartic at the same point of the full
boundary. It is not enough to record their separate discriminants
or root multiplicities. A productive next question is whether the
global source bracket constrains the multiplicities at such a
common zero; this theorem makes no sufficient-entry or existence
claim there.

## 6. Frozen exact controls

From the repository root run:

```sh
python3 -B 04-computation/planar_jc48_sep08_quartic_common_root.py
python3 -B -O 04-computation/planar_jc48_sep08_quartic_common_root.py
```

The [standalone source](../../04-computation/planar_jc48_sep08_quartic_common_root.py)
passes **568 always-active gates** in both modes, reproducing all
397 bytes of [the frozen output](planar_jc48_sep08_quartic_common_root.out).
Its explicit finite universe is `m=1..8`, `j=0..9` and `infinity`,
the 22 partitions of eight, and every split order `1<=lambda<=2k`
for `k=1,2`. The `j>=10` tail and arbitrary higher coefficients
are covered by the valuation proof, not an omitted computational
filter. All checks use exact symbolic or rational arithmetic.

| Artifact | SHA-256 |
| --- | --- |
| Source, 7,441 bytes | `df2f4ef7ae1fcca65aa9f99e9982ca26eb27768d2dad095a1402b054fbb1185b` |
| Output and both replays | `0240b334889275c2cf891fff1f995e7d4479cacab013e65000c18af8a73ebf19` |
| Semantic record | `f2bb05e16105c976140835fe0576f92c6ec3ec33a8b0f3e4e2a2bc937b1d3ec7` |

The [independent root audit](planar_jc48_sep08_quartic_common_root_audit.md)
passes the full analytic/source argument and both568-gate replays. A second
independent analytic hostile read agrees. This promotion preserves the
frozen source and output, including their pre-audit production wording.
