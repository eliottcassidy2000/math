# Independent audit of the complete binary (3,3,2) boundary stratum

**Status: PROVED / INDEPENDENT AUDIT PASS.**
The analytic proof, exact source, and independent normal/optimized replays
of [the primary](planar_jc48_sep08_three_three_two.md) are accepted. The
one requested local-text repair was applied and reread before acceptance.
The producer source and output were not edited by this referee. Root owns
the primary's status promotion and integration.

## 1. Accepted statement and its boundary

On the specified DG surface

    W=(P1_x x P1_z)\{z=x^2}, t=1/(z-x^2),
    x=1/r, t=-r^2-r^4 b_D, omega=r^2 dr wedge db_D,

take actual global sections `H in L2`, `L in L1`, and `F=H^2+L`.
If the leading binary octic of `H` has exactly three distinct roots
with multiplicities `3,3,2`, a rational mate forces the double root
to be the original infinity point and `L` to be constant. No such
`F` admits a polynomial mate of any degree. The constant-`L` rational
controls in Section 6 of the primary are genuine global sections for
every finite position `p`; they prevent strengthening the conclusion
to a universal rational-mate exclusion.

This is a complete boundary stratum on this fixed surface. Neither
an entry theorem for arbitrary Keller maps nor a solution of JC(2)
is asserted. The distinction between the original infinity and a
finite root is essential: its canonical multiplier is `r^2`.

I checked the entire primary, the entire companion source, and the
proved local suppliers used for the active and inactive branches.
The inheritance is correctly routed through
[leading exactness](planar_jc48_sep08_leading_exactness.md),
[shared-root first jets](planar_jc48_sep08_shared_roots.md),
[the M-unit Morse classification](planar_jc48_sep08_quartic_common_root.md),
[the global quadratic filtration](planar_jc48_sep08_dg_quadratic.md),
and the noncomposition/geometric-integrality route recorded in
[the all-finite (4,3,1) proof](planar_jc48_sep08_four_three_one.md).
Local pole capacity is a necessary gate, not a trial-space assertion;
the later Riemann--Roch and original-chart checks supply the missing
global information.

## 2. Locations and complete global coefficient rows

When the double root is finite, writing `sqrt(N)=u*unit` on each
local normalization sheet gives a nonzero simple-pole coefficient
for `du/sqrt(N)`. The necessary leading-exactness theorem applies
to the actual algebraic field. It excludes rational mates whether
both triples are finite or one triple lies at infinity. No change
of target coordinates is used to move a finite obstruction to the
original infinity.

With the double root at original infinity, both triples are finite.
Their nonzero separation and the leading scalar can be normalized
using the stated source/surface scalings. The remaining arbitrary
original finite location is retained by `u=x-p`. Choosing which
triple is called zero is done before this normalization; it does
not assert that a translation in `u` preserves the compactification.

For `N=u^3(u-1)^3`, I independently checked the complete rows

    P=sum A_i u^i,
    Q=(A4-1)u^2+(A3-2p A4+4p+3)u+C0,
    M=sum B_i u^i,
    R=B4 u^2+(B3-2p B4)u+E0,       0<=i<=4.

The numerator rows in original `x=u+p` coordinates are

    Q, P-2Qx^2, N-Px^2+Qx^4.

Their required degrees at most four force `deg Q<=2`, since the
terms of degrees eight and seven in the last row cannot be canceled
by its other numerator rows. Then `deg P<=4`. The coefficients of
`u^6,u^5` uniquely give the displayed quadratic and linear
coefficients of `Q`. Conversely these values make every numerator
row satisfy the full degree restriction. For `L`, the two rows
`R,M-Rx^2` similarly force and suffice for the displayed `R`.
Thus there is no unannounced coefficient restriction or missing
shifted parameter in the source universe.

Direct substitution in the original second chart gives polynomial
expressions for both sections and slopes

    partial_(b_D) H|D=2-A4,   partial_(b_D) L|D=-B4.

The source verifies all these assertions symbolically in `p`, including
an independent solve of the high-degree conditions rather than only
checking the proposed rows by substitution.

## 3. Local forms, active roots, and the repaired Morse case

At original infinity,

    r^8 N(1/r-p)=r^2(1-pr)^3(1-(p+1)r)^3.

The full canonical multiplier is retained. A normal unit is regular
by the shared-root supplier, and an M-unit is regular in multiplicity
two because the balanced condition `m=3j` cannot hold. When both
relevant units vanish, the full tangent equation is

    (a r^2+b r s+h s^2)^2+s^3(lr+es)-zeta s^4,
    a!=0.

With `s=rZ`, its constant coefficient is nonzero. If the generic
quartic had a repeated nonzero root, that root would also annihilate
`Z P0'-4P0`. This last polynomial has constant `-4a^2` and is
independent of the generic fibre coefficient. It has only finitely
many possible roots and corresponding exceptional fibre values.
The generic leading coefficient is nonzero as well, so the four
nonzero simple tangent roots exhaust the local Weierstrass degree.
The unweighted relative-form order is minus one; multiplication by
`r^2` makes all four orders one. In particular infinity contributes
no primitive pole capacity.

At a finite triple, the inherited alternatives are exhaustive: an
M-unit gives regular forms or forbidden nonzero logarithms; a shared
normal unit is regular; otherwise a rational mate requires

    P=P'=M=M'=0.

For this active case, `u=tau^2,s=tau^3 Z` gives four generic nonzero
simple determinations. The involution in the parametrization pairs
them into two normalized branches. Each has differential order
minus two and primitive pole order one, so the total capacity is
two. The nonzero leading constant and generic fibre coefficient pay
the root count even if some displayed lower coefficients vanish.

There is one pre-acceptance repair in the primary. A double root of
the balanced cubic at the other triple does **not** invariably give
a logarithm. For `m=3,j=1`, the actual Morse contact is one or two:
contact one gives a ramified unit differential, while contact two
gives the nonzero simple pole that excludes exactness. The generic
fibre parameter bounds the contact by two. The corrected primary
explicitly retains both cases and calls the three leading objects
Puiseux determinations, rather than necessarily three distinct
normalized branches. This repair follows the already-proved local
supplier; it does not change a source identity, the total primitive
capacity, or the final genus calculation. The repaired paragraph
has been reread and accepted.

## 4. Unique-active entry and constant D

If neither finite triple is active, a primitive has no poles on any
compact normalized generic component. It would therefore be constant,
contrary to its nonzero relative differential. This argument is
componentwise and does not assume generic irreducibility.

If both are active, the degree-four lower coefficient is necessarily

    M=lambda u^2(u-1)^2.

The complete inverse coefficient and field trace require `M du/N`
to be rational-exact. Its residues at zero and one are `-lambda`
and `lambda`. Thus `M=0`; the complete rows then force `L` constant.
For a nonconstant `L`, exactly one triple is active and can be
named zero. Its capacity two bounds the primitive degree on every
generic component.

If `F|D` were nonconstant, a generic transverse point of the original
divisor `D` would contribute a differential zero of order two from
`omega=r^2 dr wedge db_D`. The primitive would have local degree
three, larger than its total pole degree at most two. This proves
constant `F|D` before the later geometric-integrality step.

The square of a nonconstant affine function of `b_D` cannot be
canceled by an affine function. The slope identities therefore give
`A4=2,B4=0`. Intersecting with the active jets gives exactly

    P=2u^4+a u^3+b u^2, Q=u^2+(a+3)u+c,
    M=u^2(lambda u+mu), R=lambda u+e.

The two finite residues of `M du/N` are `-mu,mu`, so `mu=0`;
nonconstant `L` then means `lambda!=0`. The surviving differential
indeed has rational primitive `-lambda/[2(u-1)^2]`, so this residue
gate is not being overused as a complete obstruction.

At one, `M=lambda` is a unit. If `P(1)` were nonzero, the finite
triple would contribute a differential zero of order four, forcing
primitive local degree five and contradicting capacity two. Hence
`b=-a-2`. If `P(1)=0`, simple leading roots and the admissible
contact-one repeated root all give unit forms, as checked above.

## 5. Geometric integrality, genus, and every primitive pole

The noncomposition argument is valid in the full polynomial ring,
not only over a selected branch. Since the degree in `t` is four,
an outer polynomial in a nontrivial decomposition has degree two
or four. In outer degree four, its inner polynomial is linear in
`t`, and the leading coefficient would make `N^2` a fourth power;
the odd multiplicities of `N` make that impossible. In outer degree
two, complete the outer square and choose the inner leading sign
to agree with `N`. Successive coefficients of `t^4,t^3,t^2`
force the inner polynomial to be exactly `H`; the remaining term
would make `L` constant. This contradicts the residual hypothesis.
The characteristic-zero closed-polynomial criterion cited by the
primary therefore supplies a geometrically integral generic curve.
It is legitimate to work over the algebraic closure of `C(zeta)`.

On its compact normalization the differential has four order-one
zeros at original infinity and two order-minus-two poles at the
active triple. The other triple is a unit in every admissible
Morse case. The constant divisor `D` misses the generic fibre.
All other affine source points are smooth and have a unit volume
and generic relative form. Thus there are no uncounted zeros or
poles, and the canonical degree is `4-4=0`: genus one.

If `E` is the sum of the two active points, every primitive lies in
`L(E)`, and on genus one `l(E)=2`. I checked the proposed second
basis function at all points:

    J=u(u-1)^2 t+u.

* At each active point, `u~tau^2,t~tau^-3`; its leading simple
  pole is nonzero, so it is nonconstant and has the required pole
  bound.
* At the inactive triple, including the ramified contact-one case,
  `t` has order corresponding to `(u-1)^-2`; the double factor
  cancels it.
* There is no affine denominator. In particular no cancellation at
  a deleted boundary point creates an ordinary affine pole elsewhere.
* Its full numerator rows put it in the actual global `L1`. In the
  original `r,b_D` chart it is polynomial and has value `2p+2` on `D`.
* On each deleted infinity branch, `b_D=1/(rZ+O(r^2))`, direct
  substitution gives `J -> 2p+2-1/Z`. Every tangent root is nonzero,
  so all four limits are finite.

Consequently `1,J` are a complete basis, not merely two functions
with plausible local orders. Coefficients of a primitive may belong
to the algebraic closure of the generic constant field; this does
not weaken the bracket argument below.

## 6. Independent bracket and hostile checks

I independently reconstructed `J_(u,t)(F,J)` and its remainder
modulo the literal degree-four equation `F-zeta`, with

    P=u^2(2u^2+a u-a-2), Q=u^2+(a+3)u+c,
    F=(N t^2+P t+Q)^2+lambda(u^3 t+u).

The `t^3` coefficient is zero; the `t^2,t^1` coefficients are

    2c u^3(u-1)^4,
    2c u^2(u-1)^2(a+2u+2),

and the constant coefficient is

    2ac u^2-2ac u+2c^2u-2c^2+2c u^3+4c u^2-6c u
        -lambda u-2zeta u+2zeta.

A primitive `A(zeta)J+B(zeta)` has `A!=0`. The coefficient of
`t^2` forces `c=0`, after which the full remainder is

    2zeta-(2zeta+lambda)u.

It cannot be scalar on the formal generic fibre. The explicit
exceptional level `zeta=-lambda/2` gives the nonzero constant
`-lambda`; this hostile is preserved in the source and primary.
The proof uses the generic field and is not an inference from one
selected fibre. Allowing algebraic generic constants cannot cancel
the nonzero coefficient of `u`.

Finally I checked the all-position rational control directly:

    A0=u(u-1), Z0=A0 t+1,
    H=A0 Z0^2+c, G_H=(2u-1)/(A0 Z0).

Its coefficient rows satisfy the complete original globality
conditions for symbolic `p`, its leading polynomial is `A0^3`,
and its Jacobian with `G_H` equals one. Hence `H^2+e` has rational
mate `G_H/(2H)`. This is compatible with the polynomial exclusion:
when `L` is constant the polynomial bracket has the nonconstant
factor `2H`, so it cannot equal one.

## 7. Source scope, independent replays, and pins

The companion controls the entire symbolic row space, both finite
residue reductions, the original infinity coordinate and limit,
all coefficients of the decisive remainder, the special-fibre
hostile, and the actual all-p rational family. Its numerical order
and Riemann--Roch checks are arithmetic controls for the preceding
analytic arguments, not independent proofs of local exhaustion or
generic integrality. There is no finite parameter census being
promoted to an unrestricted theorem. Every gate is always active
and raises an exception on failure; Python optimization cannot
remove it.

Independent commands, run from the declared worktree, were

```sh
python3 -B 04-computation/planar_jc48_sep08_three_three_two.py > /tmp/three-three-two-audit-normal.out
python3 -B -O 04-computation/planar_jc48_sep08_three_three_two.py > /tmp/three-three-two-audit-optimized.out
```

Both completed successfully with **56 gates**. Their complete raw
outputs are byte-for-byte identical to the 268-byte frozen output.
The final repaired primary, producer source and frozen output were
read and hashed after these checks.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Repaired primary, before status promotion | 12,947 | `3373db80160074700da54760d4805ad8b2f0115c2fcafbe530be953f257cd1c0` |
| Producer source | 5,623 | `eabf9adaeeeefe68c543fab37963af1c14d75d9e39f3ee34b964963617de3ea0` |
| Frozen output and each independent replay | 268 | `37f4bd8ad522d01ccaf93a20531025e0a7b3e3645b1e91286eb6503bc458c03e` |

The semantic gate digest is
`facd85633ee6f36fb3d7851442e670d799ea602633b65301389f16ed1403dd11`.
No unresolved mathematical or source correction remains. This audit
accepts precisely the stated rational necessary condition and full
polynomial exclusion for the actual binary `(3,3,2)` stratum.
