# The infinity-octuple quartic class is excluded with its actual volume weight

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.**
The final statement excludes polynomial mates in the original
source plane. The coordinate change below need not preserve
polynomiality of a mate, and every intermediate exclusion is
typed accordingly.

## 1. Statement and the actual surface involution

Let `W=(P1_x x P1_z) minus {z=x^2}`, `t=1/(z-x^2)`. Suppose
`H=N/(z-x^2)^2` and `L=M/(z-x^2)` are global functions in
`L_2,L_1`, and `F=H^2+L`.

**Theorem.** If the complete boundary octic is the
nonzero constant `N|S=alpha`, meaning that its only zero is
the distinguished infinity point with multiplicity eight,
then no `G in C[x,t]` has `J_(x,t)(F,G)=lambda!=0`.

The product inversion `(x,z)->(1/x,1/z)` is an actual
automorphism of `P1 x P1` preserving the graph and hence `W`.
Its new source coordinates are

    u=1/x,       T=-x^2-x^4*t,
    x=1/u,       t=-u^2-u^4*T.

For the new graph coordinate `zeta=1/z`, `s=zeta-u^2`, the
transformed numerators are

    Nnew=u^4*zeta^2*Nold(1/u,1/zeta),
    Mnew=-u^2*zeta*Mold(1/u,1/zeta).

They are actual global sections, and `Nnew|S=alpha*u^8`.
Thus, writing `H,L,F` for the transformed functions, the
complete global entry is

    w=u^2*T,    v=u+u^3*T,    h=u^2+u^4*T=u*v,
    H=alpha*h^2+k0+kt*T+kxt*u*T+k2*w+k3*v+k4*h,
    L=l0+lt*T+lxt*u*T+l2*w+l3*v+l4*h.                 (1)

The full section-space proof is the one in
[the finite-octuple theorem](planar_jc48_sep08_finite_eight.md).
Its unweighted mate conclusion is **not** used: the actual
original volume form here is

    omega_old=u^2*du wedge dT.                        (2)

Rescale a proposed mate to Jacobian one. Its transformed
rational function must satisfy `dF wedge dG=omega_old`.
The transformed mate is not assumed polynomial.

On the new second chart, `u=1/r`, `T=-r^2-r^4*b`, so (2)
is `dr wedge db`: it has no pole there. On the new source
line `u=0`, it has a zero of order two. These are actual
global facts about the volume form; they must not be replaced
by the earlier unweighted boundary-zero count.

The closest mechanisms are the complete local Newton lemmas in
[M-avoidance](planar_jc48_sep08_quartic_common_root.md),
the [shared-root theorem](planar_jc48_sep08_shared_roots.md),
and the weighted field and primitive-space methods of
[octuple transport](planar_jc48_sep08_octuple_transport.md).
The canonical hostile is `F=t^4`, which has a rational mate
but no polynomial one. The corrected near miss is to apply an
unweighted source theorem after inversion. The live concepts
are the true volume weight, full branch exhaustion, local
residues, a genus-two holomorphic form, and the original
polynomial coordinate ring.

If `L` is constant, the original polynomial identity
`J(F,G)=2H*J(H,G)` already excludes polynomial mates. Assume
from now on that `L` is nonconstant.

## 2. The global exactness test with the correct weight

At the only boundary zero, use the new finite coordinates
`(u,s)`. The original generic fibre and relative form are

    E=N^2+s^3*M-c*s^4=0,
    eta=u^2*s^2*du/E_s.                              (3)

Away from this boundary point the compact generic fibre meets
no other point of `S`. Inside `W`, its relative form is
regular, including the new `D` by (2). Generic fibres are
smooth there; no generic component is contained in the fixed
zero divisor of (2). Therefore the form is nonzero on every
generic component. If (3) is holomorphic on every normalized
branch at the boundary point, rational exactness is impossible:
an exact rational differential without poles on a compact
component must vanish. If a normalized branch instead has a
simple pole with nonzero residue, exactness is also impossible.
Each local case below supplies one of these two alternatives.

If `lt!=0`, `M` is a unit and `m=8`. The unweighted M-unit
classification is regular for every possible integer normal
order `j`: the balanced equation `8=3j` never occurs. For
`j=0` the normal derivative is a unit and there are only two
local sheets. For `j>=1`, the local degree is three; the
simple low face and cancellation face when `8>3j`, or the
dominant cubic face when `8<3j`, account for every sheet.
Multiplication by `u^2` improves every order. This proves
`lt=0` without invoking the unweighted global exclusion.

If `kt!=0`, the shared-root normal-unit lemma again makes
all unweighted local forms regular. The weighted form is
regular, and `kt=0` follows. From this point onward,

    E(0,s)=(k0^2+l0-c)*s^4,                          (4)

so the exact local Weierstrass degree is four.

## 3. All lower normal jets are paid, including the balanced case

Suppose `kxt!=0`, so the first normal order is `j=1`.
The two low branches have `s` of order one and arise from

    Z^2[(kxt+k0*Z)^2+lxt*Z+(l0-c)*Z^2],  s=u*Z.

The bracket has two simple nonzero roots generically. Each
unweighted form is logarithmic, so (3) is regular of order
one. The other two determinations lie at `s` of order seven,
where the leading terms of `N` cancel. If
`n=ord(lxt*u+l2*u^2+l3*u^3+l4*u^4)`, then `1<=n<=4`.
Here `N_s` has order one and the weighted form has exponent
`(9-n)/2` in `u`. It is regular on the actual normalization,
including odd split orders. These determinations plus the
two low branches exhaust (4). Therefore `kxt=0`.

Now suppose `lxt!=0`, so `n=1`. The possible normal orders
are `j=2,3,4,6`, with the coefficient `2alpha*u^6` present
if all smaller normal coefficients vanish.

- For `j=2`, there is one low branch `s~u`, one middle
  branch `s~u^3` from the face
  `Z^2(k2^2+lxt*Z)`, and two cancellation determinations
  with `s~u^6`. The first two unweighted forms are logs,
  hence weighted-regular. The last weighted exponent is
  `5/2` in `u`, regular on the ramified normalization.

- For `j=3`, besides the low `s~u` branch, set `s=u^5*Z`.
  The remaining cubic is

      P(Z)=(alpha+k3*Z)^2+lxt*Z^3.

  It has no zero root and cannot have a triple root. Simple
  roots give weighted forms of order one. At a double root
  `q!=0`, let `Q=E/u^16`, and take the actual analytic
  critical centre of `Q` in `Z`. Its critical value `R`
  obeys `R_c=-u^4*z(u,c)^4`; thus its generic order
  `lambda` is at most four. Parametric Morse normalization
  makes (3) a unit times `u*du/sqrt(-R)`. Orders one,
  two, and three are regular after normalization. Order
  four gives two genuine simple poles with nonzero
  coefficients, so each has a nonzero residue. The simple
  cubic root and double-root group account for the other
  three determinations, including all generic contacts.

- For `j=4` or `6`, the remaining three branches at
  `s=u^5*Z` have face `alpha^2+lxt*Z^3`, with three
  simple nonzero roots. Their weighted forms have order
  one. Together with the low branch they exhaust (4).

The no-triple assertion follows directly by comparing the
coefficients of `(alpha+k3 Z)^2+lxt Z^3` with a cubic
cube: the required coefficient equalities would give
`9=12`, with nonzero `alpha,k3,lxt`. Higher jets enter the
critical value but cannot evade the generic order-four bound.
Thus every `lxt!=0` case is holomorphic or has an unavoidable
nonzero residue. It follows that `lxt=0`.

We have reached the complete reduced entry

    H=alpha*h^2+k0+k2*w+k3*v+k4*h,
    L=l0+l2*w+l3*v+l4*h.                             (5)

If `k2!=0`, the complete four-branch calculation has two low
branches with unweighted double poles; (3) makes both regular.
For the other two determinations, `s~u^6` and
`n=ord(l2*u^2+l3*u^3+l4*u^4)` is two, three, or four.
The weighted exponent is `3-n/2`, giving normalized orders
two, four, and one, respectively. Every branch is regular,
and (4) proves completeness. Therefore `k2=0` as well.

## 4. A genus-two holomorphic form or a nonzero logarithmic residue

For the remaining functions use their exact rational field:

    h=u^2+u^4*T,    v=u+u^3*T,
    u=h/v,         T=(v^2-h)*v^2/h^3,
    omega_old=(1/h)*dh wedge dv,
    w=v^2/h-1.                                      (6)

In the original source these are simply `h=-t` and `v=-x*t`.
Put `A=alpha*h^2+k4*h+k0`. The complete original function is

    F=a(h)*v^2+2b(h)*v+d(h),
    a=k3^2+l2/h,    b=k3*A+l3/2,
    d=A^2+l4*h+l0-l2.

When `a!=0`, set `Y=h(av+b)`. The exact generic equation
and relative form become

    Y^2=Q(h)=h^2[b^2+a(c-d)],
    eta=-dh/(2Y).                                   (7)

If `l2!=0`, the polynomial `Q` has degree five with leading
coefficient `-l2*alpha^2`. Its coefficient of `c` is
`B(h)=h(k3^2*h+l2)`, which has only simple roots (or just
the one simple root when `k3=0`). Generic squarefreeness is
complete: away from `B=0`, a repeated root must be a critical
point of `Q0/B`, where `Q0=Q-cB`. This rational function is
nonconstant because `deg Q0=5>deg B`; hence it has only
finitely many critical points. At a root of `B` that is also
a persistent zero of `Q0`, the derivative has nonzero `c`
coefficient `B'`, so it is simple generically. All other
roots of `B` are not roots of `Q`.

Thus the compact generic normalization is a smooth genus-two
hyperelliptic curve. The form `dh/Y` is holomorphic and
nonzero: at a finite branch point the orders of `dh` and
`Y` are both one; at infinity `ord h=-2`, `ord Y=-5`,
so the form has order two. There are no other poles.
It cannot be the derivative of a rational function. This
excludes all `l2!=0` cases, including every `k3,l3` value.

There is an elementary polynomial-only check at this boundary:
in the original variables, on `t=0`,

    F_x=0,
    F_t=-l2*x^2-(2k0*k3+l3)*x-(2k0*k4+l4).

When `l2!=0` this has a zero over `C`. The genus-two argument
is the stronger rational-mate exclusion; it is not presented
as the only way to detect a polynomial obstruction here.

If `l2=0` but `k3!=0`, put `y=k3^2*v+k3*A+l3/2`.
Then `y^2=b^2+k3^2(c-d)`, whose value at `h=0` is generically
nonzero because its coefficient of `c` is `k3^2`. The form is
`-dh/(2h*y)`. On both points over `h=0` it has a nonzero
logarithmic residue. This remains valid if the quadratic
extension becomes a constant extension after closing the
generic constants; no geometric-integrality assumption is
needed to detect those residues.

Finally let `l2=k3=0`. If `l3!=0`, the actual field is linear
in `v` and the relative form is `-dh/(l3*h)`, with nonzero
residue. If `l3=0`, then `F=P(h)` with quartic `P`, but
`h=-t` in the original source. The nonconstant factor
`P'(-t)` in the original polynomial Jacobian excludes every
polynomial mate. This final step does not assume that the
transformed rational mate is polynomial.

These cases exhaust (1), proving the theorem.

## 5. Scope and verification boundary

The original example `F=t^4`, `G=-x/(4t^3)` has Jacobian
one and belongs to the infinity-octuple class. Its mate has
poles. Thus the full conclusion must remain polynomial,
even though the noncomposite cases above exclude rational
mates. The actual surface involution preserves rational
exactness with the weight (2); it does not preserve the
old source polynomial ring.

Together with the all-finite-octuple theorem, this
excludes a global quartic whose canonical boundary octic
has just one zero anywhere on `P1`. Other multiplicity
partitions and the general quartic stratum remain outside
the claim. This is not a proof of JC(2).

## 6. Exact verification and reproduction

The new standalone [source](../../04-computation/planar_jc48_sep08_infinity_octuple.py)
and [frozen output](planar_jc48_sep08_infinity_octuple.out) pass
**141 always-active exact gates**. Ordinary and optimized Python outputs
are byte-identical. No inherited mathematical implementation is imported.

The exact universe comprises the full six-dimensional `L1` correction
space, all lower-jet strata in (1), every integer normal order supplied
by its full numerator, the complete normalized order ledger through the
balanced contact bound, and the complete remaining quadratic coefficient
family. The source checks the actual surface involution, weighted volume
on both charts, the rank-nine octic restriction map, all symbolic leading
faces and generic simple-root suppliers, and the exact full-field
hyperelliptic equation. Its named genus-two controls include a persistent
root of both `Q0` and `B`, preventing a false coprimality shortcut.
It also retains the original rational-mate hostile `t^4` and the direct
original-source criticality control for `l2!=0`.

These algebraic controls support the analytic branch-exhaustion and
compact exactness arguments; a finite sample is not substituted for either.
The no-rational conclusions are local or geometric, whereas the constant
`L` and final composite cases explicitly use the original polynomial ring.

Reproduction from the worktree root:

```bash
python3 04-computation/planar_jc48_sep08_infinity_octuple.py
python3 -O 04-computation/planar_jc48_sep08_infinity_octuple.py
```

Frozen SHA256 pins:

- Source, 12,949 bytes:
  `8fc3fcf4c719deb653970090d8ca63589a2a3649074639d59cf96ef764bde333`.
- Output, 399 bytes:
  `96805cf8e6722abd373477ed3ab6f1af62b83c50b7c18d1a0c7aba476978a33e`.
- Semantic digest:
  `a48a1cd3da3106bac97f1c9629c8063283783fd842be8aea07398be46ca89718`.

The [independent audit](planar_jc48_sep08_infinity_octuple_audit.md) accepts the complete analytic proof, source and both141-gate replays. It separately reconstructs the old global basis under inversion and the original quadratic field identities in ten exact checks. Root owns integration; no canonical ID is reserved here.
