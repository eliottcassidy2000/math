# Every finite location of the octuple quartic boundary class is excluded

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.**
The conclusion is no polynomial source mate of any degree. The
surface is fixed throughout; the proof changes its global
functions and does not assert a surface translation symmetry.

## 1. Actual global functions at an arbitrary finite point

Use `W=(P1_x x P1_z) minus {z=x^2}`, `s=z-x^2`, `t=1/s`, and
the actual chart `x=1/r`, `t=-r^2-r^4*b_D`. Let `H=N/s^2` and
`L=M/s` be global functions in `L_2,L_1`, and set `F=H^2+L`.

**Theorem.** If the complete binary boundary octic is
`N|S=alpha*(x-p)^8`, with `alpha!=0` and finite `p in C`, then
there is no `G in C[x,t]` with `J_(x,t)(F,G)=lambda!=0`.
No degree bound on the proposed mate is imposed.

Put `u=x-p` and define the following *actual global* functions:

    w=u^2*t,
    v=u+u^3*t-2p,
    h=u^2+u^4*t-2p*u=u*v.                             (1)

The numerator of `h` is `u^4+s(u^2-2p*u)` and that of `v` is
`u^3+s(u-2p)`. Substitution `s=z-(u+p)^2` places both in the
full `O(2,1)` section box. Their boundary restrictions are

    h_D=-b_D-3p^2,       v_D=0,       w_D=-1.          (2)

The complete `L_1` space has basis `1,t,ut,w,v,h`, and `h^2`
has octic `u^8`. Dividing the difference of two numerators with
the same restriction by the defining section of `S` therefore
gives the complete entry

    H=alpha*h^2+k0+kt*t+kxt*u*t+k2*w+k3*v+k4*h,
    L=l0+lt*t+lxt*u*t+l2*w+l3*v+l4*h.                 (3)

The terms `-2p` and `-2p*u` in (1) are essential. Omitting them
would use functions that generally fail to be global. The
ordinary translation of `u^2+u^4*t` has a pole `2p/r` on `D`;
the second correction in (1) cancels it. This is different from
the section-space intersection that paid
[the finite-six translation](planar_jc48_sep08_const_d_translation.md).

The closest proved mechanisms are the
[fixed finite-eight theorem](planar_jc48_sep08_finite_eight.md),
the [shared-root first-jet theorem](planar_jc48_sep08_shared_roots.md),
and the exact primitive-space operation in
[the constant-D proof](planar_jc48_sep08_const_d_quartic.md).
The live concepts are global functions, normalized branch
capacity, the original volume form, field trace, and a complete
space of allowed primitives. The corrected near miss is to
transport polynomial regularity through a rational source map.
All uses of such a map below preserve rational exactness only.

The case `p=0` is already the fixed finite-eight theorem. The
new final argument below treats `p!=0`. If `L` is constant,
`J(F,G)=2H*J(H,G)` directly excludes a polynomial mate. Assume
otherwise that a polynomial, hence rational, mate exists.
Divide it by its nonzero Jacobian constant so that the relative
exactness equations below have coefficient one.

## 2. Exact local reductions and the four-branch capacity gate

The M-avoidance theorem forces `lt=0`, since the octic has only
one boundary zero. A nonzero `kt` makes `N_s` a unit at that
point; the shared-root regularity theorem then leaves no pole
on the compact generic fibre and excludes rational exactness.
Thus `kt=0`.

Write `kappa=k0-2p*k3`, `ell=l0-2p*l3`. The actual finite
first-jet polynomial after `s=u*Z` is

    Z^2[(kxt+kappa*Z)^2+lxt*Z+(ell-c)*Z^2].           (4)

The nonzero simple-root/logarithmic-residue argument from the
fixed-point proof applies to this exact polynomial, forcing
`kxt=lxt=0`. This uses a legitimate local coordinate at the
specified point; no global coordinate-change assertion enters.
The remaining functions are

    H=alpha*h^2+k0+k2*w+k3*v+k4*h,
    L=l0+l2*w+l3*v+l4*h.                             (5)

If `k2!=0`, the full four-branch argument also survives with
its actual hypotheses unchanged. The low face after `s=u^2*Z`
is

    Z^2[k2^2+(2*kappa*k2+l2)*Z+(kappa^2+ell-c)*Z^2]. (6)

The bracket has two simple nonzero roots generically, each giving
a genuine double pole of the relative form and possible primitive
degree one. For the remaining two determinations, `s=u^6*Z`
gives `N/u^8=alpha+k2*Z+O(u)`. The exact centre `N=0` has
`Z(0)=-alpha/k2`; `N_s` has order two. If
`n=ord(l2*u^2+l3*u^3+l4*u^4)`, then `n=2,3,4`, since `L`
is nonconstant. The split has

    ord N=9+n/2,   ord E_s=11+n/2,
    eta=(unit)*u^(1-n/2)*du.

For `n=2` the two branches are regular; for `n=3` the single
normalized branch `u=tau^2` is regular; for `n=4` there are
two nonzero logarithmic residues. Higher `p` terms change none
of these orders or their nonzero leading coefficients. The
restriction `E(0,s)=(kappa^2+ell-c)*s^4` pays all four sheets.

For `n=2,3`, total possible primitive pole degree is at most
two. By (2), `F_D` has degree four with leading coefficient
`alpha^2`; a generic point of `D` on the original fibre has
relative-form order two and would require primitive local degree
three. Its own compact component cannot have more poles than
the total budget two. This is the same-component contradiction,
with no generic irreducibility assumption. Therefore `k2=0` in
every remaining polynomial-mate candidate.

## 3. The changed functions preserve an exact Poisson field

The functions (1) give the literal field identities

    u=h/v,
    t=(v^2+2p*v-h)*v^2/h^3,
    J_(u,t)(h,v)=h^3/v^2,
    omega=v^2/h^3 * dh wedge dv,
    w=(v^2+2p*v)/h-1.                               (7)

These preserve the original rational field and volume form.
They do not preserve polynomiality of a proposed mate: the
inverse has denominators. Every rational exclusion obtained
from (7) is valid for the original source; a polynomial
exclusion requires an additional argument.

Set `A(h)=alpha*h^2+k4*h+k0`. For `k2=0`, the full original
function becomes

    F=a(h)*v^2+2b(h)*v+d(h),
    a=k3^2+l2/h,
    b=k3*A+l3/2+p*l2/h,
    d=A^2+l4*h+l0-l2.                               (8)

When `a` is nonzero, the generic quadratic extension over
`K(h)`, `K=C(c)`, is irreducible: its discriminant is
`b^2+a(c-d)`, and the prime `c-(d-b^2/a)` has odd valuation
in `C(h)(c)`. This includes constant-extension cases without
assuming geometric integrality after closing the constants.

The actual form and its field trace under
`v -> -2b/a-v` are

    eta=-v^2/[2h^3(av+b)] * dh,
    Tr(eta)=2b/(a^2*h^3) * dh
       =[2k3*h*A+l3*h+2p*l2]/[h^2(k3^2*h+l2)^2] * dh. (9)

If `k3!=0`, its residue at `h=infinity` is
`-2alpha/k3^3`, which is nonzero. Trace commutes with the
relative derivative fixing `K`; an exact rational form would
have an exact rational trace with zero residues. Thus this
entire case has no rational mate, for every value of `l2`.

It remains to take `k3=0`. If `l2!=0`, the trace is

    [l3/(l2^2*h)+2p/(l2*h^2)] * dh,

whose residue forces `l3=0`. Unlike at `p=0`, this condition
need not make the original source critical. The following
primitive-space calculation supplies the missing obstruction.

## 4. The genus-two residual has no rational primitive

Assume `p*l2!=0`, `k2=k3=l3=0`, and put

    P(h)=(alpha*h^2+k4*h+k0)^2+l4*h+l0-l2,
    y=v+p.

The original fibre and relative form are exactly

    y^2=Q(h)=p^2+h(c-P(h))/l2,
    eta=-(y-p)^2/(2l2*h^2*y) * dh.                  (10)

The polynomial `Q` has degree five and leading coefficient
`-alpha^2/l2`. For generic `c` it is squarefree. A repeated
root cannot be zero because `Q(0)=p^2`; any nonzero repeated
root satisfies

    h^2*P'(h)+l2*p^2=0,

a nonzero polynomial independent of `c`. Each of its finitely
many roots determines at most one exceptional fibre value.
Thus the compact generic normalization is geometrically
integral of genus two. Its points over `h=0` are the two
distinct unramified points `Pplus=(0,p)`, `Pminus=(0,-p)`.
There is one point at infinity with `ord h=-2`, `ord y=-5`.

At `Pplus`, (10) is regular. At `Pminus` it has a double pole
with leading coefficient `2p/l2`; the residue is zero. At
infinity it has order minus four. At any other finite branch
point, `h-h0` has order two and `y` order one, so the form is
regular. These are all its poles. A rational primitive could
therefore have only a simple pole at `Pminus` and a pole of
order at most three at infinity.

The complete allowed function space has basis

    1,       h,       (y-p)/h.                       (11)

Here is a direct completeness proof, independent of a dimension
guess. Apply the hyperelliptic involution to an allowed function.
Its even part can have at most simple poles at both points over
zero and pole order at most three at infinity; as a rational
function of `h`, it has form `A0+B0*h+C0/h`. Its odd part is
`y*S(h)`. It can have at most a simple pole at `h=0`, none at
other finite points, and must decay at least as `h^-1` at
infinity because `y` has pole order five. Hence it is `D0*y/h`.
Regularity of the sum at `Pplus` forces `C0=-p*D0`, giving (11).
The three functions are independent; equivalently this is the
degree-four Riemann--Roch space on this genus-two curve.

Allow their coefficients to lie in the algebraic closure of the
generic constant field. If

    G=A0+B0*h+C0*(y-p)/h,

then comparing traces of `dG` with (10) forces
`B0=0`, `C0=1/l2`. The remaining exact identity is

    eta-(1/l2)*d[(y-p)/h]=P'(h)/(2l2^2*y) * dh.       (12)

It cannot vanish because `P` has degree four with leading
coefficient `alpha^2`. This excludes every rational mate in
the residual case. The trace alone was exact here; retaining
the full allowed primitive space supplies the decisive test.

## 5. The last linear case retains a same-fibre pole-repair obstruction

The remaining case has `k2=k3=l2=0`, so

    F=P(h)+l3*v,        deg P=4, leading coefficient alpha^2.

If `l3=0`, the nonconstant factor `P'(h)` in
`J(F,G)=P'(h)*J(h,G)` excludes every polynomial mate. This
does not require `h` to have a source critical point; for
`p!=0`, it is in fact smooth on `u=0`.

Suppose `l3!=0`. The pair `(F,h)` generates the rational source
field. For a mate normalized to Jacobian one,

    G_h|F=-(F-P(h))^2/(l3^3*h^3).

Write `P=P0+P1*h+P2*h^2+...` and `g=F-P0`. The residue is
`-(P1^2-2g*P2)/l3^3`, so rational integrability forces
`P1=P2=0`. Set `q=P-P0 in h^3*C[h]`. An exact primitive is

    G0=l3^-3 [g^2/(2h^2)+2g*integral(q/h^3)dh
                                      -integral(q^2/h^3)dh].

Substitution `g=l3*v+q` shows

    G0=1/(2l3*u^2)+a polynomial in h,v.              (13)

Every rational mate is `G0+R(F)` with `R in C(T)`: the kernel
of the nonzero field derivation on `C(F)(h)` is exactly `C(F)`.

On the original source line `E0={u=0}`, `v=-2p`, `h=0`, hence
`F_E=P0-2p*l3`. The full original source fibre factors as

    F-F_E=u*B(u,t),
    B=l3*(1+u^2*t)+q(u*v)/u,
    B|E0=l3!=0.                                      (14)

The quotient in (14) is polynomial since `q` is divisible by
`h^3`. It is nonconstant: its `t^4` coefficient is
`alpha^2*u^15`. Therefore its zero set is a nonempty source
divisor disjoint from `E0`. The primitive (13) is regular
there. To repair its double pole on `E0`, `R(F)` must have a
pole at `F_E`; that same pole produces a pole along the other
divisor `B=0`, where (13) cannot cancel it. Thus no polynomial
mate exists. This pays the changed special-fibre value and
does not import an auxiliary surface's polynomial regularity.

The cases above exhaust (3). Together with the proved `p=0`
case they establish the theorem for every finite `p`.

## 6. Hostile controls and exact scope

Inside this very class, for every `p`,

    F=h^4+h,      G=1/[3u^3(4h^3+1)]

has Jacobian one. Its mate has poles, so the main conclusion
must remain polynomial. The example `F=h^4-v` also has a
rational mate and is free of source critical points; the
same-fibre repair obstruction, rather than criticality alone,
excludes its polynomial mates.

The source-to-target map (7) preserves the full rational field,
the original fibre parameter, and its Poisson coefficient. It
loses polynomial regularity, which Section 5 restores as a
separate proof obligation. The generic trace in Section 3 can
lose geometric component data, but exactness survives. Section
4 uses the actual squarefree genus-two model and proves its
entire primitive space, rather than inferring exactness from
zero residues or from a permissible pole count.

Other octic multiplicity partitions and the distinguished
infinity-octuple case are not covered. This is not a proof of
JC(2) or of the entire quartic stratum.

## 7. Exact reproduction and verification boundary

```sh
python3 -B 04-computation/planar_jc48_sep08_octuple_transport.py
python3 -B -O 04-computation/planar_jc48_sep08_octuple_transport.py
```

The standalone source imports no inherited mathematical
implementation. It checks the complete six-column global L1
basis for symbolic `p`, the original section boxes and chart
restrictions, all three high-branch orders with the complete
derivative remainder, both rational field inverses and the volume
form, the full quadratic trace and both residue locations,
the genus-two squarefreeness supplier and every pole, the exact
primitive-space remainder, and the actual shifted same-fibre
factor. Three named smooth genus-two fibres and two genuine
inside-class rational mates are separate controls. None of these
finite controls substitutes for completeness of the analytic
cases or of the allowed primitive space.

Normal and optimized runs agree byte for byte with the frozen
output: **85 always-active exact gates, 398 output bytes**.

- [Source](../../04-computation/planar_jc48_sep08_octuple_transport.py):
  9,958 bytes, SHA256
  `b56bc24e5e79233ba74b84d3565a582dbe8d6d695c774a255fb1a175ba73fd87`.
- [Frozen output](planar_jc48_sep08_octuple_transport.out): SHA256
  `cf67f2ec2c4e68f1b9cdd2c610634ae2ddbb6a091d349484c1a1dfa97ad7aae8`.
- Semantic digest:
  `501257095e3b0f36f4657438ef9b4dd59de322b7d4c88e6978d4ab24f6497560`.

Source, output and the accepted independent analytic/source audit are frozen.

The [independent audit](planar_jc48_sep08_octuple_transport_audit.md) accepts the complete all-finite-point proof, original-source pole repair, entire source and both85-gate replays. Its independent section, field-trace and primitive-space reconstructions agree.
