# Independent audit of the all-finite octuple quartic theorem

**Status: full independent analytic/source audit PASS; normal, optimized
and frozen output agree at 85 always-active exact gates.**
This audit concerns [the complete octuple-transport proof](planar_jc48_sep08_octuple_transport.md).
The conclusion is absence of polynomial Jacobian mates of unrestricted
degree when the complete boundary octic is `alpha*(x-p)^8`, with
`alpha!=0` and finite `p`. The original DG surface is fixed. Neither
ordinary translation as a surface automorphism nor absence of all
rational mates is asserted.

## 1. The changed functions give the entire global entry

Set `u=x-p`, but keep the original defining function
`s=z-(u+p)^2` and original plane variable `t=1/s`. The actual carriers
are

    w=u^2*t, v=u+u^3*t-2p, h=u^2+u^4*t-2p*u=u*v.

Their numerators after substitution for `s` lie in the full `(2,1)`
section box. I independently formed the nine-by-fifteen restriction
matrix on all `u^i z^j`, `0<=i<=4`, `0<=j<=2`, under
`z=(u+p)^2`. It has rank nine for every `p`. Its full kernel is
`s` times the six-dimensional `(2,1)` box. The numerator of `h^2`
restricts to `u^8`; hence every numerator with the prescribed octic
is `alpha` times that numerator plus precisely such a kernel element.

Independently expressing the actual numerators of `1,t,ut,w,v,h`
in the complete `(2,1)` box gives determinant `-1`, identically in
`p`. Consequently the two six-parameter lower terms in the primary
are complete, including every normal and tangential coefficient. This
is a global section argument, not a local sufficient family.

On the actual second chart `u=1/r-p`, `t=-r^2-r^4*b_D`, direct
expansion gives

    h_D=-b_D-3p^2, v_D=0, w_D=-1.

The uncorrected function `u^2+u^4*t` instead has pole term `2p/r`.
The correction `-2pu` removes exactly that pole. This provides an
explicit hostile to a surface-translation shortcut, and confirms that
the proof is changing the global functions rather than silently moving
the graph. The complete binary octic has only the one finite zero;
its value at the distinguished point at infinity is nonzero.

The `p=0` dependency is the independently accepted
[fixed finite-eight theorem](planar_jc48_sep08_finite_eight.md), with
[its full audit](planar_jc48_sep08_finite_eight_audit.md). Its mathematical
proof and source were accepted independently. During this read, its
header still had a stale RESERVED label; the owner repaired that status
without changing the mathematics. Before freezing this audit I re-read
its PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED header. Thus this
audit does not inherit a RESERVED result as a proved dependency. The
other used local and pole-degree suppliers were also read in their
current proved versions.

## 2. The local reductions preserve their actual hypotheses

If `L` is constant, `J(F,G)=2H J(H,G)` rules out polynomial mates,
with no rational exclusion inferred. Otherwise common-root necessity
forces `lt=0`; `kt!=0` is the normal-unit regularity case at the
only possible boundary point. A compact generic component then carries
a regular, nonzero relative form, so cannot have an exact rational
primitive. This pays the first normal reduction.

The changed constant terms are exactly

    kappa=k0-2p*k3, ell=l0-2p*l3.

After `s=uZ` the tangent polynomial is

    Z^2[(kxt+kappa Z)^2+lxt Z+(ell-c)Z^2].

Its nonzero simple roots give nonzero logarithmic residues exactly as
in the fixed-point calculation. The proof uses a valid local coordinate
at the specified point and forces `kxt=lxt=0`; it does not claim
that the zero roots are already normalized branches.

For `k2!=0`, the two low branches under `s=u^2 Z` have the complete
face

    Z^2[k2^2+(2kappa*k2+l2)Z+(kappa^2+ell-c)Z^2].

The remaining quadratic has nonzero constant term and discriminant
with nonzero `c` coefficient `4k2^2`. Each actual simple branch has
a double pole of the relative form, hence at most one unit of primitive
pole degree.

I also checked the full normal coefficient, including the new terms:

    [s]N=k2*u^2+k3*u^3+k4*u^4-4alpha*p*u^5+2alpha*u^6.

Thus in the second chart `s=u^6 Z`, the exact equation `N=0`
still has the unique unit-derivative centre `Z(0)=-alpha/k2`.
For `n=ord(l2*u^2+l3*u^3+l4*u^4)`, nonconstant `L` implies
`n` is exactly 2, 3 or 4. Every `s`-dependent term of `M-cs`
has higher order there. The exact derivative comparison gives

    ord_u N=9+n/2, ord_u N_s=2,
    ord_u E_s=11+n/2,
    eta=(unit)*u^(1-n/2) du.

All other terms of `E_s` have order at least `12+n`, strictly
larger than the displayed dominant term. The `p`-dependent terms
cannot cancel these nonzero suppliers. For `n=2` there are two
regular branches; for `n=3` one normalized branch with `u=tau^2`
has regular differential; for `n=4` there are two nonzero logarithmic
residues. The generic equality
`E(0,s)=(kappa^2+ell-c)s^4` proves that the low branches and high
determinations exhaust the four-sheet local equation.

There are no other boundary intersections. For `n=2,3`, the total
possible primitive pole degree is at most two. The actual restriction
`F_D` is quartic with leading coefficient `alpha^2`, for every
remaining coefficient choice. A generic transverse point of `D` has
relative differential of order two, requiring primitive local degree
three on its own compact component. That component's pole degree is
at most the total budget two. This contradiction is componentwise;
no geometric-integrality or hidden constant-D assumption is used.

## 3. Weighted quadratic trace in the original field

Direct differentiation and inversion give

    u=h/v, t=(v^2+2pv-h)*v^2/h^3,
    J_(u,t)(h,v)=h^3/v^2,
    omega=v^2/h^3 dh wedge dv,
    w=(v^2+2pv)/h-1.

The affine coordinate `u=x-p` has Jacobian one, so this is the
original source volume. These equalities identify the full rational
field; they do not preserve polynomial regularity of an arbitrary
proposed mate.

For `k2=0`, the complete function is

    F=a(h)v^2+2b(h)v+d(h),
    a=k3^2+l2/h,
    b=k3*A+l3/2+p*l2/h,
    d=A^2+l4*h+l0-l2,
    A=alpha*h^2+k4*h+k0.

Whenever `a!=0`, the discriminant has a prime factor
`c-(d-b^2/a)` of valuation one in `C(h)(c)`. Thus the generic
quadratic extension over `K(h)`, `K=C(c)`, is a field, including
cases which enlarge the relative constant field. No geometric
integrality after closing `K` is assumed in the trace step.

The actual relative form is

    eta=-v^2/[2h^3(av+b)] dh.

I independently derived its trace from the symmetric relation
`v+v'=-2b/a`, with `v'=-2b/a-v`, obtaining

    Tr(eta)=2b/(a^2*h^3) dh
      =[2k3*h*A+l3*h+2p*l2]/[h^2(k3^2*h+l2)^2] dh.

For `k3!=0`, the infinity residue is `-2alpha/k3^3!=0`,
regardless of `l2`. Trace commutes with the relative derivation in
characteristic zero, while a rational derivative on `K(h)` has zero
residue. This therefore excludes rational mates in that entire
coefficient stratum, with constants-field degenerations retained.

For `k3=0`, `l2!=0`, the trace is
`[l3/(l2^2 h)+2p/(l2 h^2)]dh`. Its residue forces `l3=0`.
The remaining trace is exact; stopping at trace would not prove the
needed exclusion. The primary correctly proceeds to the full primitive
space instead of using a source-criticality inference that fails when
`p!=0`.

## 4. Complete genus-two model and primitive space

For `p*l2!=0` after those reductions, let

    P(h)=(alpha*h^2+k4*h+k0)^2+l4*h+l0-l2,
    y=v+p, Q(h)=p^2+h(c-P(h))/l2.

The exact original fibre and form are

    y^2=Q(h),
    eta=-(y-p)^2/(2l2*h^2*y) dh.

The polynomial `Q` has degree five, leading coefficient
`-alpha^2/l2`, and nonzero constant term `p^2`. Any repeated
finite root is nonzero and satisfies
`h^2 P'(h)+l2*p^2=0`. This nonzero polynomial has finitely many
roots, each specifying at most one exceptional value of `c`.
Therefore the generic polynomial is squarefree. Its double cover is
geometrically integral and has six simple branch points, including
infinity; Riemann--Hurwitz gives genus two. This genuine component
statement is paid here and was not needed in the earlier trace case.

The points at zero are distinct unramified points `Pplus=(0,p)`
and `Pminus=(0,-p)`. At `Pplus`, the numerator makes `eta`
regular. At `Pminus` its leading coefficient is `2p/l2` times
`dh/h^2`; expansion of `y+p^2/y-2p` shows its linear term is
zero, so the double pole has zero residue. At every other finite
branch point, `dh` cancels the simple zero of `y`, leaving a
regular differential. At the unique point at infinity,
`ord h=-2`, `ord y=-5`, `ord dh=-3`; hence `ord eta=-4`.
These exhaust the poles. A primitive could have at most a simple
pole at `Pminus` and a pole of order three at infinity.

The claimed complete space is correct:

    L(Pminus+3 infinity)=span{1,h,(y-p)/h}.

For clarity, the proof does not infer this just from its dimension.
Take even and odd parts under `y -> -y`. Each part may acquire a
simple pole at either point over zero, but no other finite pole and
no pole of order greater than three at infinity. The even part,
a rational function of `h`, is `A0+B0*h+C0/h`. The odd part is
`y*S(h)`. At other finite branch points a pole of `S` would give
an odd pole of the product, so none is allowed. At zero `S` has at
most a simple pole. At infinity, the bound combined with the order
five pole of `y` forces `S=O(1/h)`. Thus `S=D0/h`, and
regularity at `Pplus` imposes `C0=-pD0`. These arguments also
exclude cancellation-based extra basis elements. The three displayed
functions are independent. Coefficients may be taken over an algebraic
closure of the generic constant field.

After normalizing a hypothetical nonzero Jacobian to one, trace of the
derivative of `A0+B0*h+C0*(y-p)/h` forces `B0=0`, `C0=1/l2`.
Independently reducing the remaining anti-invariant part using `y^2=Q`
gives exactly

    eta-(1/l2)d[(y-p)/h]=P'(h)/(2l2^2*y) dh.

The right side cannot vanish because `P` has degree four and
`[h^3]P'=4alpha^2!=0`. This is the required rational obstruction
in the trace-zero case; zero residues and an admissible pole budget
alone would not suffice.

## 5. Original-source pole repair in the last linear case

The last possibility is `F=P(h)+l3*v`, with `P` of degree four.
For `l3=0`, the polynomial factor `P'(h)` directly excludes a
polynomial constant-Jacobian mate. This is not a claim that the
function `h` is critical on `u=0`; it has derivative `-2p` there.

For `l3!=0`, the full rational field is `C(F)(h)`. Writing
`g=F-P0` and expanding the exact derivative at fixed `F` gives
residue `-(P1^2-2gP2)/l3^3`; thus rational integrability forces
`P1=P2=0`. For `q=P-P0 in h^3 C[h]`, integration at fixed `g`
gives the primary's `G0`. Substitution `g=l3*v+q` proves

    G0=1/(2l3*u^2)+a polynomial in h,v.

In particular the pole on the original source line `E0={u=0}`
is exactly double, with nonzero coefficient. The derivation constants
are exactly `C(F)`, so every rational mate is `G0+R(F)`. No
additional arbitrary source-dependent repair is allowed.

The actual fibre value on that line is `F_E=P0-2p*l3`, and its
full source factorization is

    F-F_E=u*B(u,t),
    B=l3(1+u^2t)+q(uv)/u,
    B(0,t)=l3.

The quotient is polynomial and nonconstant; its coefficient of `t^4`
is `alpha^2*u^15`. Therefore `B=0` is a nonempty divisor in the
original affine plane, disjoint from `E0`. The primitive `G0` is
regular there. Since `F-F_E` has order exactly one at `E0`, a
repair of its double pole would require `R` to have a pole of order
two at `F_E`. That produces a pole along every component of `B=0`
which `G0` cannot cancel, regardless of multiplicities in that divisor.
This is a polynomial-mate obstruction on the original source, not an
imported regularity assertion from an auxiliary rational chart.

Together with the fixed-zero theorem, all coefficient cases and every
finite `p` are covered. The hypotheses still do not include other
octic partitions or the distinguished infinity-octuple location.

## 6. Controls and source verification

The all-p rational-mate controls are correctly inside the stated class:

    F=h^4+h, G=1/[3u^3(4h^3+1)],

and

    F=h^4-v, G=-1/(2u^2)+2v*h^2-(4/3)h^6.

Both have literal Jacobian one. The second is source-critical-free:
on `u=0` its `u` derivative is `-1`; on `v=0`, `u!=0`, its
`t` derivative is `-u^3`; elsewhere the `(h,v)` chart has nonzero
Jacobian and the `v` derivative of `h^4-v` is `-1`. They prohibit
strengthening the main conclusion to rational mates or reducing the
entire exclusion to source criticality.

I read the entire standalone source. It imports no inherited
mathematical implementation and retains every coefficient in the actual
global entry. Its named fibre checks supplement the symbolic generic
squarefreeness argument; they are not used to infer generic genus from
samples. Likewise the arithmetic dimension check in the source is not
a replacement for the complete even/odd primitive-space proof above.

Independent reproduction commands:

    python3 -B 04-computation/planar_jc48_sep08_octuple_transport.py
    python3 -B -O 04-computation/planar_jc48_sep08_octuple_transport.py

Frozen pins declared by the producer:

- [Source](../../04-computation/planar_jc48_sep08_octuple_transport.py),
  9,958 bytes, SHA256
  `b56bc24e5e79233ba74b84d3565a582dbe8d6d695c774a255fb1a175ba73fd87`.
- [Output](planar_jc48_sep08_octuple_transport.out), 398 bytes, SHA256
  `cf67f2ec2c4e68f1b9cdd2c610634ae2ddbb6a091d349484c1a1dfa97ad7aae8`.
- Semantic digest
  `501257095e3b0f36f4657438ef9b4dd59de322b7d4c88e6978d4ab24f6497560`.

Both independent runs completed successfully with **85 always-active
exact gates** and matched all 398 bytes of the frozen output. Fresh
hashes of the source and all three outputs match the pins above. The
source/output were unchanged throughout the audit. Independent full
section-matrix and symmetric-trace calculations supplement these replays.

No mathematical correction was needed. The stale fixed-zero dependency
label was corrected before acceptance. The complete candidate is accepted
for owner-controlled promotion with the polynomial-only main conclusion
and explicit rational hostiles retained.
