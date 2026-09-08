# The fixed finite-octuple quartic class has no polynomial Jacobian mate

**Status: RESERVED / analytic proof candidate with FINITE-EXACT controls; independent audit pending.**
The point is the specified finite point `(x,z)=(0,0)` of the actual
DG compactification. The main conclusion concerns polynomial
mates of unrestricted degree. Some intermediate classes exclude
rational mates; the whole class does not.

## 1. Statement and complete inherited entry

Use

    W=(P1_x x P1_z) minus S,   S={z=x^2},
    s=z-x^2,  t=1/s,
    x=1/r,  t=-r^2-r^4*b_D,
    omega=dx wedge dt=r^2 dr wedge db_D.

Let `H=N/s^2` and `L=M/s` be global functions in the proved
spaces `L_2,L_1`, and put `F=H^2+L`.

**Theorem.** If the complete binary boundary octic is

    N|S=alpha*x^8,             alpha!=0,

then no `G in C[x,t]` satisfies `J_(x,t)(F,G)=lambda!=0`.
There is no bound on the degree of `G`. No translation to a
different finite point is asserted.

Set the following actual global functions:

    w=x^2*t,     v=x(1+w),     h=x^2(1+w).

The complete `L_1` space has basis `1,t,xt,w,v,h`. Since `h^2`
has numerator `x^4*z^2`, its boundary octic is `x^8`. A numerator
with zero restriction to `S` is divisible by its defining section,
so the entire class in the theorem is

    H=alpha*h^2+K,
    K=k0+kt*t+kxt*x*t+k2*w+k3*v+k4*h,
    L=l0+lt*t+lxt*x*t+l2*w+l3*v+l4*h.                 (1)

This is the actual section-space decomposition, not a truncated
local model. The leading `t` coefficient of `H` is `alpha*x^8`,
and the leading coefficient of `F` is `alpha^2*x^16`.

The closest proved mechanisms are the
[M-avoidance theorem](planar_jc48_sep08_quartic_common_root.md),
the finite first-jet and regularity statements in
[the shared-root theorem](planar_jc48_sep08_shared_roots.md),
the componentwise [pole-degree gate](planar_jc48_sep08_pole_degree.md),
and the complete
[boundary-linear exclusion](planar_jc48_sep08_boundary_linear.md).
The corrected near miss is to confuse residue freedom with a
rational primitive, or rational integrability with a polynomial
mate. The exact hostile in Section 6 distinguishes both boundaries.
The least-used operation is field trace of the actual relative
differential, retaining its source-volume coefficient.

The live concepts are complete global sections, normalized branch
counts, componentwise pole degree, source criticality, and trace
of an exact differential. Their roles are separate: a local
residue obstruction alone does not supply the final polynomial
conclusion.

If `L` is constant, the identity
`J(F,G)=2H*J(H,G)` immediately excludes polynomial `G` because
`H` is a nonconstant polynomial. For the remaining proof we may
assume `L` nonconstant and, seeking a contradiction, that a
polynomial mate exists. It is in particular a rational mate.

## 2. The complete first-jet reduction

The only zero of the boundary octic is the finite point zero;
its value at the distinguished infinity point is nonzero.
The M-avoidance theorem therefore forces `lt=M(0,0)=0`.

If `kt!=0`, then `N_s(0,0)!=0`. The shared-root regularity lemma
makes the relative form regular on every normalized branch there,
regardless of the multiplicity of `M`. There is no other boundary
intersection of a generic compact fibre. Thus its relative form
is holomorphic on each compact component and nonzero on its
intersection with the original source plane. It cannot equal a
rational exact differential. Hence `kt=0`.

For clarity, the next use of the finite first-jet theorem can be
seen directly. With `kt=lt=0`, substitute `s=x*Z` in

    E=N^2+s^3*M-c*s^4=0,      eta=s^2*dx/E_s.

Its degree-four tangent polynomial is

    Z^2[(kxt+k0*Z)^2+lxt*Z+(l0-c)*Z^2].              (2)

If `kxt!=0`, the bracket has two simple nonzero roots for generic
`c`; each yields an actual branch and a nonzero logarithmic
residue of `eta`. If `kxt=0` but `lxt!=0`, (2) has one simple
nonzero root with the same consequence. The simple roots are
paid by the nonzero constant term in the first case and the
nonzero linear coefficient in the second; the original fibre
value supplies the generic highest coefficient. Neither argument
counts the zero roots as normalized branches.

Thus every remaining polynomial-mate candidate has exactly

    H=alpha*h^2+k0+k2*w+k3*v+k4*h,
    L=l0+l2*w+l3*v+l4*h.                             (3)

On the actual boundary `D`, these functions restrict to

    H_D=alpha*b_D^2+k0-k2-k4*b_D,
    L_D=l0-l2-l4*b_D.

In particular `F_D` has degree four with leading coefficient
`alpha^2`. A generic original fibre has points on `D` at which
the tangential derivative of `F_D` is nonzero. On the compact
component containing any such point, `eta` has order exactly
two, so a rational primitive would have local degree three.

## 3. The order-two normal coefficient is excluded by all four branches

Suppose `k2!=0`. Since `L` is nonconstant, let

    n=ord_x(l2*x^2+l3*x^3+l4*x^4) in {2,3,4}.

The first two branches are obtained by `s=x^2*Z`. The complete
leading equation is

    E/x^8=Z^2*T_c(Z)+higher powers of x,
    T_c=k2^2+(2k0*k2+l2)*Z+(k0^2+l0-c)*Z^2.          (4)

Its discriminant has nonzero `c` slope `4k2^2`; hence its two
roots are simple and nonzero for generic `c`. Each is an actual
smooth branch with parameter `x`. Since
`E_s=x^6*Z^2*T_c'(Z)+higher terms` on the branch,

    eta=(1/T_c'(Z)+O(x))*x^(-2)*dx.                  (5)

Each has a genuine double pole and possible primitive pole degree
one. Any nonzero residue already rules out exactness; even when
both residues vanish, their combined primitive budget is two.

The other two Puiseux determinations lie in the cancellation of
the leading terms of `N`. Substitute `s=x^6*Z`. Then

    N/x^8=alpha+k2*Z+O(x).

The implicit function theorem gives an exact analytic centre
`Z=Zc(x)` with `N=0`, `Zc(0)=-alpha/k2!=0`. Near this centre,
`N_s` has order two and nonzero leading coefficient `k2`.
Moreover `M-c*s` has order `n` with nonzero leading coefficient
`l_n`; its `s`-dependent terms have strictly larger order because
`n<=4<6`. Consequently the equation

    N^2=-s^3*(M-c*s)

has two split determinations with

    ord_x N=9+n/2,
    ord_x E_s=11+n/2,
    eta=(nonzero constant+...)*x^(1-n/2)*dx.          (6)

The last derivative order is exact: `2N*N_s` dominates the other
terms of `E_s`, which have order at least `12+n`. Higher normal
jets cannot cancel its nonzero leading coefficient. Equivalently,
after the analytic centre change, the split equation has a
nonzero quadratic coefficient and nonzero leading term of order
`n+2`. Thus for even `n` there are two smooth branches; for odd
`n` the two determinations give one normalized branch `x=tau^2`.

For `n=2`, both branches in (6) are regular. For `n=3`, the
normalized differential is a unit times `d tau`, also regular.
For `n=4`, both smooth branches have simple poles with nonzero
residues, so no rational primitive exists. No other local branch
is missing: `E(0,s)=(k0^2+l0-c)*s^4`, and the two branches in
(4) plus the two determinations in (6) exhaust degree four.

For `n=2,3` the total possible primitive pole degree over the
entire boundary is at most two. All generic compact components
have regular relative form inside `W`; there are no other zeros
of `N|S`. The same component as a generic point of `D` would
need degree at least three by Section 2. Its own pole degree is
no larger than the total budget two, a contradiction. No generic
irreducibility assumption or distribution of components is used.

This proves the stronger intermediate statement: **in (3), if
`k2!=0` and `L` is nonconstant, there is no rational mate.**
Together with the constant-L polynomial argument, every remaining
polynomial-mate candidate must have `k2=0`.

## 4. A linear boundary phase either has a residue or a source critical point

Assume `k2=0` but `l2!=0`. In the actual rational chart `(x,w)`
with `t=w/x^2`, write the original function as

    F(x,w)=f(w)+x*g(w)+O(x^2),
    f=k0^2+l0+l2*w,
    g=(1+w)*(2k0*k3+l3).

The volume form is `x^(-2)*dx wedge dw`, so along the original
generic fibre

    eta=-x^(-2)*dx/F_w.

The unique simple finite root of `f(w)=c` lifts to an actual
branch. Expanding its denominator gives residue

    (g'*f'-f''*g)/(f')^3=(2k0*k3+l3)/l2^2.           (7)

A rational mate forces this to vanish. But on the actual source
line `E0={x=0}` in `(x,t)`,

    F_x=2k0*k3+l3,      F_t=0.

Condition (7) therefore makes that whole source line critical,
contradicting the existence of a polynomial constant-Jacobian
mate. This step is a polynomial-mate exclusion; it does not infer
that source criticality rules out all rational mates. It leaves

    k2=l2=0.                                          (8)

## 5. The remaining nonlinear boundary case has a nonzero trace residue

Assume (8) and `k3!=0`. The functions `h,v` are birational
coordinates on the original source field:

    x=h/v,       t=(v^2-h)*v^2/h^3,
    J_(x,t)(h,v)=h^3/v^2.

Put

    epsilon=l3/(2k3),
    y=H+epsilon,
    A(h)=alpha*h^2+k4*h+k0+epsilon,
    R(h)=epsilon^2-2epsilon*A(h)+l4*h+l0.

Then `v=(y-A)/k3` and, with the original fibre value unchanged,

    F=y^2+R(h),
    omega=(y-A)^2/(k3^3*h^3) * dh wedge dy,
    eta=-(y-A)^2/(2k3^3*h^3*y) * dh.                 (9)

These are exact field and differential identities, not only
equations for an abstract curve. Let `K=C(c)`. The polynomial
`y^2-(c-R(h))` is irreducible over `K(h)`: the prime `c-R(h)`
in `C[c,h]` has odd valuation and cannot be a square in its
fraction field. This argument includes constant `R`; it makes
no assertion of geometric integrality after algebraically
closing `K` in that exceptional case.

Use the degree-two field trace for the involution `y -> -y`.
Taking the trace of the literal form (9) gives

    Tr(eta)=2*A(h)/(k3^3*h^3) * dh.                   (10)

Its residue at `h=0` is `2*alpha/k3^3`, which is nonzero.
In characteristic zero, trace commutes with the relative
derivation fixing `K`. If `eta=dG` in the original generic
field, then (10) would be `d Tr(G)` in `K(h)`. The derivative
of a rational function of `h` has zero residue at every place,
contradicting (10). This proves **no rational mate** for this
entire subfamily, including every degeneracy of `R`.

Finally suppose `k3=0` as well as (8). Then

    F=P(h)+l3*v,       P in C[h].

In the actual boundary chart `h=-b_D`, `v=-r*b_D`, so
`F=P(-b_D)-l3*r*b_D` has degree at most one in `r`.
The proved boundary-linear theorem applies to its exact global
hypothesis and excludes every polynomial mate of unrestricted
degree. If `l3=0`, the source criticality on `x=0` is also an
immediate direct control. These cases exhaust (1), proving the
theorem.

## 6. Scope, hostile controls, and the information retained

The whole theorem cannot be strengthened to rational mates.
Inside its exact finite-octuple class,

    F=h^4+h,       G=1/[3*x^3*(4*h^3+1)]

has `J(F,G)=1`. The proposed mate has poles. There is also the
source-critical-free boundary-linear example `F=h^4-v`, with
an actual rational mate but no polynomial one. Thus the proof
does not reduce the whole class to source criticality.

The earlier residue tuning for `k2!=0` and the disjoint special
fibre calculation are superseded as a proof route by Section 3:
the complete normalized branch budget excludes that entire
nonconstant-L case without tuning its residues. The trace step
then handles a different surviving coefficient, `k3`, instead
of repeating a coefficient scan.

The source of the last connection is the full actual function
field with its source volume form. The target is `K(h)` under
trace. Exactness is preserved, but generic geometric component
information may be lost; the proof does not need it. Retaining
the coefficient of the volume form supplies the nonzero residue
that a genus label, an abstract square equation, or an unweighted
involution would discard. Other octic multiplicity partitions
and other finite locations are not settled by this note.

## 7. Exact reproduction and verification boundary

```sh
python3 -B 04-computation/planar_jc48_sep08_finite_eight.py
python3 -B -O 04-computation/planar_jc48_sep08_finite_eight.py
```

The standalone source imports no inherited mathematical
implementation. Its exact universe is the complete symbolic
entry (1), both first-jet tangent cases, every normalized
`n=2,3,4` cancellation case, the low quadratic face, the full
derivative remainder, the source-criticality calculation, both
birational field inverses, the original volume form, the trace
and its constant-R boundary, the literal boundary-linear input,
and the two inside-class rational-mate hostiles. The source
checks the actual high-order suppliers before comparing orders;
it does not infer their completeness from named samples.

Normal and optimized runs agree byte for byte with the frozen
output: **92 always-active exact gates, 444 bytes**.

- [Exact source](../../04-computation/planar_jc48_sep08_finite_eight.py):
  9,242 bytes, SHA256
  `f9e86132d04caacd894651e66cbc563dbc56ec283adf4f9412334624cae9b02d`.
- [Frozen output](planar_jc48_sep08_finite_eight.out): SHA256
  `2fdb2418487c87b217877559e22167aa42f4be52e217e5d0bec2ebf2cc455bb4`.
- Semantic digest:
  `a2ba0a0d286092c9b8a153c201275c19e9ffeced53abc210a602a31681b3547f`.

The source and output are frozen. The
[independent audit](planar_jc48_sep08_finite_eight_audit.md) accepts
the complete analytic proof, full source and both exact replays,
with a separate22-check section-kernel and trace reconstruction.
No finite computation substitutes for the all-case analytic proof.
