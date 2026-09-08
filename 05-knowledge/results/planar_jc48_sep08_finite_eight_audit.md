# Independent audit of the fixed finite-octuple quartic exclusion

**Status: independent analytic/source audit PASS; normal, optimized and
frozen output agree at 92 always-active exact gates.** This audit accepts
[the finite-eight proof](planar_jc48_sep08_finite_eight.md) for its specified
point of the fixed DG compactification. Its main conclusion is absence of
polynomial Jacobian mates of unrestricted degree. It does not assert
absence of rational mates throughout that class, or transport to other
finite points.

## 1. Complete entry and proved dependencies

Write `s=z-x^2`, `t=1/s`, and use the full section spaces
`H=N/s^2 in L_2`, `L=M/s in L_1`. The assumptions are the exact binary
boundary section `N|S=alpha*x^8`, with `alpha!=0`, and `F=H^2+L`.
In particular, the value of the boundary octic at infinity is nonzero:
its only zero is the fixed finite point zero. This is not just a local
multiplicity assertion which could leave other boundary zeros unpaid.

I independently formed the restriction matrix on all fifteen monomials
`x^i z^j`, `0<=i<=4`, `0<=j<=2`. Substitution `z=x^2` has rank nine;
its kernel is the six-dimensional span of `(z-x^2)*x^i z^j`,
`0<=i<=2`, `0<=j<=1`. The particular numerator `x^4*z^2` restricts
to `x^8`. Division by `s^2` proves the complete affine description

    H=alpha*h^2+K,       K in L_1,
    w=x^2*t, v=x(1+w), h=x^2(1+w).

The six actual numerators for the basis `1,t,xt,w,v,h` are
`s,1,x,x^2,xz,x^2z`; their coefficient matrix in the full `(2,1)` box
has determinant of absolute value one. Thus both six-parameter lower
terms in the primary are complete global sections. No mate-degree
restriction or truncated jet enters this reduction. The source's exact
substitution checks both whole charts; the independent rank argument
pays completeness rather than just verifying its displayed family.

The following dependencies were read in their current proved versions:

- [Quartic common-root necessity](planar_jc48_sep08_quartic_common_root.md),
  a rational-mate obstruction when the two complete boundary sections
  have no common zero.
- [Shared-root local regularity and finite first jets](planar_jc48_sep08_shared_roots.md),
  particularly the normal-unit case and the actual nonzero tangent roots.
- [The componentwise pole-degree gate](planar_jc48_sep08_pole_degree.md),
  including rational mates and possibly reducible generic fibres.
- [The complete boundary-linear polynomial exclusion](planar_jc48_sep08_boundary_linear.md),
  whose conclusion concerns arbitrary polynomial mates in `C[x,t]`,
  rather than only mates regular on the enlarged surface.

None is a RESERVED dependency. The new finite-eight argument combines
these existing mechanisms with its complete four-branch calculation and
the trace of the literal source differential. It does not claim the
classical exactness or degree principles themselves as new.

## 2. Initial reductions and the actual boundary degree

If `L` is constant, `J(F,G)=2H J(H,G)` cannot be a nonzero constant
for polynomial `G`: `H` is nonconstant with source degree two in `t`.
This argument does not exclude rational `G` and is used only for the
polynomial conclusion.

Otherwise a proposed polynomial mate is also a rational mate. Common-root
necessity forces `lt=0` because there is just one possible boundary
intersection. If `kt!=0`, normal-unit regularity makes the relative form
holomorphic on every normalized local branch at that point. There is no
other point on `S` where a generic fibre can meet the boundary. The form
is regular inside `W` and nonzero on the open source part of a generic
component. A holomorphic exact meromorphic differential on a compact
component must vanish; hence this case has no rational mate.

After `kt=lt=0`, the complete degree-four tangent polynomial under
`s=xZ` is

    P_c(Z)=Z^2[(kxt+k0 Z)^2+lxt Z+(l0-c)Z^2].

For `kxt!=0`, the remaining quadratic has nonzero constant term and
discriminant with nonzero `c` coefficient `4kxt^2`; its two roots are
simple and nonzero generically. For `kxt=0`, `lxt!=0`, the remaining
nonzero root is simple. On any such actual branch,

    eta=(Z0^2/P_c'(Z0)+O(x))*dx/x.

It has a nonzero residue. The argument never calls its zero roots
normalized branches. It correctly forces `kxt=lxt=0`.

The surviving exact functions have

    H_D=alpha*b_D^2+k0-k2-k4*b_D,
    L_D=l0-l2-l4*b_D.

Consequently `F_D` has degree four and leading coefficient `alpha^2`
for every remaining choice of coefficients. In particular there is no
unaddressed constant-D subcase. At a generic simple root `F_D(b0)=c`,
the whole-chart volume `omega=r^2 dr wedge db_D` gives a nonzero
relative differential of order exactly two on the smooth component
through `(0,b0)`. Any rational primitive there has local degree three.
The comparison below uses that same component; it does not add degrees
of points lying in different fibres of the putative primitive.

## 3. All four determinations when k2 is nonzero

Assume `k2!=0` and `L` nonconstant. Since its two lower coefficients
have already vanished, the actual order

    n=ord_x(l2*x^2+l3*x^3+l4*x^4)

belongs to `{2,3,4}`. There is no omitted infinite-order case: that is
exactly the previously separated constant-L case.

The low substitution `s=x^2 Z` has full face

    E/x^8=Z^2 T_c(Z)+O(x),
    T_c=k2^2+(2k0*k2+l2)Z+(k0^2+l0-c)Z^2.

Its two nonzero roots are generically simple, again by the nonzero
`c` coefficient `4k2^2` in the discriminant. They give two genuine
smooth branches parameterized by `x`. The polar derivative on either
branch has leading term `x^6 Z0^2 T_c'(Z0)`, and therefore

    eta=(1/T_c'(Z0)+O(x))*x^(-2) dx.

Each pole is exactly double. Even if its residue vanishes, it supplies
only one unit of pole degree to a putative primitive.

For the remaining determinations, put `s=x^6 Z`. The exact normal
coefficient is

    N/x^8=alpha+k2 Z+O(x).

Its derivative in `Z` is a unit at the nonzero root `-alpha/k2`.
Thus the actual equation `N=0` has an analytic centre `Zc(x)`;
using that centre does not discard higher normal terms. On the nearby
fibre, `M-cs` has order precisely `n`, because `n<=4<6` and its
leading coefficient is the declared nonzero `l_n`. It follows from
`N^2=-s^3(M-cs)` that

    ord_x N=9+n/2,
    ord_x N_s=2,
    ord_x E_s=11+n/2,
    eta=(unit)*x^(1-n/2) dx.

The last derivative order is strict: all terms other than `2N N_s`
have order at least `12+n`, and cannot cancel it. In the centred
`Z` coordinate the square split has nonzero leading term of order
`n+2`. For even `n` it gives two smooth branches. For odd `n` it
gives one normalized branch with `x=tau^2`, counted twice in the
four-sheet Weierstrass projection. The orders of the actual normalized
differential are consequently

| n | Normalized branches at the high centre | Order of eta |
|---|---|---|
| 2 | two smooth branches | 0 on each |
| 3 | one branch, ramification two over x | 0 |
| 4 | two smooth branches | -1 on each |

For `n=4` both residues are nonzero, since all displayed leading
coefficients are units. For `n=2,3` these branches are regular.
Finally `E(0,s)=(k0^2+l0-c)s^4` has Weierstrass degree exactly four
at a generic value. The low two branches and high two determinations
exhaust it. This analytic count supplies completeness; a finite list
of order checks alone would not do so.

Thus either a nonzero residue already excludes exactness, or all
possible primitive poles over the entire boundary have total degree at
most two. A rational primitive has no pole inside `W`, because its
derivative is regular there. Restricting to the component meeting the
actual generic point of `D` would give

    3 <= degree(G on that component) <= 2.

This proves the stated rational-mate exclusion without assuming generic
irreducibility, geometric integrality, or a distribution of the local
branches among components.

## 4. The residue which must not be mistaken for rational exclusion

When `k2=0`, `l2!=0`, use the birational chart `t=w/x^2`.
The exact initial expansion of the original function is

    F=f(w)+x g(w)+O(x^2),
    f=k0^2+l0+l2*w,
    g=(1+w)(2k0*k3+l3).

The volume is `x^(-2) dx wedge dw`. On the actual branch through the
simple finite root `f(w)=c`, expansion of `-x^(-2)dx/F_w` gives
residue `(2k0*k3+l3)/l2^2`. Rational exactness requires its numerator
to vanish. The same numerator is `F_x` on the original source line
`x=0`, where `F_t=0` identically. This excludes a polynomial
constant-Jacobian mate. The proof explicitly does not use affine
criticality to exclude rational mates; that stronger inference would
be false. It correctly reduces only the polynomial-mate argument to
`k2=l2=0`.

## 5. Trace in the actual field, including the constant-R case

For `k3!=0`, the exact field identities are

    x=h/v,  t=(v^2-h)*v^2/h^3,
    J_(x,t)(h,v)=h^3/v^2.

The two inverse formulas are birational identities, not just maps to an
abstract model. With the primary's definitions of `epsilon,A,R` and
`y=H+epsilon`, the original function and source volume become

    F=y^2+R(h),
    omega=(y-A(h))^2/(k3^3*h^3) dh wedge dy,
    eta=-(y-A(h))^2/(2*k3^3*h^3*y) dh.

For `K=C(c)`, the element `c-R(h)` has valuation one at its prime
in `C[c,h]`; hence it is not a square in `K(h)`. The original generic
field is the degree-two field `K(h,y)` with `y^2=c-R(h)`. This remains
true if `R` is constant. In that case the relative constants may enlarge,
and the geometric generic fibre may split, but no step requires their
absence.

I independently derived the trace obstruction by parity. Every element
of this field is uniquely `G=f(h)+y*g(h)`, with `f,g in K(h)`.
The derivation fixing `K` has `y'=-R'(h)/(2y)` and commutes with the
involution `y -> -y`. Thus the even part of `dG` is `df`. The even
part of the literal differential above is

    A(h)/(k3^3*h^3) dh,

whose residue at zero is `alpha/k3^3!=0`. Equivalently the full trace
has residue `2alpha/k3^3`. A rational derivative in `K(h)` has zero
residue, a contradiction. This elementary decomposition also checks
that no geometric-integrality or constants-field hypothesis is hidden
inside the use of field trace. A nonzero Jacobian scalar just multiplies
the same nonzero residue.

If `k3=0` as well, the entire function is `F=P(h)+l3*v`. Its exact
actual boundary expression is `P(-b_D)-l3*r*b_D`, of degree at most
one in `r`. This is precisely the proved boundary-linear theorem's
input; its conclusion excludes arbitrary polynomial mates. The
constant-L and all other earlier branches together exhaust the
complete coefficient space in Section 1.

## 6. Hostiles, scope, and reproducibility

Both inside-class rational mates are retained and directly checked:

    F=h^4+h,
    G=1/[3*x^3*(4*h^3+1)],

and

    F=h^4-v,
    G=-1/(2*x^2)+2*h^3/x-(4/3)*h^6.

Each has Jacobian one. The latter first component is critical-point-free
on the original source, as checked in both charts. These controls rule
out two incorrect strengthenings: the full class does not have no
rational mates, and the proof does not reduce all of it to a source
critical point. Intermediate rational exclusions are confined to their
stated coefficient strata. Other octic partitions, other finite
locations and general JC(2) are outside this proof.

I read the full standalone source and independently ran

    python3 -B 04-computation/planar_jc48_sep08_finite_eight.py
    python3 -B -O 04-computation/planar_jc48_sep08_finite_eight.py

Both completed successfully and matched the frozen output byte for byte:
**92 always-active gates, 444 bytes**. No inherited mathematical
implementation is imported. A separate temporary reconstruction made
22 additional exact checks using the full section restriction matrix,
a unimodular actual carrier basis, parity of an arbitrary quadratic-field
primitive, and the normalized high-branch orders. It did not import the
producer. Those checks supplement the analytic proof; they are not a
claim that a finite census proves its quantifiers.

Frozen reproduction pins:

- [Source](../../04-computation/planar_jc48_sep08_finite_eight.py),
  9,242 bytes, SHA256
  `f9e86132d04caacd894651e66cbc563dbc56ec283adf4f9412334624cae9b02d`.
- [Output](planar_jc48_sep08_finite_eight.out), 444 bytes, SHA256
  `2fdb2418487c87b217877559e22167aa42f4be52e217e5d0bec2ebf2cc455bb4`.
- Semantic digest
  `a2ba0a0d286092c9b8a153c201275c19e9ffeced53abc210a602a31681b3547f`.

No source or proof correction was needed. The candidate is accepted for
status promotion by its owner with the polynomial-only main scope above.
