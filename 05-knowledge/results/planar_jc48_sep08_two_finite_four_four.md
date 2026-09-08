# Two distinct finite fourfold boundary roots exclude even a rational mate

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a rational-mate exclusion in the fixed DG global quartic section
class. The finite positions and their distinctness matter. The proof does
not identify arbitrary boundary changes with source-volume-preserving
surface automorphisms.

## 1. Actual family and complete global sections

Let

    W=(P1_x x P1_z) minus S,   S={z=x^2},
    t=1/(z-x^2),              omega=dx wedge dt.

Let `H in L2`, `L in L1` be global, `F=H^2+L`. Suppose its complete
boundary octic is

    N|S=alpha*(x-p)^4*(x-q)^4,   alpha!=0,  p!=q in C. (1)

**Theorem.** There is no `G in C(x,t)` with
`J(F,G)=lambda!=0`. In particular there is no polynomial mate of any
degree, and no assumption that `L` is nonconstant is needed.

Put `S1=p+q`, `P1=pq`, and `r=(x-p)(x-q)`. The actual global function

    Q=r^2*t+x^2-2S1*x                               (2)

has boundary numerator `r^2`. In the standard global `L1` basis it is

    Q=h-2S1*v+(S1^2+2P1)w-2S1*P1*xt+P1^2*t,
    w=x^2t,   v=x+x^3t,   h=x^2+x^4t.

The full rank-nine octic restriction map on `H^0(O(4,2))` has a
six-dimensional kernel, exactly `S` times `H^0(O(2,1))`. Thus every
entry with (1), without omitting any induced coefficient, is

    H=alpha*Q^2+K,   K in L1,   L in L1.             (3)

The `L1` restriction to the graph is the whole space of polynomials of
degree at most four; its kernel consists of constants. The actual
section-space arguments are those of the
[DG filtration](planar_jc48_sep08_dg_quadratic.md) and
[shifted octuple transport](planar_jc48_sep08_octuple_transport.md).
No translation of `W` is used in (2) or (3).

The closest local mechanisms are the `M`-unit classification in
[common-root necessity](planar_jc48_sep08_quartic_common_root.md), the
normal-unit and finite first-jet gates in
[shared roots](planar_jc48_sep08_shared_roots.md), and the exact moving-root
residue in the [4+4 proof](planar_jc48_sep08_four_four.md).
The new feature is a leading first-row term of degree four that none of
the complete global corrections can remove.

The live concepts are full global section restrictions, a regular versus
singular shared root, moving-root residues, the other root's finite
position, and the distinction between polynomial and rational exactness.
The source-to-target map restricts the original volume-relative form to
the compact generic fibre and retains each labelled residue separately.
It loses source polynomiality, so the direct residue contradiction is
stronger than a mere source critical-point obstruction here.

## 2. A verified composition reduction at two singular shared roots

Write `K=A_K(x)t+B_K(x)` and `L=A_L(x)t+B_L(x)`. The polynomials
`A_K,A_L` have degree at most four. Let `B_Q=x^2-2S1x`. In graph
coordinates `s=z-x^2`, the complete numerator is

    N=alpha*r^4+s(2alpha*r^2 B_Q+A_K)
      +s^2(alpha*B_Q^2+B_K),
    M=A_L+s B_L.                                    (4)

At either root of `r`, the term `r^2` vanishes to order two. Therefore
`N_s` and its first tangential derivative there equal `A_K` and `A_K'`.
The corresponding two jets of `M|S` are `A_L,A_L'`.

Suppose both boundary points are singular shared roots, meaning `M=0`
and `N_s=0` there. A proposed rational mate must satisfy the proved finite
first-jet gate, so `A_K,A_L` and their first derivatives vanish at both
points. As `p!=q` and their degree is at most four,

    A_K=beta*r^2,       A_L=delta*r^2.

Since the `L1` restriction kernel is constants, the entire functions are

    K=beta*Q+gamma,     L=delta*Q+epsilon.

Thus `F=P4(Q)` for a degree-four polynomial with leading coefficient
`alpha^2`. This verifies the proposed complete-section composition
reduction. It already excludes polynomial mates by the nonconstant
factor `P4'(Q)` in their source Jacobian. The argument below is stronger:
it excludes rational mates as well, and only needs one singular shared
root. Constant `L` is retained throughout.

## 3. The full boundary exhaustion for a hypothetical rational mate

At either boundary point the order of `N|S` is four. If `M` is a unit,
the proved `M`-unit local classification makes the relative form regular:
the only exceptional balance would require `4=3j`, impossible for an
integer normal order. If `N_s` is a unit, the normal-unit local lemma
makes it regular for arbitrary `M`.

The only remaining case at that point is

    M=0,   N_s=0.

If a first tangential derivative of either is nonzero, the finite
first-jet gate gives a simple pole with nonzero residue and rules out
a rational mate immediately. Hence a hypothetical mate must either
have a locally regular unit case, or have

    A_K(p)=A_K'(p)=A_L(p)=A_L'(p)=0                  (5)

at a singular shared root; the same alternative holds at `q`.

If both points are regular unit cases, all boundary branches of the
compact generic fibre have regular relative form. There are no other
boundary roots: (1) is a degree-eight polynomial with nonzero leading
coefficient, so infinity is not a boundary zero. Inside `W` the original
relative form is regular on the smooth generic fibre, including `D`,
where the original volume has its zero of order two. No generic component
is `D` or `S`; every component meets the source plane where the volume
is nonzero. Thus each compact generic component has a holomorphic,
nonzero relative form. It cannot be a rational derivative, since a
rational primitive with no differential poles is a constant on a compact
curve.

Consequently a proposed rational mate would require at least one point
with (5). We now show that this single point is impossible as well.

## 4. At one singular shared root the residue has an unavoidable leading term

Choose that root as `p`, and write

    u=x-p,   d=p-q!=0,   w=u^2t,
    v_p=u(1+w)-2p,   h_p=u^2(1+w)-2pu.

These are the actual corrected global functions for the original surface.
Condition (5) is equivalent, in the full six-dimensional global basis,
to removing only the coefficients of `t` and `ut`. Thus the whole
remaining correction space is

    K=k0+k2*w+k3*v_p+k4*h_p,
    L=l0+l2*w+l3*v_p+l4*h_p.                         (6)

It is not a finite-degree ansatz for unknown corrections; it is the full
section space with the exact two boundary jets prescribed.

The literal global function (2) becomes

    Q=Q0(w)+u*Q1(w)+u^2(w+1),
    Q0=d^2*w-p^2-2pq,
    Q1=2d*w-2q.                                    (7)

The full expansion of `H` begins `A(w)+u B(w)+O(u^2)`, where

    A=alpha*Q0^2+k2*w+k0-2p*k3,
    B=2alpha*Q0*Q1+k3(1+w)-2p*k4.

Therefore the original `F` has complete first two rows

    F=f(w)+u*g(w)+O(u^2),
    f=A^2+l2*w+l0-2p*l3,
    g=2A*B+l3(1+w)-2p*l4.                           (8)

All induced first-row corrections are included. The leading terms are

    [w^2]A=alpha*d^4,    [w^2]B=4alpha*d^3,
    [w^3]f'=4alpha^2*d^8,
    [w^4]g=8alpha^2*d^7.                            (9)

For a generic fibre value `c`, `f(w)=c` has four simple nonzero roots.
The full local Weierstrass degree is four because at `u=0` its graph
numerator is `(H(0)^2+L(0)-c)s^4`. These four branches are the actual
ones under `s=u^2/w`, not virtual moment roots.

At a generic root `w0`, the branch has
`w=w0-g(w0)u/f'(w0)+O(u^2)`. The original volume is
`u^-2 du wedge dw`; hence the restricted differential and its residue
are

    eta=-du/(u^2 F_w),
    Res eta=(g'f'-f''g)(w0)/f'(w0)^3.                (10)

A rational primitive would force the numerator to vanish at every such
root for generic `c`, equivalently the polynomial identity
`g'f'-f''g=0`. One may justify this either by varying `c`, or by
noting that a generic `w0` with `f(w0)=c` is transcendental over the
constant field. This would make `g=C f'` for a constant `C`.

But (9) makes that impossible: `g` has degree four and `f'` degree
three. More explicitly the complete residue numerator has

    degree_w(g'f'-f''g)=6,
    [w^6](g'f'-f''g)=32alpha^4*d^15!=0.              (11)

For generic `c` none of its four roots can be a zero of this fixed
nonzero polynomial, after discarding finitely many fibre values. Thus
there is an actual nonzero residue at this boundary point. This rules
out the last possible local alternative, proving the theorem.

## 5. Scope and exact hostile boundaries

The position of the second root supplies the uncancellable term in (9).
If the other root is infinity, the first-row degree may instead be three;
the all-finite/infinity [4+4 transport theorem](planar_jc48_sep08_four_four_transport.md)
uses additional critical-point arguments and only excludes polynomial mates.
Indeed

    w=x^2t,  H=w^2,  L=w,  F=w^4+w,
    G=1/[x(4w^3+1)]

has rational Jacobian one and boundary partition `0^4 infinity^4`.
Thus a rational-mate claim cannot be transported by forgetting the
actual boundary positions and volume.

Distinctness also matters to the rational statement. When the two roots
coalesce, `h=x^2+x^4t`, `H=h^2`, `L=h`, `F=h^4+h` has the rational
mate `G=1/[3x^3(4h^3+1)]`. Its octic is `x^8`, and the factor `d^15`
in (11) has disappeared.

There is a useful check on the initial composition subfamily. On a
compact generic fibre of `Q` itself, the relative form is
`-dx/[(x-p)^2(x-q)^2]`. Its residues at the two finite points are
`2/(p-q)^3` and `-2/(p-q)^3`. Their sum is zero, as required on a
compact curve, but each is nonzero. Summing residues would therefore
lose the exactness obstruction. The present proof keeps each residue.

**All-location corollary.** The independently audited finite/infinity
supplier is proved. Combining the two results gives polynomial-mate
exclusion for every boundary octic
having exactly two distinct fourfold roots anywhere on `P1`. That
supplier is not needed for the present two-finite-root proof. The stronger
rational conclusion is only asserted when both roots are finite. Other
partitions, general global quartics, and JC(2) remain outside the claim.

## 6. Exact controls and reproduction

The new standalone [source](../../04-computation/planar_jc48_sep08_two_finite_four_four.py)
and [frozen output](planar_jc48_sep08_two_finite_four_four.out) pass **60
always-active exact gates**. Normal and optimized Python outputs are
byte-identical. No inherited mathematical implementation is imported.

The exact universe is every pair of distinct finite complex points,
every nonzero `alpha`, and the full global `K,L` correction spaces.
Symbolic all-parameter identities pay the global `Q`, section boxes,
rank-nine and rank-five restrictions, two-point confluent determinant,
complete normal jets, original first rows and degree-six residue.
The actual tangent quartic is checked against `Z^4(f(1/Z)-c)`, so the
root equation and the residue computation refer to the same fibre.
The local unit/first-jet and compact exactness arguments remain analytic;
no finite root sample is used as their replacement.

Named controls include several distinct finite-point configurations,
constant `L`, both rational-mate scope hostiles, and the cancellation of
opposite individual residues. A further original-fibre check uses
`p=1,q=2,H=Q^2,L=0,F=Q^4`, fibre `F=1`, and boundary coordinate
`w0=6`, where `Q=1`. The computed residue is `-1/2`, independently
equal to the `Q`-fibre residue `-2` divided by `4Q^3=4`.

Reproduce from the worktree root:

```bash
python3 04-computation/planar_jc48_sep08_two_finite_four_four.py
python3 -O 04-computation/planar_jc48_sep08_two_finite_four_four.py
```

Frozen SHA256 pins:

- Source, 9,639 bytes:
  `2a901ed396b72a7b7eb14f6dc34f30d53e82079c454f22027b86528dca69b630`.
- Output, 423 bytes:
  `67510ca460a435fefef40ab9da962924ab7c9956f078ed88ee60be88c6e23f05`.
- Semantic digest:
  `3e1231865784dc4ad39f4ad9ad7cf13e17261b431c1d6df7e34c178c519c5d17`.

The [independent complete analytic/source audit](planar_jc48_sep08_two_finite_four_four_audit.md)
accepts all60 normal/optimized gates, complete sections and local alternatives,
the moving-root residue, constant L and split composition fibres. A separate
original-fibre computation recovers the named residue -1/2. No unproved
dependency remains in the all-location polynomial corollary.
