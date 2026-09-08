# Independent audit of the all-finite-point 4+4 transport

**Status: INDEPENDENT ANALYTIC / SOURCE / FINITE-EXACT AUDIT PASS.**
The [shifted 4+4 proof](planar_jc48_sep08_four_four_transport.md), its full standalone source, and independent normal/optimized replays are accepted without correction. The coordinating agent owns the primary's status promotion.

The new proof excludes a polynomial Jacobian mate of any degree for the complete global class with boundary octic α(x−p)⁴, αp≠0. The [distinguished-point theorem](planar_jc48_sep08_four_four.md) is currently **PROVED / INDEPENDENTLY AUDITED** and supplies p=0. Together these cover a finite point of multiplicity four and the distinguished infinity point of multiplicity four. They do not normalize an arbitrary pair of finite points, identify translation with a surface automorphism, or exclude rational mates in general.

## 1. Complete actual section space

Keep the original surface and source

    W=(P1_x×P1_z)\{z=x²},  t=1/(z−x²),
    u=x−p,  s=z−(u+p)².

The corrected functions are

    w=u²t,  v=u(1+w)−2p,
    h=u²(1+w)−2pu=uv.

Independently clearing their denominators gives the following numerators for the six L1 functions (1,t,ut,w,v,h):

    s,
    1,
    u,
    u²,
    (u−2p)z+3p²u+2p³,
    (u²−2pu)z+3p²u²+2p³u.

Each has degree at most two in u and one in z. Their coefficient determinant in that six-dimensional section box is one for every p. This verifies actual globality and completeness, including the two correction terms that a naive translation would lose.

For L2, the section box has dimension fifteen. Restriction to z=(u+p)² has rank nine: the source's ordered monomials

    1,u,u²,u³,u⁴,u³z,u⁴z,u³z²,u⁴z²

restrict to polynomials of successive degrees zero through eight with leading coefficient one. Thus the rank-nine minor is identically one, not merely nonzero at a sample p. The kernel has dimension six; multiplication of the full L1 space by s gives precisely that kernel. The base numerator αu⁴ lies in the global box. Hence every entry is exactly

    H=αw²+k0+kt t+kxt ut+k2w+k3v+k4h,
    L=l0+lt t+lxt ut+l2w+l3v+l4h.

In particular the section-space argument covers all p, while the later new critical-polynomial argument explicitly assumes p≠0. No leading coefficient of a possible mate is bounded or specified.

The exact cleared normal expressions are

    N=αu⁴+s(kt+kxt u+k2u²+k3u³+k4u⁴)
       +s²[k0−2pk3+(k3−2pk4)u+k4u²],
    M=lt+lxt u+l2u²+l3u³+l4u⁴
       +s[l0−2pl3+(l3−2pl4)u+l4u²].

These retain every shifted term. As a binary octic N|S has its only roots at u=0 and the original infinity point, both of multiplicity four. The original affine translation u=x−p has determinant one on the source plane; it is not used as an automorphism of W.

If L is constant, the nonconstant factor 2H in J(H²+L,G) immediately excludes a polynomial mate. The nonzero boundary restriction ensures that H is nonconstant. The remaining proof may therefore assume nonconstant L and normalize a hypothetical constant Jacobian to one.

## 2. Infinity with the full nonconstant unit

The actual surface inversion is (x,z)→(1/x,1/z). In its source chart

    r=1/x,  T=−x²−x⁴t,
    dx∧dt=r²dr∧dT.

The source verifies the full substituted H and L are polynomials in r,T, not just truncated series. The boundary numerator becomes

    αr⁴(1−pr)⁴.

The factor at r=0 is a unit. Its other zero at r=1/p is the original finite boundary point and is handled separately. The audit does not discard that second root or use the unweighted volume at r=0.

Here is an independent check of every local branch class. Write

    N=r⁴a(r)+sB(r)+s²C(r,s),  a(0)=α≠0,
    M=D(r)+sE(r,s),
    ℰ=N²+s³M−cs⁴,
    η=r²s²dr/ℰ_s.

Higher terms in a(r) have strictly larger weight than its leading term in every face used below. The coefficient of the generic fibre parameter c also controls the exact cancellation-centre orders, rather than only those of an approximate centre.

* If B(0)≠0, the exact analytic centre N(r,ψ(r))=0 has order four. Generically ℓ=ord(M−cψ)≤4 because ψ is independent of c with a nonzero order-four term. The unweighted relative form has exponent 2−ℓ/2 in r, which is nonnegative; the actual weight adds two. There are exactly two local determinations.
* If B(0)=0 and D(0)≠0, this is the order-four M-unit case. The balanced equation 4=3j has no integer normal order j. For j=1, one low determination is regular, and the two cancellation determinations have unweighted exponent 1/2 in r. On their actual ramified normalization r=τ² this is order two, before the additional weight. For j≥2 the dominant cubic has r=τ³, s of order eight; the unweighted form has order two. These account for the complete Weierstrass degree three.
* After both units vanish, ℰ(0,s) has generic exact degree four. If B has order one, the low face is Z² times a quadratic at s=rZ. Its two nonzero roots are simple generically because its constant term is B1²≠0 and its c-dependent discriminant cannot vanish identically. Their unweighted forms are logarithmic, so their weighted orders are one. The other two determinations are centred at s=−αr³/B1+…. There N_s has order one and 1≤ℓ=ord(M−cs)≤3 generically. Thus the weighted form has exponent (5−ℓ)/2 in r. If ℓ is even the ramified parameter change supplies the corresponding positive integral normalized order. The exact centre is independent of c, so its nonzero cubic leading term pays the upper bound even with all higher-unit coefficients present.
* If ord B≥2 and ord D=1, one low branch is weighted-regular. The other three determinations have r=τ³, s=τ⁷ times a unit and leading equation α²+D1Z³=0. The order of ℰ_s is seventeen and that of r²s²dr is twenty-two; the normalized differential has order five. Three determinations need not mean three different normalized points, and no such interpretation is used.
* When both first jets vanish, s=r²Z has the complete face

      P(Z)=(α+B2Z+C0Z²)²+D2Z³+(E0−c)Z⁴.

  Its constant term is α². Writing P=P0−cZ⁴, any repeated nonzero root satisfies ZP0'−4P0=0; that polynomial has nonzero constant term −4α². Only finitely many c can therefore be exceptional. Its generic leading coefficient is nonzero as well. The four roots are simple and nonzero, ℰ_s has order six and the weighted numerator has order six, so all four forms are regular.

This pays all coefficient degenerations at the original infinity point. The generic-fibre qualifications mean excluding finitely many c for each fixed coefficient tuple; they do not assume the coefficient tuple itself is generic. The new source independently expands the enlarged higher-unit model and checks the faces, derivative orders and c-coefficient suppliers directly.

## 3. Finite first jets and all residue coefficients

At u=0, N|S is exactly αu⁴. The same unweighted normal-unit and M-unit cases are regular; infinity has just been paid for all coefficients. On every compact generic component the relative form would then be holomorphic and nonzero. Each such component meets the original source: a fixed deleted divisor or the other boundary chart cannot be a component of a generic level. The generic fibre is smooth there in characteristic zero, and the source volume is nonzero. A rational primitive of a holomorphic differential on a compact curve has no poles and is constant, a contradiction. Thus lt=kt=0.

If kxt≠0, either simple nonzero low root of the quadratic face gives an unweighted simple pole with nonzero residue. If kxt=0 and lxt≠0, the one nonzero low root does the same. The higher shifted coefficients do not enter these leading terms. Consequently a proposed mate forces

    lt=kt=kxt=lxt=0.

In the actual chart (u,w=u²t), the volume is u⁻²du∧dw and

    F=f(w)+ug(w)+O(u²),
    A=αw²+k2w+k0−2pk3,
    f=A²+l2w+l0−2pl3,
    g=2A[k3(1+w)−2pk4]+l3(1+w)−2pl4.

For a simple generic root f(w0)=c, the actual moving root is

    w=w0−g(w0)u/f'(w0)+O(u²).

Substitution in η=−du/(u²F_w) gives residue

    (g'f'−f''g)(w0)/f'(w0)³.

Vanishing on the original generic fibres gives the polynomial identity g'f'−f''g=0. Equivalently (g/f')'=0 in C(w), hence g=Cf' for a scalar C. The inference uses a dense set of w0 and does not replace the first correction by a selected or unshifted row.

I independently expanded the complete H²+L and solved these coefficient equations sequentially. Since f' has degree three with nonzero leading coefficient, k3=0 forces C=0 and then g=0. Its value at w=0 is exactly the remaining original-source derivative F_u on u=0, while F_t vanishes there already. This gives actual affine critical points and excludes polynomial mates.

For k3≠0 the successive coefficients yield exactly

    C=k3/(2α),
    k2=2α−4αpk4/k3,
    l3=0,
    l2=−4αpl4/k3.

The entire residue polynomial becomes zero after these substitutions. Since the complete L1 basis is independent, nonconstant L is now equivalent to l4≠0. There is no omitted alternative hidden in a possible cancellation of basis functions.

## 4. Exact critical curve and source-address losses

Set a=4α/k3, q=1+w and

    z=(u²−pa)q−2pu,
    P=αq²+k3uq,
    c0=k0−α−2pk3+pak4.

The complete function, not a truncation, is

    F=(P+k4z+c0)²+l4z+(l0−l2).

Direct differentiation gives

    J(P,z)=u[2pk3−q(k3u+4αq)].

The map (P,z) is not assumed to be an invertible coordinate system. On its nonzero-q critical curve choose

    u=2p/q−aq,
    P=2pk3−3αq²,
    z=a²q³−3paq.

If H denotes P+k4z+c0 along that curve, define

    R(q)=H[k3q+2k4(p−aq²)]+l4(p−aq²).

I independently checked both partial derivatives before restricting to the curve:

    F_u=2R(q),  F_q=(4p/q²−a)R(q).

Thus any root with q≠0 and u≠0 is simultaneously critical, even when the second multiplier vanishes. No division by that multiplier is made. The actual source chart (u,t)→(u,q=1+u²t) has Jacobian u². At an allowed root, t=(q−1)/u² is finite and both original-source partial derivatives vanish. This is exactly the point where polynomiality of the mate, rather than mere rationality, is used.

Choose a nonzero square root s0 of 2p/a. Such a choice exists over C since p,a≠0. Put q=s0Q and

    B=2as0k4/k3,
    C=4+2c0/(pk3),
    Λ=2l4/(k3²s0)≠0.

The symbol Λ here denotes the primary's normalized L parameter, not the original function L. Every admissible coefficient tuple is represented, with either square-root choice. The factor relating R to

    Rbar=[C−3Q²+B(2Q³−3Q)]
           [Q+(B/2)(1−2Q²)]+Λ(1−2Q²)

is pk3²s0/2≠0. Moreover u=as0(1/Q−Q), so the excluded addresses are exactly Q=0,1,−1. There is no additional normalization denominator depending on Q.

## 5. The complete forbidden-root argument

If B≠0, Rbar has degree five and leading coefficients

    −2B²,  5B,  4B²−3.

If every root were among 0,1,−1, the fundamental theorem of algebra, including multiplicities, would give

    Rbar=−2B² Q^e(Q−1)^r(Q+1)^s,
    e+r+s=5,

with nonnegative integers. Put d=r−s and v=r+s≤5. The next two coefficients of the monic product are −d and (d²−v)/2. Hence comparison forces

    B=5/(2d),  d≠0,
    13d²=25v−100.

For v<4 the right side is negative; for v=4 it forces d=0; for v=5 it requires 13d²=25. Every case is impossible. This exhausts all C(7,2)=21 multiplicity patterns and needs no hypothesis about the remaining coefficients C or Λ.

If B=0, the polynomial is

    −3Q³−2ΛQ²+CQ+Λ.

It is genuinely cubic and zero is not a root because Λ≠0. If all roots were ±1, write it as −3(Q−1)^r(Q+1)^s with r+s=3. The constant and quadratic coefficients give

    Λ=−3(−1)^r,
    r−s=2(−1)^r.

The latter equates an odd integer with an even one, excluding all four patterns. Thus there is always at least one allowed finite root. Equations in the preceding section supply an actual source critical point for every surviving coefficient tuple. Together with the earlier cases, this proves the whole p≠0 polynomial-mate exclusion.

The proof does not depend on a Gröbner calculation, a sampled generic critical polynomial, a numerical root enclosure, or a census of mate degrees. The source's 21+4 pattern controls are an exhaustive small algebraic universe supplementary to the coefficient proof.

## 6. Positive and hostile controls

The primary's B=0,C=0,Λ=1 cubic is nonzero at all three forbidden addresses. In contrast, B=0,C=3,Λ=0 gives −3Q(Q−1)(Q+1), whose roots are all forbidden. This exactly tests the nonconstant-L condition needed in the cubic branch; it is not a counterexample to the theorem because constant L was already excluded by the polynomial factor 2H.

I additionally checked a literal original-source example independently:

    p=1/2, α=1, k0=1, k2=2, k3=4, k4=0,
    l0=l3=0, l2=−4, l4=8,
    kt=kxt=lt=lxt=0.

Then a=s0=1 and the normalized cubic is −3q³−2q²+1. Starting with the original polynomial H²+L in u,t, differentiating there, and substituting

    u=1/q−q,  t=(q−1)/u²

makes both derivative numerators divisible by that cubic. Each resulting denominator, and q(q²−1), is coprime to the cubic. This provides an exact algebraic source critical point without an approximate root or an assumption about the critical map's invertibility.

The accepted conclusion remains polynomial-only. The distinguished family already contains genuine rational-mate controls, and criticality of a polynomial first coordinate does not forbid a rational mate with a pole at the critical point. No stronger rational-mate claim is imported into the shifted proof.

## 7. Source review, replays and pins

I read the entire new standalone implementation. It imports no inherited mathematical implementation. Its global minors are literal all-p identities; the infinity block retains arbitrary higher-unit coefficients; the residue block checks all coefficients and the full reconstructed function; the critical-curve block checks both partials and every normalization factor; and the finite pattern block counts exactly 21 quintic and four cubic patterns. All gates are explicit runtime checks and remain active under Python optimization.

Independent commands run in `/tmp/math-wt-planar-jacobian-sep06`:

    python3 04-computation/planar_jc48_sep08_four_four_transport.py
    python3 -O 04-computation/planar_jc48_sep08_four_four_transport.py

Both exited successfully. Their outputs are byte-identical to the frozen output: **151 always-active gates, 420 bytes**. The semantic digest is

    ad47627c90d50ba6afc0079324c4d069634efab507670e15816c229022bdc0f8

| Audited artifact | Bytes | SHA256 |
|---|---:|---|
| Primary `.md`, pre-promotion | 11,057 | `2fba9d86230570e03db27f0c9eb90526b9e93981e469dc4b1b9093c1b80bfa6f` |
| Source `.py` | 10,516 | `4c284b826ed72bcb0e06775488f08c89a1fed1f492bbb7f118d113918d080851` |
| Frozen `.out` | 420 | `a420c162add140c10de64ff4290bf3929102b4aae72ad47d519ab5a1832a0897` |

**Final acceptance: PASS.** Every finite p is covered after combining this new nonzero-p proof with the currently proved distinguished-point theorem. The full global sections, actual weighted infinity analysis, original-fibre residue, noninvertible critical-map use, complete chart exclusions, all complex coefficient degenerations, and polynomial-versus-rational boundary have been checked. No mathematical or source correction remains.
