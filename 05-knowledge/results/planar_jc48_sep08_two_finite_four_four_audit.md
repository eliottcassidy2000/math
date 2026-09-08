# Independent audit of two distinct finite fourfold roots

**Status: INDEPENDENT ANALYTIC / SOURCE / FINITE-EXACT AUDIT PASS.**
The [two-finite-root proof](planar_jc48_sep08_two_finite_four_four.md), the complete frozen source and both fresh execution modes are accepted without correction. The coordinating agent owns promotion and integration.

The accepted statement excludes **every rational Jacobian mate** for the full fixed-DG global class H∈L2, L∈L1, F=H²+L with N|S=α(x−p)⁴(x−q)⁴, α≠0 and two distinct finite p,q. Constant L is included. It neither assumes geometric integrality of the generic fibre nor bounds a proposed mate's degree. The stronger rational assertion is not transported to a root at infinity or to coincident roots.

The [finite/infinity 4+4 supplier](planar_jc48_sep08_four_four_transport.md) and its distinguished-point dependency are now **PROVED / INDEPENDENTLY AUDITED** on disk. They are not dependencies of the two-finite rational proof. Combining their polynomial conclusion with this proof exhausts the polynomial-mate question for this exact two-distinct-fourfold-root boundary partition on P1. No surface automorphism normalizing arbitrary boundary positions is needed.

## 1. Complete global base and correction spaces

Use the original surface and source

    W=(P1_x×P1_z)\{z=x²},
    s=z−x², t=1/s, ω=dx∧dt.

Let r=(x−p)(x−q), S1=p+q and P1=pq. The new base function is exactly

    Q=r²t+x²−2S1x
      =h−2S1v+(S1²+2P1)w−2S1P1xt+P1²t,

where w=x²t, v=x+x³t and h=x²+x⁴t. This verifies Q∈L1 using the original global basis. Clearing the denominator gives sQ=r²+s(x²−2S1x), a section of bidegree at most (2,1) after substituting s=z−x². Its graph restriction is r², with nonzero degree-four leading coefficient.

The L1 restriction map has rank five onto all polynomials of degree at most four. In the basis (1,t,xt,w,v,h), the five nonconstant columns restrict to (1,x,x²,x³,x⁴), and the constant-function column restricts to zero. Its kernel is therefore exactly the constant functions.

For L2, the full section box has dimension fifteen and graph restriction has rank nine. In shifted coordinates u=x−p, the monomials

    1,u,u²,u³,u⁴,u³z,u⁴z,u³z²,u⁴z²

restrict to successive degrees zero through eight with leading coefficient one, so the rank-nine minor is identically one for every p. The six-dimensional kernel is s times the full L1 space. Since αQ² has exactly the prescribed octic numerator, every entry is

    H=αQ²+K,  K∈L1,  L∈L1.

This is an exact complete section space, not an ansatz for selected lower rows. It uses only affine source coordinates for x−p and retains the original graph and volume.

Writing K=A_K(x)t+B_K(x), L=A_L(x)t+B_L(x) and B_Q=x²−2S1x, the full cleared numerators are

    N=αr⁴+s(2αr²B_Q+A_K)
       +s²(αB_Q²+B_K),
    M=A_L+sB_L.

The polynomials A_K,A_L have degree at most four. At either root of r, the r² term and its first tangential derivative vanish. Consequently the normal value and first normal tangential jet of N are exactly the corresponding value and derivative of A_K; the two lower jets are those of A_L. No derivative of the global base Q is omitted in making this identification.

## 2. The two-shared-root composition statement is complete

If both points are singular shared roots and a rational mate is proposed, the finite first-jet obstruction forces A_K,A_L to have double zeros at p and q. Distinctness then implies

    A_K=βr²,  A_L=δr².

Their degree bounds are used here. Equivalently, the confluent value/derivative matrix on (1,x,x²,x³) has determinant (p−q)⁴≠0. I independently reconstructed this determinant. The unique degree-four direction left by all four jets is r².

The restriction-kernel result in §1 now pays the **whole functions**, not merely their boundary parts:

    K=βQ+γ, L=δQ+ε,
    F=(αQ²+βQ+γ)²+δQ+ε.

The latter is a polynomial of degree four in Q with leading coefficient α². It excludes polynomial mates via its nonconstant derivative factor, but that factor alone would not exclude rational mates. The main proof below supplies the stronger rational obstruction and needs only one singular shared root; it does not misuse this composition observation as a rational-exactness theorem.

## 3. Exhaustive local alternatives with higher-unit terms retained

I checked the precise current statements in the proved [common-root necessity](planar_jc48_sep08_quartic_common_root.md), §2, and [shared-root analysis](planar_jc48_sep08_shared_roots.md), §§2–3. They are valid for the full local analytic coefficients and every normalized generic branch. Neither lemma assumes the order-four leading term is a pure monomial with all higher coefficients zero.

Near p, with u=x−p and d=p−q≠0, the boundary term is

    N(0,u)=αu⁴(u+d)⁴.

Its order is four and its leading coefficient is αd⁴≠0. The factor (u+d)⁴ is retained as a unit; at q the same argument applies after exchanging the labels. This is a relabeling of the symmetric global data, not a transformation of the surface.

At either finite point the local fibre equation and relative form have the structure

    ℰ=N²+s³M−cs⁴,
    η=s²du/ℰ_s

up to the already fixed nonvanishing finite-coordinate multiplier. The actual finite source coordinates may be kept throughout, so no infinity weight is imported into this argument.

* If M is a unit, the M-unit classification has m=4. Its only possibly bad balance would be 4=3j, impossible for an integer j. For j=1 the simple low branch and two cancellation determinations are regular; the latter have unweighted exponent 1/2 in u and hence order two under their quadratic normalization. For j≥2, including infinite order, the dominant cubic determinations are regular. The j=0 case is also regular and can instead be covered by the normal-unit lemma.
* If N_s is a unit, the exact analytic centre N=0 has order four. At that centre ℓ=ord(M−cs) is generically at most four because the centre is independent of c and has a nonzero order-four coefficient. The form has exponent (4−ℓ)/2 in u, regular even after ramification. This covers arbitrary M. All higher-unit corrections merely change the exact centre and do not remove the c-coefficient bound.
* The remaining case is M=N_s=0. Since N|S has order four, its tangential derivative also vanishes. If either first tangential jet of N_s or of M is nonzero, the shared-root tangent face has, after s=uZ, a nonzero simple root away from Z=0. For the normal first jet it is one of the two roots of a quadratic factor; for the lower first jet alone it is the nonzero root of a linear factor. At such a branch ℰ_s has order three, s²du has order two, and the residue coefficient is a unit times Z0²/P_c'(Z0)≠0. This is a genuine finite simple pole and forbids any rational primitive. No classification of the remaining zero-face branches is needed to obtain that obstruction.

Therefore every hypothetical rational mate has, at each of the two points, either a regular unit case or the complete jet vanishing

    A_K=A_K'=A_L=A_L'=0.

The cases cover every coefficient tuple. They do not assume that the mate or tuple is generic; the generic qualifier refers only to choosing the fibre value c outside finitely many exceptional values for that fixed tuple.

## 4. Compact exactness when both points are regular

The boundary octic has degree eight and nonzero leading coefficient α, so the distinguished infinity point is not a boundary zero. The two finite roots are the only points of S met by a compact generic fibre. If both have the regular unit alternatives above, all boundary branches of that fibre have regular relative form.

Inside W the original volume is regular, with a zero of order two along the other chart divisor D. Generic fibres are smooth in characteristic zero, and restriction of ω/dF is regular there independently of whether a proposed rational G is globally regular on W. Possible fixed critical values may be discarded. A rational primitive cannot have a pole at a point where its differential is regular.

Every compact generic component meets the original source plane, where the volume is nonzero. A component contained in D would force a fixed value of F|D rather than a generic fibre; S is not a component of the fibre equation because its restriction is the nonzero octic squared. The relative form is thus nonzero on every such component. A holomorphic differential on a compact smooth component cannot be the derivative of a rational function unless it is zero: any pole of the rational function would give a differential pole of higher order, and a pole-free rational function is constant.

This is a componentwise argument and assumes no geometric-integrality or relative-constant-field simplification. A proposed G∈C(x,t) restricts to the generic fibre and to its geometric components; its denominator does not vanish identically on a component of a generic level. Consequently the both-regular alternative is impossible for a rational mate. At least one point must have the full jet vanishing, and the next computation rules out that alternative as well.

## 5. One active point and the exact original-fibre residue

Choose a point with the four vanishing jets and call it p. In the full corrected global basis,

    w=u²t, v_p=u(1+w)−2p, h_p=uv_p,

those conditions remove exactly the coefficients of t and ut from K and L. Indeed their cleared boundary restrictions have successive coefficients (kt,kxt,k2,k3,k4) and (lt,lxt,l2,l3,l4). Thus the complete remaining spaces are

    K=k0+k2w+k3v_p+k4h_p,
    L=l0+l2w+l3v_p+l4h_p.

No constant-L or k3-zero branch is removed. The exact global Q has expansion

    Q=Q0(w)+uQ1(w)+u²(w+1),
    Q0=d²w−p²−2pq, Q1=2dw−2q.

I independently substituted x=u+p, t=w/u² into the original Q and recovered all three rows. With

    A=αQ0²+k2w+k0−2pk3,
    B=2αQ0Q1+k3(1+w)−2pk4,

the complete H and F expansions give

    H=A+uB+O(u²),
    F=f+ug+O(u²),
    f=A²+l2w+l0−2pl3,
    g=2AB+l3(1+w)−2pl4.

The leading coefficients are invariant under every remaining correction:

    [w²]A=αd⁴,
    [w²]B=4αd³,
    [w³]f'=4α²d⁸,
    [w⁴]g=8α²d⁷.

The new degree-four first correction is supplied by the other **finite** root. This is the exact information that would be lost by replacing the full global base with a local monomial model.

For a generic c, f(w)=c has four distinct nonzero roots. The polynomial f has degree four and nonzero leading coefficient, so its derivative and zero-root exceptions exclude only finitely many c. These are actual branches: after the vanishing jets, the full local Weierstrass degree in s is four, and the exact tangent equation is

    [u⁸]ℰ(u,u²Z)=Z⁴(f(1/Z)−c).

Its constant term is α²d⁸ and its generic leading coefficient is H(0)²+L(0)−c. Thus the four tangent roots are precisely the reciprocal f-roots, with no omitted zero or infinite root. The local substitution s=u²/w is regular at each nonzero w0. The implicit-function branch w=w0+… is an actual normalized branch of the original compact fibre.

On it,

    w=w0−g(w0)u/f'(w0)+O(u²),
    η=−du/(u²F_w),
    Res η=(g'f'−f''g)(w0)/f'(w0)³.

This retains the induced displacement of the original root. Omitting that displacement would omit the f''g term. I independently formed the full original polynomial H²+L, extracted f and g, and recovered the same residue numerator and degrees.

If a rational primitive existed, every one of these residues would vanish for generic c. Equivalently g'f'−f''g would be identically zero and (g/f')' would vanish in C(w), forcing g=Cf' for a scalar. But g has degree four while f' has degree three. More explicitly,

    deg(g'f'−f''g)=6,
    [w⁶](g'f'−f''g)=32α⁴d¹⁵≠0.

The coefficient is 128α⁴d¹⁵−96α⁴d¹⁵ from the two products. It cannot be adjusted by any k or l coefficient, including constant L. A fixed nonzero numerator polynomial has finitely many roots; excluding their f-values ensures that **none** of the four generic roots has zero residue. Hence at least one, and in fact each of these generic local branches, contradicts a rational primitive. This closes the last local alternative without a source critical-point argument.

## 6. Hostiles and an independent original-fibre check

The finite/infinity example

    w=x²t, H=w², L=w, F=w⁴+w,
    G=1/[x(4w³+1)]

has J(F,G)=1 and boundary partition 0⁴∞⁴. The coalesced example

    h=x²+x⁴t, H=h², L=h, F=h⁴+h,
    G=1/[3x³(4h³+1)]

also has J(F,G)=1 and boundary octic x⁸. These are literal rational controls, not only formal differential models. The first prevents extending the rational conclusion by forgetting infinity's volume weight; the second confirms that distinctness is essential. In the residue supplier, d¹⁵ vanishes on the coalescence boundary exactly as it should.

For the composition subfamily, a generic fibre of Q is parametrized rationally by x, with t=(Q−B_Q)/r². Its relative form is

    −dx/[(x−p)²(x−q)²].

I independently computed its separate residues as 2/(p−q)³ and −2/(p−q)³. Their sum is zero but neither is zero. The sum therefore carries no exactness certificate; the proof correctly keeps each labelled pole separately.

The new source adds the exact literal p=1, q=2, H=Q², L=0 on F=Q⁴=1 at w0=6, where Q=1. Its residue is −1/2. This independently agrees with the Q-fibre residue −2 divided by 4Q³=4. It directly checks that the formula concerns the original F-fibre, including a geometrically split composition example, rather than assuming a primitive generic fibre. Constant L is thus positively included, not hidden in a discarded coefficient case.

## 7. Source completeness, independent replays and pins

I read the entire frozen source after final pins arrived. It imports no inherited mathematical implementation. The exact checks cover the global Q identity, actual bidegree boxes, rank-five and rank-nine restrictions, complete corrected basis, confluent two-point determinant, all normal jets, all original first rows, the inverse-coordinate tangent identity, degree-six residue, constant L and both rational hostiles. The four named finite point pairs include a complex pair. Those are supplementary controls; the universal conclusion follows from symbolic identities with leading coefficient 32α⁴(p−q)¹⁵ and the analytic boundary exhaustion.

Independent commands, run from `/tmp/math-wt-planar-jacobian-sep06` after source freeze:

    python3 04-computation/planar_jc48_sep08_two_finite_four_four.py
    python3 -O 04-computation/planar_jc48_sep08_two_finite_four_four.py

Both runs exited successfully and are byte-identical to the frozen output: **60 always-active exact gates, 423 bytes**. The semantic digest is

    3e1231865784dc4ad39f4ad9ad7cf13e17261b431c1d6df7e34c178c519c5d17

| Audited artifact | Bytes | SHA256 |
|---|---:|---|
| Primary `.md`, pre-promotion | 12,024 | `90c6525972e2568bd3b735fce344406dee3bbb99253a2190a2ffc532152186a3` |
| Source `.py` | 9,639 | `2a901ed396b72a7b7eb14f6dc34f30d53e82079c454f22027b86528dca69b630` |
| Frozen `.out` | 423 | `67510ca460a435fefef40ab9da962924ab7c9956f078ed88ee60be88c6e23f05` |

Additional independent symbolic calculations reconstructed the original Q and complete F first rows, verified the residue degree and leading coefficient, recovered the confluent determinant, and computed both individual Q-fibre residues without importing the producer.

**Final acceptance: PASS.** The two-finite-root rational exclusion, including constant L and all geometric generic components, is proved by the stated argument. The currently proved finite/infinity supplier discharges the integration condition for the all-P1 two-distinct-fourfold-root **polynomial** corollary when this primary is promoted. Other boundary partitions, unrestricted global quartics, and JC(2) remain outside the conclusion.
