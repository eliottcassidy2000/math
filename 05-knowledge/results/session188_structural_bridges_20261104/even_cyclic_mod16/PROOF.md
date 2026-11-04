# Every cyclic quotient converse has a block-uniform separator

## 1. The remaining all-even child domain modulo sixteen

Let A be any tournament of even order N>=4, with complete ordered-path-cover profile q and odd Hamilton count H=q1. Normalize every coordinate by the inverse of H modulo16.

The already proved modulo-eight profile formulas imply

    q2/H = 6 or14 (mod16),
    q3/H = 4 or12 (mod16).

For q4, forget the ordering of the four path slots. Every induced tournament has an odd Hamilton count, so

    q4/4! = S(N,4) (mod2).

The Stirling recurrence gives S(N,4)=S(N-1,3) modulo2. For even N>=4 the argument N-1 is odd and at least3, so this parity is one. Since4!=24=8 modulo16, q4=8 modulo16; multiplying by the inverse of any odd H retains8.

Similarly q5 is divisible by5!=120 and thus is either0 or8 modulo16. Normalization again leaves these possibilities unchanged. For every k>=6,16 divides k!, so qk=0 modulo16.

Consequently every actual even-order child of order at least four has normalized profile in the eight-element domain

    r(b,c,d)=(0,1,b,c,8,d,0,0,...),
    b in{6,14}, c in{4,12}, d in{0,8}.                    (1)

This is a complete upper domain, not a claim that all eight formal residue vectors are realized by tournaments. Any further relations among actual profiles only shrink it. H modulo16 has not been discarded: it is the invertible scalar used in the normalization.

## 2. Exact full physical run functional

Use three disjoint nonempty child name sets A,B,C with cyclic quotient A->B->C->A. Compare T=C3[A,B,C] and U=C3^op[A,B,C], with the same physical child interiors. For each pair(a,c) of three-bit masks, attach actual ports x->y uniformly on whole child sets.

The exact complete run-word formula is proved and independently audited in `../coarse_response_rank/PROOF.md`, Section1. To fix all conventions, put

    R=1-XYZ,
    N=[X,XY,XYZ; XYZ,Y,YZ; XZ,XYZ,Z].

The entries of N/R count all nonempty quotient run words by first and last child; quotient converse transposes N. If o=1-a,p=1-c, set

    L_a=R+1^t N a, R_p=R+p^t N1, M_(o,c)=o^t N c,
    W_(a,c)(N)=L_a(R+M_(o,c))R_p+L_c M_(p,a)R_o,
    D_(a,c)=W_(a,c)(N)-W_(a,c)(N^t).

The first term permits x before y with an empty intervening segment; the second permits y before x only with a nonempty segment. Prefixes and suffixes may be empty. Repeated visits to each child are all retained. Applying

    L(X^i Y^j Z^k)=q_i(A)q_j(B)q_k(C), q0=0,

gives the exact actual contrast

    Delta(a,c)=L(D_(a,c)/R^3).                            (2)

In particular this is the full physical two-port response, not an endpoint approximation or a quadratic/cubic truncation. Every run word uses exactly one cover coordinate from every child.

Normalize each actual child as in(1). By multilinearity,

    Delta(a,c)=H(A)H(B)H(C) L_r(D_(a,c)/R^3) (mod16).     (3)

The scalar is odd and independent of the query. Since every coordinate above index5 vanishes modulo16, the denominator expansion

    R^-3=sum_(k>=0) binom(k+2,2)(XYZ)^k

can be truncated exactly at k=5 in(3). Higher terms have a zero profile factor. No asymptotic estimate or numerical precision is involved.

## 3. Complete eight-cubed certificate

There are exactly8^3=512 ordered triples of the formal profiles(1). For every triple, the verifier evaluates all64 exact integer functionals in(3), reducing only the final residues modulo16. The complete table certifies:

- Every residue is either0 or8.
- No triple has an identically zero coarse bank.
- Exactly256 triples have24 nonzero queries,192 have28, and64 have12.

More strongly, the fixed query(a,c)=(1,3), meaning a=A,c=A union B, has residue8 for every one of the512 triples. The six queries(1,3),(1,5),(2,3),(2,6),(4,5),(4,6) are all universally nonzero on this upper domain. The checker also retains the first nonzero query in lexicographic order for every triple. This is exhaustive over the entire upper domain, including any unrealizable formal vectors. Its implication for actual tournaments is therefore one-way and valid: every actual profile triple is covered. No formal profile is being passed off as an actual tournament or a response collision.

Equation(3) turns every selected formal residue8 into the same actual residue8, since multiplying8 by an odd number leaves it unchanged modulo16. Expand the fixed quotient mask(1,3) onto every physical name of each original child. This supplies an actual separator for every triple of even-order children of order at least four, with arbitrary unequal profiles and interiors.

## 4. Local completion for all nonempty cyclic children

The independently audited `../five_profile_types/PROOF.md` supplies a modulo-eight coarse separator for every triple of nonempty children outside the all-even-order-at-least-four case. Section3 covers precisely that remaining case modulo16. Therefore:

**Theorem.** For every three nonempty tournaments on fixed disjoint physical name sets,

    F_(C3[A,B,C]) != F_(C3^op[A,B,C]),

and some block-uniform two-port query already separates them. The source and target keep the same child interiors; arbitrary profile-preserving target replacements give the same conclusion by the run functional. All child orders, orientations and unequal complete profiles are allowed. There is no exceptional child-size pattern.

The all64-query bank is a fixed finite battery on a known three-child partition. If the source complete profiles are known, the residue tables select a single separating block query. This is not reconstruction between unrelated quotients or arbitrary tournaments, and it is not a claim that the old single canonical query always works.

## 5. Root-family consequence and the ambient boundary

For a fixed source, pair a selected baseline query z with its core converse-complement Jz=(1-c,1-a). The exact full augmented-profile duality in `../monotone_dual_queries/PROOF.md` gives opposite baseline integer contrasts. Choose the member having negative source-minus-baseline contrast.

Every original-module XOR family containing the cyclic root reverses its quotient and makes every target child cover profile coordinatewise no smaller, by hereditary original-module monotonicity. Because the ports are uniform on each child, nonnegative run substitution preserves the selected negative direction. Thus the pair{z,Jz} detects every such root-containing family, without a scalar-query fallback, and one baseline-sign-selected member works for all of them. Target child profiles may be larger and unequal. This is source-relative detection, not pairwise injectivity among arbitrary orbit outputs.

This packet's new modulo-sixteen conclusion is LOCAL to the cyclic core. No arbitrary-ambient lift is inferred merely from its nonzero residue. The earlier modulo-sixteen transfer failure remains a real warning. The separate `cyclic_odd_profile_lift` theorem supplies arbitrary-ambient separation for the other nontrivial type triples, but the all-even-order-at-least-four ambient case is not settled by this packet.

## 6. Reproducible controls

The verifier retains the full512-by64 formal integer contrast/residue table and all numerator polynomials. It checks actual complete normalized profiles, literal full physical query banks on actual even children, unequal-profile controls, the entire named bidegree-at-most-(2,2) bank for C3[TT4,TT4,TT4], all zeros and cubic/quartic coefficients, and full-mask reconstruction controls. Actual root-family examples retain their literal full-query values and profile-growth checks.

The finite formal classification proves a complete algebraic upper-domain statement. Independent physical counts and an independent run-word audit are separate checks on the implementation and physical conventions; they do not replace the completeness argument in Sections1–3.
