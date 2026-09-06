# Reciprocal zero swap: an exact integer-line criterion and nonresidue obstructions

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** For the completed reciprocal permutation over an odd prime, exchanging outputs 0 and 1 has an exact divisor-and-integer-lift criterion. Its actual integer triple count is twice the number of primitive arithmetic packets defined below. A nonzero quadratic-residue discriminant is necessary for each packet, but the two standard root positions are essential.

The criterion produces three explicit infinite prime families where the swap fails even though 5 is a quadratic nonresidue. Their infinitude uses **CITED Dirichlet**; the coordinates and failures are proved directly by integer identities. This concerns a p-point permutation graph, not a 2p-point no-three-in-line construction. It gives neither an all-prime successful single-transposition construction nor a classification of all output swaps.

[Independent proof and exact audit](continuing6_20260906_reciprocal_native_audit.md) passes.

## 1. Inheritance and the operation board

The closest proved mechanism is [continuing5_20260906_wildcard_swap_compiler.md](continuing5_20260906_wildcard_swap_compiler.md), Section6, with its independent audit. Its completed reciprocal graph has exactly one integer triple, although it has (p-1)/2 modular triples. Its exact general zero-swap decoder retains native integer products and coordinate-sum carries. The canonical hostile p11 defeats every swap of zero with another output. Thus universal repair by that restricted move is already refuted.

The present operation is narrower: swap outputs 0 and 1 only. The corrected near miss is treating the discriminant-five obstruction as an iff, or treating a split quadratic as an actual integer line. The least-used sidecar is the primitive rational slope with the standard integer lifts of its modular intersection roots.

Anchor: the actual integer no-three predicate after one specified move. Niche: divisor and carry packets for the line's two conic intersections. Wildcard: fixed packets that produce arithmetic progressions of obstructions. The board has five entries:

| Object | Map or operation | Retained predicate | Lost information / required sidecar |
|---|---|---|---|
| Reciprocal permutation | Exchange outputs0,1 | One point per row and column | Actual integer triple set |
| Triple through the moved point | Reduce its slope to a/b | Its complete integer line | Primitive slope gauge |
| Line-conic intersection | Quadratic modulo p | Modular collinearity | Standard root positions in a bounded interval |
| Primitive line | Integer tuple(h,k,l,c) | Exact native triangle, bijectively | Divisibility and native upper bound |
| Fixed tuple | Vary p in a compatible progression | Explicit native failure | Prime hypothesis and CITED infinitude |

This is purposeful deterministic relabeling. The conditional-column and restart results in [overnight11_20260906_no3line_rowfreeze.md](overnight11_20260906_no3line_rowfreeze.md) and [overnight12_20260906_no3line_count_restart.md](overnight12_20260906_no3line_count_restart.md) require fresh uniform columns conditional on history; their laws are not applied to this selected move. Cooper's reciprocal example and Guy--Kelly are inherited provenance, not new literature or independence assumptions.

## 2. Exact primitive packet theorem

Let p be an odd prime. After the specified output swap, the board is

    B_p = {P=(0,1), Q=(1,0)}
          union {(x,[x^(-1)]_p): 2<=x<=p-1}.

Define N_p to be the number of integer tuples (a,b,h,k,l,c) satisfying

    a>=1, h>=5, a h=p-1,
    2<=k<l, k+l=h, 1<=c<=k-1,
    b=(h+c p)/(k l) is a positive integer,
    gcd(a,b)=1, b l<=p-1.                              (1)

**Theorem.** The tuples in (1) are in bijection with the actual integer-collinear triples through P. Their two other points are

    R=(b k,a k+1),  S=(b l,a l+1).                     (2)

Every remaining actual triple is the transpose of one of these, through Q. Consequently

    X(B_p)=2 N_p,     B_p is successful iff N_p=0.      (3)

In particular the original one-defect graph is either repaired completely or changed to an even positive number of defects. This parity statement refers to the specified swap, not to arbitrary swaps.

**Proof.** Three untouched points cannot be collinear: they lie on the nonsingular modular conic xy=1, which a line meets in at most two points. A triple containing P and Q would require an untouched point with x+y=1, impossible since both coordinates are at least2. The board is invariant under transposition, exchanging P and Q. It remains to classify triples through P.

Such a line has a unique primitive equation

    b(y-1)=a x,       a,b>0, gcd(a,b)=1.

All untouched coordinates satisfy x,y>=2, so the slope is positive. Its two points have the form (b k,a k+1),(b l,a l+1), with distinct positive integers k<l. Reducing the equation and xy=1 modulo p gives b y^2-b y-a=0. The two distinct y values therefore have sum1 modulo p. Their standard integer sum lies between4 and2p-2, so it is exactly p+1. Hence a(k+l)=p-1. Put h=k+l.

The reciprocal equation at k also gives

    b k l = h (mod p).

Write bkl=h+c p. Since bkl-h>=kl-k-l=(k-1)(l-1)-1>=-1, and p>=3, the integer c is nonnegative. Equality c=0 would force b=1/k+1/l. For k=1<l this lies strictly between1 and2; for k>=2<l it lies strictly between0 and1. Neither is an integer. Thus c>=1. The native bound b l<=p-1 gives h+c p<=k(p-1), hence c<k. Therefore k>=2, l>=3 and h>=5. All conditions in (1) follow.

Conversely suppose (1). Both points in (2) have positive native coordinates at least2. The x upper bound follows from b l<=p-1. Since l<=h-2, the y upper bound is a l+1<=p-2a<=p-2. Moreover

    b k(a k+1) = b k(p-a l) = 1 (mod p),

because a h=p-1 and bkl=h+c p. The same holds at l. Thus both are untouched reciprocal points in the standard box, and they lie on the displayed actual integer line. Primitive slope and the ordered k<l recover the tuple uniquely. Transposition gives a disjoint second family of triples, proving (3).

Two useful necessary bounds are

    a <= (p-1)/5,       1<=b<2a.                       (4)

The second follows from b h=x_R+x_S<2(p-1)=2a h. The first is sharp: at p61 the packet (a,b,h,k,l,c)=(12,11,5,2,3,1) gives the points (22,25),(33,37) through P. Thus the short-cofactor threshold five is attained, not merely a convenient search cutoff.

## 3. The equivalent character test needs a native root interval

For every positive divisor a of p-1 with a<=(p-1)/5, and every b with 1<=b<2a and gcd(a,b)=1, put

    h=(p-1)/a,  H=min(h-1,floor((p-1)/b)),
    f_(a,b)(z)=a b z^2+b z-1,
    Delta=b(b+4a) modulo p.                            (5)

There is a native P-triple with primitive slope a/b **iff** Delta is a nonzero square modulo p and the two standard roots of f_(a,b) both belong to the integer interval [1,H]. Equivalently, one standard root belongs to [h-H,H], with the nonzero-square requirement retained.

Indeed the polynomial has two distinct roots exactly under the discriminant condition, and their sum is -a^(-1)=h modulo p. If both standard roots are at most H<h, their positive sum is at most2h-2; the only possible lift of that residue is h. Conversely a native pair gives exactly these bounded roots. The interval inequalities are precisely the native coordinate bounds in (2), so the packet theorem applies. The equivalent one-root interval supplies the other root as h-k, a positive standard representative. No choice of modular square-root sign affects the test.

The equivalent polynomial in the actual x coordinate is a x^2+b x-b. Its discriminant is the same Delta, but multiplying a standard z root by b need not preserve its standard representative. This is why the z-interval and coordinate map are stated explicitly.

**Hostile to dropping the interval.** At p13 take a2,b1. Then h6,H5 and Delta9 is a nonzero square, but the z roots are7,12. The corresponding standard reciprocal points are (7,2),(12,12), and

    det(P,(7,2),(12,12))=65=5p, not0.

The actual repaired board at p13 has no integer triple. Character alone would falsely reject a successful repair. The missing coordinate is the integer lift, not a larger modular rank or another discriminant.

## 4. The discriminant-five family, and why it is not an iff

If p!=5 and 5 is a quadratic residue modulo p, the polynomial x^2+x-1 has two distinct standard roots x_1,x_2. Neither is0,1 orp-1. Their exact sum is p-1. Since x_i(x_i+1)=1 modulo p, the standard points

    P=(0,1), (x_1,x_1+1), (x_2,x_2+1)

lie on the actual line y=x+1. Transposition gives a second actual triple. Thus this character condition is sufficient for failure. At p11 the roots3,7 give the first example. At p5 the discriminant is zero and the single root2 gives only one conic point on the line, not a triple; the repaired board succeeds.

The converse is false. The first prime with character(5/p)=-1 and a failed repair is **p37**, with

    P, (24,17), (30,21),     primitive slope2/3.

The exact native packet is (2,3,18,8,10,6). The bounded minimality claim checks every smaller odd prime; it does not extrapolate an eventual failure theorem. The character-negative primes3,7,13,17,23 are genuine successful controls.

## 5. Fixed native packets give infinite nonresidue obstruction families

Here is a general construction template, allowing an unreduced displayed slope. Fix integers 2<=k<l, h=k+l and 1<=c<k. Suppose p is an odd prime satisfying

    p=1 (mod h),    k l divides h+c p,
    p(k-c)>=h+k.                                      (6)

Put a=(p-1)/h and b=(h+c p)/(k l). Formula (2) gives two untouched native reciprocal points through P. The last inequality in (6) is exactly b l<=p-1. The preceding proof works without gcd(a,b)=1 for existence; only the bijective count needs primitive normalization. Compatible reduced arithmetic progressions satisfying (6) for all sufficiently large terms therefore give explicit prime obstruction families.

Three concrete families, all with character(5/p)=-1 by quadratic reciprocity, are:

| Prime progression | (h,k,l,c) | a | b |
|---|---|---|---|
| p=37 mod360 | (18,8,10,6) | (p-1)/18 | 3(p+3)/40 |
| p=43 mod70 | (7,2,5,1) | (p-1)/7 | (p+7)/10 |
| p=97 mod120 | (8,3,5,1) | (p-1)/8 | (p+8)/15 |

Every positive term in these progressions meets the native box inequality. Each progression is coprime to its modulus. Its infinitely many prime terms follow from **CITED** [Dirichlet's 1837 theorem, English translation, arXiv0808.1408v2](https://arxiv.org/abs/0808.1408v2), whose primary theorem and opening page were checked for this use. Infinitude is the sole external input here; each actual triangle follows by the exact integer identities above.

For example write p=43+70r, r>=0. The points are

    R=(10+14r,13+20r),  S=(25+35r,31+50r),
    R_x R_y=1+p(3+4r), S_x S_y=1+p(18+25r).

All coordinates lie in2,...,p-1 and both points lie on b(y-1)=a x with a6+10r,b5+7r. Analogous coefficientwise identities for all r>=0 in the other two progressions are retained in the certificate. This verifies an unbounded family, not a bank of tested primes.

Primitive normalization is substantive for counting. At p113 the h7 template has displayed a16,b12. Dividing the slope by4 changes the primitive tuple to (a,b,h,k,l,c)=(4,3,28,8,20,4). Reusing the unreduced h7 tuple in the bijection would be wrong. Also these constructions force at least two actual triples, not exactly two: p113 has four.

## 6. Consequence, connections and stopping boundary

The new progress is an iff native decoder in primitive arithmetic packets, its equivalent character-plus-root-interval form, the sharp short-cofactor exclusion h<5, and explicit infinite failures beyond the discriminant-five character. It provides a way to generate and explain obstructions before inspecting complete boards.

The LRC connection is precise: a character or modular-root packet plays the role of a marginal overlap description, while the bounded standard roots play the role of the common native translate. The p13 modular determinant65 is the cheap hostile to forgetting that coordinate. Primitive line normalization parallels attaching one compatible owner to one actual LRC edge; the p113 factor-four change prevents counting the same geometric object in the wrong gauge. In the moments lane, split or real roots similarly do not determine the original evaluated quantity without its native normalization and interval data. These are mechanism comparisons, not theorem transfers between different predicates.

No infinite successful family among the remaining character-negative primes is proved. The divisor structure of p-1 and the locations of the associated quadratic roots both vary, so the character bit alone does not compactify the remaining question. A meaningful next test is a specified arithmetic class in which every packet in (1) can be excluded, or a second prescribed transposition that destroys all retained line packets while controlling newly created triples. Another blind prime census is not the present conclusion.

## 7. Exact evidence and reproduction

The standalone source uses no repository imports. For all odd primes through43, plus the prescribed controls61,97,113,197, it compares (i) actual primitive-slope grouping, (ii) the complete arithmetic packet enumeration, (iii) the quadratic-root interval criterion, and (iv) every literal integer triangle of the original repaired board. It checks the positive, ramified, character-only and primitive-gauge hostiles. The three infinite templates are verified by coefficientwise polynomial identities and native-box inequalities for every parameter r>=0.

Run `continuing6_20260906_reciprocal_native.py` with ordinary Python and `python -O`. Both pass **1,741 always-active exact gates**, with byte-identical raw LF output and JSON. The source writes its certificate adjacent outside the repository or to the results directory in the filed layout.

- Source SHA256: `d0981913961389520200a8215b1c2f8611ea7090bd379ba2024ceb9a6b345c11`.
- Output SHA256: `13f7a46e9f034917337233db012466b499eab3d579b8bc18955e7ae18f0279e2`.
- JSON SHA256: `1da5e84836de4aee72e0636e792a43715c9cb85eac9e621203e215ac7e228380`.

The producer was prepared in an isolated outside directory. Its independent audit passed; the continuing6 manifest pins the filed artifacts.
