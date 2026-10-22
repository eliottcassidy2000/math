# One-apex tournaments: apex classification, class counts, and a product spectral law

Session26,22 October2026. Status: PROVED by the elementary arguments below, with independent symbolic review and finite controls; final-file independent review is recorded separately. No literature-priority claim. Root supplied the classification/product-law derivations; this lane independently checked and writes their full proof. H always means the actual number of directed Hamiltonian paths, not an orbit-normalized count.

## 0. Scope and inherited interfaces

A tournament is one-apex if deleting some vertex leaves a transitive tournament; this includes transitive tournaments. Its feedback vertex number is at most1. A pointed one-apex tournament has a distinguished deletable apex. After deleting it, the remaining transitive tournament has a unique linear order, so pointed isomorphism classes on n vertices correspond exactly to binary words w=b0...b_(m-1),m=n-1, with apex x pointing to positionj iff bj=1. H(w) is computed by the exact recurrence in PROOF.md. Strong words are exactly those beginning1 and ending0. Leading0/trailing1 positions are source/sink padding and may be stripped without changing H.

The independently audited Session25 modular-spectrum dependency is:

`05-knowledge/results/session25_spectrum_wiring_20261021/session25-spectrum/PROOF.md`, Sections2--4; independent audit in sibling `session25-spectrum-review/REVIEW.md`.

It proves every parity-admissible residue modulo every positive integer occurs as actualH even among strong one-apex tournaments, using the exact reset blocks(1^L0^L)^t. This packet does not reprove or claim that modular-completeness theorem. It proves structural limitations and asymptotics for the larger one-apex family.

CUT_ADDENDUM.md proves the exact directed-cut/Ferrers formula and extrema. In particular a strong word of lengthm>=3 has H>=2^(m-2)+1, and its maximum is2^m-2^floor(m/2)-2^ceil(m/2)+3<2^m. These provide honest length cutoffs below.

An ordinary class identifies tournaments under isomorphism only. A merged class additionally identifies converse tournaments. A STRONG CORE class has no source/sink padding. The H-bounded counts used below run over strong cores of all orders. If arbitrary padding is allowed, even the single valueH1 has infinitely many classes of different orders; hence the corresponding all-order, all-class bounded-H counts are infinite. Nothing below calls those infinite counts finite.

## 1. Deletable apices and rigidity

Let D(T) be the set of vertices whose deletion leaves a transitive tournament. For a nontransitive tournament,

    D(T)=intersection of the vertex sets of all directed triangles.

Indeed a tournament is transitive iff it has no directed triangle. Deleting v kills every directed triangle exactly when v belongs to every such triangle. Thus |D(T)|<=3.

If |D(T)|=3, there is exactly one directed triangle, namely D(T). Every outside vertex uniformly dominates or is dominated by that triangle: any mixed orientation to a directed triangle creates another directed triangle, by following a transition around its cycle. Contract the triangle. A directed triangle in the quotient would lift to another triangle in T, so the quotient is transitive. Thus T's sole nontrivial strong component is this C3. In particular a strong tournament with three deletable apices is exactly C3.

Now suppose T is strong of ordern>=4 and has two distinct deletable apices x,y. The remaining vertices R form a transitive chain. Since T-x and T-y are transitive, each of y and x inserts at a cut of this same chain. Relabel x,y so x->y. Let the cuts of x,y be a,b, meaning the first a chain vertices dominate x and the first b dominate y. Cyclic triangles can only be i->x->y->i, with b<i<=a. If b>=a there are none, contradicting strongness. If b>0, the first b vertices dominate all remaining vertices and form a source prefix; if a<|R|, the final |R|-a vertices form a sink suffix. Strongness therefore forces b=0,a=|R|. Consequently

    R -> x -> y -> R,

where R is a transitive chain: the cyclic blow-up C3(1,1,n-2). Conversely this tournament is strong and has exactly two deletable apices for n>=4. Deleting x or y leaves a chain; deleting any chain vertex leaves a directed triangle because at least one chain vertex remains. Hence this is the unique strong ordinary class with two deletable apices at every n>=4. All other strong one-apex classes at those orders have one.

An automorphism fixing a deletable apex fixes the unique order of its transitive deletion and therefore is the identity. Aut(T) acts freely on D(T). Tournament automorphism groups have odd order: an involution would swap some pair and reverse its arc, impossible; Cauchy's theorem rules out even order. Thus every strong one-apex tournament of ordern>=4 is rigid, because its automorphism order divides1 or2. C3 is the separate exception, with three apices and automorphism group of order3.

## 2. Exact ordinary and selfconverse class counts

At ordern>=4 there are2^(n-3) strong words. Each rigid class contributes |D(T)| pointed words, so the unique two-apex class contributes2 and all others1. Therefore the number of ordinary strong one-apex classes is

    c_n=2^(n-3)-1, n>=4;  c_3=1.

Every nontransitive one-apex tournament has exactly one nontrivial strong component. Otherwise deleting one vertex would leave a different cyclic component. The condensation order places transitive singleton source and sink components on its two sides. Consequently an ordinary class of ordern is determined by its strong core and the numbers of leading/trailing singleton vertices. There are n-s+1 choices for a core of orders. Adding the single transitive class yields, for n>=2,

    a_n=1+sum_(s=3)^n (n-s+1)c_s
       =2^(n-1)-binom(n-1,2)-1.

The pointed converse operation on words is reverse-complement:

    w*= (1-b_(m-1))...(1-b0).

For a unique-apex class, selfconverseness is equivalent to w=w*, because any isomorphism to the converse must take the unique deletable apex to itself. Fixed strong words exist only for evenm=2r and number2^(r-1): the endpoints are forced and the other r-1 complementary pairs are free. The exceptional two-apex class C3(1,1,n-2) is selfconverse, but for n>=4 its two words10^(m-1) and1^(m-1)0 are interchanged, not fixed. Therefore the strong selfconverse counts are

    s_3=1;
    s_n=1 for even n>=4;
    s_(2r+1)=2^(r-1)+1 for r>=2.

For a nonstrong class, converse swaps source and sink padding. Selfconverseness therefore requires equal padding and a selfconverse strong core. Including the transitive class, the full fixed-order selfconverse counts are

    t_n=n/2 for even n>=2;
    t_(2r+1)=2^r+r-1 for r>=1.

These are ordinary class counts, not numbers of labelled tournaments, pointed words, paths, or converse-merged classes.

## 3. Exact binary product approximation

Use zero-based indices. Define binary integers

    A=sum_(j=0)^(m-1) bj 2^(m-1-j),
    B=sum_(i=0)^(m-1) (1-bi)2^i.

B represents the reverse-complement word. In AB/2^m, the term from indices j,i is bj(1-bi)2^(i-j-1). Terms with i=j vanish. The forward terms j<i are exactly the ordinary inversion weights in the cut formula. The reverse terms contribute

    R=sum_(i<j) bj(1-bi)2^(i-j-1),  0<=R<m/2.

The upper bound follows by summing the geometric series in i<j for each j; each such inner sum is less than1/2. If e(w) is the number of adjacent10 transitions, the exact cut formula gives

    H(w)=AB/2^m+1+e(w)-R.

Since0<=e(w)<=floor(m/2),

    |H(w)/2^m-(A/2^m)(B/2^m)| <= (m/2+1)2^(-m).       (1)

The full finite binary integers A,B are not independent. The next argument retains that distinction.

## 4. Product limit with a uniform error bound

Choose a strong word uniformly from its2^(m-2) possibilities. Let q=floor(m/2). The first q binary digits of A/2^m and B/2^m use disjoint source bit sets. Their leading digits are1 and their remaining q-1 digits are independent fair bits. Thus their truncations are independent uniform grids of mesh2^-q in[1/2,1).

Couple each grid value to an independent uniform point in its own dyadic cell, producing independent U,V uniform on[1/2,1]. Each exact binary fraction differs from the corresponding coupled variable by at most2^-q, even though its lower source digits may depend on the other prefix. Products in[0,1] differ by at most the sum of coordinate errors. Together with(1),

    |H/2^m-UV| <= 2^(1-q)+(m/2+1)2^(-m) = O(2^(-m/2)).

The product distribution has bounded density, so this gives a uniform-in-z CDF error. Direct integration gives

    F(z)=0, z<1/4;
    F(z)=4z log(4z)-4z+1, 1/4<=z<=1/2;
    F(z)=4z-4z log z-3, 1/2<=z<=1;
    F(z)=1, z>1.

The density on the two interior intervals is4log(4z) and-4log z respectively, bounded by4log2. If F_m(z)=Pr(H/2^m<=z), then

    sup_z |F_m(z)-F(z)|=O(2^(-m/2)).                    (2)

The limiting variance is55/2304, since E(U)=3/4 and E(U²)=7/12. Independently, the finite mean is exactly

    E(H)=9*2^(m-4)+1/2, m>=3.                          (3)

For(3), ordinary inversion terms contribute: mandatory endpoint pair2^(m-2), boundary pairs2^(m-2)-1, internal pairs(2^(m-2)-m+1)/4. The adjacent10 correction has mean(m+1)/4. Adding the constant1 gives(3). Formula(3) is not asserted at m2, where H=3.

## 5. All-order bounded-H word count and cancellation of scale oscillation

Let W(X) be the number of strong binary words of all lengths with H<=X, counting pointed presentations separately. It is finite by the minimum bound. For large X write X=2^k x with1<=x<2. Every word of lengthm<=k satisfies H<2^m<=X. No word of lengthm>=k+3 can satisfy H<=X, since H>=2^(m-2)+1> X. Hence exactly

    W(X)=2^(k-1)-1
         +2^(k-1)F_(k+1)(x/2)+2^k F_(k+2)(x/4).

Here x/2 is in[1/2,1) and x/4 in[1/4,1/2). Substituting the two explicit CDF branches yields the exact identity

    1+F(x/2)+2F(x/4)=2x log2.

Applying uniform error(2), both remaining count errors are O(2^(k/2))=O(sqrtX), uniformly in x. Therefore

    W(X)=(log2)X+O(sqrtX).                             (4)

The logarithm is natural. This is a count of words/pointed strong cores with multiplicity in H, not a count of distinct attained integers. The cancellation identity is why no leading dyadic oscillation remains.

## 6. Ordinary and converse-merged strong-core spectra

Let C(X) count ordinary strong one-apex classes of all orders with H<=X. The unique two-apex class at ordern>=4 has H=2^(n-2)+1. Thus the exact duplicate correction is

    C(X)=W(X)-max(0,floor(log2(X-1))-1), X>1,

where log2 here denotes base-two logarithm. It follows that

    C(X)=(log2)X+O(sqrtX).

Let S(X) count selfconverse strong-core classes with H<=X. All eligible orders satisfy n<=3+floor(log2(X-1)). The fixed-order formulas in Section2 therefore give S(X)=O(sqrtX). The number Q(X) of converse-merged strong-core classes satisfies

    Q(X)=(C(X)+S(X))/2=(log2)X/2+O(sqrtX).             (5)

No all-order nonstrong class count is substituted here: padding would make it infinite.

## 7. Positive density of missing actual counts despite modular completeness

Let A(X) count DISTINCT actual H<=X attained by any one-apex tournament, including nonstrong ones. Source/sink stripping shows that its values are exactly1 together with the values of strong cores. Since H is invariant under converse, Section6 gives

    A(X) <= Q(X)+1 <= (log2)X/2+O(sqrtX).              (6)

There is also a proof of(6) bypassing the apex-classification theorem: reverse-complement WORD orbits provide an upper bound on distinct values; their number is W(X)/2 plus half the number of fixed words. The latter is O(sqrtX) by the same length bound. Thus this strongest obstruction does not depend on class-count bookkeeping.

The number of positive odd integers <=X is X/2+O(1). If B(X) counts the odd integers <=X absent from the one-apex spectrum, then

    B(X) >= (1-log2)X/2-O(sqrtX).

In particular liminf B(X)/X >= (1-log2)/2, and the lower relative density among odd integers is at least1-log2. This is a rigorous lower-density bound; existence or value of a natural density for the attained set is not claimed. Session25 nevertheless proves that every parity-admissible residue is represented, even in a subfamily. Thus profinite/congruence completeness and a positive proportion of missing integers coexist in this mathematically natural structural family.

The explicit value13 is already absent by CUT_ADDENDUM.md. H13_WITNESS.md supplies a five-vertex strong tournament with feedback number exactly2 and actualH13, so the structural obstruction is genuine, not a forbidden value of the full spectrum.

## 8. Optional scoped Dirichlet consequence

Define Z_W(s)=sum_w H(w)^(-s), summing over strong words, and analogously Z_C,Z_Q for ordinary and merged strong-core classes. They converge absolutely for Re(s)>1. Partial summation and(4)--(5) express each as its main coefficient times s/(s-1), plus an integral of its O(sqrtX) counting remainder. That remainder integral is holomorphic on Re(s)>1/2. Hence these series continue meromorphically there, with their only pole at s=1, simple, and residues respectively log2,log2,(log2)/2.

These are presentation/class-weighted spectra. No continuation or residue is asserted for the Dirichlet series over DISTINCT H values. Including all padded classes would not even define these finite coefficients.

## 9. Reproduction, evidence types, and non-consequences

- `node verify_cut.js` reproduces cut/Ferrers, extrema, swaps and independent small Hamilton-DP controls; output `cut-output.json`.
- `node verify_asymptotic.js` exhausts strong words at lengths3..18 for exact finite means and the pointwise dyadic product error. Output `asymptotic-output.json`. Its five-threshold CDF values use floating logarithms and are diagnostics, not exact proof gates.
- Root independently checked class formulas through orders2..7; its exact receipt accompanies the broader session packet. This source does not mislabel the symbolic proof as a finite census.
- The final independent review binds this full file separately. Earlier ASYMPTOTIC_REVIEW.md records this lane's checks of root's derivations.

The binary product, Ferrers shape, apex marker, ordinary class, merged class and actual-H value are different carriers. Every quotient above has its multiplicity or loss explicitly accounted for. No full Hamiltonian spectrum closure, prime-attainment theorem, LRC result, or general merged-metagraph classification follows. The separate WIRING_REVIEW.md audits a genuine feedback-number adjacency constraint; it is not inferred merely from matching spectral counts.
