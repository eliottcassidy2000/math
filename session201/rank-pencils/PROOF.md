# Synchronized source pencils: prime exclusion and four physical probes

Session 201. **PROVED, independently audited in the stated scope.** Every rank-one
synchronized source boundary forces a composite target, without any return-word
hypothesis. Every rank-two boundary is separated by a fixed bank of four
physical sources of cover height at most six. This gives a complete
rank-pruned inverse for prime C1 targets, not a prime-existence theorem.
The full `{7,21}` spectrum remains **OPEN**.

## 1. Inheritance and the imported mechanism

The complete physical source grammar is the orbit S of t=(1,0,2,1) under
chronological letters

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

Use b=(3,1,3,1), z=(1,0,2,0) and

    Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]],
    K(v,w)=v^T Q w,
    G(s)=D+2U, H(s)=2U-d.

The synchronized orbit R consists of pairs (w,q)=(M_beta b,M_beta z),
with the SAME physical word beta. Here w is a source, while q is an
auxiliary boundary, not necessarily a source. For every v in S,

    m=K(v,w), h=K(v,q)

is the jointly physical pair (G(s),H(s)) of the source
s=t.alpha beta-dagger, where v=t.alpha and dagger reverses the word and
exchanges O/Z. Hence

    3<=h<m<2h, and m,h are odd.                         (1)

The balanced closure realizes m, and for k>=1 the source
s.O.Z^k followed by J realizes the C1 target

    n=m+L_k h, L_k=2(2^k-1).                           (2)

The earlier two-return theorem proves rank collapse for a particular
class of words. This note removes the return-word restriction and also
characterizes the surviving rank-two case. The historical normalization
hostile remains load-bearing: a primitive integral quotient of a physical
profile need not itself be physical.

The attached manuscript, C. Soule, *On the abc conjecture*,
arXiv:2611.05217v1, uses two independent sections and excludes a finite
exceptional fiber set before comparing heights. We import the question
"when do two retained measurements actually remain independent?" into
the synchronized source pencil. There is no identification of tournament
sources with line bundles or adelic lattices, and no use of the manuscript's
abc, Mordell, or Siegel-zero conclusions. See SOURCE-CONTEXT.md.

Concept board: physical joint image; rank of two boundary columns;
primitive-but-nonphysical quotient; prime target; four-source probe bank;
absorbing exceptional locus. The mathematical map below preserves exact
joint readouts, not only a scalar count or a projective direction.

## 2. Rank-one boundaries force a nontrivial integer product

In positive coordinates

    y=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E),

every physical source has y>=0 and d=e+p>=1. The bilinear matrix is

    Q_+=[[1,1,2,2],[1,1,1,1],[2,1,4,3],[2,1,3,3]].       (3)

All its entries are nonnegative, and its bottom-right two-by-two block
is entrywise at least three. Thus for any two integral vectors v,p in
this cone whose deletion coordinates are positive,

    K(v,p)>=3 d_v d_p>=3.                             (4)

The second vector in (4) need not be physical. This is an inequality on
an explicitly named cone, not an enlargement of S.

**Theorem 1 (general collapse obstruction).** Let (w,q) be in R and assume
the two columns have rank one over Q. Then there exist a primitive
integral vector p and positive integers g,a such that

    w=g p, q=a p, g odd, g>=3, g/2<a<g.              (5)

For every v in S, both m=K(v,w) and every n in (2) are composite.
More precisely,

    m=g K(v,p), n=(g+L_k a)K(v,p), K(v,p)>=3.        (6)

**Proof.** Let g be the gcd of all coordinates of w and put p=w/g.
The positive source coordinates are integral linear combinations of the
raw coordinates, so p remains in the positive cone and d_p>=1. Since
D_w is odd, g is odd. Rank one gives q=c p for a rational c. As p is
primitive, Bezout's identity implies c is an integer a. Apply (1) with
v=t. Equation (4) makes K(t,p)>0, so (1) gives 1/2<a/g<1.
In particular a>0, g cannot equal one, and odd g>=3. Equations (6)
follow by bilinearity and (4); both factors in each product exceed one.
QED.

No physical recipe for p was claimed or used. Dividing w by g is only an
integer factorization certificate. Treating p as a new physical source
would repeat the primitive-normalization error.

**Corollary 2 (prefix exclusion for primes).** A synchronized word used
in any prime target via m or (2) has rank two at EVERY prefix.

**Proof.** Multiplying both columns by the same generator cannot increase
their rank. Both remain nonzero on R by (1). Thus rank one is absorbing,
and Theorem 1 excludes a prime output at that prefix and every extension.
QED.

This strictly generalizes the two-return obstruction: no requirement that
beta return b to a scalar multiple is present. For example beta=OJOOJ has

    w=(45,9,45,9)=9(5,1,5,1),
    q=(25,5,25,5)=5(5,1,5,1),

which is rank one but w is not on the triangle ray. With v=t,k=1,
m=81,h=45,n=171=19*9. The prior two-return statement does not directly
apply to this non-return output, while Theorem 1 does.

## 3. Four small physical probes detect every surviving pencil

Fix the four actual sources

    v0=t=(1,0,2,1),
    v1=t.O=(1,0,4,1),
    v2=t.Z=(3,1,3,1),
    v3=t.Z.J=(5,2,6,3).                               (7)

Their cover heights are 2,4,3,6. The row matrix P of these vectors has
determinant -4. The response matrix is

    P Q=[[1,-1,1,0],
         [1,-1,3,-2],
         [1,0,2,0],
         [3,-1,3,0]],                                (8)

which is invertible. Define m_i=K(v_i,w), h_i=K(v_i,q) and

    Delta_i=m_0 h_i-m_i h_0, i=1,2,3.                 (9)

**Theorem 3 (physical separation and its quantitative gap).** For every
(w,q) in R the following are equivalent:

1. (w,q) has rank one;
2. Delta_1=Delta_2=Delta_3=0;
3. h/m is independent of the physical test source v in S.

If the rank is two, one of the three probes in (9) has Delta_i nonzero,
and its ratio differs from the seed ratio by at least

    |h_i/m_i-h_0/m_0| >= 1/(6 U_w^2).                (10)

**Proof.** Rank one gives (3) by proportionality and therefore (2).
Conversely m_0>0. If all Delta_i vanish, the four-dimensional vectors
(h_i) and (m_i) are proportional. Invertibility of (8) makes q and w
proportional, proving rank one. This proves the equivalence.

In rank two a nonzero Delta_i exists. All m_i,h_i are odd by (1), so
Delta_i is a nonzero even integer and |Delta_i|>=2. The inherited
physical product bound K(v,w)<=U_v U_w gives
m_0<=2U_w and m_i<=6U_w. Divide |Delta_i| by m_0 m_i to obtain (10).
QED.

Thus a single fixed small physical bank detects every loss of independence,
at any boundary height. This does not mean that four sources suffice to
hit every prescribed target. Separation of a rational slope is not
occupancy of an integer target fiber.

## 4. A complete rank-pruned inverse for prime targets

The synchronized inverse proves that a prescribed joint pair (m,h) is
physical exactly when it can be represented by v and (w,q) with both
physical cover heights at most floor(2 sqrt(m-1)). Root-change budgets
are retained by the common word and the sum of its two costs.

For a prime n>=5, use that inverse for balanced n, and use its Mersenne
branches for k>=1 with

    2^k<=(n+1)/6,
    U_v,U_w<=floor(2 sqrt(floor(n/2^k)-1)).           (11)

**Theorem 4.** At each synchronized-orbit generation step, discard any
state of rank one and do not extend it. The resulting inverse accepts
exactly the same prime targets as the unpruned inverse, with the same
minimum total J costs. The exceptional target 3 is handled separately.

**Proof.** Every discarded state and every descendant has rank one by
Corollary 2, and cannot produce any prime output. Thus no prime witness,
including a minimum-cost witness, is removed. Every retained witness is
still a literal common-word construction. The inherited cutoffs prove
the completeness and finiteness of the remaining search. QED.

One can test rank directly using the six two-by-two minors, or using
the three physical determinants in (9). The latter does not require
adding any nonphysical source to the bank.

The guard is necessary, not sufficient. The seed pair (b,z) has rank two,
but v=t.O and k=2 give m=9,h=7 and n=9+6*7=51=3*17.

## 5. Verification and scope

`verify.py` checks the exact basis determinant, positive-cone inequality,
rank persistence, the primitive factorization, and the quantitative probe
gap. It compares a complete height-bounded synchronized census with and
without the rank guard, and checks prime C1/balanced images and minimum
J costs against a complete direct source inventory. Literal graph counts
include a rank-one composite, a rank-two prime, and a rank-two composite.
Normal and optimized Python must give identical receipts.

These finite checks verify implementation and examples. Theorems 1-4
are all-order consequences of the proofs above. No all-order prime
coverage, abc theorem, or full `{7,21}` theorem is asserted.
