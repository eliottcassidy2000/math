# The square is a product: two involutions, two ports, and one small core

**PROVED:** a literal product and commuting-involution model for the balanced
two-return family; its variable-fibre growth law; its endpoint response;
and the precise weighted-core and harmonic-face maps below. **FINITE-EXACT:**
independent bounded controls in `microcosm_check.py`. The incoming square
count and its Euclid triple are inherited. The new content retains the
objects and symmetries behind those numbers, and identifies which information
their small quotients lose. No graph-minor reduction or unrestricted
tournament-count recurrence is claimed.

## 1. Inheritance and the object being counted

The closest proved mechanism is
[double_return_counts.md](../strong_blue_games_20261020f/double_return_counts.md):
for a,b>=2, let T_(a,b) be the transitive tournament on 0,...,a+b with just
0--a and a--(a+b) reversed. Its path count is

    H(T_(a,b)) = (2^(a-1)+1)(2^(b-1)+1).                    (1)

The canonical hostile in that note is the unique arc-disjoint replacement
of the standard path: it is always blue in the balanced family, even though
the marked edge has an independent black detour. The recent correction at
the front of `01-canon/MISTAKES.md` explains the first failed implication:
an abstract reordering need not be a legal path or a black move. We keep
actual path words and the positional anti-involution throughout.

The less-used sidecar is
[THM-3509, reduced-fraction harmonic K4 face](../../../01-canon/theorems/THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray.md).
It supplies an exact arithmetic decoder from a three-edge harmonic face
to Euclid's triple, without asserting that a tournament is a minor of that
carrier. MISTAKE-418 warns that raw Fibonacci parameters and primitive
normalization differ; the parameters here are consecutive, so have opposite
parity and need no content division. MISTAKE-225 warns that different sign
actions are not automatically one symmetry. The two involutions below
are explicitly defined on one common set and proved to commute.

The board is: a pointed Boolean set, reflection of actual paths, two endpoint
ports, weighted transitive chains, and the harmonic-face decoder. The anchor
is the user's product/diagonal question; the niche is the endpoint response;
the wildcard is its commuting-involution interpretation of 3-4-5.

## 2. A genuine product with a blue diagonal

Write T_k=T_(k,k), k>=2, A={0,...,k-1}, c=k, B={k+1,...,2k}, and
sigma(v)=2k-v. Let

    X_k = P({1,...,k-1}) disjoint-union {infinity},
    q_k = |X_k| = 2^(k-1)+1.

Interpret infinity as the full set A. For a pair (I,K) in X_k squared,
put J=sigma(A minus K), using this full-set interpretation, and form

    P(I,K) = I increasing, J increasing, c,
             (A minus I) increasing, (B minus J) increasing.       (2)

The incoming independent-subset theorem says (2) is a bijection onto all
directed Hamiltonian paths of T_k. In these coordinates its symmetry becomes
especially simple:

    R(P) = sigma(reverse(P)),
    R(I,K) = (K,I).                                         (3)

Both statements follow by reversing the four increasing runs in (2).
The tournament is rigid: its only repeated scores occur in the two extreme
adjacent pairs, whose internal arcs prevent swapping; all other vertices
have distinct scores. Reflection sigma is its unique anti-automorphism.
Hence a path is blue exactly when R fixes it, exactly when I=K.

Thus the q-by-q array is an actual product, its q diagonal cells are blue,
and its q(q-1) off-diagonal cells are black. Quotienting by R identifies
the black cells in pairs, indexed by unordered distinct pairs of X_k.
This is the edge set of a complete graph on the **choice states**, together
with q diagonal states. Those choice states are not the original tournament
vertices, and this set identification is not a tournament-metagraph map.

## 3. The Pythagorean numbers count specific symmetry orbits

Complement a finite subset inside {1,...,k-1}, and fix infinity. This gives
an involution kappa of X_k with exactly one fixed point. It has

    m=(q+1)/2 total orbits,
    n=(q-1)/2 nontrivial orbits,       m=n+1.                (4)

Now define a second path-set involution

    J(I,K)=(kappa(I),kappa(K)).                             (5)

J commutes with R, preserves blue and black, and has exactly one fixed
path: (infinity,infinity), the standard increasing path. Consequently

    number of nontrivial J-orbits = (H-1)/2 = 2mn = L,
    number of all J-orbits        = (H+1)/2 = m^2+n^2 = C,
    number of R-fixed paths       = q       = m^2-n^2.       (6)

This gives a set-theoretic meaning to every number in the inherited triple
`(q,(H-1)/2,(H+1)/2)`. It is not merely a substitution into q squared.

There is a finer compatibility check. R acts on the set of J-orbits.
A J-orbit is fixed by R exactly when its representative satisfies either
I=K or K=kappa(I). The diagonal and this complement-diagonal each have q
points and intersect only at (infinity,infinity). Their union therefore
has 2q-1 points, of which one is J-fixed and the others form pairs. Thus
there are **exactly q fixed J-orbits**, and

    paired R-orbits on the J-quotient = (C-q)/2 = n^2,
    total joint R,J orbits           = (C+q)/2 = m^2.        (7)

This proves the Klein-four orbit count directly. The same formulas follow
by averaging fixed-point counts `(q^2,q,1,q)`, but no general orbit-counting
theorem is needed for the argument. Equations (6) give

    (m^2-n^2)^2 + (2mn)^2 = (m^2+n^2)^2.

The parameters m,n are consecutive, hence coprime and of opposite parity;
the triple is primitive. At k=2, X has three objects, kappa has two orbits
of which one is nontrivial, and (m,n)=(2,1) gives (3,4,5). The next triples
are (5,12,13), (9,40,41), (17,144,145).

J is a permutation of the path set, not an automorphism of the tournament.
Indeed T_k is rigid and J is nontrivial. Nor have we shown that J is a black
metagraph move, that its orbit quotient preserves path-edge incidence, or
that the triple determines any graph minor. Those predicates need their
own maps and certificates.

## 4. An actual growth map and the exceptional boundary state

There is a surjection X_(k+1)->X_k: forget the last subset bit and send
infinity to infinity. Every finite state has two preimages; infinity has
one. It commutes with kappa. Applied in both coordinates, it is realized
on the **actual tournaments** by deleting vertices k and k+2 from T_(k+1)
and relabelling the surviving centre k+1 as k and later labels by subtracting
two. Formula (2) proves that deletion sends every path to its parent path.

The path fibres have sizes 4,2,1 according as zero, one, or two coordinates
are infinity. The maps commute with both involutions. Writing H=q^2 and
Delta=H-q, this proves

    q' = 2q-1,
    H' = 4H-4q+1,
    Delta' = 4Delta-2q+2.                                 (8)

A black parent has only black lifts. Each finite diagonal parent has four
lifts, two diagonal and two off-diagonal; the exceptional diagonal parent
has one blue lift. Thus exactly 2(q-1) new black paths lie over blue parents.
This explains the correction terms, rather than fitting them to a sequence.

The naive statement "every old path has four lifts" fails already at the
standard path, which has one. The singled-out state infinity is the exact
boundary information needed for the recursive counting law. The map is an
induced-subtournament deletion within this particular family; it is not a
Petersen/Kuratowski minor reduction or a universal tournament growth law.

## 5. The two-port response has a rank-one bulk and a real boundary

Let E_k(s,t) count Hamiltonian paths with start s and end t. For
0<=i,j<=k-1 define

    a_0=1,       a_i=2^(k-1-i) for 1<=i<=k-1.

Then the full rectangular block of endpoint responses satisfies

    E_k(i,2k-j)=a_i a_j.                                  (9)

Proof: these are precisely the pairs in which neither subset coordinate is
empty. Infinity produces endpoint 0 or 2k, and a nonempty subset with
minimum i has 2^(k-1-i) choices. The two choices are independent.
In particular this endpoint matrix block has rank one over the rationals.

For completeness, use formal endpoint markers U_s,V_t and put

    F(U)=U_0 + sum_(i=1)^(k-1) a_i U_i,
    F*(V)=V_(2k) + sum_(i=1)^(k-1) a_i V_(2k-i).

For each nonempty S subset {1,...,k-1}, write alpha(S)=min S and
beta(S)=max(A minus S). The *entire* endpoint polynomial is

    sum_(s,t) E_k(s,t) U_s V_t
      = F(U)F*(V) + U_0 V_k + U_k V_(2k) + U_(k+1) V_(k-1)
        + sum_(S nonempty) [U_alpha V_beta
                           + U_(2k-beta) V_(2k-alpha)].    (10)

The last terms are exactly the cases where one or both choices are empty.
This formula retains the boundary cases that a rank-one-only model loses.
The blue endpoint polynomial is

    U_0 V_(2k) + U_(k+1) V_(k-1)
       + sum_(i=1)^(k-1) a_i U_i V_(2k-i).                 (11)

Hence the endpoint fibre (i,2k-i), 1<=i<=k-1, itself contains a literal
a_i-by-a_i choice array with a_i blue diagonal cells. This is a recursive
product inside the port fibres, not just a repeated numerical ratio.

Endpoints still do not determine blue status. At k=3, endpoints (1,5)
support the following four paths:

    blue:  (1,4,6,3,0,2,5), (1,2,6,3,0,4,5);
    black: (1,6,3,0,2,4,5), (1,2,4,6,3,0,5).

The missing bit is equality of the *complete* two subset choices, not
merely equality of their minima. Similarly, boundary-hub transport in
[transport.md](../strong_blue_games_20261020f/transport.md) retains actual
endpoint pairs and cut edges. Neither H nor this endpoint count matrix
recovers those path-edge incidences without additional data.

The accompanying [quadratic port-response theorem](profile_growth.md)
uses the endpoints of **distinct components of a path cover** to compute
arbitrary one-vertex growth. Formula (10) is a one-component endpoint
response; it does not replace that theorem's two-component tensor B_1.
That distinction is the precise common mechanism: insertion consumes at
most two ports, while computing their coefficients can require the whole
parent. Our rank-one block is a special factorization of this family, not
a universal low-rank law for tournament growth.

## 6. The same five-vertex core reaches the zero-heavy boundary

The canonical maximal transitive chains of T_(a,b) are exactly

    {0}, {1,...,a-1}, {a}, {a+1,...,a+b-1}, {a+b}.          (12)

Their pair-free quotient is T_(2,2), with weights

    (1,a-1,1,b-1,1).                                     (13)

To see maximality, the quotient T_(2,2) is pair-free; each of its repeated
score pairs cannot be homogeneous, and the four possible pairs differing
in score by one have a distinguishing third vertex. Inflation of this
pair-free quotient has precisely the five displayed canonical chains.
This uses the canonical-chain mechanism of
[THM-4463](../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md).

The feedback minimum is exactly two for every a,b>=2. The reference order
has just two backward arcs. The directed triangles (0,1,a) and
(a,a+1,a+b) are edge-disjoint, so any ordering needs a backward arc from
each. This gives the matching lower bound without an enumeration.

At feedback budget 2 the inherited cap is 3. Thus the balanced k=2,3 cases
have zero heavy coordinates, while all k>=4 give the identical capped core
`(1,3,1,3,1)`, with two heavy coordinates and deficiency `2*2+1-2=3`.
This is a literal connection to the boundary of the
[fixed-deficiency grammar](../kernels_games_detours_20261020f/general_fixed_deficiency.md),
not an identification of its class count with H. As k increases, that capped
state remains fixed while H grows exponentially in the true chain lengths.

The loss is stronger than multiplicity. For T_(a,b), its score pattern
forces any anti-automorphism to be label reflection v->a+b-v: all middle
scores are distinct, and the arcs in the two repeated extreme pairs force
their images. Reflection preserves the two return arcs exactly when a=b.
Therefore

    raw blue count = 2^(a-1)+1 if a=b, and 0 otherwise.     (14)

For example T_(5,5) and T_(4,6) both have eleven vertices, feedback two,
the same capped weighted core and the same labelled score vector. Their
blue counts are 17 and 0, and their H values are 289 and 297. The actual
two chain lengths, including their equality, are the necessary sidecar.

The cap was proved sufficient for the feedback-budget predicate in
[sparse_theory.md, §3](../defect_boundary_20261020/sparse_theory.md).
It was never proved sufficient for blue status or for Hamiltonian-path
response. This example explains why a small core can encode unbounded
arithmetic exactly when its full integer weights are retained, and why
the cap cannot replace those weights for a different target predicate.

## 7. Recovering the harmonic face and the precise doubling operation

Here n=2^(k-2) and m=n+1. Under the inherited THM-3509 decoder, the
four-term additive window and its harmonic edge face are

    W=(1,n,n+1,2n+1),
    (u,v,z)=(n,n+1,n(n+1)),       vz=u(v+z).

The triple is recovered linearly:

    (u+v,2z,v+2z-u)=(q,L,C).                              (15)

This is an explicit map from the orbit parameters to a labelled K4
arithmetic carrier. It preserves the integer parameter and Euclid identity;
it forgets the original path, its endpoints, and its detour certificates.
The face equation is the familiar exact relation
`1/n = 1/(n+1) + 1/[n(n+1)]`.

The window satisfies two Fibonacci recurrence steps, but lies on the
unit-Cassini locus only when `|n+1-n^2|=1`, namely n=1 or 2. Thus the
first two windows are genuine Fibonacci windows; later members follow
binary doubling, not the golden-ratio recurrence. The equations
`n^2-n=0` and `n^2-n-2=0` prove the stated intersection exactly.

Doubling the subset dimension sends `(m,n)` to `(m+n,2n)`. Its determinant
is two, so it is not a unimodular Berggren branch operation. It induces
the exact linear triple transport

    [q']   [ 2  1 -1 ] [q]
    [L'] = [-2  2  2 ] [L].                              (16)
    [C']   [-2  1  3 ] [C]

For this matrix M, `M^t diag(1,1,-1) M = 4 diag(1,1,-1)`.
It preserves the Pythagorean cone and multiplies the quadratic form by four
away from the cone. This gives a common, verified growth operation on the
subset-product count and the arithmetic carrier. It supplies no operation
on graph minors or proof of Collatz descent.

## 8. Hostile controls and reproduction

A square H alone does not force q=sqrt(H) blue paths. The five-vertex
tournament with outgoing-neighbour bit rows `(16,1,3,7,14)` has H=9 and
exactly one blue path, `(4,3,2,1,0)`. T_(2,2) has the same H and three blue
paths. The product decoder and its reflection, not squareness alone, justify
the blue diagonal. The same-port and same-cap hostiles above test two
different, stronger quotients.

Run from the repository root:

    python 05-knowledge/results/zero_heavy_ports_20261020g/microcosm_check.py
    python -O 05-knowledge/results/zero_heavy_ports_20261020g/microcosm_check.py

The standalone checker imports no producer code. It compares direct
Hamiltonian DFS with all decoded states at k=2,...,7, tests positional
anti-automorphisms directly, builds all J and Klein-four orbits, checks the
complete endpoint matrices and five actual two-vertex deletion maps, and
recomputes homogeneous chains. It also tests all rectangular a,b=2,...,5,
all 64 weighted five-core cases with chain lengths 1,...,8, the named
same-order hostile, and arithmetic k=2,...,100. All guards are explicit
exceptions. Normal and optimized stdout and JSON agree exactly.

`microcosm_check.json` retains the literal witnesses and counts. The finite
checks supplement the all-order arguments; no finite census is promoted to
an unrestricted theorem. The research move used here is to retain incidence
through counting and compare predicates before applying a quotient.
