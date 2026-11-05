# A physical adjoint and a square-root-height inverse for C1

Session 199. **PROVED:** an exact unimodular pairing transfers every source
word to another actual source by reversing it and exchanging O with Z.
It preserves the number of root changes. The full C1 spectrum is one
positive bilinear image of physical sources. Every represented target n
has a presentation in which **both** source cover counts are at most
floor(2 sqrt(n-1)). The same bound gives an exact target-dependent inverse
for balanced sources with any prescribed bound on root changes.

**OPEN:** this inverse does not prove that every admissible odd target is
accepted. In particular the four-block coverage lemma remains open.

## 1. Physical inputs and the problem with the smaller pairing

Let S be the set of nontrivial physical rooted start-one source profiles

    s=(D,E,U,d)=(H(V),e_c(V),R_c(V),H(V-c)).

The complete physical grammar starts at t=(1,0,2,1) and uses the chronological
letters

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).                (1)

Only J changes the root. Let b=Zt=Jt=(3,1,3,1), the triangle source. Every
source word is physical, and every nontrivial physical response has such a
word. Profiles are merged only with an actual predecessor recipe retained.

The inherited full C1 readout is F(s)=D-E+U: adding a last J vertex produces
a strong tournament with exactly one Hamilton cycle. Conversely, deleting
an outdegree-one vertex from any strong C1 graph gives such a source. The
inherited balanced readout is G(s)=D+2U, realized by appending OJ; the
singleton precursor gives the separate balanced value3.

The ordinary one-entry graft pairing is

    C(v,w)=D_v D_w-E_v E_w+U_v U_w.

O and J are self-adjoint for C, but Z is not. Its exact defect is

    C(v,Zw)-C(Zv,w)
      =D_v(U_w-E_w)-(U_v-E_v)D_w
       +d_v(U_w+E_w)-(U_v+E_v)d_w.

Thus discarding d loses precisely the operation needed to reverse a full
binary block. The pairing below restores that coordinate rather than
assuming a scalar rebase exists.

## 2. Exact physical adjoint

Use column vectors in raw order (D,E,U,d), and put

    Q = [[0,0,0,1],
         [0,1,0,-1],
         [0,0,1,-1],
         [1,-1,-1,1]].                           (2)

Write K(v,w)=v^T Q w. In coordinates,

    K(v,w)=D_v d_w+d_v D_w+E_v E_w+U_v U_w
            -E_v d_w-d_v E_w-U_v d_w-d_v U_w+d_v d_w.

The matrix is symmetric and det(Q)=-1. Direct integer multiplication gives

    O^T Q=Q Z,   Z^T Q=Q O,   J^T Q=Q J,
    Q J=diag(1,-1,1,0).                         (3)

These are identities on all four coordinates, including d. The first three
are not identities for the smaller C pairing.

Define sigma(O)=Z, sigma(Z)=O, sigma(J)=J. For a chronological word w, let
w^dagger be its reversal with sigma applied to each letter. If G_w is its
column-action matrix, then

    G_w^T Q=Q G_(w^dagger),
    K(G_w v,z)=K(v,G_(w^dagger) z).              (4)

This follows by multiplying (3) in the reverse order. The dagger is an
involution and preserves the number of J letters. Every letter on its
right-hand side is an allowed physical move. Consequently (4) is a
physical source-word transfer, not an operation defined only on an
unrestricted integer cone or a signed sum of sources.

The boundary vectors give the exact readouts

    t^T Q=(1,-1,1,0),       b^T Q=(1,0,2,0),
    K(t,s)=F(s),            K(b,s)=G(s).         (5)

## 3. One positive bilinear form is the complete C1 scalar class

**Theorem 1.** The Hamilton-path spectrum of all strong tournaments with
exactly one Hamilton cycle is exactly

    {K(v,w):v,w in S}.                          (6)

Every represented value comes with a literal graph word. If v=G_alpha t
and w=G_beta t, apply the source word beta alpha^dagger to t, then append J.
By (4)-(5), its count is K(v,w), and the final J gives exactly one cycle.
The graph order is |alpha|+|beta|+3.

Conversely every strong C1 count is F(s)=K(t,s), with both t and s in S.
This proves both inclusions. It neither treats all nonnegative coordinates
as physical nor assumes an inverse closure of any response cone.

For balanced counts, let B be the response orbit starting at b under the
same letters. This is a physical subinventory: membership means that a word
from the triangle is supplied. Then the nontrivial-precursor balanced image is

    {K(v,w):v in S, w in B}.                    (7)

Indeed write v=G_alpha t and w=G_gamma b. Equation (4) changes the pair to
K(b,G_(gamma^dagger) G_alpha t). Thus the chronological precursor word is
alpha gamma^dagger, followed by the balanced closing OJ. Conversely G(s)
is obtained using w=b. Adjoin the separate balanced value3.

More precisely, let S_r and B_r consist of responses admitting such words
with at most r J letters, from t and b respectively. Then the balanced
image with at most r root changes in its precursor is exactly

    {3} union union_(i+j<=r) {K(v,w):v in S_i,w in B_j}.       (8)

The triangle seed b itself uses Z, not a charged J, so no hidden unit is
added to the budget. The proof just given retains the sum of the two J
costs. This applies to the four-block class at r=3 without enlarging it.

## 4. Positive coordinates and sharp product bounds

Use the inherited nonnegative physical coordinates

    y=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E),
    D=a+2e+p, E=e, U=a+f+2e+2p, d=e+p>=1.

Substituting into K gives the strictly positive symmetric matrix

    Q_+=[[1,1,2,2],
         [1,1,1,1],
         [2,1,4,3],
         [2,1,3,3]].                            (9)

Let u=(1,1,2,2) and delta=(0,0,1,1), so U=u y and d=delta y.
Both uu^T-Q_+ and 2Q_+-uu^T-2delta delta^T are entrywise nonnegative.
Multiplying by the two nonnegative physical vectors proves

    U_v U_w/2+d_v d_w <= K(v,w) <= U_v U_w.      (10)

In particular, if K(v,w)=n, then

    n <= U_v U_w <= 2(n-1).                     (11)

The lower factor1/2 is sharp: for the transitive sources t_2=t and
t_u=(1,0,u,1), u a power of two, K(t_2,t_u)=u+1 and the ratio to 2u tends
to1/2. The additional dd term is attained in these examples. The upper
bound is attained on every pair of balanced sources: direct substitution
gives K((h,e,h,e),(h',e',h',e'))=hh'=U_v U_w.

These bounds are physical bounds from a positive form. Their direction is
important: the earlier pairing C has a different lower bound; K is not
claimed to be at least U_v U_w.

## 5. A square-root bound on both actual sources

Every physical source letter increases U and increases it by a factor at
most2. For O this is equality; for Z use 2U-d<2U; for J use D+U<=2U.

**Theorem 2.** If n is a strong C1 count, it has a presentation K(v,w)=n
with actual v,w satisfying

    U_v,U_w <= B(n):=floor(2 sqrt(n-1)).         (12)

If n is a balanced value other than3, the same conclusion holds with
v in S and w in B. If a balanced precursor had at most r J letters, the
pair can be chosen with total J cost at most r.

**Proof.** For full C1 begin with n=K(t,G_w t). For balanced counts begin
with n=K(b,G_w t). Split the chronological word as w=alpha beta. By (4)
and symmetry the same n is

    K(G_alpha t, G_(beta^dagger) z),

where z=t in the full C1 case and z=b in the balanced case. As the split
is moved to the right by one letter, the first physical source gains one
letter, and the second loses its last physical letter. Both remain actual,
the bilinear value remains n, and the total J count remains unchanged.

Set T=sqrt(n-1). If the first source already has U>=T at the initial split,
use that split. Its height2 is at most2T because n>=3, and the other height
is at most2(n-1)/U<=2T by (11). Otherwise use the first split with first
height at least T, if one occurs. The preceding height was below T, so
the new height is less than2T. Again (11) bounds the other by2T.

If no crossing occurs, use the terminal split. The first height is below T,
and the second is the seed: height2 for full C1, height3 for balanced.
The latter case has n>=5, so3<=2T. Thus both heights are at most2T in
every case. Their integrality proves (12). All costs were retained along
the split, proving the budget assertion. QED.

The theorem gives two genuinely smaller physical inventories, not merely
the assertion that the smaller member of an arbitrary pair is small. The
source words also reconstruct the original target exactly after reversal.

## 6. An exact target-dependent inverse

For a given n, enumerate the physical orbit S through U<=B(n), retaining a
predecessor word. Test K on pairs and accept exactly when one equals n.
This is a sound and complete inverse for full C1 by Theorems1-2.

For balanced membership, also enumerate the orbit B starting at b through
the same cutoff, test the cross pairs, and allow3 separately. For the
bounded-root-change version, retain the smallest J cost for each profile
in each of the two inventories, and test only pairs with total cost<=r.
Equation (8) and Theorem2 make this an exact inverse for that restricted
class. Arbitrarily long binary blocks are included, not capped by a guessed
word length.

The enumeration is finite and complete because U strictly increases and
all coordinates are bounded by U. A priority queue ordered by U computes
the minimum J costs: every predecessor is processed before its successor.
Every accepted pair yields the explicit physical word in Section3. These
are algorithms polynomial in a crude numerical bound, not a polynomial
input-bit-time or universal-existence theorem.

The four-block coverage conjecture is now equivalently an assertion that,
for every odd n>27,781, the two inventories capped at B(n) contain a
K-pair of total J cost at most3. The source membership, cutoff, budget and
graph reconstruction are proved. The assertion that such a pair exists for
every target remains OPEN.

## 7. Verification and historical connection

The verifier checks det(Q), all adjoint identities, QJ=C, both readouts,
the positive coordinate substitution, and both product inequalities by
exact integer matrix calculations. It checks physical dagger transfers,
midpoint cutoffs and J budgets on source words, and counts literal C1 and
balanced graph constructions on small independent source pairs.

It also compares the square-root pair inverse at target cap5,001 against
a direct complete forward source bank. Full C1, balanced, and balanced
root-change budgets0,1,2,3 are compared separately. Run normal and optimized
Python serially; their receipts must agree. The optional `--target n` and
`--budget r` modes return exact physical recipes rather than a residue match.

The historical two-metric lesson was to retain both arithmetic continuation
data and an ordinary-height control. Here these become a unimodular exact
adjoint on all four response coordinates and a physical midpoint bound.
The hidden-structure lane's four binary blocks supplied the constructive
root-change budget. No Collatz, Mahler or continued-fraction conjecture is
used. The all-order result is an actual inverse and graph compiler; the
universal odd-coverage conclusion is not claimed.
