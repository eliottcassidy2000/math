# Ferromagnetic digit cuts and a physical source rearrangement theorem

Session199. **PROVED (elementary; independent audit distributed in audit/):** the retained-root
digit interaction is an explicit positive weighted cut. Its least nonzero
energy is known at every depth. More strongly, every complete fixed-partner
C1/C2/C3 Hamilton-path readout has a strict exchange law: an adjacent `OZ`
can be changed to `ZO` to lower the final path count, even before an arbitrary
source suffix including later J births. Hence exact extrema over a block have only linearly many
candidates instead of exponentially many. **OPEN:** prescribing every odd
scalar. Extrema do not supply intermediate values.

## 1. Inheritance and the actual bridge

The physical input is Session198's `bilinear-descent/PROOF.md`, Theorem1 and
(4)-(8), published at `6f3f93a3915ef2a7eba9a611f67a1bb75bd392c9`, and the
Session197 complete two-source C3 theorem. The new digit starting point is
Session198 `source-digit-blocks/PROOF.md`, published at
`527edd68877f50eaa7c0bcb0db847743cc5819d3`; the available reviewed source snapshot
is `session198-source-pair-construction/arithmetic-coverage/DIGIT_BLOCK.md`.
The companion freezes hashes of the exact local proof snapshots used.

The deliberately older sidecar is
`01-canon/theorems/THM-4465-weighted-tournament-shear-production-and-contact-kernel.md`,
Sections3-4: a scalar block state must be lifted to its quadratic contact
kernel when external contacts are nonuniform. That proof separates the
actual weighted quadratic response from an analogy about Euler evolution.
Here the external contacts are **bits in the source construction word**, not
tournament vertices. The exact map is the expansion of a word's Hamilton
response into its quadratic bit coefficients. It preserves the final scalar
and an explicit physical graph recipe. It discards the original graph once
its source profile is taken, so the source recipe remains the sidecar. No
Euler dynamics, spectrum, or weighted-tournament realization is imported.

The second older mechanism is
`01-canon/theorems/THM-3191-factorial-block-exterior-clifford-law-and-global-carry-smith-profile.md`,
Sections2 and6: cancel extension terms only on the actual transfer product,
and retain a selected-coordinate sidecar. Its cheap hostile test here is to
compare a complemented pair's **sum** with its difference. The difference
loses quadratic interaction, but the sum retains a positive cut with a growing
minimum. This is a mechanism transfer, not an identification of the factorial
and source matrices.

Concept board: physical source language; quadratic contact kernel; complement
versus positive sum; order-sensitive exchange; scalar gaps inside extrema.
The nearest hostile is a legitimate block whose odd interval has holes.

## 2. Source conventions and the cut identity

A nontrivial actual start-one source has profile

    s=(D,E,U,d)=(H(V),e_c(V),R_c(V),H(V-c)).

The inherited inequalities give `d>=1`, `E>=0`, `U>=2d`. Chronological source
letters act by

    O(s)=(D,0,2U,d),
    Z(s)=(D-E+U,d,2U-d,d).

Both are literal one-vertex constructions. The starting edge has
`tau=(1,0,2,1)`. The separately available birth
`J(s)=(D-E+U,D-E,D+U,D)` changes the root. We do not exchange J letters.

For a k-letter O/Z block, k>=1, put epsilon_i=1 at Z, otherwise0, and

    A=sum_(i=0)^(k-1) epsilon_i 2^i, B=2^k,
    R=sum_i epsilon_i 2^(k-1-i), p=epsilon_0, t=epsilon_(k-1).

The reviewed digit formula, or direct iteration, gives

    D'=D+UA-Ep-d C(epsilon), E'=dt, U'=BU-dR, d'=d,
    C(epsilon)=sum_(0<=i<j<k) w_ij epsilon_i epsilon_j,
    w_ij=2^(j-i-1)+[j=i+1].                              (1)

All w_ij are positive even integers. If `g` is the one-apex function
`g(1)=1`, `g(2l)=g(l)`, and
`g(2l-1)=g(l)+l-(l mod2)` for l>=2, the other digit expression is

    g(A+1)-1=sum_(i<j) w_ij (1-epsilon_i)epsilon_j.       (2)

Consequently

    K_k(epsilon):=g(A+1)+g(B-A)-2
      =sum_(i<j) w_ij |epsilon_i-epsilon_j|.              (3)

This is the weighted cut on the k bit positions. It is not a cut of the
source tournament. Subtracting the two versions of (2), rather than adding,
gives the inherited complement identity
`g(A+1)-g(B-A)=A-R+t-p`. Thus positive complement pairing preserves K_k;
it does not implement the cancellation in the signed difference.

## 3. The sharp complement interaction gap

For k=1 there is no nonzero cut. For k=2 the least nonzero cut is2.
For every k>=3 the least nonzero value in (3) is

    delta_k = 2^floor((k-1)/2)+2^ceil((k-1)/2).           (4)

For k>=4 equality holds exactly when one cut shore is a single central
position: `i=(k-1)/2` for odd k, or `i=k/2-1,k/2` for even k. For k=3 every
nontrivial cut has weight4.

**Proof.** An interior vertex i has weighted degree
`2^i+2^(k-1-i)`, while either endpoint has degree `2^(k-1)`. The central
singleton cuts therefore attain (4). Consider any other cut for k>=4.
If the two endpoints are separated, their edge alone has weight `2^(k-2)`.
For k>=6 this exceeds delta_k; for k=5 it equals delta_k and another positive
cut edge makes the inequality strict. For k=4 the endpoint edge and the
edge-disjoint paths through each interior vertex give lower bound
`4+min(2,2)+min(2,2)=8>6`.

Otherwise take the nonempty shore S containing neither endpoint. Each
interior i has total edge weight to the endpoints

    b_i=2^(i-1)+2^(k-2-i)+[i=1]+[i=k-2]
        >= (2^i+2^(k-1-i))/2 >= delta_k/2.

If |S|>=2, its endpoint edges already cost at least delta_k. Equality is
impossible for k=4 because the boundary corrections are positive. For k>=5,
equality would force two central positions and no additional cut edge.
At odd k there is only one central position. At even k>=6 there are other
interior positions outside this two-point S, giving an additional positive
cut edge. If |S|=1 its cut is its full degree, minimized exactly at the stated
central positions. The k=2,3 cases are direct. QED.

This all-order growing gap is about the exact complement interaction. It is
not a gap theorem for the union of all Hamilton-path counts.

## 4. Exact scalar cut optimization

Fix an actual partner and a physical integer-coefficient readout linear in the changing source:

    q(s)=alpha D+beta E+gamma U+eta d.

Substitute (1). For fixed k, the actual Hamilton counts are a quadratic
pseudo-Boolean function with pair coefficient `-alpha d w_ij`. Whenever
alpha>0, every mixed second bit difference is strictly negative. Equivalently,
put `c_ij=alpha d w_ij/2` and

    l_i=alpha U2^i-alpha E[i=0]+beta d[i=k-1]
        -gamma d2^(k-1-i)-sum_(j!=i)c_ij.

Then

    q(block(s))=alpha D+gamma BU+eta d
                  +sum_i l_i epsilon_i
                  +sum_(i<j)c_ij |epsilon_i-epsilon_j|.   (5)

This gives an exact integer source-sink cut construction: use undirected
capacity c_ij, terminal capacity max(l_i,0) from i to sink and max(-l_i,0)
from source to i, and add constant
`alpha D+gamma BU+eta d+sum_i min(l_i,0)` to the cut value. Label source-shore
vertices1. The construction has k+2 vertices and O(k^2) edges; no formal
source profile is introduced because every bit assignment is an actual word.
Forced bits can be imposed as terminal constraints. This proves exact
minimum optimization, not a claim that maximum-cut or exact-level feasibility
is solved by min-cut.

## 5. Strict physical rearrangement

**Theorem.** Suppose

    alpha>0, gamma>0, 2alpha+beta+gamma>0.                 (6)

For any nontrivial actual source, any adjacent OZ and any following O/Z
suffix, replacing that OZ by ZO strictly decreases q of the final source.
Consequently among k-letter blocks with exactly r Z letters, the unique
minimum word is Z^r O^(k-r), and the unique maximum word is O^(k-r) Z^r.
Uniqueness is among words, not graph isomorphism types.

**Proof.** Let (D,E,U,d) be the actual state immediately before the exchanged
pair. Direct subtraction gives

    OZ(s)-ZO(s)=(U+E,d,d,0).                              (7)

With no suffix, the readout difference is
`alpha(U+E)+(beta+gamma)d >= (2alpha+beta+gamma)d>0`.
For a nonempty suffix, use its parameters A,p,B in (1), applied linearly to
the difference vector in (7). The final readout difference is exactly

    alpha(U+E)+[alpha(A-p)+gamma B]d >0.                  (8)

Here A-p>=0, and B>=2. The difference vector need not itself be a source;
it is only used in the exact linear transfer. Bubble sorting all inversions
proves the unique extrema at fixed r. QED.

All complete C1/C2/C3 fixed-partner readouts satisfy (6), on either variable
source side. For C1 this is `(alpha,beta,gamma)=(1,-1,1)`. For C2/C3 use the
inherited positive coordinates

    y=(D-E-d,U-D+E-d,E,d-E), y>=0, y2+y3=d>=1.

Their positive bilinear matrices P are stated in the inherited source
proof. With the changing source on the left, let r_i denote row i of P.
Then alpha, gamma and `2alpha+beta+gamma` are respectively

    (r0-r1).y_partner, r1.y_partner,
    (r0+r2-r3).y_partner.                                (9)

For all ten matrices N,M,A,B,C,T,J0,J1,J2,J3 and their transposes,
`r1` and `r0+r2-r3` are entrywise strictly positive; `r0-r1` is nonnegative
and its last two entries are strictly positive. Hence (6) follows on every
nontrivial partner. This is an explicit 4-by-4 arithmetic certificate, printed
by the verifier. The singleton boundaries have the two readouts

    C(o,w)=5D-2E+16U-4d,
    J0(w,o)=T(o,w)=3D-3E+9U-6d,

which also satisfy (6). These are precisely the permitted singleton types.
No theorem about an arbitrary unclassified graft is asserted.

### The same exchange survives arbitrary later root changes

For the classified physical readouts, the following stronger conclusion holds:
the suffix after the exchanged OZ may be **any word in O,Z,J**. It need not
retain the root. This strengthening uses more than condition (6) alone.

Write the readout as r.y. Every fixed-partner row above belongs to the strict
polyhedral cone

    r2>r0>r1>0, r3>0, r0+r2-r3>0.                        (10a)

The extra inequality `r2>r0` follows because row2-row0 of every P and its
transpose is nonnegative, with its last two entries strictly positive.
The earlier conditions supply all other inequalities. For C1 the row is
(2,1,3,3), and the singleton rows are (21,16,36,33),(12,9,15,15).

Precomposing the readout with one physical source operation gives

    r O = (r0+r1,2r1,r0+r1+r3,2r1+r3),
    r Z = (2r0,r0+r1,r0+r2,r0+r2),
    r J = (r2,r0,r2+r3,r0+r2).                            (10b)

Each preserves (10a). In particular the last cone functional becomes,
respectively, `2r0`, `2r0`, and `r2+r3-r0`, all strictly positive. The two
strict differences in (10a) remain either a previous strict difference or r3.
Transporting a final row backwards through any suffix therefore gives another
row satisfying (6). Apply (7) to that transported row to prove the strengthened
exchange.

Consequently, for a source word with prescribed positions of J and prescribed
numbers of O and Z in each intervening block, sorting every block Z-before-O
simultaneously minimizes every classified fixed-partner C1/C2/C3 scalar.
Reverse sorting simultaneously maximizes them. This is not a claim that J can
be exchanged freely with the other letters, nor that sorting preserves H.

## 6. Closed block extrema and a sharp physical order bound

Every mixed sorted block, with r>=1 and s=k-r>=1, has

    D'=D-E+(2^r-1)U-(2^r-2)d,
    E'=0, U'=2^s[2^r(U-d)+d], d'=d.

Thus its q value is

    alpha(D-E-U+2d)+eta d+gamma2^k(U-d)
          +alpha2^r(U-d)+gamma2^(k-r)d.                  (10)

Test the k+1 possible counts r, including separately all O and all Z, to find
the global minimum. For the global maximum, test the k+1 reverse-sorted
words; an interior one has

    q(O^(k-r) Z^r(s))
     =alpha D+(alpha+gamma)2^kU
       +(beta+eta+2alpha+gamma)d
       -alpha2^(k-r)U-(alpha+gamma)2^r d.                (11)

These are exact finite algorithms with only linearly many arithmetic states,
not source membership relaxations. Their minimizers can also be localized by
balancing the two powers of two in (10)/(11).

**Corollary.** Start from tau, append exactly k>=2 retained-root vertices, and
then take the C1 closure. The minimum Hamilton-path count is exactly

    (2^floor(k/2)+1)(2^ceil(k/2)+1).                      (12)

It has an actual witness `tau . Z^floor(k/2) O^ceil(k/2) . J` (and the other
central split when k is odd). This is a strong C1 tournament on k+3 vertices.
For k=2 the all-O and all-Z words tie as well. For k>=3, the balanced split
or splits are the only minimizing words. For k=1 the minimum is5; at k=0 the
closure is the triangle with count3.

Indeed every mixed sorted word gives `(2^r+1)(2^(k-r)+1)`. The endpoints both
give `2^(k+1)+1`, and balancing the two powers minimizes the mixed expression.
The maximum at the same length k>=2 is also exact:

    (2^(floor(k/2)+1)-1)(2^(ceil(k/2)+1)-1)+2.             (13)

It is attained exactly by the balanced reverse splits
`O^ceil(k/2) Z^floor(k/2)` and, for odd k, its other balanced split. Formula
(11) makes this value
`2^(k+2)+3-2^(r+1)-2^(k-r+1)`, maximized by balancing r and k-r.
The endpoints are smaller. The maximum formula also gives the correct
boundary values3 and5 at k=0 and k=1, respectively.

This recovers a real connection to the Mersenne residual in Session198
`balanced-source-normal-form/PROOF.md`, Theorem5: that residual's formula
`(2^(r+1)-1)(2^(s+1)-1)+2` is precisely the *upper rearrangement envelope*
for retained-root words from tau. The lower envelope is its plus-sign product
`(2^r+1)(2^s+1)`. The lower witness ends in O and its final closure has a
balanced root. The upper witness has a nonempty O run followed by a nonempty
Z run, and belongs to the inherited no-balanced-root residual family.
This is a statement about those actual graphs; their scalar images may overlap
other graph families.

In particular every block target has to lie between (12) and (13). These are
not extrema over unrestricted source words containing additional J operations.

## 7. Exact coverage is still a separate problem

With tau, k=4, and the C1 readout, the entire block image is

    {25,27,31,33,37,39,45,47,51}.

The odd targets 29,35,41,43,49 lie strictly inside its extrema and are absent.
They are not claimed absent from the general spectrum. Thus the exact-level
problem cannot be replaced by min/max inequalities, cut positivity, or a
connected bit cube. In particular the rearrangement theorem does not sort a
word while preserving its count.

The direct contribution is a new, physically valid order envelope and exact
extremal compiler for every fixed-partner gate, plus the explicit interaction
gap that any complementary-digit strategy must overcome. The remaining
coverage step must use unsorted words, variable roots or partners, or another
mechanism to fill exact levels. No all-odd or cofinite theorem is claimed.

Run `python3 verify.py` and `python3 -O verify.py`. Both runs use only the
standard library, retain all exception checks, and must give equal receipts.
