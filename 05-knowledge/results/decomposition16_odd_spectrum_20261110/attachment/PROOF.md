# The complete one-birth language over the H=3 sources

**PROVED elementary refinement of inherited mechanisms; FINITE-EXACT
implementation controls. The unrestricted odd spectrum remains OPEN.**

This note supplies an explicit arithmetic compiler for every one-vertex
extension of every tournament with Hamiltonian-path count 3. It also gives
the full cut histogram, rather than only its mean. The family is a useful
controlled boundary of the proposed induction: it cannot supply a child
with count 11 at any order.

## Inheritance and scope

The complete H=3 source classification is inherited:

    T_ab = TT_a -> C3 -> TT_b,   a,b >= 0.

See [Session207, sections 1 and 6](../../session207_degree_three_metagraph_20261106/PROOF.md).
The odd-cycle collection formula is [THM-002](../../../../01-canon/theorems/THM-002-ocf.md).
The one-apex function below and its recurrence are inherited from
[THM-4875 (D)](../../../../01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md).
The full cut quadratic, its positivity and lattice are already
THM-4114/4118; the mean and the sharp order-five counterexample are in
[the prior growth obstruction](../../tournament_joint_parity_20261102/growth_obstruction.md).
The new object here is the explicit two-integer normal form and target
inverse for this entire all-order family. No external priority claim is made.

The source/target map retains the old vertex names, both flanks, the named
triangle, and the new vertex. Replacing these by its output count loses
the next attachment response. In particular, this is not an induction
state for unrestricted subsequent births.

## A1. Two binary integers and the exact count

Adjoin a vertex x. A cut bit is 1 precisely when x beats that old vertex.
Read the left flank in its dominance order as a binary word, with integer
value A. Read the **complement of the right flank in reverse order** as a
binary word, with integer value B. Leading zeroes are retained when the
original flank lengths or names matter. Write alpha=A mod 2 and beta=B mod 2;
for an empty flank the corresponding integer and bit are zero.

Let g(w) be the Hamiltonian count of a transitive chain with a new apex
having cut word w. It depends only on the integer represented by w, and

    g(0)=g(1)=1,
    g(2q+1)=g(q),
    g(2q)=g(q)+q+(q mod 2)  (q>=1).

Let k be the number of triangle vertices beaten by x. Rotating C3 shows
that k determines its cut up to an automorphism of the named triangle.
When k=1 or 2, put

    C(A,B) = 12AB+6A+3B+3g(A)+3g(B)+2alpha+beta-1.       (A1)

Then

    H(T_ab+x) = C(A,B)  if k=1,
    H(T_ab+x) = C(B,A)  if k=2.                         (A2)

For k=0 or 3 the triangle remains a module. Define

    U(A,B)=12AB+9A+3g(A)+3g(B)+3alpha-3.

Then

    H(T_ab+x)=U(A,B) if k=0, and U(B,A) if k=3.           (A3)

In particular these values are divisible by 3. A tempting stronger
claim, that one can collapse the triangle and multiply its path count
by 3, is false: paths can leave and reenter the module through the two
flanks and x. The proof and the audit hostile below retain that reentry.

### Proof of the formulas

By the odd-cycle formula, old-cycle collections contribute 3. Every other
collection consists of one odd cycle through x, possibly together with the
old triangle. Removing x from its cycle leaves an even-order directed path
in T_ab starting in the out-set of x and ending in its in-set. If this path
avoids the triangle its residual weight is 3; otherwise the residual
weight is 1. Thus, if F_0(z) counts such paths avoiding the triangle by
their number of vertices and F_1(z) counts those meeting it,

    H = 3 + 3(F_0(1)+F_0(-1)) + F_1(1)+F_1(-1).         (A4)

For a binary chain word w let

    f(w)=sum_(i<j) w_i(1-w_j) 2^(j-i-1),
    d(w)=sum_i w_i(1-w_(i+1)).

Its within-chain path polynomial takes values f(w) at 1 and d(w) at -1,
so g(w)=1+f(w)+d(w). Complementing and reversing a word preserves g.
Let L(z) count paths beginning in the left out-set and allowed to stop
anywhere in the left chain; R(z) counts paths ending in the right in-set
and allowed to begin anywhere in the right chain. Then

    L(1)=A, L(-1)=-alpha, R(1)=B, R(-1)=-beta.

Consequently

    F_0(1)+F_0(-1)=g(A)+g(B)-2+AB+alpha beta.

The triangle has 3 paths of each of lengths 1,2,3. There are k of each
length beginning in its out-set, 3-k ending in its in-set, and, for k=1,2,
one out-to-in path of length 2 and one of length 3. For k=0,3 there are
no such internal out-to-in paths. Put delta=1 for k=1,2 and delta=0
otherwise. Hence

    F_1(z)=3(z+z^2+z^3)L(z)R(z)
          +(3-k)(z+z^2+z^3)L(z)
          +k(z+z^2+z^3)R(z)+delta(z^2+z^3).

Evaluating at 1 and -1 gives

    F_1(1)+F_1(-1)=9AB-3alpha beta
          +3(3-k)A+3kB+(3-k)alpha+k beta+2delta.

Substitution into (A4) cancels the alpha beta term and yields (A1)-(A3).
This proof identifies the mixed term 12AB and the parity corrections;
neither can be inferred from H(T_ab)=3 alone.

## A2. Strong cores, a terminating inverse, and marked classes

For k=1, strip initial zero bits from the left word and final one bits from
the right word. The removed vertices are respectively an outer source
chain and an outer sink chain. The remaining graph, denoted G(A,B), is
strong and has

    4+bitlength(A)+bitlength(B)

vertices, where bitlength(0)=0. The triangle together with x is a strong
four-vertex core. If the left flank is nonempty, x beats its first vertex;
if the right flank is nonempty, its last vertex beats x. These arcs and
the forward flank arcs make every remaining vertex mutually reachable.
The k=2 graphs are the converses of these graphs after swapping the flank
coordinates. All actual one-birth graphs with mixed triangle contact are
therefore obtained from these canonical cores by outer transitive padding.

Since g is positive,

    C(A,B)>=12AB+6A+3B+5.                               (A5)

To decide whether h is in the mixed image, it suffices to test

    0<=A<=floor((h-5)/6),
    0<=B<=floor((h-5-6A)/(12A+3)).                       (A6)

For a successful pair the binary words are a concrete tournament
certificate, not merely a congruence or formal response vector. This is
an all-order target inverse for this restricted family. It need not find
a pair for an arbitrary permissible odd h.

The uniform image also has a terminating target inverse. U(0,B)=3g(B).
Any one-apex value r has an even binary representative 2q, after repeatedly
using g(2q+1)=g(q), and g(2q)>=q+1; thus q<=r-1 suffices. For A>=1,
U(A,B)>=12AB+9A+3, so it suffices to test

    A<=floor((h-3)/9), B<=floor((h-3-9A)/(12A)).

Together these inequalities decide the **entire all-order** H=3 one-birth
image at a requested h. The bounds truncate parameters by the target,
not by a speculative maximum source order.

For h not divisible by 3, every one-birth extension from an H=3 source is
mixed. Let M_h(z) count isomorphism classes of pairs (child,x), with x
marked and H(child-x)=3, by the number of vertices. Then exactly

    M_h(z) = 2/(1-z)^2 * sum_(A,B>=0:C(A,B)=h)
                             z^(4+bitlength(A)+bitlength(B)).       (A7)

To see uniqueness, delete the mark: its unique triangle and its two
ordered transitive flanks recover both words. Triangle rotations identify
the three cuts of the same cardinality. The two cardinalities k=1,2
cannot be identified by a mark-preserving isomorphism; the mark has
different out-degree into the intrinsic triangle. The factors (1-z)^-2
record the two independent outer padding lengths. Unmarking x requires
a further quotient and is not asserted in (A7).

For example C(1,0)=13 is the only solution of C(A,B)=13, so

    M_13(z)=2z^5/(1-z)^2.

This agrees with the separate SCC/deletion-core census in this packet.

## A3. The exact cut histogram and a reusable hostile

For fixed a,b, enumerate A=0..2^a-1 and B=0..2^b-1, retaining fixed-width
words. The complete
histogram of all 2^(a+b+3) equally weighted cuts is the polynomial

    sum_(A,B) [t^U(A,B)+t^U(B,A)
                         +3t^C(A,B)+3t^C(B,A)].          (A8)

The multiplicities 1,3,3,1 count actual triangle cuts; they are not
unlabelled child multiplicities. This formula recovers the order-five
histograms of TT_2->C3 and TT_1->C3->TT_1, including their identical mean
27/2 but distinct variances 279/4 and 351/4. It explains the lost
coordinate by the left/right integers and their product.

**No H=3 source at any order gives a child of count 11.** Uniform contact
would make the child count divisible by 3. In (A6), h=11 leaves only
(A,B)=(0,0),(0,1),(0,2),(1,0), with actual counts 5,9,17,13. This proves
the absence without an order cutoff. The same inverse checks any other
target exactly. A failure here is a limitation of the H=3 source lane,
not a new forbidden tournament count: 11 is attained by other sources.

**The first extra obstruction after combining H=1 and H=3 sources is 55.**
Applying the two terminating inverses above and the one-apex inverse to
the odds through 55 gives misses exactly 7,21,55. This is a FINITE-EXACT
target calculation with an all-order completeness proof, not a census
with a vertex cutoff. In particular every tournament of count 55 obeys

    H(T-v)>=5 for every vertex v.                       (A9)

The threshold is sharp. The seven-vertex tournament with rows
(54,124,56,113,96,64,5) has H=55 and deletion counts
(9,25,15,19,29,29,5). The independent deletion-core lane supplied this
witness from the inherited complete order-seven bank; this checker
recounts the graph and all deletions, also by literal permutation sums.
This is a concrete reason an induction cannot stay within the first two
parent-count lanes. It does not make 55 a globally forbidden value.

### Audit repair: a uniform module still permits reentry

The first draft incorrectly replaced U by three times the collapsed
one-apex count. The independent history audit supplied T_11 with cut
bits 11110 (x beats the left singleton and the whole triangle, and loses
to the right singleton). Its child rows are

    (30,20,24,18,32,15), H=27.

The collapsed quotient rows (6,4,8,3) have H=5. The Hamilton path
(1,2,4,5,0,3) leaves and reenters the triangle {1,2,3}; x need not be
adjacent to the triangle at both visits. The first failed implication
was "uniform module contact forces one visit". The corrected U formula
comes from the same path-polynomial proof without dropping the visits.
For k=0, if w_0 is the collapsed word, exactly

    U(A,B)-3g(w_0)=6A(B+1),

and the k=3 case swaps A,B. Divisibility by 3 survives, so the mixed
normal form, its finite inverse, and the no-11 conclusion are unaffected.
This repair is recorded before publication; no global theorem was
promoted from the faulty draft.

## Verification and next obligation

Run `python3 -B check.py` in this directory, and compare `python3 -O -B
check.py`. The checker imports no repository code. It constructs every
cut on every T_ab with a+b<=4 (1,032 actual children), compares literal
subset-DP counts to (A1)-(A3), and uses permutation enumeration on every
child of order at most 6. It separately tests 288 canonical mixed graphs
with 0<=A,B<12 and both contact cardinalities; checks all 370 mixed and
184 positive-A uniform target candidates through h=501 against direct
graph counts; and records exact low target supports.
The small support table is diagnostic, not a new global coverage record.

The next constructive obligation is a replacement rule that moves from
this completely described source family to richer actual source cores
when (A6) has no solution, while carrying their joint attachment data.
The formula gives a terminating rejection and a concrete new state when
successful; it supplies neither unlimited reentry nor the full spectrum.
