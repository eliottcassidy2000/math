# Port selection, exact factor erasure, and guarded balanced renewal

**PROVED source-preserving selector laws and an infinite branching renewal
family. FINITE-EXACT literal controls and two explicit 7p witnesses. This
family does not close the Hamilton-count spectrum.** The separate
[balanced-architecture bound](BALANCED_BOUNDARY.md) gives a stronger reason
that these renewals cannot supply all targets.

## 1. Recovery and the question that changes

The nearest proved mechanism is the [oldest-leaf reset N7](../two_port_return_induction_20261103/NATIVE_RESET_INDUCTION.md).
It erased a continuation factor three while retaining two physical names.
The hostile is its failed latest-two-leaf extrapolation: penultimate is
insufficient once an older leaf remains. The corrected near miss is treating
primitive continuation content as coverage. The least-used coordinate is
the position of a physical port inside an actual homogeneous transitive module.

We recover the [full marked block law B1](../adaptive_state_induction_20261103/BLOCK_SWITCH.md),
[universal sealing map and exact divisor R1/R3](../two_port_return_induction_20261103/ROOTED_RESETS.md),
and [THM-4145, rooted homogeneous-pair expansion, §5](../../../01-canon/theorems/THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula.md).
The latter specifically warns that child responses require a mixed two-ear
packet. The twin derivative below supplies one exact slice of that packet.
General clone counting is inherited, including the [S155 clone-slot law](../session155_spectrum_20261103/primes/README.md).

The live board has four objects: ordered clone ports and their next-birth
response; the odd gcd of all seventeen continuation queries; actual marked
reentry after a stabilized native program; and balanced rooted states on
which reset cycles can be repeated. The anchor is repeated actual source
selection; the niche is its exact factor criterion; the wildcard is the
small projective controller governing which resets remain primitive.

Write a rooted profile as z=(h,s,q,F), where h=H(T), s counts paths starting
at the root, q=H(T−root), and F counts ordered two-slot covers whose first
nonempty path ends at the root, with the second path possibly empty. Put
a=h−s and d=F−q. A birth at r adds x with r→x→T−r. The common seventeen-query
vector is X and its odd content is D=gcd(X). We use the static basis; its
dyadic equivalence to the adaptive basis is already certified in R3.

No relabeling is silently used as a physical action. Each reset below changes
explicit old arcs while fixing the retained names. At its endpoint a new
chronological source certificate is supplied. Native, temporary-birth,
deletion, and permanent-vertex clocks are separate.

## 2. An exact derivative for moving a port through a transitive module

Let u→v be a homogeneous pair in an arbitrary tournament T. Contract it to
a marked vertex w in Q. Then

    s_u(T)−s_v(T)=s_w(Q),
    F_v(T)−F_u(T)=F_w(Q),
    H(B_v(T))−H(B_u(T))=F_w(Q)+s_w(Q)>0.             (2.1)

**Proof.** Swap the two clone labels in a Hamilton path starting at v.
Their external incidences agree, and the prohibited adjacency v→u never
occurs. This injects into the paths starting at u. The omitted paths are
precisely those beginning u,v; contraction counts them by s_w(Q).

Similarly swap the labels in an ordered two-slot cover whose first path
ends at u. This injects into the covers whose first path ends at v. The
only omitted covers have their first path ending with the adjacent segment
u,v. Contracting that segment gives exactly F_w(Q). The empty second slot
is included. Subtract the two birth laws h+F−s to obtain (2.1).
Finally F_w(Q)≥H(Q−w)≥1, including the singleton convention, so the increment
is always strict. No strongness assumption is needed.

Every adjacent pair in a homogeneous transitive module is homogeneous in
the whole graph. Thus, with the graph fixed, the next-birth count strictly
increases as the chosen port moves towards the sink end of that module.
One next-birth observation distinguishes all its positions. This is a
rooted mixed-ear derivative, not a statement that the intermediate integer
counts all occur. Root and spectrum independently audited the two injections.

The probe checks both identities and the birth increment on all 285 rooted
labelled quotients of orders one through four. Literal permutation covers
independently check the rooted DP profiles on that universe.

## 3. Position determines exactly which factor can be removed

From an arbitrary nonempty rooted source (T,v), add k children at v. They
form a homogeneous transitive block. Retain v and one named child w; suppose
A children precede w in the directed order and B follow it, so k=A+B+1.
These are newer and older positions respectively. Every permutation of the
block has the same ordinary graph and H; only the retained port's position
changes the pointed source.

**Exact content theorem.** For every A,B≥0,

    D(T_with_block;v,w) = gcd(h,s,q,F)       if B=0,
                         gcd(q,F,h−s)      if B≥1.    (3.1)

Thus the oldest position has the smallest possible content, and all other
positions have exactly the same content. The first gcd divides the second.
Factor erasure is strict exactly when these two actual source gcds differ.

**Interface guard.** This theorem concerns the existing pair (v,w). The
concurrent [AGE A1/A4](../source_reset_branches_20261103/AGE.md) instead
adds a child at w, discards v, and studies the new pair (w,newborn). Its
divisor is the rooted-profile gcd at w and can depend on k. For the supplied
actual profile (179,10,41,566), our oldest pair has divisor1 at every k,
whereas resealing at that leaf gives divisor3 at k=2. The old root's name
is therefore genuine attachment data. The balanced family below satisfies
the stronger resealing guard by its explicit formulas.

The physical construction is: perform B fixed-root births, seal at the next
child w, then perform A more births at v while retaining w. The key point is
that the latter A births lose no odd divisor of the sealed source profile.
Here is the all-order inverse establishing that claim, contributed and
audited by the spectrum lane.

For an input profile (h,s,q,F), put n=A+1. After its seal and A further
fixed-root births, the four observable numbers (H,q_v,F_v,F_w) satisfy

    q_v=q,
    F_v=2^(A+1)d+q,
    H=a+q+(2^(A+1)−1)d,
    F_w=2^A a+(2·3^A−2^A)d+(2^(A+1)−1)q+s.        (3.2)

Recover q, then d, then a, then s in that order. The only divisions are by
powers of two. The four-observation determinant on the nonnegative profile
coordinates (a,q−s,s,F−h−q+s) has absolute value 2^(A+1). Each of these
four observations is itself a dyadic linear combination of the seventeen
queries, by the prior observer identities. The forward map also has only
dyadic denominators. Consequently its odd query gcd is exactly the input
profile gcd, for every A including A=0.

For B=0 this proves the first line of (3.1). For B≥1 the input to that seal is

    R^B z=(a+q+(2^B−1)d, q, q, 2^B d+q).

Its gcd is gcd(a,q,d)=gcd(h−s,q,F), because q is odd and powers of two are
units modulo every divisor in question. This proves the second line.

These are universal observer identities for actual rooted sources, including
weak ones; they do not classify which sources belong to the C3 program
language. When an initial source program retaining v is supplied, marked
reentry is constructive: use the desired reverse block order as the birth
order, retain w when it is born, and subsequently discard each new child.
That preserves both terminal names. The renewals below start at C3 and
supply complete named programs, so their entry obligation is discharged.

For the old four-vertex rooted source z=(5,2,3,6), the two gcds are 1 and 3.
For the actual balanced source z=(9,3,3,9), both are 3: no block permutation
can erase that inherited factor. This is a sharp hostile to “every module
port reset makes the source primitive.” The probe checks 196 actual
source/block/position query vectors, including these and a regular five-core.

## 4. Select one port with a linear native clock

Full reversal of a k-leaf module changes binomial(k,2) arcs. To move one
retained physical port through d positions, change only the d adjacent arcs
it crosses, leaving every other relative order unchanged. Each adjacent
swap keeps the block transitive and gives an isomorphic whole tournament;
it preserves H and |Aut|, hence also the ordinary source weight H/|Aut|.

There is a literal traced Hamilton path for each swap. Contract the module,
find a Hamilton path of the quotient by ordinary insertion, and expand its
marked vertex into the block's directed order. The selected adjacent arc
u→v appears consecutively. Apply the inherited
[one-temporary-vertex compiler, METAGRAPH §5](../adaptive_metagraph_reformulation_20261103/METAGRAPH.md):

* If v is not the last vertex of that path P, create τ beating only u,
  use Pτ, then use P with τ inserted between u and v.
* If v is last, create τ losing only to v, use τP, then the same inserted
  path. The predecessor u is not first because our sources have order at
  least three.

Both words are legal on their actual inputs. After the two native moves,
delete τ; only the old arc uv has reversed. The probe retains every temporary
cut, both path words, and the exact named deletion. Temporary names are
distinct operation names, even when their adjacency-matrix index is reused.

The d-position selector therefore uses d temporary births, 2d native moves,
d deletions, and peak order n+1. Its permanent vertex count is unchanged.
Every macro endpoint has the same H and a valid marked source program.
Intermediate temporary graphs need not have that H or belong to the adaptive
language. For newest-to-oldest selection, d=k−1, improving the previous
whole-reversal clock from k(k−1) to **2(k−1) native moves**.

The concurrent [AGE A3](../source_reset_branches_20261103/AGE.md)
independently derives this same stable rotation and clock. Its persistent-
name degree argument also proves a linear lower bound allowing temporary
vertices. These implementation results are credited as parallel work; our
twin increment, existing-pair divisor and guarded balanced renewal have
their separate stated scopes.

The source-to-source map is exact: fixed named graph and ports → explicit
adjacent arc changes → same ordinary class, altered port response → named
marked reentry. What is lost is the original chronological order of the
module's births. The replacement chronological order and slot program are
the required sidecars. Reinterpreting the physical port without changing
the graph is not this operation.

## 5. Arbitrarily large factors can be introduced and erased

Take an actual balanced rooted profile (b,t,t,b). A newborn-root birth N
gives (2b−t,b−t,b,2b). Add k≥2 children at that new root r, retaining r and
the newest physical child z. The resulting source has

    h_k=b(2^k+1),       D(newest pair)=b.             (5.1)

Move z to the oldest position by the selector of §4. Formula (3.1) gives

    D(oldest pair)=gcd(b,t).                         (5.2)

The source is primitive exactly when the original balanced pair was
coprime. This erases arbitrarily large newly introduced factors; the effect
is not limited to the original factor-three example.

For a general selected position A+B+1=k, the actual rooted leaf profile is

    h=b(2^k+1),
    s=b·2^B,
    q=b(2^(k−1)+1),
    F=b(2^(B+1)3^A+2^(A+1))−t·1_(B=0).             (5.3)

The next-birth response is h+F−s. It strictly increases with A by §2. At
the oldest leaf it is

    p_k=2b(2^k+3^(k−1))−t.                          (5.4)

The original source z=(5,2,3,6) corresponds to b=3,t=1 and recovers every
previous formula. The −t boundary term is the original s−q correction from
the marked block law; dropping it would destroy the factor-erasure proof.

**A sharp single-birth stopping boundary.** If the balanced source is
strong, then b≥3 and 0<t<b. No block position, and neither of the two
retained-root choices, supplies any of 3h_k,3h_k+2,3h_k+4 after just one birth.
For a nonoldest leaf the response is divisible by b, excluding the latter
two values. Equality with 3h_k would require

    2^B(2·3^A−1−2^(A+2))=2−2^(A+1),   B≥1.

For A=0 the left side is negative while the right is zero. For A=1 it
requires −3·2^B=−2. For A≥2 the left side is positive while the right is
negative. At the oldest leaf the residual p_k−3h_k is −b−t at k=2, and
at least 7b−t>6b for k≥3. The old root r produces
b(2^(k+1)+1)<3h_k. This proves the obstruction for all such balanced inputs,
not just finitely many observed k. Primitivity is not the missing ternary
branch theorem.

## 6. A reset can return to a new balanced actual source

After the oldest-port selection, perform N at z, creating y and marking it,
then R at y. The inherited identity R(N(h,s,q,F))=(h+2F,h,h,h+2F) gives

    b'=[3·2^k+4·3^(k−1)+1]b−2t,
    t'=(2^k+1)b.                                    (6.1)

The output rooted graph has profile (b',t',t',b'). Its retained pair is
**(y,the last newborn)**, in that order. This is a literal root-child
interface and is ready for another cycle. Its seventeen-query content is
gcd(b,t): before that last seal, the newest-root profile has gcd(h,s,F),
and in (5.3) at the oldest port gcd(h,s)=b and F≡−t modulo b.

Here is the complete controller word, with the starting root in slot 0:

    0, then k copies of 1;
    move the physical slot-1 leaf to the oldest position;
    3,1.

The six-command code is 3*root_slot+discarded_slot. The first command
replaces the root by its newborn. Each 1 retains that root and the newest
leaf. The selector fixes those two names. Command 3 replaces the old root
by the new child of the selected leaf; the final 1 retains that new root
and its last child. The endpoint's chronological source program is separately
reconstructed and replayed, including the changed leaf order.

A k-cycle has k+3 permanent adaptive births, k−1 temporary births and
deletions, and 2(k−1) native moves. During the reset its peak order is the
old order plus k+2; the final order is old order plus k+3. Thus the whole
cycle never exceeds its final permanent order. H-neutrality is asserted at
each selector macro endpoint, not during the temporary native stages.

For k=1 the same profile formula holds with no port move: (b',t')=
(11b−2t,3b), the inherited four-birth N,N,N,R boundary. The factor trap
in (5.1) requires k≥2; there is no distinct oldest/newest choice at k=1.

For coprime odd b,t, exact arithmetic gives

    gcd(b',t')=gcd(2^k+1, b(2·3^(k−1)−1)−t).        (6.2)

Indeed gcd(b',b)=1, and modulo 2^k+1 the first expression in (6.1) reduces
to twice the second argument in (6.2). This is a selection guard on actual
source parameters, not an assertion that arbitrary positive pairs are
realized.

## 7. Closed primitive cycles and an eight-state branching controller

For k=3, equation (6.1) is

    C3(b,t)=(61b−2t,9b).                            (7.1)

If gcd(b,t)=1, 3|t and 3 does not divide b, then the output is coprime and
satisfies the same guard. Any common divisor must divide 9, because it is
coprime to b; but b'≡b modulo three. Starting from the actual rooted C3
(b,t)=(3,1), one cycle gives (181,27), entering the invariant guard.
Thus k=3 can be repeated indefinitely, with actual marked reentry each time.
The balanced counts begin 3,181,10987,666949 and satisfy
b_(j+2)=61b_(j+1)−18b_j. At every cycle the provisional newest pair has
the increasing factor b_j; the selector erases it; and the final pair and
balanced rooted profile are again primitive.

This one-ray family has an exact hostile. Modulo seven its balanced pairs
cycle through (3,1),(6,6),(4,5). All six declared adaptive-stage counts

    2b−t, 3b, 5b, 9b, 34b−t, 61b−2t

are nonzero modulo seven on that cycle. Hence none of its declared birth or
return outputs lies in the 7p lane. Temporary native intermediate graphs
are not part of that declared output language.

Allow k=5 as a second branch:

    C5(b,t)=(421b−2t,33b).                          (7.2)

Under the same primitive three-divisibility guard, this remains primitive
exactly when 3b−2t is nonzero modulo eleven. Let r=t/b in P1(F_11), with
infinity allowed. The branch is forbidden exactly at r=7 and, when allowed,
sends r to zero. C3 always remains allowed and sends

    r -> 9/(6−2r).

The actual initial pair has r=4. Its C3 orbit is (4,1,5,6); after a C5
branch the C3 orbit is (0,7,3,infinity). Thus these **eight projective control
states** suffice. They are not the integer source state: b,t, the entire
literal graph, and the marked program are retained. Both first branches
from the seed are checked separately: C3 gives (181,27), while C5 gives
(1261,99); both enter the primitive three-divisibility guard.

| r | C3 successor | C5 successor |
|---|---|---|
| 4 | 1 | 0 |
| 1 | 5 | 0 |
| 5 | 6 | 0 |
| 6 | 4 | 0 |
| 0 | 7 | 0 |
| 7 | 3 | forbidden |
| 3 | infinity | 0 |
| infinity | 0 | 0 |

This proves a closed **branching actual-source family** with repeated
factor creation, erasure, and primitive reentry. It does not prove that its
integer image is dense or that an arbitrary target has an inverse branch.
The parent-derived [BALANCED_BOUNDARY](BALANCED_BOUNDARY.md) proves a global
sparsity obstruction even after allowing all block sizes.

## 8. Actual seven-prime outputs and verification limits

The short branch C3 then C5 reaches the actual balanced pair

    (b,t)=(76147,5973).

Its next initial N birth has

    H=2b−t=146321=7·20903.

The probe constructs its 18-vertex named graph and independently counts this
H by subset Hamilton DP. This demonstrates that the branching controller
escapes the fixed-C3 modulo-seven obstruction. It does not establish a new
scalar-spectrum record or prime infinitude.

Continuing with a C5 cycle gives the oldest-exit count

    H=226b−t=17203249=7·2457607.

That graph has 24 vertices. Its complete rows, native paths and temporary
cuts are recorded, and its source program is replayed. Its count here is
certified by the proved actual-source transfer (5.4), not by a 24-vertex
subset-DP run. Trial division separately proves both quotients prime.

`native_port_selector.py` imports no repository implementation. Its finite
controls are:

* 285 labelled rooted quotients for the two twin bijections and birth gap;
* 196 actual block positions for the full rooted response, strict order,
  and exact seventeen-query divisor;
* literal C3 and C5 cycles, and a second C3 cycle, with every stage count
  and incoming/outgoing rooted profile checked by subset DP (orders 9,11,15
  at the respective ends);
* each adjacent selector's temporary cut, both native path guards, unique
  changed old arc, and a complete endpoint marked source program;
* the independent 18-vertex 7p count and the exact graph/transfer distinction
  for the 24-vertex witness;
* the projective controller and the fixed-C3 modulo-seven hostile.

Reproduce from the repository root:

    python3 -I -B 05-knowledge/results/source_selection_closure_20261103/native_port_selector.py
    python3 -I -B -O 05-knowledge/results/source_selection_closure_20261103/native_port_selector.py

Both outputs must equal `native_port_selector.json`. The two lane audits
independently checked the twin injections; spectrum supplied the all-A
dyadic inverse in §3. The balanced/controller arithmetic was checked
independently by the coordinator. None of these audits is claimed to replay
all producer controls independently.

Normal and optimized receipts are byte-identical, SHA256
`414bac09a4a1ccc521b26425854ccd100a7960a3be7e5e1a84df65608e912e0e`.
Spectrum independently reviewed §§3,5–7 and the coordinator reviewed the
full note; both reported PASS without corrections. Their review preserves
the distinction between the final pair's content and the new balanced
rooted profile's gcd.

The successful connection is actual rooted source → homogeneous-block
position → exact continuation content → guarded stabilized native transport
→ actual balanced source again. It preserves physical ports and source
incidence at reset endpoints; it changes old arcs and birth chronology.
The literal path/cut ledger and reentry program restore those coordinates.
Coverage requires an operation outside the balanced renewal architecture;
the source-selection mechanism itself is now explicit and repeatable.
