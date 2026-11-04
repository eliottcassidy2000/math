# The fixed bank fails; a local seed change repairs its first selected failure

**PROVED exact inverse and scoped algorithmic bound; FINITE-EXACT bank
counterexample, local repairs, and selected continuation.** Unlimited seed
switching, prime/7p coverage, and the complete spectrum remain OPEN.

The fixed219-key bank is not closed under both ternary digits. The fourth
all-+2 descendant of170807 is

    4611815 -> 13835447 = 3*4611815+2.

The parent is in the bank's reader image, but13835447 is absent for every
bank seed, every reader length, and every outer prefix clock. A local
change to that parent's seed repairs the miss at the same23-vertex order
and the same clock4, with all old physical names retained. The next two
selected targets are accepted by the old bank again. Thus this is a
finite example of a switching procedure with an enlargement step; it is
not a uniform termination theorem for that procedure.

## Inheritance and the concept board

The closest proved mechanism is the previous
[exact reader inverse](../../zero_ready_component_returns_20261103/inverse/NOTE.md).
Its [three-round tree](../../zero_ready_component_returns_20261103/continuation/NOTE.md)
is the positive control and corrected near miss: fourteen matched branches
do not prove a fourth round. The canonical hostile is the exact219-bank
miss above. The least-used relevant sidecar is Session181's
[named defect transport](../../session181_multiport_defect_20261104/transport/PROOF.md):
full named endpoint algebra transports arbitrary operations, while target
acceptance remains a separate obligation.

The live board has five entries: fixed-bank acceptance (anchor); an exact
binary-section decoder (niche); enlarging and then shrinking the selected
seed (wildcard); physical mark/clock return; and joint attachment data.
The counterexample changes the first into a guarded enlargement problem.
The repair changes the seed and selected physical mark, but retains the
clock and old names. The later accepted pair shows that seed size need not
increase on every accepted step. The decoder decides each declared
arithmetic image without paying for a linear scan; it does not make that
image closed under ternary operations.

The finite-template obstruction in
[Session181 branches](../../session181_multiport_defect_20261104/branches/PROOF.md)
concerns fixed strong quotient templates with transitive module lengths.
Arbitrary binary reader words change their run structure and are outside
that stated template family. Conversely, neither their arbitrary length
nor the finite dimension of their arithmetic carrier proves acceptance.

## S1. An exact inverse by binary sections and local monotonicity

Write A(r) for the one-apex Hamilton count minus1, with A(0)=0. The inherited
reader arithmetic, also derived in [ARITH-1](../arithmetic/NOTE.md), gives

    A(2u+1)=A(u),        A(2u)=A(u)+u+(u mod2),
    0<=A(r)<=r.

For one actual seed with fields(h,e,t,F,p,q,P2,S3) and clock M=2^s, s>=1,
set K=t-q+Mp, d=q-e and D=N-h-MF. Here K,q>0 and0<=d<=q. Inversion is exactly

    f(r)=Kr+q A(r)+d(r mod2)=D,             r>=0,          (S1)

with r=0 representing the empty word. Positive r decodes to its unique
low-to-high complemented binary word ending0. The parity sections are

    r=2u:   (K,q,d,D) -> (2K+q,q,q,D),
    r=2u+1: (K,q,d,D) -> (2K,q,0,D-K-d).                 (S2)

These are equalities, so taking both branches neither misses nor invents
a solution. Negative residuals are impossible; residual0 has only u=0.
A common divisor of K,q,d must divide D and supplies an additional safe
pruning test.

There is no global monotonicity assumption. In a section let U=floor(D/K).
For0<=u<U,

    f(u+1)-f(u) >= K-q u-d.

Thus when K>q(U-1)+d, that section is strictly increasing on its entire
possible range[0,U], and exact integer bisection is complete. Otherwise
apply(S2) again. This is the algorithm in [inverse.py](inverse.py).

Termination follows already from doubling K in every section. A sharper
bound follows because residuals never increase and d remains in[0,q]. At
depth j, K_j>=2^j K. Once K_j^2>qD the monotonicity test must hold. Hence
there are O(1+sqrt(qD)/K) section nodes. Each leaf needs at most O(log D)
bisection evaluations, and A is evaluated in O(log D) elementary integer
steps. This is a square-root-scale algorithm up to logarithmic factors,
not a polynomial-in-log-N or finite-state acceptance theorem. Summing the
clock costs for one seed uses the geometric growth K>=2^s p, while the
number of possible clocks is O(log N).

The local test is essential. For the actual parent seed at clock4,
(K,q,d)=(1620,11,2), but

    A(510)=256, A(511)=0,
    f(510)=829016 > f(511)=827822.

More generally A(2^m-2)=2^(m-1) and A(2^m-1)=0, so the adjacent difference
is K+d-q2^(m-1), negative eventually for every fixed q>0. An ordinary
single binary search of the whole reader image would therefore be invalid.

## S2. A complete fixed-bank miss

The inherited bank consists of all labelled strong seeds through order5,
reduced by their eight fields, plus35 inherited strong order-six
representatives with every mark. It has219 keys; each retains an actual
rooted graph in [bank.json](bank.json). The counterexample uses exactly
this bank, not all tournaments or all possible seed enlargements.

For13835447 the section decoder visits56,499 states and evaluates305,060
scalar values. It finds no witness. The independent inherited interval
inverse checks1,473,275 candidate integers and also finds none. Both include
all possible clocks using h+2^s F<=N; every possible word is represented by
one of the finitely bounded r values. There is no sampled word-length cap.
The neighboring required +4 value13835449 has five bank witnesses.

No modular obstruction is asserted here: the exact integral image is
empty at the target. Nor is this a new forbidden Hamilton count:
13835447=37*179*2089 is realized by the repair below. The failed implication
is from a finite bank's successful three rounds to its closure. The strongest
survivor is the exact per-target decision procedure plus literal positive
witnesses whenever it accepts.

## S3. A two-birth seed recipe repairs the miss and keeps the actual clock

The actual H4611815 parent uses the six-seed

    rows=(38,60,24,49,33,4), root index5,
    names=(6,4,7,8,14,31),
    fields=(29,9,31,83,100,11,320,177),
    reader=001001110010, prefix=(50,51,1,2).

Its completed component fields satisfy(h,F)=(61687,284383).
Birth a seed vertex with outgoing mask37, then another with outgoing mask55,
and select old seed index1 as the new mark. The resulting strong seed is

    rows=(38,124,24,241,97,4,165,55),
    fields=(269,19,679,338,940,59,4520,5246).

Use reader1111000100 and retain clock4. This yields

    (h,F)=(594839,827538),
    h+16F=13835447.

The physical implementation promotes existing tail names21,22 into the two
extra seed slots. All six old seed names and every old induced seed arc
remain literal, and the entire prefix50,51,1,2 and its incidences remain
unchanged. Physical4 becomes the selected in-port; physical31 remains a
retained named port. All23 old names survive, with no permanent births or
old-name deletions. Exactly18 old arcs change in the full relation.
The complete target, source, physical identity map, changed arcs and seed
recipe are in [named_repair.json](named_repair.json).

The two seed births are a generative description relative to the old
six-vertex seed. They are not claimed to be permanent births in the full
23-vertex graph:21 and22 already exist there. The actual full-source change
is the recorded relation rewrite. Native legality and full joint phase
checks are recorded separately in the sibling transport lane.

### Cheap local controls identify the missing ingredient

The menus in [growth.json](growth.json) have exact word/clock quantifiers:

- One old-seed arc reversal or one arbitrary birth, allowing every mark:
  538 rooted candidates and two hits. Both require clock1 and a changed
  selected mark. Both actual strong boundary graphs are saved in the checker.
- C3 inflation of one seed vertex, every mark, clock>=2:48 rooted seeds and
  no hit.
- Two births whose first outgoing mask is11 or37, whose second is arbitrary,
  retaining strong seeds, all marks and clock>=2:2,016 rooted candidates and
  three hits. None uses old mark index5 or either newborn; those restricted
  marks account for756 failed candidates. The three surviving selected
  indices are3,4,1, with clocks3,6,4 respectively.

Thus arbitrary extra vertices alone do not repair the declared restricted
menu: switching to another existing old seed mark is essential for its
positive clock>=2 solutions. This is scoped necessity within that menu,
not a global necessity theorem for every realization of13835447.

## S4. The next selected pair switches back into the old bank

From the repaired value, the next required pair is

    41506343=3*13835447+2,
    41506345=3*13835447+4.

Each has nine witnesses in the original219-key bank. The new eight-seed
alone accepts neither target at any word or clock. This finite policy trace
therefore enlarges the bank on a miss and then switches back to older seeds;
it does not merely repeat one successful enlarged source family.

Eight +2 witnesses and one +4 witness meet the additional guards of clock>=2
and enough vertices to retain all23 old names. No witness in this bank keeps
clock4. [followup_named.json](followup_named.json) chooses actual26- and
27-vertex templates with prefix50,51, outer root39, active mark4, spare6,
and retained old mark31. They need three and four new permanent names,
and change20 and23 old arcs respectively. Full graphs and required arc
alignments are retained. This file claims actual reader templates; native
programs and independent full-interface verification are separate receipts.
No further branch layer was generated.

## Reproduction and the remaining induction obligation

Run `python3 -B check.py` and `python3 -B -O check.py`; the receipts agree
byte-for-byte. The checks compare14,448 section problems with direct integer
enumeration and2,048 address values with the independent inherited streaming
evaluator, replay all14 inherited tree targets, certify the fixed-bank miss
by the old interval inverse, and reproduce the stated local menus and
selected next-target acceptance. All validity gates remain active under -O.
Use named.py and followup_named.py to regenerate the physical templates.

A prospective recursive policy is now concrete: try the bank; on a certified
miss, generate small local descendants of an actual retained seed, allow
selected-mark switching, and require a count witness together with name,
clock and joint-interface guards. **OPEN:** prove that such enlargement
always terminates with both requested digits and that accepted outputs can
continue indefinitely. The exact integer inverse and the generic named
transport solve decision and realization of a found witness; neither proves
that the required witness exists. Bare3H+2/+4 recursion also does not by itself
cover the residue-zero lane or prove the prime/7p classification.
