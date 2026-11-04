# Conditional descent needs variable source interfaces

**PROVED logical reductions and scoped obstructions; CONDITIONAL source-bank
closure.** Full positive odd Hamiltonian spectrum coverage remains OPEN.
Assuming exact coverage is useful for separating witness selection from
transport, but it does not establish a prescribed local induction.

The inherited finite decision algorithm and conditional native compiler are
in [AXIOM](../../coverage_machine_induction_20261103/axiom/AXIOM.md).
The two-lane reduction is [THM-4094 plus its full slug](../../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
The closest obstruction is the all-order empty ready fibre at H35 in
[OBSTRUCTION O3](../../ready_return_pentagon_20261103/OBSTRUCTION.md#o3-h35-is-outside-both-ready-images-at-every-size).
The least-used sidecar here is the predecessor's *actual interface fibre*,
which can be empty even though its Hamilton count is attained.

## C1 Restricting outputs to prime and 7p lanes does not remove ready gaps

A root r of T is j-ready when F_r(T)-H(T)=j. The inherited NR reset sends
such a source of count h to count 3h+2j and balanced profile
(3h+2j,h,h,3h+2j). This describes the reset's exact input requirement;
it is not a characterization of all ways to reach that count.

**Inherited exclusion with a lane consequence.** No tournament of H35 has
a one-ready or two-ready root, at any order, by OBSTRUCTION O1--O3.
Nevertheless 35=7*5 is an allowed seven-prime input, and both

    107=3*35+2,       109=3*35+4

are prime. Hence even an unrestricted H-preserving change of source
cannot make the standard NR mechanism supply either branch from H35.
The obstruction is not cured by retaining more ports, taking more samples,
or restricting both endpoints to the designated lanes. It says nothing
against another operation realizing 107 or 109; both counts were already
in the finite attained prefix.

A second-scale control uses the incoming all-order one-ready exclusion at
12833 in [SF4](../../second_feedback_18979_20261103/PROOFS.md#sf4-positive-exterior-families-and-a-new-ready-value-gap):

    12833=41*313,       38501=3*12833+2 is prime.

Although the predecessor is outside the two lanes, the requested successor
is in the ordinary prime lane. The first case is the smaller inherited
hostile; the second demonstrates why a lane-only target restriction cannot
silently be imposed on all recursive predecessors. Primality of 107, 109
and 38501 is checked by complete trial division in `check.py`.

## C1a The same H35 source supports both branches by richer cuts

The ready obstruction is specific to the standard reset. The inherited
[complete small induced bank](../../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md)
contains the following *same named* seven-vertex source S and two induced
nine-vertex extensions. Outgoing rows use labels0,...,8, and old labels0,...,6
are unchanged in both children:

    S, H35:   [102,124,120,113,33,64,16]
    U, H107:  [486,508,504,497,417,448,16,320,64]
    V, H109:  [486,508,504,497,289,64,272,368,32].

Thus this actual source has both 3H+2 and 3H+4 extensions by two births.
These are extracted from an inherited complete bank, not newly attained
scalar values. The independent checker recounts all three graphs, all
root profiles, every attachment cut, and the old-name map. It also checks
the40 saved target representatives: two of the14 H107 classes and four
of the26 H109 classes contain an H35 subset.

A useful joint interface survives at the *same old root0* and the *same
ordered in-ports3->4*:

| Source | Root0 profile (H,s,q,F) | F-H |
|---|---|---:|
|S| (35,3,3,39) |4|
|U| (107,9,9,117) |10|
|V| (109,5,5,121) |12|

The equality q=s survives; readiness does not. None of S,U,V has a
balanced root at any of its vertices. This is an actual graded interface,
with the complete labelled relation available for the next operation.

The107 branch inserts two new vertices inside the transitive residual
chain relative to in-ports3,4 and root0. In READY's chain convention
([ui->3],[ui->4]), its word changes

    11,11,00,01  ->  11,11,00,00,00,01.

Both new vertices beat only old vertex6 among their predecessors, and
new7 beats new8. The first birth has H59. The109 branch first has H57;
its last six residual vertices include the directed triangle
old5 -> old6 -> new8 -> old5, so it leaves this particular transitive
residual chart. Its full cuts are retained in the receipt.

The all-size interpretation of the107 insertion, and its failure to
supply a closed two-branch grammar, is developed separately in the
[graded interior lane](../interior/GRADED.md). The positive finite control
here establishes why an empty ready fibre does not settle the general
source-bank premise of C5.

## C2 Exact spectrum coverage does not supply a balanced representative

A strong rooted tournament is balanced when its profile is (h,t,t,h).
The inherited defect theorem gives delta=0, hence indegree one. If a is
its unique in-neighbor, the equality q=s says that no Hamilton path of
T-r starts at a. By [ADAPTIVE 3.4](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md#3-the-defect-detects-indegree-and-the-legal-inverse),
T-r is weak, and a lies in its terminal SCC. Conversely these structural
conditions imply q=s and F=h, so they characterize strong balanced roots.

**Small all-order boundary.** Neither H11 nor H13 has a strong balanced
root. The cited sharp strong minimum is 15 at order six, so only orders
at most five need examination. For a balanced root, its weak deletion
has at most four vertices. It is either transitive, singleton=>C3, or
C3=>singleton. Adjoining the unique-in-neighbor root to these cases gives
only H in {3,5,9,15} through order five. The independent literal census
of every labelled strong tournament at these orders verifies this set,
with a second subset-DP check of all h and F values.

Thus a hypothesis demanding balanced representatives at *every* attained
value is already false on prime counts. A reformulation must allow finite
base exceptions, variable interfaces, or a return operation that changes
the scalar count. No assertion about infinitely many balanced gaps follows.

## C3 The arithmetic descent and its exact exceptions

Let L be the union of odd primes other than 7 and values 7p for odd primes
p other than 3. Apart from 3, every N in L is nonzero modulo 3. Define

    c(N)=2 if N=2 mod3, and c(N)=4 if N=1 mod3,
    h(N)=(N-c(N))/3.

For N>=5, h(N) is a positive odd integer strictly below N. The only lane
values for which h(N) is forbidden are

    N=23 with h=7,       N=67 with h=21.

Indeed the four possible outputs from forbidden predecessors are
23,25,65,67; only 23 and 67 belong to L. These need explicit base witnesses.
When h is composite, THM-4094's exact factorization using prime atoms,
7p atoms, and carries 49,63,343 constructs h entirely from smaller lane
values. Thus a well-founded induction is available at the scalar level.
It does not select a marked realization of h that admits the requested c.

More broadly, the positive odds outside {7,21} are the scalar grammar
generated by bases {1,3,23,63,67}, multiplication, and the guarded maps
h->3h+2 and h->3h+4. For an elementary proof, descend by N/3 when
3 divides N; the only forbidden quotient produces N63, a base. Otherwise
descend by (N-c)/3. Forbidden quotients yield 23,25,65,67; the first and
last are bases, while 25=5*5 and 65=5*13 use smaller generated values.
Every generator stays inside the allowed set when the target7 is guarded
out. This proves an arithmetic grammar, not a tournament construction.

## C4 A variable readiness parameter preserves a useful arithmetic guard

Rather than demanding readiness j=1 or2 at a forced predecessor, retain
an integer j>=1 and a source of actual count h with F=h+j. NR then gives

    N=3h+2j,       (N,h,h,N),       h=(N-2j)/3.       (C4)

Because F>=h by cutting Hamilton paths immediately after the root, the
parameter j is nonnegative on every actual source. For N not divisible
by3 it is positive, and the allowed residue is

    j=1 mod3 if N=2 mod3;   j=2 mod3 if N=1 mod3.

For 0<h<N, one also has

    gcd(N,h)=gcd(N,j),       provided 3 does not divide N.   (C5)

Proof: 3h=N-2j, and both2 and3 are units modulo the odd N. Thus a prime
output is automatically primitive when 1<=j<N; for a 7p output the exact
extra guard is gcd(j,7p)=1. The first alternative to the impossible
j=1 predecessor for N38501 is j=4,h=12831. This is only a necessary
arithmetic candidate: no 4-ready source at H12831 is claimed.

This changes a rigid empty-fibre demand into a finite selection problem
in j. It does not prove that any candidate is realized or reusable.

## C5 What may be assumed and then bridged back

Here is a precise **CONDITIONAL source-bank reduction**, separating the
hard local existence statement from already available transport.
For each N in L beyond the finitely supplied bases, use the h(N),c(N)
from C3. Assume the following open property:

> There is a finite actual source S_N with H(S_N)=h(N), and an actual
> induced extension U_N with H(U_N)=N. A certificate supplies all names,
> old arcs, fresh cuts, and the resulting joint attachment data.

Because N is in L, U_N has a strong component of count N. That component
contains every nontrivial strong component of S_N whose count contributes
to h(N), except that unit components may require separate name accounting.
For a clean bounded local predicate one may strengthen the displayed
property by demanding U_N strong and explicitly bounding its order by
the inherited strong-order bound. This stronger property is OPEN too;
it is not inferred from bare exact spectrum coverage.

Under the displayed property, THM-4094 and ordinary strong induction give
all lane values: construct h from smaller lane values, select the supplied
S_N in that H fibre, and follow its actual extension. A supplied S_N is
already a witness of h; the point of the induction is a reusable symbolic
rule producing these source certificates from lower-target data, rather
than an unproved universal search-acceptance assertion. This identifies
the missing constructive ingredient without claiming it has been found.

Alignment with a previously constructed named source is a separate task.
If intermediate H values may change, the inherited stabilized single-arc
compiler supplies alignment between any two specified finite tournament
sources while preserving two chosen physical port names. Map them to an
edge of the target with the same orientation, add fresh names when needed,
correct the other arcs by the two-native-move compiler, and delete excess
nonport names. Both physical port names survive throughout. Their oriented edge agrees
at the two endpoints, but a native intermediate may reverse it before the
second move restores it. Each old name that is deleted must appear in the
program's discard ledger.
If M is the larger endpoint order, at most M(M-1) native moves suffice,
with at most one temporary vertex per corrected arc. The two permanent
ports survive; the endpoint graphs and full observer data are recomputed.

This alignment is endpoint-H-neutral when both supplied sources have
count h. It is not a path of H-neutral intermediate graphs, and it does
not preserve every old arc. Demanding either stronger condition creates
a new open problem or an actual obstruction; a generic compiler cannot
be cited as its proof. In particular an induced-source obstruction does
not survive arbitrary edge-changing alignment, whereas an empty ready
fibre such as H35 survives every endpoint-H-neutral alignment.

Finally, assuming exact coverage E itself gives only the inherited finite
witness section and this generic transport. It need not give the deletion
ratio in C3, the readiness equation in C4, a balanced output, or a bounded
neutral path. To remove E, prove acceptance of a local source predicate
such as the displayed bank property, or find a different decreasing
certificate. The value-to-source section, source transport, and return
interface are three distinct obligations.

## Checks and scope

Run `python3 -B check.py` and `python3 -B -O check.py`; both reproduce the
saved receipt. The complete small graph universe is all labelled strong
tournaments at orders3,4,5, not just adaptive or two-apex objects. The
additional H35 control replays40 inherited H107/H109 representatives and
serializes two actual source-preserving branch witnesses. The
ready exclusions are inherited all-order theorems with pinned sources;
their large censuses are not rerun. Arithmetic through10000 is a control
for C3, not evidence of universal spectrum coverage. No new theorem ID,
external priority, or Lean formalization is claimed.
