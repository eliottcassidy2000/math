# Frames, actual attachments, and terminal contraction

2026-10-28. Local identifiers below are scoped to this file. PROVED means the
elementary argument given here; computations are FINITE-EXACT. These are not
reserved global theorem IDs or claims of external priority.

## N1. Free frames give an exact LCA carrier inside a fixed atom — PROVED

Let G=Aut(T) act faithfully on the vertices of a fixed finite tournament T.
A base b is an ordered tuple with trivial pointwise stabilizer. Its orbit
F=G.b is a free transitive G-set. Choose b as a temporary reference, and write
the frames at s attachment sites as g_1.b,...,g_s.b. Simultaneous change of
reference is left multiplication of all g_i. For each oriented edge ij of a
tree on the sites retain h_ij=g_i^{-1}g_j, with h_ji=h_ij^{-1}.

Multiplying h along the unique i-to-j path recovers g_i^{-1}g_j. Equivalently,
transport each branch to its LCA and compose there, preserving order in a
noncommutative group. Thus tree transport determines every joint frame orbit.
There are exactly |G|^(s-1) possibilities after fixing g_1=1. If the future
queries include all relative frame positions, all these possibilities are
distinguishable, proving minimality for this query language. A common frame
is gauge; the other alignments are data. LCA creates no information: it saves
storage only because path products are constrained by coherent transport.

For nonfree ports, pair orbits are stabilizer/double-coset data. Multiplying
these sets can introduce alternatives with incompatible intermediate
witnesses. A base removes that ambiguity; ranks, stabilizer orders and separate
pair types do not. This theorem concerns fixed copies of T with specified
frame orbits. It does not compress arbitrary paths or atom changes.

FINITE-EXACT: all automorphisms, all smaller candidate bases, and all two- and
three-distinct-port tuples were checked for C3, cyclic T5, Paley P7,
C3[C3,C3,C3], and the hostile T81/T27 pair. Minimum base sizes are respectively
1,1,2,3,3,3. P7 has ten ordered triple orbits but only eight all-pair records.
For T81/T27 all-pair records separate triples in this finite test, while root
stars merge 72 into67 and216 into201 respectively. That is no all-arity claim.

## N2. Missing joint data changes an actual attachment count — FINITE-EXACT

Use P7 on Z/7 with i->j when j-i is1,2,4. The inherited triples
a=(0,1,3), b=(0,1,5) agree on every ordered pair orbit, but have respectively
zero and one common out-neighbour. Their distinct joint orbits and literal
F8/Frobenius model were proved in the inherited native-Paley packet.

Here add three marker vertices m_j. Marker m_j beats exactly port_j among the
old vertices; every other old vertex beats m_j. Set m_0->m_1->m_2 transitively.
This fully specifies a ten-vertex tournament. Subset Hamilton DP gives
H(T_a)=4281 and H(T_b)=4315. The other five marker orders also separate them:
each pair is either4281/4315 or4287/4325. Marker roles specify the construction;
the resulting H is an uncoloured tournament invariant. Thus the witness does
not rely on treating marker colours as part of the output observation.

This is an actual finite attachment counterexample, not a statement that
these marker gadgets are a quotient-walk module substitution. The modular
composition grammar needs its own endpoint/multiplicity sidecar.

## N3. Field phases require an orientation guard — PROVED / FINITE-EXACT

Inherited [F4](../frobenius_lca_metagraph_20261028/THEORY.md): every cycle of a
tournament automorphism is odd. An even cycle has opposite vertices swapped
by a power, which would reverse their arc. Therefore on r distinct normal
basis ports cyclic shift by k is orientable only when r/gcd(r,k) is odd.
Conversely a permutation with only odd cycles has no reversed unordered-pair
orbit, so orientations can be chosen orbit by orbit.

For F64 (r=6) the allowed shifts are0,2,4; for F256 (r=8) only0. The code checks
all shifts and a reversed pair for every rejected one, using actual normal
bases over x^6+x+1 and x^8+x^4+x^3+x+1. It realizes the positive F64 action in
S=C3[C3,C3,C3]: identify the six normal-basis ports with (0,3,1,4,2,5).
The automorphism (0 1 2)(3 4 5), fixing6,7,8, sends port_j to port_(j+2).
S is strong. This intertwines the surviving Frobenius square with an actual
tournament action; it does not identify the field's alternating pairing with
the tournament relation. Frobenius may still describe external transport
between presentations when it is not an automorphism of one tournament.

## N4. Exact minimal finite-horizon unnamed carrier — PROVED recurrence

Fix a maximum order M. X_n consists of tournaments whose identity order is a
Hamilton path; the binomial(n-1,2) nonconsecutive arcs are free bits. Write
q_n(x) for ordinary isomorphism class, tau for flipping all free bits, and
b_n(x) for the actual blue fixed-point indicator under converse plus reversal
of path positions. The observation is (q_n,b_n). Allowed operations are tau,
fresh uniform Hamilton-path re-marking (reset), and ordered adjoining of one
singleton on either side, J_L,J_R, while order remains <=M.

At the top set S_M(x)=(q_M(x),b_M(x)). Downwards set

    S_n(x)=(q_n(x),b_n(x),C_(n+1)(J_L x),C_(n+1)(J_R x)),
    C_n(x)=(S_n(x),S_n(tau x)).

Merge equal tuples. This is the coarsest deterministic carrier preserving
the specified observations and transitions. Proof: complement swaps the
two components; adjoining reads the appropriate next-level component. Any
distinguishing component recursively supplies an operation word ending in
unequal observations, proving necessity by downward induction. Reset adds no
split because all states of an ordinary class have the same reset row.
Automorphisms act freely on Hamilton paths, so that row is uniform on the
H/|Aut| presentations of the class. Counted reset also has class-constant
rows; the multiplicity is |Aut| per target presentation.

At M=7, cell counts for n=1..7 are1,1,2,8,64,1019,24754, compared with
1,1,2,7,50,740,24754 before adjoining singletons. Thus even this small menu
already reconstructs every presentation through n=5. At n=6 exactly five
two-element cells remain; all other cells are singletons. Each residual pair
is separated by a context reaching order8, tested by exact uncoloured
isomorphism. Hence the carrier with horizon8 is discrete on every X_n for
n<=6. The receipt includes the five contexts and their terminal masks.
Ordinary classes are computed by legal path relocations, not supplied by
Burnside. Through n=5 every class is independently checked by all vertex
permutations. Through n=7 their numbers agree with A000568.

This theorem is deliberately typed: ordinary class is among the observations;
no individually named re-marking is available; M is finite; attachments are
singletons. It is not a universal minimum for arbitrary joint-port gadgets.
The inherited all-order sandwich decoder has a larger operation menu.

Starting from P7, exhaustive reset/tau reachability at n=7 covers all456
ordinary classes and all32768 presentations; its nine path orbits have seven
immediate complement partner classes. Those partners have trivial automorphism
groups and H in{151,153,155}. A frame inside P7 cannot describe this whole
closure without recording the changed tournament and path information.

## E1. Contract a terminal cycle stage before multiplying — PROVED

Use the inherited descending odd-cycle DP. At degree d just before processing
1-cycles, all live future gcd coordinates reduce to c, the number of already
chosen cycles. Let O(d)=d!/2^v2(d!) and let W_(d,c) be the inherited integral
weight scaled by O(d). Define P_d(t)=sum_c W_(d,c)t^c. Appending k fixed points
adds binomial(k,2)+kc pair-orbits. Its contribution to the final scaled
coefficient at n=d+k is

    O(n)/(O(d)O(k)) * 2^[binomial(k,2)-v2(k!)] * P_d(2^k).

The prefactor is integral by the inherited odd-factorial transfer proof.
Evaluate P_d at powers of two by sparse Horner shifts, then multiply by the
prefactor once for that (d,k). This is exactly the sum of the old transitions;
no histories or contributions are approximated. At k=0 take P_d(1).

For pair-free cores, the inherited substitution p_l=x^l/(1+x^l) gives excess
e>=0 at this stage. The same value P_d(2^k) is reused at n=d+k+e with prefactor
O(n)/(O(d)O(k)) times (-1)^e binomial(k+e-1,e), for k>=1. Signed weights are
retained. Divide the final total by O(n), checking exact divisibility.

This is an implementation improvement to existing unit-tail/moment formulas,
not a new Burnside formula or a new quotient of the prefix state space. It
also explains the negative control: contracting both lengths1 and3 reduces
prefix transitions but leaves more terminal polynomial groups. Fewer states
do not by themselves imply a faster evaluation.

## E2. What has, and has not, been reconstructed

Inherited reset-trace theorem E1 in the previous packet gives a(n)=tr(R_n).
The odd-cycle action on labelled tournaments evaluates the same count by the
classical Davis/Burnside formula. Our new kernel evaluates that scalar formula;
it does not generate the black/blue metagraph through its largest numeric n.
The earlier full metagraph replay through8 and independently supplied order9
bank retain their previous scopes. This packet rebuilds actual small
metagraph components to7 for the new closure test.

Strong components of a tournament form a unique ordered sequence. Consequently
the nonempty strong series is 1-1/A(x); A051337 uses constant term1 and hence
2-1/A(x). This transfers exact scalar coefficients, not attachment relations.
The equinumerosity of tournaments and automorphism-sign-even graphs similarly transfers counts;
without an operation-intertwining map it does not transfer complements, ports
or transporters.

## E3. One complete root-switch row per class — PROVED implementation rule

Use the inherited [root-switch theorem, section3](../session79_tournament_enumeration_20261027/ROOT_SWITCH.md).
The pointed state set R_N is the disjoint union of parent attachment-mask
orbits. For an unrooted child U with vertex-orbit sizes m_1,...,m_r, the entire
root-switch block is 1_r(m_1,...,m_r). Its support contains exactly the r
pointed versions of U. Thus one computed row identifies every index in that
block, along with the common row and all diagonal masses.

Algorithm: traverse pointed-state indices in their fixed order. If an index
is already assigned to a block, skip it. Otherwise construct its child,
classify its N deletion parents with the inherited exact transport, compute
the full row, and assign every supported index. Distinct discoveries have
disjoint supports. Every pointed state is eventually assigned. Choosing the
first index of each discovered support gives exactly one child representative
per ordinary class, in the same order as the old all-row traversal.

This is a direct computational use of the existing block theorem, not a new
classification theorem. Counts are used only for post-generation checking.
All attachment-mask orbits are still generated and indexed. Parent transport
is still exact. The number of computed rows is reduced from p_N to a_N, and
deletion/classification calls from N p_N to N a_N. The full matrix can be
recovered from the stored support and one row per block if later required.
The current implementation retains representatives, automorphism counts and
block-size statistics, not a full black/blue adjacency matrix.

FINITE-EXACT through9: the new producer agrees on every ordered child
representative with the inherited all-row generator through8 and with the
saved191536-member bank at9. At9 it avoids1,525,072 duplicate rows and reduces
classification calls from15,449,472 to1,723,824. Its deletion-parent cache has
186142 entries versus922030 in the saved all-row receipt. The source's9-bank
is a post-generation oracle, not an input to discovering children.

This is also the precise information connection to N4 and E1: a future
response that is already determined should be reused. In N4 the unobserved
response may force a split; in E1 a terminal weighted sum can be contracted;
here an exact complete response determines a whole equivalence block. These
are different justified operations, not a licence to merge arbitrary ports.
