# Frobenius, LCA transport and the reset trace of the actual metagraph

2026-10-28. Packet labels F1--F4, C1--C2 and E1--E3 are local, not canon
reservations. **PROVED** means an elementary argument below; **FINITE-EXACT**
means precisely the listed exhaustive universe. Inherited results are linked.
There is no claim to a first enumeration formula for tournaments or to an
efficient complete classification of arbitrary strong atoms.

**Concurrent integration.** While this packet was being checked, commits
dd1e1bbd6e and 0884018279 supplied the [normal-port/root-switch packet](../relation_tower_frobenius_ports_20261028/README.md)
and [Session86 aggregate composition theorem](../session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md).
The normal-frame minimum and all-order enumeration mechanisms are shared
or inherited, not exclusive new discoveries here. Their broader orbit-pair
tests and order-nine bank are incorporated below; the one-line form tests
and joint-parity/adaptive query refinements address additional questions.

## F1. One marked subfield line can determine the entire incidence form

Let r=2m, F=F_(2^r), K=F_(2^m), and phi(x)=x^2. Work with alternating
F2-bilinear forms on F satisfying B(phi x,phi y)=B(x,y). The distinguished
form is

    B_0(x,y)=Tr_(F/F2)(x y^(2^m)).

It is symmetric, alternating and nondegenerate. For alternation, x^(1+2^m)
belongs to K, whose absolute trace from F is zero; nondegeneracy follows
from the nondegenerate trace pairing and the invertibility of y->y^(2^m).
Every subfield line xK is totally isotropic: for a,b in K, x^(1+2^m)ab
belongs to K. Frobenius invariance follows from Tr(z^2)=Tr(z).

The inherited [P3 trace reconstruction theorem](../free_bits_signed_incidence_20261028/THEORY.md)
proves, for every m, that invariance plus isotropy of the ENTIRE subfield
spread forces B=B_0 once nondegeneracy is required. The following stronger
single-line tests are **FINITE-EXACT** at r=4,6,8:

| field | invariant alternating space dimension | nondegenerate forms | spread lines | non-base lines each determining B_0 |
|---|---:|---:|---:|---:|
| F16 | 2 | 2 | 5 | 4 |
| F64 | 3 | 2 | 9 | 8 |
| F256 | 4 | 8 | 17 | 16 |

The base line K imposes no additional equation in these three tests. Every
other single line suffices, but the mechanism matters:

- F64: six lines, each of Frobenius period six, impose two independent
  linear equations. Two period-two lines impose only one: their compatible
  linear form space still has dimension two, but only B_0 is nondegenerate.
- F256: all sixteen non-base lines have period eight and impose three
  independent equations, leaving the line {0,B_0}; nondegeneracy removes 0.
- F16: four period-four lines each impose one equation, again leaving {0,B_0}.

Thus a realization guard can eliminate a remaining linear ambiguity. It
must actually be known; “invariance” alone does not imply nondegeneracy.

**Concrete F64 witness.** In the polynomial basis modulo x^6+x+1, the
period-two line with nonzero points {3,17,18,40,43,57,58} has basis
{3,17,40}. Its isotropy constraints have rank one on the three-dimensional
invariant form space. Of the two globally nondegenerate masks 656 and
27760, only 27760 survives. Integers denote polynomial coefficient bits;
form masks list upper-triangular matrix entries lexicographically.

**Mechanism / independent audit.** In a normal basis, phi is a cyclic shift.
An invariant alternating matrix is symmetric circulant with zero diagonal,
determined by m coefficients a_1,...,a_m. Its circulant polynomial is

    c(t)=sum_(i=1)^(m-1) a_i(t^i+t^(2m-i))+a_m t^m.

It is invertible iff gcd(c(t),t^(2m)-1)=1 over F2. For r=4,8 the second
polynomial is (t+1)^r, so nondegeneracy is exactly a_m=1. For r=6 it is
(t+1)^2(t^2+t+1)^2, giving a_3=1 and a_1=a_2. This independently explains
the counts 2,2,8. `audit_fields.py` enumerates in these coordinates and
changes basis back, reproducing every nondegenerate mask and every line's
isotropy verdict. The producer instead solves the full invariance equations.

**Boundary.** Single-line reconstruction for all even degrees is OPEN here.
The known all-degree theorem requires the full spread. These are field
incidence forms, not orientations: in characteristic two an alternating
form is symmetric and cannot by itself orient tournament edges.

## F2. A determined form does not determine relative Frobenius phase

Suppose two pieces carry the same chosen normal basis e_0,...,e_(r-1),
but the second is transported by phi^k. Their full joint port matrix is

    M_k(i,j)=B_0(e_i,phi^k e_j),       k in Z/r.

**PROVED.** All r matrices are distinct. If M_k=M_l, nondegeneracy implies
phi^k e_j=phi^l e_j for every j, hence phi^k=phi^l on F and k=l mod r.
A common shift of both frames leaves every matrix unchanged. Therefore
relative phase is real observable information; absolute common phase is gauge.

For s full frames with independently arbitrary phases, simultaneous
Frobenius leaves exactly r^(s-1) joint states. A spanning tree of full
mod-r phase differences is a smallest carrier in cardinality for all these
joint matrices. Necessity follows from their separation; sufficiency follows
by summing differences along paths. Its optimally packed fixed binary storage
is ceil((s-1)log2(r)) bits, not necessarily a separate rounded code per edge.

The information needed to STORE a phase differs from the cost to READ it.
We exhaust all normal generators: 8 in F16, 24 in F64, 128 in F256.

| field | phase possibilities | binary storage bits | minimum worst-case individual incidence queries | arbitrary joint parity queries |
|---|---:|---:|---:|---:|
| F16 | 4 | 2 | 3 | 2 |
| F64 | 6 | 3 | 5 | 3 |
| F256 | 8 | 3 | 4 | 3 |

Individual queries ask one M_k(i,j), with adaptive choice allowed. The costs
hold for every normal basis in the listed finite fields when both pieces
use that SAME normal Frobenius orbit. They are not bounds for arbitrary
pairs of different normal orbits. The exhaustive
decision-tree optimization and an elementary lower bound agree:

- F64's two possible trace rows have weights one and five. Each query
  singles out just one phase; five queries are necessary in the worst case.
- F256's four possible trace rows have weights three or five. A three-query
  binary tree distinguishing eight states would need a 4/4 first split,
  which no individual entry supplies. Four fixed queries suffice.
- F16's row has weight three and similarly needs three individual queries.

**PROVED aggregate-query mechanism.** Let the available complete row be
v_k=(M_k(0,j))_j. These r vectors form an invertible circulant matrix. For
each binary digit of k, solve for a linear functional taking that digit's
value on all v_k. Exactly ceil(log2 r) such parity queries identify k, and
the binary decision-tree lower bound makes this optimal. For the selected
normal bases the parity masks are [21,24,6] at r=6 and [170,153,210] at r=8.
The receipt checks every phase explicitly.

This is a change of interface, not free access to an unavailable parity:
if only single incidence entries can be requested, computing one aggregate
may itself require several entries. If the parent already owns the required
joint parity, the corresponding phase bit needs no extra observation.

### Different normal orbits: the integrated broader test

The incoming [T2--T3 normal-port theorem](../relation_tower_frobenius_ports_20261028/THEORY.md)
allows different normal orbits for the two frames. We replayed its entire
ordered orbit-pair query universe, then optimized ADAPTIVE trees as well:

| field | ordered normal-orbit pairs | minimum fixed-entry queries: pair counts | minimum adaptive-entry queries: pair counts |
|---|---:|---|---|
| F64 | 16 | 3:8, 5:8 | 3:8, 5:8 |
| F256 | 256 | 4:192, 6:32, 7:32 | 4:224, 7:32 |

In particular, 32 F256 orbit pairs improve from six fixed probes to four
adaptive probes. Every pair admits THREE joint parity probes: the matrix
of its complete profiles is invertible because both frames are normal and
B_0 is nondegenerate. The aggregate decoder proof above applies unchanged.
The all-degree worst case r-1 for individual probes is proved in incoming
T2 by choosing dual normal bases, whose profiles are unit vectors. This
does not contradict the four-query same-orbit F256 result.

The incoming full-period-but-nonnormal hostiles were also rerun: F64 has
anchor 2 and targets 34,35 with common profile 36; F256 has anchor 2 and
targets 18,19 with common profile 221. The target pair differs by phi^(r/2)
and the anchor's conjugate span has rank r-1. Full degree alone does not
certify the injective decoder. Incoming T3 gives closure under re-rooting,
restriction and gluing over a nonempty shared port set; disjoint gluing
still needs a supplied relative phase.

## F3. LCA recovers joint matrices only at sufficient phase precision

The inherited [CA1--CA2 cyclic alignment theorem](../session72_even_graph_unifying_carriers_20261026/session72-cross-domain/THEOREM-CYCLIC-ALIGNMENT.md)
already classifies mixed exact periods d_i: a joint state modulo simultaneous
shift is determined by the coherent named differences

    delta_ij=k_j-k_i mod gcd(d_i,d_j).

There are product(d_i)/lcm(d_i) joint states. A subset of pair tests is
sufficient iff its induced graph is connected on every nonempty prime-power
layer {i:p^a divides d_i}. Connectivity of the ordinary underlying tree
alone is not sufficient.

**New explicit F64/F256 specialization, PROVED.** For any number of frames:

- If all periods divide eight, a star centred at a maximum-period frame
  is sufficient and uses the minimum number t-1 of nontrivial edges, where
  t counts frames of period greater than one (zero edges if t<=1).
- If all periods divide six and a period-six frame is present, the same
  star is sufficient and minimal. If no period-six frame is present, keep
  one tree on the period-two frames and another on the period-three frames;
  their minimum total is max(t_2-1,0)+max(t_3-1,0).

Proof: the layers for eight are nested; a maximum-period anchor belongs to
all nonempty layers. The layers for six are {2,6} and {3,6}. If a six is
present their connected subgraphs overlap and their union must connect all
nontrivial frames; a star achieves this bound. Without a six they are
disjoint and need separate trees. Apply CA2. Period-one frames carry no
phase information and require no edges.

`transport.py` independently enumerates all phase states for all 64 ordered
three-period tuples in each field and all eight observation graphs: 1,024
graph cases in total. Example joint-state counts are 6 for (2,3,6), 12 for
(2,6,6), 8 for (2,4,8), and 32 for (4,8,8).

**Hostile short anchor.** For periods (2,6,6), or (2,8,8), the states
(0,0,0) and (0,0,2) have the same two phases through the first frame but
different last-to-last phase. An LCA storing only the anchor's modulus two
cannot reconstruct the omitted port matrix. Keep its higher-precision lift.

With equal periods, let g_v be the accumulated transport from the root to v.
For u,v with lowest common ancestor w,

    delta_uv=(g_v-g_w)-(g_u-g_w) mod r.

The common ancestry cancels, so all cross-child matrices can be generated
from the edge transports without storing each pair independently. With
mixed periods this cancellation must be performed on the appropriate
prime-power layers; it is not permission to raise the precision of a
short-period edge. C6 and C8 are cyclic transport groups: binary encodings
of their elements do NOT turn group addition into bitwise XOR.

## F4. The odd-order obstruction joins Frobenius to tournament counting

**PROVED.** A permutation preserving a tournament has only odd cycles.
An even cycle of length 2h has opposite vertices swapped by its h-th
power, which reverses their directed edge. Consequently every tournament
automorphism has odd order.

On the NONZERO field points, full Frobenius has cycle counts

    F64:  1^1 2^1 3^2 6^9;
    F256: 1^1 2^1 4^3 8^30.

So neither full Frobenius action can be a tournament symmetry. On F64,
phi^2 and phi^4 each have type 1^3 3^20 and fix 2^653 tournament orientations.
Every other nonidentity Frobenius power fixes none. On F256 every
nonidentity Frobenius power fixes none. These tests concern 63 and 255
vertices, respectively; adding zero would change the exponent/count.

For odd cycle lengths l_i, the number of unordered-pair orbits is

    e=sum_i (l_i-1)/2 + sum_(i<j) gcd(l_i,l_j).

Each orbit may be oriented independently, giving 2^e invariant tournaments.
The audit follows ordered pairs until return and independently detects
whether their reverses occur, rather than reusing this exponent formula.
This links the field phase obstruction literally to E2 below; it does not
turn every field incidence model into a tournament model.

## C1. What is closed, what is minimal, and what can be inferred

The operation menu must remain explicit:

| observer and operations | established exact carrier | status |
|---|---|---|
| ordinary class; fresh uniform reset and literal complement | ([T_x],[T_tau(x)]) | inherited U2, minimal |
| joint two-port class; fresh reset and complement | pair of separately rooted endpoint classes | inherited U2, minimal |
| strong two-port state; output full class after singleton hinges on both sides | the joint rooted class is recoverable | inherited U1, all sizes |
| named legal re-marking, complement, all ordered quotient-walk attachments | aligned relations, ordered run ports and multiplicities | inherited C1, sufficient; alignment cannot have bounded size |
| full normal-basis field frames with known B_0; simultaneous Frobenius | coherent relative phases on the sufficient layer trees | F2--F3, minimal in this restricted model |
| aggregate ordinary inflation counts for every rigid palette, bounded root arity | aggregate prime vertex-action cycle index | concurrent Session86 section4, coarsest for that observer |

Sources are [unnamed U1--U2](../unnamed_carrier_dissection_20261026/README.md)
and [strong-atom C1--C5](../strong_atom_continuation_vitali_20261026/README.md).
The last row is not a replacement for the tournament row. Complement can
destroy a module tree, and arbitrary prime quotient walks revisit children.
All ordered path components, their joint endpoints and their multiplicities
still have to be assembled. Two outer endpoints alone remain insufficient.

**PROVED affine compiler.** In full arc coordinates, an affine cell C of
fixed-forward-path presentations is closed under the following specified
operations, with empty cells allowed:

1. Literal complement: translate the nonconsecutive arc bits by one.
2. A specified named legal re-marking: intersect with the equations saying
   every consecutive pair of that path is oriented forward, then apply
   the signed coordinate permutation for relabelling.
3. The blue predicate: intersect with the affine equations sigma(x)=x,
   where sigma is converse followed by position reversal.
4. Ordered join of independent cells: concatenate coordinates and fix all
   cross-arcs. Fixed quotient data and a fixed ordered walk give the same
   sort of affine wiring, provided repeated visits use the SAME variables.

Each claim follows from affine maps and intersections. A disjunction over
walks, variable quotient choices, or a reset weighted by 1/H need not remain
one affine cell. Exact execution may retain weighted mixtures and, in the
worst case, refine to singleton presentations. No polynomial-size bound is
asserted for repeated arbitrary operations. Strongness and nondegeneracy
also must not silently be treated as linear equations.

The exhaustive compiler control contains every affine subset at orders 1--4
(1,1,3,51 cells), every named path permutation, complement and blue slicing:
1,357 operation tests. It also checks 229 ordered cell products of total
order at most six. This finite check supports the symbolic argument; it
does not establish a compact all-order grammar.

For a consistent affine family Ax=b, a target parity h x is already known
iff h lies in row(A). At an LCA, store only independent new rows beyond
the child row spaces, keeping their actual coefficients and right sides.
The [inherited linear-gluing theorem](../session85_arithmetic_codes_free_bits_20261028/LINEAR_GLUING_INFORMATION.md)
gives the rank invoice and separator counting algorithm. Merely storing
the rank or separate child projections is insufficient: the even-parity
code on a,b,c projects fully to ab and bc, whose join invents one extra bit.
Keeping the parent equation a+b+c=0 recovers c from a,b.

Named legal-path guards already reconstruct presentations at every size
(inherited reconstruction theorem). Thus d free coordinates describe 2^d
distinguishable states; this is compression of coordinates, not an allowed
merger of all those states. The smallest all-order unnamed attachment
carrier beyond the established restricted menus remains OPEN.

## C2. The 81/27 hostile charges for an unverified symmetry promise

Let A=C3->1 and choose A or A^op independently at the three children of
C3. All eight resulting order-12 tournaments have H=462915. Their
automorphism orders are

    27 * |Stab_C3(the ordered three-bit decoration)|.

Only 000 and 111 have stabilizer three and Aut=81; the other six have
Aut=27. The 000 and 001 examples have exactly the same 67 LABELLED cyclic
triples. `carriers.py` recomputes their Hamilton counts by subset DP and
their automorphism groups by adjacency-checked permutation backtracking.

Knowing cyclic invariance would impose the two independent parent rows
z_0+z_1=0 and z_1+z_2=0, recovering all three decoration bits from one.
But cyclic triples and H do not establish that promise: applying it to
001 discards a legitimate hostile state. This is the tournament instance
of F1's distinction between a known constraint and an assumed constraint.

The presentation fibres have sizes H/Aut=5715 and 17145; their reset
diagonal masses are 1/5715 and 1/17145. Forgetting the colour stabilizer
would give wrong weights in E1, even though H agrees. At order four,
the two presentations with free masks 2 and 5 on a single blue loop give
Hamilton counts 11 and 13 after a singleton join followed by complement.
Thus line endpoints without the attachment's half-edge data also fail.

## E1. A000568 is the trace of the actual fresh re-marking operator

For n>=1, X_n consists of all tournaments containing the directed path
0->1->...->n-1. Its size is 2^((n-1)(n-2)/2). Let H(x) be the number of
ordered directed Hamilton paths in T_x. Define

    R_n(x,y)=#{p in HP(T_x): R_p(x)=y}/H(x).

**PROVED from the inherited reset lemma.** For an ordinary tournament
class T, each presentation y in its fibre has exactly |Aut(T)| path
preimages. Hence that fibre has c_T=H(T)/|Aut(T)| states and R_n restricts
to J_(c_T)/c_T. Every tournament has a Hamilton path: insert the last
vertex into an inductively constructed directed path at its first suitable
cut. Therefore every ordinary class occurs. It follows that

    R_n^2=R_n,     a(n)=rank(R_n)=tr(R_n)
                    =sum_(x in X_n) |Aut(T_x)|/H(T_x).       (E1)

Here a(n)=A000568(n); set a(0)=1 separately. This is an all-order identity
using the actual reset operator, not a count inferred from complement
graph cycles. In particular, tr(R_n tau) is a different observable: at
n=3 it is zero while a(3)=2; at n=4 it is 2/5 while a(4)=4.

`counting.py` independently builds every reset row through n=6 from all
legal path permutations, checks equal fibre multiplicities and equal
supported rows, and evaluates these traces exactly as rational numbers.

The metagraph convention merging converse classes also needs its orientation
sidecar. If m(n) is its vertex count and s(n) the number of self-converse
classes, then a(n)=2m(n)-s(n). At n=8: 2*3528-176=6880. Dropping which
classes are self-converse loses this lift information.

## E2. The same odd-cycle obstruction gives the full scalar sequence

**CITED classical / derivation included.** Apply Burnside to the action of
S_n on all labelled tournaments. F4 kills every permutation with an even
cycle. For an odd partition lambda=(l_1,...,l_t), let m_l be its part
multiplicities and z_lambda=product_l l^(m_l) m_l!. There are n!/z_lambda
permutations of that cycle type, each fixing 2^e orientations by F4. Thus

    a(n)=sum_(lambda partition n, all parts odd) 2^e / z_lambda.    (E2)

This is the classical Davis formula, recorded at [OEIS A000568](https://oeis.org/A000568).
The contribution here is its explicit connection to the metagraph reset
trace and the tested Frobenius transport, not a new counting formula.
`counting.py` evaluates all n=0..50 with exact rationals; the initial values
are 1,1,1,2,4,12,56,456,6880,191536,9733056. The repository already has
[production structural-DAG counting](../lca_ports_reverse_20261025/counting/README.md)
and higher-order arithmetic receipts; this session does not supersede them.

## E3. Composition as size changes: a complete recursion and its remaining cost

There is a literal extension X_(n+1)=X_n times F2^(n-1): add the next
vertex, fix n-1->n, and choose its other incident arcs. At each size,
legal re-markings build R_n; literal complement supplies the black/blue
lines. This gives an exact all-order reconstruction algorithm, exponential
in the presentation space, with E1 counting ordinary classes.

For ordered-join atoms, every tournament has a unique ordered sequence
of strongly connected components. If A(z)=sum a(n)z^n and S(z) counts
nonempty strong tournaments, then A(z)=1/(1-S(z)). Consequently

    s_n=a_n-sum_(k=1)^(n-1) s_k a_(n-k).

The strong head is 1,0,1,1,6,35,353,6008 for n=1..8. These atoms are not
the same as modular primes. Their sequence decomposition counts objects,
but a complement line may cross or destroy its decomposition; attachment
classification still needs the joint run data in C1.

The existing [Session67 enumerator](../session67_metagraph_enumeration_information_20261026/README.md)
was rebuilt and replayed exactly through n=8, matching its original census.
It uses legal right relocations on all presentations, not Burnside input.
The n=8 result includes 2,097,152 presentations, 6,880 ordinary classes,
3,528 converse-merged vertices, 2,048 blue lines and 1,046,528 black lines.
The replay is credited as a replay, not a new enumeration record.

### Integrated constructive continuation through order nine

The concurrent [T5 root-switch construction](../relation_tower_frobenius_ports_20261028/THEORY.md)
has now generated all 191536 order-nine classes from 1716608 pointed
attachment states, preserving each subset mask through the SAME deletion
parent isomorphism. We incorporated its saved bank and checked every row
is a tournament, all encodings are distinct, its class counts agree with
E2, and its entire presentation-fibre histograms agree with the actual
metagraph through eight. The order-nine construction was not rerun here;
bank type/distinct-encoding checks are not independent isomorphism proofs.
The construction's proof and its own finite run remain the source of
order-nine completeness. Its order-nine Hamilton-fibre census is not supplied.

This gives a second actual class-count projector: the incoming rerooting
operator Q_N satisfies Q_N^2=N Q_N and a(N)=tr(Q_N)/N over Q. It counts
vertex-root fibres; E1's R counts Hamilton-path presentation fibres. Both
have one block per ordinary tournament, with different masses. Neither
allows attachment masks or their shared parent alignment to be discarded.

Session86's [universal rigid-palette equivalence, section4](../session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md)
also identifies the aggregate prime cycle index as the coarsest invariant
for all ordinary inflation counts from finite rigid palettes at bounded
root arity. A base-(K+1) size palette separates every colour monomial, so
the count observer recovers the index. This resolves a precise aggregate
minimum, while individual decorated stabilizers and actual black/blue
transport still require additional data. The81/27 pair is compatible
with that theorem: summing every decoration differs from selecting one.

**Remaining research target.** Find an effective attachment grammar whose
residual states are smaller than complete presentations, preserve the
required joint run laws and transport, and remain closed under repeated
complement. F1--F3 identify provably recoverable coordinates in a controlled
subfamily. They do not prove that every prime quotient has that subfamily's
field representation. All-order completeness of unrooted strong-atom response
histograms remains OPEN despite the inherited exhaustive success through n=8.
