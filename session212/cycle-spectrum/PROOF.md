# The cycle-count spectrum is globally complete in feedback vertex number two

**Status: PROVED elementary constructions; FINITE-EXACT through order nine;
INHERITED FINITE-EXACT at order ten. No literature-priority claim.**

## 1. The question and its boundaries

Write H(T) for the number of directed Hamiltonian paths. Write C(T) for the
number of directed Hamiltonian cycles **modulo cyclic rotation**. A cycle and
its reversal are not identified; in a tournament the reversed sequence is
not a directed cycle of the same tournament. Set C(T)=0 below order three.
Let S_C(n) be the set of C(T) over n-vertex tournaments, and S_C^str(n)
the restriction to strongly connected tournaments.

For n>=3, Camion's theorem says S_C^str(n)=S_C(n) minus {0}. Zero occurs at
every positive order, for example in the transitive tournament. Singletons
are strong under the usual reachability convention but have C=0 here.

The distinction from the path spectrum is absolute: path counts are odd,
and 7 and 21 are forbidden at every order (inherited THM-1370,
`01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`). We prove
that **every nonnegative integer is a cycle count**, and every positive
integer is a cycle count of a strong tournament. This does not prove any
new case of completeness of the Hamiltonian-path spectrum.

Inheritance: the endpoint vector below was explicitly recorded in
`05-knowledge/results/chanpak_increments_audit_C_20261104.md`, section (g),
and the forced one-arc return law is the same mechanism as
`session197/lowcycle/ROOTED.md`, Proposition 1. THM-4889
(`the-hamiltonian-path-spectrum-by-order-and-by-value-exact-censuses-to-order-11-resolved-by-cyclic-triangles-least-orders-strong-windows-level-floors-and-the-zone-law.md`,
full source path in `SOURCES.json`) already records complete
(H,C) data through order ten. The contribution here is an independent route to the complete
cycle-value realization, its sharp feedback-class boundary, and the
explicit extraction/recomputation of cycle spectra. The companion
Session212 binary-comb packet independently proves the same support and
size bound, and strengthens it to gate-free realizations. These are
complementary explanations of one global theorem, not separate novelty
claims. Deleting its two apices also gives the sharp feedback-number-two
realization, jointly with its gate-free condition. No claim that the
endpoint vector or the one-arc-return identity is new is made.

## 2. A binary endpoint supply

Let U_m, m>=3, have vertices 0,...,m-1, all arcs i->j for i<j, except that
(m-1)->0 replaces 0->(m-1). Let e_j count paths ending at j. Then

    (e_0,...,e_(m-1))=(1,1,2,4,...,2^(m-3),1).        (1)
    H(U_m)=2^(m-2)+1; C(U_m)=1.                       (2)

Proof. A path that avoids the reversed arc must be the increasing path,
and ends at m-1. A path that uses it has two increasing runs: an initial
run ending at m-1, then 0 and a terminal run ending at j<m-1. For j=0
there is exactly one choice. For 1<=j<=m-2, each vertex 1,...,j-1 can be
placed independently in either run; all vertices above j must be in the
first run. This gives 2^(j-1) choices. The unique spanning cycle is the
increasing order closed by (m-1)->0. QED.

The subset sums of (1) contain every integer from zero to 2^(m-2)+1.
Indeed the weights 1,2,...,2^(m-3) alone represent 0,...,2^(m-2)-1,
and the two additional unit weights extend that interval by two.

## 3. The forced-entry gate

For any tournament U on m>=1 vertices and B subset V(U), add z and v:

    z -> every vertex of U; v -> z;
    u -> v iff u belongs to B.

Call the result G(U,B). Then

    C(G(U,B)) = sum_(u in B) e_u(U).                  (3)

Proof. The only incoming arc of z is v->z, so every Hamiltonian cycle
uses it. Root the cycle at z. The remainder is a Hamiltonian path of U
ending at some u in B, followed by u->v->z. Conversely every such path
gives one cycle. This is a bijection modulo rotation, with no extra
factor of m+2. It is valid also for m=1. QED.

Combining (1) and (3), for every n>=5,

    {0,1,...,2^(n-4)+1} is contained in S_C(n).       (4)

All positive values in (4) are attained strongly. For the constructed
examples this does not need Camion: U_m is strong, B is nonempty, and
z->U_m, U_m reaches some b in B, b->v->z gives directed return paths.
For positive outputs, the arc v->z is a strong bridge: removing it leaves z with no incoming arc.
In particular, positive cycle counts need no vertex or arc connectivity
beyond strong connectivity.

For k>=2 define

    N(k)=max(5,4+ceil(log_2(k-1))).                    (5)

There is a strong tournament with exactly k cycles at **every** order
n>=N(k), not merely at one order. Thus the global unrestricted spectrum
is N_0 and the global positive strong spectrum is N_(>0). For k=1,
the tournament U_n works at every n>=3. These are upper bounds on least
order, not exact least-order formulas and not assertions that fixed-order
spectra are intervals.

A constructive compiler chooses the largest endpoint weights greedily.
It writes an n-vertex tournament using O(n^2) output bits, with n=O(log k).
No search for a tournament is required.

## 4. The feedback-class boundary is sharp

The feedback vertex number tau(T) is the least number of vertices whose
deletion leaves an acyclic (equivalently transitive) tournament.

If tau(T)<=1, then C(T)<=1: for a feedback vertex x, deleting x from any
Hamiltonian cycle gives a Hamiltonian path of the transitive T-x, and
that path is unique. An acyclic tournament itself has no cycles.

For G(U_m,B), delete v and vertex 0 of U_m. What remains is z followed by
the transitive order 1,...,m-1. Hence tau(G)<=2. Therefore every k>=2 is
realized at every order n>=N(k) with **tau=2 exactly**, and a strong bridge.
This is the minimum feedback vertex number allowing more than one cycle.

The map U -> G(U,B) retains selected endpoint mass, not the total H(U)
alone. The sidecar is B together with the endpoint distribution. Passing
only H(U) loses which subset sums are possible. The binary vector is the
reason the selected sums have no holes. Equal cycle counts of the outputs
do not imply equal path counts or isomorphism.

## 5. Complete small-order spectra

Intervals here include every integer, unlike the odd-only path spectrum.

| n | S_C(n) | status |
|---|---|---|
| 1,2 | {0} | immediate |
| 3,4 | {0,1} | fresh complete fixed-cycle enumeration |
| 5 | [0,3] | fresh complete fixed-cycle enumeration |
| 6 | [0,5] | fresh complete fixed-cycle enumeration |
| 7 | [0,15] union {17,24} | fresh complete fixed-cycle enumeration |
| 8 | [0,54] | fresh complete fixed-cycle enumeration |
| 9 | [0,249] minus E_9 | fresh complete extension enumeration; agrees with inherited census |
| 10 | [0,961] union {971} | inherited complete (H,C) census, re-extracted and checked |

    E_9={215,219,223,226,229,232,...,242,244,...,248}.

Thus the first fixed-order gap is C=16 at n=7. All seven order-seven
holes fill at order eight. The 21 order-nine holes fill at order ten.
No all-order next-order-filling or spectrum-nesting law follows from this.
In particular, the strong windows for fixed C are eventually cofinal by
(4), unlike the bounded strong windows for fixed H from Busch's lower
bound. Early holes in those order windows are not ruled out.

The least orders of values <=249, read from these complete spectra, are:
1 at 3; 2,3 at 5; 4,5 at 6; 6..15,17,24 at 7; the remaining 16..54 at 8;
55..249 outside E_9 at 9; E_9 at 10. Values 250..961 and 971 first occur
at 10. Values 962..970 need at least 11 vertices; their exact least orders
are not determined here. The construction gives an upper bound of 14.

## 6. Reproduction and scope of evidence

`verify.py` independently checks endpoint/gate formulas, literal permutation
controls, feedback number, strongness, cycle witnesses and inherited
summaries. Ordinary and optimized Python must produce identical receipts.

`census.cpp --fixedcycle n` fixes arcs 0->1->...->n-1->0 and enumerates all
orientations of every other pair. Every strong n-tournament can be
relabeled to contain that fixed cycle, so this is a complete positive
spectrum enumeration without isomorphism software. Counts of presentations
are 1,4,32,512,16384,1048576 for n=3,...,8. They are NOT isomorphism-class
multiplicities. Add the transitive zero example to get the unrestricted set.

For n=9, `generate_classes.py` regenerates all 6880 order-eight isomorphism
classes from the singleton by every possible one-vertex extension and
exact color-refined permutation canonicalization. `census.cpp --extend8`
then checks all 6880*256=1761280 extensions. Every order-nine tournament
has an order-eight deletion represented in that list. Duplicate classes
are harmless for spectra; extension frequencies are NOT class frequencies.

The n=10 entry is extracted from THM-4889's complete 9733056-class census,
not independently re-enumerated here. Its tracked source is
`05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/data/census_10.txt`
at f9de668e6caaec77b829c13e3d87cf1ad8db61a3, SHA-256
`4d5ac84cfc46ebb4684d7531da35590880ceef98616a5d1f26df785db9aa7fd8`.
The extractor sums C-lines and verifies 9733056 total and 9355949 positive
classes. The same extraction matches our independent n=8,9 spectra.
A literal independent recount verifies the supplied maximum witness
(H,C)=(15745,971). Inherited input correctness beyond these controls is
not relabeled as fresh independent verification.

`inherited-receipt.json` records the inherited source hashes and universe
counts; `RELEASE_FILES.sha256` pins this packet. The literature search
was limited: Camion, Rédei, Busch and established low-cycle classifications
are inherited, and no priority claim is made for global count universality.
