# Fixed-apex tape kernels and polynomial numerical inversion

**Status:** PROVED deductions from the cited Busch bound and the inherited
positive tape transfer; FINITE-EXACT implementation checks. Independent audit
[PASS](../audit/FIXED_CORE_AUDIT.md). **No new scalar arithmetic gate, all-target coverage theorem, or
Chan–Pak proof is claimed.** The useful advance is an explicit fixed-core,
endpoint-aware physical compression and a substantially smaller exact search
horizon. Polynomial means in the numerical target, not its binary length.

## 1. Inheritance, claim boundary, and connection contract

A fixed labeled tournament A on k apexes and a word w over {0,...,2^k-1}
define T_A(w): its other vertices form a transitive chain in word order, and
bit i of a letter is one iff apex i beats that chain vertex. Every matrix
letter therefore costs exactly one physical vertex. The apex orientation and
labels are retained throughout. H(the empty graph)=1 by convention; an empty
word instead gives T_A(empty)=A and H(A), which need not equal one.

The positive permutation/gap transfer is inherited from Session215 and the
Session216 `inverse/PROOF.md` / `two-apex-advances/PROOF.md` family;
it is restated below, not claimed new. That prior inverse caps coordinates
at m+1 and has a state-space bound (m+2)^(k! 2·3^k). It supplies decidability,
not all-target acceptance. Session214 `surgery/PROOF.md` already uses Busch to
bound a special prime-target census. Thus global count-preserving SCC
compression and the use of Busch to make a finite search are inherited ideas.

The specific output here is the same-labeled-apex theorem for **all** targets,
an explicit word-level deletion test, its endpoint-preserving refinement,
and its fixed-k polynomial numerical-time inverse. A bounded search of
current canon and Session214/216/217 results found no prior statement of this
combined specialization; this is not a literature-priority claim.

Closest hostile: equal H is not a continuation equivalence. With apex arc
0->1, the word [3] and the empty word both have H=1, but after appending [0]
the counts are 5 and 1. Hence the new reduction is applied only to a finished
word, NEVER used to identify prefix states in the transfer search.

Concept board: physical tape words; SCC condensation; endpoint sidecars;
Busch resource cost; positive weighted automata; exact occupancy. The map is
an induced deletion of specified singleton chain SCCs. It preserves H, the
apex core, and optionally all apex endpoint labels. It loses chain labels,
chain endpoints, and continuation behaviour. The missing Chan–Pak ingredient
is still an occupancy theorem, not an efficient evaluator or a length bound.

## 2. Exact compact-core theorem

For an integer m>=1 set

    L(m) = max {2a+3b : a,b>=0 integers, 3^a 5^b <= m}.

The maximum is finite and L(1)=0. Integer comparisons compute it without
floating-point logarithms. It satisfies

    L(m) <= floor(3 log_5 m),

because 3^a >= 5^(2a/3), equivalently 27>=25.

**Theorem A.** For every fixed labeled apex tournament A and every word w
with H(T_A(w))=m, there is a subsequence u of w with

    H(T_A(u))=m,     |u|<=L(m).

The apex set, its labels, all its arcs, and all incidences of retained chain
vertices are unchanged. One obtains u by deleting exactly the chain vertices
whose SCC in the original tournament is a singleton.

**Proof.** A tournament's strong components have a transitive total order.
Every Hamiltonian path is the unique concatenation of one Hamiltonian path
from each component. Singleton components contribute a factor one, so their
deletion preserves H, and deleting singleton components cannot alter any
remaining strong component.

Every nontrivial component C meets the apex set: otherwise it would be an
induced subtournament of the transitive chain. Write c(C) for its number of
chain vertices and s(C) for its total size. Then c(C)<=s(C)-1.

The **CITED** Moon–Busch minimum is

    H(C) >= f(s(C)) = min {3^a 5^b : 2a+3b=s(C)-1}.

This is used only for nontrivial SCCs, which have at least three vertices.
Choose minimizing a_C,b_C and sum them to a,b. Multiplicativity gives
3^a5^b<=product_C H(C)=m, while

    |u| = sum_C c(C) <= sum_C(s(C)-1)=2a+3b <= L(m).

Singleton apex components remain present and cost no chain vertices. If
there are no nontrivial components, m=1 and all chain vertices are deleted.
This includes k=0. ∎

**Precise citation.** A. H. Busch, *A note on the number of Hamiltonian paths
in strong tournaments*, Electron. J. Combin. 13 (2006), N3, Corollary 1.
The inherited canonical statement and three-residue formula are in
`01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`, Setting:

    f(3q)=3·5^(q-1), f(3q+1)=5^q, f(3q+2)=9·5^(q-1).

The equivalent minimization is also recorded explicitly in THM-4861(c),
`THM-4861-the-five-fermat-primes-are-hamiltonian-path-counts-of-tournaments-with-least-orders-3-4-6-8-11-a-vertex-transitive-tournament-has-a-prime-count-only-for-c3-and-fermat-primes-are-free-vertices-of-the-prime-square-sum-graph.md`.
The classical bound is not re-proved or inferred from computation here.
All source checks used commit 38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817.

The coarse constant cannot be improved uniformly in k: take an ordinal sum
of r copies of the strong four-vertex one-apex word [1,0,0], with one apex
per block. It has m=5^r and exactly 3r indispensable chain vertices, attaining
L(m)=3r. The apex core is transitive. This is sharpness for count-preserving
induced deletion, not for the shortest arbitrary representation of 5^r. It
does not assert sharpness as m grows with k fixed, nor that every m reaches
its length cap.

## 3. Explicit word reduction and endpoint-aware version

View a mask as a subset of apexes. For position j write S=w_j,
P_j=union_{i<j} w_i, and Q_j=intersection_{i>j} w_i, with empty union empty
and empty intersection the full apex set. Then chain vertex j is a singleton
SCC **iff**

    P_j subset S subset Q_j,   and every apex in S beats every apex outside S.  (1)

To prove this, a directed triangle through j can have:

- an earlier chain vertex i and an apex in w_i minus S;
- a later chain vertex i and an apex in S minus w_i;
- two apexes a->b with a outside S and b in S.

These are exactly the three failures of (1); three chain vertices never
form a cycle. A vertex lies in a nontrivial SCC iff it lies on a directed
triangle. For the converse of this last standard fact, absence of a triangle
through v forces every in-neighbor to beat every out-neighbor, separating
N^-(v),v,N^+(v) in the condensation. This proves (1).

Precompute which apex cuts S=>A\S are valid, then prefix unions and suffix
intersections decide all deletions in O_k(|w|) time. The resulting reduction
is idempotent. This is a final-word retraction, not a semigroup congruence.

**Theorem B (endpoint sidecars).** There is also a subsequence v with
|v|<=L(m)+2 preserving H and all of the following counts, for every labeled
apex i,j:

- Hamiltonian paths starting at i;
- Hamiltonian paths ending at j;
- Hamiltonian paths starting at i and ending at j.

Delete the same singleton chain SCCs, except retain the first chain vertex
if it is a global source, and the last chain vertex if it is a global sink.
These conditions are exactly w_0=0 and w_last=2^k-1, respectively.

**Proof.** Hamilton paths traverse the SCCs in order. The deletion and
reinsertion of singleton components is a bijection. Retaining a boundary
singleton whenever the original source/sink component is a chain singleton
prevents a new apex endpoint from being exposed. When the boundary component
is nontrivial or an apex singleton it already remains. Thus the bijection
preserves whether either endpoint is each named apex, jointly as well as
separately. At most two additional vertices are kept. ∎

This does not preserve endpoint counts at deleted chain vertices, higher
path-cover coordinates, or arbitrary future contexts. Neither theorem is
permission to erase a large hidden automaton coordinate.

## 4. Exact physical transfer, including endpoints

For completeness, the inherited transfer has dimension d=k!·2·3^k. Index
coordinates by an apex permutation p and k+1 gaps before, between and after
its apexes. A gap has state 0 if empty, 1 if nonempty and its latest chain
vertex fails the next-apex boundary, and 2 if nonempty and that boundary
passes. The last gap has only states 0,2. Initially each permutation's
all-empty coordinate has weight one.

For a new mask s, choose a gap g. If it is empty and g>0, require the bit for
its preceding apex p_(g-1) to be one. A nonempty gap needs no entry test.
Set its new state to 2 if g=k or bit p_g is zero, otherwise to 1. Leave other
gaps unchanged. Sum over all choices with multiplicity; two choices reaching
the same state are both counted. These are nonnegative integer matrices M_s.

Accept a final coordinate exactly if each gap is either state 2, or state 0
with its adjacent apex arc valid when it is an interior gap. Boundary empty
gaps have no test. The final sum is H. To count starts at apex i, additionally
require first gap empty and p_0=i; for ends at j require last gap empty and
p_(k-1)=j. Combining the two tests gives joint endpoint counts. For k=0 only
H exists and equals one for every word, including the empty word.

The proof is bijective: a Hamilton path uniquely supplies its apex order and
the gap assigned to each chain vertex; vertices assigned to one gap occur in
increasing chain order. Conversely a passing assignment concatenates into
one legal Hamilton path. Thus for every actual word

    state(ws)=M_s state(w),   H(T_A(w))=final·state(w),

with the stated endpoint functionals. No formal state is asserted reachable.
Total graph order is k+|w|, hence k+L(m) after scalar compression and
k+L(m)+2 with endpoints. If nonempty tournaments are required and k=0,m=1,
use a singleton instead of the empty convention, with one extra vertex.

## 5. Polynomial numerical-time exact inverse

**Theorem C.** For fixed A on k>=1 labeled apexes, membership of a positive
integer m in {H(T_A(w))} can be decided and a witness constructed in

    O_k(m^(3k log_5 2) log m)

bit operations, with O_k((log m)^2) working bits (write log(m+1) at m=1).
The same bound holds when any specified collection of the apex start, end,
or joint endpoint counts is prescribed, using horizon L(m)+2. These counts
are integers in [0,m]; arbitrary encoded input is first validated, with an
additional cost linear in its encoding length, and impossible out-of-range
counts are rejected. The API accepts both tuple and JSON-list constraints.

Enumerate every word to the proven horizon, carrying the exact transfer
state along a depth-first tree. Theorem A or B guarantees completeness.
The alphabet size is alpha=2^k; there are O_k(alpha^L) nodes, and adding two
to L is only a fixed-k factor. Each transition costs O_k(1) integer additions.
A coordinate counts partial assignments and is at most (k+1)^|w|, so all
integers have O_k(log(m+1)) bits. DFS needs O(log(m+1)) state vectors, giving
the claimed time and working-space bounds. No arbitrary subset-DP evaluator
is used in this complexity estimate. The k=0 class is simply {1}.

Optional safe pruning rejects a prefix only when its current H exceeds m.
Hamilton-path count is monotone under adding vertices: insert a new vertex
at a legal position in every old Hamilton path, choosing deterministically;
removal recovers the old path, giving an injection. Hidden coordinates never
justify pruning. Prefix normalization would be UNSOUND by the hostile above.

For k=2 the numerical exponent is 6 log_5 2 = 2.584059..., and the returned
graph has at most 2+L(m) vertices (4+L(m) for endpoint constraints). A full
census up to M can share this same tree with L(M), but storing all witnesses
adds its output storage to the stated single-target working-space bound.

This is not polynomial time in log m, and it does not show that all odd
m except 7,21 are accepted. No enlargement of known finite coverage is
claimed. The endpoint-constrained theorem is constructive and exact even
though endpoint universality is not asserted.

## 6. Reproduction and checks

Run `python3 verify.py > VERIFY.json` and `python3 -O verify.py > VERIFY_OPT.json`
from this directory, and compare outputs. Checks use explicit exceptions,
not assert statements. `python3 automata.py` emits the compact example target
compiler output. All code is standalone; no large inherited data is copied.

The exhaustive literal universe is every labeled apex core for k=0,1,2,3
and all tape words through lengths 4,7,5,3 respectively: 7,670 cases. An
additional 90 seeded random cases use k=2,3,4 and tape lengths 4..7. Independent
vertex-subset DP checks H and apex start/end/joint endpoint counts, Floyd reachability checks
the deleted singleton SCCs, and literal reconstruction checks both normal
forms and the size bounds. Small C3/S4 ordinal products check sharp boundaries.
The complete fixed-two-apex census through m=63 returns every positive odd
value except 7,21; every returned graph receives independent subset DP. The
single-target compiler reruns 1,7,21,23,63, and six boundary/nontrivial endpoint-constrained targets. These are controls, not new coverage.
