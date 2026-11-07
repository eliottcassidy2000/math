# Hamiltonian-cycle count levels and exact inverse size sets

Session212 cycle-level lane, 2026-11-07. Status: PROVED elementary transfer
lemmas; FINITE-EXACT complete spectra through order nine, combined with those
lemmas to give exact all-size inverse sets for each count from 1 through 593.
No all-order nesting of the fixed-size cycle spectra is claimed here.

## 1. Scope and inheritance

Let c(T) count directed Hamiltonian cycles modulo cyclic rotation, with no
further automorphism quotient. Fixing a vertex as the first vertex counts each
cycle once. We use positive tournament orders and c(T)=0 at orders one and
two. The empty tournament is outside the size domain. Define

- C_n = {c(T): |T|=n}.
- E_n = {e_r(T): |T|=n, r in V(T)}, where e_r counts Hamiltonian paths ending at r.
- N_c(k) = {n>=1: k in C_n}.

Inherited mechanisms: THM-4729's closed-presentation identity is
closed(T)=n c(T)/|Aut(T)|; THM-4741 uses it in the strong-class census.
THM-4862's sole-exit extension is the converse of the endpoint closure below.
The classical unique-cycle and two-cycle classifications (Shelor/Garey;
THM-4877's later two-cycle theorem) already imply their low-level examples.
This note does not claim those results anew.

The closest hostile to path-style padding is a global source or sink: it
preserves H(T), but destroys every Hamiltonian cycle. The least-used sidecar
here is the entire endpoint-count support E_n, rather than H alone.

The independent census never reads an inherited graph bank or census table.
As a separate cross-check only, inherited class files at current main
f9de668e6caa from
04-computation/experiments/tiling_line_metagraph_n9_20261020_classes_nN.tsv
(N=6,7,8; N=9 has suffix .tsv.gz) recover the same supports using
c=closed*Aut/n. The order-nine inherited table has 191536 classes.

## 2. Exact endpoint-cut recurrence (PROVED; elementary)

For a tournament X of order m>=2 let P_X(a,b) be its number of Hamiltonian
paths starting at a and ending at b. Adjoin v with out-neighborhood A and
in-neighborhood V(X)\A. Cutting a Hamiltonian cycle at v gives exactly one
Hamiltonian path of X with first endpoint in A and last endpoint outside A.
Conversely every such path closes through v. Therefore

    c(X+v) = sum_{a in A, b notin A} P_X(a,b).
    C_(m+1) = union_{|X|=m} {sum_{a in A,b notin A}P_X(a,b): A subset V(X)}.

This is an exact all-order recurrence, not a claimed short closed form. The
endpoint matrix is essential: its scalar sum H(X) forgets the endpoint cut.
The boundary m=1 is excluded, since it has no directed two-cycle.

This recurrence supplies an exact algorithm for any fixed n. Strongness is
equivalent to c>0 for n>=3 by Camion; a non-strong tournament contributes only
zero. Do not extend that equivalence to the usual vacuously strong singleton.

## 3. Endpoint closure and count-preserving padding (PROVED)

Let X have order m>=2 and designated vertex r. Adjoin v so that r->v and
v->x for every x!=r. A Hamiltonian cycle must use the unique incoming arc
r->v; cutting there is a bijection with Hamiltonian paths of X ending at r.
Consequently

    c(X+v)=e_r(X),  and indegree(v)=1.

If any tournament T of order at least three has a vertex v of indegree one,
replace v by a transitive chain of any positive length, inheriting all its
external arcs. The unique old in-neighbor is the only outside entry into that
chain. A Hamiltonian cycle must traverse the entire chain consecutively in
its unique transitive order; otherwise it cannot return to unvisited chain
vertices. Contraction and expansion are inverse cycle bijections. Hence

    k in E_m, m>=2  =>  {m+1,m+2,...} subset N_c(k).

The converse outdegree-one version is identical after arc reversal.
This does NOT say every minimum-order witness has such a vertex. In fact
14,17,24 occur in C_7 but not E_6. We explicitly use a separate endpoint
witness and cover any preceding sizes with the complete cycle census.

For k>0 these arbitrarily large realizers remain strong, since they have a
Hamiltonian cycle. This sharply contrasts with fixed H inside strong
tournaments: Busch's exponential lower bound bounds the possible order for
each fixed H. Cycle-preserving padding generally changes H. Ordinary
source/sink path-padding instead leaves the strong class.

A useful path bridge is P_m subset E_(m+1) subset C_(m+2), where P_m is the
Hamiltonian-path count support: append a sink r to X, so e_r=H(X), then close.
All larger orders retain the same cycle count. Thus forbidden path counts
need not be forbidden cycle counts: c=7 already occurs at order seven and
c=21 at order eight.

## 4. Complete fixed-order supports (FINITE-EXACT)

    n=1,2: {0}
    n=3,4: {0,1}
    n=5:   {0,1,2,3}
    n=6:   {0,1,2,3,4,5}
    n=7:   {0,1,...,15,17,24}
    n=8:   {0,1,...,54}
    n=9:   {0,1,...,249} minus
           {215,219,223,226,229,232,233,234,235,236,237,238,
            239,240,241,242,244,245,246,247,248}.

In particular, fixed-order cycle support is neither an odd-number support
nor always an interval. The first noninterval occurs at order seven.
The positive minimum is one at every order >=3, by the closure/padding
construction. The maxima through order nine are 0,0,1,1,3,5,24,54,249.
There is no claim that an order-nine gap remains a gap at larger orders.

The support rows have 1,1,2,2,4,6,18,55,229 values respectively. Orders through
eight are complete isomorphism-class censuses. At order nine the computation
uses all 6880*256 = 1761280 extension presentations of order-eight classes,
with repetitions. Every nine-tournament is covered by deleting any vertex.
Its presentation histogram is NOT an isomorphism-class histogram.

## 5. Complete inverse-size classification for 1<=k<=593

The same enumeration finds

    {0,...,593} subset E_9,   594 notin E_9,   max E_9=661,
    C_n subset E_n for every 3<=n<=9.

These statements are FINITE-EXACT, not proposed all-order endpoint laws.
For each k in C_n the second inclusion supplies a possibly different
order-n endpoint witness. Its closure and padding realize every size >=n+1.
The original order-n cycle witness supplies n itself. The complete smaller
spectra certify minimality. For each remaining 1<=k<=593, the E_9 witness
supplies every size >=10, while the spectra through nine certify nonexistence
there. Thus the following are exact statements about ALL orders:

    N_c(0) = {1,2,3,...}.
    N_c(k) = {m_c(k),m_c(k)+1,...},  1<=k<=593,

where

    m_c(1) = 3;
    m_c(2)=m_c(3) = 5;
    m_c(4)=m_c(5) = 6;
    m_c(k)=7 for k in {6,...,15,17,24};
    m_c(k)=8 for k in {16,18,...,23,25,...,54};
    m_c(k)=9 for 55<=k<=249 except the 21 holes listed above;
    m_c(k)=10 for those 21 holes and every 250<=k<=593.

The artifact inverse_sizes_1_593.json supplies a separate rooted endpoint
witness for every count, together with its certified minimum order. For each
count, verify.py checks all intermediate sizes before its endpoint-certified
tail, rather than silently padding the original minimum-order witness.

The spectrum E_9 additionally certifies minimum order 10 and full tails for
25 counts above 593:

    595,597,599,601,603,605,607,609,613,621,627,629,631,
    633,637,639,641,643,645,649,651,653,657,659,661.

No minimal-order conclusion for the excluded values starting with 594 is
claimed by this packet. Universal realization and explicit all-k bounds are
handled by the independent binary-apex lane, not inferred from this census.

## 6. Reproducibility and validation

Commands, from this directory:

    python3 build_inputs.py
    g++ -O3 -std=c++17 -Wall -Wextra -pedantic census.cpp -o census
    ./census . > census.out
    python3 verify.py > verify.out
    python3 -O verify.py > verify_optimized.out
    cmp verify.out verify_optimized.out
    g++ -O3 -std=c++17 -Wall -Wextra -pedantic independent_labels.cpp -o independent_labels
    ./independent_labels > independent_labels.out
    python3 cut_formula.py > cut_formula.out

build_inputs.py constructs complete classes inductively. Isomorphism-invariant
outdegree/neighbor-color refinement makes ordered cells; exhaustive
permutations inside those cells provide a canonical representative. Every
extension cut is included. Counts reproduce 1,1,2,4,12,56,456,6880 at orders
one through eight. The generator therefore proves coverage inductively,
rather than trusting the known counts alone.

census.cpp uses two independent-initialization subset DPs: one for all
Hamiltonian paths (endpoint counts), one rooted at vertex zero (cycles).
All arithmetic is integer; bounds at these orders fit int. Controls include
H odd and n*c<=H for every presentation. There are no strongness or degree
filters in the enumeration.

verify.py uses uncached rooted permutation DFS, a different counting method,
to replay all 1166 saved cycle/endpoint witnesses and 1779 physical
endpoint-closure/chain-padding instances (three per k). All checks survive
Python optimization. It also checks the precise supports, the endpoint
prefix, the finite C_n subset E_n inclusions, and the inverse-size proof
conditions. An inherited census comparison runs when the optional raw class
files are present; the generator and main census do not depend on them.

independent_labels.cpp enumerates every labeled tournament through order
seven (2097152 at seven) with direct rooted permutation DFS. Its full
supports agree with the class-based DP census, including the order-seven
holes and exceptional maximum24. cut_formula.py independently computes
endpoint matrices by literal permutations and tests every cut over every
class of parent orders two through six.

No sampled row, quotient count, or isomorphism-class count is promoted to a
larger universe. Files in the dirty shared checkout were neither staged nor
modified; publication is left to the session's sole publisher.
