# Session 32 independent graph and transfer audit

Date: 2026-10-23 (UTC). Status: **PASS for the mathematical claims in the audited scope.**

The independent auditor imports no producer code. It constructs square-sum
edges by reflection in each square, reconstructs connected components by
union-find, establishes the nontrivial small path censuses with degree-obligation
edge search rather than vertex-order DFS, and glues block states in the full
25-vertex graph rather than the producer's virtual graph. Its output is
`certificate.json`; running with `python3 -O` produces an identical certificate.
All checks use explicit exceptions. The certificate pins the source documents,
programs and producer certificates by SHA-256.

## Findings at a glance

- The arrival formula, all-order degree-two classification, connectivity
  threshold, component/cycle response, endpoint-aware chain suppression and
  marked-ear correspondence are correct under the stated simple-graph and
  distinct-endpoint conventions.
- The endpoint-aware obstruction at 24 is complete. The transition to 25 is
  a compulsory marked ear, not an arrival-degree threshold.
- The local labelled signature counts are **5, 101, 149**. Forgetting endpoint
  identities while retaining their allocation and terminal matching gives
  **3, 89, 131**, with multiplicities as recorded by the producer.
- There are **1543** endpoint-budget triples and **19** degree-feasible gluings.
  Exactly **9** are disconnected; the remaining **10** are precisely the
  Hamiltonian paths, up to reversal.
- Contracting the three entire trees gives a connected 13-vertex, 20-edge
  multigraph of cycle rank 8. It has no Hamiltonian path: deleting the three
  supervertices leaves six components. Simplification removes the parallel
  pairs at C–15 and C–19 and lowers the rank to 6. The counterexample is valid.
- The integral cat-map conjugacy is exact. The rooted-forest recurrence and
  block polynomials are correct, with a genuine distinction between boundary
  matrix entries and traces.
- Gerbicz's all-order existence theorem is correctly attributed as **CITED**.
  This audit does not verify his large finite base-pair or gluing tables.

## 1. Arrival arithmetic and persistence

For the new vertex n, older neighbours correspond bijectively to squares q
with n < q <= 2n−1. Hence its degree is

    floor(sqrt(2n−1)) − floor(sqrt(n)).

The subtraction of a possible loop is therefore already built into the upper
bound. The full old-vertex degree formula is independently checked too.

For r=floor(sqrt(n)), degree two is equivalent to

    max(r², ceil(((r+2)²+1)/2)) <= n
      <= min((r+1)²−1, floor((r+3)²/2)).

The seven bands r=1,...,7 yield [13,15], [19,24], [25,32], [36,40], [49,50].
For r>=8, the polynomial r²−6r−10 is positive and increasing, so
2r²−1 >= (r+3)² and arrival degree is at least three. This proves the global
classification; checking finitely many large n is not being substituted for it.

Let r_n count the distinct old components meeting the new neighbourhood.
Adjoining the vertex first increases component count by one; its edges then
make precisely r_n mergers. Thus c_n=c_(n−1)+1−r_n and
beta_n=beta_(n−1)+d_new(n)−r_n. The square-sum producer ledger through 36
agrees entry-for-entry with independently reconstructed data. Supplementary
controls run through 64, not as an infinite proof.

The components at 12 have vertex sets of sizes 5,3,4. The new vertex 13 joins
the 5- and 4-vertex components and creates a degree-three branch at 3; 14 joins
the remaining component. The tree at 14 has three leaves. Vertex 15 then
adds one cycle rather than joining different components. These are genuinely
different structural events despite identical arrival degree two.

For n>=6, the first square above n is less than 2n; n=5 is checked directly.
Every new vertex from 5 onward therefore has an older neighbour. Once G14 is
connected, connectivity persists. The complete earlier graphs show that only
G1 is connected below 14.

## 2. Exact arrival response and endpoint-aware suppression

For an arrival vertex v with two distinct neighbours a,b, removing v from a
Hamiltonian path has two disjoint possibilities. An endpoint v leaves a
Hamiltonian path of the old graph with an exposed endpoint a or b. An internal
v leaves two spanning path components, with a and b exposed on different
components. Singleton components are allowed and cause no exception. These
operations have unique inverses, so the correspondence is exact.

If the old graph has no Hamiltonian path, the endpoint case is impossible.
Suppressing a–v–b then requires a distinguishable marked a–b edge. If a–b
already exists, its marked parallel copy must remain distinct. Replacing this
marked multigraph by a simple unmarked square-sum graph would not be valid.

For the chain kernel, put both proposed endpoints into the boundary, together
with every original vertex of degree other than two. Every internal chain
vertex is now certified to be a nonendpoint. All long chains are mandatory;
direct boundary edges are optional. Retaining loops, parallel edges and the
full chain records, the prescribed boundary degrees plus connectivity give a
connected spanning degree-at-most-two graph with exactly two leaves. It is
therefore a Hamiltonian path. Conversely every Hamiltonian path restricts to
one such selection. This establishes both completeness and uniqueness of
expansion. A pure cycle has its two chosen endpoints in the boundary, so is
not an omitted exception.

The certificate contains a G23 path ending at graph-degree-two vertex 22.
The G25 path ending at graph-degree-two vertex 8 is another direct witness
against unmarked suppression of all degree-two vertices.

## 3. The 24 obstruction and the 25 ear

The proof at 24 survives a case-by-case logical audit:

1. The leaf 18 fixes one endpoint. If the other endpoint e is not 9,
   the interior vertex 9 and leaf 18 saturate vertex 7, excluding 2–7.
   If e is neither 11 nor 22, those two interior vertices saturate 14,
   excluding 2–14. Unless e is 2, this would leave the interior vertex 2
   with only neighbour 23. Thus e belongs to {2,9,11,22}.
2. Vertex 20 is consequently interior and forces 5–20. If e is not 11,
   interior 11 also forces 5–11. If e=11, saturation at 7 forces the
   interior vertex 2 to use 14 and 23, while interior 22 uses 14. Vertex
   14 then cannot meet 11, so endpoint 11 must meet 5. In all cases 4–5
   is excluded.
3. Every vertex of 1–8–17–19–6–10–15–21–4–12–24–1 is interior. The listed
   degree-two vertices force every cycle edge except 4–12. Interior 4
   cannot use 5, so it forces 4–12 as well. A proper cycle cannot belong
   to a Hamiltonian path.

At 25 the new neighbours are 11 and 24. Nontraceability of G24 rules out an
endpoint 25, so the marked ear 11–25–24 is compulsory. The pair of square sums
changes from {25,36} at arrival 24 to {36,49} at arrival 25, without changing
arrival degree.

The independent degree-obligation search recovers the complete path sets for
n=14,15,16,17,23,24,25, using only 78 search nodes across all allowed endpoint
pairs. It confirms counts 0,1,1,1,3,0,10, respectively. All listed producer
paths through 25 are checked individually. For the other failures, explicit
component/leaf/cut certificates suffice; no redundant large census is run.

## 4. Block transfer and endpoint anonymity

Each local spanning linear forest component has m vertices and m−1 internal
edges. If it receives p selected ports and e global endpoint marks, its degree
sum gives 2(m−1)+p=2m−e, hence p+e=2. This includes a singleton component:
it can have two different ports or one port and its endpoint mark. A
component with no selected port cannot connect to the nonempty core.

The complete interface is the set of named selected boundary edges and their
pairing by components, with endpoint terminals. Retaining realization
multiplicities gives a bijection at the level of realized states, not merely a
surjective feasibility test. Expansion preserves degree and connectivity.
The producer added the explicit graph-order-at-least-two condition after this
audit identified the harmless singleton boundary case; the singleton path is
now stated separately.

For unrestricted endpoint existence, the identity of a marked endpoint
inside a block can be discarded: subsequent gluing sees only which named
port reaches that leaf. Its allocation and terminal pairing cannot be
discarded. For endpoint-conditioned questions, its identity must remain.
The anonymous multiplicity histograms are independently reconstructed:

- A: one signature of weight 1 and two of weight 2.
- B: 81 of weight 1, four of weight 2, four of weight 3.
- C: 118 of weight 1, ten of weight 2, two of weight 3, one of weight 5.

The five labelled A states and three anonymous states in the prose are exact.
The general distinction between pairing and degree data is also valid: the
four-port continuation can be either a spanning path or an edge plus a
separate cycle with the same port degrees. Global connectivity is essential.
The nine rejected degree-feasible G25 assemblies independently demonstrate
that issue in the actual graph.

## 5. Forest recurrence and integral conjugacy

A useful independent all-order derivation avoids a determinant recurrence
being guessed from its first terms. A spanning forest of P_m is a composition
m=l_1+...+l_c into consecutive component lengths. A choice of one root in each
component has weight product(l_i), and there are m−c edges. A single component
has generating series

    C(z)=sum_(l>=1) l t^(l−1) z^l = z/(1−tz)^2.

Nonempty forests are ordered sequences of these components, so

    sum_(m>=1) f_m(t) z^m
      = C(z)/(1−C(z))
      = z/[1−(1+2t)z+t²z²].

Equivalently,

    [t^k] f_m(t) = binomial(2m−k−1,k),  0<=k<=m−1.

Thus f_1=1, f_2=1+2t and
f_m=(1+2t)f_(m−1)−t²f_(m−2) for m>=3. A formal zero coefficient f_0=0
belongs to this nonempty generating function; it is not the forest count of
the empty graph, whose value is one.

The auditor independently sums rooted forest weights over every edge subset
for path orders 1,...,12. It also evaluates det(I+tL) by fraction-free integer
elimination at m distinct t-values for each such order. Since both sides have
degree at most m−1, those evaluations establish full polynomial equality in
each checked order, not just agreement at t=1. The composition argument proves
the recurrence for every order.

With A=[[2,1],[1,1]], B=[[3,−1],[1,0]], U=[[1,−1],[0,1]], direct multiplication
gives AU=UB=[[2,−1],[1,0]], and det(U)=1. The conjugacy is over SL(2,Z).
At t=1, f_m=F_(2m)=(B^(m−1))_11, whereas tr(A^m)=L_(2m). These observables
are not interchangeable: tr(A²)=7 while f_4=21.

The three G25 block polynomials are exactly

    A: 1+8t+20t²+18t³+5t⁴,
    B: 1+6t+10t²+4t³,
    C: 1+10t+36t²+56t³+35t⁴+6t⁵.

Their values are 52,21,144. The branched four-edge A block is not equivalent to
a four-edge path, whose value is 55. At 12 the genuine components P3,P4,P5
give 8·21·55=9240. Multiplying the three selected G25 block values would not
compute the full graph's forest count because boundary and core edges remain.

This audit does not validate FGST representation-theoretic claims or turn
any of these identities into a Monster, Collatz, or tournament-gap theorem.

## 6. External theorem and citation boundary

The [OEIS A090461 entry](https://oeis.org/A090461), read directly during this
audit, contains Gerbicz's January 21, 2018 attribution of Hamiltonian paths
for every n>=25 and cycles for every n>=32. The
[primary Mersenneforum discussion](https://www.mersenneforum.org/node/17331/page2)
was unavailable through direct opening but its indexed primary-author text
was accessible. It describes the consecutive-length pair invariant, shared
position parities and square-scaled shifted blocks, and distinguishes the
25 and 49 moduli. This supports the manuscript's cited scope. Neither its
large base tables nor all gluing entries were independently verified here.

The finite work in this packet therefore does not establish the cited
infinite theorem on its own. Ordinary interval Hamiltonian-path existence is
not being presented as an open conjecture.

## 7. Editorial and execution notes

The initial square-sum README said that all paths, endpoints, cut components
and arrival data were output through 36. The actual path census is through
25; only the arrival ledger runs through 36. This is a scope-wording repair,
not a mathematical error. The coordinator has corrected that wording; the
final README was reread and its updated hash is included in the certificate.

The initial `check_catmap_forests.py` uses removable Python assertions. It is
valid under ordinary Python, but running it under `-O` would disable those
checks. This independent audit supplies explicit checks that remain active
in optimized mode and does not claim that an optimized producer run tested
its assertions.

No producer source was edited by this auditor. No Git operation or repository
mutation was performed. The audit writes only its own packet.

## Reproduction

From the published packet root:

    python3 graph-audit/audit_graphs.py > graph-audit/certificate.json
    python3 -O graph-audit/audit_graphs.py > graph-audit/certificate.optimized.json
    cmp graph-audit/certificate.json graph-audit/certificate.optimized.json

The published script reads `square-sum/`, `ports/` and root-level producer files. Its recorded source pins preserve the original audit provenance; publication replay verifies the remapped files separately.
