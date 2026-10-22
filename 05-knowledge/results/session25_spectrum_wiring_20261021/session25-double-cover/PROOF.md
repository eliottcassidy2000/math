# Orientation character double cover of the fundamental merged line graph

Status: independently proved and finite-exact tested. No literature priority claim.

## Carrier and statement

Fix n >= 3 and labels 0,...,n-1. Let L have all labelled tournaments as vertices. Its edges are unordered pairs {T,U} such that U is obtained from T by retaining the arcs of a directed Hamiltonian path P and reversing every other arc. Write this operation kappa_P. Set N=binom(n,2), m=binom(n-1,2). Here m>0, so T != U. Moreover the common directed arcs of T and U are exactly P, making the path unique for an edge. Thus this definition does not forget a parallel path witness.

Let Gamma=S_n x C2 act by relabelling and optional converse. Set eta(T)=(-1)^(number of backward arcs relative to the fixed label order), chi(pi,epsilon)=sgn(pi)(-1)^(epsilon N), K=ker chi. The character is surjective since a label transposition has character -1.

Let F=L/Gamma be the quotient graph with darts retained: one edge per Gamma edge orbit, one vertex per Gamma tournament orbit, and two incidences for every loop. Let Q=L/K. Then Q -> F is an honest, unramified, two-sheet graph covering, with the two vertices above each base vertex distinguished by eta. Each edge has constant sheet voltage (-1)^m. Consequently:

* If n=0 or 3 mod 4, Q is the canonical bipartite double cover of F.
* If n=1 or 2 mod 4, Q is the disjoint union of two copies of F.

The full fundamental blue/black multiplicities are retained. In the odd-voltage case a base loop lifts to TWO parallel cross-sheet edges, not one. In the even-voltage case it lifts to one loop on each sheet.

## Proof

### 1. Character and vertex fibres

Relabelling reverses the reference direction on an odd number of unordered pairs precisely when the permutation is odd. Hence eta(gT)=chi(g)eta(T). Converse reverses all N pairs. A complement move reverses exactly m pairs, hence eta(kappa_P T)=(-1)^m eta(T).

Every tournament stabilizer in Gamma lies in K. Each Gamma orbit therefore consists of exactly two K orbits: the eta-positive and eta-negative subsets. Indeed, if gT and hT have the same eta, chi(h^{-1}g)=1, and if chi(g)=-1 the sign changes.

### 2. No edge inversions

Suppose g exchanges T and U=kappa_P T. It preserves their common directed path P as a directed arc set under the signed action.

If g is ordinary relabelling pi, it is an automorphism of the directed path, so pi is the identity. This cannot exchange distinct tournaments.

If g includes converse, pi must reverse the path order. Write P=(p_0,...,p_{n-1}). The pair {p_0,p_{n-1}} is not a path edge because n>=3. Relabelling by path reversal exchanges its endpoints, and converse reverses its direction again: their composite preserves the actual directed arc on this pair. But kappa_P reverses that arc. Thus gT cannot equal U.

Therefore every edge stabilizer fixes each endpoint, and hence lies in K. Every Gamma edge orbit splits into exactly two K edge orbits. This includes edge orbits whose endpoints become a base loop.

### 3. Covering at darts, including stabilizer issues

The quotient definitions use darts (e,v), where v is one of the two distinct labelled endpoints of e. No Gamma element exchanges the two darts of an edge, by Step 2. Consequently every base edge has exactly two dart orbits, even when its endpoints have the same base vertex.

For each base dart orbit d and a chosen K vertex above its initial base vertex, there is exactly one K dart orbit above d with that initial K vertex. Existence follows by taking a representative dart and acting by an element of the appropriate chi sign. For uniqueness, suppose two representatives of d have initial endpoints in the same K vertex orbit. First act on one by K to align their initial endpoints literally. An element of Gamma sending the first aligned dart to the second now fixes that endpoint, so lies in K. The two darts are therefore already K-equivalent.

This proves the local-star bijection; it does not assume tournament automorphism groups are trivial. It also proves that the residual Gamma/K action on Q is the free deck involution. Step 1 identifies its voltage, giving the asserted isomorphisms.

### 4. Why the base is the fundamental colored multigraph

Normalize each directed path to the fixed ordered base path. An S_n orbit of an edge is exactly a tiling line {t,kappa t}: there is no nonidentity ordinary path automorphism. The remaining signed operation is reversal of the path followed by converse, which acts on tilings as staircase reflection sigma. Thus F edges are precisely <kappa,sigma>-orbits, not raw tiling lines with merely merged vertices.

The apex pair joining the two path endpoints is fixed by sigma and toggled by kappa, so kappa sigma fixes no tiling. A line fixed by sigma is therefore pointwise sigma-fixed and is blue. A black line has its distinct sigma twin. Thus the Gamma quotient preserves each blue line and takes ONE edge for each black twin pair, exactly the fundamental black multiplicity convention of THM-4731. Colors are invariant under Gamma and pass through the covering unchanged.

## Exact one-bit holonomy

A k-step labelled complement walk from T to gamma T satisfies

    chi(gamma) = (-1)^(km).

Equivalently any closed walk in F has sheet holonomy (-1)^(km). When m is odd, odd closed walks exchange sheets and even ones preserve them; when m is even all closed walks preserve the sheet. This is only a character obstruction on the closing signed relabelling. It does not recover the permutation, the converse bit separately, or literal labelled identity holonomy. The character is independent of the choice of closing gamma because stabilizers lie in K.

For each connected base component, the odd-voltage cover is connected exactly when that base component is nonbipartite. If the base component is bipartite, the cover has two components, each isomorphic to it. Loops count as odd cycles, as required.

## Bridge transfer (exact, conditional on actual base wiring)

For m even, each lift of e is a bridge iff e is a base bridge.

For m odd, a lift of e is a bridge iff e is a base bridge AND at least one of the two components of the original connected component after removing e is bipartite. This formulation excludes unrelated connected components of F.

Proof: if e is not a bridge it lies on a cycle (a loop is a length-one cycle, parallel edges give length two). An even cycle lifts to two cycles and an odd cycle lifts to a doubled cycle, so neither lift is a bridge. If e is a bridge with sides A,B, its two lifts are the only two edges between the lifted sides. A connected side has connected canonical double cover exactly when it is nonbipartite. If both A and B are nonbipartite, their connected covers are joined by the two lifted edges: together the two edges are a cut, but neither alone is a bridge. If a side is bipartite, its two cover components each meet exactly one lifted edge, forcing both lifted edges to be bridges. This proves both directions with loops and multiplicities included.

No claim here classifies the actual bridges of the tournament base graph.

## Hostile controls and failure boundaries

* n>=3 is load-bearing. For n=2, m=0 and the purported labelled edge has T=U. Its path endpoint pair is itself the path edge, so the inversion proof and the simple labelled-edge model do not apply.
* Collapsing the two parallel edges above an odd-voltage base loop destroys the covering's local-star bijection. The n=4 census contains such a loop.
* Raw black multiplicities are twice fundamental black multiplicities, so using them gives the wrong base edge count.
* Vertex stabilizers need not be trivial; the dart proof explicitly retains their quotient effect.
* Two odd cycles joined by a bridge are a hostile to 'all base bridges lift to bridges': neither lift of that joining edge is a bridge.
* A labelled cover of L itself is NOT asserted. This is L/K -> L/Gamma, and its vertices already forget all K relabellings.

## Independent exact audit

Run `python3 verify.py` in this directory. No third-party packages are needed. It constructs all labelled tournament vertices, all Hamiltonian-path complement edges, and all signed relabellings independently of any repository engine. It computes Gamma/K vertex, edge, and dart orbit representatives, checks character covariance on every tournament/group pair, searches every edge/group pair for inversions, verifies exactly two orbit lifts, checks the local-star bijection at every base dart and lifted vertex, verifies sheet voltage, and verifies blue/fundamental-black counts directly.

Universe: n=3,4,5, all labelled tournaments with no score, regularity, self-converse, or other filter. Bridge transfer is separately checked for every graph on 1..4 fixed vertices with arbitrary loops and simple nonloop support, plus a second copy of each present edge as parallel-edge hostile control.

Counts (n; labelled V,E; base V,E; cover V,E; base loops; blue,fundamental black):

* 3; 8,6; 2,1; 4,2; 0; 1,0.
* 4; 64,96; 3,3; 6,6; 1; 2,1.
* 5; 1024,3840; 10,20; 20,40; 2; 8,12.

The bridge audit covers 2,196 graph cases and 15,975 base-edge checks, with zero failures. See verify.out for the full machine output. These finite tests supplement the general proof; they are not its justification.

## Targeted overlap audit and loss ledger

* THM-4731, `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`, is the closest inherited quotient mechanism: sigma twins and fundamental black multiplicity. It does not state this orientation-character cover.
* THM-4524, `01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md`, concerns switching gauges, Hamiltonian-path arc parity, and orientation characters for Euler-graph automorphisms. Its orientation-sign mechanism is related, but its carrier and theorem differ.
* THM-4464, `01-canon/theorems/THM-4464-checksum-projective-tournament-and-output-sheet-cocycle.md`, concerns a checksum/projective-direction output sheet, not the complement graph quotient.
* THM-4532, `01-canon/theorems/THM-4532-tournament-clocks-two-sheet-clocks-all-odd-tournaments-and-the-syracuse-clock.md`, concerns two-sheet clock tournaments and Hamiltonian-path arc counts, not this graph cover.
* `05-knowledge/results/metagraph_spectrum_torus_session16_20261020.md` records two orientation-parity-separated icosahedra for the five-vertex regular-return graph. That is a special labelled two-step/regularity carrier, not the quotient covering theorem above.

Source L/K; target L/Gamma; map forget eta; preserved predicate: edge adjacency, colored fundamental edge multiplicities, local incidences, one-bit holonomy, and the stated conditional bridge law; lost information: exact tournament labels, path labels, full relabelling/converse holonomy, and unquotiented automorphism action. Restoration sidecar for this map only: the eta sheet. Cheapest hostile is the n=4 base loop and its two parallel cross-sheet lifts. No literature-priority or full bridge-classification claim is made.
