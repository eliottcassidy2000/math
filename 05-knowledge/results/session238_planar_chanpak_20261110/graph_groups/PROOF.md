# Equivariant positive cycle carriers: the odd-group lift and its limits

Status: PROVED elementary synthesis; deterministic finite controls below. No literature-priority claim. Chan–Pak remains OPEN.

## Inheritance and precise meaning of planar

The actual recovered S190 source is `05-knowledge/results/session190_constructive_links_20261104/groups/session190-groups/PROOF.md`, Sections 1–2: finite SO(3) groups, plane point groups, and graph planarity are distinct predicates. Its odd tournament group obstruction makes a faithful SO(3) or O(2) subgroup cyclic. The earlier `thomassen_merge_hypohamiltonian_groups_20261104.md` Section 4 discusses planar Cayley groups. We do not identify either convention with planarity of a conflict graph or tournament. In particular the underlying graph of a tournament of order at least five is nonplanar.

Closest mechanisms: THM-3447 (`01-canon/theorems/THM-3447-regular-group-orbit-tournament-parity-and-involution-ties.md`) supplies the inherited odd-group orientation criterion; THM-5090 and `session237_chanpak_synthesis_20261109/cycle_tilings/REPORT.md` supply the corona construction and positive completion tax. THM-4835 (`euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments`, part 5) already records even-subgroup completion obstruction. THM-4987 (`subgroup-path-orbit-gaps-and-small-symmetry-cost`) supplies the stronger quotient-count exclusion. These are inherited, not new discoveries.

Concept board: faithful group actions; selected versus full cycle carrier; equivariant completion; positive packing tax; fixed-port state versus root motion. Anchor is a legitimate symmetry-preserving graph-to-tournament bridge. Niche is counting its completions. Wildcard is retaining noncyclic odd symmetry to demonstrate that induced carrier planarity has no SO(3) consequence.

## Theorem: free odd-group graph actions lift through the corona construction

Let a finite group A of odd order q act freely on the vertices of a nonempty finite simple graph G, preserving adjacency. Write m=|V(G)|, e=|E(G)| and o=#{v:deg(v) odd}. There is a tournament T with a faithful free A-action on its vertices, and A-equivariantly labelled directed odd cycles C_v and private triangles D_v, such that:

1. Their induced conflict graph is G corona K1, with the given action on both copies of V(G).
2. N=|V(T)|=e+5m+o <= e+6m.
3. Precisely 2^k labelled A-invariant completions of the prescribed arcs exist, where
   k=(binom(N,2)-(2e+6m+o))/q.
4. Every such completion satisfies
   H(T)=I(Omega(T),2) >= 3^m I(G,2/3)+2m.
   More precisely I(Omega(T),x) >= coefficientwise (1+x)^m I(G,x/(1+x))+mx.
5. q divides H(T), and H(T)/q is neither 7 nor 21 (the latter by inherited THM-4987).

The action need not be the full automorphism group of T. No assertion of planarity of T or Omega(T), nor equality of the displayed partition functions, is made.

### Proof

Use one vertex s_f for each graph edge f, and private vertices p_v,q_v,r_v,a_v,b_v, with one further c_v when deg(v) is odd. Let A act on edge names and on private names by its graph action. This action is free: a group element fixing an unordered graph edge either fixes its two endpoints, hence is identity by vertex freeness, or interchanges them. The latter would give an even-order permutation and cannot occur in an odd-order group. Private vertices have trivial stabilizer by hypothesis. Faithfulness follows.

For one representative v in each vertex orbit, choose an arbitrary linear ordering of its incident edges. Transport the order to av by a. Freeness makes this well-defined; no stabilizer-choice compatibility is being silently assumed. Prescribe C_v as the incident edge vertices in that order followed by the private vertices in the displayed order, closing cyclically. Its length is deg(v)+5+(deg(v) mod 2), odd and at least five. Prescribe also r_v -> p_v, making D_v=(p_v,q_v,r_v).

These prescriptions are A-invariant. Distinct C_v share at most one vertex, so their arcs do not conflict. The additional chord is not a cycle edge. Exactly sum_v (length(C_v)+1)=2e+6m+o distinct unordered pairs have prescribed orientations. The selected supports give precisely the stated corona.

On each remaining unordered-pair orbit choose an orientation at one representative and transport it. No group element reverses a pair: such an element would induce a transposition and have even order. Thus transport is consistent. The action on unordered pairs is itself free, by the same fixed-endpoints argument. There are k unassigned pair orbits, each having exactly two choices, independently. This proves both existence and the exact completion count.

The S237 packing injection applies without change. The corona polynomial is the sum over independent S in G of x^|S|(1+x)^(m-|S|). For every v there is an extra triangle in T[C_v] through a_v: a shortest directed cycle through a_v has length three, since either orientation of its first chord shortens a cycle of length above three. Such a triangle is neither D_v nor any selected long cycle. Different v give different triangles because their cycle supports meet in at most one vertex. These m additional singleton packings give mx. They can even be chosen equivariantly: choose at vertex-orbit representatives and transport, again using trivial stabilizers. All other contributions are nonnegative. Evaluation at x=2 proves the bound. Finally A acts freely on directed spanning paths, since fixing their ordered vertex list fixes every vertex. This gives q divisibility; apply THM-4987 to the actual subgroup A for the two quotient exclusions. QED.

## Sharp failure boundary and the planar-group specialization

An even-order group cannot act faithfully by automorphisms on any tournament: by Cauchy it has an involution, whose nontrivial vertex permutation swaps some pair and reverses its arc. Thus the odd-order hypothesis is sharp for the class of free actions, not a computational accident. Taking A=C_q for any odd q gives a literal cyclic geometric-group lift. For example, for odd q >= 3, a graph cycle on q vertices, with cyclic translations, lifts while retaining its rotations. Its reflection cannot lift as an ordinary tournament automorphism. Signed converse symmetries are a different problem and are not ruled out by this statement.

Conversely noncyclic odd groups also qualify. An edgeless graph on the regular A-set lifts for A=C3 x C3 (and for any finite odd A), despite the absence of a faithful SO(3) realization of that group. The carrier can even be planar. This is a direct counterexample to the inference “planar selected conflict graph implies planar geometric tournament symmetry.” There is no conflict with the original meaning of planar Cayley group: this theorem does not supply a planar Cayley graph of A.

Freeness is a sufficient condition, not a necessary condition for arbitrary equivariant realization. Without it, incident edge orders cannot in general be transported independently of stabilizer representatives. No classification for nonfree graph actions is claimed.

## What survives towards Chan–Pak, and what does not

This is a structural lift, not a value-preserving reduction. The exact missing sidecar remains the extra cycles B and their neighbor unions in the selected carrier F:

    I(Omega,x) = sum_{R independent in B} x^|R| I(F-N_F(R),x).

Every term is positive at x>0. Keeping the A-action does not cancel that tax. Dividing H by q counts genuine path orbits; dividing a partition function by an arbitrary numerical content does not manufacture a group action. The four-coordinate rooted state (h,R,e,d) from S237 closes fixed-port S/B tails but loses root motion and converse information; its equal-state hostile pair remains valid. No assertion here replaces the full root deck with a group or with a selected carrier. S236 port cones likewise concern a growth semigroup, not invertible automorphisms of a fixed tournament.

Connection ledger: source is (G, free odd A-action); target is (T, selected corona, A-action); map is equivariant partial orientation and completion; preserved predicates are induced conflicts and actual subgroup action; destroyed information when only G is kept is the remaining orientation bits and all extra-cycle packing contributions; restoration sidecar is the full cycle conflict data. Cheapest failure controls are an even regular action and a single selected long cycle, whose private completion already adds cycles.

## Verification and scope

Run `python3 04-computation/experiments/session238_graph_groups/verify.py` and the same with `python3 -O`. The verifier enumerates every undirected circulant connection set on C1,C3,C5,C7 (15 graphs), plus the empty and complete regular C3 x C3 graphs. It checks the literal arcs, all group automorphisms, free actions, induced corona supports, forced-pair and orbit counts, and existence of the extra triangles. Independent Hamilton-path DP checks the C1 and edgeless C3 cases. The `--inject-failure` control changes an output arc and must fail under both modes. No unbounded result depends on finite tests.

Targeted canon search found inherited orientation and completion obstructions, credited above; no external novelty claim is made for this combination with the S237 carrier. No count-preserving domination reduction, complete-spectrum proof, or full-automorphism-group classification results.
