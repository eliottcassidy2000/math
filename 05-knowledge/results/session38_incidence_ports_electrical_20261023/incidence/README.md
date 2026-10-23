# Twelve incidence symbols determine the order-eight port table

October 23, 2026. **PROVED:** incidence-balance and local transport identities,
conditional class-to-port determination, and the stated family-specific
minimal-alphabet implication. **FINITE-EXACT:** all 6,880 ordinary order-eight
classes and the supplied exhaustive 2^21 fixed-path presentation ledger.
**REFUTED:** completeness of the exact triangle-edge signature at every order.
**OPEN:** absolute smallest sufficient incidence signature or an all-size
bounded-profile determination rule for the metagraph.

## Inheritance and scope

The closest mechanism is `metagraph_bucket_obstruction_20261023`: J forces
order six with endpoint/color capacities; J4 separates all order-seven
classes but only 6,878 of 6,880 order-eight classes. The canonical hostile is
its reflected 2-by-2 switch preserving weak port margins. The corrected near
miss is equating a class census with replay of the actual 2^21 order-eight
ports. The least-used coordinate here is the cyclic-triangle multiplicity of
an actual outgoing edge, attached to the cyclic-triangle count of its head.

Live concepts: cyclic edge flow; associated measurements; modular alphabets;
capacity peeling; authenticated source paths; regular hostile families.
The target remains the converse-merged colored complement-line multigraph,
not the arc-flip graph, an uncolored support graph, or a loop-free graph.

## 1. A smaller local signature

For a tournament T, write c_v for the number of directed cyclic triangles
through vertex v and c_uv for the number through the unordered pair {u,v}.
Let braces below denote multisets, preserving repeated entries.

    E(T) = { { (c_w, c_vw) : v -> w } : v in V(T) },
    E_pq(T) = { { (c_w mod p, c_vw mod q) : v -> w } : v in V(T) }.

These are isomorphism invariants. No vertex labels are retained, but the two
measurements belonging to each outgoing edge remain paired, and each list
remains attached to its root before the outer multiset is taken. The exact
E recovers its root score as the row length and its root c_v as the sum of
second entries, because

    sum_{v->w} c_vw = c_v = sum_{w->v} c_vw.                 (1)

Each cyclic triangle through v has exactly one outgoing and one incoming
edge at v, proving (1). Thus an explicit root (score,c_v) tag is redundant
for E. Root c_v is NOT claimed recoverable after modular reduction.

The oriented weighted triangle graph with arc weight c_uv is Eulerian:
its incoming and outgoing weighted degrees are c_v. This is a divergence-
free incidence flow on the tournament's vertex set. It is not the electrical
Laplacian on the metagraph's tournament-class vertex set. Equating those two
carriers would discard the source/path/target map.

All counts can be computed by one O(n^3) triple pass, then sorting. Strong-
four enumeration, Hamiltonian-path counts, automorphisms, and a second
neighborhood-refinement round are unnecessary to compute E or E_43.

## 2. Exact finite boundary and limited minimality

The complete order-eight class census gives:

| Signature | Distinct classes |
|---|---:|
| J | 6,719 |
| (J,M2) | 6,853 |
| J4 | 6,878 |
| J4 plus entire unassociated c_uv histogram | 6,878 |
| E, exact paired edge data | 6,880 |
| E_43, a twelve-symbol outgoing-edge alphabet | 6,880 |
| exact c_w, but c_vw mod 2 | 6,874 |
| c_w mod 2, but exact c_vw | 6,863 |
| c_w mod 3, but exact c_vw | 6,876 |
| E with the c_w field removed | 6,749 |
| E with the c_vw field removed | 6,535 |
| E with two separate marginal lists at each root | 6,871 |
| global multiset of (c_u,c_v,c_uv) over arcs | 6,832 |

The unpaired-root signature keeps `{c_w:v->w}` and `{c_vw:v->w}` but forgets
which entries belong together. Its nine colliding pairs show that merely
retaining both marginal measurements does not retain their association.
The global edge signature loses which outgoing edges share a root.

The checker tests p=1..8 and q=1..7. E_43 is the unique coordinatewise-minimal
successful pair in this tested grid. More strongly, twelve is the minimum
alphabet size p*q for class separation within the entire positive-integer
modular family E_pq: p=1,2,3 fail even with exact edge counts (p=1 is a
coarsening of p=2), and q=1,2 fail even with exact head counts (q=1 is a
coarsening of q=2). Any product smaller than twelve has p<=3 or q<=2. E_43
attains twelve. This does NOT prove absolute minimality among arbitrary
signatures, arbitrary encodings, or arbitrary ways of choosing a sidecar.

Completeness at eight implies completeness at every lower positive order
for this same signature family. Adding one universal sink appends one
(0,0) entry to every old row and adds one empty row. It preserves signature
collisions. The new universal sink is unique, so isomorphism after padding
would imply isomorphism before padding. A smaller-order counterexample
could therefore be padded into an order-eight counterexample.

The finite census and port table use ascending numeric canonical words as
class IDs. Words are big-endian upper triangles in lexicographic i<j order,
with vertices ordered by ascending score and minimized within score blocks.
All literal collision classes and menu counts are in `receipt.json`.
`signature_ids.tsv` assigns exact first-occurrence IDs, not collision-prone
numerical hashes. Independent canonical/orbit-mass and full-mask verification
is in the separate Session38 audit packet.

## 3. A scalar repair if J4 is already retained

Define

    omega_v = sum_{v->w} binom(c_vw,2),
    Q(T) = sum_{u->v} omega_u * omega_v^2.

An omega_v term counts a strong-four subtournament whose two cyclic
triangles share an edge directed out of v, with v the tail of that unique
shared edge. Thus Q is an incidence-associated rooted-overlap moment, not
an unassociated list of motif totals.

The two surviving J4 pairs have canonical words and Q values:

    (6427778,6427781): (2881,2887),
    (27467784,27467793): (3008,3000).

Consequently (J4,Q) separates all 6,880 classes. The complete unassociated
edge histogram fails on both pairs, so retaining all its scalar moments
also fails. Q is a one-scalar repair for this particular boundary, not a
minimal bit-encoding theorem. E_43 is a cheaper structural alternative that
does not first compute J4.

## 4. Actual order-eight wiring, not only node identification

The supplied independent constructor exhausts all 2^21 fixed-path masks,
with descending path 7->...->0. Its closed-half table has total mass 2^20
and 956,339 nonzero source/target/endpoint/color entries. Our checker consumes
every one of these entries and every open endpoint capacity.

For any transported target invariant S, a row is

    (source class i, source start/end scores a,b, color beta, S(target)),

with capacity equal to its number of actual closed source ports. A column is

    (target class j, target start/end scores 8-a,6-b, color beta),

with capacity equal to its number of open target ports. An entry is allowed
exactly when its target signature, endpoint type, and color match. The
actual table supplies feasibility. Class-separating S makes each row have
one target; hence all entries are forced. Aggregating class IDs into their
converse pairs recovers the colored merged table, keeping loops and
multiplicities. This statement assumes actual transported target signatures,
not just source class signatures or target node margins.

Full replay results:

| Signature | Rows | Candidate entries | Unresolved mass after peeling |
|---|---:|---:|---:|
| J | 947,804 | 992,343 | 56,568 |
| (J,M2) | 954,479 | 961,665 | 10,780 |
| J4 | 956,171 | 956,541 | 884 |
| (J4,Q) | 956,339 | 956,339 | 0 |
| E | 956,339 | 956,339 | 0 |
| E_43 | 956,339 | 956,339 | 0 |
| E with per-root separate marginals | 955,752 | 958,189 | 4,026 |
| exact head, parity edge | 956,122 | 958,266 | 2,568 |
| parity head, exact edge | 955,547 | 960,192 | 6,262 |
| head mod 3, exact edge | 956,146 | 957,016 | 1,351 |

Peeling removes only entries forced by a sole remaining row or column
neighbor. Exhaustion proves uniqueness; failure to exhaust does NOT prove
another feasible table or another merged multigraph. The residual files
retain complete capacities, allowed cells, and native feasible entries for
the stronger residual-cycle/projected-wiring criterion in the root packet.
Family-minimal CLASS separation above must not be silently substituted for
family-minimal WIRING without that additional certificate.

## 5. All-size local incidence transport along an actual path

Normalize the directed Hamiltonian path as P=(0->1->...->n-1). Let K_P retain
path arcs and reverse every other arc. For a triple B let C_B be its cyclic
indicator, and let Delta_B=C_B(K_P(T))-C_B(T).

* If B contains no path edge, Delta_B=0: all three arcs reverse.
* If B contains two path edges, its positions are consecutive and
  Delta_B=1-2C_B: its only free diagonal reverses.
* If B contains exactly one path edge u->v and third vertex w, then

      Delta_B = A(u,w)-A(v,w).                          (2)

For (2), C_B=A(v,w)A(w,u), whereas C'_B=A(w,v)A(u,w).
Subtracting gives A(u,w)-A(v,w). Therefore

    c'_e = c_e + sum_{B containing e} Delta_B,
    c'_v = c_v + sum_{B containing v} Delta_B.           (3)

Exactly (n-2)^2 triples contain a path edge: (n-1)(n-2) path-edge/triple
incidences, minus the n-2 triples counted twice. Thus after initial triangle
incidences, an actual path complement updates all c_e,c_v in O(n^2) time.
These are all-order identities, not an extrapolation from the finite replay.

The original adjacency and actual path are required in (2). E(T) alone does
not determine E(K_P(T)). Already at n=4, on ascending path 0->1->2->3 and
free order [(0,2),(0,3),(1,3)], masks 0 and 1 have the same source E but
different partner E. Their source outgoing bit rows are [2,4,9,3] and
[6,4,8,3]. The first partner is transitive and the second has a cyclic
triangle. The structural receipt retains their complete partner signatures.
This is the smallest possible order for this particular path-dependence
failure: orders two and three have at most one path orbit per ordinary class.

## 6. A concrete all-order obstruction

On Z/15, orient x->x+d for d in either connection set

    S1={14,13,12,11,5,6,8},
    S2={14,13,12,4,10,9,8}.

Both are strongly connected because 14 generates Z/15. Every vertex in each
has score 7, c_v=28, and outgoing edge-triangle counts

    (3,3,3,4,4,5,6).

Hence their exact E, every E_pq, and all E-based weighted neighbor-refinement
rounds agree. Their rooted strong-four counts are all 184, so every ordinary
J4 neighborhood-refinement round agrees too. Yet exact subset DP gives

    H(T1)=196264905, H(T2)=196277235.

Unequal Hamiltonian-path counts prove nonisomorphism. This is an explicit
order-fifteen failure, not a claim that fifteen is the first failure.
Adding universal sinks preserves E and H, so collisions exist at every
order n>=15. The root packet gives a strong-lift continuation using the
additional independently checked ordered two-path partition totals

    F2(T1)=4857720210, F2(T2)=4857619470.

Class indistinguishability is not itself a counterfeit metagraph: the two
nodes have different H and further actual port capacities can distinguish
wiring. The conclusion is that this local signature is not an all-order
class decoder, nor can merely iterating its uniform neighbor summaries fix
this hostile. Exact path-labelled curls remain an all-order port decoder,
as proved in `metagraph_port_rigidity_20261023` and `metagraph_two_height_20261023`.

## 7. Reproduction

From this packet, using the final adjacent wiring inputs or explicit paths:

    python -B check.py --classes CLASSES --ports PORTS --capacities CAPACITIES
    python -O -B check.py --classes CLASSES --ports PORTS --capacities CAPACITIES
    python -B structural_check.py
    python -O -B structural_check.py

The main receipt contains input SHA256 values and exact signature collisions;
normal and optimized runs agree. The structural consumer independently
checks (1)--(3) over every fixed-path mask at n=3..6, the path-dependence
hostile, and both n15 witnesses including all induced-subset path counts.
No removable assertions, solver, or nonstandard package is used. Pass
`--write` only for deliberate receipt regeneration.
