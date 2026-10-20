# Independent audit: Euler, Hamiltonian, and marked torus transport

20 October 2026. **PASS, with a publication-only attribution/wording overlay.**
The producer directory remains frozen. This audit made no repository or remote
changes. The inherited shared repository was inspected at
`c402d288ac84f695e0fb355444ee67bfb7917d4f`.

## Outcome and precise qualifications

No substantive error was found in the proofs, marked slopes, counts, complete
finite certificate, or source-hypothesis application. Two prose clarifications
are recorded in `publication_overlay.diff`:

1. The concurrent `kernels_games_detours_20261020f/walks_knots.md` already proves
   directed Euler words correspond to Hamiltonian words in the directed line
   graph. Credit that inheritance; do not advertise the base bijection as new.
   The present cyclic normalization, checkerboard obstruction, marked K7
   analysis, and flux certificate remain distinct contributions of this packet.
2. The first two models disprove **nondecreasing** monotonicity of Euler count
   in raw H. They alone do not disprove nonincreasing monotonicity. All three
   models together disprove either direction. The overlay makes this explicit.

The overlay also resolves the geometric-companion link to the intended joint
publication and gives the change between its basis and this packet's basis.
Neither change alters the mathematics or any numerical artifact.

## Frozen input and replay integrity

Every one of the ten entries in the producer's `SHA256SUMS` passed before
execution. Copies alone were run under `replay/`. The normal and optimized
producer verification runs, and both certificate replays, reproduced all
recorded output and generated JSON bytes exactly. The same ten checksums
passed after those runs. The original producer files were not written.

Reproduction from this published packet directory:

```
python3 verify_euler_torus.py
python3 -O verify_euler_torus.py
python3 replay_certificate.py
python3 -O replay_certificate.py
python3 independent_audit.py .
python3 -O independent_audit.py .
```

The certificate's independent generating-function count, 2,640 with maximum
layer 713, is sufficient for its stated completeness proof. The additional
Gray-code census below was an adversarial audit, not a necessary hypothesis
or a proposed strategy for larger orders.

## Proof and convention audit

### Euler-to-Hamilton and BEST

`H` counts ordered Hamiltonian vertex words, and `E` counts Euler arc words
modulo cyclic rotation, not modulo reversal or graph automorphisms. Because
arcs are distinguishable, a cyclic Euler circuit has exactly m cuts. Balance
forces any all-edge trail to close, including the one-arc loop boundary.
Therefore H(L(D))=mE(D). Strong connectivity provides existence and excludes
isolated irrelevant vertices in the BEST application.

The cofactor of outdegree-minus-adjacency counts in-arborescences. The formula
E=tau_v product_v(d_v-1)! is correct for cyclic circuits; no extra root-degree
factor belongs here. A vertex-started count without a specified first arc
would instead introduce d_v. The source defines the same cyclic equivalence
before Theorem 1. [Creed–Cryan, pp. 1–2](https://arxiv.org/pdf/1202.2156).

The regular n>=7 parity separation, and the exceptional (H,E) values (3,1),
(15,11), are sound. The order-five uniqueness argument determines every arc.
The line graph is not a tournament: disjoint original arcs are nonadjacent.
Independent boundary controls include three labelled loops (H,E,m)=(6,2,3),
two parallel arcs each way (8,2,4), and two cycles sharing a vertex (4,1,4).
A directed two-edge open path gives (1,0,2), demonstrating why balance cannot
be removed from the cyclic-count identity.

### Homology, checkerboard, and local operations

The edge-lift label has the stated gauge change k_u-k_v; this cancels for
closed walks and for the full chain of a balanced orientation. An Euler
circuit uses every edge once, so its winding polynomial is one monomial with
coefficient E. This is not the cycle polynomial of the original tournament.
An open path retains endpoint-relative displacement, not a closed homology
class.

In the checkerboard construction, the full black-boundary chain is a boundary.
The reversed-edge chain R satisfies 2 partial R=0 and hence partial R=0 over
integral vertex chains. Thus every Euler flux is even. On a torus no such flux,
including zero, represents a simple essential curve; zero does not rule out
contractible curves. On higher-genus surfaces the homology divisibility proof
survives, but zero need not mean contractible. No smoothing theorem in S3 is
being asserted.

Reversing a directed cycle changes the full chain by -2C. The two-complement
identity J(final)=J(initial)-2P+2Q is correct in the integral edge-chain group,
with Q oriented in the intermediate tournament. Matching ordered endpoints
makes Q-P closed and restores the score vector. The independent audit checked
this edgewise and after homology for all same-endpoint presentations from the
three models: carousel 1,162, middle 1,239, Paley 1,050. It does not identify
these returns with a single original metagraph move.

### Embedding, markings, knot interpretation

The six unit residues are the six nonzero elements modulo seven, giving a
simple K7 quotient. The fourteen specified triangles have zero displacement.
Independently reconstructed oriented boundaries put each edge in two distinct
faces with opposite signs; every vertex link is one six-cycle. Thus the
connected closed orientable cellular surface has Euler characteristic zero.
The basis U=(7,0), V=(-3,1) is positively oriented and has determinant seven.
The geometric companion instead uses g1=(-3,1), g2=(1,2), so

```
g1=V, g2=U+2V;   (alpha,beta)_geometry -> (beta,alpha+2beta)_here.
```

The cocycle ((i+s-j+3t)/7,t), fluxes (6,14), (4,8), and (0,0), and every
listed cycle slope are correct. For the fixed unknotted S3 torus marking,
carousel has 16 nontrivial Hamilton-cycle knots and Paley has two. This is
conditional on the chosen ambient torus marking, not graph isomorphism.

The multiplier 2 is represented by A=[[-1,-1],[1,0]] on lattice coordinates
and M=[[2,-1],[7,-3]] on (U,V). The audit checked its action on every Paley
Hamilton cycle. In particular the slopes (1,0),(2,7),(-3,-7) form one orbit
but give an unknot and two distinct nontrivial torus knots. Hence that torus
map cannot extend to an ambient S3 homeomorphism preserving the embedding.
The scalar C3 hostile uses curved edges; its unbounded knot genus follows
from the usual positive torus-knot formula, stated in the cited
[Jabuka–Van Cott abstract](https://arxiv.org/abs/1907.12970).

H=nC+R partitions ordered Hamilton paths into closing and nonclosing ones.
The count R=21 in Paley is not a tournament with raw H=21. Its order-21 affine
automorphism group acts freely on ordered paths, so the 21 nonclosing paths
form one orbit. No automorphism quotient is substituted for raw H.

### Cited facial-flip theorem

Zhang–Qian assume strongly connected prescribed outdegrees and an orientable,
closed, cellular surface with each edge in two distinct faces. Regular
outdegree three on this K7 is strongly connected, and the embedding checks
above verify the remaining hypotheses. Their no-forbidden-face result gives
reachability iff orientation-homologous and directed distance by the normalized
integral face potential. Full-chain flux equality is equivalent here because
H1(T2;Z) has no two-torsion. [Sections 2–3, Theorem 3.4 and Corollary 3.5](https://arxiv.org/pdf/2011.07481).

The producer forest checks the **undirected** graph of reversible facial
triangles; it does not itself encode counterclockwise directions. This is
not a defect: the cited theorem supplies directed reachability, and the
independent audit additionally rebuilt the ccw-only graphs and verified all
894,540 ordered same-flux distances against the face-potential formula.
Maximum directed distance in this finite universe is 27.

## Independent finite-exact results

`independent_audit.py` imports no producer code. It visits all 2^21 orientation
masks by Gray code, rather than using the producer's degree-pruned generator
or its polynomial completeness count. Exactly the same 2,640 regular states
are obtained. Every H is independently checked against all 7! word incidence
masks. Relabelling the three named models independently gives disjoint orbits
of sizes 720, 1,680, 240 whose union is the entire regular universe.

The 25 flux fibres have size multiplicities 660:1, 273:6, 42:6, 14:6, 1:6.
Every fibre's H populations, undirected edge count, and degree histogram agree
with the producer JSON. All model slope polynomials agree. All Euler fluxes
are even, and every orientation is strongly connected.

The 273-state torus fibre at flux (0,2) has 770 undirected edges and the stated
degree histogram. Independently generated matchings of P6 disjoint union P7
have 273 states, 914 edges, maximum degree 11. These natural graphs are not
isomorphic. Cardinality alone supplies no Fibonacci mechanism for flux fibres.

## Scope and publication boundary

The line-graph counting identity, edge-cocycle sum, checkerboard obstruction,
and cycle/paired-path chain laws are general elementary proofs under their
stated hypotheses. BEST and facial-flip reachability are imported theorems,
not new discoveries. The numerical tables, finite class coverage, slope
polynomials, 25 components, directed diameters, and 273 comparison concern
only this specified labelled embedded K7 and declared finite controls.
Nothing proves all-order sector counts, a Fibonacci-sector recurrence,
ordinary tournament H-spectrum completeness, or connectivity of the original
black/merged path-complement metagraph.

The publication overlay changes only attribution, the monotonicity wording,
the basis cross-reference, and the joint-publication relative link. The new
SHA256SUMS adjusts only the proof document's digest. `audit_manifest.json`
records exact input, overlay, and audit-artifact hashes.
