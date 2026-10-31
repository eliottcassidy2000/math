# Session130: reflection, literal odd cycles, and the blue/black boundary

**Results first.** This packet proves a length-refined parity theorem for odd-cycle
packings in every self-converse tournament, and gives a smallest-order hostile
showing why even the reflection-decorated odd-cycle conflict graph cannot recover
blue degree. It also pins the inherited literal-line machinery. No publication,
theorem-number reservation, or claim of external priority is made here.

## 1. New theorem: at most one odd multiplicity survives

Let `T` be self-converse, and choose an **involutory** anti-automorphism `rho`.
A directed odd cycle is a directed cyclic vertex word modulo rotation, **not a
vertex support** and not modulo arbitrary reversal. Define

    C* = reverse(rho(C)).

This is an involution on the literal directed odd cycles of T. Put

    F_T(z3,z5,...) = sum_K product_(C in K) z_|C|,

where K runs over unordered collections of pairwise vertex-disjoint literal odd
cycles, including the empty collection. The usual OCF independence polynomial is
`F_T(x,x,...)`; OCF gives `H(T)=F_T(2,2,...)`.

Let **P_rho** be the graph whose vertices are nonfixed cycle orbits `{C,C*}`
such that `V(C)` and `V(C*)` are disjoint. Its vertices are weighted by
`z_|C|^2`; two vertices conflict iff their union supports intersect. Distinct
orbits with the same support remain distinct vertices. For a fixed cycle D,
let `P_rho[D^c]` mean the induced graph on these orbit vertices whose union
supports avoid `V(D)`. Write `I(G;weights)` for the weighted independence
polynomial. Then, coefficientwise in characteristic two,

    F_T(z) = I(P_rho; z_length^2)
             + sum_(D=D*) z_|D| I(P_rho[D^c]; z_length^2).       (1)

**PROVED, all orders.** An anti-involution has at most one fixed vertex: two
fixed vertices would force their mutual arc to reverse. It therefore has zero
fixed vertices at even order and exactly one at odd order. Every fixed odd
cycle has an invariant odd-size support, so contains that unique vertex.
Consequently two fixed odd cycles can never be disjoint. The involution `*`
pairs nonfixed cycle packings with another packing having the same length
monomial. Their contributions vanish mod 2. A fixed packing consists of complete
two-cycle orbits whose supports are disjoint, together with either zero or one
fixed cycle. These choices are exactly the right side of (1). This proves the
identity, including the conflict restrictions and multiplicities. QED.

### Coefficient consequences

For a finite vector `m=(m3,m5,...)`, let `N_m(T)` count disjoint odd-cycle
collections having precisely `m_l` cycles of each length l.

* **Every self-converse T:** if at least two entries of m are odd, then
  `N_m(T)` is even.
* **Even-order self-converse T:** if any entry of m is odd, then `N_m(T)`
  is even. Equivalently, `F_T` mod 2 is a polynomial in the squared variables.
* **Odd-order T:** all possible nonsquare sectors have the form
  `z_l B_l(z3^2,z5^2,...)`, and their coefficients are given explicitly by
  fixed length-l cycles through the unique center and the compatible paired
  cycles outside them.

For example, the number of disjoint `(3,5)` cycle pairs is even in **every**
self-converse tournament. At even order the number of disjoint `(3,3,3)`
triangle triples is also even. At odd order that latter statement fails:
`C3 -> C3 -> C3` has exactly one such triple.

This is a length-resolved extension of the inherited fixed-cycle argument, not
a new proof or discovery of OCF. Setting all variables equal yields

    I(Omega(T),x) = I(P_rho,x^2)
                    + x sum_(D=D*) I(P_rho[D^c],x^2) mod 2.    (2)

The even-order consequence `H(T)=1 mod4` was already proved in
**THM-791 / H-companion-laws-to-the-transitivity-flow** and is credited there,
not presented as new. The new content here is the full fixed-packing quotient
and its length-refined forbidden parity sectors.

### Sharp positive controls

For a chain of k directed triangles, `F_T=(1+z3)^k`. The anti-involution reverses
the component chain. Its fixed-packing polynomial is

    (1+z3^2)^(k/2)                         k even,
    (1+z3)(1+z3^2)^((k-1)/2)               k odd.

These agree mod2, and realize the optional fixed-center factor. The checker
tests k=1,...,5, including orders9,12,15 and packing degrees3,4,5. This prevents
the finite n<=7 checks from concealing a failure at the first triple channel.

## 2. New hostile: even equivariant conflict data loses blue degree

Use the owner's explorer gauge: fixed path `n -> ... -> 1`, tiles `(x,y)`
ordered by y increasing then x decreasing, bit0 follows `x -> y`.
At n5 take the reflection-fixed words **1** and **26**. Relabel vertices as
`0,1,2,3,4` along the forward base path. Their out-neighbor bit masks are

    A=(14,28,24,16,1),      B=(22,4,24,19,2).

For both, `rho(i)=4-i` is an anti-involution. Both have

    scores=(1,1,2,3,3),  H=9,  |Aut|=1,
    F_T(z3,z5)=1+3z3+z5,  Omega(T)=K4.

In each conflict graph, rho fixes one triangle and the 5-cycle, and swaps the
other two triangles. The following bijection preserves **cycle length, every
conflict edge, and the induced anti-involution**:

| A cycle | B cycle |
|---|---|
| `(0,1,2,3,4)` | `(0,4,1,2,3)` |
| `(0,2,4)` | `(1,2,3)` |
| `(0,1,4)` | `(0,2,3)` |
| `(0,3,4)` | `(1,2,4)` |

Nevertheless their blue tiling degrees are **1 and 3**, respectively. Since
both are rigid, this is also their number of palindromic Hamiltonian vertex
words. The checker directly enumerates all nine Hamilton words of each and
tests the positional-reversal anti-automorphism on every pair, independently
of the tiling census. Their fundamental-black degrees are `(9-1)/2=4` and
`(9-3)/2=3`.

Thus even

    (length-labelled odd-cycle conflict graph, induced reflection,
     score multiset, automorphism order)

does **not** determine blue degree, black degree, or black charge. These are
owner merged nodes4 and5 from session129. This is the first possible order:
the merged nodes through order4 have distinct H values, as checked by the
exhaustive small atlas; at orders1,2 there is only one tournament class.

**First missing coordinate.** In A all odd cycles share vertices `{0,4}`;
in B they share only vertex `{2}`. Pairwise nonempty intersection (K4) erases
the actual vertex-incidence hypergraph. The equivariant conflict graph retains
the parity mechanisms of section1 but does not recover Hamiltonian-path
incidence. No claim is made that just restoring these intersection sizes is
universally sufficient.

## 3. A lawful OCF-to-black-incidence bridge, with its sidecar exposed

Let K be the **literal** black multigraph on converse-merged nodes, and let Q
be its reflection-twin quotient. Keep loops and parallel edges. At a
self-converse node represented by T, set `a=|Aut(T)|`, `b=blue degree`, and
`alpha_j=#j-component odd-cycle packings`. The inherited fibre/twin laws and
OCF give exactly

    deg_Q(T) = (H(T)/a-b)/2,
    a deg_Q(T) = (1-a b)/2 + sum_(j>=1) 2^(j-1) alpha_j.       (3)

As a is odd, the black charge `chi=deg_Q mod2` is

    chi(T) = alpha_1 + (1-a b)/2 mod2.                         (4)

Here `alpha_1 mod2` can be read from the fixed literal cycles in (1).
At a non-self-converse merged node the two ordinary fibres are exchanged by
reflection, and `deg_Q=H/a` is odd, so `chi=1`.

For every set S of merged nodes, the inherited handshake/cut law therefore
becomes the explicit cycle-based congruence

    |delta_K(S)| = 2 sum_(u in S) chi(u) mod4.                 (5)

This is a proved **composition of inherited OCF and twin incidence**, not a
claim to have discovered their original cut law. The blue sidecar b is
indispensable, by section2. At n5, `alpha_1=c3+c5`; discarding 5-cycles is
already wrong, while retaining them still does not remove the blue sidecar.
Loops count twice in degree and zero in cuts. The checker tests every cut
at orders3,4,5, and every singleton cut at order6.

## 4. Inheritance: closest mechanisms and corrected near misses

Source files and exact commits/blob hashes are in `SOURCE_PINS.json`. Initial
read pin was `5431278adcd6269b5eca7756acc7a842938cd8cf`; the fresh source pin is
`5277b4505b6426f2381a01eb9186f0779566ce4e`. The latter includes the incoming
PRIMEGAME comparison below. An inherited PROVED status is not a fresh
independent audit of every proof in that file.

### Literal lines and reflection

* **THM-643 / blue-black-line-parity-structure:** `X_n=F2^m`,
  `m=binom(n-1,2)`; complement C adds the all-one word; reflection R permutes
  tile coordinates. Blue iff `Rt=t`. Literal lines are C-orbits, not full
  tournament-converse pairs. Ordinary fibre size `H/Aut` is odd. SC nodes
  have odd positive blue degree; non-SC nodes have zero blue degree.
* **THM-644 / gridsym-fiber-law-and-anti-redei**, part(a): grid symmetry iff
  path-position reversal is an anti-automorphism. Its historical conjectural
  parity label is superseded by **THM-647 / anti-redei-theorem**, which proves
  odd fixed Hamilton-path count for each involutory anti-automorphism.
* **THM-796 / three-sorted-recursive-tiling-line-node-incidence**, sections1-4:
  the faithful object is the incidence span `M_n <- X_n -> L_n`; no canonical
  complement involution exists on nodes. Endpoint recursion needs a phase bit.
  The reflection-defect exact sequence makes blue precisely zero defect and
  pairs black multiplicities even **inside each exact defect fibre**.
* **THM-4731 / tiling-line-metagraph-self-lines-by-twisted-Burnside-...**, §3:
  black R-twins have identical merged endpoints, hence K=2Q including loops.
  Its twisted Burnside computation, symmetric-core depth, and tribonacci
  partner law are additional useful sidecars, not rerun here.
* At n5 there are64 tilings,32 literal lines,8 blue and24 black, including4
  black loops. Quotienting the black twins yields12 black edges, not24;
  adding the8 blue edges gives20 quotient edges, not32.

### Colour, pure-blue atoms, and prior structural corrections

* **THM-649 / flip-parity-trilogy-affine-D-binomial-law** and
  **THM-790 / blue-parity-law-proved:** complement's triangle current is
  supported on the path legs/apex. At odd order every blue line changes
  triangle parity, so blue loops are impossible and the blue graph is
  bipartite by c3 parity. The even blue-loop guess1,2,4,8 fails: n10 has24.
* **THM-4730 / strong-pure-blue-alphabet-complete-points-of-symmetry-...:**
  the **conclusion** that strong pure-blue atoms are K1,C3,regular cyclicC5
  is current PROVED, but not by its historical odd-circulant argument.
  Nonconstant reversal-partner sums rule out affine reflection only. The
  independently audited `pure_blue_atom_classification_20261020/
  UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md` supplies the completion.
  `merged_metagraph_finite_pure_atoms_20261020/
  FINITE_ATOMS_AND_RIGID_CLOSURE.md` removes the old order62 rigid cutoff.
  Auxiliary group lemma ML is not thereby proved. Regular cyclicC5 is not
  a Paley tournament.
* Pure-blue does **not** mean leaf: regularC5 has three blue tilings, landing
  at three different nodes. Even-order pure-blue nodes do have singleton
  fibres and mixed partners; the unrestricted historical claim is false.
* **THM-4731**, §2: odd-order black self-line count `SC(n)/2` fails at n9;
  odd-order self-line endpoint rigidity fails at n11. These were conjectures,
  not dependencies of the twin theorem.
* **THM-023 / blueself-parity** is historical and should not be inherited
  indiscriminately: its sentence that blackself tilings are not self-converse
  is contradicted already by the n5 black loops at the H11 and H13 SC nodes.
  Use the corrected line conventions and THM-649's odd-order obstruction.
* Pure-blue exclusions matter in black-detour claims: `C3 -> C3` and
  `C3 -> 1 -> C3` have no black path. Active Guardrails45c points to
  `shape_ports_20261025/tournaments.md` for the repaired scope.

### Original OCF mechanism and corrections

* **THM-002 / ocf:** `H=sum alpha_k 2^k`, all orders, with vertices of Omega
  equal to **literal directed odd cycles**. Its authoritative imported proof
  is credited to the source literature in that file; we use the identity and
  make no new literature claim or claim of an independent internal OCF proof.
* **MISTAKE-023**, reiterated in the Oct26 correction of **THM-209**:
  cycle support is not a cycle. A five-vertex support can carry three distinct
  directed5-cycles. Collapsing these produces11 where H is15.
* **THM-466 / the-2-adic-digits-of-H-are-odd-cycle-collection-counts** has a
  current correction: the coefficient profile determines H, not conversely.
  `C3 -> C3` has `(1,2,1)` and a five-vertex H9 tournament plus a singleton
  has `(1,4,0)`; both evaluate to9. The full polynomial retains their difference.
* **THM-552 / c3-parity-dichotomy-self-converse** is the inherited fixed-vertex
  mechanism for triangles. **THM-791** extends it to all single odd cycles
  and even-order `H=1 mod4`. Section1 above extends it to length-refined
  disjoint collections and exhibits the admissible pair carrier explicitly.
* **MISTAKE-033/065:** tile complement is not tournament converse; converse
  read along reversed path gives reflection **without negating bits**.

### Divisor fold and PRIMEGAME: what is actually transported

* **THM-4733 / divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-...**:
  the divisor plane is an equivariant fold of the **four mirror-moved tiles**
  at n5. Summing one mirror pair gives exponent0,1,2; retaining the other two
  bits gives the p^2 q r divisor box. It preserves the Klein action and named
  gradings, not tournament isomorphism, H, line colour, or endpoint incidence.
  The two mirror-fixed tiles block an equivariant embedding into the genuine
  whole presentation cube: CR has no fixed presentation, while the divisor
  product involution has fixed points. Do not re-label this as a map onto the
  ten merged tournament nodes.
* Session129's primitive endpoint-deletion vectors classify the ten n5 merged
  nodes FINITE-EXACT. The axis rows are exactly the three pure-blue nodes at
  n5. Its encoding `7^a17^b55^c` is faithfully invertible by valuations but
  the numerical labels are chosen. The actual block products
  `55/637,539/153,585/34` are a second independently decodable dictionary.
  Neither encoding supplies the literal line pairing. At n6 a mixed sandwich
  already has an axis deletion vector.
* Incoming **THM-4854 / primegame-four-blocks-register-decomposition-...** and
  its source note are **PREAUDIT source claims for this packet**: the theorem
  contains `[[AUDIT-THM]]` and the source §7 contains only `[[AUDIT]]`, despite
  broader prose saying audited. The native identity `K_P(T)=T^op xor E(P)`
  agrees with inherited canon. To recover the32-line merged graph, quotient
  by ordinary relabelling and merge converse **vertices only**, retaining
  literal lines. A full converse edge quotient instead collapses black twins.
  Its phrase “re-marking along another Hamiltonian path is exactly a black
  move” needs the actual non-palindromic-path test: another path can be blue.
  No unqualified three-move, shortest-distance, or arithmetic classification
  is inherited from this preaudit note.

### Closest prior to an integral cycle-lattice splitting

The targeted prior search found `merged_metagraph_black_quotient_bridges_20261020/
README.md` §2: over F2, summing each twin-edge pair maps the cycle space of K
onto Q's cycle space, with one digon kernel per doubled edge, giving
`beta(K)=|E(Q)|+beta(Q)`. This is a generic doubling identity. It does **not**
contain an integral plus/minus lattice index or the full n5 graph's claimed
index `2^10`. The consulted Session75 “spectral” result is a tournament
nuclear-norm bound, not that integral graph lattice. This is a bounded
inheritance search, not a proof that no earlier record exists anywhere.

## 5. Best next structural routes

1. **Highest-value positive route:** retain original vertex supports with the
   anti-involution, and compare the fixed-packing quotient with the palindromic
   path complex. Section2 gives a cheap hostile for every proposed carrier:
   it must distinguish the two H9 tournaments, not merely H11 from H13.
2. **Coloured incidence:** carry the blue count as an honest sidecar in (3)-(5).
   The parallel reflection pairing and fixed-cycle pairing are both real
   involutions, but on different sets; a bijection between their orbit sets is
   not supplied by matching parity.
3. **Integral refinement:** compare the original literal graph with the
   black-twin quotient without discarding loops, edge orientations, or the
   index between integral invariant/anti-invariant sublattices. F2 doubling
   alone cannot pay that index.
4. **Divisor/arithmetic lane:** keep the marked endpoint-deletion span and
   exact quotient losses. Multiplicatively independent codes label the node
   table but cannot reconstruct complement adjacency without a line sidecar.

## 6. Verification and stopping boundary

`check_equivariant_ocf.py` is standard-library-only and imports no inherited
producer or data file. It checks:

* full fixed-path metagraphs n3..6, ordinary and converse-merged classes;
* direct automorphism count, Held-Karp paths, and independently enumerated
  literal cycle packings;
* scalar packing polynomial by smallest-vertex DP versus a separate explicit
  collection recursion for the length-refined polynomial;
* every reflection-fixed word through n7 (2,4,16,64,512), and128 deterministic
  reflection-fixed words at n8 (every32nd half-cube address, **not a census**);
* nonzero `(3,5)` packing coefficients at n8, the five triangle-stack controls,
  the complete multivariate mod2 quotient, and all forbidden parity sectors;
* every n5 cut congruence and the loop-aware black incidence formula;
* the explicit equivariant conflict-graph hostile and direct blue path counts.

Reproduction:

    python3 check_equivariant_ocf.py > normal.json
    python3 -O check_equivariant_ocf.py > optimized.json
    cmp normal.json optimized.json

Failure gates are explicit exceptions, active under optimization. The scope of
the theorem is all orders because of its proof; the finite checks are audits,
not an extrapolation. Higher binary digits of weighted evaluation are not
recovered merely by reducing coefficient polynomials mod2: carries matter.
The companion Hamiltonian lane's mod8 refinement is a separate statement.
