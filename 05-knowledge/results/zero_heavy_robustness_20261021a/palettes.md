# Unbounded conditional robustness with every physical port retained

**PROVED:** exact boundary repair distance and the optimal robustness of a
root-only gate composition. **FINITE-EXACT:** source gate counts and the
explicit placements/replays. This gives no new chromatic lower bound.

## 1. Inherited source and changed question

The source is the actual 16-point gate in
`kernels_games_detours_20261020f/forcing.md`, not an abstract palette name.
Its root is vertex 0 and its nine boundary ports are 7 through 15. The
three triples 7--9, 10--12, 13--15 receive distinct colours a,b,c. Write
A={a,b,c}. With six available colours, the root has precisely A as its
extension support, with 24 extensions for each member. The full induced
unit graph has 31 edges and is four-colourable.

The two-copy construction and finite-avoidance lemma in
`first_exit_complement_games_20261020/PALETTE-GLUE.md` retain physical ports
and permit an induced union meeting at the root alone. The inherited
`charge_palettes_boundary_20261020g/palettes.md` finds that every proper
one-port perturbation of the complementary two-copy bad input repairs it.
One particular perturbation is enough here: change port 7 from a to b.
For (a,b,c)=(0,1,2), its root-count vector is

    (60,84,84,8,8,8).                                    (1)

It is a proper boundary and therefore allows every root colour. For
example, after that port change a literal colouring for root 3, in vertex
order 0 through 15, is

    (3,0,4,5,4,5,3,1,0,0,1,1,1,2,2,2).

Permuting the three unused colours supplies every root colour outside A. The new
checker reconstructs the full source edge set from its exact coordinates
and recomputes (1), retaining one literal 16-vertex witness for every root
colour. No geometric port is identified merely because it has the same
prescribed colour.

The changed target is **conditional boundary robustness**, not forcing
these input equalities in every six-colouring. Repeating a fragile input
can amplify the former without making any progress on the latter.

## 2. Exact root-only repair metric

Take N>=1 copies of the gate, identifying only their actual roots, with
no other shared vertices or cross-copy edges. Retain all 9N distinct
boundary ports. In copy i prescribe the triple palette A_i, a three-subset
of the six colours, with a fixed ordering of its three elements. Let b be
this full proper boundary colouring. Define

    m_c = #{i : c belongs to A_i},
    r(b) = minimum Hamming distance from b to a proper boundary
           colouring that extends to a proper six-colouring.

Hamming distance counts changed **physical ports**, not changed palettes,
colours, copies, or quotient blocks. Then

    r(b)=N-max_c m_c.                                    (2)

More precisely, the minimum number of changed ports among extensions
whose root is c is exactly N-m_c.

**Lower bound.** Any copy with c outside A_i cannot retain all its original
nine port colours if its root is c. Each such copy requires at least one
changed port, and the port sets of different copies are disjoint. There
are N-m_c such copies. Every extension has some root colour c, so minimize
this lower bound over the six possible colours.

**Upper bound.** Choose c. In every copy whose palette excludes c, apply
the proper one-port perturbation (1), with its palette labels permuted.
Leave all other copies unchanged. Every changed copy now permits c, and
every unchanged copy already did. Their extensions combine because the
only common vertex is the root and there are no cross-copy edges. The
result differs at exactly N-m_c physical ports, proving equality.

This explicit nearest repair has exactly

    24^(m_c) * 8^(N-m_c)                                 (3)

extensions with root c. This counts internal choices, not different
boundary repairs. A repair may also admit other root colours.
The original boundary itself has 24^N extensions for each c in every A_i
and zero for all other c; in particular it is blocked iff the intersection
of the A_i is empty.

## 3. Sharp amplification, and its severe limitation

Since each palette has size three, sum_c m_c=3N. Thus max_c m_c>=ceil(N/2)
and

    r(b)<=floor(N/2).                                   (4)

This is optimal within this root-only architecture. Take floor(N/2) copies
with palette {0,1,2} and ceil(N/2) with palette {3,4,5}. Their maximum load
is ceil(N/2), attaining (4). For N=2s, the blocked boundary is unaffected by
every perturbation of fewer than s physical ports, yet an explicit s-port
repair exists. There are 30s+1 vertices, 18s ports, and 62s unit edges.
For N=1 the radius is zero, correctly: the original boundary extends.

All these graphs are still four-colourable. Fix a four-colouring in each
copy, permute its colour names so the common root agrees, and combine.
Hence arbitrarily large conditional repair radius does not imply a
six-chromatic subgraph, an equality synchronizer, or even any increase in
chromatic number. It is a sharp hostile to using an absolute robustness
threshold as a sufficient geometric-forcing certificate.

The normalized radius never exceeds 1/18 of the 9N physical input ports.
Repeating root-only copies cannot improve that limit. The full palette
load vector (m_0,...,m_5), rather than only the number of copies or the
cardinality of one pairwise intersection, is the exact robustness sidecar.

## 4. Actual geometric realization at every finite N

Keep the first gate at its original coordinates. Given the existing
finite union F and the next source gate H, rotate H about the common root
by

    rho_t=(1-3t^2+2t sqrt(-3))/(1+3t^2).

The inherited finite-avoidance proof supplies a positive integer

    t<=4(|F|-1)(|H|-1)+1=900j+1

when j copies are already present. Each nonroot cross-pair excludes at
most two t values for coincidence and two for a new unit edge. The
polynomials are nonzero because both points are nonzero. The first
surviving t therefore makes the induced unit graph exactly the desired
root-only union. This retains every existing point and its named ports;
no drawing or approximate distance replaces the exact unit-incidence test.
Induction realizes every finite N in Q(sqrt(-3),sqrt(-11)).

The checker supplies an explicit six-copy placement (91 vertices), tests
all 4,095 unordered distances exactly, verifies the 186 expected edges,
and verifies the combined four-colouring and nearest-repair six-colourings
on that actual induced graph. It also enumerates palette multisets and
independently compares the incidence formula with all subsets of repaired
copies in the declared finite universe. These computations are controls
of the all-N proof, not a search for global six-colour forcing.

## 5. Source-to-target obligation for genuine forcing

A context that forces even one equality between distinct ports in every
proper six-colouring must itself be at least six-chromatic: from a
five-colouring, recolour one selected port with a fresh sixth colour.
This inherited fresh-colour obstruction is unaffected by repetition.

For a five-colour circle-saturation attempt, the actual source object must
instead be a finite unit graph with a specified new centre and its exact
neighbour set. The target predicate is: every proper five-colouring of the
source uses all five colours on those neighbours. Adding the centre then
requires a sixth colour. A rooted gate word, a preassigned boundary, or
large repair radius supplies none of that universal quantifier. Every
boundary partition, genuine unit contact, and relative palette permutation
must survive any proposed map to that target. This is an OPEN obligation,
not a claim that the present carrier supplies it.
