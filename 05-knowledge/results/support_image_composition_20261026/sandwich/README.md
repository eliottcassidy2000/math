# Four markers recover every joint two-port tournament

2026-10-26. **PROVED: S1–S4, independently audited.** The proofs below are
self-contained; the finite controls are **FINITE-EXACT**. See the
[independent audit](../audit_sandwich/REPORT.md). The two-marker strong-core theorem is **INHERITED PROVED**.
This note proposes closure of the specific `A2*h*A2` all-core question in
[HYP-9416, Open](../../../hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md#open),
not its strong-atom response-classification or dynamic-compression questions.

## Inheritance and the changed hypothesis

The closest proved mechanism is [unnamed-carrier U1](../../unnamed_carrier_dissection_20261026/README.md#u1-one-singleton-on-each-side-recovers-every-strong-joint-port-state):
`e*h*e` recovers a **strong** core with distinct ordered ports. Its two added
vertices have outdegree one and indegree one. The only extra candidates are
the entrance and exit; deleting an incorrect candidate leaves a vertex with
no outgoing or incoming core arc. Requiring the remaining core to be strong
therefore recognizes the markers. The earlier degree-neighbour proof and
this deletion proof are inherited alternatives, not new results here.

The canonical hostile is loss of **joint** port incidence: separate entrance
and exit types do not determine a two-sided response (the inherited
`289/253` example). The corrected near miss is MISTAKE-676: a complete
unrooted attachment object was confused with a scalar path-count observer,
and a named protocol with fresh counted re-marking. The least-used sidecar
is the automorphism kernel of the added markers. Here that kernel is
nontrivial precisely at the sink/source boundary.

The live board has four items: strong-core deletion; extremal-degree chains;
joint-port composition; and a separate image-membership guard. The chosen
research move is the repository card “Controlled forgetting and unlabeled
quotients require a sidecar”: inspect the entire fibre, not merely one
successful marker choice. Source versions and exact hashes are in
[SOURCES.json](SOURCES.json).

## Convention and S1: the universal four-marker decoder

A state is `h=[C;s,t]`, where `C` is any finite tournament, `|C|=n>=2`, and
`s!=t` are ordered ports. No strongness or Hamilton-path hypothesis is made.
For two states, the inherited hinge product `h*k` keeps both interiors,
points every cross arc from the second factor to the first except the arc
`exit(h)->entrance(k)`, and retains the outer ports. It is associative
([rooted-lines R1](../../rooted_lines_frobenius_sheets_20261026/README.md#2-r1-the-joint-two-port-response-and-its-convolution)).
The singleton `e` has coinciding formal ports. Write `A2=e*e` for the
forward rooted arc and define the unrooted tournament

    E(h) = A2*h*A2 = e*e*h*e*e.

Name the left markers `a0,a1` and the right markers `b0,b1`. Their forward
chain edges are `a0->a1->s` and `t->b0->b1`. All cross arcs not among these
four consecutive-factor hinges point from later factors to earlier factors. In particular,
the right marker pair dominates the left marker pair.

**S1.** For all `n>=2`, `E(h)` and `E(k)` are isomorphic if and only if
`h` and `k` are jointly rooted-isomorphic. An unrooted input in the image
can be decoded in `O(N^2)` elementary adjacency operations, `N=n+4`, up to
joint rooted isomorphism. This is a decoder, not a general isomorphism test.

**Proof.** Both left markers have outdegree one. For each core vertex,

    d_E^+(v)=d_C^+(v)+2-[v=s]+[v=t].                    (1)

Both right markers have outdegree `n+2`. Thus the entire degree-one
set `L` is exactly `{a0,a1}`, except that it also contains `s` when `s` is
a sink of `C`.

If `s` is not a sink, the induced edge on `L` is `a0->a1`. It orders both
markers intrinsically; `s` is the unique out-neighbour of `a1`. If `s` is
a sink, `L={a0,a1,s}` is a directed triangle

    a0 -> a1 -> s -> a0.

Every vertex outside this triangle dominates all three of its vertices.
It is therefore a homogeneous set, and all three cyclic rotations extend
to automorphisms fixing its exterior. Any directed two-edge path within
it gives an equally valid interpretation of `(a0,a1,s)`.

Dually, with

    d_E^-(v)=d_C^-(v)+2+[v=s]-[v=t],                    (2)

the entire indegree-one set `R` is `{b0,b1}`, or `{t,b0,b1}` precisely when
`t` is a source of `C`. The former intrinsically orders `b0->b1` and reads
`t` as the unique in-neighbour of `b0`; the latter is the homogeneous
triangle `t->b0->b1->t`, which dominates every exterior vertex.

Since `N>=6`, the sets `L` and `R` are disjoint. Delete the two recognized
markers on each side, retaining the remaining triangle vertex as the port
when a triangle occurs. All choices give isomorphic rooted cores: each
change of choice is a cyclic rotation supported on a homogeneous triangle.
The two rotations have disjoint supports and can be made independently.
Thus every isomorphism of encoded tournaments induces a rooted isomorphism
of decoded cores. Conversely, a rooted core isomorphism extends by mapping
each named added marker to its counterpart. This proves both directions.
The construction uses degrees, adjacency and at most nine marker choices.
The two possible simultaneous triangles include `n=2`, `t->s`; nothing in
the proof requires a third core vertex. ∎

There is also a choice-free quotient description. Delete each two-vertex
degree-one chain; replace each three-vertex degree-one triangle by a single
marked vertex. The triangles are homogeneous, so the quotient arcs are
well-defined. This yields the same rooted core intrinsically.

## S2: exact marker symmetry and compositional transport

Put `epsilon=[s is a sink]+[t is a source]`. Then

    Aut(E(h)) = C3^epsilon x Aut(C;s,t).                (3)

Here the isomorphism uses the displayed encoding of `C`; intrinsically the
canonical map to the automorphism group of the decoded quotient has kernel
`C3^epsilon`. Indeed an encoded automorphism acts on the decoded core. A
port-fixing core automorphism extends while fixing all four markers. After
removing this extension, only rotations of the exceptional triangles remain.
The extension fixes each exceptional triangle pointwise, so it centralizes
those rotations: the product in (3) is **direct**, not merely semidirect.

Consequently the encoding is faithful on isomorphisms and injective on
isomorphism classes. It is fully faithful on the subcategory where neither
exception occurs; on all cores it is not full. For the reverse rooted arc
`t->s`, the source automorphism group is trivial but the encoded group has
order nine. Quotienting the image groupoid's isomorphisms by these intrinsic
rotation kernels recovers the groupoid of joint two-port states.

On encoded **isomorphism classes**, define the following actual operation:
decode both inputs, apply the rooted hinge product, then encode the result.
This is well-defined independently of every triangle choice by S1. It is
associative by the inherited rooted associativity, and

    E(h) odot E(k) = E(h*k).                            (4)

Thus the full joint-port composition law can be stored using unrooted
tournaments with four added vertices. Any isomorphism-invariant predicate
of the core or its ports can likewise be evaluated after decoding. No
strongness hypothesis is needed for intermediate results.

The same statement applies to a counted response histogram: replace each
joint key `h` by the unique ordinary key `[E(h)]`, retaining its coefficient.
The map on keys is injective, so it loses no coefficients or joint incidence.
In particular the inherited rooted response convolution transports through
(4). Forgetting keys in any other fashion still requires **summing**, not
overwriting, coefficients, as repaired in unnamed-carrier U1.

These are conjugated operations. They are not the native Hamilton-path or
hinge operations on the encoded tournament with its ports forgotten. For
the reverse rooted arc there is no entrance-to-exit Hamilton path, while
its encoding is an ordered join of two directed triangles and has nine
Hamilton paths. Fresh V1 or V2 dynamics can be transported by
decode–operate–encode; applying them directly to all paths of `E(h)` is a
different operation. No scalar `H(E(h))` completeness claim is made.

## S3: an exact recognition guard for arbitrary inputs

The promised image test is explicit:

1. Require `N>=6`; form the degree-one set `L` and indegree-one set `R`.
   Require both cardinalities to be two or three.
2. Try `a0 in L`, put `a1=N^+(a0)` and require `a1 in L`; put
   `s=N^+(a1)`. Dually try `b1 in R`, put `b0=N^-(b1)` and require
   `b0 in R`; put `t=N^-(b0)`.
3. Require four distinct markers, both ports outside the markers, and
   `s!=t`. Delete the markers and retain the ordered ports.
4. Re-encode this core and require exact equality under the proposed
   vertex correspondence. Accept iff at least one choice passes.

This recognizes the image exactly: S1 supplies choices for every encoded
input; a passing exact re-encoding is a witness in the other direction.
All accepted choices decode the same joint type by S1. There are at most
nine pairs of choices, so the `O(N^2)` bound includes the guard.

Degree cardinalities alone do not suffice. The six-vertex tournament with
out-neighbour bit rows

    (32,1,3,7,15,30)

has `L={0,1}`, `R={4,5}`. Its proposed left chain is `1->0->5` and its
proposed right chain is `0->5->4`: each proposed port is a marker on the
other side. The guard rejects it. This is the first such input in the
checker's declared labelled-mask enumeration, not a universal minimality
claim beyond the image's minimum possible order six.

This distinguishes the three Session71 obligations concretely. Joint
incidence is recovered by the ordered marker chains; transport is the
explicit product (4); realization is the separate image guard. The
[two-involution and seam-transport results](../../session71_dissection_norm_alignment_20261026/session71-root/BRIDGES.md#6-the-promised-structural-isomorphism-two-whole-side-involutions)
use the same need to retain joint alignment and then check admissibility.
They do not identify geometric dissections with these tournament images.

## S4: four added vertices are optimal in this fixed-context grammar

This is an exact limitation of the observer `h -> free(k*h*k')`, where
`k,k'` are fixed nonempty two-port hinge contexts, with distinct ports
except in a singleton. It is not a lower bound for arbitrary encodings.
If the total number of context vertices is less than four, the only
possibilities are `e,e`, `A2,e`, `A2op,e`, `e,A2`, and `e,A2op`.

All five fail already on three-vertex cores. With ports `s=0,t=2`, let

    C=(2,4,1),       X=(4,5,0),       Y=(6,0,2),

where rows encode out-neighbours. `C` is a directed triangle; `X` and `Y`
are transitive. They are distinct joint types. A complete exact
isomorphism check gives these collision fibres:

| Fixed contexts | Distinct cores in one output class |
|---|---|
| `e,e` | `C,X,Y` |
| `A2,e` | `C,X` |
| `A2op,e` | `C,X,Y` |
| `e,A2` | `C,Y` |
| `e,A2op` | `C,X,Y` |

The full two-vertex source universe has no such collision, so order three
is the first source order in this declared domain. These finite certificates
prove failure of every smaller context pair. Together with S1, they show
that four added vertices are the exact optimum for a single universally
faithful fixed two-sided hinge-context observer. This does not contradict
the inherited two-marker optimum for the restricted strong-core domain.

## Verification, scope and next question

From this directory run:

```bash
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
cmp normal.json optimized.json
```

The checker uses no external libraries or inherited producer imports.
It enumerates all rooted classes at orders two through five (respectively
`2,8,32,176`) from every labelled orientation. Its unrooted isomorphism
oracle exhausts all degree-compatible bijections; it assumes no marker
decomposition. It verifies every decoder choice and the full automorphism
formula. All `32,768` labelled order-six inputs undergo the image guard:
exactly `800=6!+6!/9` pass, accounting for both two-vertex rooted cores.
It also checks all small context collisions, choice-independent composition,
associativity and 150 relabelled larger cases through core order twenty,
forcing every sink/source boundary combination. Explicit `require` calls
remain active under `-O`.

This settles a representational bottleneck: arbitrary two-port state keys
can be encoded faithfully as ordinary tournaments in this attachment
grammar. Classification of strong atoms by response support or histograms
is a different question and remains open beyond the inherited finite
census. No unrestricted A000568 recurrence, all-order black-core
connectivity, planar Jacobian, Collatz, RH, or P-versus-NP conclusion follows.
The next structural question is an economical faithful encoding for ordered
multi-run port lists under their actual composition law; the present proof
establishes only the two-port case, not that proposed generalization.
