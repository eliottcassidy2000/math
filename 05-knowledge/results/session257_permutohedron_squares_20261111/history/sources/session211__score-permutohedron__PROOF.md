# Score permutohedra, gate-erasing cycles, and the literal 7/21 boundary

Status: PROVED elementary statements below, with independently replayable
FINITE-EXACT controls. No literature-priority claim. The classical score
zonotope and strong-component product facts are recovered and derived here,
not claimed as newly discovered. The Chan–Pak spectrum problem is OPEN.

## 0. Inheritance and connection contract

The March work recognized three distinct uses of the permutohedron:
transitive tournaments inside the tournament cube (THM-207), the score
zonotope, and the ordered-partition face structure. June's cut-side reflection
suggested comparing resolutions but did not supply a sufficient counting
coordinate. Session210 supplied a different, faithful carrier: a tournament
orients every parallel edge class of the permutohedron, and its Hamilton
paths are exactly the sinks. These objects must not be conflated.

The present bridge uses the **score projection** and the canonical closure
`Gamma(T): a -> b -> T -> a`. Its anchor is the precise geometric effect of
adding the two gate vertices. Its niche is the macro strong-component face.
Its wildcard is the smallest score-preserving local switch that erases a
gate. The concept board is: scores / tight cuts / strong components /
parallel-coherent orientations / canonical gates / directed-triangle flips.

- Source: the complete tournament, or its coherent orientation.
- Target: its labeled score point in the graphical zonotope.
- Preserved: every labeled outdegree and the ordered strong-component
  partition, recovered from the minimal face of that point.
- Lost: internal arcs, Hamilton counts, and the existence of the canonical
  gate. A directed triangle is an explicit kernel move.
- Sidecar: the full parallel-class orientation restores the tournament.
- Decisive test: two strong order-five tournaments have the exact same labeled
  score, but one has a gate and three Hamilton cycles, the other no gate and
  one Hamilton cycle. Their Hamilton-path counts also differ, 15 versus 11.

The closest proved mechanism is Session210's gate bijection. The historical
near miss is treating a score fiber as a face or treating the transitive
permutohedron's isomorphism quotient as the full tournament metagraph. The
least-used coordinate is the complete collection of **labeled cut slacks**.
SOURCES.json pins the inherited sources and historical claims separately.
The score-preserving triangle move is classical (the Ryser source recorded
in THM-1945, section1 only). Non-determination of H by scores and its first
failure at order five are already explicit in corrected THM-466(iv)-(v).
The present contribution is their exact interaction with the canonical gate,
the face-to-interior map, and literal HC=7/21 pairs inside a gated score fiber;
no novelty claim is made for those inherited mechanisms.

## 1. Scores recover the macro face exactly

Let `Z_n = sum_{i<j} [e_i,e_j]`, the score permutohedron with vertices the
permutations of `(0,1,...,n-1)`. For a tournament T on [n], choosing e_i on
an edge won by i gives its score vector `s(T)=(d_T^+(i))` in Z_n. For a
nonempty proper subset S define its Landau slack by

    L_s(S) = sum_{i in S} s_i - binom(|S|,2).

The standard facet description is `L_s(S)>=0` for every such S, with total
sum `binom(n,2)`. This description follows, for example, by sorting an
arbitrary linear functional and summing the prefix inequalities; the
corresponding permutation of `(0,...,n-1)` attains the resulting bound.
Thus no integrality assertion or tournament realization theorem is needed
for the real-polytope statements below.

**Proposition 1 (macro recovery).** Write the strong components in dominance
order as `C1 -> C2 -> ... -> Ck`. Then

    L_s(T)(S) = # {arcs of T directed from S to its complement}.       (1)

The nonempty proper tight subsets are exactly the terminal unions
`Cj union ... union Ck`, for `j=2,...,k`. Therefore the labeled score vector
recovers the ordered labeled strong-component partition. Its smallest face is

    product_{j=1}^k ( Z_(Cj) + m_j * 1_(Cj) ),
    m_j = |C_(j+1)| + ... + |Ck|,                                  (2)

in the disjoint coordinate blocks, and that face has dimension `n-k`.
In particular T is strong iff s(T) is in the relative interior of Z_n.

**Proof.** Every internal unordered pair of S contributes exactly one to its
score sum; every outgoing cross arc contributes one, proving (1). A set with
no outgoing cross arc cannot split a strong component. Since the condensation
of a tournament is a transitive total order, it must be a terminal union.
Conversely every terminal union has no outgoing cross arc. Their nested
boundaries recover the blocks. The vertices at those equalities are exactly
the total orders listing C1 first, then C2, and so on. Their score coordinates
are all internal permutation scores plus m_j, proving (2). There are exactly
k independent block-sum constraints, hence dimension n-k. QED.

A singleton is strong under the standard length-zero reachability convention;
the assertion includes its zero-dimensional permutohedron.

**Corollary 1.1 (micro factors on the macro face).** Restrict Session210's
coherent orientation to this face. Its sinks are precisely concatenations
of one Hamilton path in each component. Consequently

    sinks(minimal score face) = H(T) = product_j H(T[Cj]).           (3)

The same minimal face may contain many score points; and the same score
point may represent several different internal coherent orientations with
different H. The face is a product carrier, not a recovery of its micro data.

This is the exact point of contact between the March score zonotope and
Session210's coherent orientations. On score vertices one reads the ranking
in decreasing coordinate order; Session210's position coordinates differ by
a central reflection and a translation. That change of embedding does not
change the face graph or the sink assertion.

## 2. Adding two gate vertices moves a whole face into the interior

For a nonempty core, n>=1, put a,b after the n core coordinates and define
the affine maps

    g_t(y) = ( y_1+1, ..., y_n+1, t, n+1-t ),    0 <= t <= 1.

At t=0 this is the score map of source/sink padding `b -> T -> a`,
with `b -> a`. Its image of Z_n is the face with singleton outer blocks
`b | [n] | a`. Reversing just `b -> a` changes the scores by `e_a-e_b`,
so g_1 is exactly the score map of Gamma(T). This is the score-side version
of the separate Session211 single-parallel-class flip construction.

**Theorem 2 (face-to-interior transport).** Let n>=1. For every y in the entire closed
polytope Z_n, `g_t(y)` belongs to the relative interior of Z_(n+2) for
`0<t<=1`. Its minimum unnormalized Landau facet slack is exactly t.
At t=1 this is the integer margin 1. It is not a Euclidean-distance claim.

**Proof.** For `K subset [n]`, `k=|K|`, and `L=L_y(K)`, allowing empty/full
K with L=0, the four possible kinds of subset have slacks

    K:             L + k,
    K union {a}:   L + t,
    K union {b}:   L + n+1-t,
    K union {a,b}: L + n-k.                                      (4)

The total coordinate sum is `binom(n+2,2)`. For a nonempty proper subset,
all four expressions are positive when `0<t<=1`, and are at least t.
The singleton {a} has slack t. QED.

Thus the gate erases all proper macro faces at once: every core, including
a completely transitive one, becomes strong. The grading `n -> n+2` has a
literal affine/zonotope operation behind it. It does **not** assert the
path-count recurrence `H -> H+2`; the inherited exact recurrence is
`H(Gamma(T))=H(T)+P_T(2)`.

There is a useful distinction between this embedded copy and its affine
slice. The source/sink face at t=0 becomes a strict subpolytope of the slice
at t=1.

**Proposition 2.1 (the larger degree slice).** For n>=2, the full slice

    { x in Z_(n+2) : x_a=1, x_b=n }

is `g_1(Q_n)`, where Q_n is the permutohedron with vertex vector

    (-1, 1, 2, ..., n-2, n).

Equivalently Q_n has total sum `binom(n,2)` and all nonempty proper
Landau slacks at least -1. The original Z_n is strictly smaller than Q_n.

**Proof.** In (4) with t=1, the four inequalities at a fixed proper nonempty
K ask respectively `L>=-k`, `L>=-1`, `L>=-n`, and `L>=k-n`.
The strongest is exactly `L>=-1`. Empty and full K add no restriction after
the fixed total. The sum of the k smallest entries of the displayed vector
is `binom(k,2)-1` for `1<=k<n`, and the full sum is `binom(n,2)`.
The same sorted-functional argument as for Z_n identifies the polytope.
Its vertex with a coordinate -1 does not lie in Z_n. QED.

For n=2 the vector is `(-1,2)`; for n=1 the slice and the embedded copy are
both a single point, so the strict-containment statement is not applicable.
A condition on two scores therefore does not recover the gate or even
identify the smaller embedded score polytope. The next result is stronger:
gate status already fails to be determined *at the very same point* of that
smaller copy.

## 3. A triangle in the score kernel erases the gate

For a nonempty tournament T and a chosen vertex v, let `U(T,v)` be obtained
from Gamma(T) by reversing the directed triangle `a -> b -> v -> a`.
Only those three arcs change. Thus in U,

    a -> v -> b -> a,
    b -> (T-v) -> a,

while all core arcs, including those incident to v, are retained.

**Theorem 3 (exact gate-erasure law).** Let `n=|T|`.

1. `s(U(T,v)) = s(Gamma(T)) = (s(T)+1,1,n)`.
2. U(T,v) is strong for every n>=1.
3. With cycles counted modulo cyclic rotation,

       HC(Gamma(T)) = H(T),       HC(U(T,v)) = H(T-v).            (5)

   Here H(empty)=1 only handles the n=1 boundary.
4. For n>=3, U(T,v) has no degree-defined gate iff
   `0 < d_T^+(v) < n-1`. If v is a core source, its gate is a->v;
   if v is a core sink, its gate is v->b.

**Proof.** Each vertex of a reversed directed triangle loses one outgoing
arc and gains one, proving (1). If n>=2, for every x in T-v the paths
`x -> a -> v -> b` and `v -> b -> x` show mutual reachability with v;
a and b join the same strong component. If n=1, U is a directed triangle.

In U the only exit from a is v and the only entry into b is v. Every
Hamilton cycle therefore contains the consecutive segment `a,v,b`.
Its rest is exactly a Hamilton path of T-v from b back to a; conversely any
such path closes. The empty path gives the triangle when n=1. This gives
a literal bijection, and the original gate bijection gives the first
formula in (5).

For gate classification when n>=3, a has outdegree1 and b outdegree n.
A core vertex w has outdegree `d_T^+(w)+1`. A gate tail must be a or a core
sink; a gate head must be b or a core source. The only exit of a is v, so
it yields a gate exactly when v is a core source. A core sink w!=v has
unique exit a, which is not a gate head; if v is a core sink, its unique
exit is b and it yields a gate. There are no others. QED.

Every tournament of order n>=3 has a mixed vertex: at most one vertex
can have outdegree0 and at most one outdegree n-1. Hence *every* canonical
gate score point at orders N>=5 has both gated and nongated realizations.
Both are strong, so both have the full score permutohedron as minimal face.

The cycle count need not change: a transitive core has H(T)=H(T-v)=1.
No universal strict decrease or universal difference formula is asserted.
The nontrivial change is certified by the following smallest hostile.

### 3.1 Smallest labeled-score hostile

Use the cyclic core `T=(2,4,1)` and v=0. All rows are outgoing-neighbor
bitmasks, with vertices numbered 0,...,4 and a=3,b=4. The two tournaments are

    G = (10,12,9,16,7),     U = (18,12,9,1,14).

They have identical labeled scores `(2,2,2,1,3)`, but

    G: gate (3,4), HC=3, H=15;
    U: no gate,   HC=1, H=11.

Both are strong. Exhausting all labeled tournaments through order four
shows that gate status and H are constant on every labeled score fiber
there. Thus order five is the smallest possible labeled-score ambiguity
for either of these two observables. The H threshold was already known in
THM-466(v); its replay here is a control, not a new threshold theorem. This minimality claim is FINITE-EXACT,
not a claim inferred from these two examples alone. The complete order-five
fiber at the displayed score is recorded in RESULT.json: it has 14 labeled
realizers, H support {11,13,15}, and HC support {1,2,3}. Thus even cycle-count
parity is lost inside this single score fiber.

Consequently no function of the labeled score point, even supplied with
its minimal face and all its cut slacks, can recognize the canonical gate,
recover H, or recover HC. The destroyed coordinate is an oriented directed
triangle, not a missing scalar normalization.

## 4. The actual excluded numbers appear in the same score fibers

A directed four-cycle gives a stronger witness than the triangle: it can
expose a two-path response instead of an ordinary deletion response.

Choose distinct core vertices u,v with `v -> u` in T. In Gamma(T) reverse

    a -> b -> v -> u -> a,

and call the result W(T;u,v). Let T' be T with its one arc v->u reversed.
For a subset S containing u,v, let `H_(u,v)(T'[S])` count Hamilton paths
of T'[S] beginning at u and ending at v.

**Theorem 4 (four-cycle two-path law).** W is strong, has exactly the same
labeled score as Gamma(T), and

    HC(W) = sum_(S contains u,v) H_(u,v)(T'[S]) H(T'[V minus S]).    (6)

It has a gate iff u is a sink and v a source of the original core T; in
that case the gate is u->v.

**Proof.** Reversing any directed cycle preserves each score. In W,
`a -> u -> v -> b -> a` is a directed four-cycle; every other core vertex
x satisfies `b -> x -> a`, so W is strong. Its only exit from a is u,
and its only entry into b is v. Anchor a Hamilton cycle at a: it uniquely
splits into a u-to-v path on a subset S of the core, then b, then a path
on its complementary subset, followed by the closing arc to a. Empty
complement is permitted because b->a. All core arcs are T', and the
cross arcs needed for the two pieces are present. This proves the bijection
and (6). The degree candidates for a gate are as in Theorem3. Vertex a
cannot be a gate tail because its head u was not a core source (v->u).
A core sink other than u still has unique exit a, not a gate head. If u
was a core sink, its only exit is now v, which is a gate head exactly
when v was a core source. There are no other candidates. QED.

**Explicit score-fiber pair with HC=7.** Take the inherited order-five
H=11 hostile core from THM-466(iv), with rows

    T=(24,1,3,6,14), u=0, v=1.
    Gamma(T)=(56,33,35,38,46,64,31),
    W=(26,96,35,38,46,1,61).

Both score vectors are `(3,2,3,3,4,1,5)`. The first has gate(5,6),
HC=11 and H=111. The second is strong with no gate, HC=7 and H=91.

**Explicit score-fiber pair with HC=21.** Take

    T=(30,44,56,48,34,1), u=1, v=4.
    Gamma(T)=(94,108,120,112,98,65,128,63),
    W=(94,60,120,112,224,65,2,111).

Both score vectors are `(5,4,4,3,3,2,1,6)`. The first has gate(6,7),
HC=37 and H=389. The second is strong with no gate, HC=21 and H=293.

These counts are FINITE-EXACT, verified by anchored Hamilton-cycle DP,
independent complete cyclic-permutation enumeration, and the two-path
formula (6). The path counts are checked by subset DP and complete
permutation enumeration. No minimum-order claim is made for these two pairs.
They are not tournaments with H=7 or H=21: the exceptional counts are HC
in the ungated tournaments, while H remains an allowed odd count.

The gate spectrum's two omissions therefore disappear even if one retains
the exact labeled score, its full cut-slack function, and the macro SCC
face. This is a literal structural separation involving the two excluded
integers, not a numerical coincidence with unrelated counting objects.

## 5. What addition, multiplication, and the two gaps do not imply

There are two physical operations with different geometric effects:

- Dominance join `A -> B` puts the score point in a product face and gives
  `H(A -> B)=H(A)H(B)`.
- The canonical gate raises vertex order by two, transports an entire face
  into the interior, and has `HC(Gamma(T))=H(T)`.

Repeating the same dominance block doubles its vertex order but **squares**
its path count. For example `C3 -> C3` has six vertices and H=9, not H=6.
The odd index `q=(H-1)/2` obeys `q(A->B)=q(A)+q(B)+2q(A)q(B)`, so even
this natural quotient does not convert the product-face operation into
ordinary addition. Literal `H -> 2H` is incompatible with tournament
oddness. No general impossibility for richer multi-coordinate recursions
is claimed.

The inherited omissions 7 and 21 transfer to HC on the **gated** class via
(5). They do not transfer to all strong tournaments, to a score slice, or
to arbitrary orientations of the permutohedron. Theorems3-4 show concretely
why the score projection cannot retain that local class definition. The
full coherent orientation does retain it, because it retains every arc.
This is a structural obstruction to a proposed compression, not a new
proof of the two omissions or of spectrum completeness.

## 6. Historical repairs and limits of attribution

The dated March output `permutohedron_s202.out` contains a sound score-map
construction beside informal identifications that should not be reused.
Its score-fiber, face, and metagraph wording is repaired as follows:

1. A labeled score fiber is a set of tournaments above **one lattice point**;
   it is not a face. Proposition 1 determines the face containing that point.
2. The n! vertices corresponding to transitive tournaments collapse to one
   isomorphism class, not to the full tournament metagraph (THM-207).
3. A sorted score sequence is an orbit of lattice points, not every lattice
   point. Sorting restricts to a Weyl chamber; it does not lower the real
   dimension from n-1 to n-2.
4. General tournament Hamilton paths are not linear extensions of a poset
   formed by all its arcs: a directed triangle already defeats that typing.
   Session210's coherent sinks are the valid replacement.

These are corrections to historical descriptions, not retractions of the
actual score zonotope, transitive flip graph, or Session210 constructions.
No broad claim that all historical mistakes have been audited is made.

## 7. Verification

`python3 verify.py` emits deterministic JSON without modifying any source
file. Repeat with `python3 -O verify.py`; the bytes agree. Explicit exceptions
are used instead of assertions, and only the Python standard library is used.

- Every one of 1,099 labeled cores of orders1..5: HP subset DP versus full
  permutation scan; every cut slack versus literal outgoing-cut arcs; exact
  SCC/tight-set correspondence and SCC product law.
- All 5,405 chosen-vertex gate triangle reversals: exact scores, strongness
  for n>=2 (n=1 also follows directly), cycle counts by anchored subset DP,
  and the precise gate-loss classification for n>=3.
- Independent cyclic permutation scans on every chosen vertex for n<=4,
  plus all vertices of four declared order-five cores. No direct-permutation
  claim is made for the remaining order-seven targets.
- Exhaustive labeled-score fibers through order five, including the
  minimum-order checks and complete displayed order-five fiber.
- Every vertex of the relaxed degree slice for n=2..6, checked against all
  ambient facet inequalities, alongside all core-score interior checks.
- Every directed four-cycle choice above for every core of orders2..5:
  10,650 switches; scores, strongness, precise gate classification, and
  anchored cycle DP versus the independently enumerated two-path formula.
- Two specified seven/eight-vertex score-fiber pairs realizing HC=7,21:
  direct cyclic words, subset DP and two-path decomposition all agree;
  path counts have a separate full permutation check.
- A physical order-doubling product-face hostile `C3 -> C3`.

The all-order claims follow from the proofs. Finite controls check the
implementation and conventions, and do not establish the Chan–Pak conjecture.
