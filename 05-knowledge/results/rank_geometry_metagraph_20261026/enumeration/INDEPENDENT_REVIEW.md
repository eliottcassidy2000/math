# Independent review of canonical augmentation and clone reconstruction

2026-10-26. **ACCEPTED within the declared scopes.** Reviewer
`/root/nonconsecutive_ports` independently read [README.md](README.md),
[enumerate.py](enumerate.py), [check.py](check.py), both receipts, and
[THM-4463, homogeneous tournament contacts, §2](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md).
No proof or implementation defect was found. This is a mathematical/code
review and exact finite replay, not a formal proof-assistant certificate.
The all-order result is correctness of a terminating enumerator; the
executed census stops at eight. No smaller scalar recurrence, arbitrary-order
performance bound better than the stated coarse factorial bound, or
black-core connectivity result is promoted.

## 1. Canonical deletion really selects one unrooted class

The individualization/refinement tree is equivariant: ordered initial
color cells, outgoing-count signatures, their lexicographic order, and
the first nonsingleton cell all transport under an isomorphism. Within
that cell every possible singleton is visited. The implementation does
not discard branches by an unproved orbit rule. Its least leaf need not
be the least code over all permutations; that distinction in the proof
is correct and necessary.

All leaf orders attaining the minimum are retained. Fixing one order
`p0`, the maps `p0[i] -> p[i]` are exactly the automorphisms on the input
labels. Both directions hold: equal adjacency codes give automorphisms,
and equivariance sends the fixed leaf to a minimum leaf for every
automorphism. The transported canonical-label action in the implementation
uses `inverse_p0[p[i]]`, which is the correct conjugation. No stabilizer
or generator closure is silently assumed.

Consequently `{p[n-1]}` is one vertex orbit and transports intrinsically.
Parent plus extension-subset orbit classifies a **rooted** child exactly:
an isomorphism fixing the new vertex restricts to a parent automorphism,
and a parent automorphism carrying the subset extends to such an
isomorphism. Choosing roots in the single intrinsic deletion orbit accepts
exactly one rooted class above each unrooted class. The duplicate exception
in `enumerate_classes` tests the uniqueness conclusion instead of repairing
it by deduplication. Empty and one-vertex cases are handled literally.

The recurrence state is the canonical parent and its complete action.
The root-orbit number is not uniformly n; C3 versus TT3 is the decisive
minimal hostile. The implementation's exact subset orbit action is the
required sidecar, not merely the group order.

## 2. Pair-free cores and all clone mergers

The inherited theorem's full proof was read, not only its status. A
vertex cannot have two incoming or two outgoing homogeneous neighbors;
a directed cycle of homogeneous contacts is impossible. Along a contact
path, homogeneity forces transitivity and constant incidence to every
outside vertex. Conversely consecutive members of any transitive module
are homogeneous in the entire tournament. Thus connected components of
the pair graph are precisely the maximal transitive modules and form a
canonical partition. Homogeneous pairs themselves are not an equivalence
relation; TT3 remains an explicit hostile.

If the quotient had a homogeneous pair, the union of its two corresponding
blocks would be both homogeneous and transitive: their interblock arcs
all point in one direction, and every other block sees them alike. This
contradicts maximality. Conversely a homogeneous pair crossing two blocks
of an inflation of a pair-free core would force the two quotient vertices
to agree against every other quotient block. That would be a homogeneous
pair in the core. Internal consecutive vertices do remain homogeneous,
so the inflated blocks are exactly the canonical components. This proves
both the no-merger claim for pair-free inflations and the need to merge
arbitrary transitive block partitions first.

An isomorphism preserves this canonical partition, its block cardinalities,
and each rigid transitive internal order. It is therefore exactly a
weight-preserving core isomorphism with its unique internal lift. The
weight-orbit formula is correct with the **actual** Aut(Q) action. An
unordered multiset of weights is insufficient, as opposite weight orders
on C3 demonstrate. `weighted_key` retains both the adjacency code and the
transported ordered colors, and `weighted_core_allocations` uses one least
tuple per actual group orbit. The order-N all-one allocation still needs
all order-N pair-free cores; the stated remaining cost is honest.

An additional independent audit below enumerates every vertex subset and
tests the module predicate and absence of cyclic triples directly. Its
maximal subsets agree with the producer's pair components for every
labelled tournament through five. This is a different route to the
partition, rather than an inflation/deflation round trip alone.

## 3. Framed operations and scope

The insertion argument supplies Hamiltonian paths, including the empty
convention. The complete-word Aut(T) action is free, so chart multiplicity
is H/Aut. With ordered repeated ports the same argument uses the pointwise
stabilizer. Port occurrence vectors and `(block,internal rank)` retain
coincident endpoints and distinct ranks inside transitive blocks.

K flips exactly the nonspine arcs in the retained frame. R takes converse
and reverses the frame, which reflects the free-pair bits without another
bit complement. Their word rules and literal arc rules agree. The explicit
order-four K hostile shows why K is not a function of the unmarked class.
No assertion about connectivity of a selected blue/black move set follows
from this complete presentation inventory. These boundaries are preserved
in the prose and the code.

## 4. Replays and independent controls

Both full producer commands were independently rerun successfully:

```powershell
python -B 05-knowledge/results/rank_geometry_metagraph_20261026/enumeration/check.py
python -B -O 05-knowledge/results/rank_geometry_metagraph_20261026/enumeration/check.py
```

They each took about 12 seconds in this review and reproduced every
declared finite control through eight. There was no independent full
factorial census at order eight; the brute permutation controls stop at
five, exactly as stated in the producer.

The two committed receipts are byte-identical CRLF files, SHA256
`ec87ed60f39a3146127a1cf2454ae7e5efb53aeaaa985ef05303e8ecd95cfffd`.
After CRLF-to-LF normalization their SHA256 is
`2678b28367b78bfb478bc61a388870be82d71ad1a93f0a57d1a9511f232d16d3`.
The freshly captured normal and optimized output has that same normalized
hash. This is a newline distinction, not a difference in the census.

The following review control imports the producer only as the object
under test. Its adjacency matrix, transitive-module test, exhaustive
automorphism test, brute unmarked keys, and rooted keys are independent.
It does not import the producer checker or the reverse recurrence.
Every condition uses an exception and survives optimization.

```python
from pathlib import Path
import importlib.util, sys, itertools, json

root = Path("05-knowledge/results/rank_geometry_metagraph_20261026/enumeration")
spec = importlib.util.spec_from_file_location("audit_subject", root / "enumerate.py")
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def matrix(n, code):
    a = [[False] * n for _ in range(n)]
    for k, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        a[u][v] = bool(code >> k & 1)
        a[v][u] = not a[u][v]
    return a

def word(a, order):
    return tuple(a[u][v] for u, v in itertools.combinations(order, 2))

def brute(a):
    return min(word(a, p) for p in itertools.permutations(range(len(a))))

levels, receipts = m.enumerate_classes(6)
module_probes = module_count = deletion_probes = labelled = 0
rooted_counts = []
for n in range(6):
    brute_to_child, rooted = {}, set()
    for code in range(1 << (n * (n - 1) // 2)):
        labelled += 1
        a, vertices = matrix(n, code), set(range(n))
        modules = []
        for mask in range(1, 1 << n):
            subset = {v for v in vertices if mask >> v & 1}
            module = all(len({a[w][v] for v in subset}) == 1
                         for w in vertices - subset)
            acyclic = all(
                not (a[u][v] and a[v][w] and a[w][u]) and
                not (a[v][u] and a[u][w] and a[w][v])
                for u, v, w in itertools.combinations(subset, 3))
            module_probes += 1
            if module and acyclic:
                modules.append(frozenset(subset))
        maximal = {s for s in modules if not any(s < t for t in modules)}
        blocks, q, weights = m.clone_core(n, code)
        need(maximal == {frozenset(block) for block in blocks},
             ("maximal modules", n, code))
        module_count += len(maximal)
        need(sum(map(len, maximal)) == n and
             set().union(*maximal) == vertices, ("partition", n, code))
        qa = matrix(len(blocks), q)
        need(not any(all(qa[u][w] == qa[v][w]
                         for w in range(len(blocks)) if w not in (u, v))
                     for u, v in itertools.combinations(range(len(blocks)), 2)),
             ("pair-free", n, code))
        cf = m.canonical(n, code)
        group = [p for p in itertools.permutations(range(n))
                 if all(a[u][v] == a[p[u]][p[v]]
                        for u, v in itertools.combinations(range(n), 2))]
        need(set(group) == set(cf.automorphisms()), ("full group", n, code))
        if n:
            chosen = cf.orders[0][-1]
            need({p[chosen] for p in group} == {p[-1] for p in cf.orders},
                 ("deletion orbit", n, code))
        key = brute(a)
        need(brute_to_child.setdefault(key, cf.code) == cf.code,
             ("false split", n, code))
        for v in range(n):
            remaining = [u for u in range(n) if u != v]
            rooted.add(min(word(a, p + (v,))
                           for p in itertools.permutations(remaining)))
            deletion_probes += 1
    need(set(brute_to_child.values()) == set(levels[n]) and
         len(brute_to_child) == len(levels[n]), ("all classes", n))
    rooted_counts.append(len(rooted))
    if n:
        need(len(rooted) == receipts[n - 1]["rooted_extensions"],
             ("rooted class census", n))

# The order-six outputs are checked for literal lift, without claiming
# an independent unrestricted labelled census at six.
for n, level in enumerate(levels):
    for code in level:
        blocks, q, weights = m.clone_core(n, code)
        aa, bb = matrix(n, code), matrix(n, m.inflate(q, weights))
        order = sum(blocks, ())
        need(all(aa[order[u]][order[v]] == bb[u][v]
                 for u in range(n) for v in range(n)), ("literal lift", n, code))
print(json.dumps({
    "status": "PASS", "labelled_tournaments": labelled,
    "module_subset_probes": module_probes,
    "maximal_transitive_modules": module_count,
    "literal_root_markings": deletion_probes,
    "rooted_class_counts": rooted_counts,
    "class_counts_through6": [len(x) for x in levels],
    "accepted_lifts": sum(map(len, levels))}, sort_keys=True))
```

Extract the unique Python block and run it from the repository root, with
both normal and optimized Python. Both modes were executed in this review
and returned the identical exact output:

```json
{"accepted_lifts": 77, "class_counts_through6": [1, 1, 1, 2, 4, 12, 56], "labelled_tournaments": 1100, "literal_root_markings": 5405, "maximal_transitive_modules": 4015, "module_subset_probes": 32767, "rooted_class_counts": [0, 1, 2, 4, 12, 48], "status": "PASS"}
```

This additional control checks maximal modules directly, all labelled
deletion orbits through five, independent rooted class counts, and exact
inflation lifts of the 77 generated classes through six. It supports the
proofs above without expanding the declared large-order universe.
