# Independent audit of production reverse-potential counting

2026-10-25. **REVIEW PASS:** the proof and production implementation in
[README.md](README.md), [reverse.py](reverse.py), and [check.py](check.py)
are sound in their stated ordinary, zero-excess scope. No mathematical
gap or required producer repair was found. This review accepts recurrence
equivalence, normalized residue computation and safe pruning. It does
not claim smaller structural graphs, faster runtime, or an end-to-end
bit-complexity improvement.

## 1. Inheritance and graph semantics

Read the complete producer and the inherited reverse-potential proof
and review, integer-gauge proof, and termwise-baseline proof. The closest
mechanism is the ordinary defect transfer. Its genuine hostile is an
absolute coefficient that vanishes at low precision while its normalized
coefficient remains visible. The corrected near miss is using that
absolute support to build the graph. The less-used sidecar is stage:
a zero batch still advances the available cycle lengths.

The production constructor uses only N, stage, remaining degree and
the exact retained defect tuple. It reads no coefficient, residue or
nonzero-support test. Each generated vertex is reachable. Completing
the remaining degree with singletons gives a suffix from every vertex,
so every vertex reaches the unique terminal. Vertex IDs are assigned
layer by layer; every edge has a strictly larger target ID and a
strictly smaller stage. The reverse traversal therefore computes genuine
finite suffix minima, including early zero-remainder states and N=0.

Dropping a coordinate because its odd argument exceeds either current
stage or remaining degree is safe: neither bound increases later.
The code asserts the exact affordable mask. The stage is retained in
State; repeated(S,D) values at different stages are not merged.

The state merger is substantive. The N=26 prefixes(23) and(9,9,5)
have different cycle counts but the same stage3,S=3,D=(0,2).
For any future cycle choice, both have the same available multiplicities,
the same defect update, and the same already normalized transfer.
Appending any available suffix to either prefix is legitimate.
Thus the suffix-minimum argument survives actual merging; it does not
assume that a quotient preserving only a final total also preserves paths.

Each complete path specifies one multiplicity at every odd stage,
hence exactly one odd partition. Conversely each odd partition supplies
exactly one complete path. A merger does not create duplicate choices:
from one vertex, different k give different remaining degrees, and
the past and future stage ranges are disjoint. This is why the path
valuation can be compared to the termwise Burnside baseline.

## 2. Original edge transfer and pathwise normalization

The scalar gauge is Y=O(u)2^[S(m-b)]w. For a batch of k cycles of odd
length ell, the original cross contribution is k(ell*m-D(ell)).
The gauge ratio contributes

    S'(m+k-b')-S(m-b).

Since S-S'=ell*k, the m terms cancel. Splitting the two-part of k!
then gives exactly the implemented cost

    J=A_ell(k)-v2(k!)+k(ell*b-D(ell))+S'(k+b-b').

The inherited bounds give J>=0, including zero batches, coordinate
deletion, ell>S and the ell=1 terminal step. Empty masks and zero
factorials are handled correctly.

The factor T=O(u')/[O(u)ell^k O(k)] is an integer before modular
reduction: u'!/[u!ell^k k!] counts choosing and arranging k disjoint
ell-cycles on new positions. Its odd part is the displayed quotient.
All numerator and denominator factors are odd, so T is positive odd.
The exact replay checks that division; the modular implementation uses
the same proved quotient via an odd inverse. No exact-factorial call
occurs in residue().

The reverse minimum V satisfies V(v)<=J(v,w)+V(w), so every reduced
cost J'=J+V(w)-V(v) is nonnegative. For an individual prefix of cost c,
a minimizing suffix makes c+V(v)>=B. Hence its transformed contribution
is an odd integer times2^[c+V(v)-B], integral before merging.
All contributions at v use the same V(v), so summation commutes with
the coordinate change. The implementation starts Z_initial=1 and
never divides a truncated absolute coefficient by a power of two.

This is a gauge on the reachable path lattice, not arbitrary independent
integer assignments to the Y coordinates. The producer states that
limitation explicitly.

Along a complete path, all intermediate gauges telescope. The product
of odd factors is O(N)/product[d^k_d O(k_d)], and the total original
edge cost is e(lambda)-sum v2(k_d!). The inherited termwise proof gives
B=floor(N/2) for N>=3, with the distinct odd pair(1,N-1) providing the
even equality case. N=0,1,2 correctly have B=0. The constructor actually
computes V before checking this equality. A later carry in the sum can
increase its valuation; the N=50 hostile correctly separates25 from26.

## 3. Production evaluation and pruning

At the terminal, Z=O(N)a(N)/2^B. Multiplying by the odd inverse of O(N)
returns the desired normalized residue. Every stored recurrence residue
is reduced modulo2^K. The use of modular odd-factorial values in the
coefficient cache is valid because reduction commutes with multiplication
and inversion of odd integers. Cache keys(u,ell,k) contain every input
on which T depends; no defect or history coordinate is missing.

An edge with J'>=K has zero coefficient modulo2^K for every source,
so skipping it before computing T is valid. A source residue is tested
for zero only after its whole preceding layer has been accumulated.
Thus zero-source skipping cannot discard a later incoming contribution.
The structural graph and V were already computed and remain intact.

The unpruned route processes every edge, while the edge-only route
isolates reduced-cost pruning. In the unpruned reference, a shift>=K
is evaluated as a zero modular contribution without creating the huge
integer2^shift; this is the exact modular operation, not an omitted path.
The sum of edge updates, cost-pruned edges and source-zero edges equals
all scanned edges, with cost pruning taking priority when both apply.

The absolute comparison requests B+K bits on the same graph and targets
the same final normalized K bits. Its final residue is divisible by2^B.
Comparing to an absolute K-bit run would be invalid; the N=7 hostile
explicitly demonstrates that loss.

## 4. Counter and metadata audit

The reported structural counts include terminal vertices, all edges,
stored defect coordinates, actual mergers, stage distinctions and one
reverse relaxation per edge. The full structure and V remain stored.
Only recurrence residues use two consecutive layers. Thus reduced
residue precision is not presented as a reduction of total graph storage.

The table compares the same DAG at normalized K versus absolute B+K.
The edge-work counts and operand-bit sums are deterministic observations,
not timings. Their scopes are:

- multiply_operand_bit_sum adds bit lengths of the source residue and
  cached odd factor at each performed **edge weight update**;
- maximum product, shifted-product and modular-sum widths likewise
  describe these edge updates;
- odd-factorial setup products, odd inversions and transfer-cache entries
  have their own counters;
- the final multiplication by the inverse of O(N), coefficient-setup
  arithmetic, potential/graph integer widths, allocation overhead and
  instrumentation are not included in those edge-product widths or in
  a combined end-to-end cost;
- layer_residue_bit_sum accumulates across layers; it is not peak memory.

The code and the README's no-complexity/no-runtime claim are consistent
with this accounting. A future performance claim would need to charge
these omitted costs as well.

## 5. Reproduced producer controls

Both normal and optimized runs passed, without writing producer files:

    python 05-knowledge/results/lca_ports_reverse_20261025/counting/check.py
    python -O 05-knowledge/results/lca_ports_reverse_20261025/counting/check.py

They reproduce408 normalized comparisons,153 pruned/unpruned comparisons,
12,589 intermediate coordinate checks, and3,725 complete partition-path
checks. The samples at N80,100 and all declared hostiles pass. Exception
gates remain active under -O.

## 6. Independent controls that do not reuse the transfer derivation

The following embedded control was run normally and with -O.

First, it applies every permutation directly to unordered vertex pairs.
It checks orientation parity around each actual pair orbit, returning
zero for an inconsistent orbit and2^r otherwise. Burnside averaging then
counts all tournament isomorphism classes through order8. It does not
use the odd-cycle gcd exponent, partition formula, defect state,
factorial gauge, or producer's exact counting implementation as its
oracle. Across all46,234 permutations at orders0..8 it recovers

    1,1,1,2,4,12,56,456,6880,

and agrees with exact and normalized producer outputs.

Second, independent coefficients of product_(d odd)(1-x^d)^(-1)
count odd partitions. They equal the number of complete structural
paths for every N0..50, including3,658 paths at N50. This checks path
multiplicity independently of all costs and arithmetic weights.
The graph, edges and V are unchanged by both low-precision pruned
evaluation and a higher-precision unpruned run.

Run the code below from the repository root. It uses only standard
library code plus the production API being tested.

```python
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import sys, json
p = Path('05-knowledge/results/lca_ports_reverse_20261025/counting').resolve()
sys.path.insert(0, str(p))
from reverse import OrdinaryDAG

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def fixed_tournaments(perm):
    n = len(perm)
    unseen = set(combinations(range(n), 2))
    orbits = 0
    while unseen:
        start = min(unseen)
        pair = start
        parity = 0
        while True:
            unseen.remove(pair)
            u, v = perm[pair[0]], perm[pair[1]]
            parity ^= int(u > v)
            pair = (min(u, v), max(u, v))
            if pair == start:
                break
        if parity:
            return 0
        orbits += 1
    return 1 << orbits

observed = []
perm_count = 0
for n in range(9):
    total = 0
    for perm in permutations(range(n)):
        total += fixed_tournaments(perm)
        perm_count += 1
    count, remainder = divmod(total, factorial(n))
    need(remainder == 0, 'literal Burnside integrality')
    graph = OrdinaryDAG(n)
    need(graph.exact()[0] == count, ('direct permutation count', n))
    for K in (1, 3, 7):
        need(graph.residue(K)[0] == (count >> graph.B) % (1 << K),
             ('direct permutation normalized bits', n, K))
    observed.append(count)
need(observed == [1,1,1,2,4,12,56,456,6880], 'small counts')

# Independent formal-product coefficients count odd partitions, without
# enumerating a cycle partition or using costs/defects/transfer factors.
limit = 50
partitions = [0] * (limit+1)
partitions[0] = 1
for odd in range(1, limit+1, 2):
    for n in range(odd, limit+1):
        partitions[n] += partitions[n-odd]
for n in range(limit+1):
    graph = OrdinaryDAG(n)
    paths = [0]*len(graph.states)
    paths[0] = 1
    for v, edges in enumerate(graph.edges):
        for edge in edges:
            paths[edge.target] += paths[v]
    need(paths[graph.terminal] == partitions[n], ('complete-path count', n))
    before = (tuple(graph.states),
              tuple(tuple(row) for row in graph.edges),
              tuple(graph.potential))
    graph.residue(1)
    graph.residue(7, prune_edges=False, skip_zero_sources=False)
    after = (tuple(graph.states),
             tuple(tuple(row) for row in graph.edges),
             tuple(graph.potential))
    need(before == after, ('residues must not alter structural graph', n))
print(json.dumps({'status':'PASS','literal_permutations':perm_count,
                  'direct_counts':observed,'path_count_orders':[0,limit],
                  'odd_partition_paths_at50':partitions[50],
                  'structural_immutability_orders':[0,limit]}))
```

Both runs return PASS with the same deterministic counts. These bounded
controls supplement the all-order proof; they do not turn the observed
operation savings into an asymptotic claim.

**Disposition:** accept the production ordinary reverse-potential
integration and its stated counters. The signed clone-excess engine,
its cancellation semantics, and unrestricted state compression are
outside this implementation's scope and remain unchanged.
