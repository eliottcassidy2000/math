#!/usr/bin/env python3
"""Exact finite receipt-network grounding and boundary-energy controls.

Standard library only. No imported-paper theorem is tested, and no unbounded
Collatz orbit search is used. Run identically with or without Python -O.
"""
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations


CHECKS = Counter()


def require(test, why):
    if not test:
        raise ValueError(why)


def check(test, label):
    CHECKS[label] += 1
    if not test:
        raise RuntimeError(f"{label}, check {CHECKS[label]}")


def exact_int(x, positive=False):
    require(type(x) is int and (not positive or x > 0), "exact integer required")


def solve(matrix, rhs):
    """Exact possibly singular compatible solve; free coordinates set to zero."""
    n = len(rhs)
    require(len(matrix) == n and all(len(row) == n for row in matrix), "square system")
    require(all(type(x) is F for row in matrix for x in row)
            and all(type(x) is F for x in rhs), "Fraction system required")
    a = [list(row)+[rhs[i]] for i, row in enumerate(matrix)]
    pivots = []
    cursor = 0
    for col in range(n):
        pivot = next((i for i in range(cursor, n) if a[i][col]), None)
        if pivot is None:
            continue
        a[cursor], a[pivot] = a[pivot], a[cursor]
        q = a[cursor][col]
        a[cursor] = [x/q for x in a[cursor]]
        for i in range(n):
            if i != cursor and a[i][col]:
                q = a[i][col]
                a[i] = [x-q*y for x, y in zip(a[i], a[cursor])]
        pivots.append((cursor, col))
        cursor += 1
    require(all(any(row[:-1]) or not row[-1] for row in a), "inconsistent system")
    x = [F(0)]*n
    for row, col in pivots:
        x[col] = a[row][-1]
    return x


@dataclass(frozen=True)
class Network:
    """Abstract finite undirected network; edge orientation is bookkeeping.

    Plain integer labels do NOT authenticate Collatz facts. The Receipt
    verifier below supplies that additional interface in the actual controls.
    Parallel edges are allowed and each has its own listed conductance.
    """
    vertices: tuple
    edges: tuple

    def __post_init__(self):
        require(type(self.vertices) is tuple and type(self.edges) is tuple, "tuple network required")
        for v in self.vertices:
            exact_int(v)
        require(len(set(self.vertices)) == len(self.vertices), "distinct labels required")
        for edge in self.edges:
            require(type(edge) is tuple and len(edge) == 3, "typed edge triple required")
            a, b, c = edge
            exact_int(a); exact_int(b)
            require(a in self.vertices and b in self.vertices and a != b, "distinct listed endpoints")
            require(type(c) is F and c > 0, "positive Fraction conductance required")

    def vector(self, h):
        require(type(h) is dict and all(type(v) is int for v in h), "exact labelled vector")
        require(set(h) == set(self.vertices) and all(type(x) is F for x in h.values()),
                "complete Fraction vector required")

    def energy(self, h):
        self.vector(h)
        return sum((c*(h[a]-h[b])**2 for a, b, c in self.edges), F(0))

    def laplace(self, h):
        self.vector(h)
        out = {v: F(0) for v in self.vertices}
        for a, b, c in self.edges:
            flow = c*(h[a]-h[b])
            out[a] += flow
            out[b] -= flow
        return out

    def adjacency(self):
        out = {v: [] for v in self.vertices}
        for i, (a, b, _) in enumerate(self.edges):
            out[a].append((b, i, 1))
            out[b].append((a, i, -1))
        return out


def terminals(net, source, root):
    require(type(net) is Network, "Network required")
    exact_int(source); exact_int(root)
    require(source in net.vertices and root in net.vertices and source != root, "distinct listed terminals")


def path_flow(net, source, root):
    terminals(net, source, root)
    adjacency = net.adjacency()
    parents = {source: None}
    queue = deque([source])
    while queue:
        a = queue.popleft()
        for b, i, sign in adjacency[a]:
            if b not in parents:
                parents[b] = (a, i, sign)
                queue.append(b)
    if root not in parents:
        return None, {v: F(int(v in parents)) for v in net.vertices}
    current = root
    flow = [F(0)]*len(net.edges)
    path = [root]
    while current != source:
        a, i, sign = parents[current]
        flow[i] = F(sign)
        current = a
        path.append(current)
    return tuple(flow), tuple(reversed(path))


def boundary(net, flow):
    require(type(flow) is tuple and len(flow) == len(net.edges)
            and all(type(x) is F for x in flow), "complete Fraction edge flow required")
    out = {v: F(0) for v in net.vertices}
    for (a, b, _), x in zip(net.edges, flow):
        out[a] -= x
        out[b] += x
    return out


def flow_energy(net, flow):
    boundary(net, flow)
    return sum((x*x/c for (_, _, c), x in zip(net.edges, flow)), F(0))


def trial_capacity(net, source, root, h):
    """Exact energy correction, allowing floating disconnected components."""
    terminals(net, source, root)
    net.vector(h)
    require(h[source] == 1 and h[root] == 0, "pinned terminal values required")
    interior = [v for v in net.vertices if v not in (source, root)]
    residual = net.laplace(h)
    matrix = [[F(0) for _ in interior] for _ in interior]
    index = {v: i for i, v in enumerate(interior)}
    for a, b, c in net.edges:
        if a in index:
            matrix[index[a]][index[a]] += c
        if b in index:
            matrix[index[b]][index[b]] += c
        if a in index and b in index:
            matrix[index[a]][index[b]] -= c
            matrix[index[b]][index[a]] -= c
    y = solve(matrix, [residual[v] for v in interior])
    correction = {v: F(0) for v in net.vertices}
    correction.update(zip(interior, y))
    harmonic = {v: h[v]-correction[v] for v in net.vertices}
    defect = sum((residual[v]*correction[v] for v in interior), F(0))
    return {"capacity": net.energy(h)-defect, "trial_energy": net.energy(h),
            "dual_residual_squared": defect, "residual": residual,
            "harmonic": harmonic, "correction": correction}


def eliminate(net, vertex, protected=()):
    """One exact star-mesh elimination; proof-path provenance is an external DAG."""
    exact_int(vertex)
    require(type(protected) is tuple and all(type(v) is int for v in protected), "exact protected labels")
    require(vertex in net.vertices and vertex not in protected, "only unprotected interior may be removed")
    conductance = {}
    neighbours = {}
    for a, b, c in net.edges:
        if vertex in (a, b):
            u = b if a == vertex else a
            neighbours[u] = neighbours.get(u, F(0))+c
        else:
            key = tuple(sorted((a, b)))
            conductance[key] = conductance.get(key, F(0))+c
    total = sum(neighbours.values(), F(0))
    if total:
        for a, b in combinations(sorted(neighbours), 2):
            key = (a, b)
            conductance[key] = conductance.get(key, F(0))+neighbours[a]*neighbours[b]/total
    return Network(tuple(v for v in net.vertices if v != vertex),
                   tuple((a, b, c) for (a, b), c in sorted(conductance.items())))


def word_checked(word):
    require(type(word) is tuple and all(type(a) is int and a >= 1 for a in word), "positive exact valuation tuple")


def replay(n, word):
    exact_int(n, positive=True)
    require(n % 2 == 1, "positive odd source required")
    word_checked(word)
    x = n
    for a in word:
        require(x != 1, "ROOT padding prohibited")
        numerator = 3*x+1
        actual = (numerator & -numerator).bit_length()-1
        require(actual == a, "false valuation receipt")
        x = numerator >> a
    return x


@dataclass(frozen=True)
class Receipt:
    left: int
    right: int
    left_word: tuple
    right_word: tuple


def verify_receipt(receipt):
    require(type(receipt) is Receipt, "Receipt required")
    left = replay(receipt.left, receipt.left_word)
    right = replay(receipt.right, receipt.right_word)
    require(left == right, "unmatched actual endpoints")
    return len(receipt.left_word), len(receipt.right_word)


def then_clock(ab, cd):
    for pair in (ab, cd):
        require(type(pair) is tuple and len(pair) == 2
                and all(type(x) is int and x >= 0 for x in pair), "exact nonnegative clock pair")
    a, b = ab
    c, d = cd
    return a+max(c-b, 0), d+max(b-c, 0)


def authenticate_path(receipts, source, root=1):
    """Compose matched receipt edges and export only a bounded first-hit word."""
    exact_int(source, positive=True); exact_int(root, positive=True)
    require(root == 1 and source % 2 == 1, "odd source and actual ROOT1 required")
    require(type(receipts) is tuple, "ordered receipt tuple required")
    current, clocks = source, (0, 0)
    for item in receipts:
        require(type(item) is tuple and len(item) == 2, "typed path step required")
        rec, forward = item
        require(type(forward) is bool, "Boolean direction required")
        pair = verify_receipt(rec)
        start, end = (rec.left, rec.right) if forward else (rec.right, rec.left)
        require(current == start, "wrong labelled seam")
        clocks = then_clock(clocks, pair if forward else pair[::-1])
        current = end
    require(current == root, "path is not grounded")
    x, word = source, []
    for _ in range(clocks[0]):
        if x == 1:
            break
        numerator = 3*x+1
        a = (numerator & -numerator).bit_length()-1
        word.append(a)
        x = numerator >> a
    require(x == 1, "clock certificate did not reach ROOT")
    return clocks, tuple(word)


def export_from_receipts(receipts, source):
    """ROOT word from the stored component; None means no STORED path only."""
    require(type(receipts) is tuple, "exact receipt tuple required")
    exact_int(source, positive=True)
    require(source % 2 == 1, "odd source required")
    for rec in receipts:
        verify_receipt(rec)
    if source == 1:
        return ()
    # Zero-length geometric self-edges cannot enlarge a reachability component.
    active = tuple(rec for rec in receipts if rec.left != rec.right)
    labels = {1, source}
    for rec in active:
        labels.update((rec.left, rec.right))
    net = Network(tuple(sorted(labels)), tuple((rec.left,rec.right,F(1)) for rec in active))
    flow, path = path_flow(net, source, 1)
    if flow is None:
        return None
    steps = []
    for a, b in zip(path, path[1:]):
        i = next(i for i, rec in enumerate(active) if flow[i] and {rec.left,rec.right} == {a,b})
        steps.append((active[i], flow[i] > 0))
    return authenticate_path(tuple(steps),source)[1]


def chain(length, connected=True, conductance=F(1)):
    # Root is -1; source is0. A disconnected chain has a separate final vertex.
    vertices = (-1,)+tuple(range(length+1)) if not connected else (-1,)+tuple(range(length))
    edges = [(i, i+1, conductance) for i in range(length-1)]
    edges.append((length-1, -1 if connected else length, conductance))
    return Network(vertices, tuple(edges))


def main():
    graph_count = 0
    for n in range(2, 6):
        possible = list(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            graph_count += 1
            edges = tuple((a, b, F(1)) for i, (a, b) in enumerate(possible) if mask & (1 << i))
            net = Network(tuple(range(n)), edges)
            source, root = 0, n-1
            flow, path_or_cut = path_flow(net, source, root)
            for choice in (0, 1):
                h = {v: F((2*v+choice) % 5, 3) for v in net.vertices}
                h[source], h[root] = F(1), F(0)
                data = trial_capacity(net, source, root, h)
                cap = data["capacity"]
                check(cap >= 0 and data["dual_residual_squared"] >= 0, "nonnegative_energy")
                check(cap == net.energy(data["harmonic"]), "residual_identity")
                check(data["dual_residual_squared"] == net.energy(data["correction"]), "residual_identity")
                check(all(x == 0 for v, x in net.laplace(data["harmonic"]).items()
                          if v not in (source, root)), "harmonic_equations")
                check((cap > 0) == (flow is not None), "grounding_iff")
                if flow is not None:
                    expected = {v: F(int(v == root)-int(v == source)) for v in net.vertices}
                    check(boundary(net, flow) == expected, "unit_boundary")
                    check(cap*flow_energy(net, flow) >= 1, "Thomson_lower_certificate")
                else:
                    cut = path_or_cut
                    check(cut[source] == 1 and cut[root] == 0 and net.energy(cut) == 0, "separating_cut")
            for order in (list(range(1, n-1)), list(reversed(range(1, n-1)))):
                small = net
                for v in order:
                    small = eliminate(small, v, (source, root))
                check(small.energy({source:F(1), root:F(0)}) == cap, "Kron_boundary_preservation")

    for length in range(1, 33):
        loose = chain(length, connected=False)
        h = {-1:F(0), **{i:F(length-i, length) for i in range(length+1)}}
        data = trial_capacity(loose, 0, -1, h)
        check((data["capacity"], data["trial_energy"], data["dual_residual_squared"])
              == (0, F(1, length), F(1, length)), "small_residual_hostile")
        interior = [v for v in loose.vertices if v not in (0, -1)]
        check(max(abs(data["residual"][v]) for v in interior) == F(1, length), "small_residual_hostile")
        connected = chain(length)
        h = {-1:F(0), **{i:F(length-i, length) for i in range(length)}}
        check(trial_capacity(connected, 0, -1, h)["capacity"] == F(1, length), "refinement_resistance")
        priced = chain(length, conductance=F(length))
        check(trial_capacity(priced, 0, -1, h)["capacity"] == 1, "refinement_resistance")

    for length in range(1, 8):
        vertices = [0, 1]
        edges = []
        next_label = 2
        for j in range(length):
            path = [0]+list(range(next_label, next_label+length-1))+[1]
            vertices.extend(path[1:-1])
            next_label += length-1
            edges.extend((a, b, F(1)) for a, b in zip(path, path[1:]))
        net = Network(tuple(vertices), tuple(edges))
        h = {v:F(0) for v in vertices}; h[0] = F(1)
        check(trial_capacity(net, 0, 1, h)["capacity"] == 1, "capacity_has_no_uniform_path_deadline")
        _, path = path_flow(net, 0, 1)
        check(len(path)-1 == length, "capacity_has_no_uniform_path_deadline")

    # Formal graph only, deliberately not authenticated as a Collatz cycle.
    cycle = Network((0, 1, 2, 3), ((1,2,F(1)),(2,3,F(1)),(3,1,F(1))))
    check(all(x == 0 for x in boundary(cycle, (F(1),F(1),F(1))).values()), "cycle_has_zero_boundary")
    check(path_flow(cycle, 1, 0)[0] is None, "cycle_has_zero_boundary")

    r73 = Receipt(7, 3, (1,1,2,3), (1,))
    r35 = Receipt(3, 5, (1,), ())
    r51 = Receipt(5, 1, (4,), ())
    clocks, root_word = authenticate_path(((r73,True),(r35,True),(r51,True)), 7)
    check(clocks == (5,0) and root_word == (1,1,2,3,4), "actual_receipt_export")
    check(export_from_receipts((r51,r73,r35),7) == root_word, "actual_receipt_export")
    check(export_from_receipts((r51,r73,r35),1) == (), "actual_receipt_export")
    reverse_clocks, reverse_word = authenticate_path(((r73,False),(Receipt(7,1,root_word,()),True)), 3)
    check(reverse_clocks == (2,0) and reverse_word == (1,4), "actual_receipt_export")
    actual_net = Network((1,3,5,7), ((7,3,F(1)),(3,5,F(1)),(5,1,F(1))))
    h = {1:F(0),3:F(0),5:F(0),7:F(1)}
    check(trial_capacity(actual_net,7,1,h)["capacity"] == F(1,3), "actual_receipt_export")

    # Exactly the reported26-edge partial observation, never an unbounded search.
    x, partial_word = 4591, []
    for _ in range(26):
        numerator = 3*x+1
        a = (numerator & -numerator).bit_length()-1
        partial_word.append(a)
        x = numerator >> a
    check(x == 2717873, "retained_unresolved_prefix")
    partial = Receipt(4591, x, tuple(partial_word), ())
    check(verify_receipt(partial) == (26,0), "retained_unresolved_prefix")
    check(export_from_receipts((partial,),4591) is None, "retained_unresolved_prefix")
    pending_net = Network((1,4591,x), ((4591,x,F(1)),))
    h = {1:F(0),4591:F(1),x:F(0)}
    check(trial_capacity(pending_net,4591,1,h)["capacity"] == 0, "retained_unresolved_prefix")

    # Every fixed exact affine ROOT word admits at most one source, even though
    # its valuation cylinder contains arbitrarily many positive integer inputs.
    for depth in range(1, 33):
        u = F(1, depth+1)
        check(u/(1+u) == F(1,depth+2) and u > 0, "slow_drift_not_finite_hit")
        check((1+(1 << depth))-1 == 1 << depth, "profinite_near_ROOT_not_ROOT")
    for t in range(32):
        n = 5+32*t
        check(replay(n,(4,)) == 1+6*t, "cylinder_is_not_fixed_ROOT")

    hostile = [lambda: Network((True,1),()), lambda: Network((0,1),((0,1,1.0),)),
               lambda: Network((0,1),((False,1,F(1)),)), lambda: solve([[F(0)]],[F(1)]),
               lambda: verify_receipt(Receipt(7,9,(1,),(1,))),
               lambda: verify_receipt(Receipt(39,1,(),())),
               lambda: verify_receipt(Receipt(True,1,(),())),
               lambda: verify_receipt(Receipt(5,1,(4,2),())),
               lambda: authenticate_path(((r73,True),(r51,True)),7),
               lambda: then_clock((True,0),(0,0)),
               lambda: eliminate(actual_net,1,(1,7)),
               lambda: actual_net.energy({True:F(0),3:F(0),5:F(0),7:F(1)}),
               lambda: export_from_receipts((r35,r51),True),
               lambda: authenticate_path(([r35,True],(r51,True)),3)]
    for operation in hostile:
        try:
            operation()
        except ValueError:
            check(True, "malformed_or_unpaid_receipt")
        else:
            check(False, "malformed_or_unpaid_receipt")

    print("Finite grounding, endpoint charges and exact residual energy")
    print("Imported six paper headlines accepted; only elementary transfer controls below")
    print(f"Complete simple-graph universe: n2..5, {graph_count} graphs, fixed distinct terminals")
    print("Two rational trials per graph; both forward/reverse interior elimination orders")
    print("Disconnected lengthL chain: pointwise residual1/L, trial energy1/L, exact capacity0")
    print("Macro refinement: unit-resistance expansion has capacity1/L; resistance-preserving expansion retains1")
    print("L parallel L-edge paths: capacity1 while shortest path lengthL; L1..7")
    print(f"Authenticated7--3--5--1: clocks{clocks}, first-hit word{root_word}")
    print(f"Reverse seam3--7--1: clocks{reverse_clocks}, first-hit word{reverse_word}")
    print(f"Frozen partial4591->{x}: word{tuple(partial_word)}, isolatedROOT gives capacity0")
    print("No new ROOT coverage, global positivity, or universal deadline claimed")
    for key, value in sorted(CHECKS.items()):
        print(f"{key}: {value}")
    print(f"TOTAL: {sum(CHECKS.values())}")


if __name__ == "__main__":
    main()
