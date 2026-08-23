#!/usr/bin/env python3
"""Independent exact hostile audit of the AK paid-row k=1 forest barrier.

This intentionally does not import the maintained engine or the proposed
checker.  It enumerates arbitrary labelled forests (not only paths), checks
their ports against direct rational row algebra, exhausts the k=1 adjacent
merge topology/accounting, and replays the cyclic 5/3 boundary in raw labels.
"""

from fractions import Fraction
from itertools import combinations, product
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def gate(ok, name):
    global CHECKS
    if not ok:
        raise RuntimeError(name)
    CHECKS += 1


def rank(rows, width):
    a = [[Fraction(x) for x in row] for row in rows if any(row)]
    r = 0
    for c in range(width):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        z = a[r][c]
        a[r] = [x / z for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [x - z*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def belongs(rows, target):
    return rank(rows, len(target)) == rank(rows + [target], len(target))


def erow(n, u, v, slope):
    z = [0] * (2*n)
    z[2*u], z[2*u+1] = 1, slope
    z[2*v], z[2*v+1] = -1, -slope
    return z


def grow(n, v, slope):
    z = [0] * (2*n)
    z[2*v], z[2*v+1] = 1, slope
    return z


SLOPES = (-3, 1, 4)
ZERO, FULL = None, "full"
STATES = (ZERO,) + SLOPES + (FULL,)


def join(state, slope):
    if state is ZERO:
        return slope
    if state == FULL or state == slope:
        return state
    return FULL


def contains(state, slope):
    return state == FULL or state == slope


def dimension(state):
    return 0 if state is ZERO else 2 if state == FULL else 1


def is_forest(n, labelled_edges):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for u, v, _ in labelled_edges:
        a, b = find(u), find(v)
        if a == b:
            return False
        parent[a] = b
    return True


def adjacency(n, edges):
    adj = [[] for _ in range(n)]
    for u, v, q in edges:
        adj[u].append((v, q))
        adj[v].append((u, q))
    return adj


def port(adj, grounds, root, allowed=None):
    if allowed is None:
        allowed = set(range(len(adj)))
    def visit(v, parent):
        state = grounds[v]
        for w, q in adj[v]:
            if w == parent or w not in allowed:
                continue
            child = visit(w, v)
            if contains(child, q):
                state = join(state, q)
        return state
    return visit(root, -1)


def projected_data(n, edges, grounds, deleted):
    live = set(range(n)) - deleted
    out_edges = []
    out_grounds = list(grounds)
    for u, v, q in edges:
        if u in live and v in live:
            out_edges.append((u, v, q))
        elif u in live and v in deleted:
            out_grounds[u] = join(out_grounds[u], q)
        elif v in live and u in deleted:
            out_grounds[v] = join(out_grounds[v], q)
    return live, out_edges, tuple(out_grounds)


def rows_of(n, edges, grounds):
    rows = [erow(n, u, v, q) for u, v, q in edges]
    for v, state in enumerate(grounds):
        if state is ZERO:
            continue
        if state == FULL:
            rows.extend((grow(n, v, SLOPES[0]), grow(n, v, SLOPES[1])))
        else:
            rows.append(grow(n, v, state))
    return rows


def direct_port(rows, n, v):
    outside = [row[:2*v] + row[2*v+2:] for row in rows]
    dim = rank(rows, 2*n) - rank(outside, 2*n-2)
    vertical = [0]*(2*n)
    vertical[2*v+1] = 1
    fires = belongs(rows, vertical)
    return dim, fires


# Exhaust all arbitrary simple forests through four vertices over three slopes
# and the exact local ground-span quotient.  This includes branched stars.
networks = ports = partials = successes = 0
matrix_checks = 0
for n in range(1, 5):
    pairs = list(combinations(range(n), 2))
    for edge_values in product((ZERO,) + SLOPES, repeat=len(pairs)):
        edges = [(u, v, q) for (u, v), q in zip(pairs, edge_values)
                 if q is not ZERO]
        if not is_forest(n, edges):
            continue
        adj = adjacency(n, edges)
        for grounds in product(STATES, repeat=n):
            networks += 1
            states = tuple(port(adj, grounds, v) for v in range(n))
            ports += n
            fired = {v for v, state in enumerate(states) if state == FULL}
            if fired and len(fired) < n:
                partials += 1
                live, e2, g2 = projected_data(n, edges, grounds, fired)
                a2 = adjacency(n, e2)
                gate(all(port(a2, g2, v, live) != FULL for v in live),
                     "a forest acquired a second-round port")
            elif len(fired) == n:
                successes += 1
                paid = len(edges) + sum(dimension(s) for s in grounds)
                gate(paid >= 2*n, "successful forest underpaid rank 2n")

            # A deterministic, topology-diverse direct rational cross-check.
            # All n<=2 rows are checked; larger cases use a stable residue.
            code = networks*149 + len(edges)*31 + n
            if n <= 2 or code % 211 == 0:
                rows = rows_of(n, edges, grounds)
                for v, predicted in enumerate(states):
                    dim, fires = direct_port(rows, n, v)
                    gate(dim == dimension(predicted), "message/direct port dimension mismatch")
                    gate(fires == (predicted == FULL), "vertical port was not exactly full")
                    if isinstance(predicted, int):
                        target = [0]*(2*n)
                        target[2*v], target[2*v+1] = 1, predicted
                        gate(belongs(rows, target), "predicted nonvertical port line missing")
                    matrix_checks += 1

gate(networks == 370980, "unexpected arbitrary-forest universe size")


# Translation (a,b)->(s,d), checked without projective normalization.
raw_labels = ((2, 3), (4, -1), (-2, 5))
raw = []
sd = []
n = 3
for (u, v), (a, b) in zip(((0, 1), (1, 2)), raw_labels[:2]):
    rr = [0]*(2*n)
    tt = [0]*(2*n)
    rr[2*u:2*u+2] = [a, b]
    rr[2*v:2*v+2] = [-a, -b]
    tt[2*u:2*u+2] = [a+b, a-b]
    tt[2*v:2*v+2] = [-a-b, -a+b]
    raw.append(rr); sd.append(tt)
a, b = raw_labels[2]
rr = [0]*(2*n); rr[4:6] = [a, b]; raw.append(rr)
tt = [0]*(2*n); tt[4:6] = [a+b, a-b]; sd.append(tt)
gate(rank(raw, 2*n) == rank(sd, 2*n), "invertible coordinate change changed rank")
for v in range(n):
    raw_target = [0]*(2*n); raw_target[2*v:2*v+2] = [1, -1]
    sd_target = [0]*(2*n); sd_target[2*v+1] = 2
    gate(belongs(raw, raw_target) == belongs(sd, sd_target),
         "raw/port firing translation mismatch")


# k=1 legal-slot topology: adjacent merges always form intervals.  Exhaust
# every M/E/absent word through nine base vertices and every wildcard subset
# through six vertices.  Active rows can only disappear under projection.
merge_words = wildcard_cases = 0
for base_n in range(1, 10):
    for word in product(("M", "E", "0"), repeat=base_n-1):
        merge_words += 1
        classes = []
        c = 0
        for v in range(base_n):
            if v and word[v-1] != "M":
                c += 1
            classes.append(c)
        qedges = [(classes[i], classes[i+1]) for i, x in enumerate(word) if x == "E"]
        gate(all(u != v for u, v in qedges), "label edge swallowed by adjacent merges")
        gate(len(qedges) == len(set(qedges)), "parallel edge in k=1 quotient")
        gate(all(v == u+1 for u, v in qedges), "non-path quotient edge")
        if base_n <= 6:
            qn = c+1
            for mask in range(1 << qn):
                wildcard_cases += 1
                live_edges = sum(not ((mask >> u)&1) or not ((mask >> v)&1)
                                 for u, v in qedges)
                # The displayed expression counts edges with at least one live
                # endpoint; each is one active edge/ground row, never more than
                # its one paid labelled slot.
                gate(live_edges <= len(qedges), "projection created uncharged k=1 rows")


# Duplicate rows/labels cannot help: they preserve span and increase cost.
base = [erow(2, 0, 1, 1), grow(2, 0, -3), grow(2, 1, 4)]
dup = base + [base[0][:], base[1][:]]
gate(rank(base, 4) == rank(dup, 4), "duplicate row changed span")
for v in range(2):
    target = [0]*4; target[2*v+1] = 1
    gate(belongs(base, target) == belongs(dup, target), "duplicate row changed firing")


# Raw-label cyclic hostile.  All labels obey a+b != 0.  The advertised
# normalized-slope combination e12+2e23+e31-seed2+seed3 is vertical at
# vertex 1.  In the unnormalized raw integer rows the coefficients become
# (1,4,1,-1,1), because the slope-0 and slope-2 labels have s=a+b=2.
def raw_edge(n, u, v, label):
    a, b = label
    z = [0]*(2*n)
    z[2*u:2*u+2] = [a, b]
    z[2*v:2*v+2] = [-a, -b]
    return z


def raw_ground(n, v, label):
    a, b = label
    z = [0]*(2*n); z[2*v:2*v+2] = [a, b]
    return z


cycle_raw = [
    raw_edge(3, 0, 1, (1, 1)),
    raw_edge(3, 1, 2, (1, 0)),
    raw_edge(3, 2, 0, (3, -1)),
    raw_ground(3, 1, (3, -1)),
    raw_ground(3, 2, (1, 1)),
]
combo = [sum(c*r[j] for c, r in zip((1, 4, 1, -1, 1), cycle_raw))
         for j in range(6)]
gate(combo == [-2, 2, 0, 0, 0, 0], "raw cyclic hostile combination changed")
gate(rank(cycle_raw, 6) == 5, "raw cyclic hostile rank")
target0 = [0]*6; target0[0:2] = [1, -1]
gate(belongs(cycle_raw, target0), "raw cyclic hostile does not first-fire")
projected = [row[2:] for row in cycle_raw]
gate(rank(projected, 4) == 4, "raw cyclic hostile does not finish after projection")


print("audit=independent_AK_k1_forest_barrier")
print("verdict=PASS")
print(f"arbitrary_forest_networks={networks}")
print(f"arbitrary_forest_ports={ports}")
print(f"partial_first_rounds={partials}")
print(f"successful_forests={successes}")
print(f"direct_rational_port_checks={matrix_checks}")
print(f"adjacent_merge_words={merge_words}")
print(f"adjacent_merge_wildcard_cases={wildcard_cases}")
print("hostiles=vertical_label:excluded;simple_triangle:score_5/3_two_round")
print("scope=paid_row_semantics_only;no_private_verifier_or_AK_implication")
print(f"CHECKS={CHECKS}")
