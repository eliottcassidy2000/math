"""aud_lineage_core.py - is the S76 pair's B already in the repository?

zero_ready_component_returns_20261103/component/COMPONENT.md (2026-11-03) records the order-14 named source of the
H = 170807 target, names [1,2,4,6,7,8,14,21,22,31,32,33,40,54], rows [11134,...,894], and its canonical successive
unique-start cores; the order-7 core at mark 4 is listed with (h, e, P, F, S) = (31, 9, 644, 157, 31).
This script rebuilds that order-7 rooted core (mark 4, core {1,21,32,33,40}, tail {54}) from the stored rows, relabels it
along its unique Hamiltonian path from the mark, and compares it with THM-4872's B (and A) in the same normal form.
It also evaluates P = P(2) (all ordered two-slot covers) and S = S_v(2) (second slot nonempty and starting at v) for
A and B, to identify the packet's columns.
"""
import itertools

names = [1, 2, 4, 6, 7, 8, 14, 21, 22, 31, 32, 33, 40, 54]
rows = [11134, 632, 2682, 112, 512, 80, 528, 9087, 126, 296, 13311, 10234, 11263, 894]
n14 = len(names)
pos = {v: i for i, v in enumerate(names)}


def arc(a, b):
    return (rows[pos[a]] >> pos[b]) & 1


# the source must be a tournament
for a, b in itertools.combinations(names, 2):
    assert arc(a, b) + arc(b, a) == 1, (a, b)


def ham_paths(V, adj, start=None):
    """all Hamiltonian paths of the tournament on vertex list V (adj(a,b) = 1 iff a -> b)"""
    out = []

    def rec(path, rest):
        if not rest:
            out.append(tuple(path))
            return
        for y in list(rest):
            if adj(path[-1], y):
                rest.remove(y)
                path.append(y)
                rec(path, rest)
                path.pop()
                rest.add(y)

    for s in V:
        if start is not None and s != start:
            continue
        rec([s], set(V) - {s})
    return out


def count_paths(V, adj):
    if not V:
        return 1
    return len(ham_paths(V, adj))


def profile(V, adj, r):
    Vs = list(V)
    h = count_paths(Vs, adj)
    s = len(ham_paths(Vs, adj, start=r))
    q = count_paths([x for x in Vs if x != r], adj)
    e = sum(1 for p in ham_paths(Vs, adj) if p[-1] == r)
    F = 0
    P2 = 0
    S2 = 0
    others = [x for x in Vs if x != r]
    # ordered two-slot covers (X, Y): X, Y vertex-disjoint paths covering V, either may be empty
    for k in range(0, len(Vs) + 1):
        for X in itertools.combinations(Vs, k):
            Y = [x for x in Vs if x not in X]
            px = ham_paths(list(X), adj) if X else [()]
            py = ham_paths(Y, adj) if Y else [()]
            P2 += len(px) * len(py)
            F += sum(1 for p in px if p and p[-1] == r) * len(py)
            S2 += len(px) * sum(1 for p in py if p and p[0] == r)
    return dict(h=h, s=s, q=q, e=e, F=F, P2=P2, S2=S2)


def normal_rows(V, adj, r):
    ps = ham_paths(list(V), adj, start=r)
    assert len(ps) == 1, len(ps)
    p = ps[0]
    return [sum(1 << j for j in range(len(p)) if j != i and adj(p[i], p[j])) for i in range(len(p))]


core = [4, 1, 21, 32, 33, 40, 54]
print("170807 lineage, order-7 core at mark 4:", core)
pr = profile(core, arc, 4)
print("  profile:", pr)
print("  path-order rows:", normal_rows(core, arc, 4))


def from_rows(R):
    return lambda a, b: (R[a] >> b) & 1


A = [8, 73, 3, 4, 15, 31, 61]
B = [2, 84, 1, 7, 13, 31, 61]
for nm, R in (("A", A), ("B", B)):
    adj = from_rows(R)
    for a, b in itertools.combinations(range(7), 2):
        assert adj(a, b) + adj(b, a) == 1
    print(nm, "profile:", profile(range(7), adj, 0), "path-order rows:", normal_rows(range(7), adj, 0))
