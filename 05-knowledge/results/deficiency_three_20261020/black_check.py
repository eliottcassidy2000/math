"""Exact labelled controls for endpoint-changing R4 black detours.

No removable assertions. Run with python or python -O from the repo root.
The all-order proof and the finite universe are stated in black.md.
"""
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SPEC = spec_from_file_location("old_blue", ROOT / "05-knowledge/results/defect_transitions_20261020/blue_check.py")
B = module_from_spec(SPEC)
SPEC.loader.exec_module(B)
G, need = B.G, B.need
R4 = G.tournament(4, 0, True)
TAILS = [(0, 1, 2, 3), (1, 2, 3, 0), (2, 3, 0, 1),
         (3, 0, 1, 2), (3, 1, 2, 0)]
FORMS = {G.induced(R4, p): i for i, p in enumerate(TAILS)}
ALTERNATE = (1, 0, 0, 0, 1)


def reverse_path(t, p):
    need(B.is_path(t, p), "path reversal requires a directed spanning path")
    for u, v in zip(p, p[1:]):
        t = B.flip(t, u, v)
    return t


def edges(p):
    return {tuple(sorted(e)) for e in zip(p, p[1:])}


def template(left, middle, right, kind):
    """left/right are R4 labels 0,1,2,3; middle is an ordered opaque block."""
    a, b, c, d = right
    q = TAILS[kind]
    prefix = tuple(left[3-v] for v in reversed(q)) + tuple(middle)
    target = prefix + tuple(right[v] for v in q)
    A = prefix
    if kind == 0:
        moves = (A+(b,c,d,a), (a,d,b)+A[::-1]+(c,), (c,)+A+(a,b,d))
    elif kind == 1:
        moves = (A+(a,b,c,d), (d,b,a)+A[::-1]+(c,), (c,)+A+(b,d,a))
    elif kind == 2:
        J = tuple(left)+tuple(middle)
        moves = (J+(c,d,a,b), (a,d,b,c)+J[::-1], A+(c,b,d,a))
    elif kind == 3:
        moves = (A+(a,b,c,d), (d,c,a)+A[::-1]+(b,), (b,)+A+(d,a,c))
    else:
        moves = (A+(b,c,d,a), (b,)+A[::-1]+(a,d,c), (c,a)+A+(d,b))
    return target, moves


def verify(t, target, moves):
    need(B.is_path(t, target), "invalid target path")
    need(edges(target) == edges(moves[0]) ^ edges(moves[1]) ^ edges(moves[2]),
         "edge parity ledger fails")
    expected = reverse_path(t, target)
    colours, endpoint_sums = [], []
    for p in moves:
        need(B.is_path(t, p), "macro has invalid directed path")
        colours.append(B.blue(t, p))
        endpoint_sums.append(t[p[0]].bit_count()+t[p[-1]].bit_count())
        need(not colours[-1], "macro contains a blue move")
        t = reverse_path(t, p)
    need(t == expected, "detour does not reproduce the exact labelled target")
    return endpoint_sums


def endcap_audit():
    rows = []
    for m in range(7):
        count = 0
        for code in range(1 << len(G.free_pairs(m))):
            z = G.tournament(m, code, True)
            t = B.ordered_join([R4,z,R4])
            for kind in range(5):
                p, moves = template(tuple(range(4)), tuple(range(4,m+4)),
                                    tuple(range(m+4,m+8)), kind)
                verify(t, p, moves)
                count += 1
        rows.append((m, 1 << len(G.free_pairs(m)), count))
    return rows


def ear_data(t):
    """Standard blue path; return the guarded R4 terminal-ear macro if present."""
    n = len(t)
    if n < 8:
        return None
    tail = G.induced(t, range(n-4,n))
    if tail not in FORMS:
        return None
    kind = FORMS[tail]
    right = [None]*4
    for j,v in enumerate(TAILS[kind]):
        right[v] = n-4+j
    if not G.arc(t,n-5,right[ALTERNATE[kind]]):
        return None
    left = tuple(n-1-right[3-i] for i in range(4))
    return template(left, tuple(range(4,n-4)), right, kind)


def ear_audit():
    rows, hostile = [], None
    for size in range(4,7):
        total = tail_count = guarded = certificates = 0
        kind_counts = Counter()
        for code in range(1 << len(G.free_pairs(size))):
            sink = G.tournament(size,code,True)
            if not B.strong(sink):
                continue
            total += 1
            tail = G.induced(sink, range(size-4,size))
            tail_count += tail in FORMS
            for middle_size in range(3):
                middle = G.tournament(middle_size,0,True)
                t = B.ordered_join([G.transpose(sink),middle,sink])
                need(B.blue(t,tuple(range(len(t)))), "ear source must be blue")
                result = ear_data(t)
                if result is not None:
                    verify(t,*result)
                    certificates += 1
                    if middle_size == 0:
                        guarded += 1
                        kind_counts[FORMS[tail]] += 1
                elif size == 5 and code == 0 and middle_size == 0:
                    hostile = t
        rows.append((size,total,tail_count,guarded,certificates,
                     tuple(kind_counts[k] for k in range(5))))
    need(hostile is not None, "missing boundary hostile")
    n = len(hostile)
    right = tuple(range(n-4,n))
    left = tuple(n-1-right[3-i] for i in range(4))
    p,moves = template(left,tuple(range(4,n-4)),right,0)
    need(B.is_path(hostile,p) and B.blue(hostile,p), "hostile target changed")
    need(not B.is_path(hostile,moves[0]), "missing boundary arc must invalidate P1")
    return rows, hostile, moves[0]


def complete_small_source_audit():
    rows=[]
    for n in (8,9):
        total = eligible = 0
        for t in B.nonstrong_blue_sources(n):
            total += 1
            result = ear_data(t)
            if result:
                verify(t,*result)
                eligible += 1
        rows.append((n,total,eligible))
    return rows


def main():
    need(len(B.hamilton_paths(R4)) == 5, "R4 path inventory")
    need(set(B.hamilton_paths(R4)) == set(TAILS), "R4 path labels")
    print("PROVED-template exact controls; no connectivity conclusion")
    print("Endcap controls (middle order, all fixed-path presentations, three-black certificates):")
    for row in endcap_audit():
        print(row)
    rows,hostile,bad = ear_audit()
    print("Terminal-ear controls (SCC order, strong presentations, tail R4, guarded, certificates, five kinds):")
    for row in rows:
        print(row)
    print("Complete nonstrong blue source controls (n, sources, eligible):")
    for row in complete_small_source_audit():
        print(row)
    print("Guard hostile rows:", hostile)
    print("Guard hostile invalid first move:", bad)
    print("PASS: all exact matrices, path directions, edge XORs and three black predicates")


if __name__ == "__main__":
    main()
