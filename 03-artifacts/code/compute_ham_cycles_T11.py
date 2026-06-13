"""
Compute h_QR and h_NQR for T_11 (Paley tournament on Z/11Z).

T_11: arcs i -> j iff (j-i) mod 11 in QR = {1,3,4,5,9}

h_QR = # directed Hamiltonian cycles in T_11 \ {0, 1}  (vertices: {2,3,4,5,6,7,8,9,10})
h_NQR = # directed Hamiltonian cycles in T_11 \ {0, 2}  (vertices: {1,3,4,5,6,7,8,9,10})

A directed Hamiltonian cycle visits each vertex exactly once and returns to start.
We count by fixing the smallest vertex as the start, enumerating Hamiltonian paths
from that start back to start. Each cycle is counted exactly once this way.

This settles OPEN-Q-009 and provides evidence for/against CONJ-002.
"""

QR = {1, 3, 4, 5, 9}

def has_arc(i, j, n=11):
    """Does T_11 have arc i -> j?"""
    return (j - i) % n in QR

def count_ham_cycles(vertices):
    """
    Count directed Hamiltonian cycles in the sub-tournament induced by `vertices`.
    Strategy: fix the minimum vertex as the canonical start.
    Enumerate all Hamiltonian paths starting from it, check if last vertex
    has an arc back to start.
    """
    vlist = sorted(vertices)
    n = len(vlist)
    start = vlist[0]
    others = vlist[1:]

    # Precompute adjacency
    adj = {}
    for v in vlist:
        adj[v] = set()
        for w in vlist:
            if v != w and has_arc(v, w):
                adj[v].add(w)

    count = 0

    # Backtracking: path is current partial Hamiltonian path starting from `start`
    def backtrack(path, visited):
        nonlocal count
        if len(path) == n:
            # Check if there's an arc from last vertex back to start
            if start in adj[path[-1]]:
                count += 1
            return
        current = path[-1]
        for nxt in adj[current]:
            if nxt not in visited:
                visited.add(nxt)
                path.append(nxt)
                backtrack(path, visited)
                path.pop()
                visited.remove(nxt)

    backtrack([start], {start})
    return count

# T_11 \ {0, 1}: vertices {2, 3, 4, 5, 6, 7, 8, 9, 10}
vertices_QR = list(range(2, 11))
h_QR = count_ham_cycles(vertices_QR)

# T_11 \ {0, 2}: vertices {1, 3, 4, 5, 6, 7, 8, 9, 10}
vertices_NQR = [1] + list(range(3, 11))
h_NQR = count_ham_cycles(vertices_NQR)

print(f"h_QR  = h({{0,1}}) = {h_QR}")
print(f"h_NQR = h({{0,2}}) = {h_NQR}")
print(f"h_QR + h_NQR = {h_QR + h_NQR}")
print(f"c_9(T_11) = (55/2) * (h_QR + h_NQR) = {55 * (h_QR + h_NQR) / 2}")
print()
print(f"CONJ-002 requires h_QR + h_NQR <= 8 for H(T_11) = 4455 to hold.")
if h_QR + h_NQR <= 8:
    c9 = 55 * (h_QR + h_NQR) // 2
    print(f"  -> CONSTRAINT SATISFIED: h_QR + h_NQR = {h_QR + h_NQR} <= 8")
    print(f"  -> c_9 = {c9}, consistent with CONJ-002")
else:
    c9 = 55 * (h_QR + h_NQR) / 2
    print(f"  -> CONSTRAINT VIOLATED: h_QR + h_NQR = {h_QR + h_NQR} > 8")
    print(f"  -> c_9 = {c9}, CONJ-002 would be FALSE for p=11")
