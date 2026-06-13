"""
Directly compute H(T_11), the number of directed Hamiltonian paths in the
Paley tournament T_11 on Z/11Z (connection set QR = {1,3,4,5,9}).

This is a direct, assumption-free check of CONJ-002 (H(T_11) = 4455).
No OCF formula is assumed — just direct enumeration.

T_11 has 11 vertices. Directed Hamiltonian paths: sequences (v_1,...,v_11)
with all vertices distinct and v_i -> v_{i+1} for all i.

11! = 39,916,800 orderings to check, but with backtracking the search space
is much smaller.
"""

QR = {1, 3, 4, 5, 9}
N = 11
VERTICES = list(range(N))

def has_arc(i, j):
    return (j - i) % N in QR

# Precompute adjacency lists
adj = {v: [] for v in VERTICES}
for v in VERTICES:
    for w in VERTICES:
        if v != w and has_arc(v, w):
            adj[v].append(w)

def count_ham_paths():
    """Count all directed Hamiltonian paths in T_11."""
    count = 0
    visited = [False] * N

    def backtrack(path):
        nonlocal count
        if len(path) == N:
            count += 1
            return
        current = path[-1]
        for nxt in adj[current]:
            if not visited[nxt]:
                visited[nxt] = True
                path.append(nxt)
                backtrack(path)
                path.pop()
                visited[nxt] = False

    # By vertex-transitivity, H(T_11) = 11 * (# Ham paths starting at vertex 0)
    # But let's compute directly to be safe
    for start in VERTICES:
        visited[start] = True
        backtrack([start])
        visited[start] = False

    return count

print("Computing H(T_11) by direct enumeration...")
print("(This may take a moment for an 11-vertex tournament)")
H = count_ham_paths()
print(f"\nH(T_11) = {H}")
print(f"CONJ-002 predicts: H(T_11) = |Aut(T_11)| * 3^4 = 55 * 81 = {55*81}")
if H == 55 * 81:
    print("=> CONJ-002 CONFIRMED for p=11")
else:
    print(f"=> CONJ-002 REFUTED for p=11")
    print(f"   Difference: H(T_11) - 4455 = {H - 4455}")
    if H % 55 == 0:
        print(f"   H(T_11) / 55 = {H // 55} (check: is this 3^k for some k?)")
        k = H // 55
        import math
        if k > 0 and abs(math.log(k) / math.log(3) - round(math.log(k) / math.log(3))) < 1e-9:
            print(f"   Yes! H(T_11) = 55 * 3^{round(math.log(k)/math.log(3))}")
    else:
        print(f"   H(T_11) mod 55 = {H % 55}")

print(f"\nCross-check with OCF lower bound:")
print(f"  c_9(T_11) = 11055 (computed separately)")
print(f"  OCF lower bound: I(Omega,2) >= 1 + 2*(55+594+1320+11055) = {1 + 2*(55+594+1320+11055)}")
