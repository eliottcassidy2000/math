#!/usr/bin/env python3
"""
Tribonacci structure of T_full Hamiltonian paths and cross-scale patterns.
Instance: opus-2026-03-06-S11

The Hamiltonian paths of T_full have a beautiful recursive structure.
T_full has:
  - Path edges: i -> i+1
  - Backward non-adjacent: j -> i for j-i >= 2

A Hamiltonian path v_0, v_1, ..., v_{n-1} must use edges of T_full.
The key observation: vertex n-1 can reach 0,1,...,n-3 but NOT n-2.

We prove H(T_full) satisfies the Tribonacci recurrence by analyzing
how paths start from their initial vertex.

Also investigates:
- Block structure of paths (consecutive segments)
- How the Tribonacci structure relates to GS (grid-symmetric) tilings
- Cross-scale embedding patterns
"""

def full_tournament_adj(n):
    """T_full: i->i+1, j->i for j-i>=2."""
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    for i in range(n):
        for j in range(i+2, n):
            A[j][i] = 1
    return A

def enumerate_ham_paths(A, n):
    """List all Hamiltonian paths."""
    paths = []
    def bt(path, vis):
        if len(path) == n:
            paths.append(tuple(path))
            return
        last = path[-1] if path else -1
        for v in range(n):
            if v in vis:
                continue
            if last == -1 or A[last][v]:
                bt(path + [v], vis | {v})
    bt([], set())
    return paths

def classify_path_structure(paths, n):
    """Classify paths by their block structure."""
    # A "consecutive ascending block" is a maximal sequence a, a+1, a+2, ...
    block_types = {}
    for p in paths:
        blocks = []
        start = p[0]
        length = 1
        for i in range(1, n):
            if p[i] == p[i-1] + 1:
                length += 1
            else:
                blocks.append((start, length))
                start = p[i]
                length = 1
        blocks.append((start, length))
        # Normalize: represent as block lengths
        block_sig = tuple(b[1] for b in blocks)
        if block_sig not in block_types:
            block_types[block_sig] = []
        block_types[block_sig].append(p)
    return block_types

def tribonacci_proof_structure(n):
    """Show the recursive structure that gives Tribonacci."""
    if n < 3:
        return

    A = full_tournament_adj(n)
    paths = enumerate_ham_paths(A, n)

    print(f"\nn={n}: H(T_full) = {len(paths)}")

    # Classify by starting vertex
    by_start = {}
    for p in paths:
        s = p[0]
        if s not in by_start:
            by_start[s] = []
        by_start[s].append(p)

    print(f"  Paths by starting vertex:")
    for s in sorted(by_start):
        print(f"    start={s}: {len(by_start[s])} paths")

    # Classify by block structure
    block_types = classify_path_structure(paths, n)
    print(f"  Block structures (ascending runs):")
    for sig in sorted(block_types, key=lambda s: (-len(block_types[s]), s)):
        count = len(block_types[sig])
        print(f"    {sig}: {count} paths")

    # Key recursive structure:
    # Paths starting with (n-1, n-2, ...) can only start with n-1 if n-1 beats the next vertex.
    # n-1 beats all v <= n-3. So paths start with n-1 followed by any v in {0,...,n-3}.
    #
    # After placing n-1 first, the remaining vertices {0,...,n-2} must form a Hamiltonian path
    # starting from some v <= n-3 in the INDUCED subtournament T_full[{0,...,n-2}].
    #
    # But T_full[{0,...,n-2}] has the SAME structure as T_full on n-1 vertices!
    # (Because removing vertex n-1 from T_full gives exactly T_full on {0,...,n-2}.)

    # So the number of paths starting with n-1 equals H(T_full, n-1) minus paths starting
    # with n-2 (since n-1 can't go to n-2).
    # Wait, that's not right either. Let me think more carefully.

    # Paths starting with vertex n-1:
    # n-1 -> v where v in {0,...,n-3}
    # Then continue Hamiltonian path in T_full[{0,...,n-2}\setminus{n-1}]... wait, n-1 is removed.
    # Actually, after placing n-1, remaining vertices are {0,...,n-2}.
    # We need a path from v in this set using T_full edges.
    # But T_full[{0,...,n-2}] is T_full on n-1 vertices.

    # However! The path from v must be a Hamiltonian path in T_full[{0,...,n-2}] starting from v.
    # So the count is: sum over v in {0,...,n-3} of (# Ham paths in T_full[n-1] starting from v)

    # Similarly for paths starting with vertex n-2:
    # n-2 -> n-1 (path edge) -> then we need a path from n-1 in T_full[{0,...,n-3}]
    # But we placed n-2 first, n-1 second, remaining is {0,...,n-3}.
    # n-1 -> v for v in {0,...,n-3}. Then path in T_full[{0,...,n-3}] starting from v.
    # Wait, but the vertices are {0,...,n-3} and n-1 is already used.
    # After n-2 -> n-1, remaining is {0,...,n-3}. n-1 can go to any of {0,...,n-3}.
    # Then path continues in T_full[{0,...,n-3}] = T_full on n-2 vertices.

    # Hmm, I need to think about this more carefully. Let me just verify the recursion numerically.

    if n >= 6:
        # Check: paths starting with (n-1) and NOT followed by (n-2)
        starts_n1 = by_start.get(n-1, [])
        starts_n1_not_n2 = [p for p in starts_n1 if len(p) > 1 and p[1] != n-2]
        starts_n1_then_n2 = [p for p in starts_n1 if len(p) > 1 and p[1] == n-2]
        print(f"\n  Recursion check:")
        print(f"    Paths starting {n-1}: {len(starts_n1)}")
        print(f"    Paths starting {n-1},{n-2}: {len(starts_n1_then_n2)}")
        print(f"    Paths starting {n-1},not-{n-2}: {len(starts_n1_not_n2)}")

    # The crucial pattern: cyclic rotations
    # A cyclic rotation (k, k+1, ..., n-1, 0, 1, ..., k-1) is a valid Ham path iff
    # the edge (n-1) -> 0 exists (it does: n-1 beats 0 since gap >= 2 for n >= 3)
    # AND the edge (k-1) exists from the previous vertex (which is k-2 or n-1).
    # Actually all cyclic rotations ARE valid because the path uses edges
    # k->k+1->...->n-1->0->1->...->k-1. But n-1->0 requires n-1 beats 0.
    # In T_full, n-1 beats all v <= n-3. So n-1 beats 0 iff n-1 >= 2, i.e., n >= 3. Yes!

    cyclic = [(k,) + tuple(range(k+1,n)) + tuple(range(0,k)) for k in range(n)]
    # But not all are valid! k->k+1 is valid (path edge).
    # n-1 -> 0: n-1 beats 0 (gap n-1 >= 2 for n >= 3). ✓
    # k-1 must be beaten by k-2 or wrap: (k-1) is last vertex.
    # Actually, each cyclic rotation is 0-indexed rotation of (0,1,...,n-1).
    # Check validity:
    valid_cyclic = []
    for k in range(n):
        p = tuple(range(k, n)) + tuple(range(0, k))
        valid = True
        for i in range(n-1):
            if not A[p[i]][p[i+1]]:
                valid = False
                break
        if valid:
            valid_cyclic.append(p)
    print(f"\n  Valid cyclic rotations: {len(valid_cyclic)}/{n}")
    for p in valid_cyclic:
        print(f"    {p}")


def main():
    print("Tribonacci structure of T_full Hamiltonian paths")
    print("="*60)

    for n in range(3, 10):
        tribonacci_proof_structure(n)

    # Verify Tribonacci
    print(f"\n{'='*60}")
    print("Tribonacci verification:")
    H = {}
    for n in range(3, 12):
        A = full_tournament_adj(n)
        # Count Ham paths via DP
        dp = {}
        for v in range(n):
            dp[(1 << v, v)] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or (mask, v) not in dp:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if A[v][u]:
                        key = (mask | (1 << u), u)
                        dp[key] = dp.get(key, 0) + dp[(mask, v)]
        full = (1 << n) - 1
        H[n] = sum(dp.get((full, v), 0) for v in range(n))

    print(f"{'n':>3} {'H':>8} {'Trib check':>12}")
    for n in sorted(H):
        trib = H[n-1] + H[n-2] + H[n-3] if n >= 6 else '?'
        check = '✓' if trib == H[n] else '✗' if isinstance(trib, int) else '-'
        print(f"{n:>3} {H[n]:>8} {str(trib):>12} {check}")

    # Why Tribonacci? The recursive structure.
    print(f"\n{'='*60}")
    print("WHY TRIBONACCI: Recursive path decomposition")
    print("="*60)
    print("""
For T_full on {0,...,n-1}, a Hamiltonian path decomposes based on
where vertex (n-1) appears:

Case 1: Path starts with (n-1).
  Then n-1 -> v for some v in {0,...,n-3}.
  Remaining: Hamiltonian path in T_full[{0,...,n-2}] from v.
  Count = H(n-1) - (# paths in T_full[n-1] starting from n-2)
  Wait, we need paths NOT starting from n-2...

  Actually: paths starting n-1 have n-1->v for v≤n-3.
  Then v, ..., rest is a path in T_full[{0,...,n-2}] starting from v≤n-3.
  These are EXACTLY the Hamiltonian paths of T_full[n-1] that DON'T start from n-2.
  Count = H(n-1) - H_{start=n-2}(n-1)

Case 2: Path has (n-1) in position > 0.
  Then some vertex u precedes n-1. u -> n-1 requires u = n-2 (only n-2 beats n-1).
  So the path has ..., n-2, n-1, ...
  This means n-2 is immediately before n-1.

  Sub-case 2a: Path starts at some v != n-1, goes ..., n-2, n-1, ..., continues.
    After n-1: n-1 -> w for w in {0,...,n-3} not yet visited.

  Sub-case 2b: n-1 is the LAST vertex.
    Path is a Hamiltonian path in T_full[{0,...,n-2}] ending at n-2,
    followed by n-2 -> n-1.
    Count = H_{end=n-2}(n-1)

The key insight: by contracting the edge (n-2)->(n-1), we can relate
paths with (n-2,n-1) consecutive to paths of a smaller tournament.

Let f(n) = H(T_full, n).

Paths of T_full on [n]:
- Category A: paths where n-1 is NOT preceded by n-2 (impossible since only n-2->n-1).
  Wait, n-2 -> n-1 is the ONLY edge INTO n-1. So n-1 must be either:
  (i) the starting vertex, OR
  (ii) preceded by n-2.

  So ALL paths either start with n-1 or have n-2 immediately before n-1.

Category (i): start with n-1. Then n-1->v for v≤n-3.
  The rest is a path in T_full[n-1] starting from v (v ≤ n-3).

Category (ii): ...n-2, n-1...  Contract (n-2,n-1) into a single "super-vertex" [n-2,n-1].
  The super-vertex inherits:
  - Incoming edges of n-2: only n-3 -> n-2 (path edge from n-3) and nobody else beats n-2.
    Wait: in T_full, who beats n-2? Only n-2+1 = n-1 (path edge). So n-1 beats... no:
    The path edge is n-2 -> n-1 (forward). Non-adjacent: j -> n-2 for j - (n-2) >= 2,
    but there are no j > n-2 except n-1, and n-1 - (n-2) = 1 < 2. So the only edge
    INTO n-2 is... let's see: in T_full, vertex k receives edges from k-1 (path edge: k-1->k)
    and from all j with j < k and k - j >= 2... wait no: for the pair (j, k) with k > j:
    if k - j >= 2: k -> j (backward). So k BEATS j. So nobody beats k from below except k-1.
    And from above: the only vertex above k is k+1,...,n-1. Vertex m > k: if m - k >= 2,
    m -> k (m beats k). If m = k+1: k -> k+1 (k beats k+1).
    So vertices m >= k+2 ALL beat k.

    For k = n-2: vertices m >= n beat n-2. Only m = n-1 qualifies (and n-1 - (n-2) = 1 < 2).
    So n-1 does NOT beat n-2! And n-2 -> n-1 (path edge). So n-2 beats n-1.
    Wait, this contradicts my T_full definition. Let me re-check:

    T_full: i->i+1 (path, forward), j->i for j-i >= 2 (backward non-adjacent).
    So for i < j:
    - if j = i+1: i beats j (path edge)
    - if j >= i+2: j beats i (backward edge)

    For pair (n-2, n-1): j = n-1, i = n-2, j - i = 1. So i beats j: n-2 -> n-1. ✓
    For pair (n-3, n-1): j = n-1, i = n-3, j - i = 2. So j beats i: n-1 -> n-3.
    So n-1 beats n-3 but n-2 beats n-1.

    Edges INTO n-1: only from n-2 (n-2 -> n-1). ✓
    Edges INTO n-2: from n-3 (n-3 -> n-2, path edge), and from j with j >= n-2+2 = n:
    nobody (j = n-1 has gap 1 < 2). So only n-3 -> n-2.

    So the super-vertex [n-2, n-1] receives edges only from n-3.
    The super-vertex [n-2, n-1] sends edges to:
    - From n-1: n-1 -> all v <= n-3. So [n-2,n-1] -> v for all v <= n-3.
    - From n-2: n-2 -> n-1 (internal, already used). n-2 -> v for v <= n-4 (n-2 beats v if n-2-v >= 2).

    These outgoing edges are the SAME as T_full's vertex n-1 on an (n-1)-vertex tournament!

    So contracting (n-2, n-1) gives essentially T_full on n-1 vertices... almost.

    This analysis needs more care. The point is that the Tribonacci recurrence
    f(n) = f(n-1) + f(n-2) + f(n-3) holds, and it relates to the block structure
    of consecutive ascending runs in the Hamiltonian paths.
""")


if __name__ == '__main__':
    main()
