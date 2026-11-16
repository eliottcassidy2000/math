"""Collatz tree (rooted at 1; children of m: 2m, and (m-1)/3 when m = 4 mod 6, m > 4) built one node at a time.
Orders: BFS with even-first / odd-first ties, DFS, and numeric order of labels restricted to a fixed BFS prefix.
Each order gives a tree-reversal tournament (new node beats only its parent).  H via ./treeH."""
import subprocess, sys
def children(m):
    c = [2 * m]
    if m % 6 == 4 and m > 4: c.append((m - 1) // 3)
    return c
def bfs(N, pref):
    order = [1]; par = {1: None}; frontier = [1]
    while len(order) < N:
        nxt = []
        for m in frontier:
            for c in children(m):
                par[c] = m; nxt.append(c)
        if pref == 'even': nxt.sort(key=lambda c: (c % 2, c))
        elif pref == 'odd': nxt.sort(key=lambda c: (1 - c % 2, c))
        elif pref == 'num': nxt.sort()
        for c in nxt:
            if len(order) < N: order.append(c)
        frontier = nxt
    return order, par
def dfs(N):
    order = []; par = {1: None}; stack = [1]
    while stack and len(order) < N:
        m = stack.pop(); order.append(m)
        for c in sorted(children(m), reverse=True):   # visit smaller child first
            par[c] = m; stack.append(c)
    return order, par
def H_of_orders(orders):
    lines = []
    for order, par in orders:
        idx = {v: i for i, v in enumerate(order)}
        ps = [idx[par[v]] for v in order[1:]]
        lines.append(f"{len(order)} " + " ".join(map(str, ps)))
    out = subprocess.run(['./w46_treeH'], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout.split()
    return [int(x) for x in out]
if __name__ == '__main__':
    NM = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    for pref in ('even', 'odd'):
        orders = [bfs(N, pref) for N in range(1, NM + 1)]
        print(f"BFS {pref}-first labels:", orders[-1][0])
        print("   H:", H_of_orders(orders))
    orders = [dfs(N) for N in range(1, NM + 1)]
    print("DFS labels:", orders[-1][0]); print("   H:", H_of_orders(orders))
    chain = [(list(range(N)), {i: (i - 1 if i else None) for i in range(N)}) for N in range(1, NM + 1)]
    print("single doubling row (chain):", H_of_orders(chain))
