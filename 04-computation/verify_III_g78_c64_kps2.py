"""
kind-pasteur-2026-06-09-S2 : adversarial follow-up. Resolve the G78-class-1 C64
disagreement: sibling's count_cycles_capped(adj, 64, cap=0) said EXISTS (0s);
my exhaustive randomized DFS (verify_III_dyadic_hunt_kps2) found none.

Method: re-run THEIR counter on the class-1 adjacency; then re-implement their
exact enumeration with path recording to extract the witness cycle; verify the
witness explicitly (64 distinct vertices, consecutive adjacency, closure).
Also instrument MY finder to report whether it exhausted or hit budget.

Output -> 05-knowledge/results/verify_III_g78_c64_kps2.out
"""
import sys, os, time
from collections import deque
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/verify_III_g78_c64_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

from verify_III_dyadic_hunt_kps2 import (G78_CLASS1, parse_adj, bfs_dist,
                                         find_cycle_len, verify_cycle)
from dyadic_cycle_checker_kps2 import count_cycles_len, bfs_dist_restricted
from dyadic_gap_hunt_kps2 import count_cycles_capped

A = G78_CLASS1

def their_counter_with_witness(adj, L):
    """EXACT copy of the sibling's count_cycles_capped DFS logic, but recording
    the path so the first counted cycle is returned as an explicit vertex list."""
    n = len(adj)
    for s in range(n):
        dist = bfs_dist_restricted(adj, s, n)
        on_path = [False] * n
        on_path[s] = True
        nbrs_s = [w for w in adj[s] if w > s and dist[w] <= L - 1]
        for v1 in nbrs_s:
            on_path[v1] = True
            stack = [(v1, 1, iter(adj[v1]))]
            path = [v1]
            while stack:
                u, depth, it = stack[-1]
                advanced = False
                for w in it:
                    if depth == L - 1:
                        if w == s and u > v1:
                            return [s] + path
                        continue
                    if w <= s or on_path[w]:
                        continue
                    if dist[w] > L - depth - 1:
                        continue
                    on_path[w] = True
                    stack.append((w, depth + 1, iter(adj[w])))
                    path.append(w)
                    advanced = True
                    break
                if not advanced:
                    on_path[u] = False
                    stack.pop()
                    path.pop()
            on_path[v1] = False
        on_path[s] = False
    return None

def main():
    log("=" * 90)
    log("G78 class 1: C64 disagreement resolution")
    log("=" * 90)

    t = time.time()
    r = count_cycles_capped(A, 64, cap=0, budget=80_000_000)
    log(f"sibling count_cycles_capped(A,64,cap=0): returned {r}  ({time.time()-t:.1f}s)")

    t = time.time()
    w = their_counter_with_witness(A, 64)
    log(f"witness extraction via their exact DFS logic: "
        f"{'FOUND' if w else 'NONE'}  ({time.time()-t:.1f}s)")
    if w:
        ok_len = len(w) == 64
        ok_distinct = len(set(w)) == len(w)
        ok_adj = all(w[(i + 1) % len(w)] in A[w[i]] for i in range(len(w)))
        log(f"  witness: {w}")
        log(f"  len==64: {ok_len}   all distinct: {ok_distinct}   closed cycle in A: {ok_adj}")
        if ok_len and ok_distinct and ok_adj:
            log("  => C64 EXISTS in G78 class 1. MY find_cycle_len was the buggy/limited one.")
        else:
            log("  => THEIR counter produced an invalid 'cycle' -- checker bug confirmed!")

    # my finder, instrumented
    log("")
    for seed in (3, 7, 11):
        t = time.time()
        res = find_cycle_len(A, 64, seed=seed, budget=200_000_000)
        kind = ("budget" if res == "budget" else
                ("found+valid" if verify_cycle(A, res, 64) else
                 ("found+INVALID" if res else "exhausted->None")))
        log(f"my find_cycle_len(A,64,seed={seed}): {kind}  ({time.time()-t:.1f}s)")
        if isinstance(res, list):
            log(f"  witness: {res}")
            break
    save()

if __name__ == "__main__":
    main()
