"""Replay 'WORD <hex rows> : <path>; <path>; ...' lines (geodesic_search_audit output) in pure Python: the start must be
regular, the word must have exactly c = (n^2-1)/8 paths, each a directed Hamilton path of the current state; the native
move reverses every arc not on the path; L must drop by one per move and the final state must be transitive.
Also converts gentourng upper-triangle lines to hex rows:  python3 replay_words.py --tohex n < gentourng_lines
Usage: python3 replay_words.py n < search_output
"""
import sys
if sys.argv[1] == "--tohex":
    n = int(sys.argv[2])
    for line in sys.stdin:
        s = line.strip()
        if not s or s[0] not in "01": continue
        o = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[k] == "1": o[i] |= 1 << j
                else: o[j] |= 1 << i
                k += 1
        print(" ".join(f"{x:x}" for x in o))
    sys.exit(0)
n = int(sys.argv[1]); m = (n - 1) // 2; c = (n * n - 1) // 8
def L(A):
    s = sorted(sum(r) for r in A); return sum(abs(s[i] - i) for i in range(n)) // 2
ok_count = bad = 0
for line in sys.stdin:
    if not line.startswith("WORD"): continue
    head, tail = line[4:].split(":", 1)
    rows = [int(x, 16) for x in head.split()]
    A = [[(rows[i] >> j) & 1 for j in range(n)] for i in range(n)]
    good = all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(n) if i != j) and all(sum(r) == m for r in A)
    paths = [tuple(int(x) for x in p.split(",")) for p in tail.split("[")[0].split(";") if p.strip()]
    good &= len(paths) == c
    Ls = [L(A)]
    for Q in paths:
        if sorted(Q) != list(range(n)) or not all(A[Q[i]][Q[i + 1]] for i in range(n - 1)): good = False; break
        on = {(Q[i], Q[i + 1]) for i in range(n - 1)}
        A = [[0 if i == j else (1 if (i, j) in on else 0 if (j, i) in on else 1 - A[i][j]) for j in range(n)] for i in range(n)]
        Ls.append(L(A))
    good &= Ls == list(range(c, -1, -1)) and sorted(sum(r) for r in A) == list(range(n))
    ok_count += good; bad += not good
print(f"n={n}: words replayed {ok_count + bad}; valid (regular start, {c} moves, Hamilton paths at every stage, L = {c}..0, transitive end): {ok_count}; invalid: {bad}")
