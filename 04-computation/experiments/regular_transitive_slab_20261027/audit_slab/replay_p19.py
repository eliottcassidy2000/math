"""Independent replay of the 45-move Paley P19 word of THM-4819 (C) (paley19_geodesic.out, line 'word: ...').

Start: P19 on Z_19 with i -> j iff j - i is a nonzero square mod 19 (squares computed here, not read from the file).
For every listed path: it must be a permutation of 0..18 and a directed Hamilton path of the CURRENT tournament; the native
move reverses every arc not on the path.  Checks L = (1/2) sum |s_(i) - i| after every move (THM-4810: |dL| <= 1, so a
45-move word to L = 0 must drop L by exactly one each time), that the final tournament is transitive (acyclic: its
scores are 0..18), and reports the colour of every move (black / blue) for information.
Usage: python3 replay_p19.py path/to/paley19_geodesic.out
"""
import sys
n = 19
squares = {(x * x) % n for x in range(1, n)}
A = [[1 if i != j and (j - i) % n in squares else 0 for j in range(n)] for i in range(n)]
assert all(sum(A[i]) == 9 for i in range(n)) and all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(n) if i != j)
word_line = [l for l in open(sys.argv[1]) if l.startswith("word:")]
assert len(word_line) == 1
paths = [tuple(int(x) for x in p.split(",")) for p in word_line[0][len("word:"):].strip().split(";") if p.strip()]
def L(A):
    s = sorted(sum(r) for r in A); return sum(abs(s[i] - i) for i in range(n)) // 2
def blue(A, Q):
    sig = {Q[i]: Q[n - 1 - i] for i in range(n)}
    return all(A[x][y] == A[sig[y]][sig[x]] for x in range(n) for y in range(n) if x != y)
Ls = [L(A)]; colours = []
for k, Q in enumerate(paths):
    assert sorted(Q) == list(range(n)), f"move {k+1}: not a permutation"
    assert all(A[Q[i]][Q[i + 1]] for i in range(n - 1)), f"move {k+1}: not a directed Hamilton path of the current state"
    colours.append("blue" if blue(A, Q) else "black")
    on = {(Q[i], Q[i + 1]) for i in range(n - 1)}
    B = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if (i, j) in on: B[i][j] = 1
            elif (j, i) in on: B[i][j] = 0
            else: B[i][j] = 1 - A[i][j]
    A = B; Ls.append(L(A))
scores = sorted(sum(r) for r in A)
acyclic = scores == list(range(n))
print(f"moves {len(paths)}; every path a directed Hamilton path at its stage: True; L sequence {Ls}")
print(f"L drops by exactly one per move: {all(Ls[i] - Ls[i+1] == 1 for i in range(len(paths)))}; final transitive: {acyclic}")
print(f"colours: {colours.count('black')} black, {colours.count('blue')} blue")
