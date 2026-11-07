"""Literal witnesses for the root-change depth table (kind-pasteur-2026-11-07-S14), Python stdlib only.

For each target N: (1) an exhaustive word search (letters O, Z, J; F strictly increases along every extension, so
pruning at F > N is exact) finds the minimal number of root changes J among all words w with F(t.w) = N;
(2) a witness word w of that depth is closed by J and built as a literal tournament with Session200's builder
(new vertex beaten by every old vertex for O; for Z and J it beats the current root; J moves the root);
(3) bitmask dynamic programming counts Hamiltonian paths H and Hamiltonian cycles C of that tournament, which must be
(N, 1): the strong C = 1 readout.  Every check raises RuntimeError on failure; nothing uses assert.

usage: python3 literal_witnesses.py > LITERAL_WITNESSES.txt
"""
import sys

T = (1, 0, 2, 1)
TARGETS = [1877, 1897, 667, 3523, 373, 1811, 2149, 2999, 4411, 4547, 7129, 17383, 19951, 22091, 27781, 205, 281, 313]
# 1877, 1897 depth 7; 667, 3523 depth 5; 373 ... 27781 depth 4; 205, 281, 313: F_<=2 misses that F_<=3 covers.
C1_MISSES = [7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583,
             637, 821, 1067, 1339, 1861, 2647]


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def step(s, c):
    D, E, U, d = s
    if c == 'O':
        return (D, 0, 2 * U, d)
    if c == 'Z':
        return (D - E + U, d, 2 * U - d, d)
    return (D - E + U, D - E, D + U, D)


def F(s):
    return s[0] - s[1] + s[2]


def min_depth_word(N, kmax=None):
    """Minimal #J over words w with F(t.w) = N, and one witness of that depth (exhaustive, exact pruning)."""
    best = [None, None]
    stack = [(T, '', 0)]
    while stack:
        s, w, j = stack.pop()
        f = F(s)
        if f > N:
            continue
        if f == N and (best[0] is None or j < best[0] or (j == best[0] and len(w) < len(best[1]))):
            best = [j, w]
        for c in 'OZJ':
            jj = j + (c == 'J')
            if kmax is None or jj <= kmax:
                q = step(s, c)
                need(F(q) > f, 'F must increase strictly')
                stack.append((q, w + c, jj))
    return best


def graph(w):
    rows, root = [2, 0], 0
    for c in w:
        n = len(rows)
        rows.append(0)
        for v in range(n):
            if c != 'O' and v == root:
                rows[n] |= 1 << v
            else:
                rows[v] |= 1 << n
        if c == 'J':
            root = n
    return rows


def ham_counts(rows):
    """(Hamiltonian paths, Hamiltonian cycles) by subset DP; cycles fixed to start at vertex 0."""
    n = len(rows)
    full = (1 << n) - 1
    paths = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        paths[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = paths[S]
        for v in range(n):
            x = row[v]
            if x:
                out = rows[v] & ~S
                while out:
                    b = out & -out
                    u = b.bit_length() - 1
                    paths[S | b][u] += x
                    out ^= b
    H = sum(paths[full])
    cyc = [[0] * n for _ in range(1 << n)]
    cyc[1][0] = 1
    for S in range(1, 1 << n, 2):
        row = cyc[S]
        for v in range(n):
            x = row[v]
            if x:
                out = rows[v] & ~S
                while out:
                    b = out & -out
                    u = b.bit_length() - 1
                    cyc[S | b][u] += x
                    out ^= b
    C = sum(cyc[full][v] for v in range(1, n) if rows[v] & 1)
    return H, C


def main():
    print('target N, minimal root changes, witness word w (closed by J), vertices, literal (H, C)')
    for N in TARGETS:
        k, w = min_depth_word(N)
        need(k is not None, 'no word for %d' % N)
        rows = graph(w + 'J')
        H, C = ham_counts(rows)
        need((H, C) == (N, 1), 'literal check failed for %d: %s' % (N, (H, C)))
        print('%d %d %s %d (%d, %d)' % (N, k, w, len(rows), H, C))
    for N in C1_MISSES:
        k, w = min_depth_word(N)
        need(k is None, 'C1 miss %d unexpectedly has a word %s' % (N, w))
    print('the 29 C1 misses have no word with any number of root changes (exhaustive: F increases by >= 2 per letter)')
    print('RESULT PASS')


if __name__ == '__main__':
    sys.setrecursionlimit(10000)
    main()
