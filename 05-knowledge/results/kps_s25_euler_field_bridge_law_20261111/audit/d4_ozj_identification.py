"""d4_ozj_identification.py - the family of addendum D (iterated direct joins of one-apex blocks at their apexes,
'apex-path tournaments') IS the O/Z/J root-change word family of THM-4904 / HYP-9494 (Sessions 197-200).

THM-4904, Setting: words over {O, Z, J} act on t = (D, E, U, d), from the rooted edge (1, 0, 2, 1):
    O = (D, 0, 2U, d),   Z = (D - E + U, d, 2U - d, d),   J = (D - E + U, D - E, D + U, D);
the value of a word is its final D.

Dictionary (audit): J = a new apex (one-apex block with the empty word), O = a spine letter 1 (the apex beats the
new spine vertex), Z = a spine letter 0; the rooted edge = the block with word '1'; (D, E, U) = (A, d, W) =
(H, o - e, o + e) of README B4.

Checks
  Z1  all O/Z/J words of length <= 9: D equals the literal Hamiltonian-path count of the apex-path tournament with
      blocks ('1' + letters before the first J), (letters between J's), ... (subset DP of tt.py), and (D, E, U)
      equals (H, o - e, o + e) at the last apex
  Z2  all O/Z/J words of length <= 15: D equals the first coordinate of M_(w_k) ... M_(w_0) (1, 0, 0)^T (README B4)
  Z3  value sets: {D of O/Z/J words} and {first coordinates of the block semigroup}, both closed below 2^11 by plain
      breadth-first search on states, have the same odd misses (the 28 numbers of runs/direct_closure_11.txt)
"""
import os
import sys
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tt import invariants, g_fun, word_PQ                    # noqa: E402
from b2_one_apex_chains import build                         # noqa: E402


def ozj(word, t=(1, 0, 2, 1)):
    D, E, U, d = t
    for c in word:
        if c == "O":
            D, E, U, d = D, 0, 2 * U, d
        elif c == "Z":
            D, E, U, d = D - E + U, d, 2 * U - d, d
        else:
            D, E, U, d = D - E + U, D - E, D + U, D
    return D, E, U, d


def blocks_of(word):
    parts = word.split("J")
    words = []
    for j, part in enumerate(parts):
        letters = tuple(1 if c == "O" else 0 for c in part)
        words.append(((1,) + letters) if j == 0 else letters)
    return words


def Mw(w):
    n = len(w)
    P, Q = word_PQ(w)
    return ((g_fun(P), Q + 1, -(1 if Q % 2 == 0 else 0)),
            (P + 1, 1 << n, 0),
            ((1 if P % 2 == 0 else 0), 0, -(1 if n == 0 else 0)))


def apply(M, v):
    return tuple(sum(M[r][c] * v[c] for c in range(3)) for r in range(3))


def z1():
    tot = bad = 0
    for L in range(0, 10):
        for word in product("OZJ", repeat=L):
            word = "".join(word)
            words = blocks_of(word)
            nv = sum(len(w) + 1 for w in words)
            if nv > 12:
                continue
            D, E, U, d = ozj(word)
            n, adj, apexes = build(words, [0] * (len(words) - 1))
            inv = invariants(n, adj, apexes[-1])
            tot += 1
            bad += ((D, E, U) != (inv["H"], inv["o"] - inv["e"], inv["o"] + inv["e"]))
    print("Z1  O/Z/J words of length <= 9 (at most 12 vertices): %d words, (D, E, U) = (H, o - e, o + e) of the literal apex-path tournament: %d failures" % (tot, bad))


def z2():
    tot = bad = 0
    for L in range(0, 13):
        for word in product("OZJ", repeat=L):
            word = "".join(word)
            v = (1, 0, 0)
            for w in blocks_of(word):
                v = apply(Mw(w), v)
            D, E, U, d = ozj(word)
            tot += 1
            bad += (v != (D, U, E))
    print("Z2  O/Z/J words of length <= 12: %d words, state (D, U, E) = M_(w_k)...M_(w_0)(1,0,0)^T: %d failures" % (tot, bad))


def z3(K=11):
    X = 1 << K
    # O/Z/J closure on states (D, E, U, d).  New values are F = D - E + U (after Z or J); O keeps D.  With the
    # invariants D - E >= d >= 1, D <= U, U >= 2d, E >= 0 (Session200; asserted below on every state) F never
    # decreases along a word (O: D + 2U; Z: D - E + 3U - 2d; J: D + 2U), so a state with F > X is dead: exact rule.
    seen = set()
    stack = [(1, 0, 2, 1)]
    seen.add(stack[0])
    hit = {1}
    while stack:
        D, E, U, d = stack.pop()
        assert E >= 0 and d >= 1 and D - E >= d and D <= U and U >= 2 * d
        for nxt in ((D, 0, 2 * U, d), (D - E + U, d, 2 * U - d, d), (D - E + U, D - E, D + U, D)):
            if nxt[0] > X:
                continue
            hit.add(nxt[0])
            if nxt[0] - nxt[1] + nxt[2] > X:
                continue
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    miss_ozj = [n for n in range(1, X + 1, 2) if n not in hit]
    print("Z3  O/Z/J closure below 2^%d: %d states, odd misses: %d %s" % (K, len(seen), len(miss_ozj), miss_ozj))
    path = os.path.join(HERE, "..", "runs", "direct_closure_%d.txt" % K)
    with open(path) as f:
        lines = f.read().splitlines()
    lst = [int(s) for s in lines[1].split("[")[1].split("]")[0].split(",")]
    print("    identical to the misses of the direct joins (runs/direct_closure_%d.txt): %s" % (K, lst == miss_ozj))
    s197 = [7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583, 637,
            821, 1067, 1339, 1861, 2647]
    print("    Session197's 29 values (C1_COMPLETE.md; THM-4884 (2b); HYP-9494 (D)) below 2^%d: identical: %s" % (K, [v for v in s197 if v <= X] == miss_ozj))


def z4(K=12):
    """the pure chains of README B2-B3 (a transitive chain of i >= 1 vertices before every new apex) are the O/Z/J
    words in which every J is immediately preceded by O (the rooted edge counts as ending in O).  Also: OJ = JZ as
    maps on states (Session205's JZ = OJ, THM-4888 (D)1), which is the normal form 'leading letters 0 are absorbed
    in the exponent'."""
    X = 1 << K
    bad = 0
    for D in range(1, 6):
        for E in range(0, 4):
            for U in range(1, 9):
                for d in range(1, 4):
                    oj = ozj("OJ", (D, E, U, d))
                    jz = ozj("JZ", (D, E, U, d))
                    bad += (oj != jz)
    print("Z4  OJ = JZ as maps on (D, E, U, d) (1440 integer states): %d failures" % bad)
    seen = set()
    start = ((1, 0, 2, 1), True)
    stack = [start]
    seen.add(start)
    hit = {1}
    while stack:
        (D, E, U, d), lastO = stack.pop()
        moves = [((D, 0, 2 * U, d), True), ((D - E + U, d, 2 * U - d, d), False)]
        if lastO:
            moves.append(((D - E + U, D - E, D + U, D), False))
        for nxt, flag in moves:
            if nxt[0] > X:
                continue
            hit.add(nxt[0])
            if nxt[0] - nxt[1] + nxt[2] > X:
                continue
            if (nxt, flag) not in seen:
                seen.add((nxt, flag))
                stack.append((nxt, flag))
    miss = [n for n in range(1, X + 1, 2) if n not in hit]
    path = os.path.join(HERE, "b3_misses_closure_%d.txt" % K)
    ref = None
    if os.path.exists(path):
        with open(path) as f:
            ref = [int(s) for s in f.read().splitlines()[1].split()]
    print("    O/Z/J words with every J preceded by O, closure below 2^%d: %d states, odd misses: %d; identical to the misses of the pure chain semigroup (b3_misses_closure_%d.txt): %s"
          % (K, len(seen), len(miss), K, ref == miss))


def ham_cycles(n, adj):
    """number of directed Hamiltonian cycles (n >= 3), by the subset DP from vertex 0"""
    full = (1 << n) - 1
    dp = {(1, 0): 1}
    for mask in range(1, full + 1, 2):
        for last in range(n):
            c = dp.get((mask, last))
            if not c:
                continue
            nxt = adj[last] & ~mask
            w = 0
            while nxt:
                if nxt & 1:
                    key = (mask | (1 << w), w)
                    dp[key] = dp.get(key, 0) + c
                nxt >>= 1
                w += 1
    return sum(c for (mask, last), c in dp.items() if mask == full and (adj[last] & 1))


def z5():
    """apex-path tournaments have at most one Hamiltonian cycle (hc(J) = s_a(A) e_b(B), THM-4866 (B)); the strong ones
    are uniquely Hamiltonian: this is why the family is the 'C = 1' family of Sessions 197-200."""
    from tt import H
    tot = strong = bad = 0
    for L in range(1, 9):
        for word in product("OZJ", repeat=L):
            words = blocks_of("".join(word))
            n, adj, apexes = build(words, [0] * (len(words) - 1))
            if n < 3 or n > 10:
                continue
            hc = ham_cycles(n, adj)
            tot += 1
            strong += (hc >= 1)
            bad += (hc > 1)
    print("Z5  O/Z/J words of length <= 8 (3 to 10 vertices): %d tournaments, %d with a Hamiltonian cycle (strong), %d with more than one Hamiltonian cycle" % (tot, strong, bad))


if __name__ == "__main__":
    z1()
    z2()
    z3()
    z4()
    z5()
