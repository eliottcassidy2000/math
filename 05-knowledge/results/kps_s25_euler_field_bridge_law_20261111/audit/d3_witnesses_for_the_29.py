"""d3_witnesses_for_the_29.py - the 29 odd numbers missed by the direct joins of one-apex blocks (addendum D) are,
except 7 and 21, Hamiltonian-path counts of OTHER tournaments: explicit witnesses, each recounted by the subset
dynamic programme of tt.py on the literal adjacency matrix.

Search (auditor's own):
  S1  all labeled tournaments with at most 6 vertices (plain witnesses; also shows 7 and 21 do not occur there)
  S2  T = (X, u) => Y_w with the arc between u and the apex of Y_w reversed, X any rooted tournament with at most 5
      vertices, Y_w a one-apex tournament with word w of length <= 9 (candidates by the law D1, then literal recount)
Witness format: n ; adjacency rows as bitmasks (bit j of row i set iff i -> j).
"""
from itertools import product

from tt import tables, invariants, ordinal_sum, reverse_arc, one_apex, all_tournaments, g_fun, word_PQ, H

TARGETS = [7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583,
           637, 821, 1067, 1339, 1861, 2647]


def main():
    found = {}
    spectrum = {}
    for n in range(1, 7):
        vals = set()
        for adj in all_tournaments(n):
            h = H(n, adj)
            vals.add(h)
            if h in TARGETS and h not in found:
                found[h] = ("plain tournament of order %d" % n, n, list(adj))
        spectrum[n] = sorted(vals)
    print("S1  H-spectra of all labeled tournaments: " + "; ".join("n=%d: %s" % (n, spectrum[n]) for n in range(1, 6)))
    print("    n=6:", spectrum[6])
    print("    7 or 21 attained with at most 6 vertices:", any(v in (7, 21) for n in spectrum for v in spectrum[n]))

    bank = {}
    for n in range(1, 6):
        for adj in all_tournaments(n):
            tabs = tables(n, adj)
            for u in range(n):
                inv = invariants(n, adj, u, tabs)
                key = (inv["H"], inv["o"], inv["e"])
                if key not in bank:
                    bank[key] = (n, list(adj), u)
    print("S2  bank: %d distinct rooted states (A, o, e) of rooted tournaments with at most 5 vertices" % len(bank))
    words = [w for L in range(0, 10) for w in product((0, 1), repeat=L)]
    wdata = []
    for w in words:
        P, Q = word_PQ(w)
        wdata.append((w, g_fun(P), 1 + Q // 2, (Q + 1) // 2))
    todo = set(t for t in TARGETS if t not in found)
    for (A, o, e), (n, adj, u) in sorted(bank.items()):
        if not todo:
            break
        for (w, g, O, E) in wdata:
            N = A * g + 2 * (o * E + e * O)
            if N in todo:
                m, a2, offs = ordinal_sum([(n, adj), one_apex(w)])
                reverse_arc(a2, offs[0] + u, offs[1])
                if H(m, a2) == N:
                    found[N] = ("join of a rooted %d-vertex tournament (root %d) with the one-apex block %s" % (n, u, "".join(map(str, w)) or "(empty word)"), m, a2)
                    todo.discard(N)
    # S3: direct join of two rooted tournaments with at most 6 vertices each (arc between the roots reversed):
    #     H = H(X) H(Y) + 2 (o_X Em_Y + e_X Om_Y),  Om / Em = W-_Y(v) split by the parity of the number of vertices
    if todo:
        left, right = {}, {}
        for n in range(1, 7):
            for adj in all_tournaments(n):
                tabs = tables(n, adj)
                for u in range(n):
                    inv = invariants(n, adj, u, tabs)
                    kl = (inv["H"], inv["o"], inv["e"])
                    kr = (inv["H"], inv["Om"], inv["Em"])
                    if kl not in left:
                        left[kl] = (n, list(adj), u)
                    if kr not in right:
                        right[kr] = (n, list(adj), u)
        print("S3  banks of rooted tournaments with at most 6 vertices: %d left states (H, o, e), %d right states (H, O-, E-)" % (len(left), len(right)))
        for (A, o, e), (n, adj, u) in sorted(left.items()):
            if not todo:
                break
            for (HY, Om, Em), (m2, adj2, v) in right.items():
                N = A * HY + 2 * (o * Em + e * Om)
                if N in todo:
                    m, a2, offs = ordinal_sum([(n, adj), (m2, adj2)])
                    reverse_arc(a2, offs[0] + u, offs[1] + v)
                    if H(m, a2) == N:
                        found[N] = ("join X => Y of rooted tournaments of orders %d and %d (roots %d and %d, arc between the roots reversed)" % (n, m2, u, v), m, a2)
                        todo.discard(N)
    print()
    ok = 0
    for t in TARGETS:
        if t in found:
            desc, n, adj = found[t]
            recount = H(n, adj)
            ok += (recount == t)
            print("H = %4d : %s ; n = %d ; rows = %s ; literal recount = %d" % (t, desc, n, adj, recount))
        else:
            print("H = %4d : no witness found in this search%s" % (t, " (expected: 7 and 21 are not path counts, THM-1370)" if t in (7, 21) else ""))
    print()
    print("witnesses found and recounted: %d of %d targets other than 7 and 21" % (ok, len(TARGETS) - 2))


if __name__ == "__main__":
    main()
