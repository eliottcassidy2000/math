"""b6_xe_identity.py - the dictionary of README B1 / B4: 'x_e = o - e'.

For a rooted tournament (X, u): o (e) = sum of H(X - P) over the paths P starting at u with an odd (even) number of
vertices.  Claim used by the note: o - e is the third coordinate of the rooted triple of THM-4866, i.e. the number
of Hamiltonian paths of X that END at u.  Checked on every labeled rooted tournament with at most 6 vertices
(and the dual statement for the split of W-).
"""
from tt import tables, invariants, all_tournaments


def main():
    tot = bad = bad2 = 0
    for n in range(1, 7):
        for adj in all_tournaments(n):
            tabs = tables(n, adj)
            hsub, st, en = tabs
            full = (1 << n) - 1
            for u in range(n):
                inv = invariants(n, adj, u, tabs)
                tot += 1
                bad += (inv["o"] - inv["e"] != int(en[full, u]))
                bad2 += (inv["Om"] - inv["Em"] != int(st[full, u]))
    print("rooted labeled tournaments with at most 6 vertices: %d; o - e = number of Hamiltonian paths ending at the root: %d failures; "
          "(odd - even part of W-) = number of Hamiltonian paths starting at the root: %d failures" % (tot, bad, bad2))


if __name__ == "__main__":
    main()
