"""cd_orbit.py - the group of sign changes that preserves the parity of every zero-sum 4-set (S23).

Sign patterns on the nonzero vectors of Z_2^n = tournaments; x_{uv} in F_2 (u < v).  Two families of moves:
  switching at a vertex u            : flips every arc at u                    y_{vw} = s_v + s_w
  line coboundary at a point f       : flips every arc {v, w} with v + w = f   y_{vw} = t_(v+w)
(rescaling the basis, e_u -> -e_u, is the product of the switching at u and the coboundary at u... no: it flips
 the arcs {u, w} and the arcs {v, w} with v + w = u; it is s = t = indicator of u.)
Claim checked here for n = 3, 4: the moves span exactly the kernel of the parity map (all homogeneous equations),
of dimension 2 (2^n - 1) - (n + 1); so two tournaments have the same violation pattern iff they differ by moves."""
from cd_more import F2System, quads


def rank(vectors):
    piv = {}
    for v in vectors:
        while v:
            hb = v.bit_length() - 1
            if hb in piv:
                v ^= piv[hb]
            else:
                piv[hb] = v
                break
    return len(piv), piv


for n in (3, 4):
    S = F2System(n)
    N = 1 << n
    nv = len(S.var)
    moves = []
    for u in range(1, N):
        moves.append(sum(1 << k for (a, b), k in S.var.items() if a == u or b == u))
    for f in range(1, N):
        moves.append(sum(1 << k for (a, b), k in S.var.items() if a ^ b == f))
    rk_moves, _ = rank(moves)
    rows = []
    for q in quads(n):
        rows += [m for (m, c) in S.quad_rows(q)]
    rk_eq, _ = rank(rows)
    killed = all(bin(mv & r).count("1") % 2 == 0 for mv in moves for r in rows)
    print("n = %d: unknowns %d, rank of the parity equations %d, nullity %d; moves span dimension %d (expected %d); every move preserves every equation: %s"
          % (n, nv, rk_eq, nv - rk_eq, rk_moves, 2 * (N - 1) - (n + 1), killed))
