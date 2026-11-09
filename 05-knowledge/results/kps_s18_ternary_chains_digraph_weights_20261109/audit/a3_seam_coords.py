#!/usr/bin/env python3
"""AUDIT: is SEAM.md's (H, S, R) at x = 2 a fixed linear image of THM-5070's (h, E, O)?  (PROOF.md: "At x = 2,
SEAM.md's (H, S, R) transfer is the case treated here, in other coordinates, by OCF.")

SEAM glues the converse way, so for a port chain (T, u) the SEAM vector of T^op at exit u is
  H = #Hamiltonian paths of T, S = #Hamiltonian paths of T ENDING at u, R = sum over paths P of T from u of H(T - V(P)).
By OCF, H = h(2) and R = E(2) + O(2). The test: fit S = a h + b E + c O on random port chains and check every chain.
"""
import itertools
import random
import sys
from fractions import Fraction

sys.path.insert(0, __import__("os").path.dirname(__file__))
from a1_blocks import state_polys, ev, glue, classes  # noqa: E402


def ham_end(n, A, u):
    tot = endu = 0
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[i + 1]] for i in range(n - 1)):
            tot += 1
            endu += p[-1] == u
    return tot, endu


def main():
    rng = random.Random(7)
    blocks = [(B, w, up) for k in range(1, 4) for B in classes(k) for w in range(k) for up in range(k)]
    rows = []
    for t in range(60):
        T, u = [[0]], 0
        target = rng.randint(1, 8)
        while len(T) < target:
            B, w, up = blocks[rng.randrange(len(blocks))]
            if len(T) + len(B) > 8:
                break
            T, u = glue(T, u, B, w, up)
        h, E, O = (ev(p, 2) for p in state_polys(len(T), T, u))
        H, S = ham_end(len(T), T, u)
        assert H == h
        rows.append((h, E, O, S, len(T)))
    # least-squares-free exact fit on three independent rows, then check all
    import itertools as it
    fit = None
    for a, b, c in it.combinations(rows, 3):
        M = [[Fraction(r[0]), Fraction(r[1]), Fraction(r[2])] for r in (a, b, c)]
        det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
               + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        if det:
            # Cramer
            ys = [Fraction(r[3]) for r in (a, b, c)]
            sol = []
            for col in range(3):
                Mc = [row[:] for row in M]
                for i in range(3):
                    Mc[i][col] = ys[i]
                d = (Mc[0][0] * (Mc[1][1] * Mc[2][2] - Mc[1][2] * Mc[2][1]) - Mc[0][1] * (Mc[1][0] * Mc[2][2] - Mc[1][2] * Mc[2][0])
                     + Mc[0][2] * (Mc[1][0] * Mc[2][1] - Mc[1][1] * Mc[2][0]))
                sol.append(d / det)
            fit = sol
            break
    bad = sum(1 for r in rows if fit[0] * r[0] + fit[1] * r[1] + fit[2] * r[2] != r[3])
    print(f"60 random port chains of order <= 8: best exact fit S = {fit[0]} h + {fit[1]} E + {fit[2]} O fails on {bad} chains")
    print("sample rows (h, E, O, S, order):", rows[:8])


if __name__ == "__main__":
    main()
