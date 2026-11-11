#!/usr/bin/env python3
"""The path form of conjecture (bK) (section 5.3): numbers p_j of directed paths on j vertices, modulo 4, for doubly
regular tournaments, against C(n, j) (n = 7 mod 8) and C(n, j) + 2 [j odd, j >= 3] C(m, (j-1)/2) (n = 3 mod 8)."""
from math import comb
from paley_tools import paley_out
from cayley_dickson import tournament
from local_forms import endpoint_counts
from drt_zoo import gs_solutions, goethals_seidel, descendant
def predicted(n, j):
    m = (n - 1) // 2
    return (comb(n, j) + (2 * comb(m, (j - 1) // 2) if (n % 8 == 3 and j % 2 == 1 and j >= 3) else 0)) % 4
if __name__ == "__main__":
    a, b, c, d = next(gs_solutions(5)); D19 = descendant(goethals_seidel(a, b, c, d), 0)
    for name, out in [("P_3", paley_out(3)), ("P_7", paley_out(7)), ("P_11", paley_out(11)), ("CD_4", tournament(4)[0]),
                      ("P_19", paley_out(19)), ("D_19", D19)]:
        n = len(out); E = endpoint_counts(out, n); got = [sum(E[k]) % 4 for k in range(1, n + 1)]
        pred = [predicted(n, j) for j in range(1, n + 1)]
        print("%-5s n = %2d: p_j mod 4 (j = 1..n): %s  predicted: %s  %s   (C(n, j) mod 4: %s)" % (
            name, n, ''.join(map(str, got)), ''.join(map(str, pred)), "ok" if got == pred else "DIFFERENT",
            ''.join(str(comb(n, j) % 4) for j in range(1, n + 1))))
