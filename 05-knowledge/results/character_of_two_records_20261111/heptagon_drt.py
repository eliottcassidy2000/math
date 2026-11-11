#!/usr/bin/env python3
"""The number of 7-cycles of a doubly regular tournament of order n = 4t + 3 (section 5.4 of the note).

For every tournament:  tr(A^7) = 7 c_7 + 7 (N_1 + N_2 + N_3),  N_j = number of rotation classes of closed 7-walks
with exactly j repeated vertices, and  sum_u t_u f_u = N_1 + 2 N_2 + 3 N_3  (t_u, f_u = numbers of 3- and 4-cycles
through u).  For a doubly regular tournament:
    tr(A^7) = (2t+1) ((2t+1)^6 + 7q^3 - 14q^2 + 7q - 1),  q = t + 1,
    sum_u t_u f_u = 2 n t (t+1)^2 (2t+1)^2,   N_2 = (3/2) n (n-1) t^2 (t+1),   N_3 = 3 t c_3,   c_3 = n (n^2 - 1)/24,
    c_7 = tr(A^7)/7 - sum_u t_u f_u + N_2 + 2 N_3.
Checks: the formula against direct counts (cycle_lengths.c, modulo 256; exact for P_7); integrality and evenness for
t = 0..40 (evenness at 8 consecutive t proves it for all t, the polynomial having degree 7).
Usage: heptagon_drt.py CYC      (CYC = compiled cycle_lengths.c)"""
import sys, json
from fractions import Fraction
from paley_tools import paley_out
from cayley_dickson import tournament
from drt_zoo import gs_solutions, goethals_seidel, descendant, is_drt
from cycle_parities import cycles
def c7(t):
    n = 4 * t + 3; q = t + 1
    tr7 = (2 * t + 1) * ((2 * t + 1) ** 6 + 7 * q ** 3 - 14 * q ** 2 + 7 * q - 1)
    c3 = Fraction(n * (n * n - 1), 24)
    stf = 2 * n * t * (t + 1) ** 2 * (2 * t + 1) ** 2
    N2 = Fraction(3 * n * (n - 1) * t * t * (t + 1), 2); N3 = 3 * t * c3
    return Fraction(tr7, 7) - stf + N2 + 2 * N3
if __name__ == "__main__":
    exe = sys.argv[1]
    vals = [c7(t) for t in range(0, 41)]
    print("integral for t = 0..40:", all(v.denominator == 1 for v in vals), "; even for t = 0..40:", all(int(v) % 2 == 0 for v in vals))
    print("c_7 for n = 3, 7, 11, 15, 19, 23, 27, 31:", [int(v) for v in vals[:8]])
    a, b, c, d = next(gs_solutions(5)); D19 = descendant(goethals_seidel(a, b, c, d), 0)
    a, b, c, d = next(gs_solutions(7)); G27 = descendant(goethals_seidel(a, b, c, d), 0)
    for name, out in [("P_3", paley_out(3)), ("P_7", paley_out(7)), ("P_11", paley_out(11)), ("CD_4", tournament(4)[0]), ("P_19", paley_out(19)),
                      ("D_19", D19), ("P_23", paley_out(23)), ("P_27", paley_out(27)), ("a Goethals-Seidel class of order 27", G27), ("P_31", None)]:
        if out is None: continue
        assert is_drt(out) or len(out) == 3
        n = len(out); t = (n - 3) // 4; got = cycles(exe, out).get(7, 0) if n >= 7 else 0
        print("%-36s n = %2d: 7-cycles (mod 256) counted %3d, formula %d = %3d (mod 256): %s" % (name, n, got, int(c7(t)), int(c7(t)) % 256, "ok" if got == int(c7(t)) % 256 else "MISMATCH"))
