"""AUDIT: the number K_q of Hamiltonian paths of the Paley tournament P_q that start with the arc 0 -> 1, modulo 256,
q = 7, 11, 19, 23, 27 (run F / R8 of the note). K_q = number of Hamiltonian paths of P_q minus the vertex 0 that
start at the vertex 1 (aud_ham8.c).   usage: python3 aud_paley27.py /path/aud_ham8"""
import subprocess, sys, time
from aud_common import gf_field, paley

HAM8 = sys.argv[1]


def K_mod256(q):
    A = paley(q)
    assert A[0][1] == 1
    idx = [v for v in range(q) if v != 0]
    B = [[A[a][b] for b in idx] for a in idx]
    text = str(len(B)) + "\n" + "\n".join("".join(map(str, r)) for r in B) + "\n"
    return int(subprocess.run([HAM8, "0"], input=text, capture_output=True, text=True, check=True).stdout)


known = {7: 9, 11: 1729, 19: 6857869865, 23: 62293308207033}
for q in (7, 11, 19, 23, 27):
    t = time.time()
    k = K_mod256(q)
    e = q * (q - 1) // 2
    line = "q = %2d: paths starting 0 -> 1, mod 256: %3d" % (q, k)
    if q in known:
        line += " (exact value %d = %d mod 256: %s)" % (known[q], known[q] % 256, known[q] % 256 == k)
    line += "; H(P_q) = q(q-1)/2 * K = %d (mod 256), %d (mod 8); q(q-1)/2 = %d (mod 8)   [%.0f s]" % ((e * k) % 256, (e * k) % 8, e % 8, time.time() - t)
    print(line)
