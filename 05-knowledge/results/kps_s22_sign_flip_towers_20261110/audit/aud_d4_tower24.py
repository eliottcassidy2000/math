"""aud_d4_tower24.py - D4: build the four-closer tower tournaments of order 24 explicitly (start (TT_n, {1}),
levels (1, 8-n, 0), (1, 4, 0), (2, 2, 1), closing tail t = 1) and count their Hamiltonian paths by LITERAL
generation of every path (aud_pairstate mode 3: depth-first extension along arcs, no tables, no formula).
Expected: 2^24 - 1 = 16777215.   usage: python aud_d4_tower24.py n1 n2 ..."""
import subprocess, os, sys, time
from aud_core import transitive, close, is_tournament

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "aud_pairstate.exe")
for n in map(int, sys.argv[1:]):
    A = transitive(n)
    O = [0]
    for (d, i, m) in ((1, 8 - n, 0), (1, 4, 0), (2, 2, 1)):
        a = i + d + m
        nn = len(A)
        A = close(A, O, transitive(a))
        O = [nn + i - 1, nn + i - 1 + d]
    T = close(A, O, transitive(1))
    N = len(T)
    assert is_tournament(T) and N == 24
    line = " ".join(["3", str(N)] + [str(sum(1 << j for j in range(N) if T[i][j])) for i in range(N)] + ["0"])
    t0 = time.time()
    out = subprocess.run([EXE], input=line + "\n", capture_output=True, text=True).stdout.strip()
    H = int(out.split()[1])
    print("n = %d: order %d, H (literal generation) = %d ; 2^24 - 1 = %d ; %s  (%.0fs)"
          % (n, N, H, 2 ** 24 - 1, "PASS" if H == 2 ** 24 - 1 else "FAIL", time.time() - t0))
    sys.stdout.flush()
