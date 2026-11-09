"""Evidence for HYP-9540 (A) at n = 3k + 1: the run-window staircase realises every tight colouring.

For a tight colouring w of [0, n-3] with maximal runs [p, q] of 1s, the run-window tiling reverses exactly the tiles
inside the vertex windows [p, q + 2]. At n = 3k + 2 this is THM-4951 (3b)'s multi-staircase. The script checks, for
every tight colouring, whether the tiling has triple word w and is tight (nu(t) + nu(t + 1) = floor(n/3)).

Expected:
- It works for every tight colouring at n = 3k + 1, k <= 5, and at n = 6.
- It fails at n = 9 and 12 on 1 and 3 words, each containing 1110111: an isolated 0 between two 1-runs of length >= 2.
  The follow-up audit showed this is exactly the failure set for every n = 3k <= 18. An isolated 0 next to a 1-run
  of length 1, as in 1010, is harmless. The failing words are still realised by other tight lines (THM-4951 (4)).
usage: python run_window_check.py
"""
import os
import importlib.util

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("sc", os.path.join(here, "staircase_check.py"))
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)


def P(xs):
    c, last = 0, -10
    for x in xs:
        if x - last >= 3:
            c += 1
            last = x
    return c


def tight(L):
    k = -(-L // 3)
    out = []
    for w in range(1 << L):
        if not (w & 1):
            continue
        A = [i for i in range(L) if w >> i & 1]
        B = [i for i in range(L) if not w >> i & 1]
        if P(A) + P(B) == k:
            out.append(''.join('1' if w >> i & 1 else '0' for i in range(L)))
    return out


def main():
    for n in (4, 7, 10, 13, 16, 6, 9, 12):
        k, L = n // 3, n - 2
        good, bad = 0, []
        for w in tight(L):
            runs, i = [], 0
            while i < L:
                if w[i] == '1':
                    s = i
                    while i < L and w[i] == '1':
                        i += 1
                    runs.append((s, i - 1))
                else:
                    i += 1
            rev = sc.inside([(p, q + 2) for (p, q) in runs], n)
            t, u = sc.tournament(n, rev), sc.tournament(n, sc.complement(n, rev))
            a, b = sc.nu(t, n), sc.nu(u, n)
            word = ''.join('1' if (x, x + 2) in rev else '0' for x in range(L))
            if a + b == k and word == w:
                good += 1
            else:
                bad.append((w, a, b))
        print(f"n = {n}: tight colourings {good + len(bad)}, run-window tight with the right word {good}, "
              f"fails {len(bad)} {bad[:4]}")


if __name__ == "__main__":
    main()
