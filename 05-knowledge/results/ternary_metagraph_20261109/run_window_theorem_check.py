"""Brute-force checks of THM-4951 (3e) and (3f) (added after the main audit).

(3f) Run-window theorem. Let w be a colouring of the positions [0, n-3] in which every maximal run (of 0s or of 1s)
has length >= 2. Let t_w reverse exactly the tiles inside the windows [p, q + 2] of its 1-runs [p, q]. Then:
- the triple word of t_w is w;
- nu(t_w) = P(A) and nu(t_w + 1) = P(B), with P(X) = sum over the runs of X of ceil(length / 3);
- F is a feedback vertex set of t_w + 1 with |F| = P(B). F is the union over the 0-runs [p0, q0] of the points
  min(p0 + 3i + 2, q0 + 1), for i >= 0 with p0 + 3i <= q0.
(3e) At n = 3k + 1:
- the tight colourings are exactly the block words for a decomposition of [0, 3k - 2] into consecutive blocks of
  length 3 except one block of length 2;
- their number up to complement is (k + 1) 2^(k-2).
Checks: (3f) for every colouring with all runs >= 2 at n = 5..13; (3e) by brute force for k = 1..6.
usage: python run_window_theorem_check.py
"""
import os
import itertools
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


def runs_of(w):
    out, i, L = [], 0, len(w)
    while i < L:
        s = i
        while i < L and w[i] == w[s]:
            i += 1
        out.append((w[s], s, i - 1))
    return out


def main():
    bad = 0
    total = 0
    for n in range(5, 14):
        L = n - 2
        for bits in itertools.product('01', repeat=L):
            w = ''.join(bits)
            R = runs_of(w)
            if any(q - p + 1 < 2 for (_, p, q) in R):
                continue
            total += 1
            A = [i for i in range(L) if w[i] == '1']
            B = [i for i in range(L) if w[i] == '0']
            wins = [(p, q + 2) for (ch, p, q) in R if ch == '1']
            rev = sc.inside(wins, n)
            t, u = sc.tournament(n, rev), sc.tournament(n, sc.complement(n, rev))
            word = ''.join('1' if (x, x + 2) in rev else '0' for x in range(L))
            F = set()
            for (ch, p0, q0) in R:
                if ch == '0':
                    i = 0
                    while p0 + 3 * i <= q0:
                        F.add(min(p0 + 3 * i + 2, q0 + 1))
                        i += 1
            alive = ((1 << n) - 1) & ~sum(1 << f for f in F)
            sumA = sum(-(-(q - p + 1) // 3) for (ch, p, q) in R if ch == '1')
            sumB = sum(-(-(q - p + 1) // 3) for (ch, p, q) in R if ch == '0')
            ok = (word == w and P(A) == sumA and P(B) == sumB and sc.nu(t, n) == sumA and sc.nu(u, n) == sumB
                  and len(F) == sumB and sc.acyclic(u, n, alive))
            if not ok:
                bad += 1
                if bad <= 5:
                    print("FAIL (3f)", n, w, sc.nu(t, n), sc.nu(u, n), sumA, sumB, sorted(F))
    print(f"(3f) run-window theorem: {total} colourings with all runs >= 2 at n = 5..13; failures {bad}")
    bad = 0
    for k in range(1, 7):
        L = 3 * k - 1
        tight = set()
        for bits in itertools.product('01', repeat=L):
            w = ''.join(bits)
            A = [i for i in range(L) if w[i] == '1']
            B = [i for i in range(L) if w[i] == '0']
            if P(A) + P(B) == k:
                tight.add(w)
        blockw = set()
        for s in range(k):
            sizes = [3] * k
            sizes[s] = 2
            for beta in itertools.product('01', repeat=k):
                blockw.add(''.join(beta[b] * sizes[b] for b in range(k)))
        half = len([w for w in tight if w[0] == '1'])
        ok = tight == blockw and half == (k + 1) * 2 ** k // 4
        bad += not ok
        print(f"(3e) k = {k}, n = {3 * k + 1}: tight colourings {len(tight)}, block words with one short block "
              f"{len(blockw)}, equal {tight == blockw}; up to complement {half} vs (k+1) 2^(k-2) = {(k + 1) * 2 ** k / 4}")
    print(f"(3e) failures {bad}")


if __name__ == "__main__":
    main()
