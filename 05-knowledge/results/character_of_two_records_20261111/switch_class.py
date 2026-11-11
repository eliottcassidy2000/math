#!/usr/bin/env python3
"""Run H1, H2 (and the exploratory in-sample runs): H mod 4 on the switching class of T^+ = (a new source) + T for a
doubly regular tournament T.  The members of the class are the tournaments whose +-1 skew matrix S satisfies
S^2 = -(N-1) I and which are switching-equivalent to T^+.

  switch_class.py insample EXACT            P_3, P_7 (all members), P_11, CD_4, P_19 (random members), exact H
  switch_class.py h1 LAYERED [members]      N = 24: classes of P_23^+ and of D(P_11)^+, H mod 256
  switch_class.py h2 EXACT [members]        N = 20: class of the non-Paley doubly regular tournament of order 19
  switch_class.py general                   control: H mod 4 is not a switching invariant in general (orders 4-9)

EXACT / LAYERED = compiled hampaths_exact.c / hampaths_layered.c."""
import sys, random
from collections import Counter
from paley_tools import paley_out
from cayley_dickson import tournament
from drt_zoo import exact_H, layered_H, gs_solutions, goethals_seidel, descendant, double, paley_matrix, is_drt, find_map

def plus(out):
    n = len(out); return [((1 << (n + 1)) - 2)] + [(x << 1) for x in out]      # vertex 0 = the new source
def switch(out, X):
    n = len(out); o2 = [0] * n
    for i in range(n):
        for j in range(n):
            if i == j: continue
            a = (out[i] >> j) & 1
            if ((X >> i) & 1) != ((X >> j) & 1): a ^= 1
            if a: o2[i] |= 1 << j
    return o2
def is_sh(out):
    n = len(out); S = [[0 if i == j else (1 if (out[i] >> j) & 1 else -1) for j in range(n)] for i in range(n)]
    return all(sum(S[i][k] * S[k][j] for k in range(n)) == (-(n - 1) if i == j else 0) for i in range(n) for j in range(n))
def report(name, hs, modulus_note=""):
    print("%s: %d members; H mod 4: %s; H mod 8: %s%s" % (name, len(hs), sorted(Counter(h % 4 for h in hs).items()),
          sorted(Counter(h % 8 for h in hs).items()), modulus_note)); sys.stdout.flush()

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "insample":
        exact = sys.argv[2]; random.seed(5)
        for name, out in [("P_3", paley_out(3)), ("P_7", paley_out(7)), ("P_11", paley_out(11)), ("CD_4", tournament(4)[0]), ("P_19", paley_out(19))]:
            T = plus(out); N = len(T); assert is_sh(T)
            Xs = list(range(1 << (N - 1))) if N <= 8 else [random.randrange(1 << (N - 1)) for _ in range(60 if N <= 16 else 30)]
            tours = [switch(T, X) for X in Xs]; assert all(is_sh(u) for u in tours[:5])
            report("class of %s^+ (N = %d)" % (name, N), exact_H(exact, tours))
    elif mode == "h1":
        layered = sys.argv[2]; k = int(sys.argv[3]) if len(sys.argv) > 3 else 40; random.seed(24)
        K = double(paley_matrix(11)); D = descendant(K, 0); assert is_drt(D)
        ok = True
        for name, out in [("P_23", paley_out(23)), ("D(P_11)", D)]:
            T = plus(out); N = len(T); assert is_sh(T)
            Xs = [random.randrange(1 << (N - 1)) for _ in range(k)]
            hs = []
            for X in Xs:
                u = switch(T, X); hs.append(layered_H(layered, u))
            report("class of %s^+ (N = %d)" % (name, N), hs, " (H known mod 256)")
            ok = ok and all(h % 4 == 1 for h in hs)
        print("R11 (H1: H = 1 mod 4): %s" % ("HOLDS" if ok else "FAILS"))
    elif mode == "h2":
        exact = sys.argv[2]; k = int(sys.argv[3]) if len(sys.argv) > 3 else 40; random.seed(20)
        a, b, c, d = next(gs_solutions(5)); M = goethals_seidel(a, b, c, d); G = descendant(M, 0); assert is_drt(G)
        assert find_map(G, paley_out(19)) is None
        T = plus(G); N = len(T); assert is_sh(T)
        tours = [switch(T, random.randrange(1 << (N - 1))) for _ in range(k)]
        hs = exact_H(exact, tours); report("class of (non-Paley doubly regular tournament of order 19)^+ (N = 20)", hs)
        print("R11 (H2: H = 3 mod 4): %s" % ("HOLDS" if all(h % 4 == 3 for h in hs) else "FAILS"))
    elif mode == "general":
        from fold import ham
        random.seed(2)
        for n in range(4, 10):
            same = diff = 0
            for _ in range(200 if n <= 8 else 60):
                out = [0] * n
                for i in range(n):
                    for j in range(i + 1, n):
                        if random.random() < 0.5: out[i] |= 1 << j
                        else: out[j] |= 1 << i
                X = random.randrange(1, (1 << n) - 1)
                if ham(out, n) % 4 == ham(switch(out, X), n) % 4: same += 1
                else: diff += 1
            print("order %d, random tournament and random switching set: H mod 4 unchanged in %d cases, changed in %d" % (n, same, diff))
