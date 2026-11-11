#!/usr/bin/env python3
"""H mod 4 on a switching class is a quadratic function of the switching set (observation; section 5.3).

For a tournament U and a vertex set X let U_X be U switched at X, and f_U(X) = (H(U_X) - H(U))/2 mod 2.
Observation: f_U is a polynomial of degree at most 2 over F_2 in the indicator of X; equivalently, for any three
vertices u, v, w the eight tournaments U_S, S subset {u, v, w}, satisfy  sum_S H(U_S) = 0 (mod 4).

  switch_quadratic.py anf                  algebraic degree of f_U: all labelled tournaments of order 4, 5, 6;
                                           random tournaments of order 7..10
  switch_quadratic.py third EXACT          sum over the 8 switchings at subsets of a random triple, mod 4: random
                                           tournaments of order 11..20 and members of skew Hadamard classes (N = 12,
                                           16, 20)
  switch_quadratic.py third24 LAYERED JSON the same at N = 24 (P_23, the double of P_11, three annealed classes)

EXACT / LAYERED = compiled hampaths_exact.c / hampaths_layered.c."""
import sys, random, json
from fold import ham
from local_forms import rand_tournament, all_tournaments
from switch_class import plus, switch

def anf_degree(f, k):
    g = f[:]
    for i in range(k):
        for m in range(1 << k):
            if (m >> i) & 1: g[m] ^= g[m ^ (1 << i)]
    return max((bin(m).count('1') for m in range(1 << k) if g[m]), default=-1)
def degree_of(out):
    n = len(out); k = n - 1; h0 = ham(out, n)          # the last vertex is never switched (X and its complement agree)
    return anf_degree([((ham(switch(out, X), n) - h0) // 2) % 2 for X in range(1 << k)], k)

if __name__ == "__main__":
    mode = sys.argv[1]; rng = random.Random(20261111)
    if mode == "anf":
        for n in (4, 5, 6):
            degs = {}
            for o in all_tournaments(n):
                d = degree_of(o); degs[d] = degs.get(d, 0) + 1
            print("order %d, all %d labelled tournaments: algebraic degree of f_U: %s  (-1 = identically zero)" % (n, sum(degs.values()), sorted(degs.items())))
            sys.stdout.flush()
        for n, tr in ((7, 60), (8, 40), (9, 12), (10, 4)):
            degs = {}
            for _ in range(tr):
                d = degree_of(rand_tournament(n, rng)); degs[d] = degs.get(d, 0) + 1
            print("order %d, %d random tournaments: algebraic degree of f_U: %s" % (n, tr, sorted(degs.items()))); sys.stdout.flush()
    elif mode == "third":
        from drt_zoo import exact_H, gs_solutions, goethals_seidel, descendant
        from paley_tools import paley_out
        from cayley_dickson import tournament
        exe = sys.argv[2]
        def third(U, trials):
            N = len(U); tours = []
            for _ in range(trials):
                X = rng.getrandbits(N); V = switch(U, X); u, v, w = rng.sample(range(N), 3)
                tours += [switch(V, sum((1 << x) for x, b in zip((u, v, w), (S & 1, S >> 1 & 1, S >> 2 & 1)) if b)) for S in range(8)]
            hs = exact_H(exe, tours)
            return sum(1 for i in range(0, len(hs), 8) if sum(hs[i:i + 8]) % 4)
        for n in (11, 12, 13, 14, 15, 16, 18, 20):
            tr = 40 if n <= 16 else 10; bad = 0
            for _ in range(tr): bad += third(rand_tournament(n, rng), 1)
            print("order %d: %d random (tournament, triple): sum of H over the 8 switchings != 0 (mod 4) in %d" % (n, tr, bad)); sys.stdout.flush()
        a, b, c, d = next(gs_solutions(5)); D19 = descendant(goethals_seidel(a, b, c, d), 0)
        for name, T in [("P_11", paley_out(11)), ("CD_4", tournament(4)[0]), ("P_19", paley_out(19)), ("D_19", D19)]:
            U = plus(T); tr = 30 if len(U) <= 16 else 10
            print("class of (source) + %s (N = %d): %d random (member, triple): violations %d" % (name, len(U), tr, third(U, tr))); sys.stdout.flush()
    elif mode == "third24":
        from drt_zoo import layered_H, descendant, double, paley_matrix
        from paley_tools import paley_out
        from cycle_parities import code_to_out
        lay = sys.argv[2]; d = json.load(open(sys.argv[3]))
        rows = [("P_23", paley_out(23)), ("double of P_11", descendant(double(paley_matrix(11)), 0))] + [(rec['label'], code_to_out(rec['code'], 23)) for rec in d[:3]]
        for name, T in rows:
            U = plus(T); N = 24; bad = 0; tr = 6
            for _ in range(tr):
                V = switch(U, rng.getrandbits(N)); u, v, w = rng.sample(range(N), 3)
                tot = sum(layered_H(lay, switch(V, sum((1 << x) for x, b in zip((u, v, w), (S & 1, S >> 1 & 1, S >> 2 & 1)) if b))) for S in range(8))
                bad += (tot % 4 != 0)
            print("class of (source) + %s (N = 24): %d random (member, triple): sum of H over the 8 switchings != 0 (mod 4) in %d" % (name, tr, bad)); sys.stdout.flush()
