"""
metagraph_anti_redei_test_opus_S139.py

TEST OF THE FIBER LAW + THE ANTI-REDEI CONJECTURE (opus-S139, companion to the fiber census).

PROVED (2 lines, in the reflection): a tiling t is grid-symmetric iff rho(i) = n+1-i is an
ANTI-AUTOMORPHISM of T_t (the base path is rho-anti-symmetric; on tiles the condition is
bit(sigma tile) = bit(tile)).

CONJECTURE (C1, the fiber law): for every class C,
    g(C) * |Aut(C)| = H_anti(C) := #{ Hamiltonian paths P = (p_1..p_n) of T :
                                       the map p_j -> p_{n+1-j} is an anti-automorphism of T }
(the exact anti-symmetric analog of LEM-003's N(C) * |Aut| = H).

CONSEQUENCES if true:
  * N(C) odd  = Redei (H odd) + odd |Aut|                     [PARITY LAW half 1 — proved]
  * g(C) odd on transpose-self classes  <=  H_anti odd there  [ANTI-REDEI — new parity thm]
  * pure-black = non-transpose-self (no anti-auto at all)     [proved by the 2-line lemma]
  * mixed vs pure-blue: b(C) = (H - H_anti)/|Aut| = 0 iff every Ham path is anti-symmetric.

This script verifies (C1) exactly at n = 4, 5, 6 (every class: enumerate Ham paths of a
representative, count anti-symmetric ones, compare with g*|Aut|), and verifies H_anti
oddness on all transpose-self classes.
"""
import sys, time
from itertools import permutations
from collections import Counter

sys.path.insert(0, ".")
from metagraph_fiber_allocation_opus_S139 import tiles_of, sigma_map, arcs_of_tiling, canon

def census_with_reps(n):
    T = tiles_of(n); m = len(T)
    sig = sigma_map(n, T)
    perms = list(permutations(range(n)))
    cls_of = {}; rep = {}
    for bits in range(1 << m):
        a = arcs_of_tiling(n, T, bits)
        c = canon(a, n, perms)
        cls_of[bits] = c
        if c not in rep: rep[c] = a
    def gridsym(bits):
        return all(((bits >> i) & 1) == ((bits >> s) & 1) for i, s in enumerate(sig))
    from collections import Counter as Ctr
    N = Ctr(cls_of.values())
    g = Ctr(cls_of[b] for b in range(1 << m) if gridsym(b))
    return T, cls_of, rep, N, g, perms

def ham_paths(a, n):
    """All Hamiltonian paths of tournament a (as vertex sequences)."""
    paths = []
    def ext(seq, used):
        if len(seq) == n:
            paths.append(tuple(seq)); return
        for v in range(n):
            if not used[v] and a[seq[-1]][v]:
                used[v] = True; seq.append(v)
                ext(seq, used)
                seq.pop(); used[v] = False
    for s in range(n):
        used = [False] * n; used[s] = True
        ext([s], used)
    return paths

def is_antiauto(a, n, f):
    return all(a[f[u]][f[v]] == a[v][u] for u in range(n) for v in range(n) if u != v)

def aut_order(a, n, perms):
    return sum(1 for p in perms if all(a[p[u]][p[v]] == a[u][v]
               for u in range(n) for v in range(n) if u != v))

def main():
    for n in (4, 5, 6):
        t0 = time.time()
        T, cls_of, rep, N, g, perms = census_with_reps(n)
        ok = bad = 0
        antired_ok = True
        for c, a in rep.items():
            H = ham_paths(a, n)
            Hanti = 0
            for P in H:
                f = [0] * n
                for j, v in enumerate(P):
                    f[v] = P[n - 1 - j]
                if is_antiauto(a, n, f):
                    Hanti += 1
            au = aut_order(a, n, perms)
            gc = g.get(c, 0)
            if gc * au == Hanti:
                ok += 1
            else:
                bad += 1
                print(f"   n={n} class {c}: g*|Aut| = {gc}*{au} = {gc*au}  !=  H_anti = {Hanti}")
            if gc > 0 and Hanti % 2 == 0:
                antired_ok = False
                print(f"   n={n} class {c}: H_anti = {Hanti} EVEN on a transpose-self class!")
            # sanity: N*|Aut| = H (LEM-003)
            assert N[c] * au == len(H), (n, c, N[c], au, len(H))
        print(f"n={n}: (C1) g*|Aut| == H_anti on {ok}/{ok+bad} classes"
              f"{'  *** FAILURES ***' if bad else ''};  LEM-003 sanity all pass;"
              f"  anti-Redei (H_anti odd on selfop) {'HOLDS' if antired_ok else 'FAILS'}"
              f"  [{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
