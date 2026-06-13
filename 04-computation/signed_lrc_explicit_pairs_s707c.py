"""Extract EXPLICIT colliding cut-pairs + verify the silent-condition formula.
monad-explorer-2026-06-06-S707.

Silent condition (derived): flipping H_q = {C/q, 2C/q, .., ((q-1)/2)C/q} takes cut eps to
eps' = eps with H_q negated.  At freq t with q|t, Phi unchanged.  At q-nmid-t the move is
silent iff Phi_{eps'}(t) = +- Phi_eps(t).  We:
  (1) print explicit colliding pairs for small composite C (runner sign strings),
  (2) for each colliding pair, show the per-frequency sign pattern s_t = Phi_eps'(t)/Phi_eps(t),
  (3) confirm s_t in {+1,-1} and s_t = +1 exactly when q|t (CHECK the clean reflection claim).
"""
import math
from itertools import product
from collections import defaultdict


def folded_spectrum(V, eps, C):
    out = []
    n1 = len(V)
    for i in range(n1):
        for j in range(i + 1, n1):
            f = (eps[i] * V[i] - eps[j] * V[j]) % C
            out.append(min(f, C - f))
    return tuple(sorted(out))


def Phi(V, eps, C, t):
    return sum(eps[i] * math.sin(2 * math.pi * t * V[i] / C) for i in range(len(V)))


def is_prime(C):
    return C > 1 and all(C % d for d in range(2, int(C**0.5) + 1))


def epsstr(eps):
    return "".join("+" if e == 1 else "-" for e in eps)


for n in range(5, 15):
    C = 2 * n - 1
    if is_prime(C):
        continue
    V = list(range(1, n))
    n1 = len(V)
    fold_map = defaultdict(list)
    for bits in product([0, 1], repeat=n1 - 1):
        eps = (1,) + tuple(1 if b else -1 for b in bits)
        fold_map[folded_spectrum(V, eps, C)].append(eps)
    nontrivial = [g for g in fold_map.values() if len(g) > 1]
    print(f"\n===== C={C} (n={n}) : {len(nontrivial)} nontrivial classes =====")
    shown = 0
    for g in nontrivial:
        if shown >= 4:
            print(f"   ... ({len(nontrivial)} classes total)")
            break
        # show class
        print(f"  class size {len(g)}:")
        base = g[0]
        for other in g:
            T = [V[i] for i in range(n1) if base[i] != other[i]]
            print(f"     {epsstr(other)}   (flip vs first: runners {T})")
        # per-frequency sign for the first nontrivial pair
        e0, e1 = g[0], g[1]
        Tset = set(V[i] for i in range(n1) if e0[i] != e1[i])
        q_candidates = sorted(set(C // math.gcd(C, x) for x in Tset))
        signs = []
        for t in range(1, n1 + 1):
            a, b = Phi(V, e0, C, t), Phi(V, e1, C, t)
            if abs(a) < 1e-9 and abs(b) < 1e-9:
                signs.append(0)
            elif abs(a) < 1e-9 or abs(b) < 1e-9:
                signs.append(9)  # one zero one not -> would break collision (shouldn't happen)
            else:
                signs.append(1 if (b / a) > 0 else -1)
        # compare with q | t pattern for each candidate q
        print(f"     flip-set runners={sorted(Tset)} subgroup-order q in {q_candidates}")
        print(f"     per-freq sign s_t (t=1..{n1}): {signs}")
        for q in q_candidates:
            qpat = [1 if (t % q == 0) else -1 for t in range(1, n1 + 1)]
            match = all((s == 0) or (s == qp) for s, qp in zip(signs, qpat))
            print(f"        q={q}: s_t == (+1 iff q|t) ? {match}   [+1-iff-q|t pattern={qpat}]")
        shown += 1
