"""Finite reachability and per-visit merge probabilities for the T-clock relation chain (lane 'groups').

Box B(J, C) = {(j, N) : |j| <= J, |N| <= C 3^|j|}  (N = c 3^max(0,-j); for j >= 0 this is |c/3^j| <= C, for j < 0 it
is |c| <= C, i.e. the reversed relation's normalized translation).  Part A (FINITE-EXACT): for every state of the box,
the exact probability of reaching (0,0) within m fair-coin steps (dynamic programming over all 2^m words, with
identical states merged).  Part B (NUMERICAL): along unmerged paths of the M1 chain (x = y + 1), the frequency of
visits to the box, the law of c at visits to j = +-1, and the hit rate of the merge precursors (1,1), (-1,-1/3).
Usage: python3 tclock_box.py [J C m]   defaults 1 10 24
"""
import sys, random, math
from collections import defaultdict, Counter
from tclock_chain import step

args = [int(a) for a in sys.argv[1:4]]
J, C, M = (args + [1, 10, 24][len(args):])[:3]


def box_states(J, C):
    out = []
    for j in range(-J, J + 1):
        R = C * 3 ** abs(j)
        for N in range(-R, R + 1):
            out.append((j, N))
    return out


def absorb_probs(J1, C1, J2, C2):
    """FINITE-EXACT: for every state of B(J1,C1), the exact probability that the chain reaches (0,0) before leaving
    the larger box B(J2,C2) (absorbing chain on a finite state space; sparse linear solve in exact float64,
    cross-checked by value iteration).  This is a rigorous lower bound for the merge probability."""
    import numpy as np
    from scipy.sparse import lil_matrix, identity
    from scipy.sparse.linalg import spsolve
    states = box_states(J2, C2)
    idx = {s: i for i, s in enumerate(states)}
    n = len(states)
    A = lil_matrix((n, n))
    b = np.zeros(n)
    for s, i in idx.items():
        if s == (0, 0):
            continue
        for p in (0, 1):
            t = step(s[0], s[1], p)
            if t == (0, 0):
                b[i] += 0.5
            elif t in idx:
                A[i, idx[t]] += 0.5
            # else: exits the box (counted as failure)
    h = spsolve((identity(n) - A.tocsr()).tocsc(), b)
    # value iteration cross-check (monotone from 0)
    Ac = A.tocsr(); v = np.zeros(n)
    for _ in range(4000):
        v = Ac @ v + b
    small = box_states(J1, C1)
    res = [(h[idx[s]], s) for s in small if s != (0, 0)]
    return res, float(np.max(np.abs(v - h))), n


def partA():
    J2, C2 = J + 2, 3 * C
    print('Part A (FINITE-EXACT): exact P(merge before leaving B(%d,%d)) from every state of B(%d,%d)' % (J2, C2, J, C))
    res, err, n = absorb_probs(J, C, J2, C2)
    res.sort()
    print('  finite chain: %d states; value-iteration vs direct solve max difference %.2e' % (n, err))
    print('  min over B(%d,%d) = %.5f at %s;  five smallest: %s' % (J, C, res[0][0], res[0][1],
          ', '.join('%s: %.4f' % (s, p) for p, s in res[:5])))
    print('  mean over the box = %.4f;  states with probability 0: %d' % (sum(p for p, _ in res) / len(res),
          sum(1 for p, _ in res if p == 0)))
    for st in [(0, 1), (1, 2), (1, 1), (0, -1), (1, 0), (2, 1)]:
        hit = [p for p, s in res if s == st]
        if hit:
            print('    from %s: %.5f' % (st, hit[0]))
    return res


def partB(seed=11, npaths=3000, smax=30000):
    print('\nPart B (NUMERICAL): M1 chain (x = y + 1); visits of unmerged paths to j = +-1, the law of c there,')
    print('  and how often the precursor states (1,1) and (-1,-1/3) are hit (each hit merges with probability 1/2)')
    rng = random.Random(seed)
    visits1 = Counter(); visitsm1 = Counter()
    nvis = 0; nbox = 0; steps_alive = 0; prec_hits = 0
    for _ in range(npaths):
        j, N = 0, 1
        for s in range(smax):
            if j == 0 and N == 0:
                break
            steps_alive += 1
            if abs(j) <= J and abs(N) <= C * 3 ** abs(j):
                nbox += 1
            if j == 1:
                visits1[N] += 1; nvis += 1
            elif j == -1:
                visitsm1[N] += 1; nvis += 1
            if (j, N) in ((1, 1), (-1, -1)):
                prec_hits += 1
            j, N = step(j, N, rng.getrandbits(1))
    print('  alive steps %d; visits to j = +-1: %d (%.4f per alive step); box visits %.4f per alive step'
          % (steps_alive, nvis, nvis / steps_alive, nbox / steps_alive))
    print('  precursor hits per visit to j = +-1: %.4f' % (prec_hits / nvis))
    print('  most frequent c at j = 1 :', ', '.join('%d:%.3f' % (c, k / sum(visits1.values())) for c, k in visits1.most_common(8)))
    print('  most frequent N at j = -1 (c = N/3):', ', '.join('%d:%.3f' % (c, k / sum(visitsm1.values())) for c, k in visitsm1.most_common(8)))
    big = sum(k for c, k in visits1.items() if abs(c) > 3 * C) / max(1, sum(visits1.values()))
    print('  share of j = 1 visits with |c| > %d: %.4f' % (3 * C, big))
    # tail of |c|/3 at j = 1 visits
    tot = sum(visits1.values())
    for u in (1, 3, 10, 30, 100, 300, 1000):
        sh = sum(k for c, k in visits1.items() if abs(c) / 3 > u) / tot
        print('    P(|c/3| > %5d | j = 1) = %.5f   u * P = %.3f' % (u, sh, u * sh))


if __name__ == '__main__':
    partA()
    partB()
