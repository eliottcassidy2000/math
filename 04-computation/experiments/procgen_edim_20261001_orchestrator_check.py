#!/usr/bin/env python3
"""procgen_edim_20261001_orchestrator_check.py -- the orchestrator's independent audit of the edim lane (2026-10-01).

Written from the definition (Allikvere, arXiv:2608.09983): for an edge uv of Q_d and a vertex s,
d(uv, s) = min(d(u,s), d(v,s)); S is edge-multiset resolving iff the histograms H_e(r) = #{s in S: d(e,s) = r}
are pairwise distinct over all edges. The lane's code was not read; only DATA was taken from it (the explicit
Q7..Q12 sets, parsed from the CERT literal of procgen_edim_20261001_run.py, and the deposited orbit list).

Part I  (C engine procgen_edim_20261001_orchestrator_check.c, a THIRD exhaustive search with its own reduction:
         split by bit 5 only, |A| >= |B|, A an Aut(Q5)-orbit minimum, every B enumerated; no imbalance normal form):
         no resolving set of size k <= 14; the resolving 15-sets found reduce to the deposited 229 orbits.
Part II (Python): the 229 orbit representatives are resolving, are orbit minima, are pairwise inequivalent,
         have trivial stabilizer; the paper's set lies in one of the orbits; the k = 14 near miss has defect 1;
         the explicit Q7..Q12 sets are resolving; the counting bound (L4) and the entropy bound table (L5).
usage: python3 procgen_edim_20261001_orchestrator_check.py [--skip-search]
"""
import ast
import itertools
import math
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', '..', '05-knowledge', 'results')
OKS = []


def ok(cond, msg):
    OKS.append(bool(cond))
    print(('[OK] ' if cond else '[FAIL] ') + msg, flush=True)


def edges(d):
    return [(u, u ^ (1 << i)) for u in range(1 << d) for i in range(d) if not u >> i & 1]


def hist_matrix(d, S):
    """rows = edges, columns = r = 0..d-1, entries H_e(r)"""
    E = np.array(edges(d), dtype=np.int64)
    S = np.array(sorted(S), dtype=np.int64)
    pc = np.vectorize(lambda x: bin(int(x)).count('1'))
    du = pc(E[:, 0:1] ^ S[None, :])
    dv = pc(E[:, 1:2] ^ S[None, :])
    dm = np.minimum(du, dv)
    H = np.zeros((len(E), d), dtype=np.int64)
    for r in range(d):
        H[:, r] = (dm == r).sum(axis=1)
    return H


def n_distinct(H):
    return len({tuple(row) for row in H.tolist()})


def mask_to_set(m):
    return [v for v in range(64) if m >> v & 1]


def aut_q6():
    """all 46080 automorphisms of Q6 as a (46080, 64) array of vertex images"""
    perms = []
    for pi in itertools.permutations(range(6)):
        base = np.array([sum(1 << pi[i] for i in range(6) if x >> i & 1) for x in range(64)], dtype=np.int64)
        for t in range(64):
            perms.append(base ^ t)
    return np.array(perms, dtype=np.int64)


def canon(P, S):
    imgs = P[:, S]                                       # (46080, |S|)
    masks = np.zeros(P.shape[0], dtype=np.uint64)
    for j in range(imgs.shape[1]):
        masks |= (np.uint64(1) << imgs[:, j].astype(np.uint64))
    return int(masks.min()), int((masks == masks.min()).sum())


def search_part():
    with tempfile.TemporaryDirectory() as tdir:
        b = os.path.join(tdir, 'ea')
        subprocess.run(['cc', '-O3', '-o', b, os.path.join(HERE, 'procgen_edim_20261001_orchestrator_check.c')], check=True)
        found = os.path.join(tdir, 'found.txt')
        t0 = time.time()
        r = subprocess.run(['nice', b, '1', '15', found], capture_output=True, text=True, check=True)
        print(r.stdout.strip(), flush=True)
        print(f'  (search time {time.time() - t0:.0f} s)', flush=True)
        ok('MATCH Burnside' in r.stdout, 'Aut(Q5)-orbit representatives of a-subsets (a <= 15) match the Burnside counts')
        for k in range(1, 16):
            m = re.search(rf'^k={k} expected leaves=(\d+)$', r.stdout, re.M)
            m2 = re.search(rf'^k={k} leaves=(\d+) fullchecks=\d+ resolving_found=(\d+)$', r.stdout, re.M)
            ok(m and m2 and m.group(1) == m2.group(1), f'k={k}: enumeration complete (leaves = sum_a #reps(a) C(32, k-a) = {m.group(1) if m else "?"})')
            if k <= 14:
                ok(m2 and m2.group(2) == '0', f'k={k}: no edge-multiset resolving set')
            else:
                ok(m2 and int(m2.group(2)) > 0, f'k=15: {m2.group(2) if m2 else "?"} resolving leaves found')
        with open(found) as f:
            masks = [int(line.split()[1], 16) for line in f if line.startswith('k=15')]
        return masks


def main():
    t0 = time.time()
    P = aut_q6()
    ok(P.shape == (46080, 64) and all(len(set(row)) == 64 for row in P[::997].tolist()), '|Aut(Q6)| = 46080 vertex permutations')
    # every automorphism preserves Hamming distance (spot check)
    pc = lambda x: bin(int(x)).count('1')
    ok(all(pc(P[g, x] ^ P[g, y]) == pc(x ^ y) for g in range(0, 46080, 4093) for x in range(64) for y in range(0, 64, 7)), 'automorphisms preserve distance')

    found = None
    if '--skip-search' not in sys.argv:
        print('==== Part I: third exhaustive search (orchestrator C engine) ====', flush=True)
        found = search_part()

    print('==== Part II: certificates, orbits, bounds ====', flush=True)
    with open(os.path.join(RES, 'procgen_edim_20261001_q6_resolving15_orbits.txt')) as f:
        reps = [int(l.strip(), 16) for l in f if l.strip() and not l.startswith('#')]
    ok(len(reps) == 229, f'deposited orbit list has {len(reps)} representatives')
    allres = all(n_distinct(hist_matrix(6, mask_to_set(m))) == 192 and bin(m).count('1') == 15 for m in reps)
    ok(allres, 'all 229 representatives are 15-sets that resolve all 192 edges (own checker)')
    cans = [canon(P, mask_to_set(m)) for m in reps]
    ok(all(c == m for (c, _), m in zip(cans, reps)), 'each representative is the minimum mask of its Aut(Q6)-orbit')
    ok(len({c for c, _ in cans}) == 229, 'the 229 orbits are pairwise distinct')
    ok(all(st == 1 for _, st in cans), 'every orbit has trivial stabilizer (the minimum is attained by exactly one automorphism)')
    paper = 0x02283022a042a00a
    ok(sorted(mask_to_set(paper)) == [1, 3, 13, 15, 17, 22, 29, 31, 33, 37, 44, 45, 51, 53, 57], "paper's set decoded")
    ok(n_distinct(hist_matrix(6, mask_to_set(paper))) == 192, "paper's 15-set (Table 1) is resolving")
    ok(canon(P, mask_to_set(paper))[0] in set(reps), "paper's set lies in one of the 229 orbits")
    if found is not None:
        fc = {canon(P, mask_to_set(m))[0] for m in found}
        ok(fc == set(reps), f'the {len(found)} resolving 15-sets found by the third search reduce to exactly the deposited 229 orbits ({len(fc)} orbits)')
    # near miss at k = 14
    nm = 0x000001810690226d
    H = hist_matrix(6, mask_to_set(nm))
    E = edges(6)
    groups = {}
    for i, row in enumerate(H.tolist()):
        groups.setdefault(tuple(row), []).append(E[i])
    coll = [(k, v) for k, v in groups.items() if len(v) > 1]
    ok(bin(nm).count('1') == 14 and n_distinct(H) == 191 and len(coll) == 1 and coll[0][0] == (1, 1, 5, 5, 1, 1)
       and sorted(coll[0][1]) == [(24, 26), (37, 39)], 'k=14 near miss 0x000001810690226d: defect 1, the one collision is {24,26} ~ {37,39} with histogram (1,1,5,5,1,1)')
    # explicit Q7..Q12 sets (data parsed from the lane's runner, checked by this code)
    src = open(os.path.join(HERE, 'procgen_edim_20261001_run.py')).read()
    i = src.index('CERT = {')
    j = src.index('}', i)
    cert = ast.literal_eval(src[i + len('CERT = '):j + 1])
    for d in sorted(cert):
        S = cert[d]
        ne = d * 2 ** (d - 1)
        good = len(set(S)) == len(S) and all(0 <= s < 2 ** d for s in S) and n_distinct(hist_matrix(d, S)) == ne
        ok(good, f'Q{d}: explicit set of size {len(S)} resolves all {ne} edges  =>  edim_m(Q{d}) <= {len(S)}')
    # L4 counting bound for d = 6
    ok(all(192 - 12 * m > math.comb(m + 3, 3) for m in range(1, 7)) and not (192 - 12 * 7 > math.comb(10, 3)),
       'L4: 192 - 12m > C(m+3,3) exactly for m <= 6, so edim_m(Q6) >= 7')
    # L5 entropy bound table: least m with log2(d 2^(d-1)) <= sum_{r=0}^{d-2} g(m C(d-1,r)/2^(d-1))
    def g(mu):
        return 0.0 if mu <= 0 else (mu + 1) * math.log2(mu + 1) - mu * math.log2(mu)
    def least_m(d):
        lhs = math.log2(d) + (d - 1)
        # p_r = C(d-1, r) / 2^(d-1) in floating point via lgamma (the right side is increasing in m: binary search)
        pr = [math.exp(math.lgamma(d) - math.lgamma(r + 1) - math.lgamma(d - r) - (d - 1) * math.log(2)) for r in range(d - 1)]
        rhs = lambda m: sum(g(m * p) for p in pr)
        lo, hi = 1, 2
        while rhs(hi) < lhs:
            hi *= 2
        while lo < hi:
            mid = (lo + hi) // 2
            if rhs(mid) >= lhs:
                hi = mid
            else:
                lo = mid + 1
        return lo
    table = {6: 4, 7: 5, 8: 5, 10: 6, 12: 8, 16: 11, 20: 15, 32: 28, 64: 81, 128: 275, 256: 1159, 1024: 50116}
    got = {d: least_m(d) for d in table}
    ok(got == table, f'L5 entropy bound table reproduced: {got}')
    # L2 antipodal reversal and L1-type sanity on random sets
    rng = np.random.default_rng(20261001)
    good = True
    for _ in range(20):
        S = sorted(rng.choice(64, size=15, replace=False).tolist())
        Hm = hist_matrix(6, S)
        idx = {e: i for i, e in enumerate(E)}
        for i, (u, v) in enumerate(E):
            a, b = 63 ^ u, 63 ^ v
            jj = idx[(min(a, b), max(a, b))]
            if list(Hm[jj]) != list(Hm[i][::-1]):
                good = False
    ok(good, 'L2: H_{antipode(e)}(r) = H_e(5 - r) on 20 random 15-sets')
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')


if __name__ == '__main__':
    main()
