#!/usr/bin/env python3
"""
procgen_drift_20260926_ihs.py -- exact flip distance of q n + 1 to class (i) by an implicit hitting set (IHS),
drift lane 2026-09-26 (session collatz-procgen-20260922).

delta_k(q) = least number of odd residues mod 2^k at which a class-(i) strategy differs from q n + 1.

Method (exact):
  * y_v in {0,1} for every odd residue v (y_v = 1: sigma(v) = -1).  Objective: min sum y_v.
  * NO-GOODS.  A cycle C of G_sigma is a cycle of G_sigma' whenever sigma' agrees with sigma on the odd nodes of
    C (edges out of even nodes never depend on sigma).  So an expanding cycle C of G_sigma gives the valid clause
    "sigma' differs from sigma at some odd node of C".  Each no-good is stored as (cycle, minus-signed odd nodes of
    the cycle) and is re-checkable on its own (closed walk under those signs, q^a > 2^p).
  * The optimum of the hitting-set problem over the no-goods found so far is a LOWER bound; a hitting set whose
    G_sigma has no expanding cycle is class (i), so the first such optimum is optimal.  Any class-(i) set gives an
    UPPER bound; when the solver proves "no hitting set of size <= UB-1", UB is optimal.
  * Seeds: every expanding cycle of the all-plus graph (de Bruijn B(2,k)) of length <= P.
  * Each round: solve the hitting set exactly (core-guided MaxSAT RC2, single thread), harvest node-disjoint expanding cycles of the
    optimum (C engine), and run greedy repair + pruning phases (more no-goods, better UB).
  * Every reported optimum is re-verified class (i) by an edge-checked potential and (k <= 13) by exact Karp.
"""
import os
import sys
import json
import gzip
import time
import random
import itertools

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_drift_20260926_lib as L  # noqa: E402


def base_expanding_cycles(k, q, P):
    """all expanding simple cycles of the all-plus graph G_0 of length p <= P: for each primitive necklace w of
    length p with q^(ones) > 2^p (and all k-windows of w^inf distinct), the node list of residues whose parity
    k-word is the window of w^inf at positions 0..p-1."""
    table = {}
    for r in range(1 << k):
        table[tuple(L.parity_word(r, k, q))] = r
    out = []
    seen = set()
    for p in range(1, P + 1):
        for ones in range(p + 1):
            if not q ** ones > 2 ** p:
                continue
            for pos in itertools.combinations(range(p), ones):
                w = [0] * p
                for i in pos:
                    w[i] = 1
                rots = [tuple(w[i:] + w[:i]) for i in range(p)]
                if tuple(w) != min(rots) or len(set(rots)) != p:
                    continue
                ext = w * ((k // p) + 2)
                cyc = [table[tuple(ext[i:i + k])] for i in range(p)]
                if len(set(cyc)) != p:
                    continue
                key = tuple(sorted(set(v for v in cyc if v & 1)))
                if key not in seen:
                    seen.add(key)
                    out.append(cyc)
    return out


def nogood_ok(k, q, cyc, minus_odd):
    """the stored no-good is genuine: cyc is a closed walk of G_sigma for every sigma with sigma = - exactly on
    minus_odd among the odd nodes of cyc, and q^a > 2^p."""
    N = 1 << k
    H = N >> 1
    p = len(cyc)
    for i in range(p):
        s, t = cyc[i], cyc[(i + 1) % p]
        if not (0 <= s < N and 0 <= t < N):
            return False
        if s & 1:
            u = (q * s + (-1 if s in minus_odd else 1)) // 2
        else:
            u = s // 2
        if t % H != u % H:
            return False
    a = sum(1 for v in cyc if v & 1)
    return q ** a > 2 ** p


class IHS:
    def __init__(self, k, q, P=None, seed=0, log=print):
        self.k, self.q = k, q
        self.N = 1 << k
        self.odd = list(range(1, self.N, 2))
        self.F = L.best_lower_approx(self.N, q)
        self.nogoods = []
        self.clauses = []
        self.seen = set()
        self.counts = np.zeros(self.N, dtype=np.int64)
        self.rng = random.Random(seed)
        self.log = log
        P = P if P is not None else k + 2
        for cyc in base_expanding_cycles(k, q, P):
            self.add(cyc, frozenset())
        self.nseed = len(self.clauses)
        self.UB = None
        self.best = None
        self.LB = 0

    def add(self, cyc, minus_odd):
        lits = []
        for v in sorted(set(u for u in cyc if u & 1)):
            lits.append((v, 0 if v in minus_odd else 1))   # (v, value that breaks the no-good): y_v must be != recorded
        key = tuple(lits)
        if key in self.seen:
            return False
        self.seen.add(key)
        self.nogoods.append((list(cyc), frozenset(minus_odd)))
        self.clauses.append(lits)
        for v, _ in lits:
            self.counts[v] += 1
        return True

    # ------------------------------------------------------------------ hitting set (core-guided MaxSAT, RC2)
    def solve_hs(self, cutoff=None, hint=None, time_limit=3600.0):
        """exact min sum y subject to the no-good clauses (pysat RC2 with Glucose4; exact optimum).
        Returns ('optimal', set) or ('infeasible', None); with cutoff = c, returns ('cutoff', None) when the optimum
        exceeds c (so the true optimum is > c).  (time_limit is not enforced inside RC2.)"""
        from pysat.examples.rc2 import RC2
        from pysat.formula import WCNF
        w = WCNF()
        var = {v: i + 1 for i, v in enumerate(self.odd)}
        for lits in self.clauses:
            w.append([var[v] if val == 1 else -var[v] for v, val in lits])
        for v in self.odd:
            w.append([-var[v]], weight=1)
        with RC2(w, solver='g4', adapt=True, exhaust=True, minz=True) as rc2:
            m = rc2.compute()
            if m is None:
                return 'infeasible', None
            sol = set(v for v in self.odd if m[var[v] - 1] > 0)
            L.check(len(sol) == rc2.cost, "RC2 cost mismatch")
        # independent check that sol satisfies every clause
        for lits in self.clauses:
            L.check(any((v in sol) == (val == 1) for v, val in lits), "RC2 solution violates a no-good")
        if cutoff is not None and len(sol) > cutoff:
            return 'cutoff', None
        return 'optimal', sol

    # ------------------------------------------------------------------ harvesting
    def harvest(self, flipset, M=60):
        fl = L.flip_array(self.k, flipset)
        cycs = L.disjoint_cycles(self.k, self.q, fl, M, self.F)
        added = 0
        for cyc in cycs:
            mo = frozenset(v for v in cyc if (v & 1) and fl[v])
            L.check(nogood_ok(self.k, self.q, cyc, mo), "harvested no-good not genuine")
            added += self.add(cyc, mo)
        return cycs, added

    def repair(self, flipset, max_steps=60):
        cur = set(flipset)
        added = 0
        for _ in range(max_steps):
            cycs, a = self.harvest(cur)
            added += a
            if not cycs:
                return cur, added
            for cyc in cycs:
                cand = [v for v in cyc if (v & 1) and v not in cur]
                if cand:
                    v = max(cand, key=lambda u: (self.counts[u], self.rng.random()))
                    cur.add(v)
                else:
                    cur.discard(max((v for v in cyc if v & 1), key=lambda u: (-self.counts[u], self.rng.random())))
        return None, added

    def prune(self, flipset, rounds=2):
        """greedy: try to restore sigma = + at each flipped residue (random order), incremental engine"""
        fl = L.flip_array(self.k, flipset)
        E = L.Engine(self.k, self.q, fl, self.F)
        L.check(E.solve(), "prune start not class (i)")
        for _ in range(rounds):
            order = [v for v in self.odd if fl[v]]
            self.rng.shuffle(order)
            changed = False
            for v in order:
                if E.toggle(v):
                    changed = True
            fl = E.flip()
            if not changed:
                break
        return set(int(v) for v in np.nonzero(fl)[0])

    def verify_class_i(self, flipset):
        fl = L.flip_array(self.k, flipset)
        st = L.classify(self.k, self.q, fl, self.F)
        return st[0] == 'I'

    def offer_ub(self, flipset):
        if flipset is None:
            return False
        L.check(self.verify_class_i(flipset), "offered UB set not class (i)")
        if self.UB is None or len(flipset) < self.UB:
            self.UB = len(flipset)
            self.best = sorted(flipset)
            return True
        return False

    # ------------------------------------------------------------------ main loop (incremental RC2)
    def short_nogoods(self, flipset, P, cap=3000):
        fl = L.flip_array(self.k, flipset)
        Fp = L.best_lower_approx(max(P, 1), self.q)   # exact threshold for cycles of length <= P
        E = L.Engine(self.k, self.q, fl, Fp)
        added = 0
        for cyc in E.short_cycles(P, cap, needflip=True):
            mo = frozenset(v for v in cyc if (v & 1) and fl[v])
            L.check(nogood_ok(self.k, self.q, cyc, mo), "short cycle no-good not genuine")
            added += self.add(cyc, mo)
        return added

    def run(self, time_limit=3600.0, ckpt=None, verbose=True, shortP=None, repair_every=5, ub_hook=None,
            max_iter=None):
        from pysat.examples.rc2 import RC2
        from pysat.formula import WCNF
        t0 = time.time()
        var = {v: i + 1 for i, v in enumerate(self.odd)}
        w = WCNF()
        for lits in self.clauses:
            w.append([var[v] if val == 1 else -var[v] for v, val in lits])
        for v in self.odd:
            w.append([-var[v]], weight=1)
        nfed = len(self.clauses)
        shortP = shortP if shortP is not None else self.k + 4
        it = 0
        with RC2(w, solver='g4', adapt=True, exhaust=True, minz=True) as rc2:
            while True:
                it += 1
                if time.time() - t0 > time_limit or (max_iter is not None and it > max_iter):
                    return 'timeout'
                m = rc2.compute()
                if m is None:
                    return 'infeasible'
                sol = set(v for v in self.odd if m[var[v] - 1] > 0)
                for lits in self.clauses:
                    L.check(any((v in sol) == (val == 1) for v, val in lits), "RC2 solution violates a no-good")
                self.LB = max(self.LB, len(sol))
                if self.UB is not None and self.LB >= self.UB:
                    return 'optimal'
                added = self.short_nogoods(sol, shortP)
                cycs, a2 = self.harvest(sol)
                added += a2
                if not cycs and added == 0:
                    self.offer_ub(sol)
                    self.LB = len(sol)
                    return 'optimal'
                if ub_hook is not None and it % repair_every == 1:
                    cand = ub_hook(self, sol)
                    if cand is not None:
                        self.offer_ub(cand)
                for lits in self.clauses[nfed:]:
                    rc2.add_clause([var[v] if val == 1 else -var[v] for v, val in lits])
                nfed = len(self.clauses)
                if ckpt:
                    self.save(ckpt)
                if verbose:
                    self.log(f"    k={self.k} q={self.q} it {it}: LB {self.LB} UB {self.UB} no-goods {len(self.clauses)} "
                             f"(seed {self.nseed}, +{added}) {time.time() - t0:.0f}s")
                if self.UB is not None and self.LB >= self.UB:
                    return 'optimal'

    def save(self, path):
        data = {'k': self.k, 'q': self.q, 'LB': self.LB, 'UB': self.UB, 'best': self.best,
                'nogoods': [[c, sorted(m)] for c, m in self.nogoods[self.nseed:]]}
        with open(path + '.tmp', 'w') as f:
            json.dump(data, f)
        os.replace(path + '.tmp', path)

    def load(self, path):
        with open(path) as f:
            data = json.load(f)
        for c, m in data['nogoods']:
            L.check(nogood_ok(self.k, self.q, c, frozenset(m)), "checkpoint no-good fails")
            self.add(c, frozenset(m))
        if data.get('best'):
            self.offer_ub(set(data['best']))
        self.LB = max(self.LB, data.get('LB') or 0)
        return data


def export_certificate(S, path, status):
    out = {'k': S.k, 'q': S.q, 'status': status, 'LB': S.LB, 'UB': S.UB, 'best': S.best,
           'seed_P': None,
           'note': 'no-goods = expanding cycles (node list) with the minus-signed odd nodes at the time they were '
                   'found; the seeds (expanding cycles of the all-plus graph of length <= k+2) are regenerated by the '
                   'verifier',
           'nogoods': [[c, sorted(m)] for c, m in S.nogoods[S.nseed:]]}
    with gzip.open(path, 'wt') as f:
        json.dump(out, f, separators=(',', ':'))


def verify_certificate(path, time_limit=3600.0):
    """independent re-check: every stored no-good is re-derived (closed walk + expanding); the seeds are
    regenerated; the hitting-set problem below the claimed bound is re-solved (LOWER bound); the stored best set
    is re-verified class (i) (UPPER bound).  Returns (LB_certified, UB)."""
    with gzip.open(path, 'rt') as f:
        data = json.load(f)
    k, q = data['k'], data['q']
    S = IHS(k, q, log=lambda *a: None)
    for c, m in data['nogoods']:
        L.check(nogood_ok(k, q, c, frozenset(m)), "stored no-good is not a genuine expanding cycle")
        S.add(c, frozenset(m))
    for c, m in S.nogoods:
        L.check(nogood_ok(k, q, c, frozenset(m)), "seed no-good is not a genuine expanding cycle")
    LB = data['LB']
    st, _ = S.solve_hs(cutoff=LB - 1, time_limit=time_limit)
    L.check(st == 'cutoff', "a hitting set below the certified lower bound exists")
    UB = None
    if data.get('best'):
        L.check(S.verify_class_i(set(data['best'])), "stored best set is not class (i)")
        UB = len(data['best'])
    return LB, UB, len(S.clauses)


if __name__ == '__main__':
    # python3 procgen_drift_20260926_ihs.py q k [time=SECONDS] [ckpt]
    q = int(sys.argv[1])
    k = int(sys.argv[2])
    tl = 3600.0
    ck = False
    for a in sys.argv[3:]:
        if a.startswith('time='):
            tl = float(a[5:])
        elif a == 'ckpt':
            ck = True
    S = IHS(k, q, log=lambda *a: print(*a, flush=True))
    path = os.path.join(L.SCR, f"ihs_q{q}_k{k}.json") if ck else None
    if path and os.path.exists(path):
        S.load(path)
        print('resumed', len(S.clauses) - S.nseed, 'no-goods; UB', S.UB, flush=True)
    t = time.time()
    st = S.run(time_limit=tl, ckpt=path)
    print(st, 'LB', S.LB, 'UB', S.UB, 'best', S.best, f'{time.time() - t:.0f}s', flush=True)

