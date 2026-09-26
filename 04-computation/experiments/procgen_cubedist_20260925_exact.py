#!/usr/bin/env python3
"""
procgen_cubedist_20260925_exact.py -- exact delta_k (least number of sign changes of the all-plus strategy
reaching class (i)) by an implicit-hitting-set (IHS) loop; cube-distance lane 2026-09-25 (HYP-9138).

Method (exact; the logic of the cube lane's lazy MaxSAT, with a MIP hitting-set solver):
  * x_r in {0,1} for every odd residue r mod 2^k (x_r = 1: sigma(r) differs from the base strategy).
    Objective: min sum x_r.
  * A cycle C of G_sigma is a cycle of G_sigma' whenever sigma' agrees with sigma on the odd nodes of C (edges out
    of even nodes never depend on sigma).  So an expanding cycle C of G_sigma gives the valid no-good
    "sigma' must differ from sigma at some odd node of C".  Every no-good is a necessary condition for class (i),
    so the hitting-set optimum is a LOWER bound; a hitting set whose G_sigma has no expanding cycle is class (i),
    so the first optimum without expanding cycles is optimal.
  * A no-good is stored as (cycle node list, the odd nodes of the cycle where sigma = -); it is re-checkable
    on its own: the cycle must be a closed walk under those signs and 3^a > 2^p (mul^a > 2^p).
  * Seeds (all-plus base only): every expanding cycle of the all-plus graph G_0 (= de Bruijn B(2,k) in parity-word
    coordinates) of length <= P.
  * Each round: solve the hitting set exactly (HiGHS MIP; integer data), then run greedy repair phases from the
    optimum and from a few perturbations of it, harvesting node-disjoint expanding cycles (C engine; threshold
    q*/r* = best lower approximation of log_mul 2 with denominator <= 2^k, exact for simple cycles).
  * Every reported optimum is re-verified by exact Karp (max cycle density < log_mul 2, compared as mul^a < 2^p).
  mul = 3: Collatz sheet (base 'plus' = Collatz, base 'minus' = 3n-1); mul = 5: DRIFT control (base 'plus' = 5n+1).
"""
import os
import sys
import json
import time
import random
import itertools
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_cubedist_20260925_lib as L  # noqa: E402


def base_minus(k, base):
    """odd residues with sigma = - in the base strategy"""
    return set() if base == 'plus' else set(range(1, 1 << k, 2))


def de_bruijn_expanding_cycles(k, P, mul=3):
    """all expanding cycles of the all-plus graph G_0 of length p <= P, as ordered node lists: for each primitive
    necklace w of length p with mul^(ones) > 2^p, the nodes are the residues whose parity k-word is the window of
    w^infinity at positions 0..p-1 (consecutive windows are joined by edges of G_0)."""
    K = 1 << k
    table = {}
    for r in range(K):
        table[tuple(L.parity_word(r, k, mul))] = r
    out = []
    seen = set()
    for p in range(1, P + 1):
        for ones in range(p + 1):
            if not mul ** ones > 2 ** p:
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
                key = tuple(sorted(set(v for v in cyc if v & 1)))
                if key not in seen:
                    seen.add(key)
                    out.append(cyc)
    return out


def nogood_ok(k, mul, cyc, minus_odd):
    """the stored no-good is genuine: cyc is a closed walk of G_sigma for every sigma with sigma = - exactly on
    minus_odd among the odd nodes of cyc, and it is expanding."""
    K = 1 << k
    H = K >> 1
    p = len(cyc)
    for i in range(p):
        s, t = cyc[i], cyc[(i + 1) % p]
        if s & 1:
            u = (mul * s + (-1 if s in minus_odd else 1)) // 2
        else:
            u = s // 2
        if t % H != u % H:
            return False
    a = sum(1 for v in cyc if v & 1)
    return mul ** a > 2 ** p


def clause_of(nogood, base_min, idx):
    """no-good (cyc, minus_odd) -> (pos, neg) over x-indices: some odd v of cyc must have its sign changed
    relative to the recorded one.  x_v = 1 means sigma(v) != base(v)."""
    cyc, minus_odd = nogood
    pos, neg = [], []
    for v in sorted(set(u for u in cyc if u & 1)):
        rec_minus = v in minus_odd
        rec_x = 1 if (rec_minus != (v in base_min)) else 0
        (pos if rec_x == 0 else neg).append(idx[v])
    return pos, neg


def solve_hitting(nvar, clauses, time_limit=3000, cutoff=None, threads=2):
    """min sum x s.t. every clause (pos, neg): sum_{pos} x + sum_{neg} (1 - x) >= 1  (HiGHS MIP, integer data,
    at most 2 threads).  cutoff = c: only solutions with objective <= c are sought; 'cutoff' is returned when none
    exists (so the true optimum is > c).  Returns (x, status), status in {'optimal', 'infeasible', 'cutoff',
    'timeout'}."""
    import highspy
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    h.setOptionValue('threads', threads)
    h.setOptionValue('time_limit', float(time_limit))
    h.setOptionValue('mip_rel_gap', 0.0)
    inf = highspy.kHighsInf
    for _ in range(nvar):
        h.addVar(0, 1)
    h.changeColsIntegrality(nvar, np.arange(nvar, dtype=np.int32),
                            np.array([highspy.HighsVarType.kInteger] * nvar))
    h.changeColsCost(nvar, np.arange(nvar, dtype=np.int32), np.ones(nvar))
    for pos, neg in clauses:
        idx = list(pos) + list(neg)
        vals = [1.0] * len(pos) + [-1.0] * len(neg)
        h.addRow(1 - len(neg), inf, len(idx), np.array(idx, dtype=np.int32), np.array(vals))
    if cutoff is not None:
        # sum x <= cutoff as an explicit row (exact, integer)
        h.addRow(-inf, float(cutoff), nvar, np.arange(nvar, dtype=np.int32), np.ones(nvar))
    h.run()
    st = h.getModelStatus()
    if st == highspy.HighsModelStatus.kOptimal:
        x = np.round(np.array(h.getSolution().col_value[:nvar])).astype(int)
        return x, 'optimal'
    if st == highspy.HighsModelStatus.kInfeasible:
        return None, ('cutoff' if cutoff is not None else 'infeasible')
    return None, 'timeout'


class IHS:
    def __init__(self, k, mul=3, base='plus', P=None, M=60, seed=0):
        self.k, self.mul, self.base = k, mul, base
        self.res = list(range(1, 1 << k, 2))
        self.idx = {r: i for i, r in enumerate(self.res)}
        self.base_min = base_minus(k, base)
        self.thr = L.expanding_threshold(k, mul)
        self.M = M
        self.nogoods = []          # (cyc, frozenset minus_odd)
        self.clauses = []
        self.seen = set()
        self.counts = {}
        self.rng = random.Random(seed)
        self.shortP = 0
        if base == 'plus':
            for cyc in de_bruijn_expanding_cycles(k, P if P is not None else k + 3, mul):
                self.add(cyc, frozenset())
        self.nseed = len(self.clauses)

    def add(self, cyc, minus_odd):
        cl = clause_of((cyc, minus_odd), self.base_min, self.idx)
        key = (tuple(cl[0]), tuple(cl[1]))
        if key in self.seen:
            return False
        self.seen.add(key)
        self.nogoods.append((list(cyc), frozenset(minus_odd)))
        self.clauses.append(cl)
        for j in cl[0]:
            self.counts[j] = self.counts.get(j, 0) + 1
        return True

    def minus_of(self, changed):
        return self.base_min.symmetric_difference(changed)

    def cycles(self, changed):
        minus = self.minus_of(changed)
        cycs = L.disjoint_cycles(self.k, L.mask_of(minus), self.thr.numerator, self.thr.denominator, self.M, self.mul)
        out = []
        for cyc in cycs:
            mo = frozenset(v for v in cyc if v & 1 and v in minus)
            L.check(nogood_ok(self.k, self.mul, cyc, mo), "engine returned a non-expanding or broken cycle")
            out.append((cyc, mo))
        return out

    def short(self, changed, P, cap=4000):
        """all short expanding cycles through changed-sign nodes, as no-goods"""
        minus = self.minus_of(changed)
        added = 0
        for cyc in L.short_cycles(self.k, L.mask_of(minus), self.thr.numerator, self.thr.denominator, P, cap,
                                  True, self.mul):
            mo = frozenset(v for v in cyc if v & 1 and v in minus)
            L.check(nogood_ok(self.k, self.mul, cyc, mo), "short cycle is not a genuine expanding cycle")
            added += self.add(cyc, mo)
        return added

    def repair(self, changed, max_steps=40):
        cur = set(changed)
        added = 0
        for _ in range(max_steps):
            if self.shortP:
                added += self.short(cur, self.shortP)
            cyc_list = self.cycles(cur)
            if not cyc_list:
                return cur, added
            for cyc, mo in cyc_list:
                added += self.add(cyc, mo)
            for cyc, mo in cyc_list:
                cand = [v for v in cyc if v & 1 and v not in cur]
                if cand:
                    v = max(cand, key=lambda u: (self.counts.get(self.idx[u], 0), u))
                    cur.add(v)
                else:
                    cur.discard(max(v for v in cyc if v & 1))
        return None, added

    def save(self, path, LB, UB, best):
        data = {'k': self.k, 'mul': self.mul, 'base': self.base, 'LB': LB, 'UB': UB, 'best': best,
                'nogoods': [[c, sorted(m)] for c, m in self.nogoods[self.nseed:]]}
        json.dump(data, open(path + '.tmp', 'w'))
        os.replace(path + '.tmp', path)

    def load(self, path):
        data = json.load(open(path))
        for c, m in data['nogoods']:
            L.check(nogood_ok(self.k, self.mul, c, frozenset(m)), "checkpoint no-good fails")
            self.add(c, frozenset(m))
        return data


def exact_delta(k, mul=3, base='plus', P=None, M=60, verbose=True, time_limit=200000, log=print,
                ub_set=None, max_repair=40, extra_repairs=0, checkpoint=False, shortP=0):
    t0 = time.time()
    S = IHS(k, mul, base, P, M)
    S.shortP = shortP
    ckpt = os.path.join(L.SCR, f"ihs_k{k}_m{mul}_{base}.json") if checkpoint else None
    UB, best = None, None
    if ckpt and os.path.exists(ckpt):
        data = S.load(ckpt)
        if data.get('best') is not None:
            ub_set = ub_set if (ub_set is not None and len(ub_set) <= len(data['best'])) else data['best']
        if verbose:
            log(f"    resumed {len(S.clauses) - S.nseed} no-goods from {ckpt}")
    if ub_set is not None:
        minus = S.minus_of(set(ub_set))
        d = L.karp_density(k, L.mask_of(minus), mul)
        L.check(mul ** d.numerator < 2 ** d.denominator, "supplied upper-bound set is not class (i)")
        UB, best = len(ub_set), sorted(ub_set)
    it = 0
    while True:
        it += 1
        x, st = solve_hitting(len(S.res), S.clauses, time_limit=max(60, time_limit - (time.time() - t0)),
                              cutoff=(UB - 1 if UB is not None else None))
        if st == 'cutoff':
            # no hitting set of size <= UB - 1: the verified class-(i) set of size UB is optimal
            minus = S.minus_of(set(best))
            d = L.karp_density(k, L.mask_of(minus), mul)
            L.check(mul ** d.numerator < 2 ** d.denominator, "IHS optimum is not class (i) by Karp")
            if ckpt:
                S.save(ckpt, UB, UB, best)
            return {'status': 'optimal', 'k': k, 'delta': UB, 'changed': sorted(best), 'rho_max': d,
                    'iterations': it, 'clauses': len(S.clauses), 'seed': S.nseed,
                    'seconds': time.time() - t0, 'minus_set': sorted(minus), 'ihs': S}
        if st == 'infeasible':
            return {'status': 'infeasible', 'k': k, 'iterations': it, 'clauses': len(S.clauses), 'ihs': S}
        if st != 'optimal':
            return {'status': 'timeout', 'k': k, 'iterations': it, 'clauses': len(S.clauses), 'UB': UB, 'ihs': S}
        lb = int(x.sum())
        chosen = [S.res[i] for i in range(len(S.res)) if x[i]]
        if UB is not None and lb >= UB:
            chosen = best
        else:
            sol, added = S.repair(chosen, max_repair)
            if sol is not None and (UB is None or len(sol) < UB):
                UB, best = len(sol), sorted(sol)
            for _ in range(extra_repairs):
                start = set(chosen)
                if start:
                    start.discard(S.rng.choice(sorted(start)))
                sol2, a2 = S.repair(start, max_repair)
                added += a2
                if sol2 is not None and (UB is None or len(sol2) < UB):
                    UB, best = len(sol2), sorted(sol2)
            if ckpt:
                S.save(ckpt, lb, UB, best)
            if verbose:
                log(f"    k={k} mul={mul} base={base} it {it}: LB {lb}, UB {UB}, no-goods {len(S.clauses)} "
                    f"(seed {S.nseed}, +{added}), {time.time() - t0:.0f}s")
            if not (UB is not None and lb >= UB):
                continue
            chosen = best
        minus = S.minus_of(set(chosen))
        d = L.karp_density(k, L.mask_of(minus), mul)
        L.check(mul ** d.numerator < 2 ** d.denominator, "IHS optimum is not class (i) by Karp")
        return {'status': 'optimal', 'k': k, 'delta': len(chosen), 'changed': sorted(chosen),
                'rho_max': d, 'iterations': it, 'clauses': len(S.clauses), 'seed': S.nseed,
                'seconds': time.time() - t0, 'minus_set': sorted(minus), 'ihs': S}


def verify_certificate_file(path):
    """Independent re-check of a stored delta_k certificate: every no-good is re-derived (closed walk + expanding),
    the hitting-set MIP over them is re-solved (LOWER bound), and the stored optimum is re-verified class (i) by
    Karp (UPPER bound)."""
    data = json.load(open(path))
    k, mul, base = data['k'], data['mul'], data['base']
    S = IHS(k, mul, base)
    for c, m in data['nogoods']:
        L.check(nogood_ok(k, mul, c, frozenset(m)), "stored no-good is not a genuine expanding cycle")
        S.add(c, frozenset(m))
    best = data['best']
    x, st = solve_hitting(len(S.res), S.clauses, cutoff=len(best) - 1)
    L.check(st == 'cutoff', "a hitting set smaller than the stored optimum exists: certificate incomplete")
    lb = len(best)
    minus = S.minus_of(set(best))
    d = L.karp_density(k, L.mask_of(minus), mul)
    L.check(mul ** d.numerator < 2 ** d.denominator, "stored optimum is not class (i)")
    return lb, len(best), d, len(S.clauses), S.nseed


def export_certificate(k, path, mul=3, base='plus'):
    """write the checkpoint of a (stopped) long IHS run as a gzipped certificate: the stored no-goods (seeds are
    regenerated by the verifier), the certified lower bound LB, and the best class-(i) set found."""
    import gzip
    ck = os.path.join(L.SCR, f"ihs_k{k}_m{mul}_{base}.json")
    d = json.load(open(ck))
    out = {'k': k, 'mul': mul, 'base': base, 'LB': d['LB'], 'best': d['best'],
           'note': 'no-goods = expanding cycles (node list) with the minus-signed odd nodes at the time they were '
                   'found; seeds (all expanding de Bruijn cycles of length <= k+3) are regenerated by the verifier',
           'nogoods': d['nogoods']}
    with gzip.open(path, 'wt') as f:
        json.dump(out, f, separators=(',', ':'))
    return d['LB'], len(d['best'] or [])


if __name__ == '__main__':
    # examples:
    #   python3 procgen_cubedist_20260925_exact.py 2 3 4 5 6 7 8 9          (exact delta_k)
    #   python3 procgen_cubedist_20260925_exact.py 10 ckpt short=14 ub=prune time=5400
    #        (the long k = 10 run: checkpointed in scratch/procgen_cubedist/, resumable, stopped by the time limit)
    #   python3 procgen_cubedist_20260925_exact.py export 10   (checkpoint -> 05-knowledge/results/..._k10_certificate.json.gz)
    mul, base, extra, ck, shortP, ub, tl = 3, 'plus', 0, False, 0, None, 200000
    ks = []
    args = sys.argv[1:]
    if args and args[0] == 'export':
        k = int(args[1])
        res_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '05-knowledge', 'results'))
        path = os.path.join(res_dir, f"procgen_cubedist_20260925_k{k}_certificate.json.gz")
        print(export_certificate(k, path), path)
        sys.exit(0)
    for a in args:
        if a.startswith('mul='):
            mul = int(a[4:])
        elif a.startswith('base='):
            base = a[5:]
        elif a.startswith('extra='):
            extra = int(a[6:])
        elif a.startswith('short='):
            shortP = int(a[6:])
        elif a.startswith('time='):
            tl = float(a[5:])
        elif a == 'ub=prune':
            ub = 'prune'
        elif a == 'ckpt':
            ck = True
        else:
            ks.append(int(a))
    for k in ks:
        ub_set = None
        if ub == 'prune':
            import procgen_cubedist_20260925_prune as PR
            best, stats, nb = PR.prune(k, ('asc', 'desc', 'margin_hi', 'rand1'))
            ub_set = sorted(best)
            print('pruned upper bound', len(ub_set), stats, flush=True)
        r = exact_delta(k, mul, base, extra_repairs=extra, checkpoint=ck, shortP=shortP, ub_set=ub_set,
                        time_limit=tl)
        print({kk: v for kk, v in r.items() if kk != 'ihs'}, flush=True)
