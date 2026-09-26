#!/usr/bin/env python3
"""
procgen_mykk_20260926_search.py -- the (long) exact searches of the expanding-cycle feedback lane.

Usage:  python3 procgen_mykk_20260926_search.py INSTANCE [INSTANCE ...]
  INSTANCE = log3:K | log5:K | log3odd:K | log5odd:K | gt:A/B:K | ge:A/B:K | stepfun:K | packs
Each instance writes scratch/procgen_mykk/results/<name>.json (value, optimal set R, lower-bound certificate
cycles in word encoding, optional packing).  The runner (procgen_mykk_20260926_run.py) never trusts these
files: it re-verifies every certificate (Bellman-Ford + potential for the upper set, RC2 or a disjoint packing
for the lower bound).  The certificates used by the runner are frozen into procgen_mykk_20260926_certs.py by
`python3 procgen_mykk_20260926_search.py freeze`.
"""
import os
import sys
import json
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from procgen_mykk_20260926_lib import *  # noqa

RES = os.path.join(SCR, 'results')
os.makedirs(RES, exist_ok=True)


def log(msg):
    print(msg, flush=True)


def parse(inst):
    kind, rest = inst.split(':', 1)
    if kind in ('log3', 'log5', 'log3odd', 'log5odd'):
        q = int(kind[3])
        return dict(name=inst.replace(':', '_').replace('/', '-'), k=int(rest), thr=('log', q),
                    odd=kind.endswith('odd'))
    if kind in ('gt', 'ge', 'gtodd', 'geodd'):
        frac, k = rest.split(':')
        a, b = frac.split('/')
        return dict(name=inst.replace(':', '_').replace('/', '-'), k=int(k), thr=(kind[:2], int(a), int(b)),
                    odd=kind.endswith('odd'))
    raise SystemExit("bad instance " + inst)


def make_thr(spec):
    if spec[0] == 'log':
        return Thr('log', q=spec[1])
    return Thr(spec[0], spec[1], spec[2])


def run_instance(inst, per_round=1000, pack_time=120):
    opts = {}
    if '@' in inst:
        parts = inst.split('@')
        inst = parts[0]
        for o in parts[1:]:
            key, val = o.split('=')
            opts[key] = int(val)
    P = parse(inst)
    k, odd = P['k'], P['odd']
    thr = make_thr(P['thr'])
    out = os.path.join(RES, P['name'] + '.json')
    if os.path.exists(out):
        log("[skip] %s exists" % out)
        return
    t0 = time.time()
    ck = os.path.join(SCR, 'ckpt_' + P['name'] + '.json')
    seeds = None
    # seed an odd-only run with the cycles of the corresponding all-node run, and vice versa
    other = os.path.join(RES, (P['name'].replace('odd', '') if odd else
                               P['name'].replace('log3_', 'log3odd_').replace('log5_', 'log5odd_')) + '.json')
    inc = None
    if os.path.exists(other) and other != out:
        with open(other) as f:
            od = json.load(f)
        seeds = [word_to_cycle(e, k) for e in od['cert']]
        if not odd and od.get('R'):
            inc = od['R']          # an optimal odd set is feasible for the all-node problem
    lmax = min(2 * k + 2, 22 if k >= 11 else 20)
    if 'lmax' in opts:
        lmax = opts['lmax']
    r = ihs_fvs(k, thr, odd_only=odd, lmax=lmax, per_round=per_round, log=log, verbose=True, ckpt=ck,
                seed_cycles=seeds, mip_time=opts.get('mip'), time_limit=opts.get('time'), incumbent=inc)
    n = 1 << k
    allowed = [v for v in range(n) if (v & 1 or not odd)]
    if not r['optimal']:
        # bounds only: record the best feasible set, the MIP lower bound, and an exact LP-dual fractional
        # packing (a weaker but solver-free lower bound); the runner re-verifies both kinds.
        y, yval, lpv = lp_dual_packing(r['cycles'], allowed)
        keep = [i for i, yc in enumerate(y) if yc > 0]
        res = dict(inst=inst, k=k, thr=list(P['thr']), odd=odd, value=None, lb=r['lb'], ub=r['ub'],
                   R=r['R'], cert=[cycle_to_word(c) for c in r['cycles']], packing=[],
                   frac=[[cycle_to_word(r['cycles'][i]), str(y[i])] for i in keep], frac_value=str(yval),
                   rounds=r['rounds'], ncycles=len(r['cycles']), time=time.time() - t0, rss=peak_rss_mb())
        with open(out.replace('.json', '_bounds.json'), 'w') as f:
            json.dump(res, f)
        log("  %s: NOT finished: %d <= FVS <= %s (%.1fs)" % (inst, r['lb'], r['ub'], time.time() - t0))
        return
    cert = shrink_certificate(r['cycles'], allowed, r['value'])
    log("  %s: value %d, certificate %d cycles (from %d), %.1fs" % (inst, r['value'], len(cert),
                                                                   len(r['cycles']), time.time() - t0))
    # packing attempt: certificate cycles + short pool cycles
    pools = pool_cycles(k, thr, min(k + 6, 18))
    pool = [list(map(int, row)) for l in sorted(pools) for row in pools[l]]
    cand = cert + pool
    nu, ch, opt = max_packing_cpsat(cand, time_limit=pack_time)
    packing = [cand[i] for i in ch]
    log("  packing: %d (optimal within candidates: %s), lengths %s" % (nu, opt, sorted(len(c) for c in packing)))
    res = dict(inst=inst, k=k, thr=list(P['thr']), odd=odd, value=r['value'], R=r['R'],
               cert=[cycle_to_word(c) for c in cert],
               packing=[cycle_to_word(c) for c in packing], packing_opt_in_pool=opt,
               rounds=r['rounds'], ncycles=len(r['cycles']), time=time.time() - t0, rss=peak_rss_mb())
    with open(out + '.tmp', 'w') as f:
        json.dump(res, f)
    os.replace(out + '.tmp', out)
    log("  wrote %s (%.1fs, rss %.0f MB)" % (out, time.time() - t0, peak_rss_mb()))


def run_stepfun(k):
    out = os.path.join(RES, 'stepfun_%d.json' % k)
    if os.path.exists(out):
        log("[skip] %s exists" % out)
        return
    t0 = time.time()
    segs = step_function(k, log=log)
    js = [dict(top=str(g['top']), bottom=str(g['bottom']), V=g['V'], R=g['R'],
               cert=[cycle_to_word(c) for c in g['cert']]) for g in segs]
    with open(out, 'w') as f:
        json.dump(dict(k=k, segs=js, time=time.time() - t0), f)
    log("  wrote %s (%.1fs)" % (out, time.time() - t0))


def run_cover(q, k):
    """certified fractional cover (upper bound for nu* = tau*) for the log_q 2 instance of size k"""
    out = os.path.join(RES, 'cover_%d_%d.json' % (q, k))
    if os.path.exists(out):
        log("[skip] %s exists" % out)
        return
    t0 = time.time()
    src = os.path.join(RES, 'log%d_%d.json' % (q, k))
    with open(src) as f:
        d = json.load(f)
    thr = Thr('log', q=q)
    cert = [word_to_cycle(e, k) for e in d['cert']]
    val, units, D, cyc = tau_star(k, thr, cert)
    res = dict(q=q, k=k, D=D, units=[[int(v), int(u)] for v, u in sorted(units.items())],
               value=str(val), fvs=d['value'], nu_found=len(d['packing']), time=time.time() - t0)
    with open(out, 'w') as f:
        json.dump(res, f)
    log("  cover q=%d k=%d: %s (%.4f), FVS %d, packing %d  (%.1fs)" % (q, k, val, float(val), d['value'],
                                                                      len(d['packing']), time.time() - t0))


def bounds_from_ckpt(inst):
    """write results/<name>_bounds.json from a checkpoint of an interrupted run (cycles, best set, MIP bound)"""
    P = parse(inst)
    k, odd = P['k'], P['odd']
    ck = os.path.join(SCR, 'ckpt_' + P['name'] + '.json')
    with open(ck) as f:
        d = json.load(f)
    thr = make_thr(P['thr'])
    cycles = [list(c) for c in necklace_constraints(k, thr)] + d['cycles']
    R = d['best_R']
    cyc, _ = has_expanding_cycle(k, thr, R)
    check(cyc is None, "checkpoint best set infeasible")
    n = 1 << k
    allowed = [v for v in range(n) if (v & 1 or not odd)]
    res = dict(inst=inst, k=k, thr=list(P['thr']), odd=odd, value=None, lb=d.get('lb'), ub=len(R), R=R,
               cert=[cycle_to_word(c) for c in cycles], packing=[], note='from checkpoint of an interrupted run')
    out = os.path.join(RES, P['name'] + '_bounds.json')
    with open(out, 'w') as f:
        json.dump(res, f)
    log("  wrote %s: %s <= FVS <= %d" % (out, d.get('lb'), len(R)))


def freeze():
    """pack every results/*.json into procgen_mykk_20260926_certs.py (base64 gzip json strings)"""
    data = {}
    for fn in sorted(os.listdir(RES)):
        if fn.endswith('.json'):
            with open(os.path.join(RES, fn)) as f:
                data[fn[:-5]] = json.load(f)
    path = os.path.join(HERE, 'procgen_mykk_20260926_certs.py')
    with open(path, 'w') as f:
        f.write('"""procgen_mykk_20260926_certs.py -- frozen search results (generated by '
                'procgen_mykk_20260926_search.py freeze).\nEvery entry is re-verified by the runner; nothing here is '
                'trusted."""\n')
        f.write('DATA = %r\n' % pack_json(data))
    log("froze %d result files into %s (%d bytes)" % (len(data), path, os.path.getsize(path)))


if __name__ == '__main__':
    for inst in sys.argv[1:]:
        if inst == 'freeze':
            freeze()
        elif inst.startswith('stepfun:'):
            run_stepfun(int(inst.split(':')[1]))
        elif inst.startswith('bounds:'):
            bounds_from_ckpt(inst[len('bounds:'):])
        elif inst.startswith('cover:'):
            _, qq, kk = inst.split(':')
            run_cover(int(qq), int(kk))
        else:
            run_instance(inst)
