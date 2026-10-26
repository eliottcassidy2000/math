#!/usr/bin/env python3
"""A6 orbit to large depth with the certified accelerator; per-depth near misses of the two exact targets (S56).

Targets at depth d (sweep starts J(0,0,W) with 24*2^(d-1) - 3(d-1) - 16 < W <= 24*2^d - 3d - 16):
   T_d = 24*2^d - 3d - 16   (the sweep halts: last W of the depth window)
   N_d = 18*2^d - 3d - 4    (the sweep never returns: certified translated cycler)
For every depth we record the orbit points immediately below and above each target (exact, also through level-2
jumps, whose visited points are b + k*D + e_i), and the sweep index of the depth entry.
Usage: python a6_orbit_deep_20261025.py [seconds] [max_depth]
"""
import sys, os, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a6_accel_20261025 as acc
from a6_fast_20261025 import depth_of


def targets(d):
    return {'N': 18 * (1 << d) - 3 * d - 4, 'T': 24 * (1 << d) - 3 * d - 16}


class Tracker:
    def __init__(self):
        self.rows = {}
        self.cur_d = None

    def _row(self, d):
        if d not in self.rows:
            self.rows[d] = {'N_below': None, 'N_above': None, 'T_below': None, 'T_above': None, 'entry': None}
        return self.rows[d]

    def point(self, W, nsw):
        d = depth_of(0, W)
        r = self._row(d)
        if r['entry'] is None:
            r['entry'] = (nsw, W)
        for d2 in (d - 1, d, d + 1):
            if d2 < 1:
                continue
            r2 = self._row(d2)
            for k, X in targets(d2).items():
                if W <= X:
                    g = X - W
                    if r2[k + '_below'] is None or g < r2[k + '_below']:
                        r2[k + '_below'] = g
                else:
                    g = W - X
                    if r2[k + '_above'] is None or g < r2[k + '_above']:
                        r2[k + '_above'] = g

    def progression(self, b, D, offs, k0, k1, nsw):
        """Orbit points b + k*D + e for k0 <= k <= k1 (e in offs)."""
        lo = b + k0 * D + min(offs)
        hi = b + k1 * D + max(offs)
        dlo, dhi = depth_of(0, max(0, min(lo, hi))), depth_of(0, max(lo, hi))
        for d2 in range(max(1, dlo - 1), dhi + 2):
            r2 = self._row(d2)
            for k, X in targets(d2).items():
                for e in offs:
                    # largest k in [k0,k1] with b + kD + e <= X (D > 0 assumed) and the next one
                    if D <= 0:
                        continue
                    kk = (X - b - e) // D
                    for kc in (kk, kk + 1):
                        if k0 <= kc <= k1:
                            Wp = b + kc * D + e
                            if Wp <= X:
                                g = X - Wp
                                if r2[k + '_below'] is None or g < r2[k + '_below']:
                                    r2[k + '_below'] = g
                            else:
                                g = Wp - X
                                if r2[k + '_above'] is None or g < r2[k + '_above']:
                                    r2[k + '_above'] = g


def run(seconds, max_depth):
    tr = Tracker()
    t0 = time.time()
    Ws = [0]
    toks = []
    occ = {}
    nsw = 0
    tr.point(0, 0)
    last_print = 0
    status = 'time'
    info = None
    gram, pmax = 2, 3000
    computed = 0
    cur_depth = 0
    prof = []
    while time.time() - t0 < seconds:
        W = Ws[-1]
        if depth_of(0, W) > max_depth:
            status = 'depth'
            break
        try:
            W2, t = acc.sweep(W)
        except acc.Halt as h:
            status, info = 'HALT', {'W': W, 'tokens': str(h.args[0])[:2000]}
            break
        except acc.Capture as c:
            b, D, cyc, E = c.args[0]
            status, info = 'CAPTURE', {'W': W, 'base': b, 'translation': D, 'labels': cyc, 'E': E}
            break
        nsw += 1
        computed += 1
        toks.append(t)
        Ws.append(W2)
        tr.point(W2, nsw)
        dW2 = depth_of(0, W2)
        if dW2 > cur_depth:
            cur_depth = dW2
            prof.append((dW2, round(time.time() - t0, 1), nsw, computed, acc.STATS['phases'], acc.STATS['l2_jumps']))
            print('depth', dW2, 'time', prof[-1][1], 'sweeps %.3e' % nsw, 'computed sweeps', computed,
                  'phases', acc.STATS['phases'], 'l2 jumps', acc.STATS['l2_jumps'], flush=True)
        n = len(toks)
        if len(Ws) > 20000:
            del Ws[:10000]
            del toks[:10000]
            occ = {}
            continue
        if n < 2 * gram:
            continue
        key = tuple(toks[-gram:])
        lst = occ.setdefault(key, [])
        found = False
        for i in reversed(lst):
            p = n - i
            if p > pmax or 2 * p > n:
                break
            D = Ws[n] - Ws[n - p]
            if D != Ws[n - p] - Ws[n - 2 * p] or D == 0:
                continue
            if toks[n - p:n] != toks[n - 2 * p:n - p]:
                continue
            cyc = tuple(toks[n - p:n])
            b = Ws[n - 2 * p]
            maxbits = max(abs(x).bit_length() for x in Ws[n - 2 * p:n + 1])
            E = max(1000, 3 * maxbits + 300)

            def valid(K):
                return acc._run_sweeps(b + K * D, cyc) == b + (K + 1) * D
            Kmax = acc._max_valid(valid, E)
            if Kmax is None:
                status, info = 'INFINITE-SWEEP-CYCLE', {'base': b, 'D': D, 'p': p, 'E': E}
                break
            m = Kmax + 1
            offs = [Ws[n - 2 * p + i] - b for i in range(p)]
            tr.progression(b, D, offs, 2, m - 1, nsw)
            new = b + m * D
            tr.point(new, nsw + m * p - 2 * p)
            nsw += m * p - 2 * p
            del Ws[n - 2 * p + 1:]
            del toks[n - 2 * p:]
            toks.append(('cyc2', len(cyc), m))
            Ws.append(new)
            occ = {}
            found = True
            break
        if status not in ('time', 'depth'):
            break
        if not found:
            lst.append(n)
            if len(lst) > 8:
                del lst[0]
        if time.time() - last_print > 60:
            last_print = time.time()
            print(f'[{time.time() - t0:.0f}s] sweeps {nsw:.3e} depth {depth_of(0, Ws[-1])} W bits {Ws[-1].bit_length()}',
                  flush=True)
    run.prof = prof
    return status, info, tr, nsw, Ws[-1]


if __name__ == '__main__':
    secs = float(sys.argv[1]) if len(sys.argv) > 1 else 600
    maxd = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9
    status, info, tr, nsw, W = run(secs, maxd)
    print('status', status, 'sweeps', nsw, 'final depth', depth_of(0, W), 'W bits', W.bit_length())
    if info:
        print('info', info)
    print('depth  entry_sweep  N_below N_above  T_below T_above(overshoot)')
    out = []
    for d in sorted(tr.rows):
        r = tr.rows[d]
        print(d, r['entry'][0] if r['entry'] else None, r['N_below'], r['N_above'], r['T_below'], r['T_above'], flush=True)
        out.append({'d': d, **{k: (v if not isinstance(v, tuple) else list(v)) for k, v in r.items()}})
    json.dump({'status': status, 'info': str(info), 'sweeps': nsw, 'final_W': str(W), 'rows': out, 'profile': run.prof},
              open(os.path.join(HERE, 'a6_orbit_deep.json'), 'w'), indent=1, default=str)
