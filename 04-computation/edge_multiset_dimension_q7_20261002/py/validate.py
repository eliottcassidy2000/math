"""validate.py -- cross-checks for completed driver runs.

usage: python3 validate.py D RUNDIR k1 [k2 ...]
For each k:
  * coverage: re-checks that the chunks of every a tile [0, nrep) of the representative file exactly;
  * leaf counts per (k, a): C search == independent numpy DP (leafcount_dp.py); a-values without a
    representative file must have been excluded by the weight lemma (and are re-derived here);
  * every RESOLVING set printed is re-verified by the pure-Python checker (pycheck.py) and
    canonicalised under the full Aut(Q_D) (src/canon); orbits and stabiliser orders are reported.
Writes RUNDIR/validation_k<k>.json and prints a summary line per k.
"""
import sys, os, json, subprocess
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import leafcount_dp, pycheck

def wmin(n, a):
    s, r = 0, a
    for w in range(n + 1):
        t = min(r, comb(n, w)); s += t * w; r -= t
        if r == 0: return s

def main():
    D = int(sys.argv[1]); rundir = sys.argv[2]; ks = [int(x) for x in sys.argv[3:]]
    n = D - 1
    allrec = [json.loads(l) for l in open(os.path.join(rundir, 'chunks.jsonl')) if l.strip()]
    for k in ks:
        recs = [r for r in allrec if r.get('D') == D and r.get('K') == k]
        out = dict(D=D, k=k, per_a={}, ok=True, problems=[])
        resolving = []
        for a in range((k + 1) // 2, k + 1):
            b = k - a
            rf = leafcount_dp.repfile(n, a)
            ra = [r for r in recs if r['a'] == a]
            if not os.path.exists(rf):
                t = a - 2 * b
                excl = t > 0 and wmin(n, a) > n * b
                good = excl and any(r.get('excluded') for r in ra)
                out['per_a'][a] = dict(excluded_by_weight_lemma=excl, wmin=wmin(n, a), bound=n * b)
                if not good: out['ok'] = False; out['problems'].append('a=%d: no rep file and not excluded' % a)
                continue
            nrep = os.path.getsize(rf) // 8
            cov = sorted((r['first'], r['last']) for r in ra if 'first' in r)
            pos = 0; tiled = True
            for f, l in cov:
                if f != pos: tiled = False
                pos = l
            tiled = tiled and pos == nrep
            c_leaves = sum(r['leaves'] for r in ra if 'first' in r)
            _, dp_leaves = leafcount_dp.leafcount(D, k, a, check_symmetry=True)
            found = sum(r['found'] for r in ra if 'first' in r)
            for r in ra: resolving += r.get('resolving', [])
            out['per_a'][a] = dict(nrep=nrep, chunks=len(cov), tiled=tiled, leaves_C=c_leaves, leaves_DP=dp_leaves,
                                   found=found, sec=round(sum(r['sec'] for r in ra if 'first' in r), 2))
            if not tiled: out['ok'] = False; out['problems'].append('a=%d coverage' % a)
            if c_leaves != dp_leaves: out['ok'] = False; out['problems'].append('a=%d leaves C %d != DP %d' % (a, c_leaves, dp_leaves))
        sets = [sorted(int(x) for x in l.split('set=')[1].split(',')) for l in resolving]
        py_ok = all(len(S) == k and pycheck.is_resolving(D, S) for S in sets)
        out['resolving_leaves'] = len(sets); out['python_reverified'] = py_ok
        if not py_ok: out['ok'] = False; out['problems'].append('python re-verification failed')
        if sets:
            inp = '\n'.join(' '.join(map(str, S)) for S in sets) + '\n'
            co = subprocess.run([os.path.join(ROOT, 'src', 'canon'), str(D)], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
            co = [l.split() for l in co if l.strip()]
            orbits = sorted(set(x[0] for x in co)); stabs = sorted(set(int(x[1]) for x in co))
            out['orbits'] = len(orbits); out['stabilizer_orders'] = stabs; out['orbit_list'] = orbits
        out['leaves_total_C'] = sum(v.get('leaves_C', 0) for v in out['per_a'].values())
        out['leaves_total_DP'] = sum(v.get('leaves_DP', 0) for v in out['per_a'].values())
        out['cpu_sec'] = round(sum(v.get('sec', 0) for v in out['per_a'].values()), 1)
        with open(os.path.join(rundir, 'validation_k%d.json' % k), 'w') as f: json.dump(out, f, indent=1, default=str)
        print('VALIDATE D=%d k=%d ok=%s leaves C=%d DP=%d resolving_leaves=%d python_ok=%s orbits=%s stab=%s cpu=%.1fs problems=%s' % (
            D, k, out['ok'], out['leaves_total_C'], out['leaves_total_DP'], len(sets), py_ok, out.get('orbits', 0),
            out.get('stabilizer_orders', []), out['cpu_sec'], out['problems']), flush=True)

if __name__ == '__main__':
    main()
