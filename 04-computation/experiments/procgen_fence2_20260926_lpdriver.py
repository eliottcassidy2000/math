"""procgen_fence2_20260926_lpdriver.py -- cutting-plane driver for the typed piece-potential LPs of
lane "fence2" (session collatz-procgen-20260922, 2026-09-26).  Writes lp_<tag>.pkl and cert_<tag>.json
into the CURRENT directory (run it from a scratch directory).

usage: [env] python3 -u procgen_fence2_20260926_lpdriver.py NL NA kmax iters tag lib alphabet rmax [init.pkl ...]
env:   LPCLS=2 (TypedLP2, non-separable 0 pieces) or 1 (TypedLP); RESUME=1 (import all field rows of the
       init pickles); NCUTS (cuts per type word); CENTER=prox|ref; PURGE (slack threshold); STARTS
       (pricing starts in full rounds); KMAXR (max convex corners of words with a reflex corner).
The certificates of the note were produced by (in order): the coarse separable run (cert_sep), then
TypedLP2 runs n1 -> n1b -> n1s -> n1t -> main2 (cert_main), nr1 -> nr2 (cert_nonconvex), and the
pure-pinwheel runs pin2 -> pin3 -> pin4 (pincert); see scratch/procgen_fence2/NOTES.md.
"""
import sys, time, pickle, warnings
warnings.filterwarnings('ignore')
import os as _os
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import numpy as np
import procgen_fence2_20260926_lp as LPm
NL, NA, kmax, iters, tag, lib = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6])
alpha = sys.argv[7]
rmax = int(sys.argv[8])
inits = sys.argv[9:]
fan = 1 if 'E' not in alpha else 4
conv = 2 if 'E' not in alpha else 6
t0 = time.time()
import os
cls = LPm.TypedLP2 if os.environ.get('LPCLS', '2') == '2' else LPm.TypedLP
L = cls(NL=NL, NA=NA, kmax=kmax, fan_max=fan, conv_max=conv, seed=2, alphabet=alpha, rmax=rmax, kmax_r=(int(os.environ['KMAXR']) if 'KMAXR' in os.environ else None))
print('patterns', len(L.pats), 'vars', L.nv, flush=True)
n0 = 0
for p in inits:
    D = pickle.load(open(p, 'rb'))
    duals = D.get('duals', None)
    resume = os.environ.get('RESUME', '0') == '1'
    for r, (pat, th, ell) in D.get('geom', {}).items():
        if resume or (duals is not None and r < len(duals) and duals[r] > 1e-9):
            try:
                if L.add_field(pat, th, ell):
                    n0 += 1
            except Exception:
                pass
    for pat, zs in D.get('pool', {}).items():
        L.pool.setdefault(pat, []).extend(zs[-2:])
print('imported rows', n0, flush=True)
print('library rows', L.seed_library(lib) if os.environ.get('RESUME', '0') != '1' else 0, flush=True)
L.ncuts = int(os.environ.get('NCUTS', '3'))
L.center_mode = os.environ.get('CENTER', 'ref')
def save():
    val = L.solve()
    pickle.dump(dict(x=L.x, rows=L.rows, duals=L.duals, NL=NL, NA=NA, kmax=kmax, val=val, pool=L.pool,
                     alphabet=alpha, rmax=rmax, geom=getattr(L, 'geom', {})), open(f'lp_{tag}.pkl', 'wb'))
    LPm.export_cert(L, f'cert_{tag}.json')
    return val
done = 0
if os.environ.get('RESUME', '0') == '1':
    v0 = L.solve()
    print('resumed value', v0, 'purged', L.purge(float(os.environ.get('PURGE', '0.02'))), 'rows', len(L.rows), flush=True)
while done < iters:
    h = L.run(iters=5, starts=int(os.environ.get('STARTS', '4')), full_every=5)
    done += 5
    v = save()
    print('saved', done, v, 'time', round(time.time() - t0), 'lastworst', h[-1][3], flush=True)
    npg = L.purge(float(os.environ.get('PURGE', '0.02')))
    print('purged', npg, 'rows', len(L.rows), 'center_fail', getattr(L, 'center_fail', 0), flush=True)
    if len(h) < 5:   # converged (no violated rows in a FULL round)
        print('converged', flush=True)
        break
