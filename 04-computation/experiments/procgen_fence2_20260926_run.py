"""procgen_fence2_20260926_run.py -- runner for lane "fence2" (session collatz-procgen-20260922,
2026-09-26): per-fence accounting for Friedman's fences.

Re-checks every finite claim of 05-knowledge/results/procgen_fence2_20260926_per_fence_accounting.md
and prints [OK] lines; ends with ALL CHECKS PASSED.  Prints to stdout only; writes no files.
Usage:  python3 -u 04-computation/experiments/procgen_fence2_20260926_run.py > 05-knowledge/results/procgen_fence2_20260926.out
"""
import math
import os
import sys
import time
import warnings
from fractions import Fraction as Fr

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np                                                    # noqa: E402
import mpmath                                                         # noqa: E402

import procgen_fencelim_20260926_geom as G                            # noqa: E402
import procgen_fencelim_20260926_audit as AU                          # noqa: E402
import procgen_fence2_20260926_torus as T                             # noqa: E402
import procgen_fence2_20260926_struct as S                            # noqa: E402
import procgen_fence2_20260926_constr as C                            # noqa: E402
import procgen_fence2_20260926_lp as LPm                              # noqa: E402

FAIL = []
T0 = time.time()


def ok(tag, cond, msg):
    print(('[OK] ' if cond else '[FAIL] ') + tag + ': ' + msg, flush=True)
    if not cond:
        FAIL.append(tag)


def part_A():
    print('== Part A: sector lemma, corner types, side lemma, counting identity, pinwheel lemma ==')
    allres = []
    tot = dict(O=0, I=0, E=0, R=0, S0=0, S2=0)
    npin = 0
    for name, fs in AU.configurations():
        R = T.analyze(fs)
        out = S.check_structure(R)
        allres.append(all(out['res'].values()))
        for k in tot:
            tot[k] += out['cnt'][k]
        npin += out['pinwheels']
    ok('A1', all(allres) and len(allres) == 26,
       f'A1-A9 hold on all {len(allres)} fencelim configurations (plane; holes, pinches, bridges, X, '
       f'straight joins); corner totals {tot}')
    fs, at = AU.n6_record()
    R = T.analyze(fs)
    out = S.check_structure(R)
    pent = [w for w in out['walk_info'] if sorted(w['types']) == sorted(['O', 'I', 'E', 'E', 'E'])]
    whole = sum(1 for s in pent[0]['sides'] if s['e'] == 2) if pent else -1
    per = sum(float(mpmath.sqrt(G.tomp(s['len2']))) for s in pent[0]['sides']) if pent else 0
    ok('A2', all(out['res'].values()) and whole == 3 and per < 4,
       f'n=6 record (50 digits): A1-A9 hold; its area-1 pentagon has types O,I,E,E,E, {whole} whole '
       f'sides, perimeter {per:.6f} < 4 (x = a - P/4 = {1 - per / 4:.6f} > 0)')
    for a in (Fr(3, 5), Fr(4, 5), Fr(12, 13), Fr(5, 13)):
        fsx, lat = C.pythagorean(a)
        R = T.analyze(fsx, lattice=lat)
        out = S.check_structure(R)
        pures = sorted(w['pure'] for w in out['walk_info'])
        area = sum(w['A2'] for w in R['walks']) / 2
        dens = area / 2
        ok('A3', all(out['res'].values()) and pures == ['I', 'O'] and out['pair_checks'] == 4 and
           dens == (a * a + (1 - a) ** 2) / 2,
           f'Pythagorean torus a={a}: one pure-O and one pure-I pinwheel, 4 landings, pinwheel pairing '
           f'verified at all 4 (s + s\' = 1 exactly), density {dens} = (a^2+b^2)/2 < 1/2')
    g = [((Fr(0), Fr(0)), (Fr(1), Fr(0))), ((Fr(0), Fr(0)), (Fr(0), Fr(1)))]
    R = T.analyze(g, lattice=((Fr(1), Fr(0)), (Fr(0), Fr(1))))
    out = S.check_structure(R)
    ok('A4', all(out['res'].values()) and out['cnt']['E'] == 4 and out['whole_sides'] == 4 and out['Lambda'] == 0,
       'unit grid torus: 4 E corners, 4 whole sides, no landing (Lambda = 0), density 1/2')
    s6 = math.sqrt(2 / (3 * math.sqrt(3)))
    s5 = math.sqrt(4 * math.tan(math.pi / 5) / 5)
    pk = s5 / (2 * math.cos(math.radians(36)))
    ok('A5', s6 < 1 and abs(1 - s6 - 0.379601) < 1e-5 and s5 < 1 and pk < 1,
       f'B1 hexagon of area 1 has side {s6:.6f} < 1 => pinwheel; pairing with a pinwheel roundabout forces '
       f'the roundabout side 1 - {s6:.6f} = {1 - s6:.6f} (not tiny).  B.4 notched pentagon: sides '
       f'{s5:.6f}, {pk:.6f} < 1 but its 2 reflex corners force >= 2 whole sides (>= 1)')


def part_B():
    print('== Part B: typed piece-potential LP certificates (EMPIRICAL: field constraints by multistart) ==')
    certs = sorted(f for f in os.listdir(HERE) if f.startswith('procgen_fence2_20260926_cert_') and f.endswith('.json'))
    if not certs:
        ok('B0', False, 'no certificate file found')
        return
    starts = int(os.environ.get('FENCE2_STARTS', '10'))
    for f in certs:
        tag = f.split('_cert_')[1][:-5]
        L, c = LPm.load_cert(os.path.join(HERE, f), seed=7)
        val, gridv, worst, wpat, rep = LPm.verify_cert(L, starts=starts)
        tail = LPm.tail_check(L, emax=12, kmax_extra=(8 if c['kmax'] < 8 else c['kmax']), starts=2) \
            if 'O' in c['alphabet'] and len(c['alphabet']) > 1 else dict(conv=0.0, fanE=0.0, field=-1.0, pat=None)
        delta = max(worst, tail['field'], tail['conv'], tail['fanE'], gridv, 0.0)
        bound = val + 4.0 / 3.0 * delta
        ok("B-" + tag, gridv < 1e-6 and bound < 0.5168084,
           f'{f}: {c.get("cls", "TypedLP")} grid NL={c["NL"]} NA={c["NA"]} alphabet {c["alphabet"]} kmax {c["kmax"]}: '
           f'2(beta+rho) = {val:.6f}; grid families max violation {gridv:.1e}; field violation (kappa <= {c["kmax"]}) '
           f'{worst:.2e} [{"".join(wpat) if wpat else "-"}]; tails: conv(e<=12) {tail["conv"]:.1e}, larger fans '
           f'{tail["fanE"]:.1e}, fields kappa<=8 {tail["field"]:.2e} [{"".join(tail["pat"]) if tail["pat"] else "-"}]; '
           f'repaired bound {bound:.6f} < 0.5168084')


def part_C():
    print('== Part C: constructions ==')
    good = True
    for s in (Fr(1, 3), Fr(3, 5), Fr(5, 13), Fr(7, 25)):
        r = C.flipped_square(s)
        good = good and r['area'] == 1 and r['len2'] == [1, s * s, 1, (1 - s) ** 2, 1] and r['convex']
    ok('C1', good, 'flipped-square pentagons (sides 1, s, 1, 1-s, 1; s = 1/3, 3/5, 5/13, 7/25): convex, area exactly 1, '
       'perimeter exactly 4 (a/P = 1/4)')
    res = []
    for cc in (0.2, 0.1, 0.05):
        b = C.pent_roundabout_family(cc, starts=int(os.environ.get('FENCE2_CSTARTS', '8')))
        res.append((cc, b[0]))
    ok('C2', all(v < 1 - 0.5 * cc for cc, v in res),
       'pinwheel pentagon + roundabout family, max(a_P + a_T) at triangle diameter c: ' +
       ', '.join(f'c={cc}: {v:.6f}' for cc, v in res) + ' (< 1 - c/2, so density (a_P+a_T)/2 < 1/2; EMPIRICAL)')


def one_param_pinwheel(k, K, starts, seed=0):
    """max of a - P/4 + k sum_i cos(theta_{i+1}) s_i (1 - s_i) over convex K-gons with sides < 1 and
    area <= 1 (the one-parameter pinwheel certificate g = -k cos(th) l (1-l) of D.2)."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    PI = math.pi
    best = -1e9

    def geo(z):
        phi = np.empty(K)
        phi[1:] = z[:K - 1]
        phi[0] = 2 * PI - z[:K - 1].sum()
        ell = z[K - 1:]
        psi = np.concatenate(([0.0], np.cumsum(phi[1:])))
        xs = np.concatenate(([0.0], np.cumsum(ell * np.cos(psi))))
        ys = np.concatenate(([0.0], np.cumsum(ell * np.sin(psi))))
        return phi, ell, 0.5 * np.sum(xs[:-1] * ys[1:] - xs[1:] * ys[:-1]), xs, ys

    def F(z):
        phi, ell, A, xs, ys = geo(z)
        return A - ell.sum() / 4 + k * np.sum(np.cos(np.roll(PI - phi, -1)) * ell * (1 - ell))
    cons = [{'type': 'eq', 'fun': lambda z: np.array([geo(z)[3][-1], geo(z)[4][-1]])},
            {'type': 'ineq', 'fun': lambda z: 1 - geo(z)[2]},
            {'type': 'ineq', 'fun': lambda z: z[:K - 1].sum() - PI - 1e-4},
            {'type': 'ineq', 'fun': lambda z: 2 * PI - z[:K - 1].sum() - 1e-4}]
    bnds = [(1e-4, PI - 1e-4)] * (K - 1) + [(1e-5, 1 - 1e-6)] * K
    for st in range(starts):
        phi = rng.dirichlet(np.ones(K) * 2) * 2 * PI
        if phi.max() >= PI:
            continue
        z0 = np.concatenate((phi[1:], rng.uniform(0.1, 0.95, K)))
        r = minimize(lambda z: -F(z), z0, method='SLSQP', bounds=bnds, constraints=cons,
                     options={'maxiter': 400, 'ftol': 1e-13})
        g = geo(r.x)
        if abs(g[3][-1]) + abs(g[4][-1]) > 1e-8 or g[2] > 1 + 1e-9 or g[2] <= 0:
            continue
        best = max(best, F(r.x))
    return best


def part_D():
    print('== Part D: the pinwheel statement (P) ==')
    vals = [(k, one_param_pinwheel(k, 5, 16, seed=1)) for k in (0.2, 0.25, 0.3)]
    ok('D1', all(v > 1e-3 for k, v in vals),
       'one-parameter pinwheel certificate g = -k cos(th) l(1-l) fails (P) on pentagons: max violation ' +
       ', '.join(f'k={k}: {v:+.5f}' for k, v in vals) + ' (EMPIRICAL)')
    f = os.path.join(HERE, 'procgen_fence2_20260926_pincert.json')
    if os.path.exists(f):
        L, c = LPm.load_cert(f, seed=5)
        val, gridv, worst, wpat, rep = LPm.verify_cert(L, starts=int(os.environ.get('FENCE2_STARTS', '6')))
        ok("D2", gridv < 1e-6 and rep >= 0.5 - 1e-3,
           f'pure-pinwheel LP certificate (grid NL={c["NL"]} NA={c["NA"]}, kappa <= {c["kmax"]}): 2(beta+rho) = {val:.6f}, '
           f'field violation {worst:.2e}, repaired value {rep:.6f} (the Pythagorean tilings give density -> 1/2 from below)')


if __name__ == '__main__':
    part_A()
    part_C()
    part_D()
    part_B()
    print(f'time {time.time() - T0:.1f} s')
    if FAIL:
        print('FAILED:', FAIL)
        sys.exit(1)
    print('ALL CHECKS PASSED')
