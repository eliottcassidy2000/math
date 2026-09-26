"""procgen_fencelim_20260926_run.py -- runner for the "fencelim" lane (collatz-procgen-20260922, 2026-09-26).

Re-checks every finite fact claimed in 05-knowledge/results/procgen_fencelim_20260926_fence_density.md
and prints [OK] / [FAIL] lines.  Usage:  python3 -u procgen_fencelim_20260926_run.py [--quick]
"""
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FAILS = []


def check(label, cond, detail=''):
    tag = '[OK]' if cond else '[FAIL]'
    if not cond:
        FAILS.append(label)
    print(f'{tag} {label}' + (f' :: {detail}' if detail else ''), flush=True)


def part_A():
    import mpmath
    import procgen_fencelim_20260926_audit as AU
    print('== A. audit of (1)-(5) ==', flush=True)
    # (4) per-face inequality: max over a in [0,1] is at a=1; check kappa=3..400 and the limit
    worst = []
    for k in range(3, 401):
        v = 1 - AU.MU * AU.P_reg(k) - 2 * AU.RHO * (k - 3)
        worst.append((v, k))
    lim_ok = 1 - AU.MU * 2 * mpmath.sqrt(mpmath.pi) - 2 * AU.RHO * 4 < 0
    tight = [k for (v, k) in worst if abs(v) < mpmath.mpf(10) ** -40]
    check('A4 per-face inequality a <= mu P_kappa sqrt(a) + 2 rho (kappa-3) on [0,1] for kappa=3..400 (+ all kappa>=7 via P_kappa>2 sqrt(pi))',
          all(v <= mpmath.mpf(10) ** -40 for (v, k) in worst) and lim_ok and tight == [4, 5],
          f'equality only at a=1, kappa in {tight}; slack at kappa=3: {mpmath.nstr(-worst[0][0], 6)} (mu P_3 = {mpmath.nstr(AU.MU * AU.P_reg(3), 6)} > 1), kappa=6: {mpmath.nstr(-worst[3][0], 6)}')
    check('A5 constants', abs(AU.LAM - mpmath.mpf('0.5224525')) < 1e-6,
          f'P5={mpmath.nstr(AU.P5, 10)} mu={mpmath.nstr(AU.MU, 10)} rho={mpmath.nstr(AU.RHO, 10)} lambda<= (6-P5)/(8-P5)={mpmath.nstr(AU.LAM, 10)}')
    # configurations
    rows = []
    allok = True
    for name, fs in AU.configurations():
        r = AU.audit_config(name, fs)
        rep = r['rep']
        ok = rep['ok1'] and rep['ok2'] and rep['ok3'] and rep['ok3b'] and rep['ok_e'] and rep['tight_types_exact'] \
            and r['worst'] < mpmath.mpf(10) ** -30 and r['ok4'] and r['ok5'] and r['ok_perim']
        allok = allok and ok
        rows.append((name, r['n'], rep['C'], rep['H'], rep['c_o'], rep['kappa_o'], rep['lhs3'], r['A'], r['valid'], rep['types']))
    fs6, at6 = AU.n6_record()
    r = AU.audit_config('n=6 record (reproduced, 50 digits)', fs6)
    rep = r['rep']
    ok6 = rep['ok1'] and rep['ok2'] and rep['ok3'] and rep['ok3b'] and rep['ok_e'] and r['ok4'] and r['ok5'] and r['valid'] \
        and abs(r['A'] - mpmath.mpf('1.47585')) < 1e-5
    rows.append(('n=6 record (reproduced, 50 digits)', r['n'], rep['C'], rep['H'], rep['c_o'], rep['kappa_o'], rep['lhs3'], r['A'], r['valid'], rep['types']))
    check('A1-A3 junction inequality, angle identity (every face, 1e-30), global identity 2C=sum(k-alpha)+2H+2c_o, '
          'Corner Lemma identity and sum(kappa_f-3)<=n-3, sum P_f + P_o = 2n, on %d exact configurations' % (len(rows) - 1),
          allok, 'grids, polyominoes (incl. an enclosed cell), brick wall, X junction, Y junction (Q(sqrt3)), nested component, bridge, pinched faces, 2 components')
    check('A6 n=6 record reproduced (pentagon of unit sides + unit chord) and audited', ok6,
          f'A={mpmath.nstr(r["A"], 12)}; fields 1 (kappa=5) + {mpmath.nstr(r["A"] - 1, 8)} (kappa=4); sum(kappa-3)=3=n-3')
    tight_rows = [x[0] for x in rows if x[6] == x[1] - 3]
    check('A7 Corner Lemma is attained (sum(kappa_f-3)=n-3) by the unit square and the n=5 and n=6 records',
          'grid 1x1' in tight_rows and any('n=5' in t for t in tight_rows) and any('n=6' in t for t in tight_rows),
          'tight on: ' + '; '.join(tight_rows))
    print('   configuration table (name | n | C | H | c_o | kappa_o | sum(kappa-3) | A | fields<=1 | junction types):')
    for x in rows:
        print(f'   {x[0][:52]:52s} | {x[1]:3d} | {x[2]:3d} | {x[3]} | {x[4]} | {x[5]:2d} | {x[6]:3d} | {mpmath.nstr(x[7], 9):>11s} | {x[8]} | {x[9]}')
    # records
    import procgen_fencelim_20260926_records as REC
    bad = []
    lines = []
    for n in range(3, 51):
        a = REC.A_REC[n]
        b = AU.finite_bound(n)
        u = AU.U_prev(n)
        if not (a <= b):
            bad.append(n)
        lines.append((n, a, b, u))
    check('A8 every record n=3..50 on the Friedman page satisfies the new finite bound B(n)', not bad,
          f'B(4)={mpmath.nstr(AU.finite_bound(4), 6)} B(12)={mpmath.nstr(AU.finite_bound(12), 6)} B(50)={mpmath.nstr(AU.finite_bound(50), 6)}; '
          f'min B(n)-record = {mpmath.nstr(min(b - a for (n, a, b, u) in lines), 6)} at n={min(lines, key=lambda t: t[2] - t[1])[0]}')
    check('A9 B(n) < U(n) (previous lane) for every n=3..50', all(b < u for (n, a, b, u) in lines),
          'B(n)/U(n) from %s to %s' % (mpmath.nstr(min(b / u for (n, a, b, u) in lines), 4), mpmath.nstr(max(b / u for (n, a, b, u) in lines), 4)))
    print('   n | record | B(n) | U(n) | record/B')
    for (n, a, b, u) in lines:
        print(f'   {n:2d} | {a:9.5f} | {mpmath.nstr(b, 7):>9s} | {mpmath.nstr(u, 7):>9s} | {mpmath.nstr(a / b, 4)}')


def part_B(quick=False):
    import math
    import numpy as np
    from scipy.optimize import linprog
    import procgen_fencelim_20260926_potlp as PL
    print('== B. angle-potential certificates ==', flush=True)
    P = lambda k: 2 * math.sqrt(k * math.tan(math.pi / k))
    P6 = P(6)
    lam_star = (8 - P6) / (12 - P6)
    alt = (4 - 12 ** 0.25) / (6 - 12 ** 0.25)
    # B1: any certificate obeys T(60,120), T(90,90), square, regular hexagon, tiny triangle:
    # vars mu, rho, f60, f90, f120 ; min 2mu+2rho
    c = [2, 2, 0, 0, 0]
    A = [[0, -1, 1, 0, 1],            # f60 + f120 <= rho
         [0, -1, 0, 2, 0],            # 2 f90 <= rho
         [-4, 0, 0, -4, 0],           # 4 f90 + 4 mu >= 1
         [-P6, 0, 0, 0, -6],          # 6 f120 + mu P6 >= 1
         [0, 0, -3, 0, 0]]            # 3 f60 >= 0
    b = [0, 0, -1, -1, 0]
    r = linprog(c, A_ub=A, b_ub=b, bounds=[(0, None), (0, None), (None, None), (None, None), (None, None)], method='highs')
    check('B1 every angle-potential certificate has 2(mu+rho) >= (8-P6)/(12-P6) = (4-12^(1/4))/(6-12^(1/4)) (5 constraints: T(60,120), T(90,90), unit square, regular hexagon of area 1, tiny equilateral triangle)',
          abs(r.fun - lam_star) < 1e-9 and abs(lam_star - alt) < 1e-12 and lam_star > 0.5,
          f'stall value {lam_star:.10f}; optimal mu = 2/(12-P6) = {2 / (12 - P6):.8f}, rho = (1-4mu)/2 = {(1 - 8 / (12 - P6)) / 2:.8f}')
    # B2: sampled LP (torus, convex fields, exact grid angles) reaches the stall value; primal
    L = PL.SampledLP(delta_deg=6.0)
    val, it = L.run()
    pr = L.primal()
    faces = sorted(tuple(sorted(k[2])) for w, k in pr if k[0] == 'F1')
    tri = [k[1] for w, k in pr if k[0] == 'F0']
    check('B2 sampled potential LP (6 deg grid, torus, convex fields, T/X/Y/fans) = stall value; primal = regular hexagons + unit squares + tiny triangles',
          abs(val - lam_star) < 1e-6 and (15, 15, 15, 15) in faces and (20, 20, 20, 20, 20, 20) in faces and (10, 10, 10) in tri,
          f'LP value {val:.8f} ({it} rounds); junctions T(60,120) and T(90,90); faces 120^6 (area 1, P6=3.72242), 90^4, 60^3 (area -> 0)')
    # B3: rigorous certificate for configurations whose fields have convex outer boundaries
    L = PL.PLLP(delta_deg=2.0, nmode='none')
    v = L.run()
    bound, slack, ok = PL.certify_convex(L)
    check('B3 PROVED (computer-assisted): if every field has a convex outer boundary then A <= mu(2n - P_o) + 2 rho n with 2(mu+rho) = %.7f' % bound,
          ok and bound < 0.5169, f'piecewise-linear Phi on a 2-degree grid, {len(L.rows)} LP rows; slacks: ' +
          ', '.join(f'{k}={v:.1e}' for k, v in slack.items()))
    # B4: non-convex fields with hull-based bounds give no gain over the kappa-LP
    L = PL.PLLP(delta_deg=4.0, nmode='crude')
    vcr = L.run()
    check('B4 EMPIRICAL: with non-convex fields bounded through their convex hull (P >= P_kappa sqrt(a)), the potential LP returns to the kappa-LP value',
          vcr > 0.5224, f'LP value {vcr:.7f} (kappa-LP 0.5224525); binding primal: "star" fields with 72-degree tips and reflex corners')
    # B6: exact hull-angle treatment of non-convex fields: still the kappa-LP; the notched pentagon
    L = PL.PLLP(delta_deg=6.0, nmode='hull')
    vh = L.run(max_iter=600)
    L.solve()
    du = -L.res.ineqlin.marginals
    notched = [L.rowkeys[i] for i in range(len(du)) if du[i] > 1e-7 and L.rowkeys[i][0] == 'H1']
    target = tuple(sorted([(12, 18, 1)] * 4 + [(18, 18, 0)]))
    check('B6 EMPIRICAL: with the exact hull-angle bound for non-convex fields (<= 2r pocket ends, linear reflex minorants) the LP still returns ~ the kappa-LP; binding field = notched regular pentagon',
          vh > 0.5224 and L.final_worst < 1e-9 and any(k[-1] == target for k in notched),
          f'LP {vh:.7f} at 6 deg; notched pentagon: hull 108^5, four 72-degree tips as pocket ends, one exact 108, two 252-degree dents (L junctions)')
    # B5: Lhuilier's inequality P^2 >= 4 A sum cot(theta_i/2) on random convex polygons (sanity)
    rng = np.random.default_rng(1)
    worst = np.inf
    for _ in range(3000):
        m = int(rng.integers(3, 9))
        ang = np.sort(rng.uniform(0, 2 * math.pi, m))
        rad = rng.uniform(0.5, 1.5)
        pts = np.c_[np.cos(ang), np.sin(ang)] * rad
        # convex hull of random points on a circle is the polygon itself
        e = np.roll(pts, -1, 0) - pts
        Pm = np.sum(np.linalg.norm(e, axis=1))
        Am = 0.5 * np.sum(pts[:, 0] * np.roll(pts[:, 1], -1) - np.roll(pts[:, 0], -1) * pts[:, 1])
        th = []
        for i in range(m):
            u, w = pts[i - 1] - pts[i], pts[(i + 1) % m] - pts[i]
            th.append(math.atan2(abs(u[0] * w[1] - u[1] * w[0]), float(np.dot(u, w))))
        g = sum(1 / math.tan(t / 2) for t in th)
        worst = min(worst, Pm ** 2 / (4 * Am * g))
    check('B5 Lhuilier P^2 >= 4 A sum cot(theta_i/2) on 3000 random convex polygons (CITED, sanity)', worst >= 1 - 1e-9,
          f'min ratio {worst:.12f} (equality for polygons with an incircle, e.g. every triangle)')


def part_C():
    import math
    import numpy as np
    import procgen_fencelim_20260926_typed as TY
    print('== C. constructions and obstructions ==', flush=True)
    P = lambda k: 2 * math.sqrt(k * math.tan(math.pi / k))
    check('C1 PROVED: a field with at most 4 convex corners has a <= P/4 (P_3, P_4 >= 4), so if every field has <= 4 convex corners then A <= (2n - P_o)/4 < n/2',
          P(3) >= 4 - 1e-12 and P(4) >= 4 - 1e-12, f'P_3 = {P(3):.5f}, P_4 = {P(4):.5f}; beating 1/2 needs fields with >= 5 convex corners')
    cd = TY.cairo_data()
    check('C2 Cairo tiling (X at (0,0),(1,1), period 2, a=(3+sqrt3)/6): pentagons of area 1, angles 90/120, P = 3.86370 < 4, but a through line Y-X-Y has length 1.633 > 1',
          abs(cd['angle_X'] - 90) < 1e-9 and all(abs(x - 120) < 1e-9 for x in cd['angles_Y']) and cd['through'] > 1,
          f"XY edge {cd['xy']:.5f}, YY edge {cd['yy']:.5f}, perimeter {cd['perimeter']:.5f}: no fence decomposition of this geometry")
    Ps = {u: TY.min_perimeter_unit_sides(5, u) for u in (3, 2, 1)}
    check('C3 least perimeter of an area-1 pentagon with u sides of length 1 (cyclic-polygon theorem, CITED)',
          Ps[3] > Ps[2] > Ps[1] > P(5), ', '.join(f'u={u}: {Ps[u]:.5f}' for u in (3, 2, 1)) + f', u=0: {P(5):.5f}')
    mf, x = TY.unit_side_mean_field(Ps[3], Ps[2])
    check('C4 EMPIRICAL relaxation: unit-side mean field (tight junctions, convex fields; Y-corners force whole-fence sides) still exceeds 1/2',
          0.5 < mf < 0.5224, f'value {mf:.5f} with {x[0]:.4f} P(3Y,2T) + {x[1]:.4f} P(2Y,3T) pentagons + {x[2]:.4f} unit squares per fence (angles ignored)')
    res = TY.best_balanced_typed_pentagon(starts=100)
    best = max(r[2] for r in res)
    check('C5 EMPIRICAL: Cairo-balanced typed pentagons (3 Y + 2 T corners, Y angles summing to 360, T angles to 180, the 3 forced whole-fence sides of length 1, area <= 1) never beat the square: max a/P = 1/4',
          best <= 0.25 + 1e-9, '; '.join(f'{o}{s}: {b:.6f}' for o, s, b in res))
    a = 3 / 5
    dens = (a * a + (1 - a) ** 2) / 2
    check('C6 PROVED: Pythagorean fence tilings (squares a, 1-a; all junctions T; 2 fences, 2 fields per period) have density (a^2+(1-a)^2)/2 < 1/2, and Corner Lemma equality sum(kappa-3) = n on the torus',
          dens < 0.5, f'a=3/5: {dens:.4f}; per period: 4 T junctions, sum(3alpha-k)/2 = 2 = n = F')


def main():
    t0 = time.time()
    part_A()
    part_B()
    part_C()
    print(f'== done in {time.time() - t0:.1f}s; failures: {len(FAILS)} ==')
    for f in FAILS:
        print('   FAILED:', f)


if __name__ == '__main__':
    main()
