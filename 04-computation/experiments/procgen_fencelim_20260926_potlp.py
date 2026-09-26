"""procgen_fencelim_20260926_potlp.py -- Task B: angle-potential certificates for the fence density.

Session collatz-procgen-20260922, lane "fencelim" (2026-09-26).

Certificate family.  Prices mu (perimeter) and rho (fence end), and a corner potential
    Phi(theta) = c0 + c1 theta + D[k]   for theta in I_k = [k delta, (k+1) delta), 0 <= k < N = pi/delta,
    Phi(pi)    = 0,
    Phi(psi)   = S[l]                   for psi - pi in (l delta, (l+1) delta]  (reflex corners).
The linear part c0 + c1 theta is summed exactly (angle sums at junctions and faces are fixed), so
only the correction D carries discretisation error.  Constraints, for all real angle tuples
(interval compatibility is taken with closed intervals, hence conservatively):
  (J) every junction j:  sum over its sectors of Phi <= rho e_j.  Families: fan(m) = m stems on one
      side of a through fence (m+1 convex sectors, sum pi, bound rho m; m=1 is T; X = two T's);
      conv(e) = e ends, all sectors convex (sum 2 pi; e=3 is Y); refl(e) = e ends with one reflex
      sector (e=2 is L); straight(e) is implied by fan(e-2).
  (F) every convex field with corners theta_i:  sum Phi >= max(0, 1 - 2 mu sqrt(g)), g = sum cot(theta_i/2)
      (Lhuilier: P^2 >= 4 a g for a convex polygon with these angles; a <= 1 and convexity in sqrt(a)).
      Used conservatively:  g >= G_lo = sum_k min_{I_k}(cot(t/2) - sig t) + sig (kappa-2) pi,
      sqrt >= L = piecewise-linear interpolant of sqrt on [G_A, G0], and mu >= 1/(2 sqrt(G0)).
  (N) every non-convex field (kappa convex corners):  sum Phi + sum Phi_reflex >= max(0, 1 - mu P_kappa)
      (hull bound), using Phi_reflex(psi) >= beta (psi - pi) with beta >= c1.
  (O) Phi(theta) >= a + c1 theta (all convex theta), Phi_reflex(psi) >= c1 (psi - pi), a + c1 pi >= 0,
      c1 >= 0: then every outer-face walk and every hole walk has sum Phi >= 0.
  Tails (kappa > kappa_max, many stems/ends) are closed by linear minorant/majorant rows.
Then A <= mu (2n - P_o) + 2 rho n, hence lambda <= 2 (mu + rho).
The LP is solved by cutting planes; separation is exact dynamic programming over interval multisets.
"""
import math
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

NEG = -1e18
POS = 1e18


def dp_table(w, cmax, smax, maximize=False):
    """D[c][s] = best sum of c items (item k: weight w[k], size k) of total size s; choice tables."""
    N = len(w)
    fill = NEG if maximize else POS
    D = [np.full(smax + 1, fill)]
    D[0][0] = 0.0
    CH = [None]
    for c in range(1, cmax + 1):
        prev = D[-1]
        cur = np.full(smax + 1, fill)
        ch = np.full(smax + 1, -1, dtype=np.int32)
        for k in range(min(N, smax + 1)):
            cand = prev[:smax + 1 - k] + w[k]
            tgt = cur[k:]
            better = (cand > tgt) if maximize else (cand < tgt)
            tgt[better] = cand[better]
            ch[k:][better] = k
        D.append(cur)
        CH.append(ch)
    return D, CH


def backtrack(CH, c, s):
    items = []
    while c > 0:
        k = int(CH[c][s])
        items.append(k)
        s -= k
        c -= 1
    return items


def h_min(k, delta, sig):
    """min over t in [k delta, (k+1) delta] of cot(t/2) - sig t  (convex in t)."""
    lo, hi = k * delta, (k + 1) * delta
    f = lambda t: (1 / math.tan(t / 2) if t > 0 else 1e9) - sig * t
    # stationary point: -1/(2 sin^2(t/2)) = sig  (sig < 0)
    cands = [hi]
    if lo > 0:
        cands.append(lo)
    if sig < 0:
        s2 = -1 / (2 * sig)
        if 0 < s2 < 1:
            t = 2 * math.asin(math.sqrt(s2))
            if lo <= t <= hi:
                cands.append(t)
    return min(f(t) for t in cands)


class PotLP:
    def __init__(self, delta_deg=1.0, kappa_max=10, fan_max=10, conv_max=9, refl_max=8,
                 G_A=2.0, G0=5.0, sig=-0.8, chords=None, outer=True, verbose=False):
        self.dd = delta_deg
        self.N = int(round(180 / delta_deg))
        assert abs(self.N * delta_deg - 180) < 1e-9
        self.U = self.N
        self.delta = math.pi / self.N
        self.kappa_max, self.fan_max, self.conv_max, self.refl_max = kappa_max, fan_max, conv_max, refl_max
        self.G_A, self.G0, self.sig = G_A, G0, sig
        if chords is None:
            chords = [G_A, 2.6, 3.0, 3.2, 3.35, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.25, 4.5, G0]
        self.chords = chords
        self.outer = outer
        self.verbose = verbose
        N = self.N
        self.iD = 0
        self.iS = N
        names = ['mu', 'rho', 'c0', 'c1', 'a', 'beta', 'Dmax', 'Smax']
        for i, nm in enumerate(names):
            setattr(self, 'i' + nm, 2 * N + i)
        self.nv = 2 * N + len(names)
        self.hmin = np.array([h_min(k, self.delta, sig) for k in range(N)])
        self.rows = []
        self.keys = set()
        self.P = [None, None, None] + [2 * math.sqrt(k * math.tan(math.pi / k)) for k in range(3, 400)]

    def add(self, coefs, rhs, key=None):
        if key is not None:
            if key in self.keys:
                return False
            self.keys.add(key)
        self.rows.append((coefs, rhs))
        return True

    # --------------------------------------------------------------------------- static rows
    def static_rows(self):
        N, d = self.N, self.delta
        pi = math.pi
        self.add({self.imu: -1.0}, -1 / (2 * math.sqrt(self.G0)), key=('mu_min',))
        # (N) reflex minorant: S[l] >= beta (l+1) delta, beta >= c1
        self.add({self.ic1: 1.0, self.ibeta: -1.0}, 0.0, key=('beta>=c1',))
        for l in range(N):
            self.add({self.ibeta: (l + 1) * d, self.iS + l: -1.0}, 0.0, key=('Smin', l))
        if self.outer:
            # (O) D[k] >= a - c0 ; a + c1 pi >= 0 ; c1 >= 0 ; S[l] >= c1 (l+1) delta (implied by beta>=c1)
            for k in range(N):
                self.add({self.ia: 1.0, self.ic0: -1.0, self.iD + k: -1.0}, 0.0, key=('Omin', k))
            self.add({self.ia: -1.0, self.ic1: -pi}, 0.0, key=('Oab',))
            self.add({self.ic1: -1.0}, 0.0, key=('c1>=0',))
            # face tail kappa > kappa_max: (K+1)(a + c1 pi) - 2 pi c1 >= max(0, 1 - 2 mu sqrt(pi))
            K1 = self.kappa_max + 1
            self.add({self.ia: -K1, self.ic1: -(K1 - 2) * pi}, 0.0, key=('Ftail0',))
            self.add({self.ia: -K1, self.ic1: -(K1 - 2) * pi, self.imu: -2 * math.sqrt(pi)}, -1.0, key=('Ftail1',))
        # junction tails: D[k] <= Dmax, S[l] <= Smax, c1 >= 0 (needed for the refl tail)
        for k in range(N):
            self.add({self.iD + k: 1.0, self.iDmax: -1.0}, 0.0, key=('Dmax', k))
            self.add({self.iS + k: 1.0, self.iSmax: -1.0}, 0.0, key=('Smax', k))
        self.add({self.ic0: 1.0, self.iDmax: 1.0, self.irho: -1.0}, 0.0, key=('c0+Dmax<=rho',))
        M = self.fan_max + 1
        self.add({self.ic0: M + 1, self.iDmax: M + 1, self.ic1: pi, self.irho: -M}, 0.0, key=('Jtail_fan',))
        E = self.conv_max + 1
        self.add({self.ic0: E, self.iDmax: E, self.ic1: 2 * pi, self.irho: -E}, 0.0, key=('Jtail_conv',))
        R = self.refl_max + 1
        self.add({self.ic0: R - 1, self.iDmax: R - 1, self.ic1: pi, self.iSmax: 1.0, self.irho: -R}, 0.0, key=('Jtail_refl',))
        if not self.outer:
            self.add({self.ic1: -1.0}, 0.0, key=('c1>=0',))

    def initial(self):
        self.static_rows()
        N, U = self.N, self.U
        pi = math.pi
        for k in range(N):
            for l in range(k, N):
                if U - 2 <= k + l <= U:
                    self.add({self.ic0: 2.0, self.ic1: pi, self.iD + k: 1.0, self.iD + l: 1.0 if l != k else 2.0,
                              self.irho: -1.0} if l != k else
                             {self.ic0: 2.0, self.ic1: pi, self.iD + k: 2.0, self.irho: -1.0}, 0.0, key=('fan', 1, (k, l)))

    # --------------------------------------------------------------------------- LP
    def solve(self):
        data, ri, ci, b = [], [], [], []
        for r, (coefs, rhs) in enumerate(self.rows):
            for v, c in coefs.items():
                ri.append(r)
                ci.append(v)
                data.append(c)
            b.append(rhs)
        A = csr_matrix((data, (ri, ci)), shape=(len(self.rows), self.nv))
        c = np.zeros(self.nv)
        c[self.imu] = 2.0
        c[self.irho] = 2.0
        bounds = [(-3, 3)] * (2 * self.N) + [(0, 1), (0, 1), (-3, 3), (-3, 3), (-3, 3), (-3, 3), (-3, 3), (-3, 3)]
        res = linprog(c, A_ub=A, b_ub=np.array(b), bounds=bounds, method='highs')
        if res.status != 0:
            raise RuntimeError(res.message)
        self.x = res.x
        return res.fun

    def get(self, name):
        return self.x[getattr(self, 'i' + name)]

    # --------------------------------------------------------------------------- separation
    def jrow(self, items, extra, rhs_mult_rho, lin_c0, lin_c1):
        coefs = {self.ic0: lin_c0, self.ic1: lin_c1, self.irho: -float(rhs_mult_rho)}
        for k in items:
            coefs[self.iD + k] = coefs.get(self.iD + k, 0) + 1.0
        for (v, c) in extra:
            coefs[v] = coefs.get(v, 0) + c
        return coefs

    def sep_junctions(self, tol=1e-10):
        N, U = self.N, self.U
        pi = math.pi
        Dv = self.x[self.iD:self.iD + N]
        Sv = self.x[self.iS:self.iS + N]
        rho, c0, c1 = self.get('rho'), self.get('c0'), self.get('c1')
        cuts, worst = 0, 0.0
        cmax = max(self.fan_max + 1, self.conv_max)
        D, CH = dp_table(Dv, cmax, 2 * U, maximize=True)
        for m in range(1, self.fan_max + 1):
            c = m + 1
            lo, hi = max(0, U - c), U
            s = int(np.argmax(D[c][lo:hi + 1])) + lo
            v = D[c][s] + c * c0 + c1 * pi - rho * m
            worst = max(worst, v)
            if v > tol:
                items = backtrack(CH, c, s)
                if self.add(self.jrow(items, [], m, c, pi), 0.0, key=('fan', m, tuple(sorted(items)))):
                    cuts += 1
        for e in range(3, self.conv_max + 1):
            lo, hi = max(0, 2 * U - e), 2 * U
            s = int(np.argmax(D[e][lo:hi + 1])) + lo
            v = D[e][s] + e * c0 + 2 * pi * c1 - rho * e
            worst = max(worst, v)
            if v > tol:
                items = backtrack(CH, e, s)
                if self.add(self.jrow(items, [], e, e, 2 * pi), 0.0, key=('conv', e, tuple(sorted(items)))):
                    cuts += 1
        for e in range(2, self.refl_max + 1):
            c = e - 1
            best = (NEG, None)
            for l in range(N):
                lo, hi = max(0, U - l - e), U - l
                if hi < 0:
                    continue
                s = int(np.argmax(D[c][lo:hi + 1])) + lo
                # convex sectors sum to 2pi - psi in [pi - (l+1) delta, pi - l delta)
                lin = max(c1 * (pi - (l + 1) * self.delta), c1 * (pi - l * self.delta))
                v = D[c][s] + c * c0 + lin + Sv[l] - rho * e
                if v > best[0]:
                    best = (v, (l, s))
            v, (l, s) = best
            worst = max(worst, v)
            if v > tol:
                items = backtrack(CH, c, s)
                for T in (pi - (l + 1) * self.delta, pi - l * self.delta):
                    if self.add(self.jrow(items, [(self.iS + l, 1.0)], e, c, T), 0.0, key=('refl', e, l, T, tuple(sorted(items)))):
                        cuts += 1
        return cuts, worst

    def chord_lines(self):
        out = []
        for g1, g2 in zip(self.chords[:-1], self.chords[1:]):
            s = (math.sqrt(g2) - math.sqrt(g1)) / (g2 - g1)
            out.append((math.sqrt(g1) - s * g1, s))
        return out

    def sep_faces(self, tol=1e-10, per_family=4):
        N, U = self.N, self.U
        pi = math.pi
        Dv = self.x[self.iD:self.iD + N]
        mu, c0, c1 = self.get('mu'), self.get('c0'), self.get('c1')
        cuts, worst = 0, 0.0
        smax = (self.kappa_max - 2) * U
        fams = [('F0', Dv.copy(), None)]
        for i, (al, sl) in enumerate(self.chord_lines()):
            fams.append(('F1', Dv + 2 * mu * sl * self.hmin, (i, al, sl)))
        for (name, w, info) in fams:
            D, CH = dp_table(w, self.kappa_max, smax, maximize=False)
            found = []
            for kap in range(3, self.kappa_max + 1):
                lo, hi = max(0, (kap - 2) * U - kap), (kap - 2) * U
                s = int(np.argmin(D[kap][lo:hi + 1])) + lo
                lin = kap * c0 + c1 * (kap - 2) * pi
                if name == 'F0':
                    viol = -(D[kap][s] + lin)
                else:
                    i, al, sl = info
                    viol = 1 - (D[kap][s] + lin + 2 * mu * (al + sl * self.sig * (kap - 2) * pi))
                worst = max(worst, viol)
                if viol > tol:
                    found.append((viol, kap, s))
            found.sort(reverse=True)
            for (viol, kap, s) in found[:per_family]:
                items = backtrack(CH, kap, s)
                coefs = {self.ic0: -float(kap), self.ic1: -(kap - 2) * pi}
                for k in items:
                    coefs[self.iD + k] = coefs.get(self.iD + k, 0) - 1.0
                if name == 'F0':
                    if self.add(coefs, 0.0, key=('F0', tuple(sorted(items)))):
                        cuts += 1
                else:
                    i, al, sl = info
                    G = float(sum(self.hmin[k] for k in items)) + self.sig * (kap - 2) * pi
                    coefs[self.imu] = -2 * (al + sl * G)
                    if self.add(coefs, -1.0, key=('F1', i, tuple(sorted(items)))):
                        cuts += 1
        c2, w2 = self.sep_nonconvex(tol)
        return cuts + c2, max(worst, w2)

    def sep_nonconvex(self, tol):
        """(N): kappa c0 + c1 (sum theta) + sum D + sum_reflex >= ..., with sum theta = (kappa-2)pi - E,
        reflex part >= beta E >= c1 E; so LHS >= kappa c0 + c1 (kappa-2) pi + sum D, over sum_k <= (kappa-2)U."""
        N, U = self.N, self.U
        pi = math.pi
        Dv = self.x[self.iD:self.iD + N]
        mu, c0, c1 = self.get('mu'), self.get('c0'), self.get('c1')
        cuts, worst = 0, 0.0
        smax = (self.kappa_max - 2) * U
        D, CH = dp_table(Dv, self.kappa_max, smax, maximize=False)
        for kap in range(3, self.kappa_max + 1):
            target = (kap - 2) * U
            s = int(np.argmin(D[kap][:target + 1]))
            v = D[kap][s] + kap * c0 + c1 * (kap - 2) * pi
            for rhs, tag in ((0.0, 'N0'), (1 - mu * self.P[kap], 'N1')):
                viol = rhs - v
                worst = max(worst, viol)
                if viol > tol:
                    items = backtrack(CH, kap, s)
                    coefs = {self.ic0: -float(kap), self.ic1: -(kap - 2) * pi}
                    for k in items:
                        coefs[self.iD + k] = coefs.get(self.iD + k, 0) - 1.0
                    if tag == 'N1':
                        coefs[self.imu] = -self.P[kap]
                        rr = -1.0
                    else:
                        rr = 0.0
                    if self.add(coefs, rr, key=(tag, kap, tuple(sorted(items)))):
                        cuts += 1
        return cuts, worst

    def run(self, max_iter=500, tol=1e-10):
        self.initial()
        hist = []
        val = None
        for it in range(max_iter):
            val = self.solve()
            cj, wj = self.sep_junctions(tol)
            cf, wf = self.sep_faces(tol)
            hist.append((it, val, cj, cf, wj, wf))
            if self.verbose and (it % 10 == 0 or (cj == 0 and cf == 0)):
                print(f'  it {it}: bound {val:.7f} cuts J{cj} F{cf} worst J{wj:.2e} F{wf:.2e} rows {len(self.rows)}', flush=True)
            if cj == 0 and cf == 0:
                break
        self.final_worst = (wj, wf)
        return val, hist


class SampledLP:
    """Exploration (not a certificate): point-sampled angles theta = k*delta (k=1..N-1), exact angle
    sums, convex fields only, junctions fan(m)/conv(e) with exact sums, torus (no outer face).
    Its optimum lies between the best grid-angle fractional tiling and the continuous certificate
    optimum, so it locates where the angle-potential family stalls."""

    def __init__(self, delta_deg=1.0, kappa_max=8, fan_max=6, conv_max=6, gs=None):
        self.N = int(round(180 / delta_deg))
        self.U = self.N
        self.delta = math.pi / self.N
        self.kappa_max, self.fan_max, self.conv_max = kappa_max, fan_max, conv_max
        N = self.N
        # variables: Phi[1..N-1] at indices 0..N-2, mu, rho
        self.nphi = N - 1
        self.imu, self.irho = N - 1, N
        self.nv = N + 1
        self.cot = np.array([1 / math.tan(k * self.delta / 2) for k in range(1, N)])
        if gs is None:
            gs = list(np.linspace(2.9, 5.2, 116))
        self.lines = []
        for g1, g2 in zip(gs[:-1], gs[1:]):
            s = (math.sqrt(g2) - math.sqrt(g1)) / (g2 - g1)
            self.lines.append((math.sqrt(g1) - s * g1, s))
        self.rows, self.keys, self.rowkeys = [], set(), []

    def add(self, coefs, rhs, key):
        if key in self.keys:
            return False
        self.keys.add(key)
        self.rows.append((coefs, rhs))
        self.rowkeys.append(key)
        return True

    def primal(self, thr=1e-9):
        """re-solve and return the fractional tiling: (dual weight, row key) of the active rows."""
        self.solve()
        du = -self.res.ineqlin.marginals
        return sorted([(float(du[i]), self.rowkeys[i]) for i in range(len(du)) if du[i] > thr], reverse=True)

    def solve(self):
        data, ri, ci, b = [], [], [], []
        for r, (coefs, rhs) in enumerate(self.rows):
            for v, c in coefs.items():
                ri.append(r); ci.append(v); data.append(c)
            b.append(rhs)
        A = csr_matrix((data, (ri, ci)), shape=(len(self.rows), self.nv))
        c = np.zeros(self.nv)
        c[self.imu] = 2.0
        c[self.irho] = 2.0
        bounds = [(-10, 10)] * self.nphi + [(0.2236, 1), (0, 1)]
        res = linprog(c, A_ub=A, b_ub=np.array(b), bounds=bounds, method='highs')
        if res.status != 0:
            raise RuntimeError(res.message)
        self.x = res.x
        self.res = res
        return res.fun

    def sep(self, tol=1e-10):
        N, U = self.N, self.U
        phi = self.x[:self.nphi]
        mu, rho = self.x[self.imu], self.x[self.irho]
        # DP over items k=1..N-1 : shift sizes by 1 (item index i=k-1 has size k)
        w_pad = np.concatenate([[NEG], phi])       # size 0 item forbidden
        cuts, worst = 0, 0.0
        cmax = max(self.fan_max + 1, self.conv_max)
        D, CH = dp_table(w_pad, cmax, 2 * U, maximize=True)
        for m in range(1, self.fan_max + 1):
            v = D[m + 1][U] - rho * m
            worst = max(worst, v)
            if v > tol:
                items = backtrack(CH, m + 1, U)
                coefs = {self.irho: -float(m)}
                for k in items:
                    coefs[k - 1] = coefs.get(k - 1, 0) + 1.0
                cuts += self.add(coefs, 0.0, ('fan', tuple(sorted(items))))
        for e in range(3, self.conv_max + 1):
            v = D[e][2 * U] - rho * e
            worst = max(worst, v)
            if v > tol:
                items = backtrack(CH, e, 2 * U)
                coefs = {self.irho: -float(e)}
                for k in items:
                    coefs[k - 1] = coefs.get(k - 1, 0) + 1.0
                cuts += self.add(coefs, 0.0, ('conv', tuple(sorted(items))))
        smax = (self.kappa_max - 2) * U
        wpos = np.concatenate([[POS], phi])
        fams = [(None, wpos)]
        cot_pad = np.concatenate([[0.0], self.cot])
        for (al, sl) in self.lines:
            fams.append(((al, sl), wpos + 2 * mu * sl * cot_pad))
        for info, w in fams:
            D, CH = dp_table(w, self.kappa_max, smax, maximize=False)
            for kap in range(3, self.kappa_max + 1):
                t = (kap - 2) * U
                v = D[kap][t]
                if v >= POS / 2:
                    continue
                if info is None:
                    viol = -v
                else:
                    al, sl = info
                    viol = 1 - v - 2 * mu * al
                worst = max(worst, viol)
                if viol > tol:
                    items = backtrack(CH, kap, t)
                    coefs = {}
                    for k in items:
                        coefs[k - 1] = coefs.get(k - 1, 0) - 1.0
                    if info is None:
                        cuts += self.add(coefs, 0.0, ('F0', tuple(sorted(items))))
                    else:
                        al, sl = info
                        g = float(sum(self.cot[k - 1] for k in items))
                        coefs[self.imu] = -2 * (al + sl * g)
                        cuts += self.add(coefs, -1.0, ('F1', al, tuple(sorted(items))))
        return cuts, worst

    def run(self, max_iter=600, tol=1e-10, verbose=False):
        it = 0
        # seed: T pairs
        for k in range(1, self.U):
            l = self.U - k
            if k <= l:
                coefs = {self.irho: -1.0}
                coefs[k - 1] = coefs.get(k - 1, 0) + 1.0
                coefs[l - 1] = coefs.get(l - 1, 0) + 1.0
                self.add(coefs, 0.0, ('fan', (k, l)))
        for it in range(max_iter):
            val = self.solve()
            cuts, worst = self.sep(tol)
            if verbose and it % 20 == 0:
                print(f'   it {it} val {val:.7f} cuts {cuts} worst {worst:.2e} rows {len(self.rows)}', flush=True)
            if cuts == 0:
                break
        return val, it


# =====================================================================================================
# Rigorous certificate with a piecewise-LINEAR potential (grid nodes k*delta), see the note, section B.
# =====================================================================================================
def cot2(t):
    return math.cos(t / 2) / math.sin(t / 2)


def cot2_dd(t):
    """second derivative of cot(t/2): (1/2) cos(t/2)/sin(t/2)^3 (decreasing on (0, pi))."""
    s = math.sin(t / 2)
    return 0.5 * math.cos(t / 2) / s ** 3


class PLLP:
    """Phi piecewise linear with nodes PHI[k] = Phi(k delta), k = 0..N (k=N is the limit pi-);
    reflex potential R[l] = Phi(pi + l delta), l = 0..N (limits pi+ and 2pi-).  Phi(pi) = 0.
    Junction constraints reduce exactly to node tuples (vertices of cell-hyperplane sections are nodes).
    Face constraints reduce to node tuples with the Taylor correction E (concavity on cells)."""

    def __init__(self, delta_deg=2.0, kappa_max=10, fan_max=10, conv_max=9, refl_max=8,
                 G_A=1.0, G0=5.0, lam=0.8, kmin_deg=22.0, chords=None, verbose=False, nmode='crude'):
        self.nmode = nmode
        self.dd = delta_deg
        self.N = N = int(round(180 / delta_deg))
        assert abs(N * delta_deg - 180) < 1e-9
        self.delta = d = math.pi / N
        self.kappa_max, self.fan_max, self.conv_max, self.refl_max = kappa_max, fan_max, conv_max, refl_max
        self.G_A, self.G0, self.lam = G_A, G0, lam
        self.kmin = int(math.floor(2 * math.atan(1 / G0) / d + 1e-12))
        assert cot2(self.kmin * d) >= G0            # an angle < kmin delta forces g > G0
        if chords is None:
            chords = [G_A, 2.6, 3.0, 3.2, 3.35, 3.45, 3.55, 3.65, 3.75, 3.85, 3.95, 4.05, 4.2, 4.4, 4.7, G0]
        self.chords = chords
        self.lines = []
        for g1, g2 in zip(chords[:-1], chords[1:]):
            s = (math.sqrt(g2) - math.sqrt(g1)) / (g2 - g1)
            self.lines.append((math.sqrt(g1) - s * g1, s))
        self.verbose = verbose
        # corrected cot and h_lambda node values
        self.chat = np.full(N + 1, np.inf)
        for k in range(max(1, self.kmin), N + 1):
            E = 0.5 * (d / 2) ** 2 * cot2_dd((k - 1) * d) if k >= 2 else np.inf
            self.chat[k] = cot2(k * d) - E
        self.theta_star = 2 * math.asin(1 / math.sqrt(2 * lam))
        m_lam = cot2(self.theta_star) - lam * (math.pi - self.theta_star)
        self.m_lam = m_lam
        self.hhat = np.zeros(N + 1)
        for k in range(N + 1):
            t = k * d
            h = m_lam if t <= self.theta_star else cot2(t) - lam * (math.pi - t)
            if (k + 1) * d > self.theta_star:
                E = 0.5 * (d / 2) ** 2 * cot2_dd(max((k - 1) * d, self.theta_star))
            else:
                E = 0.0
            self.hhat[k] = h - E
        # variables
        self.iF = 0
        self.iR = N + 1
        names = ['mu', 'rho', 'a', 'c', 'beta', 'A', 'B', 'Smax', 'gam', 'g1', 'b1', 'g2', 'b2', 'g3', 'b3']
        for i, nm in enumerate(names):
            setattr(self, 'i' + nm, 2 * N + 2 + i)
        self.nv = 2 * N + 2 + len(names)
        self.rows, self.keys, self.rowkeys = [], set(), []
        self.P = [None, None, None] + [2 * math.sqrt(k * math.tan(math.pi / k)) for k in range(3, 400)]

    def add(self, coefs, rhs, key):
        if key in self.keys:
            return False
        self.keys.add(key)
        self.rows.append((coefs, rhs))
        self.rowkeys.append(key)
        return True

    def get(self, nm):
        return self.x[getattr(self, 'i' + nm)]

    def static_rows(self):
        N, d, pi = self.N, self.delta, math.pi
        self.add({self.imu: -1.0}, -1 / (2 * math.sqrt(self.G0)), ('mu_min',))
        if getattr(self, 'bare', False):
            return
        for k in range(N + 1):
            self.add({self.ia: 1.0, self.ic: k * d, self.iF + k: -1.0}, 0.0, ('Omin', k))       # Phi >= a + c theta
            self.add({self.ic: k * d, self.iR + k: -1.0}, 0.0, ('ORmin', k))                   # R >= c (psi - pi)
            self.add({self.igam: 1.0, self.ibeta: k * d, self.iR + k: -1.0}, 0.0, ('NRmin', k))  # R >= gam + beta (psi - pi)
            self.add({self.iF + k: 1.0, self.iA: -1.0, self.iB: -k * d}, 0.0, ('Maj', k))     # Phi <= A + B theta
            self.add({self.iR + k: 1.0, self.iSmax: -1.0}, 0.0, ('MajR', k))
        for (ig, ib, rr) in ((self.ig1, self.ib1, 1), (self.ig2, self.ib2, 2), (self.ig3, self.ib3, 3)):
            for k in range(N + 1):
                self.add({ig: 1.0 / rr, ib: k * d, self.iR + k: -1.0}, 0.0, ('Rmin_r', rr, k))
        self.add({self.ia: -1.0, self.ic: -pi}, 0.0, ('a+cpi>=0',))
        K1 = self.kappa_max + 1
        self.add({self.ia: -K1, self.ic: -(K1 - 2) * pi}, 0.0, ('Ftail0',))
        self.add({self.ia: -K1, self.ic: -(K1 - 2) * pi, self.imu: -2 * math.sqrt(pi)}, -1.0, ('Ftail1',))
        self.add({self.iA: 1.0, self.irho: -1.0}, 0.0, ('A<=rho',))
        M = self.fan_max + 1
        self.add({self.iA: M + 1, self.iB: pi, self.irho: -M}, 0.0, ('Jtail_fan',))
        E = self.conv_max + 1
        self.add({self.iA: E, self.iB: 2 * pi, self.irho: -E}, 0.0, ('Jtail_conv',))
        R = self.refl_max + 1
        self.add({self.iA: R - 1, self.iB: pi, self.iSmax: 1.0, self.irho: -R}, 0.0, ('Jtail_refl',))

    def solve(self):
        data, ri, ci, b = [], [], [], []
        for r, (coefs, rhs) in enumerate(self.rows):
            for v, c in coefs.items():
                ri.append(r); ci.append(v); data.append(c)
            b.append(rhs)
        A = csr_matrix((data, (ri, ci)), shape=(len(self.rows), self.nv))
        c = np.zeros(self.nv)
        c[self.imu] = 2.0
        c[self.irho] = 2.0
        bounds = [(-3, 3)] * (2 * self.N + 2) + [(0, 1), (0, 1), (-3, 3), (0, 3), (0, 3), (-3, 3), (0, 3), (-3, 3), (0, 3)] + [(0, 3)] * 6
        res = linprog(c, A_ub=A, b_ub=np.array(b), bounds=bounds, method='highs')
        if res.status != 0:
            raise RuntimeError(res.message)
        self.x, self.res = res.x, res
        return res.fun

    def _row(self, items, sign, extra=None):
        coefs = {} if extra is None else dict(extra)
        for k in items:
            coefs[self.iF + k] = coefs.get(self.iF + k, 0) + sign
        return coefs

    def sep(self, tol=1e-10, per_family=3):
        N, d, pi = self.N, self.delta, math.pi
        F = self.x[self.iF:self.iF + N + 1]
        R = self.x[self.iR:self.iR + N + 1]
        mu, rho, beta, gam = self.get('mu'), self.get('rho'), self.get('beta'), self.get('gam')
        cuts, worst = 0, 0.0
        # ---- junctions
        cmax = max(self.fan_max + 1, self.conv_max)
        D, CH = dp_table(F, cmax, 2 * N, maximize=True)
        for m in range(1, self.fan_max + 1):
            v = D[m + 1][N] - rho * m
            worst = max(worst, v)
            if v > tol:
                it = backtrack(CH, m + 1, N)
                cuts += self.add(self._row(it, 1.0, {self.irho: -float(m)}), 0.0, ('fan', tuple(sorted(it))))
        for e in range(3, self.conv_max + 1):
            v = D[e][2 * N] - rho * e
            worst = max(worst, v)
            if v > tol:
                it = backtrack(CH, e, 2 * N)
                cuts += self.add(self._row(it, 1.0, {self.irho: -float(e)}), 0.0, ('conv', tuple(sorted(it))))
        for e in range(2, self.refl_max + 1):
            best = (NEG, None)
            for l in range(N + 1):
                v = D[e - 1][N - l] + R[l]
                if v > best[0]:
                    best = (v, l)
            v, l = best
            v -= rho * e
            worst = max(worst, v)
            if v > tol:
                it = backtrack(CH, e - 1, N - l)
                cuts += self.add(self._row(it, 1.0, {self.iR + l: 1.0, self.irho: -float(e)}), 0.0, ('refl', l, tuple(sorted(it))))
        # ---- convex faces
        smax = (self.kappa_max - 2) * N
        fams = [('F0', F.copy(), None)]
        Fm = F.copy()
        Fm[:self.kmin] = POS
        for j, (al, sl) in enumerate(self.lines):
            w = Fm + 2 * mu * sl * np.where(np.isfinite(self.chat), self.chat, 0.0)
            w[:self.kmin] = POS
            fams.append(('F1', w, (j, al, sl)))
        for name, w, info in fams:
            D, CH = dp_table(w, self.kappa_max, smax, maximize=False)
            found = []
            for kap in range(3, self.kappa_max + 1):
                t = (kap - 2) * N
                v = D[kap][t]
                if v >= POS / 2:
                    continue
                viol = -v if name == 'F0' else 1 - v - 2 * mu * info[1]
                worst = max(worst, viol)
                if viol > tol:
                    found.append((viol, kap))
            found.sort(reverse=True)
            for viol, kap in found[:per_family]:
                it = backtrack(CH, kap, (kap - 2) * N)
                if name == 'F0':
                    cuts += self.add(self._row(it, -1.0), 0.0, ('F0', tuple(sorted(it))))
                else:
                    j, al, sl = info
                    G = float(sum(self.chat[k] for k in it))
                    cuts += self.add(self._row(it, -1.0, {self.imu: -2 * (al + sl * G)}), -1.0, ('F1', j, tuple(sorted(it))))
        # ---- non-convex faces (hull + Lagrangian bound): sum_k (F_k - beta k d) + beta (kap-2) pi [+ 2 mu L(G)]
        if self.nmode == 'hull':
            c2, w2 = sep_hull(self, tol)
            return cuts + c2, max(worst, w2)
        if self.nmode == 'pocket2':
            c2, w2 = sep_pocket2(self, tol)
            return cuts + c2, max(worst, w2)
        if self.nmode == 'pocket':
            if not hasattr(self, '_pocket'):
                self._pocket = pocket_setup(self, self.lam)
            c2, w2 = sep_pocket(self, self.lam, tol)
            return cuts + c2, max(worst, w2)
        base = F - beta * d * np.arange(N + 1)
        fams = [('N0', base.copy(), None)] if self.nmode != 'none' else []
        if self.nmode == 'none':
            pass
        elif self.nmode == 'lagr':
            for j, (al, sl) in enumerate(self.lines):
                fams.append(('N1', base + 2 * mu * sl * self.hhat, (j, al, sl)))
        else:
            fams.append(('NC', base.copy(), None))
        for name, w, info in fams:
            D, CH = dp_table(w, self.kappa_max, smax, maximize=False)
            found = []
            for kap in range(3, self.kappa_max + 1):
                t = (kap - 2) * N
                s = int(np.argmin(D[kap][:t + 1]))
                v = D[kap][s] + beta * (kap - 2) * pi + (gam if self.nmode != 'nogam' else 0.0)
                if name == 'N0':
                    viol = -v
                elif name == 'NC':
                    viol = 1 - mu * self.P[kap] - v
                else:
                    j, al, sl = info
                    viol = 1 - v - 2 * mu * (al + sl * 2 * pi * self.lam)
                worst = max(worst, viol)
                if viol > tol:
                    found.append((viol, kap, s))
            found.sort(reverse=True)
            for viol, kap, s in found[:per_family]:
                it = backtrack(CH, kap, s)
                ex = {self.ibeta: -(kap - 2) * pi + d * sum(it), self.igam: -1.0}
                if name == 'N0':
                    cuts += self.add(self._row(it, -1.0, ex), 0.0, ('N0', tuple(sorted(it))))
                elif name == 'NC':
                    ex[self.imu] = -self.P[kap]
                    cuts += self.add(self._row(it, -1.0, ex), -1.0, ('NC', kap, tuple(sorted(it))))
                else:
                    j, al, sl = info
                    G = float(sum(self.hhat[k] for k in it)) + 2 * pi * self.lam
                    ex[self.imu] = -2 * (al + sl * G)
                    cuts += self.add(self._row(it, -1.0, ex), -1.0, ('N1', j, tuple(sorted(it))))
        return cuts, worst

    def run(self, max_iter=800, tol=1e-10):
        self.static_rows()
        N = self.N
        for k in range(N + 1):   # seed T rows
            l = N - k
            if k <= l:
                self.add(self._row([k, l], 1.0, {self.irho: -1.0}), 0.0, ('fan', tuple(sorted([k, l]))))
        val = None
        for it in range(max_iter):
            val = self.solve()
            cuts, worst = self.sep(tol)
            if self.verbose and it % 25 == 0:
                print(f'   it {it} val {val:.7f} cuts {cuts} worst {worst:.2e} rows {len(self.rows)}', flush=True)
            if cuts == 0:
                break
        self.iters = it
        self.final_worst = worst
        return val


def dp3(wtypes, flexflag, cmax, smax, fmax):
    """min-plus DP over corners with a type per corner: D[c][s][f] = min sum of weights of c corners,
    node-size sum s, f = number of 'flexible' corners.  wtypes: list of arrays (one per type)."""
    D = [np.full((smax + 1, fmax + 1), POS)]
    D[0][0, 0] = 0.0
    CH = [None]
    for c in range(1, cmax + 1):
        prev = D[-1]
        cur = np.full((smax + 1, fmax + 1), POS)
        chk = np.full((smax + 1, fmax + 1), -1, dtype=np.int32)
        cht = np.full((smax + 1, fmax + 1), -1, dtype=np.int8)
        for ti, w in enumerate(wtypes):
            df = 1 if flexflag[ti] else 0
            for k in range(min(len(w), smax + 1)):
                if not np.isfinite(w[k]) or w[k] >= POS / 2:
                    continue
                cand = prev[:smax + 1 - k, :fmax + 1 - df] + w[k]
                tgt = cur[k:, df:]
                better = cand < tgt
                tgt[better] = cand[better]
                chk[k:, df:][better] = k
                cht[k:, df:][better] = ti
        D.append(cur)
        CH.append((chk, cht))
    return D, CH


def back3(CH, flexflag, c, s, f):
    items = []
    while c > 0:
        chk, cht = CH[c]
        k, ti = int(chk[s, f]), int(cht[s, f])
        items.append((k, ti))
        s -= k
        f -= 1 if flexflag[ti] else 0
        c -= 1
    return items


def pocket_setup(L, lam):
    """node values (with Taylor corrections) of the three corner types for the pocket bound."""
    N, d = L.N, L.delta
    ts = 2 * math.asin(1 / math.sqrt(2 * lam))
    m = cot2(ts) - lam * (math.pi - ts)
    ve = np.full(N + 1, np.inf)
    vh = np.zeros(N + 1)
    for k in range(N + 1):
        t = k * d
        if k >= max(2, L.kmin):
            Ee = 0.5 * (d / 2) ** 2 * cot2_dd((k - 1) * d)
            ve[k] = cot2(t) - lam * (math.pi - t) - Ee
        h = m if t <= ts else cot2(t) - lam * (math.pi - t)
        Eh = 0.5 * (d / 2) ** 2 * cot2_dd(max((k - 1) * d, ts)) if (k + 1) * d > ts else 0.0
        vh[k] = h - Eh
    return ve, vh, np.zeros(N + 1)


def reflex_tables(R, N, rmax, tmax):
    """Rex[r][t] = min sum R over exactly r reflex nodes with node sum t (r < rmax);
    Rall[t] = min over any r >= 1 (unbounded knapsack), with choices for backtracking."""
    Rex = [np.full(tmax + 1, POS)]
    Rex[0][0] = 0.0
    CHx = [None]
    for r in range(1, rmax):
        prev = Rex[-1]
        cur = np.full(tmax + 1, POS)
        ch = np.full(tmax + 1, -1, dtype=np.int32)
        for l in range(min(N, tmax) + 1):
            cand = prev[:tmax + 1 - l] + R[l]
            tgt = cur[l:]
            better = cand < tgt
            tgt[better] = cand[better]
            ch[l:][better] = l
        Rex.append(cur)
        CHx.append(ch)
    # unbounded: A[t] = min over >=1 items
    A = np.full(tmax + 1, POS)
    cha = np.full(tmax + 1, -1, dtype=np.int32)
    for t in range(tmax + 1):
        best, arg = POS, -1
        for l in range(0, min(N, t) + 1):
            rest = 0.0 if t - l == 0 else A[t - l]
            if l == 0 and t == 0:
                v = R[0]
            elif l == 0:
                continue  # zero-excess items only help if R[0] < 0 (excluded: R >= 0)
            else:
                v = R[l] + (0.0 if t - l == 0 else min(A[t - l], POS))
            if v < best:
                best, arg = v, l
        A[t], cha[t] = best, arg
    return Rex, CHx, A, cha


def sep_pocket(L, lam, tol=1e-10, per_family=3):
    """non-convex fields (r >= 1 reflex corners): hull + pocket bound, see note section B."""
    N, d, pi = L.N, L.delta, math.pi
    F = L.x[L.iF:L.iF + N + 1]
    R = np.maximum(L.x[L.iR:L.iR + N + 1], 0.0)   # rows R >= gam + beta l d >= 0 keep R >= 0
    mu = L.get('mu')
    ve, vh, v0 = L._pocket
    km = L.kappa_max
    smax = (km - 2) * N
    rmax = (km + 1) // 2
    Rex, CHx, Rall, cha = reflex_tables(R, N, rmax, smax)
    flexflag = [False, True, False]
    cuts, worst = 0, 0.0
    fams = [('P0', None)] + [('P1', (j, al, sl)) for j, (al, sl) in enumerate(L.lines)]
    for name, info in fams:
        if name == 'P0':
            wt = [F.copy(), np.full(N + 1, POS), np.full(N + 1, POS)]
            const = 0.0
        else:
            j, al, sl = info
            wt = [F + 2 * mu * sl * ve, F + 2 * mu * sl * vh, F + 2 * mu * sl * v0]
            wt = [np.where(np.isfinite(w), w, POS) for w in wt]
            const = 2 * mu * (al + sl * 2 * pi * lam)
        D, CH = dp3(wt, flexflag, km, smax, km)
        found = []
        for kap in range(3, km + 1):
            target = (kap - 2) * N
            Dk = D[kap]
            cum = np.minimum.accumulate(Dk, axis=1)       # min over f <= F
            for r in range(1, rmax + 1):
                Fa = min(2 * r, kap)
                col = cum[:target + 1, Fa]
                if r < rmax:
                    Rt = Rex[r][:target + 1][::-1]
                else:
                    Rt = Rall[:target + 1][::-1]
                tot = col + Rt
                s = int(np.argmin(tot))
                v = tot[s] + const
                viol = (0.0 - v) if name == 'P0' else (1.0 - v)
                worst = max(worst, viol)
                if viol > tol:
                    found.append((viol, kap, r, s, Fa))
        found.sort(reverse=True)
        for viol, kap, r, s, Fa in found[:per_family]:
            target = (kap - 2) * N
            f = int(np.argmin(D[kap][s, :Fa + 1]))
            items = back3(CH, flexflag, kap, s, f)
            t = target - s
            ritems = []
            if r < rmax:
                rr, tt = r, t
                while rr > 0:
                    l = int(CHx[rr][tt]); ritems.append(l); tt -= l; rr -= 1
            else:
                tt = t
                while True:
                    l = int(cha[tt]); ritems.append(l); tt -= l
                    if tt <= 0:
                        break
            coefs = {}
            for (k, ti) in items:
                coefs[L.iF + k] = coefs.get(L.iF + k, 0) - 1.0
            for l in ritems:
                coefs[L.iR + l] = coefs.get(L.iR + l, 0) - 1.0
            if name == 'P0':
                cuts += L.add(coefs, 0.0, ('P0', tuple(sorted(items)), tuple(sorted(ritems))))
            else:
                j, al, sl = info
                G = sum((ve, vh, v0)[ti][k] for (k, ti) in items) + 2 * pi * lam
                coefs[L.imu] = -2 * (al + sl * G)
                cuts += L.add(coefs, -1.0, ('P1', j, tuple(sorted(items)), tuple(sorted(ritems))))
    return cuts, worst


def certify_convex(L, eps=2e-9):
    """Re-verify a PLLP solution (nmode='none': fields with convex outer boundary) after a safety
    bump: Phi += eps, A += eps, mu += 4 eps, rho += (fan_max+2) eps.  Float64 dynamic programming with
    margin; every check must pass with slack >= 0 (the float error of these sums is < 1e-13).
    Returns (certified bound 2(mu+rho), dict of worst slacks)."""
    N, d, pi = L.N, L.delta, math.pi
    F = L.x[L.iF:L.iF + N + 1] + eps
    R = L.x[L.iR:L.iR + N + 1] + eps
    mu = L.get('mu') + 4 * eps
    rho = L.get('rho') + (max(L.fan_max, L.conv_max, L.refl_max) + 3) * eps
    a, c = L.get('a'), L.get('c')
    A, B, Smax = L.get('A') + 2 * eps, L.get('B'), L.get('Smax') + 2 * eps
    slack = {}
    # junctions
    cmax = max(L.fan_max + 1, L.conv_max)
    D, _ = dp_table(F, cmax, 2 * N, maximize=True)
    slack['fan'] = min(rho * m - D[m + 1][N] for m in range(1, L.fan_max + 1))
    slack['conv'] = min(rho * e - D[e][2 * N] for e in range(3, L.conv_max + 1))
    slack['refl'] = min(rho * e - max(D[e - 1][N - l] + R[l] for l in range(N + 1)) for e in range(2, L.refl_max + 1))
    # tails and linear bounds (at nodes; piecewise linear => everywhere)
    slack['major'] = min(min(A + B * k * d - F[k] for k in range(N + 1)), min(Smax - R[l] for l in range(N + 1)))
    slack['A<=rho'] = rho - A
    slack['Jtail'] = min(rho * (L.fan_max + 1) - (L.fan_max + 2) * A - B * pi,
                         rho * (L.conv_max + 1) - (L.conv_max + 1) * A - 2 * pi * B,
                         rho * (L.refl_max + 1) - L.refl_max * A - B * pi - Smax)
    slack['O'] = min(min(F[k] - a - c * k * d for k in range(N + 1)), min(R[l] - c * l * d for l in range(N + 1)), a + c * pi, c)
    K1 = L.kappa_max + 1
    slack['Ftail'] = min(K1 * a + (K1 - 2) * pi * c, K1 * a + (K1 - 2) * pi * c - (1 - 2 * mu * math.sqrt(pi)))
    slack['mu_min'] = mu - 1 / (2 * math.sqrt(L.G0))
    # faces
    smax = (L.kappa_max - 2) * N
    D0, _ = dp_table(F, L.kappa_max, smax, maximize=False)
    slack['F0'] = min(D0[k][(k - 2) * N] for k in range(3, L.kappa_max + 1))
    worst = np.inf
    Fm = F.copy()
    Fm[:L.kmin] = POS
    for (al, sl) in L.lines:
        w = Fm + 2 * mu * sl * np.where(np.isfinite(L.chat), L.chat, 0.0)
        w[:L.kmin] = POS
        Dj, _ = dp_table(w, L.kappa_max, smax, maximize=False)
        for k in range(3, L.kappa_max + 1):
            v = Dj[k][(k - 2) * N]
            if v < POS / 2:
                worst = min(worst, v + 2 * mu * al - 1)
    slack['F1'] = worst
    ok = all(v >= 0 for v in slack.values())
    return 2 * (mu + rho), slack, ok


# ---------------------------------------------------------------------------------------------
# Non-convex fields, budget-dual pocket bound (see note B.4):  for a field with convex corners
# theta_v and reflex excess E = (kappa-2) pi - sum theta_v, the hull satisfies
#   g_hull >= G_lam = sum_v f_role(v)(theta_v) - lam (kappa-2) pi   for every lam >= 0,
# roles: exact hull vertex  f = cot(t/2) + lam t;  pocket end (<= 2r of them)  f = min_{u in [t,pi]} cot(u/2) + lam u;
# pocket interior  f = lam pi.   lam may depend on the block of the convex node sum s.
# ---------------------------------------------------------------------------------------------
def role_values(L, lam):
    N, d = L.N, L.delta
    ts = 2 * math.asin(1 / math.sqrt(2 * lam)) if lam >= 0.5 else math.pi
    ve = np.full(N + 1, np.inf)
    vp = np.zeros(N + 1)
    vi = np.full(N + 1, lam * math.pi)
    for k in range(N + 1):
        t = k * d
        if k >= max(2, L.kmin):
            ve[k] = cot2(t) + lam * t - 0.5 * (d / 2) ** 2 * cot2_dd((k - 1) * d)
        u = max(t, ts)
        val = (cot2(u) if u < math.pi else 0.0) + lam * u
        if (k + 1) * d > ts and k >= 1:
            E = 0.5 * (d / 2) ** 2 * cot2_dd(max((k - 1) * d, ts))
        else:
            E = 0.0
        vp[k] = val - E
    return ve, vp, vi


def face_G(roles, ks, r, lam, kap):
    """min over role assignments (at most 2r pocket ends) of sum role values - lam (kap-2) pi."""
    ve, vp, vi = roles
    base, gains = 0.0, []
    for k in ks:
        m = min(ve[k], vi[k])
        base += m
        gains.append(max(0.0, m - vp[k]))
    gains.sort(reverse=True)
    return base - sum(gains[:2 * r]) - lam * (kap - 2) * math.pi


def sep_pocket2(L, tol=1e-10, per_family=3, lams=(0.05, 0.12, 0.25, 0.4, 0.55, 0.7, 0.9, 1.2, 1.7, 2.5, 4.0), block=None):
    N, d, pi = L.N, L.delta, math.pi
    F = L.x[L.iF:L.iF + N + 1]
    R = np.maximum(L.x[L.iR:L.iR + N + 1], 0.0)
    mu = L.get('mu')
    km = L.kappa_max
    smax = (km - 2) * N
    rmax = (km + 1) // 2
    if block is None:
        block = max(2, N // 6)
    Rex, CHx, Rall, cha = reflex_tables(R, N, rmax, smax)
    flexflag = [False, True, False]
    if not hasattr(L, '_roles'):
        L._roles = {lam: role_values(L, lam) for lam in lams}
    cuts, worst = 0, 0.0
    # N0: sum Phi + sum R >= 0 (roles irrelevant)
    D0, CH0 = dp_table(F, km, smax, maximize=False)
    for kap in range(3, km + 1):
        target = (kap - 2) * N
        tot = D0[kap][:target + 1] + Rall[:target + 1][::-1]
        s = int(np.argmin(tot))
        viol = -tot[s]
        worst = max(worst, viol)
        if viol > tol:
            it = backtrack(CH0, kap, s)
            t = target - s
            ritems = []
            tt = t
            while tt > 0:
                l = int(cha[tt]); ritems.append(l); tt -= l
            if t == 0:
                ritems = [0]
            coefs = {}
            for k in it:
                coefs[L.iF + k] = coefs.get(L.iF + k, 0) - 1.0
            for l in ritems:
                coefs[L.iR + l] = coefs.get(L.iR + l, 0) - 1.0
            cuts += L.add(coefs, 0.0, ('Q0', tuple(sorted(it)), tuple(sorted(ritems))))
    # chords
    for j, (al, sl) in enumerate(L.lines):
        tabs = {}
        for lam in lams:
            ve, vp, vi = L._roles[lam]
            wt = [F + 2 * mu * sl * ve, F + 2 * mu * sl * vp, F + 2 * mu * sl * vi]
            wt = [np.where(np.isfinite(w), w, POS) for w in wt]
            tabs[lam] = dp3(wt, flexflag, km, smax, km)
        for kap in range(3, km + 1):
            target = (kap - 2) * N
            for r in range(1, rmax + 1):
                Fa = min(2 * r, kap)
                Rt = (Rex[r] if r < rmax else Rall)[:target + 1][::-1]
                vals = {}
                for lam in lams:
                    D, CH = tabs[lam]
                    cum = np.minimum.accumulate(D[kap], axis=1)[:target + 1, Fa]
                    vals[lam] = cum + Rt + 2 * mu * (al - sl * lam * (kap - 2) * pi)
                # blocks of the cell base sum s0; vertices have s in [s0, s0 + kap]
                worst_blk = None
                for b0 in range(0, target + 1, block):
                    lo, hi = b0, min(target, b0 + block - 1 + kap)
                    best_lam, best_v, best_s = None, -np.inf, None
                    for lam in lams:
                        seg = vals[lam][lo:hi + 1]
                        s_rel = int(np.argmin(seg))
                        if seg[s_rel] > best_v:
                            best_lam, best_v, best_s = lam, seg[s_rel], lo + s_rel
                    viol = 1.0 - best_v
                    worst = max(worst, viol)
                    if viol > tol and (worst_blk is None or viol > worst_blk[0]):
                        worst_blk = (viol, best_lam, best_s)
                if worst_blk is None:
                    continue
                viol, lam, s = worst_blk
                D, CH = tabs[lam]
                f = int(np.argmin(D[kap][s, :Fa + 1]))
                if D[kap][s, f] >= POS / 2:
                    continue
                items = back3(CH, flexflag, kap, s, f)
                t = target - s
                ritems = []
                if r < rmax:
                    rr, tt = r, t
                    while rr > 0:
                        l = int(CHx[rr][tt]); ritems.append(l); tt -= l; rr -= 1
                else:
                    tt = t
                    while tt > 0:
                        l = int(cha[tt]); ritems.append(l); tt -= l
                    if t == 0:
                        ritems = [0]
                ks = [k for (k, ti) in items]
                G = max(math.pi, max(face_G(L._roles[lm], ks, r, lm, kap) for lm in lams))
                coefs = {}
                for k in ks:
                    coefs[L.iF + k] = coefs.get(L.iF + k, 0) - 1.0
                for l in ritems:
                    coefs[L.iR + l] = coefs.get(L.iR + l, 0) - 1.0
                coefs[L.imu] = -2 * (al + sl * G)
                cuts += L.add(coefs, -1.0, ('Q1', j, tuple(sorted(ks)), tuple(sorted(ritems))))
    return cuts, worst


# ---------------------------------------------------------------------------------------------
# Non-convex fields, exact hull-angle formulation (note B.4).  A field with convex corners theta_v,
# r >= 1 reflex corners (excess E) has hull angles phi_v: phi_v = theta_v (hull vertex that is not
# a pocket end), phi_v in [theta_v, pi) (pocket end; at most 2r of them), phi_v = pi (corner inside
# a pocket); sum phi_v = (kappa-2) pi, E = sum (phi_v - theta_v), g_hull = sum cot(phi_v/2).
# With Phi_reflex(psi) >= gam + beta (psi - pi):  sum Phi + sum Phi_reflex >= sum Phi(theta_v)
# + r gam + beta sum (phi_v - theta_v).  The constraint matrix (phi-theta differences + one sum
# row + boxes) is totally unimodular, so node tuples suffice.
# ---------------------------------------------------------------------------------------------
def sep_hull(L, tol=1e-10):
    cuts, worst = 0, 0.0
    km = L.kappa_max
    rmax = (km + 1) // 2
    for (ig, ib, rs) in ((L.ig1, L.ib1, [1]), (L.ig2, L.ib2, [2]), (L.ig3, L.ib3, list(range(3, rmax + 1)))):
        c, w = _sep_hull_r(L, tol, ig, ib, rs)
        cuts += c
        worst = max(worst, w)
    return cuts, worst


def _sep_hull_r(L, tol, ig, ib, rlist):
    N, d, pi = L.N, L.delta, math.pi
    F = L.x[L.iF:L.iF + N + 1]
    mu, beta, gam0 = L.get('mu'), L.x[ib], L.x[ig]
    km = L.kappa_max
    target_max = (km - 2) * N
    ks = np.arange(N + 1)
    base = F - beta * d * ks                        # Phi_k - beta k d
    pre = np.minimum.accumulate(base)               # min_{k <= k'} (Phi_k - beta k d)
    argpre = np.zeros(N + 1, dtype=int)
    for kp in range(1, N + 1):
        argpre[kp] = argpre[kp - 1] if base[argpre[kp - 1]] <= base[kp] else kp
    kI = int(np.argmin(base))
    flexflag = [False, True, False]
    cuts, worst = 0, 0.0
    chat = np.where(np.isfinite(L.chat), L.chat, 0.0)
    fams = [('H0', None)] + [('H1', (j, al, sl)) for j, (al, sl) in enumerate(L.lines)]
    for name, info in fams:
        wE = F.copy()
        wP = pre + beta * d * ks
        wI = np.full(N + 1, POS)
        wI[N] = base[kI] + beta * d * N
        if name == 'H1':
            j, al, sl = info
            wE = wE + 2 * mu * sl * chat
            wP = wP + 2 * mu * sl * chat
            wE[:L.kmin] = POS
            wP[:L.kmin] = POS
        D, CH = dp3([wE, wP, wI], flexflag, km, target_max, km)
        for kap in range(3, km + 1):
            t = (kap - 2) * N
            cum = np.minimum.accumulate(D[kap][t, :])
            for r in rlist:
                Fa = min(2 * r, kap)
                gcoef = 1.0 if r <= 3 else r / 3.0      # r >= 3 corners carry at least (r/3) g3
                v = cum[Fa] + gcoef * gam0
                if name == 'H0':
                    viol = -v
                else:
                    viol = 1 - v - 2 * mu * info[1]
                worst = max(worst, viol)
                if viol > tol:
                    f = int(np.argmin(D[kap][t, :Fa + 1]))
                    items = back3(CH, flexflag, kap, t, f)
                    coefs = {ig: -gcoef}
                    Esum = 0.0
                    G = 0.0
                    full = []
                    for (kp, ti) in items:
                        if ti == 0:
                            k = kp
                        elif ti == 1:
                            k = int(argpre[kp])
                        else:
                            k = kI
                        coefs[L.iF + k] = coefs.get(L.iF + k, 0) - 1.0
                        full.append((k, kp, ti))
                        Esum += (kp - k) * d
                        if ti != 2 and name == 'H1':
                            G += chat[kp]
                    coefs[ib] = coefs.get(ib, 0) - Esum
                    if name == 'H0':
                        cuts += L.add(coefs, 0.0, ('H0', ig, r, tuple(sorted(full))))
                    else:
                        coefs[L.imu] = -2 * (info[1] + info[2] * G)
                        cuts += L.add(coefs, -1.0, ('H1', ig, info[0], r, tuple(sorted(full))))
    return cuts, worst
