"""procgen_fence2_20260926_lp.py -- the typed piece-potential LPs of lane "fence2" (Part B),
session collatz-procgen-20260922, 2026-09-26.

Model (EMPIRICAL relaxation): periodic configurations (torus), no straight joins; fields convex
(corner types O / I / E) or, with alphabet 'OIER', with reflex corners R.  Every face side is one
fence-side piece.  A certificate (dual solution) consists of
  P1(l, th)   potential of a 1a piece of length l ending at an O corner of angle th
              (= potential of a 1b piece of length l starting at an I corner of angle th: mirror);
  0 pieces    TypedLP:  P0(l) + T0(th_I) + T0(th_O)   (separable)
              TypedLP2: P0N(l, th_I, th_O), symmetric in the two angles (mirror)
  P2          potential of a whole piece (length exactly 1);
  PE(th), PR(psi)  potentials of E corners and reflex corners;  rho per fence end;  beta per side.
Constraints:
  fields   a <= sum of the potentials of its sides and E/R corners      (all typed polygons)
  sides    sum over the pieces of a fence side (and the E sectors of fans at its landings)
           <= beta + rho * (#stems landing on it)                           (all partitions of [0,1])
  t=0      sum PE (+ PR) over the sectors <= rho e                        (conv(e), L, refl(e))
Weak duality: density <= 2 (beta + rho).
Potentials are grid functions (l step 1/NL, angle step pi/NA) with multilinear interpolation.  The
side and junction constraints are encoded exactly in the LP by auxiliary DP variables (vertex
argument); grid_check re-verifies them independently.  Field constraints are separated by multistart
SLSQP over typed polygons (price_pattern) with a random library, L1 centering and row purging: the
field part is heuristic, hence EMPIRICAL.  export_cert / load_cert / verify_cert / tail_check handle
certificates (JSON).
"""
import itertools
import math

import numpy as np
from scipy.optimize import linprog, minimize
from scipy.sparse import coo_matrix

PI = math.pi


def patterns(kmin=3, kmax=7, alphabet='OIE', rmax=0, lmax=8, kmax_r=None):
    """cyclic words over O,I,E (and R: reflex) up to rotation and mirror (reverse + swap O<->I);
    kmin..kmax convex letters and at most rmax letters R, total length <= lmax."""
    sw = {'O': 'I', 'I': 'O', 'E': 'E', 'R': 'R'}
    seen = set()
    out = []
    conv = ''.join(c for c in alphabet if c != 'R')
    words = []
    for k in range(kmin, kmax + 1):
        for w in itertools.product(conv, repeat=k):
            words.append(w)
    if 'R' in alphabet:
        base = list(words)
        for w in base:
            if kmax_r is not None and len(w) > kmax_r:
                continue
            for r in range(1, rmax + 1):
                if len(w) + r > lmax:
                    continue
                for pos in itertools.combinations_with_replacement(range(len(w)), r):
                    ww = list(w)
                    for p in sorted(pos, reverse=True):
                        ww.insert(p + 1, 'R')
                    words.append(tuple(ww))
    for w in words:
        if True:
            k = len(w)
            rots = [w[i:] + w[:i] for i in range(k)]
            mir = tuple(sw[c] for c in reversed(w))
            rots += [mir[i:] + mir[:i] for i in range(k)]
            can = min(rots)
            if can in seen:
                continue
            seen.add(can)
            out.append(can)
    return out


def side_types(pat):
    k = len(pat)
    o = {'O': 1, 'I': 0, 'E': 1, 'R': 1}
    io = {'O': 0, 'I': 1, 'E': 1, 'R': 1}
    st = []
    for i in range(k):
        a, b = pat[i], pat[(i + 1) % k]
        st.append({(1, 1): '2', (1, 0): '1a', (0, 1): '1b', (0, 0): '0'}[(o[a], io[b])])
    return st


class TypedLP:
    def __init__(self, NL=20, NA=18, kmax=7, fan_max=5, conv_max=8, bound=3.0, seed=0, chain_rows=12, alphabet='OIE', rmax=0, refl_max=4, kmax_r=None):
        self.NL, self.NA = NL, NA
        self.chain_rows = chain_rows
        self.kmax, self.fan_max, self.conv_max = kmax, fan_max, conv_max
        idx = 0
        self.iP1 = np.arange(idx, idx + (NL + 1) * (NA + 1)).reshape(NL + 1, NA + 1)
        idx += (NL + 1) * (NA + 1)
        self.iP0 = np.arange(idx, idx + NL + 1)
        idx += NL + 1
        self.iT0 = np.arange(idx, idx + NA + 1)
        idx += NA + 1
        self.iE = np.arange(idx, idx + NA + 1)
        idx += NA + 1
        self.iR = np.arange(idx, idx + NA + 1)    # reflex potential at pi + j*pi/NA
        idx += NA + 1
        self.iP2 = idx
        self.irho = idx + 1
        self.ibeta = idx + 2
        idx += 3
        self.cmax = max(conv_max, fan_max - 1)
        self.iG = np.arange(idx, idx + self.cmax * (2 * NA + 1)).reshape(self.cmax, 2 * NA + 1)
        idx += self.cmax * (2 * NA + 1)          # iG[c-1, s] = G[c][s]
        self.iK = np.arange(idx, idx + NA + 1)
        idx += NA + 1
        self.iU = np.arange(idx, idx + NA + 1)
        idx += NA + 1
        self.iAS = np.arange(idx, idx + NL + 1)
        idx += NL + 1
        self.iAM = idx
        idx += 1
        self.iZ = np.arange(idx, idx + NL + 1)
        idx += NL + 1
        self.nv = idx
        self.bound = bound
        self.rows = []      # (cols, vals, rhs, kind)
        self.keys = set()
        self.rng = np.random.default_rng(seed)
        self.rmax, self.refl_max = rmax, refl_max
        self.kmax_r = kmax_r
        self.pats = patterns(3, kmax, alphabet, rmax, kmax_r=kmax_r)
        self.alphabet = alphabet
        self.pool = {}      # pattern -> list of x vectors (warm starts)
        self.x = None
        self.static_rows()

    # ---------------------------------------------------------------- rows
    def add(self, coef, rhs, kind, key=None):
        if key is not None:
            if key in self.keys:
                return False
            self.keys.add(key)
        cols = np.array(list(coef.keys()), dtype=int)
        vals = np.array([coef[c] for c in cols], dtype=float)
        self.rows.append((cols, vals, float(rhs), kind))
        return True

    @staticmethod
    def acc(coef, var, w):
        if w != 0.0:
            coef[var] = coef.get(var, 0.0) + w

    def static_rows(self):
        """exact LP encoding (auxiliary DP variables) of all fence-side, fan and t=0 junction
        constraints on the grid; only the field constraints need separation."""
        NL, NA, cm = self.NL, self.NA, self.cmax
        A = self.acc

        def G(c, s):
            return self.iG[c - 1, s]
        # G[1][s] >= PE[s]
        for s in range(NA + 1):
            c = {}
            A(c, self.iE[s], 1.0)
            A(c, G(1, s), -1.0)
            self.add(c, 0.0, 'aux')
        for cc in range(2, cm + 1):
            for s in range(2 * NA + 1):
                for j in range(0, min(NA, s) + 1):
                    if s - j > (cc - 1) * NA:
                        continue
                    c = {}
                    A(c, G(cc - 1, s - j), 1.0)
                    A(c, self.iE[j], 1.0)
                    A(c, G(cc, s), -1.0)
                    self.add(c, 0.0, 'aux')
        # conv(e): G[e][2NA] <= e rho
        for e in range(3, self.conv_max + 1):
            c = {}
            A(c, G(e, 2 * NA), 1.0)
            A(c, self.irho, -float(e))
            self.add(c, 0.0, 'conv')
        # L and refl(e): (e-1) E sectors with sum s*delta and one reflex sector pi + (NA-s)*delta
        if 'R' in self.alphabet:
            for e in range(2, self.refl_max + 1):
                for s in range(NA + 1):
                    c = {}
                    A(c, G(e - 1, s), 1.0)
                    A(c, self.iR[NA - s], 1.0)
                    A(c, self.irho, -float(e))
                    self.add(c, 0.0, 'refl')
        # K[s] >= G[c][s] - c rho (c = 1..fan_max-1), K[0] >= 0
        c = {}
        A(c, self.iK[0], -1.0)
        self.add(c, 0.0, 'aux')
        for cc in range(1, self.fan_max):
            for s in range(NA + 1):
                c = {}
                A(c, G(cc, s), 1.0)
                A(c, self.irho, -float(cc))
                A(c, self.iK[s], -1.0)
                self.add(c, 0.0, 'aux')
        if self.side_rows_hook():
            return
        # U[r] >= K[r - jI] + T0[jI]
        for r in range(NA + 1):
            for jI in range(r + 1):
                c = {}
                A(c, self.iK[r - jI], 1.0)
                A(c, self.iT0[jI], 1.0)
                A(c, self.iU[r], -1.0)
                self.add(c, 0.0, 'aux')
        # AS[i] >= P1[i, jO] + U[NA - jO] - rho ;  AM >= T0[jO] + U[NA - jO] - rho
        for jO in range(NA + 1):
            for i in range(NL + 1):
                c = {}
                A(c, self.iP1[i, jO], 1.0)
                A(c, self.iU[NA - jO], 1.0)
                A(c, self.irho, -1.0)
                A(c, self.iAS[i], -1.0)
                self.add(c, 0.0, 'aux')
            c = {}
            A(c, self.iT0[jO], 1.0)
            A(c, self.iU[NA - jO], 1.0)
            A(c, self.irho, -1.0)
            A(c, self.iAM, -1.0)
            self.add(c, 0.0, 'aux')
        # k = 1 fence sides (T, X sides and fans)
        for i in range(NL + 1):
            for jO in range(NA + 1):
                for jI in range(NA + 1 - jO):
                    c = {}
                    A(c, self.iP1[i, jO], 1.0)
                    A(c, self.iK[NA - jO - jI], 1.0)
                    A(c, self.iP1[NL - i, jI], 1.0)
                    A(c, self.irho, -1.0)
                    A(c, self.ibeta, -1.0)
                    self.add(c, 0.0, 'side1')
        # chains k >= 2
        for i0 in range(NL + 1):
            for i1 in range(NL + 1 - i0):
                c = {}
                A(c, self.iAS[i0], 1.0)
                A(c, self.iP0[i1], 1.0)
                A(c, self.iZ[i0 + i1], -1.0)
                self.add(c, 0.0, 'aux')
        for s in range(NL + 1):
            for t in range(NL + 1 - s):
                c = {}
                A(c, self.iZ[s], 1.0)
                A(c, self.iAM, 1.0)
                A(c, self.iP0[t], 1.0)
                A(c, self.iZ[s + t], -1.0)
                self.add(c, 0.0, 'aux')
            c = {}
            A(c, self.iZ[s], 1.0)
            A(c, self.iAS[NL - s], 1.0)
            A(c, self.ibeta, -1.0)
            self.add(c, 0.0, 'chain')
        c = {}
        A(c, self.iP2, 1.0)
        A(c, self.ibeta, -1.0)
        self.add(c, 0.0, 'side0')
        # the unit square with E corners (grid) keeps the LP bounded
        if 'E' in self.alphabet:
            self.add_field(('E', 'E', 'E', 'E'), [PI / 2] * 4, [1.0] * 4)
        else:
            # a near-unit O-pinwheel square keeps the O/I-only LP bounded
            self.add_field(('O', 'O', 'O', 'O'), [PI / 2] * 4, [0.999] * 4)

    def side_rows_hook(self):
        return False

    # ---------------------------------------------------------------- interpolation
    def w1(self, l, th):
        """bilinear weights of P1 at (l, th): list of (var, w)."""
        NL, NA = self.NL, self.NA
        u = min(max(l * NL, 0.0), NL - 1e-12)
        v = min(max(th / PI * NA, 0.0), NA - 1e-12)
        i, j = int(u), int(v)
        fu, fv = u - i, v - j
        return [(self.iP1[i, j], (1 - fu) * (1 - fv)), (self.iP1[i + 1, j], fu * (1 - fv)),
                (self.iP1[i, j + 1], (1 - fu) * fv), (self.iP1[i + 1, j + 1], fu * fv)]

    def wlin(self, arr, t, N):
        u = min(max(t * N, 0.0), N - 1e-12)
        i = int(u)
        f = u - i
        return [(arr[i], 1 - f), (arr[i + 1], f)]

    def w0(self, l, th_i, th_o):
        """weights of a 0-piece (separable model): P0(l) + T0(th_I) + T0(th_O)."""
        return self.wP0(l) + self.wT(self.iT0, th_i) + self.wT(self.iT0, th_o)

    def zero_piece(self, V, l, th_i, th_o):
        NL, NA = self.NL, self.NA
        v0, du = interp1g(V['P0'], l * NL)
        v1, d1 = interp1g(V['T0'], th_i / PI * NA)
        v2, d2 = interp1g(V['T0'], th_o / PI * NA)
        return v0 + v1 + v2, du * NL, d1 * NA / PI, d2 * NA / PI

    def wP0(self, l):
        return self.wlin(self.iP0, l, self.NL)

    def wT(self, arr, th):
        return self.wlin(arr, th / PI, self.NA)

    def field_coef(self, pat, th, ell):
        st = side_types(pat)
        k = len(pat)
        c = {}
        for i in range(k):
            t = st[i]
            if t == '2':
                self.acc(c, self.iP2, 1.0)
            elif t == '1a':
                for v, w in self.w1(ell[i], th[(i + 1) % k]):
                    self.acc(c, v, w)
            elif t == '1b':
                for v, w in self.w1(ell[i], th[i]):
                    self.acc(c, v, w)
            else:
                for v, w in self.w0(ell[i], th[i], th[(i + 1) % k]):
                    self.acc(c, v, w)
        for i in range(k):
            if pat[i] == 'E':
                for v, w in self.wT(self.iE, th[i]):
                    self.acc(c, v, w)
            elif pat[i] == 'R':
                for v, w in self.wT(self.iR, th[i] - PI):
                    self.acc(c, v, w)
        return c

    def add_field(self, pat, th, ell):
        area = polygon_area(th, ell)
        c = self.field_coef(pat, th, ell)
        neg = {v: -w for v, w in c.items()}
        key = ('F', pat, tuple(np.round(th, 9)), tuple(np.round(ell, 9)))
        ok = self.add(neg, -area, 'field:' + ''.join(pat), key=key)
        if ok:
            if not hasattr(self, 'geom'):
                self.geom = {}
            self.geom[len(self.rows) - 1] = (tuple(pat), np.array(th, float), np.array(ell, float))
        return ok

    # ---------------------------------------------------------------- LP solve
    def ref_point(self):
        """centering target: the 'perimeter only' potentials (P = l/4 per piece, beta = 1/4)."""
        x = np.zeros(self.nv)
        NL, NA = self.NL, self.NA
        for i in range(NL + 1):
            x[self.iP1[i, :]] = 0.25 * i / NL
            x[self.iP0[i]] = 0.25 * i / NL
        x[self.iP2] = 0.25
        x[self.ibeta] = 0.25
        return x

    def solve(self, center=True, slack=1e-7):
        m = len(self.rows)
        I, J, V, b = [], [], [], np.zeros(m)
        for r, (cols, vals, rhs, kind) in enumerate(self.rows):
            I.extend([r] * len(cols))
            J.extend(cols.tolist())
            V.extend(vals.tolist())
            b[r] = rhs
        A = coo_matrix((V, (I, J)), shape=(m, self.nv)).tocsr()
        cvec = np.zeros(self.nv)
        cvec[self.ibeta] = 2.0
        cvec[self.irho] = 2.0
        bnds = [(-self.bound, self.bound)] * self.nv
        res = linprog(cvec, A_ub=A, b_ub=b, bounds=bnds, method='highs')
        for meth in ('highs-ipm', 'highs-ds'):
            if res.status == 0:
                break
            res = linprog(cvec, A_ub=A, b_ub=b, bounds=bnds, method=meth)
        if res.status != 0:
            raise RuntimeError(res.message)
        val = res.fun
        self.duals = -res.ineqlin.marginals
        self.res = res
        self.x = res.x
        if center:
            # phase 2: stay optimal, minimise the L1 distance of the potentials to the reference
            ref = self.ref_point()
            if getattr(self, 'center_mode', 'ref') == 'prox' and getattr(self, 'prev_x', None) is not None \
                    and len(self.prev_x) == self.nv:
                ref = self.prev_x.copy()
            pot = self.pot_idx() if hasattr(self, 'pot_idx') else np.concatenate((self.iP1.ravel(), self.iP0, self.iT0, self.iE, [self.iP2]))
            npot = len(pot)
            nv2 = self.nv + npot
            from scipy.sparse import hstack, vstack, csr_matrix, identity
            Z = csr_matrix((m, npot))
            A1 = hstack([A, Z])
            # x_p - d <= ref ; -x_p - d <= -ref
            Sel = csr_matrix((np.ones(npot), (np.arange(npot), pot)), shape=(npot, self.nv))
            Id = identity(npot, format='csr')
            A2 = hstack([Sel, -Id])
            A3 = hstack([-Sel, -Id])
            obj_row = np.zeros(nv2)
            obj_row[self.ibeta] = 2.0
            obj_row[self.irho] = 2.0
            Aall = vstack([A1, A2, A3, csr_matrix(obj_row)]).tocsr()
            ball = np.concatenate((b, ref[pot], -ref[pot], [val + slack]))
            c2 = np.concatenate((np.zeros(self.nv), np.ones(npot)))
            bnds2 = bnds + [(0, None)] * npot
            r2 = linprog(c2, A_ub=Aall, b_ub=ball, bounds=bnds2, method='highs')
            if r2.status != 0 and getattr(self, 'center_mode', 'ref') == 'prox':
                # fall back to centering on the reference point
                ball2 = np.concatenate((b, self.ref_point()[pot], -self.ref_point()[pot], [val + slack]))
                r2 = linprog(c2, A_ub=Aall, b_ub=ball2, bounds=bnds2, method='highs')
            if r2.status == 0:
                self.x = r2.x[:self.nv]
            else:
                self.center_fail = getattr(self, 'center_fail', 0) + 1
        self.prev_x = self.x.copy()
        return val

    def vals(self):
        x = self.x
        return dict(P1=x[self.iP1], P0=x[self.iP0], T0=x[self.iT0], PE=x[self.iE], PR=x[self.iR], P2=x[self.iP2],
                    rho=x[self.irho], beta=x[self.ibeta])

    # ---------------------------------------------------------------- grid separations
    def FE_table(self, cmax):
        """FE[c][s] = max sum of c E sectors (grid angles 0..NA each) with total s (0..2NA), and argmax."""
        NA = self.NA
        PE = self.x[self.iE]
        S = 2 * NA
        FE = np.full((cmax + 1, S + 1), -np.inf)
        CH = np.zeros((cmax + 1, S + 1), dtype=int)
        FE[0][0] = 0.0
        for c in range(1, cmax + 1):
            for s in range(S + 1):
                best, arg = -np.inf, -1
                for j in range(0, min(NA, s) + 1):
                    v = FE[c - 1][s - j] + PE[j]
                    if v > best:
                        best, arg = v, j
                FE[c][s] = best
                CH[c][s] = arg
        return FE, CH

    @staticmethod
    def back(CH, c, s):
        out = []
        while c > 0:
            j = CH[c][s]
            out.append(j)
            s -= j
            c -= 1
        return out

    def landing_best(self, FE, CH, left, right):
        """max over m>=1 stems and angles jO + sum(E) + jI = NA of left[jO] + FE + right[jI] - m rho.
        left/right: arrays over angle index.  Returns value and (m, jO, Es, jI)."""
        NA = self.NA
        rho = self.x[self.irho]
        best, arg = -np.inf, None
        for m in range(1, self.fan_max + 1):
            for jO in range(NA + 1):
                for jI in range(NA + 1 - jO):
                    rest = NA - jO - jI
                    if m == 1 and rest != 0:
                        continue
                    fe = FE[m - 1][rest]
                    if fe == -np.inf:
                        continue
                    v = left[jO] + fe + right[jI] - m * rho
                    if v > best:
                        best, arg = v, (m, jO, self.back(CH, m - 1, rest), jI)
        return best, arg

    def sep_grid(self, tol=1e-9):
        NL, NA = self.NL, self.NA
        x = self.x
        P1 = x[self.iP1]
        P0 = x[self.iP0]
        T0 = x[self.iT0]
        rho, beta = x[self.irho], x[self.ibeta]
        FE, CH = self.FE_table(max(self.conv_max, self.fan_max))
        added = 0
        # conv(e) junctions
        for e in range(3, self.conv_max + 1):
            v = FE[e][2 * NA] - e * rho
            if v > tol:
                js = self.back(CH, e, 2 * NA)
                c = {}
                for j in js:
                    self.acc(c, self.iE[j], 1.0)
                self.acc(c, self.irho, -float(e))
                if self.add(c, 0.0, 'conv', key=('conv', tuple(sorted(js)))):
                    added += 1
        # k = 1 sides with fans (m >= 2); m = 1 rows are static
        for i in range(NL + 1):
            best, arg = self.landing_best(FE, CH, P1[i], P1[NL - i])
            if best - beta > tol and arg[0] >= 2:
                m, jO, Es, jI = arg
                c = {}
                self.acc(c, self.iP1[i, jO], 1.0)
                self.acc(c, self.iP1[NL - i, jI], 1.0)
                for j in Es:
                    self.acc(c, self.iE[j], 1.0)
                self.acc(c, self.ibeta, -1.0)
                self.acc(c, self.irho, -float(m))
                if self.add(c, 0.0, 'fan1', key=('fan1', i, m, jO, tuple(sorted(Es)), jI)):
                    added += 1
        # chains k >= 2
        Ast, Aen = [], []
        for i in range(NL + 1):
            Ast.append(self.landing_best(FE, CH, P1[i], T0))
            Aen.append(self.landing_best(FE, CH, T0, P1[i]))
        Amid = self.landing_best(FE, CH, T0, T0)
        # V[s]: best partial chain after a 0-piece, total length index s; track structure
        V = np.full(NL + 1, -np.inf)
        VA = [None] * (NL + 1)
        for i0 in range(NL + 1):
            for i1 in range(NL + 1 - i0):
                v = Ast[i0][0] + P0[i1]
                if v > V[i0 + i1]:
                    V[i0 + i1] = v
                    VA[i0 + i1] = [('start', i0, Ast[i0][1]), ('p0', i1)]
        for it in range(8):
            changed = False
            for s in range(NL + 1):
                if V[s] == -np.inf:
                    continue
                for t in range(NL + 1 - s):
                    v = V[s] + Amid[0] + P0[t]
                    if v > V[s + t] + 1e-15:
                        V[s + t] = v
                        VA[s + t] = VA[s] + [('mid', Amid[1]), ('p0', t)]
                        changed = True
            if not changed:
                break
        best = -np.inf
        cands = []
        for s in range(NL + 1):
            if V[s] == -np.inf:
                continue
            ie = NL - s
            v = V[s] + Aen[ie][0]
            best = max(best, v)
            if v - beta > tol:
                cands.append((v, VA[s] + [('end', ie, Aen[ie][1])]))
        cands.sort(key=lambda t: -t[0])
        for v, barg in cands[:self.chain_rows]:
            if len(barg) >= 40:
                continue
            c = {}
            nst = 0
            for item in barg:
                if item[0] == 'start':
                    _, i0, (m, jO, Es, jI) = item
                    self.acc(c, self.iP1[i0, jO], 1.0)
                    self.acc(c, self.iT0[jI], 1.0)
                elif item[0] == 'end':
                    _, ie, (m, jO, Es, jI) = item
                    self.acc(c, self.iT0[jO], 1.0)
                    self.acc(c, self.iP1[ie, jI], 1.0)
                elif item[0] == 'mid':
                    (m, jO, Es, jI) = item[1]
                    self.acc(c, self.iT0[jO], 1.0)
                    self.acc(c, self.iT0[jI], 1.0)
                else:
                    self.acc(c, self.iP0[item[1]], 1.0)
                    continue
                for j in Es:
                    self.acc(c, self.iE[j], 1.0)
                nst += m
            self.acc(c, self.ibeta, -1.0)
            self.acc(c, self.irho, -float(nst))
            key = ('chain', tuple(sorted((k, round(v, 12)) for k, v in c.items())))
            if self.add(c, 0.0, 'chain', key=key):
                added += 1
        worst = max([FE[e][2 * NA] - e * rho for e in range(3, self.conv_max + 1)] or [-np.inf])
        return added, dict(conv=worst, chain=best - beta)

    def grid_check(self):
        """independent DP evaluation of the worst fence-side / fan / chain / conv violation."""
        NL, NA = self.NL, self.NA
        x = self.x
        P1, P0, T0 = x[self.iP1], x[self.iP0], x[self.iT0]
        rho, beta = x[self.irho], x[self.ibeta]
        FE, CH = self.FE_table(max(self.conv_max, self.fan_max))
        conv = max([FE[e][2 * NA] - e * rho for e in range(3, self.conv_max + 1)] or [-np.inf])
        side1 = max(self.landing_best(FE, CH, P1[i], P1[NL - i])[0] for i in range(NL + 1)) - beta
        Ast = [self.landing_best(FE, CH, P1[i], T0)[0] for i in range(NL + 1)]
        Amid = self.landing_best(FE, CH, T0, T0)[0]
        V = np.full(NL + 1, -np.inf)
        for i0 in range(NL + 1):
            for i1 in range(NL + 1 - i0):
                V[i0 + i1] = max(V[i0 + i1], Ast[i0] + P0[i1])
        for it in range(NL + 2):
            for s in range(NL + 1):
                for t in range(NL + 1 - s):
                    V[s + t] = max(V[s + t], V[s] + Amid + P0[t])
        chain = max(V[s] + Ast[NL - s] for s in range(NL + 1)) - beta
        refl = -np.inf
        if 'R' in self.alphabet:
            PR = x[self.iR]
            for e in range(2, self.refl_max + 1):
                for s in range(NA + 1):
                    refl = max(refl, FE[e - 1][s] + PR[NA - s] - e * rho)
        return dict(conv=max(conv, refl), chain=max(chain, side1, x[self.iP2] - beta))

    def seed_library(self, M=40):
        """add field rows for M random closed typed polygons per pattern (area <= 1)."""
        added = 0
        for pat in self.pats:
            k = len(pat)
            st = side_types(pat)
            free = [i for i in range(k) if st[i] != '2']
            if len(free) < 2:
                continue
            got = 0
            for trial in range(40 * M):
                if got >= M:
                    break
                phi = sample_ext(self.rng, pat, 2.0)
                if phi is None:
                    continue
                psi = np.concatenate(([0.0], np.cumsum(phi[1:])))
                u = np.stack((np.cos(psi), np.sin(psi)), axis=1)
                ell = np.ones(k)
                a2, b2 = self.rng.choice(free, 2, replace=False)
                others = [i for i in free if i not in (a2, b2)]
                ell[others] = self.rng.uniform(0.02, 0.98, len(others))
                rhs = -sum(ell[i] * u[i] for i in range(k) if i not in (a2, b2))
                Mx = np.stack((u[a2], u[b2]), axis=1)
                if abs(np.linalg.det(Mx)) < 1e-9:
                    continue
                sol = np.linalg.solve(Mx, rhs)
                ell[a2], ell[b2] = sol
                if ell[free].min() <= 1e-3 or ell[free].max() >= 1 - 1e-3:
                    continue
                if not free or len(free) == k:
                    # pinwheel: scale up towards area 1 while all sides stay < 1
                    A0 = polygon_area(PI - phi, ell)
                    if A0 <= 0:
                        continue
                    s = min(1 / math.sqrt(A0), 0.999 / ell.max())
                    ell = ell * s
                A = polygon_area(PI - phi, ell)
                if A <= 0 or A > 1:
                    continue
                if 'R' in pat and not is_simple(PI - phi, ell):
                    continue
                if self.add_field(pat, PI - phi, ell):
                    added += 1
                    got += 1
        return added

    # ---------------------------------------------------------------- field pricing
    def pot_eval(self, pat, th, ell, V):
        st = side_types(pat)
        k = len(pat)
        s = 0.0
        for i in range(k):
            t = st[i]
            if t == '2':
                s += V['P2']
            elif t == '1a':
                s += interp2(V['P1'], ell[i] * self.NL, th[(i + 1) % k] / PI * self.NA)
            elif t == '1b':
                s += interp2(V['P1'], ell[i] * self.NL, th[i] / PI * self.NA)
            else:
                s += self.zero_piece(V, ell[i], th[i], th[(i + 1) % k])[0]
        for i in range(k):
            if pat[i] == 'E':
                s += interp1(V['PE'], th[i] / PI * self.NA)
            elif pat[i] == 'R':
                s += interp1(V['PR'], (th[i] - PI) / PI * self.NA)
        return s

    def price_pattern(self, pat, V, starts=6, tol=1e-7):
        k = len(pat)
        st = side_types(pat)
        free = [i for i in range(k) if st[i] != '2']
        nf = len(free)
        eps_a = 1e-3
        NL, NA = self.NL, self.NA
        isE = [c == 'E' for c in pat]
        isR = [c == 'R' for c in pat]

        def unpack(z):
            phi = np.empty(k)
            phi[1:] = z[:k - 1]
            phi[0] = 2 * PI - z[:k - 1].sum()
            ell = np.ones(k)
            ell[free] = z[k - 1:]
            return phi, ell

        def geom(z):
            phi, ell = unpack(z)
            psi = np.concatenate(([0.0], np.cumsum(phi[1:])))
            return phi, ell, psi

        def closure(z):
            phi, ell, psi = geom(z)
            return np.array([np.dot(ell, np.cos(psi)), np.dot(ell, np.sin(psi))])

        def closure_jac(z):
            phi, ell, psi = geom(z)
            J = np.zeros((2, k - 1 + nf))
            s, c = np.sin(psi), np.cos(psi)
            for m in range(1, k):
                J[0, m - 1] = -np.dot(ell[m:], s[m:])
                J[1, m - 1] = np.dot(ell[m:], c[m:])
            J[0, k - 1:] = c[free]
            J[1, k - 1:] = s[free]
            return J

        def area_grad(z):
            phi, ell, psi = geom(z)
            D = psi[:, None] - psi[None, :]
            M = np.tril(np.sin(D), -1)
            A = 0.5 * ell @ M @ ell
            dl = 0.5 * (M + M.T) @ ell
            W = np.tril(np.cos(D), -1) * np.outer(ell, ell)
            dphi = np.zeros(k - 1)
            for m in range(1, k):
                dphi[m - 1] = 0.5 * W[m:, :m].sum()
            g = np.concatenate((dphi, dl[free]))
            return A, g

        def pot_grad(z):
            phi, ell = unpack(z)
            th = PI - phi
            dth = np.zeros(k)
            dl = np.zeros(k)
            s = 0.0
            for i in range(k):
                t = st[i]
                if t == '2':
                    s += V['P2']
                elif t == '1a' or t == '1b':
                    j = (i + 1) % k if t == '1a' else i
                    val, du, dv = interp2g(V['P1'], ell[i] * NL, th[j] / PI * NA)
                    s += val
                    dl[i] += du * NL
                    dth[j] += dv * NA / PI
                else:
                    j2 = (i + 1) % k
                    val, du, d1, d2 = self.zero_piece(V, ell[i], th[i], th[j2])
                    s += val
                    dl[i] += du
                    dth[i] += d1
                    dth[j2] += d2
            for i in range(k):
                if isE[i]:
                    val, dv = interp1g(V['PE'], th[i] / PI * NA)
                    s += val
                    dth[i] += dv * NA / PI
                elif isR[i]:
                    val, dv = interp1g(V['PR'], (th[i] - PI) / PI * NA)
                    s += val
                    dth[i] += dv * NA / PI
            # theta_m = pi - phi_m (m >= 1), theta_0 = sum_{m>=1} phi_m - pi
            dphi = -dth[1:] + dth[0]
            return s, np.concatenate((dphi, dl[free]))

        def obj(z):
            A, gA = area_grad(z)
            p, gp = pot_grad(z)
            return -(A - p), -(gA - gp)

        ones = np.concatenate((np.ones(k - 1), np.zeros(nf)))
        # exterior angle phi_0 = 2pi - sum(z[:k-1]) must lie in (lo0, hi0)
        lo0, hi0 = ((-PI + eps_a, -eps_a) if isR[0] else (eps_a, PI - eps_a))
        cons = [{'type': 'eq', 'fun': closure, 'jac': closure_jac},
                {'type': 'ineq', 'fun': lambda z: 1.0 - area_grad(z)[0], 'jac': lambda z: -area_grad(z)[1]},
                {'type': 'ineq', 'fun': lambda z: z[:k - 1].sum() - (2 * PI - hi0), 'jac': lambda z: ones},
                {'type': 'ineq', 'fun': lambda z: (2 * PI - lo0) - z[:k - 1].sum(), 'jac': lambda z: -ones}]
        bnds = [((-PI + eps_a, -eps_a) if isR[i] else (eps_a, PI - eps_a)) for i in range(1, k)] + \
            [(1e-4, 1 - 1e-4)] * nf
        best = None
        inits = list(self.pool.get(pat, []))
        for _ in range(starts):
            phi = sample_ext(self.rng, pat, 3.0)
            if phi is None:
                continue
            ell0 = self.rng.uniform(0.2, 0.95, nf)
            inits.append(np.concatenate((phi[1:], ell0)))
        for z0 in inits:
            try:
                with np.errstate(all='ignore'):
                    r = minimize(obj, z0, jac=True, method='SLSQP', bounds=bnds, constraints=cons,
                                 options={'maxiter': 200, 'ftol': 1e-12})
            except Exception:
                continue
            z = np.clip(r.x, [b[0] for b in bnds], [b[1] for b in bnds])
            if np.max(np.abs(closure(z))) > 1e-8:
                continue
            A = area_grad(z)[0]
            if A > 1 + 1e-9 or A <= 0:
                continue
            phi, ell = unpack(z)
            if any((phi[i] >= 0 or phi[i] <= -PI) if isR[i] else (phi[i] <= 0 or phi[i] >= PI) for i in range(k)):
                continue
            if any(isR) and not is_simple(PI - phi, ell):
                continue
            val = -obj(z)[0]
            if val > 1e-7:
                if not hasattr(self, '_extra'):
                    self._extra = []
                self._extra.append((val, z.copy(), PI - phi, ell.copy()))
            if best is None or val > best[0]:
                best = (val, z, PI - phi, ell)
        if best is not None:
            lst = self.pool.setdefault(pat, [])
            lst.append(best[1])
            if len(lst) > 3:
                lst.pop(0)
        return best

    def sep_fields(self, starts=4, tol=1e-7, pats=None):
        V = self.vals()
        added = 0
        worst = -np.inf
        wpat = None
        if not hasattr(self, 'last'):
            self.last = {}
        for pat in (pats if pats is not None else self.pats):
            b = self.price_pattern(pat, V, starts=starts)
            if b is None:
                continue
            val, z, th, ell = b
            self.last[pat] = val
            if val > worst:
                worst, wpat = val, pat
            if val > tol:
                if self.add_field(pat, th, ell):
                    added += 1
            # further distinct violated local optima of this pattern (up to ncuts - 1)
            ex = sorted(getattr(self, '_extra', []), key=lambda t: -t[0])
            self._extra = []
            nadd = 0
            for v2, z2, th2, ell2 in ex:
                if nadd >= getattr(self, 'ncuts', 1) - 1:
                    break
                if np.max(np.abs(ell2 - ell)) + np.max(np.abs(th2 - th)) < 1e-3:
                    continue
                if v2 > tol and self.add_field(pat, th2, ell2):
                    added += 1
                    nadd += 1
        return added, worst, wpat

    def purge(self, slack_min=0.02):
        """drop field rows (never the grid rows, never the first field row) whose slack at the
        current LP solution exceeds slack_min and whose dual weight is zero; keys are released so
        that a purged row can be re-added later.  Returns the number of rows removed."""
        x = self.x
        geom = getattr(self, 'geom', {})
        duals = getattr(self, 'duals', None)
        keep = []
        first_field = None
        for i, (cols, vals, rhs, kind) in enumerate(self.rows):
            if not kind.startswith('field'):
                keep.append(True)
                continue
            if first_field is None:
                first_field = i
                keep.append(True)
                continue
            slack = rhs - float(np.dot(vals, x[cols]))
            w = duals[i] if duals is not None and i < len(duals) else 0.0
            keep.append(not (slack > slack_min and w <= 1e-12))
        newrows, newgeom, removed = [], {}, 0
        for i, r in enumerate(self.rows):
            if keep[i]:
                if i in geom:
                    newgeom[len(newrows)] = geom[i]
                newrows.append(r)
            else:
                removed += 1
                if i in geom:
                    pat, th, ell = geom[i]
                    self.keys.discard(('F', pat, tuple(np.round(th, 9)), tuple(np.round(ell, 9))))
        self.rows = newrows
        self.geom = newgeom
        self.duals = None
        return removed

    def run(self, iters=60, starts=4, verbose=True, tol=1e-7, full_every=5, act=0.03, min_iters=0):
        hist = []
        for it in range(iters):
            val = self.solve()
            ag, gw = 0, self.grid_check()
            full = (it % full_every == 0)
            if full or not hasattr(self, 'last'):
                pats = None
            else:
                pats = [p for p in self.pats if self.last.get(p, 1.0) > -act]
            af, fw, fp = self.sep_fields(starts=starts if full else 1, tol=tol, pats=pats)
            hist.append((val, ag, af, fw))
            if verbose:
                print(f'it {it:3d} bound {val:.7f} rows {len(self.rows)} grid+{ag} fields+{af} '
                      f'{"FULL" if full else len(pats)} worst {fw:.2e} {"".join(fp) if fp else ""} '
                      f'chain {gw["chain"]:.2e} conv {gw["conv"]:.2e}', flush=True)
            if ag == 0 and af == 0 and full and it >= min_iters:
                break
        return hist


def polygon_area(th, ell):
    k = len(th)
    phi = PI - np.asarray(th)
    psi = np.cumsum(np.concatenate(([0.0], phi[1:])))
    xs = np.concatenate(([0.0], np.cumsum(ell * np.cos(psi))))
    ys = np.concatenate(([0.0], np.cumsum(ell * np.sin(psi))))
    return 0.5 * float(np.sum(xs[:-1] * ys[1:] - xs[1:] * ys[:-1]))


def interp1(arr, u):
    N = len(arr) - 1
    u = min(max(u, 0.0), N - 1e-12)
    i = int(u)
    f = u - i
    return (1 - f) * arr[i] + f * arr[i + 1]


def interp2(arr, u, v):
    NL, NA = arr.shape[0] - 1, arr.shape[1] - 1
    u = min(max(u, 0.0), NL - 1e-12)
    v = min(max(v, 0.0), NA - 1e-12)
    i, j = int(u), int(v)
    fu, fv = u - i, v - j
    return ((1 - fu) * (1 - fv) * arr[i, j] + fu * (1 - fv) * arr[i + 1, j] +
            (1 - fu) * fv * arr[i, j + 1] + fu * fv * arr[i + 1, j + 1])


def interp1g(arr, u):
    N = len(arr) - 1
    u = min(max(u, 0.0), N - 1e-12)
    i = int(u)
    f = u - i
    return (1 - f) * arr[i] + f * arr[i + 1], arr[i + 1] - arr[i]


def interp2g(arr, u, v):
    NL, NA = arr.shape[0] - 1, arr.shape[1] - 1
    u = min(max(u, 0.0), NL - 1e-12)
    v = min(max(v, 0.0), NA - 1e-12)
    i, j = int(u), int(v)
    fu, fv = u - i, v - j
    a, b, c, d = arr[i, j], arr[i + 1, j], arr[i, j + 1], arr[i + 1, j + 1]
    val = (1 - fu) * (1 - fv) * a + fu * (1 - fv) * b + (1 - fu) * fv * c + fu * fv * d
    du = (1 - fv) * (b - a) + fv * (d - c)
    dv = (1 - fu) * (c - a) + fu * (d - b)
    return val, du, dv


def sample_ext(rng, pat, conc):
    """random exterior angles (phi = pi - theta) with the sign pattern of pat (R: reflex), sum 2pi."""
    k = len(pat)
    isR = np.array([c == 'R' for c in pat])
    nr = int(isR.sum())
    phi = np.empty(k)
    neg = rng.uniform(0.05, 0.6 * PI, nr) if nr else np.zeros(0)
    tot = 2 * PI + neg.sum()
    pos = rng.dirichlet(np.ones(k - nr) * conc) * tot
    if pos.max() >= PI - 1e-3 or pos.min() < 1e-3:
        return None
    phi[isR] = -neg
    phi[~isR] = pos
    return phi


def is_simple(th, ell):
    """the closed polygon with interior angles th and sides ell has no self-intersection."""
    k = len(th)
    phi = PI - np.asarray(th)
    psi = np.concatenate(([0.0], np.cumsum(phi[1:])))
    P = np.zeros((k + 1, 2))
    P[1:, 0] = np.cumsum(ell * np.cos(psi))
    P[1:, 1] = np.cumsum(ell * np.sin(psi))

    def cr(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    for i in range(k):
        for j in range(i + 1, k):
            if j == i + 1 or (i == 0 and j == k - 1):
                continue
            a, b, c, d = P[i], P[i + 1], P[j], P[j + 1]
            d1, d2, d3, d4 = cr(c, d, a), cr(c, d, b), cr(a, b, c), cr(a, b, d)
            if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)):
                return False
    return True


# ------------------------------------------------------------------------------ certificates
def export_cert(L, path, extra_polys=()):
    """save the potentials (grid functions) and warm-start polygons as JSON."""
    import json
    x = L.x
    polys = []
    geom = getattr(L, 'geom', {})
    for r, w in enumerate(getattr(L, 'duals', [])):
        if w > 1e-7 and r in geom:
            pat, th, ell = geom[r]
            polys.append([''.join(pat), list(map(float, th)), list(map(float, ell)), float(w)])
    for pat, zs in L.pool.items():
        pass
    for p in extra_polys:
        polys.append(p)
    cert = dict(NL=L.NL, NA=L.NA, kmax=L.kmax, fan_max=L.fan_max, conv_max=L.conv_max,
                alphabet=L.alphabet, rmax=L.rmax, refl_max=L.refl_max, kmax_r=getattr(L, 'kmax_r', None),
                P1=x[L.iP1].tolist(), P0=x[L.iP0].tolist(), T0=x[L.iT0].tolist(), PE=x[L.iE].tolist(),
                PR=x[L.iR].tolist(), P2=float(x[L.iP2]), rho=float(x[L.irho]), beta=float(x[L.ibeta]),
                polys=polys, cls=type(L).__name__)
    if hasattr(L, 'iP0N'):
        cert['P0N'] = x[L.iP0N].tolist()
    with open(path, 'w') as f:
        json.dump(cert, f)
    return cert


def load_cert(path, seed=0):
    import json
    with open(path) as f:
        c = json.load(f)
    cls = TypedLP2 if c.get('cls') == 'TypedLP2' else TypedLP
    L = cls(NL=c['NL'], NA=c['NA'], kmax=c['kmax'], fan_max=c['fan_max'], conv_max=c['conv_max'],
            seed=seed, alphabet=c['alphabet'], rmax=c['rmax'], refl_max=c['refl_max'], kmax_r=c.get('kmax_r'))
    x = np.zeros(L.nv)
    x[L.iP1] = np.array(c['P1'])
    x[L.iP0] = np.array(c['P0'])
    x[L.iT0] = np.array(c['T0'])
    x[L.iE] = np.array(c['PE'])
    x[L.iR] = np.array(c['PR'])
    x[L.iP2] = c['P2']
    x[L.irho] = c['rho']
    x[L.ibeta] = c['beta']
    if 'P0N' in c:
        x[L.iP0N] = np.array(c['P0N'])
    L.x = x
    for pat, th, ell, w in c['polys']:
        z = np.concatenate(((PI - np.array(th))[1:], np.array(ell)[[i for i, t in enumerate(side_types(tuple(pat))) if t != '2']]))
        L.pool.setdefault(tuple(pat), []).append(z)
    return L, c


def verify_cert(L, starts=10, extra_seeds=0):
    """independent re-check: exact DP on the grid families; multistart pricing on all field patterns.
    Returns (value 2(beta+rho), grid violation, max field violation, worst pattern, repaired bound)."""
    val = 2 * (L.x[L.ibeta] + L.x[L.irho])
    g = L.grid_check()
    gridv = max(g['conv'], g['chain'])
    V = L.vals()
    worst, wpat = -np.inf, None
    for pat in L.pats:
        b = L.price_pattern(pat, V, starts=starts)
        if b is not None and b[0] > worst:
            worst, wpat = b[0], pat
    rep = val + 4.0 / 3.0 * max(worst, 0.0) + 4.0 / 3.0 * max(gridv, 0.0)
    return val, gridv, worst, wpat, rep


def tail_check(L, emax=12, kmax_extra=7, starts=4):
    """tails outside the LP: conv(e) and fans with up to emax ends/stems (exact DP on the grid), and
    field patterns with kappa = kmax+1 .. kmax_extra convex corners (multistart pricing)."""
    NA = L.NA
    rho = L.x[L.irho]
    FE, CH = L.FE_table(emax)
    conv = max(FE[e][2 * NA] - e * rho for e in range(3, emax + 1))
    # fans: the LP used K[r] = max_{c <= fan_max-1} (FE[c][r] - c rho) for the E sectors of a landing;
    # fanE = max_r (K over c <= emax-1) - (K over c <= fan_max-1): the extra gain of larger fans
    Kc = np.array([max([0.0 if r == 0 else -np.inf] + [FE[c][r] - c * rho for c in range(1, L.fan_max)]) for r in range(NA + 1)])
    Ka = np.array([max([0.0 if r == 0 else -np.inf] + [FE[c][r] - c * rho for c in range(1, emax)]) for r in range(NA + 1)])
    with np.errstate(invalid='ignore'):
        dif = np.where(np.isfinite(Ka), Ka - np.where(np.isfinite(Kc), Kc, -1e9), -np.inf)
    fanE = float(np.max(dif)) if L.fan_max > 1 else 0.0
    V = L.vals()
    worst, wpat = -np.inf, None
    for pat in patterns(L.kmax + 1, kmax_extra, L.alphabet, 0):
        b = L.price_pattern(pat, V, starts=starts)
        if b is not None and b[0] > worst:
            worst, wpat = b[0], pat
    return dict(conv=conv, fanE=fanE, field=worst, pat=wpat)


class TypedLP2(TypedLP):
    """typed piece-potential LP with a NON-separable 0-piece potential P0N[l, th_I, th_O]
    (length and both end angles of a 0 piece linked).  Fence sides are encoded by the chain DP
    Z[s][jO] (partial side ending with a piece whose end O-angle is jO, total length s) and
    Y[s][jI] (after a landing whose I-angle is jI)."""

    def side_rows_hook(self):
        NL, NA = self.NL, self.NA
        A = self.acc
        idx = self.nv
        self.iP0N = np.arange(idx, idx + (NL + 1) * (NA + 1) * (NA + 1)).reshape(NL + 1, NA + 1, NA + 1)
        idx += (NL + 1) * (NA + 1) * (NA + 1)
        self.iY = np.arange(idx, idx + (NL + 1) * (NA + 1)).reshape(NL + 1, NA + 1)
        idx += (NL + 1) * (NA + 1)
        self.iZ2 = np.arange(idx, idx + (NL + 1) * (NA + 1)).reshape(NL + 1, NA + 1)
        idx += (NL + 1) * (NA + 1)
        self.nv = idx
        for i in range(NL + 1):
            for jO in range(NA + 1):
                c = {}
                A(c, self.iP1[i, jO], 1.0)
                A(c, self.iZ2[i, jO], -1.0)
                self.add(c, 0.0, 'aux')
        for s in range(NL + 1):
            for jO in range(NA + 1):
                for jI in range(NA + 1 - jO):
                    c = {}
                    A(c, self.iZ2[s, jO], 1.0)
                    A(c, self.iK[NA - jO - jI], 1.0)
                    A(c, self.irho, -1.0)
                    A(c, self.iY[s, jI], -1.0)
                    self.add(c, 0.0, 'aux')
        for s in range(NL + 1):
            for t in range(NL + 1 - s):
                for jI in range(NA + 1):
                    for jO in range(NA + 1):
                        c = {}
                        A(c, self.iY[s, jI], 1.0)
                        A(c, self.iP0N[t, jI, jO], 1.0)
                        A(c, self.iZ2[s + t, jO], -1.0)
                        self.add(c, 0.0, 'aux')
        for s in range(NL + 1):
            for jI in range(NA + 1):
                c = {}
                A(c, self.iY[s, jI], 1.0)
                A(c, self.iP1[NL - s, jI], 1.0)
                A(c, self.ibeta, -1.0)
                self.add(c, 0.0, 'side')
        c = {}
        A(c, self.iP2, 1.0)
        A(c, self.ibeta, -1.0)
        self.add(c, 0.0, 'side0')
        # mirror symmetry: P0N[l, a, b] = P0N[l, b, a].  The mirror image of a certificate is a
        # certificate (the configuration class is mirror invariant), so symmetric certificates lose
        # nothing, and with them a pattern and its mirror image give the same field constraint
        # (patterns are enumerated up to mirror) and the chain constraints are mirror invariant.
        for t in range(NL + 1):
            for a in range(NA + 1):
                for b in range(a + 1, NA + 1):
                    for sg in (1.0, -1.0):
                        c = {}
                        A(c, self.iP0N[t, a, b], sg)
                        A(c, self.iP0N[t, b, a], -sg)
                        self.add(c, 0.0, 'sym')
        if 'E' in self.alphabet:
            self.add_field(('E', 'E', 'E', 'E'), [PI / 2] * 4, [1.0] * 4)
        else:
            self.add_field(('O', 'O', 'O', 'O'), [PI / 2] * 4, [0.999] * 4)
        return True

    def pot_idx(self):
        return np.concatenate((self.iP1.ravel(), self.iP0N.ravel(), self.iE, self.iR, [self.iP2]))

    def ref_point(self):
        x = np.zeros(self.nv)
        NL = self.NL
        for i in range(NL + 1):
            x[self.iP1[i, :]] = 0.25 * i / NL
            x[self.iP0N[i, :, :]] = 0.25 * i / NL
        x[self.iP2] = 0.25
        x[self.ibeta] = 0.25
        return x

    def vals(self):
        V = TypedLP.vals(self)
        V['P0N'] = self.x[self.iP0N]
        return V

    def w0(self, l, th_i, th_o):
        NL, NA = self.NL, self.NA
        u = min(max(l * NL, 0.0), NL - 1e-12)
        v = min(max(th_i / PI * NA, 0.0), NA - 1e-12)
        w = min(max(th_o / PI * NA, 0.0), NA - 1e-12)
        i, j, k = int(u), int(v), int(w)
        fu, fv, fw = u - i, v - j, w - k
        out = []
        for a, wa in ((i, 1 - fu), (i + 1, fu)):
            for b, wb in ((j, 1 - fv), (j + 1, fv)):
                for cc, wc in ((k, 1 - fw), (k + 1, fw)):
                    out.append((self.iP0N[a, b, cc], wa * wb * wc))
        return out

    def zero_piece(self, V, l, th_i, th_o):
        NL, NA = self.NL, self.NA
        T = V['P0N']
        u = min(max(l * NL, 0.0), NL - 1e-12)
        v = min(max(th_i / PI * NA, 0.0), NA - 1e-12)
        w = min(max(th_o / PI * NA, 0.0), NA - 1e-12)
        i, j, k = int(u), int(v), int(w)
        fu, fv, fw = u - i, v - j, w - k
        C = T[i:i + 2, j:j + 2, k:k + 2]
        wu = np.array([1 - fu, fu])
        wv = np.array([1 - fv, fv])
        ww = np.array([1 - fw, fw])
        val = np.einsum('abc,a,b,c->', C, wu, wv, ww)
        du = np.einsum('abc,a,b,c->', C, np.array([-1.0, 1.0]), wv, ww) * NL
        dv = np.einsum('abc,a,b,c->', C, wu, np.array([-1.0, 1.0]), ww) * NA / PI
        dw = np.einsum('abc,a,b,c->', C, wu, wv, np.array([-1.0, 1.0])) * NA / PI
        return val, du, dv, dw

    def grid_check(self):
        NL, NA = self.NL, self.NA
        x = self.x
        P1, P0N = x[self.iP1], x[self.iP0N]
        rho, beta = x[self.irho], x[self.ibeta]
        FE, CH = self.FE_table(max(self.conv_max, self.fan_max))
        conv = max([FE[e][2 * NA] - e * rho for e in range(3, self.conv_max + 1)] or [-np.inf])
        # K[r] = max_c (FE[c][r] - c rho), c = 0..fan_max-1
        K = np.full(NA + 1, -np.inf)
        K[0] = 0.0
        for cc in range(1, self.fan_max):
            K = np.maximum(K, FE[cc][:NA + 1] - cc * rho)
        Z = P1.copy()                               # Z[s][jO]
        Y = np.full((NL + 1, NA + 1), -np.inf)
        for it in range(NL + 3):
            Ynew = np.full((NL + 1, NA + 1), -np.inf)
            for jO in range(NA + 1):
                for jI in range(NA + 1 - jO):
                    Ynew[:, jI] = np.maximum(Ynew[:, jI], Z[:, jO] + K[NA - jO - jI] - rho)
            Znew = P1.copy()
            for s in range(NL + 1):
                for t in range(NL + 1 - s):
                    # Znew[s+t][jO] >= max_jI Y[s][jI] + P0N[t][jI][jO]
                    Znew[s + t] = np.maximum(Znew[s + t], np.max(Ynew[s][:, None] + P0N[t], axis=0))
            if np.allclose(Znew, Z) and np.allclose(np.nan_to_num(Ynew, neginf=-1e9), np.nan_to_num(Y, neginf=-1e9)):
                Y, Z = Ynew, Znew
                break
            Y, Z = Ynew, Znew
        side = max(np.max(Y[s] + P1[NL - s]) for s in range(NL + 1)) - beta
        refl = -np.inf
        if 'R' in self.alphabet:
            PR = x[self.iR]
            for e in range(2, self.refl_max + 1):
                for s in range(NA + 1):
                    refl = max(refl, FE[e - 1][s] + PR[NA - s] - e * rho)
        return dict(conv=max(conv, refl), chain=max(side, x[self.iP2] - beta))
